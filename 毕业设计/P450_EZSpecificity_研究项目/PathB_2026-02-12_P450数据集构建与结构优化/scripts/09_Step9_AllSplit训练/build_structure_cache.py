"""
Offline structure cache builder for EZSpecificity Step 9 (multiprocessing).

Iterates valid samples from the original StructureDataset (transform=None),
runs BuildStructureCacheData to extract split real/knn partitions, and writes
them to LMDB.  One LMDB per dataset_id.

SAFETY:
- Processes one dataset at a time (sequential LMDB writes).
- Windows LMDB preallocates full map_size — per-dataset sizing + disk check.
- Workers only READ source LMDBs; main process is the sole LMDB writer.
- OMP/MKL/torch threads capped at 1 per process to prevent CPU oversubscription.
- Bounded queues to prevent RAM explosion.

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe "../毕业设计/.../build_structure_cache.py" --split all
"""
from __future__ import annotations

# --- Cap BLAS/OMP threads BEFORE numpy/torch import ---
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import argparse
import math
import multiprocessing as mp
import pickle
import queue
import shutil
import sys
import time
import warnings

import lmdb
import yaml
from easydict import EasyDict
from tqdm.auto import tqdm

# ---------------------------------------------------------------------------
# 1. Monkey-patch lmdb.open for large databases (read-only opens only)
# ---------------------------------------------------------------------------
_ORIGINAL_LMDB_OPEN = lmdb.open
_READ_MAP_SIZE = 200 * 1024 ** 3


def _patched_lmdb_open(path, **kwargs):
    if kwargs.get("readonly", False) or not kwargs.get("create", True):
        kwargs["map_size"] = max(kwargs.get("map_size", 0), _READ_MAP_SIZE)
    return _ORIGINAL_LMDB_OPEN(path, **kwargs)


lmdb.open = _patched_lmdb_open

# ---------------------------------------------------------------------------
# 2. Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", "src"))
PATHB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))

for d in (SCRIPT_DIR, SRC_DIR):
    if d not in sys.path:
        sys.path.insert(0, d)

from cache_utils import BuildStructureCacheData  # noqa: E402
from Datasets.Structure.structure import StructureDataset  # noqa: E402
from Datasets.utils import read_datasets  # noqa: E402


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Build structure cache LMDBs")
    p.add_argument("--config", default=os.path.join(SCRIPT_DIR, "train_allsplit_config.yml"))
    p.add_argument("--output-dir", default=None, help="Cache output directory")
    p.add_argument("--split", choices=["train", "val", "test", "all"], default="all")
    p.add_argument("--commit-interval", type=int, default=500)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--num-workers", type=int, default=2,
                   help="Worker processes for parallel structure building (default: 2)")
    p.add_argument("--chunk-size", type=int, default=4,
                   help="Samples per worker task (default: 4)")
    p.add_argument("--queue-size", type=int, default=4,
                   help="Max in-flight result batches (default: 4)")
    return p.parse_args()


def load_config(path: str) -> EasyDict:
    with open(path, "r", encoding="utf-8") as f:
        return EasyDict(yaml.safe_load(f))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _get_split_df(config, split: str):
    import pandas as pd
    frames = []
    if split in ("train", "all"):
        frames.append(read_datasets(config.data.train_data_path))
    if split in ("val", "all"):
        frames.append(read_datasets(config.data.val_data_path))
    if split in ("test", "all"):
        frames.append(read_datasets(config.data.test_data_path))
    return pd.concat(frames, ignore_index=True).reset_index(drop=True)


def _fmt_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} PB"


def _default_cache_dir() -> str:
    return os.path.normpath(os.path.join(SRC_DIR, "..", "data", "step9_structure_cache"))


def _dataset_names(config) -> list[str]:
    paths = config.data.structure_processed_path
    if isinstance(paths, str):
        paths = [paths]
    names = []
    for i, p in enumerate(paths):
        try:
            name = os.path.basename(os.path.dirname(os.path.dirname(p)))
            if not name:
                name = f"dataset{i}"
        except Exception:
            name = f"dataset{i}"
        names.append(name)
    return names


def _close_handles(dataset: StructureDataset):
    for attr in ("complex_dbs", "ligand_dbs"):
        dbs = getattr(dataset, attr, None)
        if dbs is not None:
            for db in dbs:
                if db is not None and hasattr(db, "close"):
                    try:
                        db.close()
                    except Exception:
                        pass
            setattr(dataset, attr, None)


def _scan_existing_keys(path: str) -> set[bytes]:
    if not os.path.exists(path):
        return set()
    env = _ORIGINAL_LMDB_OPEN(path, map_size=_READ_MAP_SIZE, create=False,
                               subdir=False, readonly=True, lock=False,
                               readahead=False, meminit=False)
    try:
        with env.begin(write=False) as txn:
            keys = set(txn.cursor().iternext(values=False))
        return keys
    finally:
        env.close()


def _chunked(items, chunk_size):
    chunk = []
    for item in items:
        chunk.append(item)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


# ---------------------------------------------------------------------------
# Size estimation (single-process, runs before workers start)
# ---------------------------------------------------------------------------
def estimate_size(dataset, builder, n_preview=100):
    indices = dataset.valid_idx[:min(n_preview, len(dataset.valid_idx))]
    sizes = []
    for idx in tqdm(indices, desc="Size estimate", leave=False):
        try:
            data = dataset.getitem_with_real_idx(idx)
            if str(getattr(data, "str_tag", "complex")) != "complex":
                continue
            payload = builder(data)
            blob = pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)
            sizes.append(len(blob))
        except Exception as e:
            warnings.warn(f"[Estimate] idx={idx}: {e}")

    if sizes:
        import numpy as np
        mean_sz = np.mean(sizes)
        total_est = mean_sz * len(dataset.valid_idx)
        print(f"[Estimate] {len(sizes)} samples, mean={_fmt_bytes(mean_sz)}, "
              f"projected={_fmt_bytes(total_est)}")
    return sizes


# ---------------------------------------------------------------------------
# Disk safety
# ---------------------------------------------------------------------------
GIB = 1024 ** 3
SAFETY_MARGIN = 10 * GIB
MIN_MAP_SIZE = 1 * GIB
GROWTH_FACTOR = 1.4


def _round_up_gib(num_bytes: int) -> int:
    return max(MIN_MAP_SIZE, int(math.ceil(num_bytes / GIB) * GIB))


def _disk_free(path: str) -> int:
    return int(shutil.disk_usage(path).free)


def _estimated_map_size(num_unique_keys: int, avg_bytes: int) -> int:
    raw = int(math.ceil(num_unique_keys * avg_bytes * GROWTH_FACTOR))
    return _round_up_gib(max(raw, MIN_MAP_SIZE))


# ---------------------------------------------------------------------------
# Worker process (module-level for Windows spawn picklability)
# ---------------------------------------------------------------------------
def _worker_process(worker_id, config_path, split, k, max_substrate_length,
                    task_queue, result_queue, stop_event):
    """Worker: create own StructureDataset, process chunks, return (key, blob)."""
    import torch
    torch.set_num_threads(1)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass

    dataset = None
    try:
        config = load_config(config_path)
        df = _get_split_df(config, split)
        dataset = StructureDataset(config=config, df=df, transform=None, is_train=False)
        builder = BuildStructureCacheData(k=k, max_substrate_length=max_substrate_length)
        dock_indices = dataset.df["Dock Index"].values

        while not stop_event.is_set():
            try:
                task = task_queue.get(timeout=2.0)
            except queue.Empty:
                continue
            if task is None:  # poison pill
                break

            batch_id, real_indices = task
            records = []
            errors = []

            for real_idx in real_indices:
                if stop_event.is_set():
                    break
                try:
                    data = dataset.getitem_with_real_idx(int(real_idx))
                    if str(getattr(data, "str_tag", "complex")) != "complex":
                        continue
                    payload = builder(data)
                    blob = pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)
                    key = str(int(dock_indices[int(real_idx)])).encode()
                    records.append((key, blob))
                except Exception as e:
                    errors.append((int(real_idx), repr(e)))

            result_queue.put((batch_id, records, errors))

    except Exception as e:
        try:
            result_queue.put(("__FATAL__", worker_id, repr(e)))
        except Exception:
            pass
    finally:
        if dataset is not None:
            _close_handles(dataset)
        try:
            result_queue.cancel_join_thread()
        except Exception:
            pass


# ---------------------------------------------------------------------------
# Parallel build for one dataset
# ---------------------------------------------------------------------------
def _build_dataset_parallel(*, config_path, split, ds_name,
                            pending_real_indices, output_path,
                            resume_key_set, map_size, commit_interval,
                            num_workers, chunk_size, queue_size, k,
                            max_substrate_length):
    """Build one dataset's LMDB using worker processes + single main writer."""
    if not pending_real_indices:
        return {"writes": 0, "errors": 0, "dup": 0}

    ctx = mp.get_context("spawn")
    num_workers = max(1, min(num_workers, len(pending_real_indices)))
    task_maxsize = max(num_workers * 2, 4)

    task_queue = ctx.Queue(maxsize=task_maxsize)
    result_queue = ctx.Queue(maxsize=max(2, queue_size))
    stop_event = ctx.Event()

    workers = []
    for wid in range(num_workers):
        p = ctx.Process(
            target=_worker_process,
            args=(wid, config_path, split, k, max_substrate_length,
                  task_queue, result_queue, stop_event),
            daemon=True,
        )
        p.start()
        workers.append(p)

    env = None
    txn = None
    seen = set(resume_key_set)
    writes = 0
    errors = 0
    dup = 0

    chunks = list(_chunked(pending_real_indices, chunk_size))
    total_items = len(pending_real_indices)
    next_chunk = 0
    in_flight = 0

    try:
        env = _ORIGINAL_LMDB_OPEN(
            output_path, map_size=map_size,
            create=True, subdir=False, readonly=False,
        )
        txn = env.begin(write=True)

        # Prime the task queue
        while next_chunk < len(chunks) and in_flight < task_maxsize:
            task_queue.put((next_chunk, chunks[next_chunk]))
            next_chunk += 1
            in_flight += 1

        pbar = tqdm(total=total_items, desc=f"Building {ds_name}",
                    unit="sample", smoothing=0.1)
        no_progress_count = 0
        MAX_NO_PROGRESS = 5  # 5 × 60s = 5 min without any result → abort

        while in_flight > 0:
            try:
                result = result_queue.get(timeout=60.0)
                no_progress_count = 0  # reset on any result
            except queue.Empty:
                no_progress_count += 1
                # Check for dead workers
                dead = [p for p in workers
                        if not p.is_alive() and p.exitcode not in (0, None)]
                if dead:
                    raise RuntimeError(
                        f"Worker crashed: "
                        + ", ".join(f"pid={p.pid} exit={p.exitcode}" for p in dead))
                if no_progress_count >= MAX_NO_PROGRESS:
                    raise RuntimeError(
                        f"No progress for {MAX_NO_PROGRESS} min. "
                        f"Workers may be stuck. in_flight={in_flight}")
                continue

            if result and result[0] == "__FATAL__":
                _, wid, msg = result
                raise RuntimeError(f"Worker {wid} fatal error: {msg}")

            batch_id, records, batch_errors = result
            in_flight -= 1

            for real_idx, msg in batch_errors:
                errors += 1
                if errors <= 10:
                    warnings.warn(f"[Build:{ds_name}] idx={real_idx}: {msg}")

            batch_wrote = 0
            for key, blob in records:
                if key in seen:
                    dup += 1
                    continue
                txn.put(key, blob)
                seen.add(key)
                writes += 1
                batch_wrote += 1

                if writes % commit_interval == 0:
                    txn.commit()
                    txn = env.begin(write=True)

            pbar.update(len(records) + len(batch_errors))

            # Dispatch more work
            while next_chunk < len(chunks) and in_flight < task_maxsize:
                task_queue.put((next_chunk, chunks[next_chunk]))
                next_chunk += 1
                in_flight += 1

        pbar.close()

    finally:
        # Signal workers to stop
        stop_event.set()
        for _ in workers:
            try:
                task_queue.put_nowait(None)
            except Exception:
                pass

        if txn is not None:
            try:
                txn.commit()
            except Exception:
                pass
        if env is not None:
            try:
                env.sync()
            except Exception:
                pass
            try:
                env.close()
            except Exception:
                pass

        for p in workers:
            p.join(timeout=15)
            if p.is_alive():
                p.terminate()
                p.join(timeout=5)

        try:
            task_queue.close()
            result_queue.close()
        except Exception:
            pass

    return {"writes": writes, "errors": errors, "dup": dup}


# ---------------------------------------------------------------------------
# Main build — one dataset at a time, parallel workers per dataset
# ---------------------------------------------------------------------------
def build_cache(config, config_path, split, output_dir, commit_interval,
                overwrite, num_workers, chunk_size, queue_size):
    if commit_interval <= 0:
        raise ValueError(f"commit_interval must be > 0, got {commit_interval}")
    cache_dir = output_dir or _default_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)

    ds_names = _dataset_names(config)
    n_datasets = len(ds_names)
    k = int(config.transform.k)
    max_substrate_length = int(config.data.max_substrate_length)

    # Output paths + resume scan
    output_paths = []
    resume_keys = [set() for _ in range(n_datasets)]
    for i, name in enumerate(ds_names):
        path = os.path.join(cache_dir, f"complex_cache_{name}.lmdb")
        output_paths.append(path)
        if os.path.exists(path):
            if overwrite:
                os.remove(path)
                lock_path = path + "-lock"
                if os.path.exists(lock_path):
                    os.remove(lock_path)
            else:
                resume_keys[i] = _scan_existing_keys(path)
                if resume_keys[i]:
                    print(f"[Resume] dataset {i} ({name}): {len(resume_keys[i])} existing entries")

    print(f"[Build] Cache directory: {cache_dir}")
    for i, p in enumerate(output_paths):
        print(f"  dataset {i} ({ds_names[i]}): {os.path.basename(p)}")

    # Load data (main process only, for metadata extraction)
    print(f"[Build] Loading DataFrame (split={split})...")
    df = _get_split_df(config, split)
    print(f"[Build] DataFrame: {len(df)} rows")

    print("[Build] Initializing StructureDataset for metadata + size estimate...")
    t0 = time.time()
    dataset = StructureDataset(config=config, df=df, transform=None, is_train=False)
    n_valid = len(dataset.valid_idx)
    print(f"[Build] Ready in {time.time()-t0:.1f}s, valid={n_valid}")

    builder = BuildStructureCacheData(k=k, max_substrate_length=max_substrate_length)

    # Size estimate
    preview_sizes = estimate_size(dataset, builder)
    if not preview_sizes:
        _close_handles(dataset)
        raise RuntimeError("Size estimate produced no complex payloads.")
    avg_sample_bytes = int(sum(preview_sizes) / len(preview_sizes))
    print(f"[Estimate] Average payload: {_fmt_bytes(avg_sample_bytes)} "
          f"(from {len(preview_sizes)} preview samples)")

    # Extract metadata before closing main dataset
    valid_idx = list(dataset.valid_idx)
    dataset_ids = dataset.df["dataset_id"].values.copy()
    dock_indices = dataset.df["Dock Index"].values.copy()

    # Close main dataset to free RAM before spawning workers
    _close_handles(dataset)
    del dataset, builder, df
    print("[Build] Main dataset closed (RAM freed for workers)")

    # Organize valid rows by dataset_id
    rows_by_ds = [[] for _ in range(n_datasets)]
    target_keys_by_ds = [set() for _ in range(n_datasets)]
    for real_idx in valid_idx:
        ds_id = int(dataset_ids[real_idx])
        dock_idx = int(dock_indices[real_idx])
        rows_by_ds[ds_id].append(real_idx)
        target_keys_by_ds[ds_id].add(str(dock_idx).encode())

    # Print plan
    print(f"\n[Plan] Sequential per-dataset, {num_workers} workers per dataset")
    for ds_id, name in enumerate(ds_names):
        n_unique = len(target_keys_by_ds[ds_id])
        est_ms = _estimated_map_size(n_unique, avg_sample_bytes)
        print(f"  [{ds_id}] {name}: rows={len(rows_by_ds[ds_id])}, "
              f"unique_docks={n_unique}, resume={len(resume_keys[ds_id])}, "
              f"map_size={_fmt_bytes(est_ms)}")

    total_new = 0
    total_err = 0
    total_dup = 0
    total_skipped_resume = 0
    overall_start = time.time()

    for ds_id, name in enumerate(ds_names):
        ds_rows = rows_by_ds[ds_id]
        target_keys = target_keys_by_ds[ds_id]
        output_path = output_paths[ds_id]

        if not ds_rows:
            print(f"\n[Build] Dataset {ds_id} ({name}): no valid rows, skipping.")
            continue

        # Pre-deduplicate: skip resume keys + intra-dataset duplicates
        pending = []
        scheduled_keys = set()
        ds_skipped_resume = 0
        ds_pre_dup = 0
        for real_idx in ds_rows:
            key = str(int(dock_indices[real_idx])).encode()
            if key in resume_keys[ds_id]:
                ds_skipped_resume += 1
                continue
            if key in scheduled_keys:
                ds_pre_dup += 1
                continue
            scheduled_keys.add(key)
            pending.append(real_idx)

        # Compute map_size
        existing_size = os.path.getsize(output_path) if os.path.exists(output_path) else 0
        map_size = max(_estimated_map_size(len(target_keys), avg_sample_bytes), existing_size)
        map_size = _round_up_gib(map_size)

        # Disk safety check
        free_before = _disk_free(cache_dir)
        additional_needed = max(0, map_size - existing_size)
        if free_before < additional_needed + SAFETY_MARGIN:
            raise RuntimeError(
                f"[Safety] Not enough disk for dataset {ds_id} ({name}). "
                f"free={_fmt_bytes(free_before)}, need={_fmt_bytes(additional_needed)}, "
                f"safety_margin={_fmt_bytes(SAFETY_MARGIN)}, map_size={_fmt_bytes(map_size)}")

        print(f"\n[Build] Dataset {ds_id} ({name})")
        print(f"  unique_docks={len(target_keys)}, pre_existing={len(resume_keys[ds_id])}, "
              f"pending_new={len(pending)}, pre_dup={ds_pre_dup}")
        print(f"  map_size={_fmt_bytes(map_size)}, existing={_fmt_bytes(existing_size)}, "
              f"free={_fmt_bytes(free_before)}")

        ds_start = time.time()

        stats = _build_dataset_parallel(
            config_path=config_path,
            split=split,
            ds_name=name,
            pending_real_indices=pending,
            output_path=output_path,
            resume_key_set=resume_keys[ds_id],
            map_size=map_size,
            commit_interval=commit_interval,
            num_workers=num_workers,
            chunk_size=chunk_size,
            queue_size=queue_size,
            k=k,
            max_substrate_length=max_substrate_length,
        )

        # Per-dataset verification
        final_keys = _scan_existing_keys(output_path)
        expected = len(target_keys)
        actual = len(final_keys)
        status = "OK" if actual == expected else "MISMATCH"
        ds_elapsed = time.time() - ds_start
        file_size = os.path.getsize(output_path) if os.path.exists(output_path) else 0

        print(f"  finished in {ds_elapsed/60:.1f} min")
        print(f"  new_writes={stats['writes']}, errors={stats['errors']}, "
              f"dup={stats['dup'] + ds_pre_dup}, skipped_resume={ds_skipped_resume}")
        print(f"  verify: expected={expected}, actual={actual} [{status}]")
        print(f"  file_size={_fmt_bytes(file_size)}")

        if actual != expected:
            raise RuntimeError(
                f"[Verify] Dataset {ds_id} ({name}) count mismatch: "
                f"expected={expected}, actual={actual}")

        total_new += stats["writes"]
        total_err += stats["errors"]
        total_dup += stats["dup"] + ds_pre_dup
        total_skipped_resume += ds_skipped_resume

    elapsed = time.time() - overall_start
    total_pre = sum(len(k) for k in resume_keys)
    print(f"\n[Build] Complete in {elapsed/60:.1f} min")
    print(f"  new_writes={total_new}, errors={total_err}, duplicates={total_dup}, "
          f"skipped_resume={total_skipped_resume}, pre_existing={total_pre}")

    print("\n[Verify] Final LMDB entry counts:")
    for ds_id, name in enumerate(ds_names):
        path = output_paths[ds_id]
        if not os.path.exists(path):
            print(f"  {name}: missing")
            continue
        actual = len(_scan_existing_keys(path))
        expected = len(target_keys_by_ds[ds_id])
        status = "OK" if actual == expected else "MISMATCH"
        size = os.path.getsize(path)
        print(f"  {name}: expected={expected}, actual={actual} [{status}], size={_fmt_bytes(size)}")

    if total_err > 0:
        print(f"\n[WARNING] {total_err} errors occurred during build.")


def main():
    args = parse_args()
    config_path = args.config
    print(f"[Config] {config_path}")
    config = load_config(config_path)
    print(f"[Config] split={args.split}, commit={args.commit_interval}, "
          f"overwrite={args.overwrite}")
    print(f"[MP] num_workers={args.num_workers}, chunk_size={args.chunk_size}, "
          f"queue_size={args.queue_size}")
    print("[Safety] Sequential datasets + bounded queues + thread caps")

    build_cache(
        config=config,
        config_path=config_path,
        split=args.split,
        output_dir=args.output_dir,
        commit_interval=args.commit_interval,
        overwrite=args.overwrite,
        num_workers=args.num_workers,
        chunk_size=args.chunk_size,
        queue_size=args.queue_size,
    )


if __name__ == "__main__":
    main()
