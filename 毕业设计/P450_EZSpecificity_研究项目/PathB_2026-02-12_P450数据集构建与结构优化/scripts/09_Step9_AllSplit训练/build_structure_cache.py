"""
Offline structure cache builder for EZSpecificity Step 9.

Iterates valid samples from the original StructureDataset (transform=None),
runs BuildStructureCacheData to extract split real/knn partitions, and writes
them to LMDB. One LMDB per dataset_id.

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe "../毕业设计/.../scripts/09_Step9_AllSplit训练/build_structure_cache.py" --split all
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import time
import warnings

import lmdb
import yaml
from easydict import EasyDict
from tqdm.auto import tqdm

# ---------------------------------------------------------------------------
# 1. Monkey-patch lmdb.open for large databases
# ---------------------------------------------------------------------------
_ORIGINAL_LMDB_OPEN = lmdb.open
_READ_MAP_SIZE = 256 * 1024 ** 3
_WRITE_MAP_SIZE = 200 * 1024 ** 3  # 512GB fails on Windows; 200GB is safe


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
    # Use ASCII-only path to avoid LMDB Unicode issues on Windows
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


# ---------------------------------------------------------------------------
# Size estimation
# ---------------------------------------------------------------------------
def estimate_size(dataset, builder, n_preview=100):
    indices = dataset.valid_idx[:min(n_preview, len(dataset.valid_idx))]
    sizes = []
    preview = {}
    for idx in tqdm(indices, desc="Size estimate", leave=False):
        try:
            data = dataset.getitem_with_real_idx(idx)
            if str(getattr(data, "str_tag", "complex")) != "complex":
                continue
            payload = builder(data)
            blob = pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)
            sizes.append(len(blob))
            ds_id = int(dataset.df.loc[idx, "dataset_id"])
            preview[idx] = (ds_id, blob)
        except Exception as e:
            warnings.warn(f"[Estimate] idx={idx}: {e}")

    if sizes:
        import numpy as np
        mean_sz = np.mean(sizes)
        total_est = mean_sz * len(dataset.valid_idx)
        print(f"[Estimate] {len(sizes)} samples, mean={_fmt_bytes(mean_sz)}, "
              f"projected={_fmt_bytes(total_est)}")
    return preview


# ---------------------------------------------------------------------------
# Main build
# ---------------------------------------------------------------------------
def _scan_existing_keys(path: str) -> set[bytes]:
    """Scan an existing LMDB for all keys (for resume support)."""
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


def build_cache(config, split, output_dir, commit_interval, overwrite):
    cache_dir = output_dir or _default_cache_dir()
    os.makedirs(cache_dir, exist_ok=True)

    ds_names = _dataset_names(config)
    n_datasets = len(ds_names)

    # Output paths: one LMDB per dataset_id
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
                # Resume mode: scan existing keys
                resume_keys[i] = _scan_existing_keys(path)
                if resume_keys[i]:
                    print(f"[Resume] dataset {i} ({name}): {len(resume_keys[i])} existing entries")

    print(f"[Build] Cache directory: {cache_dir}")
    for i, p in enumerate(output_paths):
        print(f"  dataset {i} ({ds_names[i]}): {os.path.basename(p)}")

    # Load data
    print(f"[Build] Loading DataFrame (split={split})...")
    df = _get_split_df(config, split)
    print(f"[Build] DataFrame: {len(df)} rows")

    # StructureDataset with NO transform — gives raw data with str_tag/sample_weight
    print("[Build] Initializing StructureDataset (transform=None)...")
    t0 = time.time()
    dataset = StructureDataset(config=config, df=df, transform=None, is_train=False)
    n_valid = len(dataset.valid_idx)
    print(f"[Build] Ready in {time.time()-t0:.1f}s, valid={n_valid}")

    builder = BuildStructureCacheData(
        k=int(config.transform.k),
        max_substrate_length=int(config.data.max_substrate_length),
    )

    # Size estimate
    preview = estimate_size(dataset, builder)

    # Open output LMDBs (create=True works for both new and existing)
    envs = []
    for path in output_paths:
        envs.append(_ORIGINAL_LMDB_OPEN(
            path, map_size=_WRITE_MAP_SIZE, create=True, subdir=False, readonly=False,
        ))
    txns = [env.begin(write=True) for env in envs]
    writes = [0] * n_datasets
    # Initialize seen sets from resume keys
    seen = [set(rk) for rk in resume_keys]
    n_resumed = sum(len(s) for s in seen)
    if n_resumed > 0:
        print(f"[Resume] Skipping {n_resumed} already-cached entries across all datasets")

    ok, err, dup, skip, skipped_resume = 0, 0, 0, 0, 0
    t_start = time.time()

    try:
        for real_idx in tqdm(dataset.valid_idx, desc=f"Building cache ({split})"):
            try:
                ds_id = int(dataset.df.loc[real_idx, "dataset_id"])
                dock_idx = int(dataset.df.loc[real_idx, "Dock Index"])
                key = str(dock_idx).encode()

                if key in seen[ds_id]:
                    # Distinguish resumed vs true duplicate
                    if key in resume_keys[ds_id]:
                        skipped_resume += 1
                    else:
                        dup += 1
                    continue

                # Reuse preview if available
                if real_idx in preview:
                    blob = preview[real_idx][1]
                else:
                    data = dataset.getitem_with_real_idx(real_idx)
                    if str(getattr(data, "str_tag", "complex")) != "complex":
                        skip += 1
                        continue
                    payload = builder(data)
                    blob = pickle.dumps(payload, protocol=pickle.HIGHEST_PROTOCOL)

                txns[ds_id].put(key, blob)
                seen[ds_id].add(key)
                writes[ds_id] += 1
                ok += 1

                if writes[ds_id] % commit_interval == 0:
                    txns[ds_id].commit()
                    txns[ds_id] = envs[ds_id].begin(write=True)

            except KeyboardInterrupt:
                raise
            except Exception as e:
                err += 1
                if err <= 10:
                    warnings.warn(f"[Build] idx={real_idx}: {e}")
    finally:
        for i in range(n_datasets):
            try:
                txns[i].commit()
            except Exception:
                pass
            try:
                envs[i].sync()
                envs[i].close()
            except Exception:
                pass
        _close_handles(dataset)

    elapsed = time.time() - t_start
    print(f"\n[Build] Complete in {elapsed/60:.1f} min")
    print(f"  new_writes={ok}, errors={err}, duplicates={dup}, non_complex={skip}, "
          f"skipped_resume={skipped_resume}, pre_existing={n_resumed}")
    for i, name in enumerate(ds_names):
        sz = os.path.getsize(output_paths[i]) if os.path.exists(output_paths[i]) else 0
        print(f"  {name}: {writes[i]} new + {len(resume_keys[i])} resumed, {_fmt_bytes(sz)}")

    # Post-build count verification
    print("\n[Verify] Checking final LMDB entry counts...")
    all_ok = True
    for i, name in enumerate(ds_names):
        if not os.path.exists(output_paths[i]):
            continue
        final_keys = _scan_existing_keys(output_paths[i])
        expected = writes[i] + len(resume_keys[i])
        actual = len(final_keys)
        status = "OK" if actual == expected else "MISMATCH"
        if actual != expected:
            all_ok = False
        print(f"  {name}: expected={expected}, actual={actual} [{status}]")
    if all_ok:
        print("[Verify] All counts match.")
    else:
        print("[Verify] WARNING: Count mismatches detected!")

    if err > 0:
        print(f"\n[WARNING] {err} errors occurred during build. Inspect and consider rebuilding affected samples.")


def main():
    args = parse_args()
    print(f"[Config] {args.config}")
    config = load_config(args.config)
    print(f"[Config] split={args.split}, commit={args.commit_interval}, overwrite={args.overwrite}")

    build_cache(
        config=config,
        split=args.split,
        output_dir=args.output_dir,
        commit_interval=args.commit_interval,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
