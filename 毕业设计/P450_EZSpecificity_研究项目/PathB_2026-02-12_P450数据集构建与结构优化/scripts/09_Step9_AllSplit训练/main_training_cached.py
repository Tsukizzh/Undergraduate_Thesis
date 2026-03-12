"""
Step 9: Train EZSpecificity from scratch with CACHED structure features.

Drop-in replacement for main_training.py. Uses:
- CachedStructureSequenceDataset (sequence live + structure cached)
- RebuildComplexEdgeAttr (runtime edge reconstruction, supports fixed/legacy_bug mode)

Fixes the edge-attr/edge-index ordering bug in original EdgeConnection (transforms.py:130-147)
when run with --edge-mode fixed (default).

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe "../毕业设计/.../scripts/09_Step9_AllSplit训练/main_training_cached.py"
    D:/anaconda3/envs/torch/python.exe "../毕业设计/.../scripts/09_Step9_AllSplit训练/main_training_cached.py" --edge-mode legacy_bug
"""
from __future__ import annotations

import argparse
import os
import random as _random
import sys
import time
import warnings
from collections import defaultdict
from queue import Queue
from threading import Thread

import lmdb
import yaml

warnings.filterwarnings("ignore")

# ============================================================
# 1. Monkey-patch lmdb.open for large databases
# ============================================================
_original_lmdb_open = lmdb.open


def _patched_lmdb_open(path, **kwargs):
    if kwargs.get("readonly", False) or not kwargs.get("create", True):
        kwargs["map_size"] = 256 * 1024 ** 3
    return _original_lmdb_open(path, **kwargs)


lmdb.open = _patched_lmdb_open

# ============================================================
# 2. Paths
# ============================================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", "src"))
PATHB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))

for d in (SCRIPT_DIR, SRC_DIR):
    if d not in sys.path:
        sys.path.insert(0, d)

CONFIG_PATH = os.path.join(SCRIPT_DIR, "train_allsplit_config.yml")
LOG_DIR = os.path.join(PATHB_DIR, "logs", "09_Step9_AllSplit训练")
CKPT_DIR = os.path.join(PATHB_DIR, "results", "09_Step9_AllSplit训练", "checkpoints")
CACHE_DIR = os.path.normpath(os.path.join(SRC_DIR, "..", "data", "step9_structure_cache"))

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(CKPT_DIR, exist_ok=True)

# ============================================================
# 3. Imports
# ============================================================
import torch
import pytorch_lightning as pl
from easydict import EasyDict
from pytorch_lightning.callbacks import EarlyStopping, LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger
from torch_geometric.data import DataLoader

from cache_utils import CachedStructureSequenceDataset, RebuildComplexEdgeAttr
from Datasets.utils import read_datasets
from Models.ss import SS

# TF32 + cuDNN
torch.set_float32_matmul_precision("high")
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
torch.backends.cudnn.benchmark = True

# Non-blocking H2D transfer
_orig_transfer = SS.transfer_batch_to_device


def _fast_transfer(self, batch, device, dataloader_idx):
    return batch.to(device, non_blocking=True)


SS.transfer_batch_to_device = _fast_transfer

NUM_WORKERS = 1


# ============================================================
# 4. CLI
# ============================================================
def parse_args():
    p = argparse.ArgumentParser(description="Train EZSpecificity with cached structure data")
    p.add_argument("--config", default=CONFIG_PATH)
    p.add_argument("--edge-mode", choices=["fixed", "legacy_bug"], default="fixed",
                   help="Edge ordering mode (fixed=bug-corrected, legacy_bug=original)")
    p.add_argument("--cache-dir", default=CACHE_DIR)
    p.add_argument("--num-workers", type=int, default=NUM_WORKERS)
    p.add_argument("--benchmark", action="store_true",
                   help="Benchmark mode: disable logging/ckpt/val, run limited batches")
    p.add_argument("--benchmark-batches", type=int, default=500,
                   help="Total train batches in benchmark mode (includes warmup)")
    p.add_argument("--benchmark-warmup", type=int, default=50,
                   help="Warmup batches to skip before measuring speed")
    p.add_argument("--no-prefetch-wrapper", action="store_true",
                   help="Disable BackgroundPrefetchLoader (for A/B testing)")
    p.add_argument("--resume", type=str, default=None,
                   help="Resume from checkpoint path (or 'last' to auto-find last.ckpt)")
    return p.parse_args()


def load_config(path: str) -> EasyDict:
    with open(path, "r", encoding="utf-8") as f:
        return EasyDict(yaml.safe_load(f))


# ============================================================
# 5. Cache path resolution
# ============================================================
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


def resolve_cache_paths(config, cache_dir: str) -> list[str]:
    names = _dataset_names(config)
    paths = [os.path.join(cache_dir, f"complex_cache_{n}.lmdb") for n in names]
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(
            f"Cache LMDB(s) not found:\n  " + "\n  ".join(missing) +
            "\nRun build_structure_cache.py first."
        )
    return paths


# ============================================================
# 6. BackgroundPrefetchLoader
# ============================================================
class BackgroundPrefetchLoader:
    def __init__(self, loader, max_prefetch=4):
        self.loader = loader
        self.max_prefetch = max_prefetch

    def __len__(self):
        return len(self.loader)

    def __getattr__(self, name):
        return getattr(self.loader, name)

    def __iter__(self):
        q = Queue(maxsize=self.max_prefetch)
        sentinel = object()
        errors = []

        def producer():
            try:
                for batch in self.loader:
                    q.put(batch)
            except Exception as e:
                errors.append(e)
            finally:
                q.put(sentinel)

        Thread(target=producer, daemon=True).start()
        while True:
            item = q.get()
            if item is sentinel:
                if errors:
                    raise errors[0]
                break
            yield item


# ============================================================
# 6b. BlockShuffleSampler — locality-aware sampling for large LMDB caches
# ============================================================
class BlockShuffleSampler(torch.utils.data.Sampler):
    """Locality-aware sampler: groups samples by LMDB file and key order,
    then shuffles blocks of contiguous samples to preserve SSD read locality
    while maintaining training randomness.
    """

    def __init__(self, dataset, block_size=256, batch_size=8, seed=0):
        self.n = len(dataset)
        self.block_size = block_size
        self.batch_size = batch_size
        self.seed = seed
        self._epoch = 0
        self._blocks = self._build_blocks(dataset)

    def _build_blocks(self, dataset):
        df = dataset.df
        valid_idx = dataset.valid_idx
        ds_col = df["dataset_id"].values
        dock_col = df["Dock Index"].values

        per_dataset = defaultdict(list)
        for item_idx, global_idx in enumerate(valid_idx):
            ds_id = int(ds_col[global_idx])
            key_bytes = str(int(dock_col[global_idx])).encode("ascii")
            per_dataset[ds_id].append((key_bytes, item_idx))

        blocks = []
        for ds_id in sorted(per_dataset):
            items = per_dataset[ds_id]
            items.sort(key=lambda pair: pair[0])  # LMDB byte-key order
            ordered = [item_idx for _, item_idx in items]
            for start in range(0, len(ordered), self.block_size):
                blocks.append(ordered[start:start + self.block_size])
        return blocks

    def __iter__(self):
        rng = _random.Random(self.seed + self._epoch)
        self._epoch += 1
        blocks = list(self._blocks)
        rng.shuffle(blocks)
        for block in blocks:
            if len(block) > self.batch_size:
                chunks = [block[i:i + self.batch_size]
                          for i in range(0, len(block), self.batch_size)]
                rng.shuffle(chunks)
                for chunk in chunks:
                    yield from chunk
            else:
                yield from block

    def __len__(self):
        return self.n


# ============================================================
# 7. DataModule
# ============================================================
class CachedTrainingDataModule(pl.LightningDataModule):
    def __init__(self, config, cache_paths, edge_mode, num_workers,
                 use_prefetch_wrapper=True):
        super().__init__()
        self.config = config
        self.cache_paths = cache_paths
        self.edge_mode = edge_mode
        self.num_workers = num_workers
        self.batch_size = config.data.batch_size
        self.use_prefetch_wrapper = use_prefetch_wrapper

        print("[DataModule] Loading CSVs...")
        t0 = time.time()
        self.train_df = read_datasets(config.data.train_data_path)
        self.val_df = read_datasets(config.data.val_data_path)
        self.test_df = read_datasets(config.data.test_data_path)
        print(f"[DataModule] CSVs loaded in {time.time()-t0:.1f}s")
        print(f"[DataModule] Train: {len(self.train_df)}, Val: {len(self.val_df)}, Test: {len(self.test_df)}")

    def _transform(self, is_train: bool) -> RebuildComplexEdgeAttr:
        return RebuildComplexEdgeAttr(
            mode=self.edge_mode,
            is_train=is_train,
            dist_noise=self.config.transform.dist_noise if is_train else False,
            cutoff=self.config.transform.cutoff,
            num_r_gaussian=self.config.transform.num_r_gaussian,
        )

    def _make_loader(self, df, is_train, shuffle=False, use_block_sampler=False):
        dataset = CachedStructureSequenceDataset(
            df=df, config=self.config,
            complex_cache_paths=self.cache_paths,
            structure_transform=self._transform(is_train),
            is_train=is_train,
        )
        sampler = None
        if use_block_sampler:
            sampler = BlockShuffleSampler(
                dataset, block_size=256, batch_size=self.batch_size,
                seed=int(getattr(self.config, "seed", 0)),
            )
            shuffle = False
            print(f"[DataModule] BlockShuffleSampler: {len(sampler._blocks)} blocks of <={sampler.block_size}")
        kwargs = dict(
            batch_size=self.batch_size, shuffle=shuffle, sampler=sampler,
            num_workers=self.num_workers, pin_memory=True,
            follow_batch=["ligand_index"],
        )
        if self.num_workers > 0:
            kwargs.update(persistent_workers=True, prefetch_factor=4)
        return dataset, DataLoader(dataset, **kwargs)

    def train_dataloader(self):
        print("[DataModule] Building train dataset...")
        t0 = time.time()
        data, loader = self._make_loader(self.train_df, is_train=True, shuffle=True)
        print(f"[DataModule] Train built in {time.time()-t0:.1f}s, valid={len(data)}")
        if self.use_prefetch_wrapper:
            return BackgroundPrefetchLoader(loader)
        return loader

    def val_dataloader(self):
        print("[DataModule] Building val dataset...")
        t0 = time.time()
        data, loader = self._make_loader(self.val_df, is_train=False)
        print(f"[DataModule] Val built in {time.time()-t0:.1f}s, valid={len(data)}")
        return loader

    def test_dataloader(self):
        print("[DataModule] Building test dataset...")
        t0 = time.time()
        data, loader = self._make_loader(self.test_df, is_train=False)
        print(f"[DataModule] Test built in {time.time()-t0:.1f}s, valid={len(data)}")
        return loader


# ============================================================
# 8. Main
# ============================================================
def main():
    args = parse_args()
    config = load_config(args.config)
    cache_paths = resolve_cache_paths(config, args.cache_dir)

    is_benchmark = args.benchmark
    use_prefetch = not args.no_prefetch_wrapper

    print("=" * 70)
    if is_benchmark:
        print("Step 9: BENCHMARK MODE (cached structure, all_split fold 0)")
        print(f"  Batches: {args.benchmark_batches} (warmup: {args.benchmark_warmup})")
        print(f"  Prefetch wrapper: {'ON' if use_prefetch else 'OFF'}")
    else:
        print("Step 9: Training EZSpecificity (cached structure, all_split fold 0)")
    print(f"  Edge mode: {args.edge_mode}")
    print("=" * 70)
    print(f"Config: {args.config}")
    print(f"GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
    print(f"CUDA: {torch.version.cuda}, PyTorch: {torch.__version__}")
    print(f"Batch: {config.data.batch_size} x {config.training.accumulate_grad_batches} "
          f"= {config.data.batch_size * config.training.accumulate_grad_batches}")
    print(f"Workers: {args.num_workers}")
    print("Cache LMDBs:")
    for i, p in enumerate(cache_paths):
        print(f"  [{i}] {p}")
    print()

    pl.seed_everything(config.seed)

    # DataModule
    print("[1/4] Initializing DataModule...")
    dm = CachedTrainingDataModule(
        config, cache_paths, args.edge_mode, args.num_workers,
        use_prefetch_wrapper=use_prefetch,
    )

    # Model
    print("[2/4] Initializing Model (random weights)...")
    model = SS(config)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[2/4] Trainable parameters: {n_params:,}")

    if is_benchmark:
        # Validate benchmark args
        assert args.benchmark_batches > 0, "--benchmark-batches must be > 0"
        assert args.benchmark_warmup >= 0, "--benchmark-warmup must be >= 0"
        assert args.benchmark_batches > args.benchmark_warmup, \
            f"--benchmark-batches ({args.benchmark_batches}) must be > --benchmark-warmup ({args.benchmark_warmup})"

        # Benchmark mode: no logging, no checkpointing, no validation
        # Disable ReduceLROnPlateau (needs val metric which is unavailable)
        _orig_configure_optimizers = model.configure_optimizers

        def _bench_configure_optimizers():
            result = _orig_configure_optimizers()
            if isinstance(result, tuple) and len(result) == 2:
                return result[0]  # return only optimizers, drop schedulers
            return result

        model.configure_optimizers = _bench_configure_optimizers

        print("[3/4] Initializing Trainer (BENCHMARK)...")
        trainer = pl.Trainer(
            benchmark=True, max_epochs=1,
            accelerator="gpu", devices=1, precision=16,
            gradient_clip_val=config.training.gradient_clip_val,
            accumulate_grad_batches=config.training.accumulate_grad_batches,
            limit_train_batches=args.benchmark_batches,
            limit_val_batches=0,
            num_sanity_val_steps=0,
            callbacks=[],
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=True, log_every_n_steps=50,
        )

        warmup_n = args.benchmark_warmup

        class BenchmarkCallback(pl.Callback):
            def __init__(self):
                self._measure_start = None
                self._measure_end = None
                self._n_measured = 0

            def on_train_batch_start(self, trainer, pl_module, batch, batch_idx):
                if batch_idx == warmup_n and self._measure_start is None:
                    torch.cuda.synchronize()
                    self._measure_start = time.perf_counter()

            def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
                if batch_idx == warmup_n - 1 and warmup_n > 0:
                    # End of warmup: sync GPU and start timer
                    torch.cuda.synchronize()
                    self._measure_start = time.perf_counter()
                elif batch_idx >= warmup_n:
                    self._n_measured += 1

            def on_train_end(self, trainer, pl_module):
                torch.cuda.synchronize()
                self._measure_end = time.perf_counter()

        bench_cb = BenchmarkCallback()
        trainer.callbacks.append(bench_cb)

        print(f"[4/4] Starting benchmark ({args.benchmark_batches} batches, "
              f"{warmup_n} warmup + {args.benchmark_batches - warmup_n} measured)...")
        print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")

        trainer.fit(model, datamodule=dm)

        # Compute results
        if bench_cb._measure_start is not None and bench_cb._n_measured > 0:
            elapsed_measure = bench_cb._measure_end - bench_cb._measure_start
            n_measured = bench_cb._n_measured
            its = n_measured / elapsed_measure
            print()
            print("=" * 70)
            print("BENCHMARK RESULTS")
            print(f"  Warmup batches: {warmup_n}")
            print(f"  Measured batches: {n_measured}")
            print(f"  Measured time: {elapsed_measure:.1f}s")
            print(f"  Speed: {its:.2f} it/s ({1/its:.3f} s/batch)")
            print(f"  Workers: {args.num_workers}, Prefetch wrapper: {'ON' if use_prefetch else 'OFF'}")
            print(f"  Method 1 baseline: 4.3 it/s")
            print(f"  Speedup: {its/4.3:.2f}x")
            print("=" * 70)
        else:
            print("[WARN] Not enough measured batches for timing.")

    else:
        # Full training mode
        lr_monitor = LearningRateMonitor(logging_interval="step")
        ckpt_cb = ModelCheckpoint(
            dirpath=CKPT_DIR,
            filename=f"allsplit-fold0-{args.edge_mode}-epoch{{epoch:02d}}-aupr{{aupr/val:.4f}}",
            monitor="aupr/val", mode="max",
            save_top_k=3, save_last=True, verbose=True,
            auto_insert_metric_name=False,
        )
        early_stop_cb = EarlyStopping(
            monitor="aupr/val", mode="max", patience=15, verbose=True,
        )
        logger = TensorBoardLogger(save_dir=LOG_DIR, name=f"allsplit_fold0_{args.edge_mode}")

        print("[3/4] Initializing Trainer...")
        trainer = pl.Trainer(
            benchmark=True, max_epochs=50,
            accelerator="gpu", devices=1, precision=16,
            gradient_clip_val=config.training.gradient_clip_val,
            accumulate_grad_batches=config.training.accumulate_grad_batches,
            check_val_every_n_epoch=config.training.val_frequency,
            num_sanity_val_steps=2,
            callbacks=[lr_monitor, ckpt_cb, early_stop_cb],
            logger=logger,
            enable_progress_bar=True, log_every_n_steps=50,
        )

        # Resume support
        ckpt_path = None
        if args.resume:
            if args.resume.lower() == "last":
                last_ckpt = os.path.join(CKPT_DIR, "last.ckpt")
                if os.path.exists(last_ckpt):
                    ckpt_path = last_ckpt
                    print(f"[4/4] Resuming from: {ckpt_path}")
                else:
                    print(f"[WARN] last.ckpt not found at {last_ckpt}, training from scratch")
            else:
                if os.path.exists(args.resume):
                    ckpt_path = args.resume
                    print(f"[4/4] Resuming from: {ckpt_path}")
                else:
                    raise FileNotFoundError(f"Checkpoint not found: {args.resume}")

        print("[4/4] Starting training...")
        print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        t0 = time.time()

        trainer.fit(model, datamodule=dm, ckpt_path=ckpt_path)

        elapsed = time.time() - t0
        print()
        print("=" * 70)
        print("Training complete!")
        print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Time: {elapsed:.0f}s = {elapsed/60:.1f}min = {elapsed/3600:.2f}h")
        n_batches = trainer.num_training_batches
        n_epochs = trainer.current_epoch  # PL 1.9: already = completed epochs count
        if n_batches > 0 and n_epochs > 0:
            total_batches = n_batches * n_epochs
            print(f"Batches/epoch: {n_batches}, Epochs: {n_epochs}, "
                  f"Total batches: {total_batches}, avg {elapsed/total_batches:.3f}s/batch")
        print(f"Best ckpt: {ckpt_cb.best_model_path}")
        print(f"Best aupr/val: {ckpt_cb.best_model_score}")
        print("=" * 70)


if __name__ == "__main__":
    main()
