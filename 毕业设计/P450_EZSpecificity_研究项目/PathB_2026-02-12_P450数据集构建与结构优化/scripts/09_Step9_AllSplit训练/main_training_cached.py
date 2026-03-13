"""
Step 9: Train EZSpecificity from scratch with CACHED structure features.

Drop-in replacement for main_training.py. Uses:
- CachedStructureSequenceDataset (sequence live + structure cached)
- RebuildComplexEdgeAttr (runtime edge reconstruction, supports fixed/legacy_bug mode)

Fixes the edge-attr/edge-index ordering bug in original EdgeConnection (transforms.py:130-147)
when run with --edge-mode fixed (default).

Checkpoint selection & early stopping: monitor AUC-ROC (auc/val).
Metrics CSV: epoch, lr, grad_norm, train_loss, val_loss, auc, aupr, per-family (7x2).
Plots: auto-generated after training (3x2 training dynamics + per-family bar chart).

=== HOW TO START TRAINING ===

  # Fresh start:
  cd D:/EZSpecificity_Project/src
  D:/anaconda3/envs/torch/python.exe "../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/main_training_cached.py" --edge-mode fixed --num-workers 2 --no-prefetch-wrapper

  # Resume from last checkpoint:
  cd D:/EZSpecificity_Project/src
  D:/anaconda3/envs/torch/python.exe "../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/main_training_cached.py" --edge-mode fixed --num-workers 2 --no-prefetch-wrapper --resume last

  # Stop training:
  taskkill /IM python.exe /F

=== HOW TO EVALUATE ===

  cd D:/EZSpecificity_Project/src
  D:/anaconda3/envs/torch/python.exe "../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/eval_checkpoints.py"

=== OUTPUT FILES ===

  results/09_Step9_AllSplit训练/
  ├── checkpoints/           # Model checkpoints (top-3 by AUC + last.ckpt)
  ├── metrics_history.csv    # Per-epoch metrics (written during training, append-mode)
  └── eval_epochXX/          # Evaluation figures (created by eval_checkpoints.py)
      ├── fig1_training_dynamics.png
      ├── fig2_best_checkpoint_analysis.png
      ├── fig3_family_heatmap.png
      └── metrics_history.csv

  Data persistence: metrics_history.csv is flushed after EVERY validation epoch.
  Killing training mid-epoch only loses the current epoch's in-memory accumulators.
  All completed epochs are safe on disk.
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
from pytorch_lightning.callbacks import Callback, EarlyStopping, LearningRateMonitor, ModelCheckpoint
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
# 8. Metrics CSV Logger Callback
# ============================================================
METRICS_DIR = os.path.join(PATHB_DIR, "results", "09_Step9_AllSplit训练")
METRICS_CSV = os.path.join(METRICS_DIR, "metrics_history.csv")

# Family name mapping (tag index → family name, from config data sources)
FAMILY_NAMES = {
    "0": "brenda", "1": "Duf", "2": "Esterase", "3": "Gt_acceptor",
    "4": "Nitrilase", "5": "Phosphatase", "6": "Thiolase",
}


def _metric_val(m, key):
    """Extract float from callback_metrics, return None if missing/nan."""
    v = m.get(key)
    if v is None:
        return None
    v = v.item() if hasattr(v, "item") else float(v)
    import math
    return v if not math.isnan(v) else None


class MetricsCSVLogger(Callback):
    """Log comprehensive metrics to CSV after each validation epoch."""

    HEADER = (
        "epoch,timestamp,lr,grad_norm,"
        "train_loss,val_loss,"
        "val_auc,val_aupr,"
        "val_auc_0,val_auc_1,val_auc_2,val_auc_3,val_auc_4,val_auc_5,val_auc_6,"
        "val_aupr_0,val_aupr_1,val_aupr_2,val_aupr_3,val_aupr_4,val_aupr_5,val_aupr_6"
    )

    def __init__(self, csv_path: str = METRICS_CSV):
        super().__init__()
        self.csv_path = csv_path
        self._train_loss_sum = 0.0
        self._train_loss_count = 0
        self._grad_norm_sum = 0.0
        self._grad_norm_count = 0
        self._ensure_csv_header(csv_path)

    @staticmethod
    def _ensure_csv_header(csv_path):
        """Create CSV if missing, or migrate old CSV if header changed."""
        expected_cols = MetricsCSVLogger.HEADER.replace("\n", "").split(",")
        if not os.path.exists(csv_path):
            with open(csv_path, "w", encoding="utf-8") as f:
                f.write(MetricsCSVLogger.HEADER + "\n")
            return
        with open(csv_path, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()
        existing_cols = first_line.split(",")
        if existing_cols == expected_cols:
            return
        # Migrate: re-read all rows, re-write with new header
        import csv as csv_mod
        rows = []
        with open(csv_path, "r", encoding="utf-8") as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                rows.append(row)
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            writer = csv_mod.DictWriter(f, fieldnames=expected_cols, extrasaction="ignore")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        print(f"[MetricsCSV] Migrated {len(rows)} rows to new schema ({csv_path})")

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
        # SS.training_step returns a bare tensor, so handle all cases
        import torch
        if isinstance(outputs, torch.Tensor):
            loss = outputs
        elif isinstance(outputs, dict):
            loss = outputs.get("loss")
        else:
            loss = getattr(outputs, "loss", None)
        if loss is not None:
            self._train_loss_sum += loss.item()
            self._train_loss_count += 1

    def on_before_optimizer_step(self, trainer, pl_module, optimizer, opt_idx=0):
        import torch
        grads = [p.grad.detach().flatten() for p in pl_module.parameters() if p.grad is not None]
        if grads:
            total_norm = torch.cat(grads).norm(2).item()  # single GPU sync
            self._grad_norm_sum += total_norm
            self._grad_norm_count += 1

    def on_validation_epoch_end(self, trainer, pl_module):
        if trainer.sanity_checking:
            return
        m = trainer.callback_metrics
        epoch = trainer.current_epoch
        ts = time.strftime("%Y-%m-%d %H:%M:%S")

        # LR from optimizer
        lr = ""
        if trainer.optimizers:
            opt = trainer.optimizers[0]
            lr = f"{opt.param_groups[0]['lr']:.2e}"

        # Gradient norm (accumulated average)
        grad_norm = ""
        if self._grad_norm_count > 0:
            grad_norm = f"{self._grad_norm_sum / self._grad_norm_count:.4f}"
            self._grad_norm_sum = 0.0
            self._grad_norm_count = 0

        # Train loss (accumulated)
        train_loss = ""
        if self._train_loss_count > 0:
            train_loss = f"{self._train_loss_sum / self._train_loss_count:.6f}"
            self._train_loss_sum = 0.0
            self._train_loss_count = 0

        # Val metrics
        val_loss = _metric_val(m, "sp_loss/val")
        val_auc = _metric_val(m, "auc/val")
        val_aupr = _metric_val(m, "aupr/val")

        # Per-family metrics
        family_aucs = [_metric_val(m, f"auc/{i}/val") for i in range(7)]
        family_auprs = [_metric_val(m, f"aupr/{i}/val") for i in range(7)]

        def fmt(v):
            return f"{v:.6f}" if v is not None else ""

        row = (
            f"{epoch},{ts},{lr},{grad_norm},"
            f"{train_loss},{fmt(val_loss)},"
            f"{fmt(val_auc)},{fmt(val_aupr)},"
            + ",".join(fmt(a) for a in family_aucs) + ","
            + ",".join(fmt(a) for a in family_auprs)
        )
        with open(self.csv_path, "a", encoding="utf-8") as f:
            f.write(row + "\n")

        # Console summary
        fam_aucs_valid = [a for a in family_aucs if a is not None]
        macro_auc = sum(fam_aucs_valid) / len(fam_aucs_valid) if fam_aucs_valid else 0
        worst_auc = min(fam_aucs_valid) if fam_aucs_valid else 0
        print(f"[Metrics] ep={epoch} auc={fmt(val_auc)} aupr={fmt(val_aupr)} "
              f"macro_auc={macro_auc:.4f} worst_auc={worst_auc:.4f} "
              f"train_loss={train_loss} val_loss={fmt(val_loss)} grad_norm={grad_norm}")


def plot_training_curves(csv_path: str = METRICS_CSV):
    """Generate comprehensive training report figures."""
    try:
        import csv as csv_mod
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        rows = []
        with open(csv_path, encoding="utf-8") as f:
            reader = csv_mod.DictReader(f)
            for row in reader:
                rows.append(row)

        if len(rows) < 2:
            print("[Plot] Not enough data points.")
            return

        def col_floats(key):
            return [float(r[key]) if r.get(key) else np.nan for r in rows]

        epochs = [int(r["epoch"]) for r in rows]
        val_auc = col_floats("val_auc")
        val_aupr = col_floats("val_aupr")
        val_loss = col_floats("val_loss")
        train_loss = col_floats("train_loss")
        lr_vals = col_floats("lr")
        grad_norm = col_floats("grad_norm")

        # Per-family AUC
        fam_aucs = {}
        for i in range(7):
            name = FAMILY_NAMES.get(str(i), str(i))
            fam_aucs[name] = col_floats(f"val_auc_{i}")

        # ---- Figure 1: Core Learning Curves (3x2) ----
        fig, axes = plt.subplots(3, 2, figsize=(14, 15))

        # 1a: Train/Val Loss
        ax = axes[0, 0]
        if not all(np.isnan(train_loss)):
            ax.plot(epochs, train_loss, "o-", label="Train", color="#FF5722", markersize=4)
        if not all(np.isnan(val_loss)):
            ax.plot(epochs, val_loss, "s-", label="Val", color="#2196F3", markersize=4)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Loss"); ax.set_title("Train / Val Loss")
        ax.legend(); ax.grid(True, alpha=0.3)

        # 1b: AUC + AUPR
        ax = axes[0, 1]
        ax.plot(epochs, val_auc, "o-", label="AUC-ROC", color="#2196F3", markersize=4)
        ax.plot(epochs, val_aupr, "s-", label="AUPR", color="#4CAF50", markersize=4)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Score"); ax.set_title("Val AUC-ROC & AUPR")
        ax.legend(); ax.grid(True, alpha=0.3)
        bi = int(np.nanargmax(val_auc))
        ax.annotate(f"Best AUC: {val_auc[bi]:.4f} (ep{epochs[bi]})",
                    xy=(epochs[bi], val_auc[bi]), xytext=(0, 12),
                    textcoords="offset points", arrowprops=dict(arrowstyle="->", color="red"),
                    color="red", fontweight="bold", ha="center", fontsize=9)

        # 1c: Learning Rate
        ax = axes[1, 0]
        if not all(np.isnan(lr_vals)):
            ax.plot(epochs, lr_vals, "D-", color="#FF9800", markersize=4, linewidth=2)
            ax.set_yscale("log")
            ax.set_xlabel("Epoch"); ax.set_ylabel("Learning Rate")
            ax.set_title("Learning Rate Schedule")
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "No LR data (pre-upgrade epochs)", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            ax.set_title("Learning Rate Schedule")

        # 1d: Gradient Norm
        ax = axes[1, 1]
        if not all(np.isnan(grad_norm)):
            ax.plot(epochs, grad_norm, "^-", color="#9C27B0", markersize=4, linewidth=2)
            ax.set_xlabel("Epoch"); ax.set_ylabel("Avg Gradient L2 Norm")
            ax.set_title("Gradient Norm (per epoch avg)")
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "No gradient norm data (pre-upgrade epochs)", ha="center",
                    va="center", transform=ax.transAxes, fontsize=12, color="gray")
            ax.set_title("Gradient Norm")

        # 1e: Generalization gap
        ax = axes[2, 0]
        if not all(np.isnan(train_loss)) and not all(np.isnan(val_loss)):
            gap = [v - t if not (np.isnan(v) or np.isnan(t)) else np.nan
                   for t, v in zip(train_loss, val_loss)]
            ax.plot(epochs, gap, "D-", color="#795548", markersize=4)
            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
            ax.set_xlabel("Epoch"); ax.set_ylabel("Val Loss - Train Loss")
            ax.set_title("Generalization Gap")
            ax.grid(True, alpha=0.3)
        else:
            ax.text(0.5, 0.5, "No train loss data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            ax.set_title("Generalization Gap")

        # 1f: Per-Family AUC over epochs
        ax = axes[2, 1]
        colors = ["#E53935", "#1E88E5", "#43A047", "#FB8C00", "#8E24AA", "#00ACC1", "#6D4C41"]
        macro_aucs = []
        for idx, (name, vals) in enumerate(fam_aucs.items()):
            valid = [v for v in vals if not np.isnan(v)]
            if valid:
                ax.plot(epochs, vals, "-", label=name, color=colors[idx % 7],
                        alpha=0.7, linewidth=1.5, markersize=3)
        for ei in range(len(epochs)):
            fam_vals = [fam_aucs[n][ei] for n in fam_aucs
                        if not np.isnan(fam_aucs[n][ei])]
            macro_aucs.append(np.mean(fam_vals) if fam_vals else np.nan)
        ax.plot(epochs, macro_aucs, "k-", linewidth=2.5, label="Macro Avg", zorder=10)
        ax.set_xlabel("Epoch"); ax.set_ylabel("AUC-ROC"); ax.set_title("Per-Family AUC Trend")
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

        plt.suptitle("Step 9: AllSplit Training Report", fontsize=16, fontweight="bold")
        plt.tight_layout()
        path1 = os.path.join(os.path.dirname(csv_path), "fig1_training_dynamics.png")
        plt.savefig(path1, dpi=150, bbox_inches="tight")
        print(f"[Plot] Fig1 saved: {path1}")
        plt.close()

        # ---- Figure 2: Per-Family Bar Chart (final epoch) ----
        fig, ax = plt.subplots(figsize=(12, 6))
        last = rows[-1]
        names, aucs_bar, auprs_bar = [], [], []
        for i in range(7):
            name = FAMILY_NAMES.get(str(i), str(i))
            a = float(last.get(f"val_auc_{i}", "nan"))
            p = float(last.get(f"val_aupr_{i}", "nan"))
            if not (np.isnan(a) and np.isnan(p)):
                names.append(name)
                aucs_bar.append(a if not np.isnan(a) else 0)
                auprs_bar.append(p if not np.isnan(p) else 0)

        if names:
            x = np.arange(len(names))
            w = 0.35
            bars1 = ax.bar(x - w/2, aucs_bar, w, label="AUC-ROC", color="#2196F3", alpha=0.8)
            bars2 = ax.bar(x + w/2, auprs_bar, w, label="AUPR", color="#4CAF50", alpha=0.8)
            ax.set_xticks(x); ax.set_xticklabels(names, rotation=30, ha="right")
            ax.set_ylabel("Score"); ax.set_title(f"Per-Family Performance (Epoch {last['epoch']})")
            ax.legend(); ax.grid(True, alpha=0.3, axis="y")
            # Value labels
            for bar in bars1:
                h = bar.get_height()
                if h > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, h + 0.01, f"{h:.3f}",
                            ha="center", va="bottom", fontsize=8)
            for bar in bars2:
                h = bar.get_height()
                if h > 0:
                    ax.text(bar.get_x() + bar.get_width()/2, h + 0.01, f"{h:.3f}",
                            ha="center", va="bottom", fontsize=8)
            # Overall lines
            oa = float(last.get("val_auc", "nan"))
            op = float(last.get("val_aupr", "nan"))
            if not np.isnan(oa):
                ax.axhline(y=oa, color="#2196F3", linestyle="--", alpha=0.5, label=f"Overall AUC={oa:.3f}")
            if not np.isnan(op):
                ax.axhline(y=op, color="#4CAF50", linestyle="--", alpha=0.5, label=f"Overall AUPR={op:.3f}")
            ax.legend(fontsize=9)

        plt.tight_layout()
        path2 = os.path.join(os.path.dirname(csv_path), "fig2_family_performance.png")
        plt.savefig(path2, dpi=150, bbox_inches="tight")
        print(f"[Plot] Fig2 saved: {path2}")
        plt.close()

    except Exception as e:
        import traceback
        print(f"[Plot] Failed: {e}")
        traceback.print_exc()


# ============================================================
# 9. Main
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
            filename=f"allsplit-fold0-{args.edge_mode}-epoch{{epoch:02d}}-auc{{auc/val:.4f}}",
            monitor="auc/val", mode="max",
            save_top_k=3, save_last=True, verbose=True,
            auto_insert_metric_name=False,
        )
        early_stop_cb = EarlyStopping(
            monitor="auc/val", mode="max", patience=15, verbose=True,
        )
        logger = TensorBoardLogger(save_dir=LOG_DIR, name=f"allsplit_fold0_{args.edge_mode}")
        metrics_csv_cb = MetricsCSVLogger()

        print("[3/4] Initializing Trainer...")
        trainer = pl.Trainer(
            benchmark=True, max_epochs=50,
            accelerator="gpu", devices=1, precision=16,
            gradient_clip_val=config.training.gradient_clip_val,
            accumulate_grad_batches=config.training.accumulate_grad_batches,
            check_val_every_n_epoch=config.training.val_frequency,
            num_sanity_val_steps=2,
            callbacks=[lr_monitor, ckpt_cb, early_stop_cb, metrics_csv_cb],
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
        print(f"Best auc/val: {ckpt_cb.best_model_score}")
        print("=" * 70)

        # ---- Auto post-training analysis ----
        # 1) Plot training curves from CSV (loss, AUC, LR, grad_norm, per-family)
        plot_training_curves()

        # 2) Detailed analysis of best checkpoint (PR/ROC curves, confusion matrix, etc.)
        best_path = ckpt_cb.best_model_path
        if best_path and os.path.exists(best_path):
            print("\n[Post-training] Running detailed eval on best checkpoint...")
            try:
                _run_best_checkpoint_analysis(model_cls=SS, config=config, dm=dm,
                                               ckpt_path=best_path, trainer=trainer)
            except Exception as e:
                import traceback
                print(f"[Post-training] Analysis failed: {e}")
                traceback.print_exc()


def _run_best_checkpoint_analysis(model_cls, config, dm, ckpt_path, trainer):
    """Run validation on best checkpoint and generate detailed figures (PR/ROC, confusion matrix)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn import metrics as skmetrics

    # Run validation to get raw logits
    model = model_cls(config)
    trainer_eval = pl.Trainer(
        accelerator="gpu", devices=1, precision=16,
        enable_progress_bar=True, logger=False,
    )
    trainer_eval.validate(model, datamodule=dm, ckpt_path=ckpt_path)

    logits = getattr(model, "logits", None)
    labels = getattr(model, "labels", None)
    if logits is None or labels is None:
        print("[Post-training] Model did not store logits/labels, skipping detailed analysis.")
        return

    import numpy as np
    logits = np.array(logits).ravel()
    labels = np.array(labels).ravel()
    probs = 1.0 / (1.0 + np.exp(-logits))

    # Extract epoch from checkpoint
    ckpt_data = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    epoch = ckpt_data.get("epoch", -1)

    output_dir = os.path.join(METRICS_DIR, f"eval_epoch{epoch:02d}")
    os.makedirs(output_dir, exist_ok=True)

    # ---- Figure: Best Checkpoint Analysis (3x2) ----
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))

    # Score distribution
    ax = axes[0, 0]
    pos_scores, neg_scores = probs[labels == 1], probs[labels == 0]
    bins = np.linspace(0, 1, 50)
    ax.hist(neg_scores, bins=bins, alpha=0.6, label=f"Neg (n={len(neg_scores)})",
            color="#2196F3", density=True)
    ax.hist(pos_scores, bins=bins, alpha=0.6, label=f"Pos (n={len(pos_scores)})",
            color="#E53935", density=True)
    ax.set_xlabel("Predicted Probability"); ax.set_ylabel("Density")
    ax.set_title("Score Distribution"); ax.legend(); ax.grid(True, alpha=0.3)

    # PR Curve
    ax = axes[0, 1]
    precision, recall, _ = skmetrics.precision_recall_curve(labels, probs)
    aupr_val = skmetrics.average_precision_score(labels, probs)
    prevalence = labels.mean()
    ax.plot(recall, precision, color="#4CAF50", linewidth=2, label=f"AUPR={aupr_val:.4f}")
    ax.axhline(y=prevalence, color="gray", linestyle="--", alpha=0.7,
                label=f"Prevalence={prevalence:.3f}")
    ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
    ax.set_title("Precision-Recall Curve"); ax.legend(); ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # ROC Curve
    ax = axes[1, 0]
    fpr, tpr, _ = skmetrics.roc_curve(labels, probs)
    auroc = skmetrics.auc(fpr, tpr)
    ax.plot(fpr, tpr, color="#2196F3", linewidth=2, label=f"AUC={auroc:.4f}")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
    ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
    ax.set_title("ROC Curve"); ax.legend(); ax.grid(True, alpha=0.3)

    # Threshold sweep
    ax = axes[1, 1]
    thresholds = np.linspace(0.01, 0.99, 200)
    f1s, precs, recs = [], [], []
    for t in thresholds:
        preds = (probs >= t).astype(int)
        tp = ((preds == 1) & (labels == 1)).sum()
        fp = ((preds == 1) & (labels == 0)).sum()
        fn = ((preds == 0) & (labels == 1)).sum()
        p = tp / (tp + fp) if (tp + fp) > 0 else 0
        r = tp / (tp + fn) if (tp + fn) > 0 else 0
        f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
        f1s.append(f1); precs.append(p); recs.append(r)
    ax.plot(thresholds, precs, label="Precision", color="#FF9800", linewidth=1.5)
    ax.plot(thresholds, recs, label="Recall", color="#03A9F4", linewidth=1.5)
    ax.plot(thresholds, f1s, label="F1", color="#4CAF50", linewidth=2)
    best_f1_idx = np.argmax(f1s)
    best_t = thresholds[best_f1_idx]
    ax.axvline(x=best_t, color="red", linestyle="--", alpha=0.5,
                label=f"Best F1={f1s[best_f1_idx]:.3f} @t={best_t:.2f}")
    ax.set_xlabel("Threshold"); ax.set_ylabel("Score")
    ax.set_title("Threshold Sweep"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

    # Confusion matrix (best F1 threshold)
    ax = axes[2, 0]
    preds = (probs >= best_t).astype(int)
    cm = skmetrics.confusion_matrix(labels, preds, labels=[0, 1])
    im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Neg", "Pos"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Neg", "Pos"])
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title(f"Confusion Matrix (threshold={best_t:.2f})")
    for i_cm in range(2):
        for j_cm in range(2):
            color = "white" if cm[i_cm, j_cm] > cm.max() / 2 else "black"
            ax.text(j_cm, i_cm, f"{cm[i_cm, j_cm]}", ha="center", va="center",
                    fontsize=16, fontweight="bold", color=color)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # Confusion matrix (0.5 threshold)
    ax = axes[2, 1]
    preds_05 = (probs >= 0.5).astype(int)
    cm05 = skmetrics.confusion_matrix(labels, preds_05, labels=[0, 1])
    im2 = ax.imshow(cm05, interpolation="nearest", cmap="Oranges")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Neg", "Pos"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Neg", "Pos"])
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title("Confusion Matrix (threshold=0.50)")
    for i_cm in range(2):
        for j_cm in range(2):
            color = "white" if cm05[i_cm, j_cm] > cm05.max() / 2 else "black"
            ax.text(j_cm, i_cm, f"{cm05[i_cm, j_cm]}", ha="center", va="center",
                    fontsize=16, fontweight="bold", color=color)
    plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)

    plt.suptitle(f"Best Checkpoint Analysis (Epoch {epoch})", fontsize=16, fontweight="bold")
    plt.tight_layout()
    p = os.path.join(output_dir, "fig_best_checkpoint_analysis.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[Post-training] Saved: {p}")


if __name__ == "__main__":
    main()
