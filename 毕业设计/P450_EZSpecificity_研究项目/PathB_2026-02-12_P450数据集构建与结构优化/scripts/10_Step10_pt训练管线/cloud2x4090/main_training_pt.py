"""
Step 10: Train EZSpecificity from .pt shard cache (PtCacheDataset).

Standalone script — NO LMDB, NO monkey-patching.
Uses PtCacheDataset for all data loading (enzyme, substrate, structure).

Usage:
  # Fresh start:
  python main_training_pt.py --config train_pt_config.yml --cache-dir /path/to/allsplit_pt_cache

  # Resume:
  python main_training_pt.py --config train_pt_config.yml --cache-dir /path/to/allsplit_pt_cache --resume last
"""
from __future__ import annotations

import argparse
import csv as csv_mod
import math
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Paths — ensure src/ is importable
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", "src"))
PATHB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))

for d in (SCRIPT_DIR, SRC_DIR):
    if d not in sys.path:
        sys.path.insert(0, d)

RESULTS_DIR = os.path.join(PATHB_DIR, "results", "10_Step10_pt训练")
CKPT_DIR = os.path.join(RESULTS_DIR, "checkpoints")
LOG_DIR = os.path.join(PATHB_DIR, "logs", "10_Step10_pt训练")
# METRICS_CSV is set per-run in main() based on edge_mode, avoiding cross-run mixing
METRICS_CSV = None  # placeholder

# ---------------------------------------------------------------------------
# Torch setup (before any torch import in submodules)
# ---------------------------------------------------------------------------
import torch
import pytorch_lightning as pl
import yaml
from easydict import EasyDict
from pytorch_lightning.callbacks import Callback, EarlyStopping, LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger
from torch_geometric.loader import DataLoader

from pt_dataset import PtCacheDataset
from Models.ss import SS

torch.set_float32_matmul_precision("high")
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
torch.backends.cudnn.benchmark = True

# Fast H2D transfer patch
_orig_transfer = SS.transfer_batch_to_device


def _fast_transfer(self, batch, device, dataloader_idx):
    return batch.to(device, non_blocking=True)


SS.transfer_batch_to_device = _fast_transfer

# ---------------------------------------------------------------------------
# Family mapping
# ---------------------------------------------------------------------------
FAMILY_NAMES = {
    "0": "brenda", "1": "Duf", "2": "Esterase", "3": "Gt_acceptor",
    "4": "Nitrilase", "5": "Phosphatase", "6": "Thiolase",
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description="Train EZSpecificity with .pt shard cache")
    p.add_argument("--config", required=True, help="YAML config file")
    p.add_argument("--cache-dir", required=True, help="Path to allsplit_pt_cache directory")
    p.add_argument("--edge-mode", choices=["fixed", "legacy_bug"], default="fixed")
    p.add_argument("--num-workers", type=int, default=0,
                   help="DataLoader workers. 0=main process only (safest). "
                        "Each worker caches ~2GB of shards, so memory = (1+N)*2GB.")
    p.add_argument("--batch-size", type=int, default=None,
                   help="Override config batch size (default: from config)")
    p.add_argument("--resume", type=str, default=None,
                   help="'last' or path to checkpoint")
    p.add_argument("--max-epochs", type=int, default=50)
    p.add_argument("--devices", type=int, default=1,
                   help="Number of GPUs (1=single, 2+=DDP)")
    p.add_argument("--shutdown", action="store_true",
                   help="Auto shutdown server after training completes (for cloud instances)")
    p.add_argument("--preload", action="store_true",
                   help="Preload all .pt samples into RAM (needs ~28GB, eliminates disk I/O)")
    return p.parse_args()


def load_config(path: str) -> EasyDict:
    with open(path, "r", encoding="utf-8") as f:
        return EasyDict(yaml.safe_load(f))


# ---------------------------------------------------------------------------
# DataModule
# ---------------------------------------------------------------------------
class PtTrainingDataModule(pl.LightningDataModule):
    def __init__(self, cache_dir: str, config: EasyDict,
                 edge_mode: str, batch_size: int, num_workers: int,
                 preload: bool = False):
        super().__init__()
        self.cache_dir = cache_dir
        self.config = config
        self.edge_mode = edge_mode
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.preload = preload

    def _make_dataset(self, split: str, is_train: bool) -> PtCacheDataset:
        return PtCacheDataset(
            cache_dir=self.cache_dir,
            split=split,
            edge_mode=self.edge_mode,
            dist_noise=self.config.transform.dist_noise if is_train else False,
            cutoff=self.config.transform.cutoff,
            num_r_gaussian=self.config.transform.num_r_gaussian,
            max_enzyme_len=self.config.data.max_enzyme_length,
            max_substrate_len=self.config.data.max_substrate_length,
            preload=self.preload,
        )

    def _make_loader(self, split: str, is_train: bool, shuffle: bool = False):
        ds = self._make_dataset(split, is_train)
        kw = dict(
            batch_size=self.batch_size, shuffle=shuffle,
            num_workers=self.num_workers, pin_memory=(self.num_workers > 0),
            follow_batch=["ligand_index"],
        )
        if self.num_workers > 0:
            kw.update(persistent_workers=True, prefetch_factor=4)
        return ds, DataLoader(ds, **kw)

    def train_dataloader(self):
        t0 = time.time()
        ds, loader = self._make_loader("train", is_train=True, shuffle=True)
        print(f"[DataModule] Train: {len(ds)} samples, built in {time.time()-t0:.1f}s")
        return loader

    def val_dataloader(self):
        t0 = time.time()
        ds, loader = self._make_loader("val", is_train=False)
        print(f"[DataModule] Val: {len(ds)} samples, built in {time.time()-t0:.1f}s")
        return loader

    def test_dataloader(self):
        t0 = time.time()
        ds, loader = self._make_loader("test", is_train=False)
        print(f"[DataModule] Test: {len(ds)} samples, built in {time.time()-t0:.1f}s")
        return loader


# ---------------------------------------------------------------------------
# MetricsCSVLogger
# ---------------------------------------------------------------------------
def _metric_val(m, key):
    v = m.get(key)
    if v is None:
        return None
    v = v.item() if hasattr(v, "item") else float(v)
    return v if not math.isnan(v) else None


class MetricsCSVLogger(Callback):
    HEADER = (
        "epoch,timestamp,lr,grad_norm,"
        "train_loss,val_loss,"
        "val_auc,val_aupr,"
        "val_auc_0,val_auc_1,val_auc_2,val_auc_3,val_auc_4,val_auc_5,val_auc_6,"
        "val_aupr_0,val_aupr_1,val_aupr_2,val_aupr_3,val_aupr_4,val_aupr_5,val_aupr_6"
    )

    def __init__(self, csv_path: str):
        super().__init__()
        self.csv_path = csv_path
        self._train_loss_sum = 0.0
        self._train_loss_count = 0
        self._grad_norm_sum = 0.0
        self._grad_norm_count = 0
        self._ensure_csv_header(csv_path)

    @staticmethod
    def _ensure_csv_header(csv_path):
        expected_cols = MetricsCSVLogger.HEADER.replace("\n", "").split(",")
        if not os.path.exists(csv_path):
            os.makedirs(os.path.dirname(csv_path), exist_ok=True)
            with open(csv_path, "w", encoding="utf-8") as f:
                f.write(MetricsCSVLogger.HEADER + "\n")
            return
        with open(csv_path, "r", encoding="utf-8") as f:
            first_line = f.readline().strip()
        existing_cols = first_line.split(",")
        if existing_cols == expected_cols:
            return
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
        print(f"[MetricsCSV] Migrated {len(rows)} rows to new schema")

    def on_train_batch_end(self, trainer, pl_module, outputs, batch, batch_idx):
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
        # O(1) memory: accumulate per-parameter norm, no giant cat tensor
        total_norm_sq = 0.0
        for p in pl_module.parameters():
            if p.grad is not None:
                total_norm_sq += p.grad.detach().norm(2).item() ** 2
        if total_norm_sq > 0:
            self._grad_norm_sum += total_norm_sq ** 0.5
            self._grad_norm_count += 1

    def on_validation_epoch_end(self, trainer, pl_module):
        if trainer.sanity_checking:
            return
        if trainer.global_rank != 0:
            return
        m = trainer.callback_metrics
        epoch = trainer.current_epoch
        ts = time.strftime("%Y-%m-%d %H:%M:%S")

        lr = ""
        if trainer.optimizers:
            lr = f"{trainer.optimizers[0].param_groups[0]['lr']:.2e}"

        grad_norm = ""
        if self._grad_norm_count > 0:
            grad_norm = f"{self._grad_norm_sum / self._grad_norm_count:.4f}"
            self._grad_norm_sum = 0.0
            self._grad_norm_count = 0

        train_loss = ""
        if self._train_loss_count > 0:
            train_loss = f"{self._train_loss_sum / self._train_loss_count:.6f}"
            self._train_loss_sum = 0.0
            self._train_loss_count = 0

        val_loss = _metric_val(m, "sp_loss/val")
        val_auc = _metric_val(m, "auc/val")
        val_aupr = _metric_val(m, "aupr/val")
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

        fam_aucs_valid = [a for a in family_aucs if a is not None]
        macro_auc = sum(fam_aucs_valid) / len(fam_aucs_valid) if fam_aucs_valid else 0
        worst_auc = min(fam_aucs_valid) if fam_aucs_valid else 0
        print(f"[Metrics] ep={epoch} auc={fmt(val_auc)} aupr={fmt(val_aupr)} "
              f"macro_auc={macro_auc:.4f} worst_auc={worst_auc:.4f} "
              f"train_loss={train_loss} val_loss={fmt(val_loss)} grad_norm={grad_norm}")

        plot_training_curves(self.csv_path)


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def plot_training_curves(csv_path: str):
    try:
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

        fam_aucs = {}
        for i in range(7):
            fam_aucs[FAMILY_NAMES.get(str(i), str(i))] = col_floats(f"val_auc_{i}")

        fig, axes = plt.subplots(3, 2, figsize=(14, 15))

        # Train/Val Loss
        ax = axes[0, 0]
        if not all(np.isnan(train_loss)):
            ax.plot(epochs, train_loss, "o-", label="Train", color="#FF5722", markersize=4)
        if not all(np.isnan(val_loss)):
            ax.plot(epochs, val_loss, "s-", label="Val", color="#2196F3", markersize=4)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Loss"); ax.set_title("Train / Val Loss")
        ax.legend(); ax.grid(True, alpha=0.3)

        # AUC + AUPR
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

        # LR
        ax = axes[1, 0]
        if not all(np.isnan(lr_vals)):
            ax.plot(epochs, lr_vals, "D-", color="#FF9800", markersize=4, linewidth=2)
            ax.set_yscale("log")
        ax.set_xlabel("Epoch"); ax.set_ylabel("Learning Rate")
        ax.set_title("Learning Rate Schedule"); ax.grid(True, alpha=0.3)

        # Grad norm
        ax = axes[1, 1]
        if not all(np.isnan(grad_norm)):
            ax.plot(epochs, grad_norm, "^-", color="#9C27B0", markersize=4, linewidth=2)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Avg Gradient L2 Norm")
        ax.set_title("Gradient Norm"); ax.grid(True, alpha=0.3)

        # Generalization gap
        ax = axes[2, 0]
        if not all(np.isnan(train_loss)) and not all(np.isnan(val_loss)):
            gap = [v - t if not (np.isnan(v) or np.isnan(t)) else np.nan
                   for t, v in zip(train_loss, val_loss)]
            ax.plot(epochs, gap, "D-", color="#795548", markersize=4)
            ax.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Val - Train Loss")
        ax.set_title("Generalization Gap"); ax.grid(True, alpha=0.3)

        # Per-family AUC
        ax = axes[2, 1]
        colors = ["#E53935", "#1E88E5", "#43A047", "#FB8C00", "#8E24AA", "#00ACC1", "#6D4C41"]
        macro_aucs = []
        for idx, (name, vals) in enumerate(fam_aucs.items()):
            if any(not np.isnan(v) for v in vals):
                ax.plot(epochs, vals, "-", label=name, color=colors[idx % 7],
                        alpha=0.7, linewidth=1.5)
        for ei in range(len(epochs)):
            fam_vals = [fam_aucs[n][ei] for n in fam_aucs if not np.isnan(fam_aucs[n][ei])]
            macro_aucs.append(np.mean(fam_vals) if fam_vals else np.nan)
        ax.plot(epochs, macro_aucs, "k-", linewidth=2.5, label="Macro Avg", zorder=10)
        ax.set_xlabel("Epoch"); ax.set_ylabel("AUC-ROC"); ax.set_title("Per-Family AUC Trend")
        ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

        plt.suptitle("Step 10: .pt Cache Training Report", fontsize=16, fontweight="bold")
        plt.tight_layout()
        out_path = os.path.join(os.path.dirname(csv_path), "fig1_training_dynamics.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close()

        # Per-family bar chart (last epoch)
        fig, ax = plt.subplots(figsize=(12, 6))
        last = rows[-1]
        names, aucs_bar, auprs_bar = [], [], []
        for i in range(7):
            name = FAMILY_NAMES.get(str(i), str(i))
            a = float(last.get(f"val_auc_{i}", "nan") or "nan")
            p = float(last.get(f"val_aupr_{i}", "nan") or "nan")
            if not (np.isnan(a) and np.isnan(p)):
                names.append(name)
                aucs_bar.append(a if not np.isnan(a) else 0)
                auprs_bar.append(p if not np.isnan(p) else 0)
        if names:
            x = np.arange(len(names)); w = 0.35
            bars1 = ax.bar(x - w/2, aucs_bar, w, label="AUC-ROC", color="#2196F3", alpha=0.8)
            bars2 = ax.bar(x + w/2, auprs_bar, w, label="AUPR", color="#4CAF50", alpha=0.8)
            ax.set_xticks(x); ax.set_xticklabels(names, rotation=30, ha="right")
            ax.set_ylabel("Score"); ax.set_title(f"Per-Family Performance (Epoch {last['epoch']})")
            ax.legend(); ax.grid(True, alpha=0.3, axis="y")
            for bars in (bars1, bars2):
                for bar in bars:
                    h = bar.get_height()
                    if h > 0:
                        ax.text(bar.get_x() + bar.get_width()/2, h + 0.01, f"{h:.3f}",
                                ha="center", va="bottom", fontsize=8)
            for key, color in [("val_auc", "#2196F3"), ("val_aupr", "#4CAF50")]:
                v = float(last.get(key, "nan") or "nan")
                if not np.isnan(v):
                    ax.axhline(y=v, color=color, linestyle="--", alpha=0.5,
                               label=f"Overall {key.split('_')[1].upper()}={v:.3f}")
            ax.legend(fontsize=9)
        plt.tight_layout()
        out_path2 = os.path.join(os.path.dirname(csv_path), "fig2_family_performance.png")
        plt.savefig(out_path2, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"[Plot] Saved: {out_path}, {out_path2}")

    except Exception as e:
        import traceback
        print(f"[Plot] Failed: {e}")
        traceback.print_exc()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    config = load_config(args.config)

    batch_size = args.batch_size or config.data.batch_size
    accum = config.training.accumulate_grad_batches
    seed = getattr(config, "seed", 3407)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    os.makedirs(CKPT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)

    # Validate cache dir
    cache_dir = os.path.abspath(args.cache_dir)
    for required in ["train", "val", "enzymes", "substrates"]:
        p = os.path.join(cache_dir, required)
        if not os.path.isdir(p):
            print(f"[ERROR] Missing directory: {p}")
            print(f"  Expected .pt cache structure under: {cache_dir}")
            sys.exit(1)

    # Banner
    print("=" * 70)
    print("Step 10: EZSpecificity Training (.pt shard cache)")
    print(f"  Edge mode   : {args.edge_mode}")
    print(f"  Cache dir   : {cache_dir}")
    print(f"  Config      : {args.config}")
    print(f"  Batch       : {batch_size} x {accum} = {batch_size * accum} effective")
    print(f"  Workers     : {args.num_workers}")
    print(f"  Max epochs  : {args.max_epochs}")
    print(f"  Seed        : {seed}")
    n_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 0
    gpu_info = f"{torch.cuda.get_device_name(0)} x{args.devices}" if n_gpus > 0 else "CPU"
    print(f"  GPU         : {gpu_info} (devices={args.devices})")
    print(f"  CUDA        : {torch.version.cuda}, PyTorch: {torch.__version__}")
    print(f"  Results     : {RESULTS_DIR}")
    print(f"  Checkpoints : {CKPT_DIR}")
    print("=" * 70)

    pl.seed_everything(seed)

    # DataModule
    print("[1/4] Initializing DataModule...")
    dm = PtTrainingDataModule(
        cache_dir=cache_dir, config=config,
        edge_mode=args.edge_mode, batch_size=batch_size,
        num_workers=args.num_workers,
        preload=args.preload,
    )

    # Model
    print("[2/4] Initializing Model (random weights)...")
    model = SS(config)
    # Patch self.log: always sync_dist=False because our gather patch
    # already computes global metrics on every rank — no need for Lightning sync.
    _orig_log = model.log
    def _nosync_log(name, value, **kw):
        kw["sync_dist"] = False
        return _orig_log(name, value, **kw)
    model.log = _nosync_log

    # Patch validation_epoch_end: gather outputs from all DDP ranks before computing AUC
    _orig_val_epoch_end = model.validation_epoch_end
    def _ddp_val_epoch_end(outputs):
        import numpy as np
        import torch.distributed as dist
        if dist.is_initialized() and dist.get_world_size() > 1:
            # Each rank has its own subset of validation outputs
            # Gather all logits/labels to every rank so AUC is computed on full data
            local_logits = np.concatenate([a[0] for a in outputs], axis=0).ravel()
            local_labels = np.concatenate([a[1] for a in outputs], axis=0).ravel()
            local_logits_t = torch.tensor(local_logits, dtype=torch.float32).cuda()
            local_labels_t = torch.tensor(local_labels, dtype=torch.float32).cuda()
            # Gather sizes (ranks may have slightly different counts)
            local_n = torch.tensor([local_logits_t.shape[0]], device=local_logits_t.device)
            world = dist.get_world_size()
            all_n = [torch.zeros_like(local_n) for _ in range(world)]
            dist.all_gather(all_n, local_n)
            max_n = max(n.item() for n in all_n)
            # Pad and gather
            pad_logits = torch.zeros(max_n, device=local_logits_t.device)
            pad_labels = torch.zeros(max_n, device=local_labels_t.device)
            pad_logits[:local_logits_t.shape[0]] = local_logits_t
            pad_labels[:local_labels_t.shape[0]] = local_labels_t
            gathered_logits = [torch.zeros_like(pad_logits) for _ in range(world)]
            gathered_labels = [torch.zeros_like(pad_labels) for _ in range(world)]
            dist.all_gather(gathered_logits, pad_logits)
            dist.all_gather(gathered_labels, pad_labels)
            # Reconstruct full outputs: merge into a single fake "output" tuple
            all_logits = torch.cat([g[:all_n[i].item()] for i, g in enumerate(gathered_logits)])
            all_labels = torch.cat([g[:all_n[i].item()] for i, g in enumerate(gathered_labels)])
            # Also gather tags (list of tuples of tag_0, tag_1 strings)
            # Tags are needed for per-family AUC; gather via pickle
            import pickle
            local_tags = []
            for out in outputs:
                for t0, t1 in zip(out[2][0], out[2][1]):
                    local_tags.append((t0, t1))
            tag_bytes = pickle.dumps(local_tags)
            tag_tensor = torch.tensor(list(tag_bytes), dtype=torch.uint8).cuda()
            tag_sizes = [torch.zeros(1, dtype=torch.long, device=tag_tensor.device) for _ in range(world)]
            dist.all_gather(tag_sizes, torch.tensor([len(tag_bytes)], dtype=torch.long, device=tag_tensor.device))
            max_tag = max(s.item() for s in tag_sizes)
            pad_tags = torch.zeros(max_tag, dtype=torch.uint8, device=tag_tensor.device)
            pad_tags[:len(tag_bytes)] = tag_tensor
            gathered_tags = [torch.zeros_like(pad_tags) for _ in range(world)]
            dist.all_gather(gathered_tags, pad_tags)
            all_tags = []
            for i in range(world):
                t = bytes(gathered_tags[i][:tag_sizes[i].item()].cpu().tolist())
                all_tags.extend(pickle.loads(t))
            # Build merged outputs as single entry
            tag0_list = [t[0] for t in all_tags]
            tag1_list = [t[1] for t in all_tags]
            merged_outputs = [(all_logits.cpu().numpy(), all_labels.cpu().numpy(), (tag0_list, tag1_list))]
            _orig_val_epoch_end(merged_outputs)
        else:
            _orig_val_epoch_end(outputs)
    model.validation_epoch_end = _ddp_val_epoch_end
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Trainable parameters: {n_params:,}")

    # Callbacks
    lr_monitor = LearningRateMonitor(logging_interval="step")
    ckpt_cb = ModelCheckpoint(
        dirpath=os.path.join(CKPT_DIR, args.edge_mode),
        filename=f"pt-{args.edge_mode}-ep{{epoch:02d}}-auc{{auc/val:.4f}}",
        monitor="auc/val", mode="max",
        save_top_k=3, save_last=True, verbose=True,
        auto_insert_metric_name=False,
    )
    early_stop_cb = EarlyStopping(
        monitor="auc/val", mode="max", patience=15, verbose=True,
    )
    # Per-run metrics CSV (edge_mode in filename to avoid cross-run mixing)
    run_metrics_csv = os.path.join(RESULTS_DIR, f"metrics_{args.edge_mode}.csv")
    metrics_csv_cb = MetricsCSVLogger(csv_path=run_metrics_csv)
    logger = TensorBoardLogger(save_dir=LOG_DIR, name=f"pt_{args.edge_mode}")

    # Trainer
    print("[3/4] Initializing Trainer...")
    n_devices = args.devices if torch.cuda.is_available() else 1
    strategy = pl.strategies.DDPStrategy(
        find_unused_parameters=False,
        gradient_as_bucket_view=True,
    ) if n_devices > 1 else None
    trainer = pl.Trainer(
        max_epochs=args.max_epochs,
        accelerator="gpu" if torch.cuda.is_available() else "cpu",
        devices=n_devices,
        strategy=strategy,
        precision=16 if torch.cuda.is_available() else 32,
        gradient_clip_val=config.training.gradient_clip_val,
        accumulate_grad_batches=accum,
        check_val_every_n_epoch=config.training.val_frequency,
        num_sanity_val_steps=2,
        callbacks=[lr_monitor, ckpt_cb, early_stop_cb, metrics_csv_cb],
        logger=logger,
        enable_progress_bar=True, log_every_n_steps=50,
    )

    # Resume
    ckpt_path = None
    if args.resume:
        if args.resume.lower() == "last":
            last_ckpt = os.path.join(CKPT_DIR, args.edge_mode, "last.ckpt")
            if os.path.exists(last_ckpt):
                ckpt_path = last_ckpt
                print(f"[4/4] Resuming from: {ckpt_path}")
            else:
                print(f"[WARN] last.ckpt not found at {last_ckpt}, training from scratch")
        elif os.path.exists(args.resume):
            ckpt_path = args.resume
            print(f"[4/4] Resuming from: {ckpt_path}")
        else:
            print(f"[ERROR] Checkpoint not found: {args.resume}")
            sys.exit(1)

    # Train
    print("[4/4] Starting training...")
    print(f"  Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    t0 = time.time()

    trainer.fit(model, datamodule=dm, ckpt_path=ckpt_path)

    elapsed = time.time() - t0
    print()
    print("=" * 70)
    print("Training complete!")
    print(f"  End          : {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Time         : {elapsed:.0f}s = {elapsed/60:.1f}min = {elapsed/3600:.2f}h")
    n_batches = trainer.num_training_batches
    n_epochs = trainer.current_epoch
    if n_batches > 0 and n_epochs > 0:
        total = n_batches * n_epochs
        print(f"  Batches/ep   : {n_batches}, Epochs: {n_epochs}, "
              f"Total: {total}, avg {elapsed/total:.3f}s/batch")
    print(f"  Best ckpt    : {ckpt_cb.best_model_path}")
    print(f"  Best auc/val : {ckpt_cb.best_model_score}")
    print("=" * 70)

    plot_training_curves(run_metrics_csv)

    # Auto-shutdown for cloud instances
    if args.shutdown:
        print("\n[AUTO-SHUTDOWN] Training complete. Shutting down in 30s... (Ctrl+C to cancel)")
        try:
            time.sleep(30)
            os.system("shutdown -h now")
        except KeyboardInterrupt:
            print("[AUTO-SHUTDOWN] Cancelled.")


if __name__ == "__main__":
    main()
