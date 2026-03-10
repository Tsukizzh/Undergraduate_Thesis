"""
Step 9: Train EZSpecificity from scratch with all_split (unknown enzyme+substrate).
BRENDA + 7 small families, fold 0. Full training pipeline.

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe ../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/main_training.py
"""

import os
import sys
import time
from queue import Queue
from threading import Thread

import lmdb
import yaml
import warnings

warnings.filterwarnings("ignore")

# ============================================================
# 1. Monkey-patch lmdb.open: increase map_size for large DBs
# ============================================================
_original_lmdb_open = lmdb.open

def _patched_lmdb_open(path, **kwargs):
    if kwargs.get('readonly', False) or not kwargs.get('create', True):
        kwargs['map_size'] = 128 * (1024 ** 3)
    return _original_lmdb_open(path, **kwargs)

lmdb.open = _patched_lmdb_open

# ============================================================
# 2. Setup paths
# ============================================================
SRC_DIR = "D:/EZSpecificity_Project/src"
sys.path.insert(0, SRC_DIR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PATHB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
CONFIG_PATH = os.path.join(SCRIPT_DIR, "train_allsplit_config.yml")
LOG_DIR = os.path.join(PATHB_DIR, "logs", "09_Step9_AllSplit训练")
CKPT_DIR = os.path.join(PATHB_DIR, "results", "09_Step9_AllSplit训练", "checkpoints")

os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(CKPT_DIR, exist_ok=True)

# ============================================================
# 3. Load config
# ============================================================
from easydict import EasyDict

with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
    config = EasyDict(yaml.safe_load(f))

# ============================================================
# 4. Imports (after sys.path is set)
# ============================================================
import torch
import pytorch_lightning as pl
from pytorch_lightning.callbacks import LearningRateMonitor, ModelCheckpoint
from pytorch_lightning.loggers import TensorBoardLogger

from torch_geometric.data import DataLoader
from torch_geometric.transforms import Compose

from Datasets.data_representer import get_representer
import Datasets.Structure.transforms as utils_trans
from Datasets.utils import read_datasets
from Models.ss import SS

# TF32 + cuDNN benchmark
torch.set_float32_matmul_precision("high")
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True
torch.backends.cudnn.benchmark = True

# ============================================================
# 4b. Monkey-patch: non-blocking H2D transfer
# ============================================================
_original_transfer = SS.transfer_batch_to_device
def _fast_transfer(self, batch, device, dataloader_idx):
    return batch.to(device, non_blocking=True)
SS.transfer_batch_to_device = _fast_transfer


class BackgroundPrefetchLoader:
    """Wraps a DataLoader to prefetch batches in a background thread."""
    def __init__(self, loader, max_prefetch=2):
        self.loader = loader
        self.max_prefetch = max_prefetch

    def __len__(self):
        return len(self.loader)

    def __getattr__(self, name):
        return getattr(self.loader, name)

    def __iter__(self):
        queue = Queue(maxsize=self.max_prefetch)
        sentinel = object()
        errors = []

        def _producer():
            try:
                for batch in self.loader:
                    queue.put(batch)
            except Exception as exc:
                errors.append(exc)
            finally:
                queue.put(sentinel)

        thread = Thread(target=_producer, daemon=True)
        thread.start()

        while True:
            item = queue.get()
            if item is sentinel:
                if errors:
                    raise errors[0]
                break
            yield item


# ============================================================
# 5. LMDB handle management for Windows multiprocessing
# ============================================================
NUM_WORKERS = 2
print(f"[Config] num_workers={NUM_WORKERS}")

def _close_lmdb_handles(dataset):
    """Close LMDB handles after valid_idx is built, enabling pickle for spawn workers."""
    seq_db = dataset.sequence_db
    for attr in ['grover_dbs', 'enzyme_dbs', 'reaction_dbs']:
        dbs = getattr(seq_db, attr, None)
        if dbs is not None:
            for db in dbs:
                if hasattr(db, 'close'):
                    db.close()
            setattr(seq_db, attr, None)

    str_db = dataset.structure_db
    for attr in ['complex_dbs', 'ligand_dbs']:
        dbs = getattr(str_db, attr, None)
        if dbs is not None:
            for db in dbs:
                if db is not None and hasattr(db, 'close'):
                    db.close()
            setattr(str_db, attr, None)


# ============================================================
# 6. Custom DataModule
# ============================================================
class TrainingDataModule(pl.LightningDataModule):
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.batch_size = config.data.batch_size

        print(f"[DataModule] Loading CSVs...")
        t0 = time.time()
        self.train_df = read_datasets(config.data.train_data_path)
        self.val_df = read_datasets(config.data.val_data_path)
        self.test_df = read_datasets(config.data.test_data_path)
        print(f"[DataModule] CSVs loaded in {time.time()-t0:.1f}s")
        print(f"[DataModule] Train: {len(self.train_df)}, Val: {len(self.val_df)}, Test: {len(self.test_df)}")

        self.train_transform = self._get_transform(is_train=True)
        self.val_transform = self._get_transform(is_train=False)

    def _get_transform(self, is_train=True):
        return Compose([
            utils_trans.FeaturizeProteinAtom(),
            utils_trans.FeaturizeLigandAtom(),
            utils_trans.EdgeConnection(
                dist_noise=self.config.transform.dist_noise if is_train else False,
                cutoff=self.config.transform.cutoff,
                num_r_gaussian=self.config.transform.num_r_gaussian,
                k=self.config.transform.k
            )
        ])

    def _make_loader(self, df, transform, is_train, shuffle=False):
        data = get_representer(
            df=df, config=self.config,
            transform=transform, is_train=is_train
        )
        if NUM_WORKERS > 0:
            _close_lmdb_handles(data)
        loader_kwargs = dict(
            batch_size=self.batch_size,
            shuffle=shuffle,
            num_workers=NUM_WORKERS, pin_memory=True,
            follow_batch=['ligand_index']
        )
        if NUM_WORKERS > 0:
            loader_kwargs.update(persistent_workers=True, prefetch_factor=2)
        return data, DataLoader(data, **loader_kwargs)

    def train_dataloader(self):
        print(f"[DataModule] Building train dataset...")
        t0 = time.time()
        data, loader = self._make_loader(
            self.train_df, self.train_transform, is_train=True, shuffle=True
        )
        elapsed = time.time() - t0
        print(f"[DataModule] Train dataset built in {elapsed:.1f}s, valid samples: {len(data)}")
        return BackgroundPrefetchLoader(loader)

    def val_dataloader(self):
        print(f"[DataModule] Building val dataset...")
        t0 = time.time()
        data, loader = self._make_loader(
            self.val_df, self.val_transform, is_train=False
        )
        elapsed = time.time() - t0
        print(f"[DataModule] Val dataset built in {elapsed:.1f}s, valid samples: {len(data)}")
        return loader

    def test_dataloader(self):
        print(f"[DataModule] Building test dataset...")
        t0 = time.time()
        data, loader = self._make_loader(
            self.test_df, self.val_transform, is_train=False
        )
        elapsed = time.time() - t0
        print(f"[DataModule] Test dataset built in {elapsed:.1f}s, valid samples: {len(data)}")
        return loader


# ============================================================
# 7. Main
# ============================================================
def main():
    print("=" * 70)
    print("Step 9: Training EZSpecificity (all_split, fold 0, BRENDA + 7 families)")
    print("  >> Full pipeline: train + val + checkpoint + scheduler")
    print("=" * 70)
    print(f"Config: {CONFIG_PATH}")
    print(f"GPU: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")
    print(f"CUDA: {torch.version.cuda}")
    print(f"PyTorch: {torch.__version__}")
    print(f"Batch size: {config.data.batch_size}, Accumulate: {config.training.accumulate_grad_batches}")
    print(f"Effective batch: {config.data.batch_size * config.training.accumulate_grad_batches}")
    print()

    pl.seed_everything(config.seed)

    # DataModule
    print("[1/4] Initializing DataModule...")
    dm = TrainingDataModule(config)

    # Model
    print("[2/4] Initializing Model (random weights)...")
    model = SS(config)
    param_count = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[2/4] Trainable parameters: {param_count:,}")

    # Callbacks
    lr_monitor = LearningRateMonitor(logging_interval='step')
    checkpoint_cb = ModelCheckpoint(
        dirpath=CKPT_DIR,
        filename='allsplit-fold0-epoch{epoch:02d}-aupr{aupr/val:.4f}',
        monitor='aupr/val',
        mode='max',
        save_top_k=3,
        save_last=True,
        verbose=True,
        auto_insert_metric_name=False,
    )
    logger = TensorBoardLogger(save_dir=LOG_DIR, name='allsplit_fold0')

    # Trainer
    print("[3/4] Initializing Trainer...")
    trainer = pl.Trainer(
        benchmark=True,
        max_epochs=50,
        accelerator='gpu',
        devices=1,
        precision=16,
        gradient_clip_val=config.training.gradient_clip_val,
        accumulate_grad_batches=config.training.accumulate_grad_batches,
        check_val_every_n_epoch=config.training.val_frequency,
        num_sanity_val_steps=2,
        callbacks=[lr_monitor, checkpoint_cb],
        logger=logger,
        enable_progress_bar=True,
        log_every_n_steps=50,
    )

    # Train
    print("[4/4] Starting training...")
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    t_start = time.time()

    trainer.fit(model, datamodule=dm)

    t_end = time.time()
    elapsed = t_end - t_start
    print()
    print("=" * 70)
    print(f"Training complete!")
    print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Total time: {elapsed:.1f}s = {elapsed/60:.1f}min = {elapsed/3600:.2f}h")
    n_batches = trainer.num_training_batches
    print(f"Batches per epoch: {n_batches}")
    if n_batches > 0:
        print(f"Avg sec/batch: {elapsed/n_batches:.3f}")
    print(f"Best checkpoint: {checkpoint_cb.best_model_path}")
    print(f"Best aupr/val: {checkpoint_cb.best_model_score}")
    print("=" * 70)


if __name__ == '__main__':
    main()
