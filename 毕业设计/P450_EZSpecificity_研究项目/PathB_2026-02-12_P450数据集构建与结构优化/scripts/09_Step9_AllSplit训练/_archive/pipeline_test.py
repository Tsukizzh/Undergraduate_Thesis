"""
Pipeline test: verify full training loop with a few batches.
Uses Google Drive for enzyme_features.lmdb if local copy is incomplete.
"""

import os
import sys
import time
import lmdb
import yaml
import warnings

warnings.filterwarnings("ignore")

# Monkey-patch lmdb.open
_original_lmdb_open = lmdb.open
def _patched_lmdb_open(path, **kwargs):
    if kwargs.get('readonly', False) or not kwargs.get('create', True):
        kwargs['map_size'] = 128 * (1024 ** 3)
    return _original_lmdb_open(path, **kwargs)
lmdb.open = _patched_lmdb_open

SRC_DIR = "D:/EZSpecificity_Project/src"
sys.path.insert(0, SRC_DIR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(SCRIPT_DIR, "train_allsplit_config.yml")

from easydict import EasyDict

with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
    config = EasyDict(yaml.safe_load(f))

# Check if local enzyme file is complete (~51GB)
LOCAL_ENZYME = config.data.enzyme_lmdb_path
GDRIVE_ENZYME = "G:/.shortcut-targets-by-id/173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr/ESIBank/brenda/enzyme_features.lmdb"

local_size = os.path.getsize(LOCAL_ENZYME) if os.path.exists(LOCAL_ENZYME) else 0
if local_size < 50 * (1024**3):  # less than 50GB = incomplete
    print(f"[WARN] Local enzyme file is {local_size/1024**3:.1f}GB (incomplete). Using Google Drive version.")
    config.data.enzyme_lmdb_path = GDRIVE_ENZYME
else:
    print(f"[OK] Local enzyme file is {local_size/1024**3:.1f}GB (complete).")

# Same for grover
LOCAL_GROVER = config.data.grover_path
GDRIVE_GROVER = "G:/.shortcut-targets-by-id/173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr/ESIBank/brenda/grover_fingerprint.lmdb"

local_grover_size = os.path.getsize(LOCAL_GROVER) if os.path.exists(LOCAL_GROVER) else 0
if local_grover_size < 8 * (1024**3):  # less than 8GB = incomplete
    print(f"[WARN] Local grover file is {local_grover_size/1024**3:.1f}GB (incomplete). Using Google Drive version.")
    config.data.grover_path = GDRIVE_GROVER
else:
    print(f"[OK] Local grover file is {local_grover_size/1024**3:.1f}GB (complete).")

import torch
import pytorch_lightning as pl
from pytorch_lightning.loggers import TensorBoardLogger

from torch_geometric.data import DataLoader
from torch_geometric.transforms import Compose

from Datasets.data_representer import get_representer
import Datasets.Structure.transforms as utils_trans
from Datasets.utils import read_datasets
from Models.ss import SS


class TrainingDataModule(pl.LightningDataModule):
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.batch_size = config.data.batch_size

        self.train_df = read_datasets(config.data.train_data_path)
        self.val_df = read_datasets(config.data.val_data_path)
        self.test_df = read_datasets(config.data.test_data_path)

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

    def train_dataloader(self):
        data = get_representer(
            df=self.train_df, config=self.config,
            transform=self.train_transform, is_train=True
        )
        return DataLoader(
            data, batch_size=self.batch_size, shuffle=True,
            num_workers=0, pin_memory=True,
            follow_batch=['ligand_index']
        )

    def val_dataloader(self):
        data = get_representer(
            df=self.val_df, config=self.config,
            transform=self.val_transform, is_train=False
        )
        return DataLoader(
            data, batch_size=self.batch_size, shuffle=False,
            num_workers=0, pin_memory=True,
            follow_batch=['ligand_index']
        )


def main():
    print("=" * 60)
    print("Pipeline Test: 20 train batches + 5 val batches")
    print("=" * 60)

    pl.seed_everything(config.seed)

    print("[1/3] Building DataModule...")
    t0 = time.time()
    dm = TrainingDataModule(config)
    print(f"  DataModule init: {time.time()-t0:.1f}s")

    print("[2/3] Initializing Model...")
    model = SS(config)

    LOG_DIR = os.path.join(SCRIPT_DIR, "..", "..", "sessions", "09_Step9_AllSplit训练", "logs")
    os.makedirs(LOG_DIR, exist_ok=True)
    logger = TensorBoardLogger(save_dir=LOG_DIR, name='pipeline_test')

    print("[3/3] Running 20 train + 5 val batches...")
    trainer = pl.Trainer(
        max_epochs=1,
        accelerator='gpu',
        devices=1,
        precision=16,
        gradient_clip_val=config.training.gradient_clip_val,
        accumulate_grad_batches=config.training.accumulate_grad_batches,
        limit_train_batches=20,
        limit_val_batches=5,
        num_sanity_val_steps=0,
        logger=logger,
        enable_progress_bar=True,
        log_every_n_steps=5,
    )

    t_start = time.time()
    trainer.fit(model, datamodule=dm)
    elapsed = time.time() - t_start

    print()
    print("=" * 60)
    print(f"Pipeline test complete! Time: {elapsed:.1f}s")
    print(f"  20 train batches: ~{elapsed:.1f}s total")
    print(f"  Per batch: ~{elapsed/20:.2f}s")
    print(f"  Estimated full epoch ({len(dm.train_df)}÷{config.data.batch_size}≈{len(dm.train_df)//config.data.batch_size} batches):")
    est_epoch = (elapsed / 20) * (len(dm.train_df) // config.data.batch_size)
    print(f"    ~{est_epoch:.0f}s = {est_epoch/60:.1f}min = {est_epoch/3600:.2f}h")
    print("=" * 60)


if __name__ == '__main__':
    main()
