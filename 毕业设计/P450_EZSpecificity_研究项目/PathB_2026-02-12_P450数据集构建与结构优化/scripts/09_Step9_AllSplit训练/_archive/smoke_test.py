"""
Smoke test: verify config loading, CSV column mapping, model init.
Does NOT require all feature files to be present.
"""

import os
import sys
import yaml
import time

SRC_DIR = "D:/EZSpecificity_Project/src"
sys.path.insert(0, SRC_DIR)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_PATH = os.path.join(SCRIPT_DIR, "train_allsplit_config.yml")

from easydict import EasyDict

print("[1] Loading config...")
with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
    config = EasyDict(yaml.safe_load(f))

print(f"  tag: {config.data.tag}")
print(f"  batch_size: {config.data.batch_size}")
print(f"  accumulate: {config.training.accumulate_grad_batches}")
print(f"  effective batch: {config.data.batch_size * config.training.accumulate_grad_batches}")
print(f"  sample_weight: {config.data.sample_weight}")
print(f"  high_quality_id_path: {config.data.high_quality_id_path} (type={type(config.data.high_quality_id_path).__name__})")
print()

print("[2] Testing CSV reading...")
from Datasets.utils import read_datasets, get_paths
train_df = read_datasets(config.data.train_data_path)
print(f"  Train rows: {len(train_df)}")
print(f"  Columns: {list(train_df.columns)}")
print(f"  Required columns present: ", end="")
required = ['Enzyme Index', 'Substrate Index', 'Label', 'Dock Index', 'dataset_id', 'tag']
missing = [c for c in required if c not in train_df.columns]
if missing:
    print(f"MISSING: {missing}")
else:
    print("ALL OK")
print(f"  First row: {dict(train_df.iloc[0])}")
print(f"  dataset_id values: {sorted(train_df['dataset_id'].unique())}")
print()

print("[3] Testing get_paths on feature paths...")
for name in ['enzyme_lmdb_path', 'reaction_lmdb_path', 'grover_path', 'morgan_path',
             'structure_processed_path', 'sequence_processed_path', 'high_quality_id_path']:
    val = getattr(config.data, name)
    paths = get_paths(val)
    exists = [os.path.exists(p) if p != 'None' else 'SKIP(None)' for p in paths]
    print(f"  {name}: {paths} -> exists={exists}")
print()

print("[4] Testing model init...")
import torch
print(f"  CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"  GPU: {torch.cuda.get_device_name(0)}")
    print(f"  VRAM: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.1f} GB")

from Models.ss import SS
model = SS(config)
param_count = sum(p.numel() for p in model.parameters() if p.requires_grad)
print(f"  Model created: {param_count:,} trainable parameters")
print()

print("[5] Testing high_quality_id_path handling...")
from Datasets.utils import get_paths
hq_paths = get_paths(config.data.high_quality_id_path)
print(f"  get_paths result: {hq_paths}")
for p in hq_paths:
    if p == 'None':
        print(f"  '{p}' is string 'None' -> will be caught by try/except in structure.py -> OK")
    elif os.path.exists(p):
        print(f"  '{p}' exists -> will be loaded")
    else:
        print(f"  '{p}' does not exist -> will fail in open()")
print()

print("=" * 50)
print("SMOKE TEST COMPLETE")
print("=" * 50)
