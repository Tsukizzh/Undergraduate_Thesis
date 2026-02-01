"""
Step 9: Windows-Safe EZSpecificity Inference Script
=====================================================
This script runs inference on the P450 test dataset using the pre-trained model.

Key Windows Compatibility Fixes:
1. No torch.multiprocessing.set_sharing_strategy('file_descriptor') - Linux only
2. No DDP strategy - use pure forward pass instead of trainer.test()
3. num_workers=0 - avoid Windows multiprocessing issues
4. full_data=False - we don't have str_features.lmdb

Usage:
    D:/anaconda3/envs/torch/python.exe step9_inference.py

Output:
    data/09_Step9_模型推理/predictions.csv
"""

import sys
import os
from pathlib import Path

# ============================================================================
# Path Configuration
# ============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
PATHA_ROOT = SCRIPT_DIR.parent.parent  # PathA_2026-01-08_模型评估测试集构建
PROJECT_ROOT = PATHA_ROOT.parent.parent  # EZSpecificity_Project
SRC_DIR = PROJECT_ROOT / "src"

# Add src to sys.path for imports
sys.path.insert(0, str(SRC_DIR))

# ============================================================================
# Feature File Paths
# ============================================================================
DATA_ROOT = PATHA_ROOT / "data"

# Input files
TEST_CSV = DATA_ROOT / "04_Step4_格式修正后数据" / "data.csv"
ENZYME_LMDB = DATA_ROOT / "05_Step5_ESM酶嵌入" / "enzyme_features.lmdb"
REACTION_LMDB = DATA_ROOT / "06_Step6_反应图特征" / "reaction_features.lmdb"
GROVER_LMDB = DATA_ROOT / "07_Step7_分子指纹" / "grover_fingerprint.lmdb"
MORGAN_NPY = DATA_ROOT / "07_Step7_分子指纹" / "morgan_fingerprint.npy"
STRUCTURE_LMDB = DATA_ROOT / "08_Step8_结构特征" / "structure_features.lmdb"
HIGH_QUALITY_IDS = DATA_ROOT / "08_Step8_结构特征" / "high_quality_id.txt"

# Model files
CONFIG_YAML = PROJECT_ROOT / "saved_model" / "model" / "run_0" / "complete-full-random-all-0-complex.yml"
CHECKPOINT = PROJECT_ROOT / "saved_model" / "model" / "run_0" / "models" / "best-checkpoint.ckpt"

# Output
OUTPUT_DIR = DATA_ROOT / "09_Step9_模型推理"
OUTPUT_CSV = OUTPUT_DIR / "predictions.csv"

# ============================================================================
# Imports (after sys.path setup)
# ============================================================================
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import torch
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from utils import load_config


def validate_inputs():
    """Check that all required input files exist."""
    required_files = [
        ("Test CSV", TEST_CSV),
        ("Enzyme LMDB", ENZYME_LMDB),
        ("Reaction LMDB", REACTION_LMDB),
        ("GROVER LMDB", GROVER_LMDB),
        ("Morgan NPY", MORGAN_NPY),
        ("Structure LMDB", STRUCTURE_LMDB),
        ("High Quality IDs", HIGH_QUALITY_IDS),
        ("Config YAML", CONFIG_YAML),
        ("Checkpoint", CHECKPOINT),
    ]

    missing = []
    for name, path in required_files:
        if not path.exists():
            missing.append(f"  - {name}: {path}")

    if missing:
        raise FileNotFoundError(
            "Missing required files:\n" + "\n".join(missing)
        )

    print("All input files validated successfully.")


def build_config():
    """Build inference config with Windows-safe overrides."""
    # Load base config from training
    config = load_config(str(CONFIG_YAML))

    # Windows-safe overrides
    config.num_cpus = 0  # num_workers=0 for DataLoader
    config.num_gpus = 1

    # Override data paths for P450 test set
    config.data.tag = "p450_inference"
    config.data.representer = "structure_sequence"

    # All paths point to test CSV (we only use test_dataloader)
    config.data.train_data_path = str(TEST_CSV)
    config.data.val_data_path = str(TEST_CSV)
    config.data.test_data_path = str(TEST_CSV)

    # Feature paths (single dataset, not list)
    config.data.enzyme_lmdb_path = str(ENZYME_LMDB)
    config.data.reaction_lmdb_path = str(REACTION_LMDB)
    config.data.grover_path = str(GROVER_LMDB)
    config.data.morgan_path = str(MORGAN_NPY)
    config.data.structure_processed_path = str(STRUCTURE_LMDB)
    config.data.high_quality_id_path = str(HIGH_QUALITY_IDS)

    # Critical: full_data=False (we don't have str_features.lmdb)
    config.data.full_data = False

    # Placeholder for sequence_processed_path (not used when full_data=False)
    config.data.sequence_processed_path = "str_features.lmdb"

    # Keep other settings from training config
    config.data.batch_size = 16
    config.data.max_substrate_length = 280
    config.data.max_enzyme_length = 1450
    config.data.features = ["morgan", 1024, "grover_mean", 4885]
    config.data.atom_features = ["grover", 2400]
    config.data.fake_sequence_ratio = 0
    config.data.sample_weight = [1.0, 1.0]

    return config


def sigmoid_np(x):
    """Numerically stable sigmoid for numpy arrays."""
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    neg_exp = np.exp(x[~pos])
    out[~pos] = neg_exp / (1.0 + neg_exp)
    return out.astype(np.float32)


@torch.no_grad()
def run_inference(config, device):
    """Run inference using pure forward pass (Windows-safe, no DDP)."""
    # Late imports to ensure sys.path is set
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    print(f"\n{'='*60}")
    print("Preparing data...")
    dm = Singledataset(config)
    test_loader = dm.test_dataloader()  # Also sets dm.test_prediction_df

    n_samples = len(dm.test_prediction_df)
    print(f"Test samples: {n_samples}")
    print(f"Batch size: {config.data.batch_size}")
    print(f"Total batches: {(n_samples + config.data.batch_size - 1) // config.data.batch_size}")

    print(f"\n{'='*60}")
    print(f"Loading model from: {CHECKPOINT}")
    model = SS.load_from_checkpoint(
        str(CHECKPOINT),
        config=config,
        map_location="cpu"
    )
    model.eval()
    model.to(device)
    print(f"Model loaded and moved to {device}")

    print(f"\n{'='*60}")
    print("Running inference...")

    logits_chunks = []
    n_processed = 0

    try:
        from tqdm import tqdm
        iterator = tqdm(test_loader, desc="Inference", unit="batch")
    except ImportError:
        iterator = test_loader
        print("(tqdm not available, showing progress manually)")

    for batch_idx, batch in enumerate(iterator):
        batch = batch.to(device)
        logits, _tags = model(batch)  # logits: [B, 1]
        logits = logits.squeeze(-1).detach().float().cpu().numpy().ravel()
        logits_chunks.append(logits)
        n_processed += len(logits)

        if not hasattr(iterator, 'set_description'):  # No tqdm
            if (batch_idx + 1) % 10 == 0:
                print(f"  Processed {n_processed}/{n_samples} samples...")

    print(f"Inference complete: {n_processed} samples processed")

    # Aggregate results
    logits_all = np.concatenate(logits_chunks) if logits_chunks else np.zeros((0,), dtype=np.float32)
    pred_df = dm.test_prediction_df.copy()

    # Safety check
    assert len(pred_df) == len(logits_all), \
        f"Length mismatch: DataFrame={len(pred_df)}, logits={len(logits_all)}"

    return pred_df, logits_all


def main():
    print("="*60)
    print("Step 9: EZSpecificity Model Inference")
    print("="*60)

    # Validate inputs
    print("\n[1/5] Validating input files...")
    validate_inputs()

    # Setup device
    print("\n[2/5] Setting up device...")
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        print(f"  CUDA available: {torch.version.cuda}")
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
        print(f"  GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    else:
        device = torch.device("cpu")
        print("  CUDA not available, using CPU")

    # Build config
    print("\n[3/5] Building configuration...")
    config = build_config()
    print(f"  full_data: {config.data.full_data}")
    print(f"  batch_size: {config.data.batch_size}")
    print(f"  num_workers: {config.num_cpus}")

    # Run inference
    print("\n[4/5] Running inference...")
    pred_df, logits = run_inference(config, device)

    # Save results
    print("\n[5/5] Saving results...")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Add score columns
    pred_df["logit"] = logits.astype(np.float32)
    pred_df["prob"] = sigmoid_np(logits)
    pred_df["score"] = pred_df["logit"]  # Default score is raw logit

    # Ensure required columns are present and ordered
    required_cols = ["Dock Index", "Enzyme Index", "Substrate Index", "Label"]
    extra_cols = ["score", "logit", "prob"]

    # Verify required columns exist
    for col in required_cols:
        if col not in pred_df.columns:
            raise ValueError(f"Missing required column: {col}")

    # Reorder columns
    final_cols = required_cols + extra_cols
    pred_df = pred_df[final_cols]

    # Save
    pred_df.to_csv(OUTPUT_CSV, index=False)
    print(f"  Saved: {OUTPUT_CSV}")
    print(f"  Total predictions: {len(pred_df)}")

    # Summary statistics
    print(f"\n{'='*60}")
    print("Summary Statistics:")
    print(f"  Predictions: {len(pred_df)}")

    if len(logits) > 0:
        print(f"  Logit range: [{logits.min():.4f}, {logits.max():.4f}]")
        print(f"  Prob range: [{pred_df['prob'].min():.4f}, {pred_df['prob'].max():.4f}]")

        # Label distribution in predictions
        if "Label" in pred_df.columns:
            label_counts = pred_df["Label"].value_counts()
            print(f"  Label=1 (positive): {label_counts.get(1, 0)}")
            print(f"  Label=0 (negative): {label_counts.get(0, 0)}")

        # Score distribution
        print(f"\n  Score percentiles:")
        for p in [10, 25, 50, 75, 90]:
            print(f"    P{p}: {np.percentile(logits, p):.4f}")
    else:
        print("  WARNING: No valid predictions generated!")

    print(f"\n{'='*60}")
    print("Inference completed successfully!")
    print(f"Output: {OUTPUT_CSV}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
