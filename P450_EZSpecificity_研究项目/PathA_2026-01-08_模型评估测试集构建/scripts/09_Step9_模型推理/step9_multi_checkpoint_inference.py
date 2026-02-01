"""
Multi-Checkpoint Inference Script
==================================
This script runs inference using different model checkpoints.
Supports: best-checkpoint.ckpt, best-checkpoint-v1.ckpt, v2, v3, v4

Usage:
    python step9_multi_checkpoint_inference.py --checkpoint best-checkpoint-v1

    Available checkpoints:
    - best-checkpoint (default, same as best-checkpoint.ckpt)
    - best-checkpoint-v1
    - best-checkpoint-v2
    - best-checkpoint-v3
    - best-checkpoint-v4
"""

import sys
import os
import argparse
from pathlib import Path
from datetime import datetime

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
# Available Checkpoints
# ============================================================================
CHECKPOINT_DIR = PROJECT_ROOT / "saved_model" / "model" / "run_0" / "models"
AVAILABLE_CHECKPOINTS = {
    "best-checkpoint": "best-checkpoint.ckpt",
    "best-checkpoint-v1": "best-checkpoint-v1.ckpt",
    "best-checkpoint-v2": "best-checkpoint-v2.ckpt",
    "best-checkpoint-v3": "best-checkpoint-v3.ckpt",
    "best-checkpoint-v4": "best-checkpoint-v4.ckpt",
}

# ============================================================================
# Feature File Paths
# ============================================================================
DATA_ROOT = PATHA_ROOT / "data"

# Input files (same for all runs)
TEST_CSV = DATA_ROOT / "04_Step4_格式修正后数据" / "data.csv"
ENZYME_LMDB = DATA_ROOT / "05_Step5_ESM酶嵌入" / "enzyme_features.lmdb"
REACTION_LMDB = DATA_ROOT / "06_Step6_反应图特征" / "reaction_features.lmdb"
GROVER_LMDB = DATA_ROOT / "07_Step7_分子指纹" / "grover_fingerprint.lmdb"
MORGAN_NPY = DATA_ROOT / "07_Step7_分子指纹" / "morgan_fingerprint.npy"
STRUCTURE_LMDB = DATA_ROOT / "08_Step8_结构特征" / "structure_features.lmdb"
HIGH_QUALITY_IDS = DATA_ROOT / "08_Step8_结构特征" / "high_quality_id.txt"
CONFIG_YAML = PROJECT_ROOT / "saved_model" / "model" / "run_0" / "complete-full-random-all-0-complex.yml"

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


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Run inference with different checkpoints")
    parser.add_argument(
        "--checkpoint", "-c",
        type=str,
        default="best-checkpoint",
        choices=list(AVAILABLE_CHECKPOINTS.keys()),
        help="Checkpoint to use for inference"
    )
    parser.add_argument(
        "--run-id", "-r",
        type=str,
        default=None,
        help="Custom run ID (default: auto-generated based on checkpoint and timestamp)"
    )
    return parser.parse_args()


def validate_inputs(checkpoint_path):
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
        ("Checkpoint", checkpoint_path),
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
    config = load_config(str(CONFIG_YAML))

    # Windows-safe overrides
    config.num_cpus = 0
    config.num_gpus = 1

    # Override data paths
    config.data.tag = "p450_inference"
    config.data.representer = "structure_sequence"
    config.data.train_data_path = str(TEST_CSV)
    config.data.val_data_path = str(TEST_CSV)
    config.data.test_data_path = str(TEST_CSV)
    config.data.enzyme_lmdb_path = str(ENZYME_LMDB)
    config.data.reaction_lmdb_path = str(REACTION_LMDB)
    config.data.grover_path = str(GROVER_LMDB)
    config.data.morgan_path = str(MORGAN_NPY)
    config.data.structure_processed_path = str(STRUCTURE_LMDB)
    config.data.high_quality_id_path = str(HIGH_QUALITY_IDS)
    config.data.full_data = False
    config.data.sequence_processed_path = "str_features.lmdb"
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
def run_inference(config, checkpoint_path, device):
    """Run inference using pure forward pass."""
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    print(f"\n{'='*60}")
    print("Preparing data...")
    dm = Singledataset(config)
    test_loader = dm.test_dataloader()

    n_samples = len(dm.test_prediction_df)
    print(f"Test samples: {n_samples}")

    print(f"\n{'='*60}")
    print(f"Loading model from: {checkpoint_path}")
    model = SS.load_from_checkpoint(
        str(checkpoint_path),
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

    for batch_idx, batch in enumerate(iterator):
        batch = batch.to(device)
        logits, _tags = model(batch)
        logits = logits.squeeze(-1).detach().float().cpu().numpy().ravel()
        logits_chunks.append(logits)
        n_processed += len(logits)

    print(f"Inference complete: {n_processed} samples processed")

    logits_all = np.concatenate(logits_chunks) if logits_chunks else np.zeros((0,), dtype=np.float32)
    pred_df = dm.test_prediction_df.copy()

    assert len(pred_df) == len(logits_all), \
        f"Length mismatch: DataFrame={len(pred_df)}, logits={len(logits_all)}"

    return pred_df, logits_all


def main():
    args = parse_args()

    # Determine checkpoint
    checkpoint_name = args.checkpoint
    checkpoint_file = AVAILABLE_CHECKPOINTS[checkpoint_name]
    checkpoint_path = CHECKPOINT_DIR / checkpoint_file

    # Generate run ID
    if args.run_id:
        run_id = args.run_id
    else:
        # Find next available run number
        output_base = DATA_ROOT / "09_Step9_模型推理"
        existing_runs = [d.name for d in output_base.iterdir() if d.is_dir() and d.name.startswith("Run")]
        run_numbers = []
        for r in existing_runs:
            try:
                num = int(r.split("_")[0].replace("Run", ""))
                run_numbers.append(num)
            except:
                pass
        next_num = max(run_numbers, default=0) + 1
        run_id = f"Run{next_num:02d}_{checkpoint_name}"

    print("="*60)
    print("Multi-Checkpoint Inference")
    print("="*60)
    print(f"\nCheckpoint: {checkpoint_name}")
    print(f"Checkpoint file: {checkpoint_file}")
    print(f"Run ID: {run_id}")

    # Setup output directories
    OUTPUT_DIR_STEP9 = DATA_ROOT / "09_Step9_模型推理" / run_id
    OUTPUT_CSV = OUTPUT_DIR_STEP9 / "predictions.csv"

    # Validate inputs
    print("\n[1/5] Validating input files...")
    validate_inputs(checkpoint_path)

    # Setup device
    print("\n[2/5] Setting up device...")
    if torch.cuda.is_available():
        device = torch.device("cuda:0")
        print(f"  GPU: {torch.cuda.get_device_name(0)}")
    else:
        device = torch.device("cpu")
        print("  Using CPU")

    # Build config
    print("\n[3/5] Building configuration...")
    config = build_config()

    # Run inference
    print("\n[4/5] Running inference...")
    pred_df, logits = run_inference(config, checkpoint_path, device)

    # Save results
    print("\n[5/5] Saving results...")
    OUTPUT_DIR_STEP9.mkdir(parents=True, exist_ok=True)

    pred_df["logit"] = logits.astype(np.float32)
    pred_df["prob"] = sigmoid_np(logits)
    pred_df["score"] = pred_df["logit"]

    required_cols = ["Dock Index", "Enzyme Index", "Substrate Index", "Label"]
    extra_cols = ["score", "logit", "prob"]
    final_cols = required_cols + extra_cols
    pred_df = pred_df[final_cols]

    pred_df.to_csv(OUTPUT_CSV, index=False)
    print(f"  Saved: {OUTPUT_CSV}")

    # Create README
    readme_content = f"""# {run_id}

## 基本信息

| 项目 | 值 |
|------|-----|
| **运行编号** | {run_id} |
| **检查点文件** | `{checkpoint_file}` |
| **检查点路径** | `saved_model/model/run_0/models/{checkpoint_file}` |
| **运行日期** | {datetime.now().strftime('%Y-%m-%d %H:%M')} |

## 输出文件

- `predictions.csv`: {len(pred_df)}条预测结果

## 统计摘要

| 指标 | 值 |
|------|-----|
| 预测数量 | {len(pred_df)} |
| Logit范围 | [{logits.min():.4f}, {logits.max():.4f}] |
| Prob范围 | [{pred_df['prob'].min():.4f}, {pred_df['prob'].max():.4f}] |
| Logit均值 | {logits.mean():.4f} |
| Logit标准差 | {logits.std():.4f} |
"""

    readme_path = OUTPUT_DIR_STEP9 / "README.md"
    with open(readme_path, "w", encoding="utf-8") as f:
        f.write(readme_content)
    print(f"  Saved: {readme_path}")

    # Summary
    print(f"\n{'='*60}")
    print("Summary Statistics:")
    print(f"  Predictions: {len(pred_df)}")
    print(f"  Logit range: [{logits.min():.4f}, {logits.max():.4f}]")
    print(f"  Prob range: [{pred_df['prob'].min():.4f}, {pred_df['prob'].max():.4f}]")
    print(f"\n{'='*60}")
    print(f"Output: {OUTPUT_DIR_STEP9}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
