"""
PathB Step 2 - Inference for one factorial experiment.

Runs EZSpecificity inference on a single EXP directory and writes
predictions.csv to --output_dir (defaults to --experiment_dir).

Windows-safe: no DDP, no file_descriptor sharing, num_workers=0.

Usage:
    python step9_inference.py \
        --experiment_dir  <path_to_EXP* in data/> \
        --shared_features <path_to_shared/features> \
        --shared_datasets <path_to_shared/datasets> \
        --checkpoint_dir  <path_to_run_0> \
        --output_dir      <path_to_EXP* in results/>  # optional
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

REQUIRED_OUTPUT_COLS = ["Dock Index", "Enzyme Index", "Substrate Index", "Label"]


def _find_project_root(start: Path) -> Path:
    p = start.resolve()
    while p != p.parent:
        if (p / "src").is_dir() and (p / "saved_model").is_dir():
            return p
        p = p.parent
    raise RuntimeError(f"Cannot find project root from {start}")


def _resolve_dataset_csv(shared_datasets: Path) -> Path:
    for candidate in [shared_datasets / "B6_v1" / "data.csv",
                      shared_datasets / "data.csv"]:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(f"data.csv not found under {shared_datasets}")


def _resolve_checkpoint(checkpoint_dir: Path) -> Tuple[Path, Path]:
    ckpt = checkpoint_dir / "models" / "best-checkpoint.ckpt"
    if not ckpt.is_file():
        candidates = sorted(checkpoint_dir.rglob("best-checkpoint.ckpt"))
        if len(candidates) != 1:
            raise FileNotFoundError(f"Cannot resolve best-checkpoint.ckpt under {checkpoint_dir}")
        ckpt = candidates[0]

    yaml_path = checkpoint_dir / "complete-full-random-all-0-complex.yml"
    if not yaml_path.is_file():
        yamls = sorted(checkpoint_dir.glob("*.yml")) + sorted(checkpoint_dir.glob("*.yaml"))
        if not yamls:
            raise FileNotFoundError(f"No config YAML in {checkpoint_dir}")
        yaml_path = yamls[0]

    return yaml_path, ckpt


def _validate_paths(required: List[Tuple[str, Path]]) -> None:
    missing = [f"  - {name}: {path}" for name, path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(missing))


def _sigmoid_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out.astype(np.float32)


def _build_config(
    config_yaml, dataset_csv, enzyme_lmdb, reaction_lmdb,
    grover_lmdb, morgan_npy, structure_lmdb, high_quality_ids,
    experiment_name,
):
    from utils import load_config

    config = load_config(str(config_yaml))

    config.num_cpus = 0
    config.num_gpus = 1

    config.data.tag = f"pathb_{experiment_name}_inference"
    config.data.representer = "structure_sequence"

    config.data.train_data_path = str(dataset_csv)
    config.data.val_data_path = str(dataset_csv)
    config.data.test_data_path = str(dataset_csv)

    config.data.enzyme_lmdb_path = str(enzyme_lmdb)
    config.data.reaction_lmdb_path = str(reaction_lmdb)
    config.data.grover_path = str(grover_lmdb)
    config.data.morgan_path = str(morgan_npy)
    config.data.structure_processed_path = str(structure_lmdb)
    config.data.high_quality_id_path = str(high_quality_ids)

    config.data.full_data = False
    config.data.sequence_processed_path = "str_features.lmdb"

    config.data.batch_size = 16
    config.data.sample_weight = [1.0, 1.0]
    config.data.fake_sequence_ratio = 0
    config.data.max_substrate_length = 280
    config.data.max_enzyme_length = 1450
    config.data.features = ["morgan", 1024, "grover_mean", 4885]
    config.data.atom_features = ["grover", 2400]

    return config


@torch.no_grad()
def _run_inference(config, checkpoint_path: Path, device: torch.device):
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    dm = Singledataset(config)
    test_loader = dm.test_dataloader()

    model = SS.load_from_checkpoint(str(checkpoint_path), config=config, map_location="cpu")
    model.eval()
    model.to(device)

    logits_chunks: List[np.ndarray] = []
    n_processed = 0
    n_total = len(dm.test_prediction_df)

    for batch_idx, batch in enumerate(test_loader):
        batch = batch.to(device)
        logits, _tags = model(batch)
        logits_np = logits.squeeze(-1).detach().float().cpu().numpy().ravel()
        logits_chunks.append(logits_np)
        n_processed += len(logits_np)
        if (batch_idx + 1) % 10 == 0:
            print(f"  Processed {n_processed}/{n_total} samples...")

    logits_all = (
        np.concatenate(logits_chunks).astype(np.float32)
        if logits_chunks
        else np.zeros((0,), dtype=np.float32)
    )
    pred_df = dm.test_prediction_df.copy()
    assert len(pred_df) == len(logits_all), \
        f"Length mismatch: df={len(pred_df)}, logits={len(logits_all)}"

    return pred_df, logits_all


def parse_args():
    p = argparse.ArgumentParser(description="PathB Step 2 inference for one experiment.")
    p.add_argument("--experiment_dir", required=True)
    p.add_argument("--shared_features", required=True)
    p.add_argument("--shared_datasets", required=True)
    p.add_argument("--checkpoint_dir", required=True)
    p.add_argument("--output_dir", default=None,
                   help="Where to write predictions.csv (default: same as experiment_dir)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    project_root = _find_project_root(Path(__file__).resolve().parent)
    sys.path.insert(0, str(project_root / "src"))

    experiment_dir = Path(args.experiment_dir).resolve()
    shared_features = Path(args.shared_features).resolve()
    shared_datasets = Path(args.shared_datasets).resolve()
    checkpoint_dir = Path(args.checkpoint_dir).resolve()
    output_dir = Path(args.output_dir).resolve() if args.output_dir else experiment_dir

    dataset_csv = _resolve_dataset_csv(shared_datasets)
    config_yaml, checkpoint_path = _resolve_checkpoint(checkpoint_dir)

    structure_lmdb = experiment_dir / "structure_features.lmdb"
    high_quality_ids = experiment_dir / "high_quality_id.txt"
    enzyme_lmdb = shared_features / "enzyme_features.lmdb"
    reaction_lmdb = shared_features / "reaction_features.lmdb"
    grover_lmdb = shared_features / "grover_fingerprint.lmdb"
    morgan_npy = shared_features / "morgan_fingerprint.npy"

    _validate_paths([
        ("dataset_csv", dataset_csv),
        ("enzyme_features.lmdb", enzyme_lmdb),
        ("reaction_features.lmdb", reaction_lmdb),
        ("grover_fingerprint.lmdb", grover_lmdb),
        ("morgan_fingerprint.npy", morgan_npy),
        ("structure_features.lmdb", structure_lmdb),
        ("high_quality_id.txt", high_quality_ids),
        ("config_yaml", config_yaml),
        ("best-checkpoint.ckpt", checkpoint_path),
    ])

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    exp_name = experiment_dir.name
    print(f"[{exp_name}] Device: {device}")

    config = _build_config(
        config_yaml, dataset_csv, enzyme_lmdb, reaction_lmdb,
        grover_lmdb, morgan_npy, structure_lmdb, high_quality_ids,
        exp_name,
    )

    print(f"[{exp_name}] Running inference...")
    pred_df, logits = _run_inference(config, checkpoint_path, device)

    pred_df["logit"] = logits.astype(np.float32)
    pred_df["prob"] = _sigmoid_np(logits)
    pred_df["score"] = pred_df["logit"]

    for col in REQUIRED_OUTPUT_COLS:
        assert col in pred_df.columns, f"Missing column: {col}"

    output_dir.mkdir(parents=True, exist_ok=True)
    out_csv = output_dir / "predictions.csv"
    pred_df[REQUIRED_OUTPUT_COLS + ["score", "logit", "prob"]].to_csv(out_csv, index=False)

    print(f"[{exp_name}] Saved: {out_csv}")
    print(f"[{exp_name}] Samples: {len(pred_df)}, "
          f"prob range=({pred_df['prob'].min():.4f}, {pred_df['prob'].max():.4f})")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
