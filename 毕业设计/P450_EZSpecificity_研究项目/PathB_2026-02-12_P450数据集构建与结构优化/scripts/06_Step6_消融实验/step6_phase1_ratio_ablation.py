"""
PathB Step 6 — Phase 1: Ratio Ablation Experiments.

Isolates the effect of positive:negative class ratio on AUC-ROC by
subsampling Step 5's Vina-docked dataset at 1:1, 1:3, and 1:9 ratios,
while keeping structure source (Vina) and negative type (random) constant.

Experiments:
  ABL-01: 265 pos + 265 neg  (1:1)
  ABL-02: 265 pos + 795 neg  (1:3)
  ABL-03: 265 pos + 2501 neg (1:9.44 — all negatives, reproduces Step 5)

Design choices (per Codex review):
  - Nested subsets: ABL-01 negatives ⊂ ABL-02 negatives ⊂ ABL-03
  - Multi-seed repeats: 5 seeds for ABL-01/02; ABL-03 uses all negatives
  - Bootstrap CI: 2000 iterations for AUC-ROC confidence intervals
  - Reuse: Step 5's structure_features.lmdb + high_quality_id.txt (no regeneration)

Usage:
    python step6_phase1_ratio_ablation.py \\
        --step5_data      <results/04_Step4_批量对接/data.csv> \\
        --step5_lmdb      <data/05_Step5_重构评估/structure_features.lmdb> \\
        --step5_hqid      <data/05_Step5_重构评估/high_quality_id.txt> \\
        --shared_features <data/00_shared/features> \\
        --checkpoint_dir  <saved_model/model/run_0> \\
        --data_dir        <data/06_Step6_消融实验> \\
        --results_dir     <results/06_Step6_消融实验> \\
        [--exp01_predictions <results/02_Step2_因子实验/EXP01_.../predictions.csv>] \\
        [--step5_predictions <results/05_Step5_重构评估/predictions.csv>] \\
        [--seeds 42,123,456,789,2026]
"""

from __future__ import annotations

import argparse
import csv
import sys
import warnings
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

# ─── Constants ──────────────────────────────────────────────────────────
DEFAULT_SEEDS = [42, 123, 456, 789, 2026]
BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 42

PHASE1_EXPERIMENTS = {
    "ABL-01": {"neg_ratio": 1, "label": "Vina 1:1 (random neg)"},
    "ABL-02": {"neg_ratio": 3, "label": "Vina 1:3 (random neg)"},
    "ABL-03": {"neg_ratio": None, "label": "Vina 1:9.44 (all neg, ≈Step 5)"},
}

REQUIRED_COLS = ["Dock Index", "Enzyme Index", "Substrate Index", "Label"]


# ═══════════════════════════════════════════════════════════════════════
#  1. DATA SUBSAMPLING
# ═══════════════════════════════════════════════════════════════════════

def subsample_nested(
    df: pd.DataFrame,
    ratios: List[int],
    seed: int,
) -> Dict[int, pd.DataFrame]:
    """Subsample negatives at multiple ratios using nested subsets.

    For a given seed, negatives are shuffled once. Then:
      ratio=1 → first N negatives
      ratio=3 → first 3N negatives (superset of ratio=1)
      ratio=9 → first 9N or all negatives

    Returns dict mapping ratio → subsampled DataFrame.
    """
    pos = df[df["Label"] == 1].copy()
    neg = df[df["Label"] == 0].copy()
    n_pos = len(pos)

    rng = np.random.RandomState(seed)
    neg_shuffled = neg.sample(frac=1, random_state=rng)

    result = {}
    for ratio in sorted(ratios):
        n_neg = min(n_pos * ratio, len(neg_shuffled))
        neg_subset = neg_shuffled.head(n_neg)
        combined = pd.concat([pos, neg_subset], ignore_index=True)
        combined = combined.sort_values("Dock Index").reset_index(drop=True)
        result[ratio] = combined

    return result


# ═══════════════════════════════════════════════════════════════════════
#  2. INFERENCE
# ═══════════════════════════════════════════════════════════════════════

def find_project_root(start: Path) -> Path:
    p = start.resolve()
    while p != p.parent:
        if (p / "src").is_dir() and (p / "saved_model").is_dir():
            return p
        p = p.parent
    raise RuntimeError(f"Cannot find project root from {start}")


def resolve_checkpoint(checkpoint_dir: Path) -> Tuple[Path, Path]:
    ckpt = checkpoint_dir / "models" / "best-checkpoint.ckpt"
    if not ckpt.is_file():
        candidates = sorted(checkpoint_dir.rglob("best-checkpoint.ckpt"))
        if len(candidates) != 1:
            raise FileNotFoundError(f"Cannot resolve checkpoint under {checkpoint_dir}")
        ckpt = candidates[0]

    yaml_path = checkpoint_dir / "complete-full-random-all-0-complex.yml"
    if not yaml_path.is_file():
        yamls = sorted(checkpoint_dir.glob("*.yml")) + sorted(checkpoint_dir.glob("*.yaml"))
        if not yamls:
            raise FileNotFoundError(f"No config YAML in {checkpoint_dir}")
        yaml_path = yamls[0]

    return yaml_path, ckpt


def sigmoid_np(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out.astype(np.float32)


def build_config(
    config_yaml: Path,
    dataset_csv: Path,
    enzyme_lmdb: Path,
    reaction_lmdb: Path,
    grover_lmdb: Path,
    morgan_npy: Path,
    structure_lmdb: Path,
    high_quality_ids: Path,
    experiment_name: str,
):
    from utils import load_config

    config = load_config(str(config_yaml))
    config.num_cpus = 0
    config.num_gpus = 1

    config.data.tag = f"ablation_{experiment_name}"
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
def run_inference(
    config,
    checkpoint_path: Path,
    device: torch.device,
) -> Tuple[pd.DataFrame, np.ndarray]:
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

    logits_all = (
        np.concatenate(logits_chunks).astype(np.float32)
        if logits_chunks
        else np.zeros((0,), dtype=np.float32)
    )
    pred_df = dm.test_prediction_df.copy()
    assert len(pred_df) == len(logits_all), \
        f"Length mismatch: df={len(pred_df)}, logits={len(logits_all)}"

    return pred_df, logits_all


def run_single_experiment(
    dataset_csv: Path,
    structure_lmdb: Path,
    high_quality_ids: Path,
    shared_features: Path,
    checkpoint_dir: Path,
    experiment_name: str,
    device: torch.device,
) -> pd.DataFrame:
    """Run inference on a single data.csv and return predictions DataFrame."""
    config_yaml, checkpoint_path = resolve_checkpoint(checkpoint_dir)

    config = build_config(
        config_yaml=config_yaml,
        dataset_csv=dataset_csv,
        enzyme_lmdb=shared_features / "enzyme_features.lmdb",
        reaction_lmdb=shared_features / "reaction_features.lmdb",
        grover_lmdb=shared_features / "grover_fingerprint.lmdb",
        morgan_npy=shared_features / "morgan_fingerprint.npy",
        structure_lmdb=structure_lmdb,
        high_quality_ids=high_quality_ids,
        experiment_name=experiment_name,
    )

    pred_df, logits = run_inference(config, checkpoint_path, device)
    pred_df["score"] = logits.astype(np.float32)
    pred_df["logit"] = logits.astype(np.float32)
    pred_df["prob"] = sigmoid_np(logits)

    return pred_df


# ═══════════════════════════════════════════════════════════════════════
#  3. ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def bootstrap_auc(labels, scores, n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED):
    from sklearn.metrics import roc_auc_score
    n = len(labels)
    if n == 0:
        return float("nan"), float("nan")
    rng = np.random.RandomState(seed)
    vals = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)
        bl, bs = labels[idx], scores[idx]
        if bl.sum() == 0 or bl.sum() == n:
            vals[i] = np.nan
        else:
            vals[i] = roc_auc_score(bl, bs)
    vals = vals[~np.isnan(vals)]
    if vals.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def compute_metrics(labels, scores, name=""):
    from sklearn.metrics import roc_auc_score, average_precision_score
    labels = np.asarray(labels).astype(int)
    scores = np.asarray(scores, dtype=float)

    n_total = len(labels)
    n_pos = int(labels.sum())
    n_neg = n_total - n_pos

    if n_total == 0 or n_pos == 0 or n_neg == 0:
        auc_roc = float("nan")
        auc_pr = float("nan")
        ci_lo, ci_hi = float("nan"), float("nan")
    else:
        auc_roc = roc_auc_score(labels, scores)
        auc_pr = average_precision_score(labels, scores)
        ci_lo, ci_hi = bootstrap_auc(labels, scores)

    score_pos = float(scores[labels == 1].mean()) if n_pos > 0 else float("nan")
    score_neg = float(scores[labels == 0].mean()) if n_neg > 0 else float("nan")
    return {
        "name": name,
        "n_total": n_total,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "ratio": f"1:{n_neg / max(n_pos, 1):.1f}",
        "auc_roc": round(auc_roc, 4),
        "auc_roc_ci_lo": round(ci_lo, 4),
        "auc_roc_ci_hi": round(ci_hi, 4),
        "auc_pr": round(auc_pr, 4),
        "score_mean_pos": round(float(score_pos), 4),
        "score_mean_neg": round(float(score_neg), 4),
        "score_separation": round(float(score_pos - score_neg), 4),
    }


def load_predictions(csv_path: Path):
    df = pd.read_csv(csv_path)
    missing = {"Label", "score"} - set(df.columns)
    if missing:
        raise ValueError(f"{csv_path} missing required columns: {sorted(missing)}")
    labels = df["Label"].values.astype(int)
    scores = df["score"].values.astype(float)
    return labels, scores


def parse_seeds(seed_arg: str) -> List[int]:
    seeds: List[int] = []
    for token in seed_arg.split(","):
        token = token.strip()
        if not token:
            continue
        seeds.append(int(token))
    if not seeds:
        raise ValueError("--seeds produced an empty list")
    return seeds


def validate_reused_predictions(pred_csv: Path, expected_df: pd.DataFrame) -> None:
    """Validate that reused predictions match the expected dataset."""
    pred_df = pd.read_csv(pred_csv, usecols=REQUIRED_COLS)
    expected = expected_df[REQUIRED_COLS]

    if len(pred_df) != len(expected):
        raise ValueError(
            f"Reused predictions row count mismatch: got {len(pred_df)}, "
            f"expected {len(expected)}"
        )
    if set(pred_df["Dock Index"].astype(int)) != set(expected["Dock Index"].astype(int)):
        raise ValueError(
            "Reused predictions Dock Index set does not match expected subset"
        )


# ═══════════════════════════════════════════════════════════════════════
#  4. VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def plot_phase1_roc(
    all_metrics: List[dict],
    all_labels_scores: Dict[str, Tuple[np.ndarray, np.ndarray]],
    output_path: Path,
    reference_datasets: Optional[Dict[str, Tuple[np.ndarray, np.ndarray]]] = None,
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))
    colors = {"ABL-01": "#2196F3", "ABL-02": "#FF9800", "ABL-03": "#4CAF50"}

    for name, (labels, scores) in sorted(all_labels_scores.items()):
        if len(np.unique(labels)) < 2:
            print(f"  [WARN] Skipping ROC for {name}: single-class labels")
            continue
        fpr, tpr, _ = roc_curve(labels, scores)
        auc = roc_auc_score(labels, scores)
        color = colors.get(name, "gray")
        ax.plot(fpr, tpr, color=color, linewidth=2,
                label=f"{name} (AUC={auc:.4f})")

    if reference_datasets:
        ref_styles = {"EXP01": ("green", "--"), "Step 5": ("red", ":")}
        for name, (labels, scores) in reference_datasets.items():
            if len(np.unique(labels)) < 2:
                continue
            fpr, tpr, _ = roc_curve(labels, scores)
            auc = roc_auc_score(labels, scores)
            color, ls = ref_styles.get(name, ("gray", "-."))
            ax.plot(fpr, tpr, color=color, linestyle=ls, linewidth=1.5,
                    label=f"{name} (AUC={auc:.4f})")

    ax.plot([0, 1], [0, 1], "k--", linewidth=1, alpha=0.4, label="Random (0.5)")
    ax.set_xlabel("False Positive Rate", fontsize=12)
    ax.set_ylabel("True Positive Rate", fontsize=12)
    ax.set_title("Phase 1: Ratio Ablation — ROC Curves", fontsize=14)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  ROC plot saved: {output_path}")


def plot_auc_vs_ratio(all_metrics: List[dict], output_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ratios = []
    aucs = []
    ci_los = []
    ci_his = []
    labels = []

    for m in all_metrics:
        if m["name"].startswith("ABL-"):
            neg_count = m["n_neg"]
            pos_count = m["n_pos"]
            ratio_val = neg_count / max(pos_count, 1)
            ratios.append(ratio_val)
            aucs.append(m["auc_roc"])
            ci_los.append(m["auc_roc_ci_lo"])
            ci_his.append(m["auc_roc_ci_hi"])
            labels.append(m["name"])

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    yerr_lo = [a - lo for a, lo in zip(aucs, ci_los)]
    yerr_hi = [hi - a for a, hi in zip(aucs, ci_his)]
    ax.errorbar(ratios, aucs, yerr=[yerr_lo, yerr_hi],
                fmt="o-", color="#2196F3", capsize=5, markersize=8, linewidth=2)
    for r, a, lab in zip(ratios, aucs, labels):
        ax.annotate(f"{lab}\n{a:.4f}", (r, a), textcoords="offset points",
                    xytext=(10, 10), fontsize=9)

    ax.axhline(0.5, color="red", linestyle="--", alpha=0.5, label="Random (0.5)")
    ax.set_xlabel("Negative:Positive Ratio", fontsize=12)
    ax.set_ylabel("AUC-ROC", fontsize=12)
    ax.set_title("Phase 1: AUC-ROC vs Class Ratio", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"  AUC vs ratio plot saved: {output_path}")


def plot_score_distributions(
    all_labels_scores: Dict[str, Tuple[np.ndarray, np.ndarray]],
    output_path: Path,
):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(all_labels_scores)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, (name, (labels, scores)) in zip(axes, sorted(all_labels_scores.items())):
        pos_scores = scores[labels == 1]
        neg_scores = scores[labels == 0]
        lo = min(scores.min(), -15)
        hi = max(scores.max(), 5)
        bins = np.linspace(lo, hi, 50)
        ax.hist(pos_scores, bins=bins, alpha=0.6, color="blue", density=True,
                label=f"Pos (n={len(pos_scores)})")
        ax.hist(neg_scores, bins=bins, alpha=0.6, color="red", density=True,
                label=f"Neg (n={len(neg_scores)})")
        ax.set_xlabel("Score (logit)")
        ax.set_title(name, fontsize=11)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    fig.suptitle("Phase 1: Score Distributions", fontsize=13, y=1.02)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Score distribution plot saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
#  5. REPORT GENERATION
# ═══════════════════════════════════════════════════════════════════════

def write_phase1_report(
    all_metrics: List[dict],
    seed_metrics: Dict[str, List[dict]],
    output_path: Path,
):
    lines = [
        "# Phase 1: Ratio Ablation Report",
        "",
        f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Branch**: pathb-ablation",
        "",
        "## Experimental Design",
        "",
        "Isolate the effect of positive:negative ratio on AUC-ROC,",
        "holding structure source (Vina) and negative type (random) constant.",
        "",
        "| Experiment | Ratio | Pos | Neg | Structure | Neg Type |",
        "|-----------|-------|-----|-----|-----------|----------|",
    ]
    for m in all_metrics:
        if m["name"].startswith("ABL-"):
            lines.append(
                f"| {m['name']} | {m['ratio']} | {m['n_pos']} | {m['n_neg']} "
                f"| Vina | Random |"
            )
    lines.append("")

    lines.append("## Results")
    lines.append("")
    lines.append("| Experiment | AUC-ROC | 95% CI | AUC-PR | Score Sep |")
    lines.append("|-----------|---------|--------|--------|----------|")
    for m in all_metrics:
        lines.append(
            f"| {m['name']} | **{m['auc_roc']:.4f}** | "
            f"[{m['auc_roc_ci_lo']:.4f}, {m['auc_roc_ci_hi']:.4f}] | "
            f"{m['auc_pr']:.4f} | {m['score_separation']:.4f} |"
        )
    lines.append("")

    # Multi-seed statistics
    for exp_name, seeds_m in sorted(seed_metrics.items()):
        if len(seeds_m) > 1:
            aucs = [m["auc_roc"] for m in seeds_m]
            lines.append(f"### {exp_name} Multi-Seed Results ({len(seeds_m)} seeds)")
            lines.append("")
            lines.append(f"- Mean AUC-ROC: **{np.mean(aucs):.4f}** ± {np.std(aucs):.4f}")
            lines.append(f"- Range: [{min(aucs):.4f}, {max(aucs):.4f}]")
            lines.append(f"- Seeds: {[m.get('seed', '?') for m in seeds_m]}")
            lines.append(f"- AUCs: {[round(a, 4) for a in aucs]}")
            lines.append("")

    # Interpretation
    abl01_m = next((m for m in all_metrics if m["name"] == "ABL-01"), None)
    abl03_m = next((m for m in all_metrics if m["name"] == "ABL-03"), None)
    exp01_m = next((m for m in all_metrics if m["name"] == "EXP01"), None)

    lines.append("## Interpretation")
    lines.append("")
    if abl01_m and abl03_m:
        delta = abl01_m["auc_roc"] - abl03_m["auc_roc"]
        lines.append(f"**Ratio Effect**: ABL-01 (1:1) - ABL-03 (1:9) = "
                      f"{delta:+.4f}")
        lines.append("")
        if abs(delta) < 0.05:
            lines.append("Ratio has **minimal effect** on AUC-ROC. "
                          "The performance drop from EXP01 to Step 5 is NOT primarily "
                          "due to class imbalance.")
            lines.append("")
            lines.append("**Next step**: Phase 2 — isolate structure source vs "
                          "negative type effect.")
        elif delta > 0.10:
            lines.append("Ratio has a **substantial positive effect** on AUC-ROC. "
                          "Balanced classes improve model discrimination.")
            lines.append("")
            if exp01_m:
                gap_to_exp01 = exp01_m["auc_roc"] - abl01_m["auc_roc"]
                lines.append(f"**Remaining gap** to EXP01: {gap_to_exp01:.4f}")
                if gap_to_exp01 > 0.05:
                    lines.append("Ratio alone does NOT fully explain the EXP01-Step5 gap. "
                                  "Structure source and/or negative type also contribute.")
                    lines.append("")
                    lines.append("**Next step**: Phase 2 — further ablation needed.")
                else:
                    lines.append("Ratio appears to be the **primary factor**. "
                                  "The gap to EXP01 is small.")
        else:
            lines.append(f"Ratio has a **moderate effect** ({delta:+.4f}). "
                          "It is a contributing factor but not the sole cause.")
            lines.append("")
            lines.append("**Next step**: Phase 2 — isolate structure source effect.")

    lines.append("")
    lines.append("## Files Generated")
    lines.append("")
    lines.append("- `phase1_metrics.csv` — All experiment metrics")
    lines.append("- `phase1_roc.png` — ROC curve comparison")
    lines.append("- `phase1_auc_vs_ratio.png` — AUC vs ratio plot")
    lines.append("- `phase1_score_dist.png` — Score distributions")
    lines.append("- `ABL-*/predictions*.csv` — Per-experiment predictions")
    lines.append("")

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"  Report saved: {output_path}")


# ═══════════════════════════════════════════════════════════════════════
#  6. MAIN
# ═══════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description="Step 6 Phase 1: Ratio Ablation")
    p.add_argument("--step5_data", required=True,
                   help="Path to Step 4/5 data.csv (2766 rows)")
    p.add_argument("--step5_lmdb", required=True,
                   help="Path to Step 5 structure_features.lmdb")
    p.add_argument("--step5_hqid", required=True,
                   help="Path to Step 5 high_quality_id.txt")
    p.add_argument("--shared_features", required=True,
                   help="Path to data/00_shared/features/")
    p.add_argument("--checkpoint_dir", required=True,
                   help="Path to saved_model/model/run_0/")
    p.add_argument("--data_dir", required=True,
                   help="Path to data/06_Step6_消融实验/")
    p.add_argument("--results_dir", required=True,
                   help="Path to results/06_Step6_消融实验/")
    p.add_argument("--exp01_predictions", default="",
                   help="EXP01 predictions.csv for reference comparison")
    p.add_argument("--step5_predictions", default="",
                   help="Step 5 predictions.csv (if available, skip ABL-03 inference)")
    p.add_argument("--seeds", default=",".join(str(s) for s in DEFAULT_SEEDS),
                   help="Comma-separated random seeds (default: 42,123,456,789,2026)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    project_root = find_project_root(Path(__file__).resolve().parent)
    sys.path.insert(0, str(project_root / "src"))

    step5_data = Path(args.step5_data).resolve()
    step5_lmdb = Path(args.step5_lmdb).resolve()
    step5_hqid = Path(args.step5_hqid).resolve()
    shared_features = Path(args.shared_features).resolve()
    checkpoint_dir = Path(args.checkpoint_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    results_dir = Path(args.results_dir).resolve()
    seeds = parse_seeds(args.seeds)

    # Validate inputs
    for name, path in [
        ("step5_data", step5_data),
        ("step5_lmdb", step5_lmdb),
        ("step5_hqid", step5_hqid),
        ("shared_features", shared_features),
        ("checkpoint_dir", checkpoint_dir),
    ]:
        if not path.exists():
            print(f"[ERROR] {name} not found: {path}")
            return 1

    # Load Step 5 data
    print(f"Loading Step 5 data: {step5_data}")
    df_full = pd.read_csv(step5_data)
    n_pos = (df_full["Label"] == 1).sum()
    n_neg = (df_full["Label"] == 0).sum()
    print(f"  Total: {len(df_full)} rows ({n_pos} pos, {n_neg} neg, ratio 1:{n_neg/n_pos:.1f})")

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Precompute nested subsets for all seeds (Codex recommendation:
    # one shuffle per seed, then slice at each ratio → guarantees nesting)
    sampled_ratios = sorted(
        cfg["neg_ratio"]
        for cfg in PHASE1_EXPERIMENTS.values()
        if cfg["neg_ratio"] is not None
    )
    print(f"Precomputing nested subsets for ratios {sampled_ratios} x {len(seeds)} seeds...")
    subsets_by_seed = {
        seed: subsample_nested(df_full, sampled_ratios, seed) for seed in seeds
    }

    # ── Phase 1 experiments ──────────────────────────────────────────
    print("\n" + "=" * 60)
    print("PHASE 1: RATIO ABLATION")
    print("=" * 60)

    all_metrics: List[dict] = []
    all_labels_scores: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    seed_metrics: Dict[str, List[dict]] = defaultdict(list)

    for exp_name, exp_config in PHASE1_EXPERIMENTS.items():
        neg_ratio = exp_config["neg_ratio"]
        exp_label = exp_config["label"]

        exp_data_dir = data_dir / f"{exp_name}_vina_1to{neg_ratio or 'all'}"
        exp_results_dir = results_dir / f"{exp_name}_vina_1to{neg_ratio or 'all'}"
        exp_data_dir.mkdir(parents=True, exist_ok=True)
        exp_results_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n── {exp_name}: {exp_label} ──")

        if neg_ratio is None:
            run_seeds = [seeds[0]]  # ABL-03: all negatives, no sampling randomness
        else:
            run_seeds = seeds

        exp_aucs = []
        exp_seed_payloads: List[Tuple[np.ndarray, np.ndarray]] = []

        for seed in run_seeds:
            seed_suffix = f"_seed{seed}" if len(run_seeds) > 1 else ""

            # Subsample (from precomputed nested cache)
            if neg_ratio is not None:
                df_sub = subsets_by_seed[seed][neg_ratio]
            else:
                df_sub = df_full.copy()

            # Save subsampled data.csv
            data_csv_path = exp_data_dir / f"data{seed_suffix}.csv"
            df_sub[REQUIRED_COLS].to_csv(data_csv_path, index=False)
            n_p = (df_sub["Label"] == 1).sum()
            n_n = (df_sub["Label"] == 0).sum()
            print(f"  Seed {seed}: {len(df_sub)} rows ({n_p} pos, {n_n} neg)")

            # Check if Step 5 predictions can be reused for ABL-03
            pred_csv_path = exp_results_dir / f"predictions{seed_suffix}.csv"
            if neg_ratio is None and args.step5_predictions and Path(args.step5_predictions).exists():
                print(f"  Reusing Step 5 predictions: {args.step5_predictions}")
                import shutil
                shutil.copy2(args.step5_predictions, pred_csv_path)
                validate_reused_predictions(pred_csv_path, df_sub)
            else:
                # Run inference
                print(f"  Running inference...")
                pred_df = run_single_experiment(
                    dataset_csv=data_csv_path,
                    structure_lmdb=step5_lmdb,
                    high_quality_ids=step5_hqid,
                    shared_features=shared_features,
                    checkpoint_dir=checkpoint_dir,
                    experiment_name=f"{exp_name}_s{seed}",
                    device=device,
                )
                pred_df[REQUIRED_COLS + ["score", "logit", "prob"]].to_csv(
                    pred_csv_path, index=False
                )

            # Compute metrics (on post-filter predictions — may differ from pre-filter counts)
            labels, scores = load_predictions(pred_csv_path)
            actual_pos = int((labels == 1).sum())
            actual_neg = int((labels == 0).sum())
            if actual_pos != int(n_p) or actual_neg != int(n_n):
                print(
                    f"  [WARN] Post-filter count drift (seed={seed}): "
                    f"expected ({int(n_p)} pos, {int(n_n)} neg) -> "
                    f"actual ({actual_pos} pos, {actual_neg} neg)"
                )
            m = compute_metrics(labels, scores, f"{exp_name} (seed={seed})")
            m["seed"] = seed
            seed_metrics[exp_name].append(m)
            exp_aucs.append(m["auc_roc"])
            exp_seed_payloads.append((labels, scores))
            print(f"  AUC-ROC = {m['auc_roc']:.4f} [{m['auc_roc_ci_lo']:.4f}, {m['auc_roc_ci_hi']:.4f}]")

        # Aggregate multi-seed results
        mean_auc = np.mean(exp_aucs)
        std_auc = np.std(exp_aucs) if len(exp_aucs) > 1 else 0.0
        print(f"  {'Mean ' if len(exp_aucs) > 1 else ''}AUC-ROC: {mean_auc:.4f}"
              f"{f' ± {std_auc:.4f}' if std_auc > 0 else ''}")

        # Use median-seed as representative (consistent for both metrics and plots)
        valid_aucs = [(i, a) for i, a in enumerate(exp_aucs) if not np.isnan(a)]
        if valid_aucs:
            med = np.median([a for _, a in valid_aucs])
            median_idx = min(valid_aucs, key=lambda x: abs(x[1] - med))[0]
        else:
            median_idx = 0
        representative_m = seed_metrics[exp_name][median_idx].copy()
        representative_m["name"] = exp_name
        all_metrics.append(representative_m)
        rep_labels, rep_scores = exp_seed_payloads[median_idx]
        all_labels_scores[exp_name] = (rep_labels, rep_scores)

        # Save per-experiment metrics
        metrics_csv = exp_results_dir / "metrics.csv"
        with open(metrics_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(seed_metrics[exp_name][0].keys()))
            writer.writeheader()
            writer.writerows(seed_metrics[exp_name])
        print(f"  Metrics saved: {metrics_csv}")

    # ── Reference datasets ───────────────────────────────────────────
    ref_datasets: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}

    if args.exp01_predictions and Path(args.exp01_predictions).exists():
        e1_labels, e1_scores = load_predictions(Path(args.exp01_predictions))
        e1_m = compute_metrics(e1_labels, e1_scores, "EXP01")
        all_metrics.append(e1_m)
        ref_datasets["EXP01"] = (e1_labels, e1_scores)
        print(f"\nReference EXP01: AUC={e1_m['auc_roc']:.4f}")

    if args.step5_predictions and Path(args.step5_predictions).exists():
        s5_labels, s5_scores = load_predictions(Path(args.step5_predictions))
        s5_m = compute_metrics(s5_labels, s5_scores, "Step 5")
        ref_datasets["Step 5"] = (s5_labels, s5_scores)
        print(f"Reference Step 5: AUC={s5_m['auc_roc']:.4f}")

    # ── Comparative analysis ─────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("COMPARATIVE ANALYSIS")
    print("=" * 60)

    analysis_dir = results_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Metrics CSV
    metrics_csv = analysis_dir / "phase1_metrics.csv"
    with open(metrics_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(all_metrics[0].keys()))
        writer.writeheader()
        writer.writerows(all_metrics)
    print(f"  Metrics: {metrics_csv}")

    # Plots
    plot_phase1_roc(all_metrics, all_labels_scores, analysis_dir / "phase1_roc.png",
                    reference_datasets=ref_datasets)
    plot_auc_vs_ratio(all_metrics, analysis_dir / "phase1_auc_vs_ratio.png")
    plot_score_distributions(all_labels_scores, analysis_dir / "phase1_score_dist.png")

    # Report
    write_phase1_report(all_metrics, dict(seed_metrics), analysis_dir / "phase1_report.md")

    # ── Summary ──────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("PHASE 1 COMPLETE")
    print("=" * 60)
    print("\nResults summary:")
    for m in all_metrics:
        ci = f"[{m['auc_roc_ci_lo']:.4f}, {m['auc_roc_ci_hi']:.4f}]"
        print(f"  {m['name']:12s}  AUC-ROC = {m['auc_roc']:.4f}  {ci}  "
              f"({m['n_pos']} pos, {m['n_neg']} neg)")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\n[ERROR] {exc}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
