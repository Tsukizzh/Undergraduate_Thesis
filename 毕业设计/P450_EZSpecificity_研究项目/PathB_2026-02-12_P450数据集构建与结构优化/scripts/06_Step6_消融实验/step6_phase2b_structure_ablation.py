"""
Phase 2b: Structure-source ablation — ABL-04 bridge experiment.

Fills the missing cell in the 2x2 (Structure x NegType) matrix:
  - EXP01:  Crystal + Inhibitor  → AUC-ROC = 0.7115  (Step 2)
  - ABL-03: Vina    + Random     → AUC-ROC = 0.5170  (Phase 1 / Step 5)
  - ABL-04: Vina    + Inhibitor  → ???               (THIS experiment)
  - (Crystal + Random: impossible — no crystal structures for random pairs)

Gap decomposition:
  Total gap      = EXP01 - ABL-03 = 0.1945
  Structure eff  = EXP01 - ABL-04  (holding neg type = inhibitor)
  Neg-type eff   = ABL-04 - ABL-03 (holding structure = Vina)
  Check: Structure eff + Neg-type eff = Total gap

Pipeline:
  1. Merge Step 5 LMDB (pos entries) + Phase 2a LMDB (inhib neg entries)
  2. Create ABL-04 data.csv (Step 5 positives + Phase 2a inhibitor negatives)
  3. Run inference on merged data
  4. Matched cohort analysis: EXP01 vs ABL-04 on overlapping pairs
  5. 2x2 decomposition report

Usage:
    python step6_phase2b_structure_ablation.py \
        --step5_data        <results/04_Step4_批量对接/data.csv> \
        --step5_lmdb        <data/05_Step5_重构评估/structure_features.lmdb> \
        --step5_hqid        <data/05_Step5_重构评估/high_quality_id.txt> \
        --phase2a_dir       <data/06_Step6_消融实验/phase2_inhibitor_docking> \
        --exp01_predictions <results/02_Step2_因子实验/EXP01_.../predictions.csv> \
        --shared_features   <data/00_shared/features> \
        --checkpoint_dir    <saved_model/model/run_0> \
        --data_dir          <data/06_Step6_消融实验> \
        --results_dir       <results/06_Step6_消融实验> \
        [--phase1_metrics   <results/06_Step6_消融实验/analysis/phase1_metrics.csv>]
"""

from __future__ import annotations

import argparse
import struct
import sys
import warnings
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import torch
from rdkit import RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 42


# ═══════════════════════════════════════════════════════════════════════
#  1. LMDB MERGING
# ═══════════════════════════════════════════════════════════════════════

def merge_lmdbs(lmdb_paths: List[Path], output_path: Path) -> int:
    """Merge multiple structure_features.lmdb into one.

    Keys are Dock Index strings (e.g., b"42", b"4001"). No collision expected
    because Dock Index ranges are non-overlapping.

    Returns total number of entries written.
    """
    import lmdb

    output_path.parent.mkdir(parents=True, exist_ok=True)
    out_lmdb_dir = str(output_path)

    all_keys: Dict[bytes, bytes] = {}
    for lp in lmdb_paths:
        lp_str = str(lp)
        # LMDB on Windows with Unicode paths: use subdir=False for .lmdb files
        is_file = lp.is_file() or lp_str.endswith(".lmdb")
        env = lmdb.open(lp_str, readonly=True, lock=False, subdir=not is_file,
                        map_size=1 * 1024 * 1024 * 1024)
        with env.begin() as txn:
            cursor = txn.cursor()
            for key, value in cursor:
                if key in all_keys:
                    raise ValueError(
                        f"Duplicate key {key!r} found in {lp.name} — "
                        f"Dock Index collision between LMDBs"
                    )
                all_keys[key] = value
        env.close()
        print(f"  Read {lp.name}: {len(all_keys)} total entries so far")

    # Write merged LMDB
    map_size = sum(len(k) + len(v) for k, v in all_keys.items()) * 2 + 10 * 1024 * 1024
    is_out_file = out_lmdb_dir.endswith(".lmdb")
    env_out = lmdb.open(out_lmdb_dir, map_size=map_size, subdir=not is_out_file)
    with env_out.begin(write=True) as txn:
        for key, value in sorted(all_keys.items(), key=lambda kv: int(kv[0])):
            txn.put(key, value)
    env_out.close()
    print(f"  Merged LMDB: {len(all_keys)} entries → {output_path}")
    return len(all_keys)


def merge_high_quality_ids(id_paths: List[Path], output_path: Path) -> List[str]:
    """Merge high_quality_id.txt files (union)."""
    all_ids = set()
    for p in id_paths:
        if p.exists():
            ids = [l.strip() for l in p.read_text().strip().split("\n") if l.strip()]
            all_ids.update(ids)
            print(f"  Read {p.name}: {len(ids)} IDs")
    sorted_ids = sorted(all_ids, key=lambda x: int(x))
    output_path.write_text("\n".join(sorted_ids) + "\n")
    print(f"  Merged HQ IDs: {len(sorted_ids)} → {output_path}")
    return sorted_ids


# ═══════════════════════════════════════════════════════════════════════
#  2. INFERENCE (reused from Phase 1)
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
def run_inference(config, checkpoint_path: Path, device: torch.device):
    from Datasets.brenda import Singledataset
    from Models.ss import SS
    dm = Singledataset(config)
    test_loader = dm.test_dataloader()
    model = SS.load_from_checkpoint(str(checkpoint_path), config=config, map_location="cpu")
    model.eval()
    model.to(device)
    logits_chunks = []
    for batch in test_loader:
        batch = batch.to(device)
        logits, _ = model(batch)
        logits_chunks.append(logits.squeeze(-1).detach().float().cpu().numpy().ravel())
    logits_all = (
        np.concatenate(logits_chunks).astype(np.float32)
        if logits_chunks else np.zeros((0,), dtype=np.float32)
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
#  3. METRICS
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


def bootstrap_auc_diff(
    labels_a, scores_a, labels_b, scores_b,
    n_boot=BOOTSTRAP_N, seed=BOOTSTRAP_SEED,
):
    """Bootstrap the difference in AUC-ROC between two paired experiments.

    Both must have matching row order (matched cohort).
    Returns (mean_diff, ci_lo, ci_hi).
    """
    from sklearn.metrics import roc_auc_score
    n = len(labels_a)
    assert n == len(labels_b), "Matched cohort must have same length"
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.RandomState(seed)
    diffs = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.randint(0, n, size=n)
        la, sa = labels_a[idx], scores_a[idx]
        lb, sb = labels_b[idx], scores_b[idx]
        if la.sum() == 0 or la.sum() == n or lb.sum() == 0 or lb.sum() == n:
            diffs[i] = np.nan
        else:
            diffs[i] = roc_auc_score(la, sa) - roc_auc_score(lb, sb)
    diffs = diffs[~np.isnan(diffs)]
    if diffs.size == 0:
        return float("nan"), float("nan"), float("nan")
    return float(np.mean(diffs)), float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))


def compute_metrics(labels, scores, name="", seed=None):
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
    m = {
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
    if seed is not None:
        m["seed"] = seed
    return m


# ═══════════════════════════════════════════════════════════════════════
#  4. MATCHED COHORT ANALYSIS
# ═══════════════════════════════════════════════════════════════════════

def matched_cohort_analysis(
    exp01_df: pd.DataFrame,
    abl04_df: pd.DataFrame,
) -> dict:
    """Match EXP01 vs ABL-04 on (Enzyme Index, Substrate Index) pairs.

    Returns dict with matched metrics, AUC differences, and bootstrap CIs.
    """
    exp01_df = exp01_df.copy()
    abl04_df = abl04_df.copy()

    exp01_df["pair_key"] = (
        exp01_df["Enzyme Index"].astype(str) + "_" +
        exp01_df["Substrate Index"].astype(str)
    )
    abl04_df["pair_key"] = (
        abl04_df["Enzyme Index"].astype(str) + "_" +
        abl04_df["Substrate Index"].astype(str)
    )

    # Verify pair_key uniqueness (each experiment should have at most one row per pair)
    exp01_dups = exp01_df["pair_key"].duplicated().sum()
    abl04_dups = abl04_df["pair_key"].duplicated().sum()
    if exp01_dups > 0 or abl04_dups > 0:
        print(f"  [WARN] Duplicate pair_keys: EXP01={exp01_dups}, ABL-04={abl04_dups}")
        exp01_df = exp01_df.drop_duplicates(subset="pair_key", keep="first")
        abl04_df = abl04_df.drop_duplicates(subset="pair_key", keep="first")

    common_keys = set(exp01_df["pair_key"]) & set(abl04_df["pair_key"])
    print(f"  Matched cohort: {len(common_keys)} common pairs "
          f"(EXP01 has {len(exp01_df)}, ABL-04 has {len(abl04_df)})")

    exp01_matched = (
        exp01_df[exp01_df["pair_key"].isin(common_keys)]
        .sort_values("pair_key").reset_index(drop=True)
    )
    abl04_matched = (
        abl04_df[abl04_df["pair_key"].isin(common_keys)]
        .sort_values("pair_key").reset_index(drop=True)
    )

    # Verify labels match
    assert (exp01_matched["Label"].values == abl04_matched["Label"].values).all(), \
        "Label mismatch in matched cohort — pair_key collision?"

    labels = exp01_matched["Label"].values.astype(int)
    exp01_scores = exp01_matched["score"].values.astype(float)
    abl04_scores = abl04_matched["score"].values.astype(float)

    exp01_m = compute_metrics(labels, exp01_scores, name="EXP01 (matched)")
    abl04_m = compute_metrics(labels, abl04_scores, name="ABL-04 (matched)")

    _, ci_lo, ci_hi = bootstrap_auc_diff(
        labels, exp01_scores, labels, abl04_scores
    )
    # Use observed diff (not bootstrap mean) as point estimate
    observed_diff = exp01_m["auc_roc"] - abl04_m["auc_roc"]

    return {
        "n_matched": len(common_keys),
        "n_pos": int(labels.sum()),
        "n_neg": len(labels) - int(labels.sum()),
        "exp01_metrics": exp01_m,
        "abl04_metrics": abl04_m,
        "auc_diff_mean": round(observed_diff, 4),
        "auc_diff_ci_lo": round(ci_lo, 4),
        "auc_diff_ci_hi": round(ci_hi, 4),
        "exp01_matched_df": exp01_matched,
        "abl04_matched_df": abl04_matched,
    }


# ═══════════════════════════════════════════════════════════════════════
#  5. VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════

def plot_phase2_roc(all_metrics, all_labels_scores, output_path):
    """ROC curves for EXP01, ABL-03, ABL-04 (and matched variants)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.metrics import roc_curve, roc_auc_score

    fig, ax = plt.subplots(1, 1, figsize=(8, 7))
    colors = {"EXP01": "#2196F3", "ABL-03": "#4CAF50", "ABL-04": "#FF5722",
              "EXP01 (matched)": "#90CAF9", "ABL-04 (matched)": "#FFAB91"}

    for name, (labels, scores) in sorted(all_labels_scores.items()):
        if len(np.unique(labels)) < 2:
            print(f"  [WARN] Skipping ROC for {name}: single-class labels")
            continue
        fpr, tpr, _ = roc_curve(labels, scores)
        auc_val = roc_auc_score(labels, scores)
        color = colors.get(name, None)
        ls = "--" if "(matched)" in name else "-"
        ax.plot(fpr, tpr, label=f"{name} (AUC={auc_val:.4f})", color=color, ls=ls)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.3, label="Random")
    ax.set_xlabel("FPR")
    ax.set_ylabel("TPR")
    ax.set_title("Phase 2: Structure Source Ablation — ROC Curves")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150)
    plt.close(fig)
    print(f"  ROC plot → {output_path}")


def plot_2x2_heatmap(matrix_data, output_path):
    """2x2 heatmap of Structure x NegType AUC-ROC values."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    # matrix_data: dict with keys (struct, neg_type) → auc
    labels_struct = ["Crystal", "Vina"]
    labels_neg = ["Inhibitor", "Random"]
    matrix = np.full((2, 2), np.nan)
    for (s, n), v in matrix_data.items():
        si = labels_struct.index(s) if s in labels_struct else -1
        ni = labels_neg.index(n) if n in labels_neg else -1
        if si >= 0 and ni >= 0:
            matrix[si, ni] = v

    im = ax.imshow(matrix, cmap="RdYlGn", vmin=0.45, vmax=0.75, aspect="auto")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels_neg)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(labels_struct)
    ax.set_xlabel("Negative Type")
    ax.set_ylabel("Structure Source")
    ax.set_title("2×2 Ablation: AUC-ROC")

    for i in range(2):
        for j in range(2):
            val = matrix[i, j]
            if np.isnan(val):
                ax.text(j, i, "N/A\n(impossible)", ha="center", va="center",
                        fontsize=11, color="gray")
            else:
                ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                        fontsize=14, fontweight="bold")

    fig.colorbar(im, ax=ax, shrink=0.8, label="AUC-ROC")
    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150)
    plt.close(fig)
    print(f"  Heatmap → {output_path}")


def plot_score_distributions(all_labels_scores, output_path):
    """Score distributions for EXP01, ABL-03, ABL-04."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = [n for n in sorted(all_labels_scores.keys()) if "(matched)" not in n]
    n_plots = len(names)
    fig, axes = plt.subplots(1, n_plots, figsize=(5 * n_plots, 4), sharey=True)
    if n_plots == 1:
        axes = [axes]
    colors_pos = "#2196F3"
    colors_neg = "#FF5722"

    for ax, name in zip(axes, names):
        labels, scores = all_labels_scores[name]
        ax.hist(scores[labels == 1], bins=30, alpha=0.6, color=colors_pos,
                label=f"Pos (n={int(labels.sum())})", density=True)
        ax.hist(scores[labels == 0], bins=30, alpha=0.6, color=colors_neg,
                label=f"Neg (n={int((1 - labels).sum())})", density=True)
        ax.set_xlabel("Score (logit)")
        ax.set_title(name)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Density")
    fig.suptitle("Phase 2: Score Distributions", fontsize=13)
    fig.tight_layout()
    fig.savefig(str(output_path), dpi=150)
    plt.close(fig)
    print(f"  Score dist → {output_path}")


# ═══════════════════════════════════════════════════════════════════════
#  6. REPORT
# ═══════════════════════════════════════════════════════════════════════

def write_phase2_report(
    report_path: Path,
    all_metrics: List[dict],
    matched: dict,
    decomposition: dict,
):
    lines = [
        "# Phase 2: Structure-Source Ablation Report",
        "",
        f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Branch**: pathb-ablation",
        "",
        "## 2×2 Matrix Design",
        "",
        "| | Inhibitor Neg | Random Neg |",
        "|---|---|---|",
        f"| **Crystal** | EXP01 = **{decomposition['exp01_auc']:.4f}** | N/A (impossible) |",
        f"| **Vina** | ABL-04 = **{decomposition['abl04_auc']:.4f}** | ABL-03 = **{decomposition['abl03_auc']:.4f}** |",
        "",
        "## All Experiment Results",
        "",
        "| Experiment | AUC-ROC | 95% CI | AUC-PR | Pos | Neg | Score Sep |",
        "|-----------|---------|--------|--------|-----|-----|----------|",
    ]
    for m in all_metrics:
        lines.append(
            f"| {m['name']} | **{m['auc_roc']:.4f}** | "
            f"[{m['auc_roc_ci_lo']:.4f}, {m['auc_roc_ci_hi']:.4f}] | "
            f"{m['auc_pr']:.4f} | {m['n_pos']} | {m['n_neg']} | "
            f"{m['score_separation']:.4f} |"
        )

    lines.extend([
        "",
        "## Gap Decomposition",
        "",
        f"**Total gap** = EXP01 - ABL-03 = {decomposition['total_gap']:.4f}",
        "",
        f"**Structure effect** = EXP01 - ABL-04 = {decomposition['structure_effect']:.4f}",
        f"  (Crystal → Vina, holding neg type = inhibitor)",
        "",
        f"**Neg-type effect** = ABL-04 - ABL-03 = {decomposition['negtype_effect']:.4f}",
        f"  (Inhibitor → Random, holding structure = Vina)",
        "",
        f"**Consistency check**: {decomposition['structure_effect']:.4f} + "
        f"{decomposition['negtype_effect']:.4f} = "
        f"{decomposition['structure_effect'] + decomposition['negtype_effect']:.4f} "
        f"(total gap = {decomposition['total_gap']:.4f})",
        "",
        "### Relative Contributions",
        "",
    ])

    total = abs(decomposition["total_gap"]) if decomposition["total_gap"] != 0 else 1
    struct_pct = decomposition["structure_effect"] / total * 100
    neg_pct = decomposition["negtype_effect"] / total * 100
    lines.extend([
        f"- Structure source: {struct_pct:.1f}% of total gap",
        f"- Negative type: {neg_pct:.1f}% of total gap",
        "",
    ])

    # Interpretation
    lines.extend([
        "## Interpretation",
        "",
    ])
    if abs(decomposition["negtype_effect"]) > abs(decomposition["structure_effect"]) * 2:
        lines.append(
            "**Negative type is the dominant factor.** Switching from inhibitor to random "
            "negatives causes most of the performance drop. This confirms that the task "
            "difficulty (inhibitor vs random) is the primary driver, not structure quality."
        )
    elif abs(decomposition["structure_effect"]) > abs(decomposition["negtype_effect"]) * 2:
        lines.append(
            "**Structure source is the dominant factor.** Switching from crystal to Vina "
            "structures causes most of the performance drop. Improving docking quality "
            "would be the most impactful intervention."
        )
    else:
        lines.append(
            "**Both factors contribute comparably.** Neither structure source nor negative "
            "type alone explains the gap. Both need to be addressed for significant improvement."
        )

    # Matched cohort
    lines.extend([
        "",
        "## Matched Cohort Analysis (EXP01 vs ABL-04)",
        "",
        f"- Matched pairs: {matched['n_matched']} "
        f"({matched['n_pos']} pos + {matched['n_neg']} neg)",
        f"- EXP01 AUC (matched): {matched['exp01_metrics']['auc_roc']:.4f}",
        f"- ABL-04 AUC (matched): {matched['abl04_metrics']['auc_roc']:.4f}",
        f"- AUC difference: {matched['auc_diff_mean']:.4f} "
        f"[{matched['auc_diff_ci_lo']:.4f}, {matched['auc_diff_ci_hi']:.4f}]",
        "",
        "The matched cohort controls for enzyme-substrate pair composition, "
        "isolating the pure structure-source effect.",
    ])

    lines.extend([
        "",
        "## Files Generated",
        "",
        "- `phase2_metrics.csv` — Experiment metrics",
        "- `phase2_roc.png` — ROC curve comparison",
        "- `phase2_heatmap.png` — 2×2 ablation heatmap",
        "- `phase2_score_dist.png` — Score distributions",
        "- `ABL-04_*/predictions.csv` — ABL-04 predictions",
        "",
    ])

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"  Report → {report_path}")


# ═══════════════════════════════════════════════════════════════════════
#  7. MAIN
# ═══════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description="Phase 2b: Structure-source ablation")
    p.add_argument("--step5_data", required=True,
                   help="Step 5 data.csv (2766 rows: 265 pos + 2501 neg)")
    p.add_argument("--step5_lmdb", required=True,
                   help="Step 5 structure_features.lmdb")
    p.add_argument("--step5_hqid", required=True,
                   help="Step 5 high_quality_id.txt")
    p.add_argument("--phase2a_dir", required=True,
                   help="Phase 2a output dir (inhibitor docking)")
    p.add_argument("--exp01_predictions", required=True,
                   help="EXP01 predictions.csv (Crystal+Inhibitor)")
    p.add_argument("--shared_features", required=True,
                   help="Shared features dir (enzyme, reaction, grover, morgan)")
    p.add_argument("--checkpoint_dir", required=True,
                   help="Model checkpoint dir")
    p.add_argument("--data_dir", required=True,
                   help="Intermediate data output dir")
    p.add_argument("--results_dir", required=True,
                   help="Final results output dir")
    p.add_argument("--phase1_metrics", default="",
                   help="Phase 1 metrics CSV (for ABL-03 reference)")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    step5_data = Path(args.step5_data).resolve()
    step5_lmdb = Path(args.step5_lmdb).resolve()
    step5_hqid = Path(args.step5_hqid).resolve()
    phase2a_dir = Path(args.phase2a_dir).resolve()
    exp01_pred_path = Path(args.exp01_predictions).resolve()
    shared_features = Path(args.shared_features).resolve()
    checkpoint_dir = Path(args.checkpoint_dir).resolve()
    data_dir = Path(args.data_dir).resolve()
    results_dir = Path(args.results_dir).resolve()

    # Validate inputs
    for p, label in [
        (step5_data, "Step 5 data.csv"),
        (step5_lmdb, "Step 5 LMDB"),
        (step5_hqid, "Step 5 high_quality_id.txt"),
        (exp01_pred_path, "EXP01 predictions"),
    ]:
        if not p.exists():
            print(f"[ERROR] {label} not found: {p}")
            return 1

    phase2a_lmdb = phase2a_dir / "structure_features.lmdb"
    phase2a_hqid = phase2a_dir / "high_quality_id.txt"
    phase2a_data = phase2a_dir / "inhibitor_data.csv"
    for p, label in [
        (phase2a_lmdb, "Phase 2a LMDB"),
        (phase2a_hqid, "Phase 2a high_quality_id.txt"),
        (phase2a_data, "Phase 2a inhibitor_data.csv"),
    ]:
        if not p.exists():
            print(f"[ERROR] {label} not found: {p}")
            print("  Have you run step6_phase2a_dock_inhibitors.py first?")
            return 1

    # Setup project root for model imports
    project_root = find_project_root(Path(__file__))
    src_dir = project_root / "src"
    if str(src_dir) not in sys.path:
        sys.path.insert(0, str(src_dir))

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Output dirs
    abl04_data_dir = data_dir / "ABL-04_vina_inhibitor"
    abl04_results_dir = results_dir / "ABL-04_vina_inhibitor"
    analysis_dir = results_dir / "analysis"
    for d in [abl04_data_dir, abl04_results_dir, analysis_dir]:
        d.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Merge LMDBs ──────────────────────────────────────
    print("\n" + "=" * 60)
    print("STEP 1: Merge structure feature LMDBs")
    print("=" * 60)

    merged_lmdb = abl04_data_dir / "structure_features.lmdb"
    merged_hqid = abl04_data_dir / "high_quality_id.txt"

    if merged_lmdb.exists():
        print(f"  Merged LMDB already exists: {merged_lmdb}")
        print(f"  [NOTE] Delete {merged_lmdb} to force re-merge if inputs changed")
    else:
        merge_lmdbs([step5_lmdb, phase2a_lmdb], merged_lmdb)

    hq_ids = merge_high_quality_ids([step5_hqid, phase2a_hqid], merged_hqid)
    print(f"  Total high-quality IDs: {len(hq_ids)}")

    # ── Step 2: Create ABL-04 data.csv ────────────────────────────
    print("\n" + "=" * 60)
    print("STEP 2: Create ABL-04 data.csv")
    print("=" * 60)

    step5_df = pd.read_csv(step5_data)
    pos_df = step5_df[step5_df["Label"] == 1].copy()
    inhib_df = pd.read_csv(phase2a_data)

    abl04_csv = abl04_data_dir / "data.csv"
    combined = pd.concat([pos_df, inhib_df], ignore_index=True)
    combined = combined.sort_values("Dock Index").reset_index(drop=True)
    combined.to_csv(abl04_csv, index=False)
    print(f"  ABL-04 data.csv: {len(pos_df)} pos + {len(inhib_df)} neg = "
          f"{len(combined)} total → {abl04_csv}")

    # Verify Dock Indices exist in high_quality_id set
    hq_set = set(int(x) for x in hq_ids)
    data_docks = set(combined["Dock Index"].astype(int))
    in_hq = data_docks & hq_set
    coverage = len(in_hq) / len(data_docks) if data_docks else 0
    print(f"  Dock Index coverage: {len(in_hq)}/{len(data_docks)} "
          f"({coverage:.1%}) in high_quality_id.txt")
    if len(in_hq) == 0:
        print("[ERROR] No Dock Indices in HQ set. Check LMDB merge.")
        return 1
    if coverage < 0.8:
        print(f"  [WARN] Low coverage ({coverage:.1%}). Some pairs may be missing.")

    # ── Step 3: Run ABL-04 inference ──────────────────────────────
    print("\n" + "=" * 60)
    print("STEP 3: Run ABL-04 inference")
    print("=" * 60)

    pred_csv = abl04_results_dir / "predictions.csv"
    if pred_csv.exists():
        print(f"  Predictions already exist: {pred_csv}")
        abl04_pred = pd.read_csv(pred_csv)
        if len(abl04_pred) != len(combined):
            print(f"  [WARN] Cached predictions ({len(abl04_pred)}) != "
                  f"current data ({len(combined)}). Re-running inference.")
            pred_csv.unlink()
    if not pred_csv.exists():
        abl04_pred = run_single_experiment(
            dataset_csv=abl04_csv,
            structure_lmdb=merged_lmdb,
            high_quality_ids=merged_hqid,
            shared_features=shared_features,
            checkpoint_dir=checkpoint_dir,
            experiment_name="ABL-04",
            device=device,
        )
        abl04_pred.to_csv(pred_csv, index=False)
        print(f"  Predictions saved: {len(abl04_pred)} rows → {pred_csv}")

    # Compute ABL-04 metrics
    abl04_labels = abl04_pred["Label"].values.astype(int)
    abl04_scores = abl04_pred["score"].values.astype(float)
    abl04_metrics = compute_metrics(abl04_labels, abl04_scores, name="ABL-04")
    print(f"  ABL-04 AUC-ROC: {abl04_metrics['auc_roc']:.4f} "
          f"[{abl04_metrics['auc_roc_ci_lo']:.4f}, {abl04_metrics['auc_roc_ci_hi']:.4f}]")

    # ── Step 4: Load reference experiments ────────────────────────
    print("\n" + "=" * 60)
    print("STEP 4: Load reference experiments")
    print("=" * 60)

    # EXP01
    exp01_pred = pd.read_csv(exp01_pred_path)
    exp01_labels = exp01_pred["Label"].values.astype(int)
    exp01_scores = exp01_pred["score"].values.astype(float)
    exp01_metrics = compute_metrics(exp01_labels, exp01_scores, name="EXP01")
    print(f"  EXP01 AUC-ROC: {exp01_metrics['auc_roc']:.4f} ({len(exp01_pred)} samples)")

    # ABL-03 (from Phase 1 metrics or Step 5)
    abl03_auc = None
    if args.phase1_metrics and Path(args.phase1_metrics).exists():
        p1 = pd.read_csv(args.phase1_metrics)
        abl03_row = p1[p1["name"] == "ABL-03"]
        if not abl03_row.empty:
            abl03_auc = float(abl03_row["auc_roc"].iloc[0])
            print(f"  ABL-03 AUC-ROC: {abl03_auc:.4f} (from Phase 1 metrics)")

    if abl03_auc is None:
        abl03_auc = 0.5170
        print(f"  ABL-03 AUC-ROC: {abl03_auc:.4f} (hardcoded from Phase 1)")

    # ── Step 5: Matched cohort analysis ───────────────────────────
    print("\n" + "=" * 60)
    print("STEP 5: Matched cohort analysis (EXP01 vs ABL-04)")
    print("=" * 60)

    matched = matched_cohort_analysis(exp01_pred, abl04_pred)
    print(f"  EXP01 AUC (matched): {matched['exp01_metrics']['auc_roc']:.4f}")
    print(f"  ABL-04 AUC (matched): {matched['abl04_metrics']['auc_roc']:.4f}")
    print(f"  AUC diff: {matched['auc_diff_mean']:.4f} "
          f"[{matched['auc_diff_ci_lo']:.4f}, {matched['auc_diff_ci_hi']:.4f}]")

    # ── Step 6: 2×2 Decomposition ─────────────────────────────────
    print("\n" + "=" * 60)
    print("STEP 6: 2×2 Gap Decomposition")
    print("=" * 60)

    decomposition = {
        "exp01_auc": exp01_metrics["auc_roc"],
        "abl04_auc": abl04_metrics["auc_roc"],
        "abl03_auc": abl03_auc,
        "total_gap": exp01_metrics["auc_roc"] - abl03_auc,
        "structure_effect": exp01_metrics["auc_roc"] - abl04_metrics["auc_roc"],
        "negtype_effect": abl04_metrics["auc_roc"] - abl03_auc,
    }
    print(f"  Total gap (EXP01 - ABL-03): {decomposition['total_gap']:.4f}")
    print(f"  Structure effect (EXP01 - ABL-04): {decomposition['structure_effect']:.4f}")
    print(f"  Neg-type effect (ABL-04 - ABL-03): {decomposition['negtype_effect']:.4f}")

    # ── Step 7: Visualization ─────────────────────────────────────
    print("\n" + "=" * 60)
    print("STEP 7: Generate plots")
    print("=" * 60)

    all_labels_scores = {
        "EXP01": (exp01_labels, exp01_scores),
        "ABL-04": (abl04_labels, abl04_scores),
    }
    if "exp01_matched_df" in matched:
        m_exp01 = matched["exp01_matched_df"]
        m_abl04 = matched["abl04_matched_df"]
        all_labels_scores["EXP01 (matched)"] = (
            m_exp01["Label"].values.astype(int),
            m_exp01["score"].values.astype(float),
        )
        all_labels_scores["ABL-04 (matched)"] = (
            m_abl04["Label"].values.astype(int),
            m_abl04["score"].values.astype(float),
        )

    plot_phase2_roc(
        [exp01_metrics, abl04_metrics],
        all_labels_scores,
        analysis_dir / "phase2_roc.png",
    )

    plot_2x2_heatmap(
        {("Crystal", "Inhibitor"): exp01_metrics["auc_roc"],
         ("Vina", "Inhibitor"): abl04_metrics["auc_roc"],
         ("Vina", "Random"): abl03_auc,
         ("Crystal", "Random"): float("nan")},
        analysis_dir / "phase2_heatmap.png",
    )

    plot_score_distributions(
        all_labels_scores,
        analysis_dir / "phase2_score_dist.png",
    )

    # ── Step 8: Save metrics and report ───────────────────────────
    print("\n" + "=" * 60)
    print("STEP 8: Save metrics and report")
    print("=" * 60)

    all_metrics = [exp01_metrics, abl04_metrics]
    metrics_csv = analysis_dir / "phase2_metrics.csv"
    pd.DataFrame(all_metrics).to_csv(metrics_csv, index=False)
    print(f"  Metrics → {metrics_csv}")

    write_phase2_report(
        analysis_dir / "phase2_report.md",
        all_metrics, matched, decomposition,
    )

    # ── Summary ───────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("PHASE 2B COMPLETE")
    print(f"{'=' * 60}")
    print(f"  ABL-04 AUC-ROC: {abl04_metrics['auc_roc']:.4f}")
    print(f"  Structure effect: {decomposition['structure_effect']:.4f}")
    print(f"  Neg-type effect: {decomposition['negtype_effect']:.4f}")
    total = abs(decomposition["total_gap"]) if decomposition["total_gap"] != 0 else 1
    print(f"  Structure %: {decomposition['structure_effect'] / total * 100:.1f}%")
    print(f"  Neg-type %: {decomposition['negtype_effect'] / total * 100:.1f}%")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\n[ERROR] {exc}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
