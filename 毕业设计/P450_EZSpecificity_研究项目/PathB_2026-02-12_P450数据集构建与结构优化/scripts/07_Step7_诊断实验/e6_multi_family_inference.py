"""
E6 Expansion: Multi-Family Inference + Per-Substrate Promiscuity Analysis
Tests: Can the model distinguish positive from negative across different enzyme families?
Method: Run trained model inference on 6 small families' test sets, compute AUC per family,
then analyze AUC vs per-substrate promiscuity (within-family and cross-family).
"""
from __future__ import annotations

import sys
import os
import warnings
import json
import tempfile
import shutil

import numpy as np
import pandas as pd
import torch
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

# --- Paths ---
PROJECT_ROOT = Path(r"D:\EZSpecificity_Project")
PATHB = PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目" / "PathB_2026-02-12_P450数据集构建与结构优化"
G_SMALL_FAMILY = Path(r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\small_family")
LOCAL_TMP = PROJECT_ROOT / "tmp_lmdb"
OUTPUT_DIR = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E6_expansion_多家族推理"
CHECKPOINT_DIR = PROJECT_ROOT / "saved_model" / "model" / "run_0"

# P450 predictions for comparison
P450_PREDICTIONS = PATHB / "results" / "05_Step5_重构评估" / "predictions.csv"

FAMILIES = {
    "Duf": {"dir_name": "Duf"},
    "Esterase": {"dir_name": "Esterase"},
    "Gt_acceptor": {"dir_name": "Gt_acceptor"},
    "Nitrilase": {"dir_name": "Nitrilase"},
    "Phosphatase": {"dir_name": "Phosphatase"},
    "Thiolase": {"dir_name": "Thiolase"},
}

# Column mapping: small family CSV → model expected
COL_MAP = {
    "reaction": "Substrate Index",
    "enzyme": "Enzyme Index",
    "label": "Label",
    "structure_index": "Dock Index",
}


def _find_config_yaml():
    yaml_path = CHECKPOINT_DIR / "complete-full-random-all-0-complex.yml"
    if yaml_path.is_file():
        return yaml_path
    raise FileNotFoundError(f"No config YAML in {CHECKPOINT_DIR}")


def _find_checkpoint():
    ckpt = CHECKPOINT_DIR / "models" / "best-checkpoint.ckpt"
    if not ckpt.is_file():
        candidates = sorted(CHECKPOINT_DIR.rglob("best-checkpoint.ckpt"))
        if not candidates:
            raise FileNotFoundError(f"No best-checkpoint.ckpt under {CHECKPOINT_DIR}")
        ckpt = candidates[0]
    return ckpt


def prepare_family_csv(family_name):
    """Prepare a data CSV with renamed columns for model compatibility."""
    csv_path = LOCAL_TMP / f"{family_name}_testing_datas_0.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"Missing {csv_path}")

    df = pd.read_csv(csv_path)
    df = df.rename(columns=COL_MAP)
    for col in ["Substrate Index", "Enzyme Index", "Dock Index", "Label"]:
        if col in df.columns:
            df[col] = df[col].astype(int)
    return df


def build_config_for_family(config_yaml, family_name, tmp_csv_path):
    """Build a config for single-family inference."""
    sys.path.insert(0, str(PROJECT_ROOT / "src"))
    from utils import load_config

    config = load_config(str(config_yaml))

    config.num_cpus = 0
    config.num_gpus = 1

    config.data.tag = f"e6_{family_name}_inference"
    config.data.representer = "structure_sequence"

    config.data.train_data_path = str(tmp_csv_path)
    config.data.val_data_path = str(tmp_csv_path)
    config.data.test_data_path = str(tmp_csv_path)

    fam_dir = G_SMALL_FAMILY / FAMILIES[family_name]["dir_name"]

    config.data.enzyme_lmdb_path = str(fam_dir / "enzyme_features.lmdb")
    config.data.reaction_lmdb_path = str(fam_dir / "reaction_features.lmdb")
    config.data.grover_path = str(fam_dir / "grover_fingerprint.lmdb")
    config.data.morgan_path = str(fam_dir / "morgan_fingerprint.npy")
    config.data.structure_processed_path = str(fam_dir / "structure" / "af2" / "structure_features.lmdb")
    config.data.sequence_processed_path = "nonexistent_str_features.lmdb"
    config.data.high_quality_id_path = "nonexistent_hq.txt"

    config.data.full_data = False
    config.data.batch_size = 16
    config.data.sample_weight = [1.0, 1.0]
    config.data.fake_sequence_ratio = 0
    config.data.max_substrate_length = 280
    config.data.max_enzyme_length = 1450
    config.data.features = ["morgan", 1024, "grover_mean", 4885]
    config.data.atom_features = ["grover", 2400]

    return config


def _sigmoid_np(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out.astype(np.float32)


@torch.no_grad()
def run_inference(config, checkpoint_path, device):
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    dm = Singledataset(config)
    test_loader = dm.test_dataloader()

    model = SS.load_from_checkpoint(str(checkpoint_path), config=config, map_location="cpu")
    model.eval()
    model.to(device)

    logits_chunks = []
    n_processed = 0
    n_total = len(dm.test_prediction_df)

    for batch_idx, batch in enumerate(test_loader):
        batch = batch.to(device)
        logits, _tags = model(batch)
        logits_np = logits.squeeze(-1).detach().float().cpu().numpy().ravel()
        logits_chunks.append(logits_np)
        n_processed += len(logits_np)
        if (batch_idx + 1) % 20 == 0:
            print(f"    Batch {batch_idx+1}: {n_processed}/{n_total}")

    logits_all = np.concatenate(logits_chunks).astype(np.float32) if logits_chunks else np.zeros(0, dtype=np.float32)
    pred_df = dm.test_prediction_df.copy()
    assert len(pred_df) == len(logits_all), f"Mismatch: df={len(pred_df)}, logits={len(logits_all)}"

    pred_df["logit"] = logits_all
    pred_df["prob"] = _sigmoid_np(logits_all)
    return pred_df


def compute_promiscuity(family_name):
    """Compute per-substrate promiscuity from full data.csv (not just test set)."""
    full_csv = LOCAL_TMP / f"{family_name}_data.csv"
    df = pd.read_csv(full_csv)
    pos = df[df["label"] == 1]
    promiscuity = pos.groupby("reaction")["enzyme"].nunique().reset_index()
    promiscuity.columns = ["Substrate Index", "promiscuity"]
    return promiscuity


def compute_per_substrate_auc(pred_df, min_pos=2, min_neg=2):
    """Compute AUC for each substrate with sufficient positive and negative samples."""
    results = []
    for sub_idx, group in pred_df.groupby("Substrate Index"):
        n_pos = (group["Label"] == 1).sum()
        n_neg = (group["Label"] == 0).sum()
        if n_pos >= min_pos and n_neg >= min_neg:
            auc = roc_auc_score(group["Label"], group["logit"])
            results.append({"Substrate Index": sub_idx, "auc": auc, "n_pos": n_pos, "n_neg": n_neg, "n_total": len(group)})
    return pd.DataFrame(results) if results else pd.DataFrame(columns=["Substrate Index", "auc", "n_pos", "n_neg", "n_total"])


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

    config_yaml = _find_config_yaml()
    checkpoint_path = _find_checkpoint()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    print(f"Config: {config_yaml}")
    print(f"Checkpoint: {checkpoint_path}")

    # --- Run inference for each family ---
    all_family_results = {}
    all_pred_dfs = {}

    for fam_name in FAMILIES:
        print(f"\n{'='*60}")
        print(f"  Family: {fam_name}")
        print(f"{'='*60}")

        # Prepare CSV with renamed columns
        df = prepare_family_csv(fam_name)
        print(f"  Test samples: {len(df)} (pos={int((df['Label']==1).sum())}, neg={int((df['Label']==0).sum())})")

        # Save temp CSV
        tmp_csv = OUTPUT_DIR / f"{fam_name}_test_renamed.csv"
        df.to_csv(tmp_csv, index=False)

        # Build config
        config = build_config_for_family(config_yaml, fam_name, tmp_csv)

        # Run inference
        print(f"  Running inference...")
        pred_df = run_inference(config, checkpoint_path, device)
        print(f"  Valid predictions: {len(pred_df)}")

        # Save predictions
        pred_csv = OUTPUT_DIR / f"{fam_name}_predictions.csv"
        pred_df.to_csv(pred_csv, index=False)

        # Compute metrics
        if len(pred_df) > 0 and pred_df["Label"].nunique() == 2:
            auc_roc = roc_auc_score(pred_df["Label"], pred_df["logit"])
            auc_pr = average_precision_score(pred_df["Label"], pred_df["logit"])
            prevalence = pred_df["Label"].mean()
        else:
            auc_roc = auc_pr = prevalence = float("nan")

        # Compute promiscuity
        promiscuity = compute_promiscuity(fam_name)
        mean_prom = promiscuity["promiscuity"].mean()
        median_prom = promiscuity["promiscuity"].median()

        # Per-substrate AUC
        per_sub_auc = compute_per_substrate_auc(pred_df)
        if len(per_sub_auc) > 0:
            per_sub_auc = per_sub_auc.merge(promiscuity, on="Substrate Index", how="left")
            per_sub_auc["family"] = fam_name

        result = {
            "family": fam_name,
            "n_test_input": len(df),
            "n_valid_predictions": len(pred_df),
            "n_pos": int(pred_df["Label"].sum()) if len(pred_df) > 0 else 0,
            "n_neg": int((pred_df["Label"] == 0).sum()) if len(pred_df) > 0 else 0,
            "auc_roc": float(auc_roc),
            "auc_pr": float(auc_pr),
            "prevalence": float(prevalence),
            "promiscuity_mean": float(mean_prom),
            "promiscuity_median": float(median_prom),
            "n_substrates_with_auc": len(per_sub_auc),
            "per_substrate_auc_mean": float(per_sub_auc["auc"].mean()) if len(per_sub_auc) > 0 else float("nan"),
            "logit_mean": float(pred_df["logit"].mean()) if len(pred_df) > 0 else float("nan"),
            "logit_std": float(pred_df["logit"].std()) if len(pred_df) > 0 else float("nan"),
        }
        all_family_results[fam_name] = result
        all_pred_dfs[fam_name] = per_sub_auc

        print(f"  AUC-ROC: {auc_roc:.4f}, AUC-PR: {auc_pr:.4f}")
        print(f"  Promiscuity: mean={mean_prom:.1f}, median={median_prom:.0f}")
        print(f"  Per-substrate AUC (n={len(per_sub_auc)}): mean={result['per_substrate_auc_mean']:.4f}" if len(per_sub_auc) > 0 else "  No per-substrate AUC")

    # --- Add P450 for comparison ---
    print(f"\n{'='*60}")
    print(f"  P450 (from Step 5)")
    print(f"{'='*60}")
    p450_df = pd.read_csv(P450_PREDICTIONS)
    p450_auc = roc_auc_score(p450_df["Label"], p450_df["logit"])
    p450_aupr = average_precision_score(p450_df["Label"], p450_df["logit"])
    # P450 promiscuity ~ 1.0 by definition (our Direction A dataset)
    all_family_results["P450"] = {
        "family": "P450",
        "n_test_input": len(p450_df),
        "n_valid_predictions": len(p450_df),
        "n_pos": int(p450_df["Label"].sum()),
        "n_neg": int((p450_df["Label"] == 0).sum()),
        "auc_roc": float(p450_auc),
        "auc_pr": float(p450_aupr),
        "prevalence": float(p450_df["Label"].mean()),
        "promiscuity_mean": 1.0,
        "promiscuity_median": 1.0,
        "n_substrates_with_auc": 0,
        "per_substrate_auc_mean": float("nan"),
        "logit_mean": float(p450_df["logit"].mean()),
        "logit_std": float(p450_df["logit"].std()),
    }
    print(f"  AUC-ROC: {p450_auc:.4f}, AUC-PR: {p450_aupr:.4f}")

    # --- Summary Table ---
    print(f"\n\n{'='*80}")
    print(f"{'Family':<15} {'N':>6} {'Valid':>6} {'AUC-ROC':>9} {'AUC-PR':>9} {'Prom_mean':>10} {'Prom_med':>9}")
    print("-" * 80)
    for fam in list(FAMILIES.keys()) + ["P450"]:
        r = all_family_results[fam]
        print(f"{r['family']:<15} {r['n_test_input']:>6} {r['n_valid_predictions']:>6} "
              f"{r['auc_roc']:>9.4f} {r['auc_pr']:>9.4f} "
              f"{r['promiscuity_mean']:>10.1f} {r['promiscuity_median']:>9.0f}")

    # --- Combine per-substrate AUCs across families ---
    all_per_sub = pd.concat([df for df in all_pred_dfs.values() if len(df) > 0], ignore_index=True)
    if len(all_per_sub) > 0:
        all_per_sub.to_csv(OUTPUT_DIR / "all_per_substrate_auc.csv", index=False)
        print(f"\nPer-substrate AUC data: {len(all_per_sub)} substrates across families")

    # --- Visualization ---
    _plot_family_comparison(all_family_results, OUTPUT_DIR)
    if len(all_per_sub) > 0:
        _plot_promiscuity_vs_auc(all_per_sub, all_family_results, OUTPUT_DIR)
        _plot_within_family_promiscuity(all_per_sub, OUTPUT_DIR)

    # --- Save results ---
    with open(OUTPUT_DIR / "e6_expansion_results.json", "w", encoding="utf-8") as f:
        json.dump(all_family_results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


def _plot_family_comparison(results, output_dir):
    """Bar chart comparing AUC-ROC across families, sorted by promiscuity."""
    families = sorted(results.keys(), key=lambda f: results[f]["promiscuity_mean"], reverse=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # AUC-ROC bar chart
    ax = axes[0]
    aucs = [results[f]["auc_roc"] for f in families]
    proms = [results[f]["promiscuity_mean"] for f in families]
    colors = ["crimson" if f == "P450" else "steelblue" for f in families]
    bars = ax.bar(range(len(families)), aucs, color=colors, alpha=0.7)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5, label="Random")
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels([f"{f}\n(prom={p:.0f})" for f, p in zip(families, proms)], fontsize=8, rotation=0)
    ax.set_ylabel("AUC-ROC")
    ax.set_title("AUC-ROC by Family (sorted by promiscuity)")
    ax.set_ylim(0.3, 1.0)
    ax.legend()
    for bar, auc in zip(bars, aucs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, f"{auc:.3f}",
                ha="center", fontsize=8, fontweight="bold")

    # Scatter: promiscuity vs AUC
    ax = axes[1]
    for f in families:
        r = results[f]
        color = "crimson" if f == "P450" else "steelblue"
        marker = "D" if f == "P450" else "o"
        ax.scatter(r["promiscuity_mean"], r["auc_roc"], s=100, c=color, marker=marker, zorder=5)
        ax.annotate(f, (r["promiscuity_mean"], r["auc_roc"]), textcoords="offset points",
                    xytext=(5, 5), fontsize=8, fontweight="bold" if f == "P450" else "normal")
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Mean Substrate Promiscuity")
    ax.set_ylabel("AUC-ROC")
    ax.set_title("Family AUC-ROC vs Mean Promiscuity")

    plt.tight_layout()
    plt.savefig(output_dir / "e6_family_comparison.png", dpi=150)
    plt.close()


def _plot_promiscuity_vs_auc(per_sub_df, family_results, output_dir):
    """Scatter plot: per-substrate promiscuity vs per-substrate AUC (cross-family)."""
    fig, ax = plt.subplots(figsize=(10, 7))

    cmap = plt.get_cmap("tab10")
    families = sorted(per_sub_df["family"].unique())

    for i, fam in enumerate(families):
        fam_data = per_sub_df[per_sub_df["family"] == fam]
        color = cmap(i % 10)
        family_auc = family_results[fam]["auc_roc"]
        ax.scatter(fam_data["promiscuity"], fam_data["auc"],
                   s=30, alpha=0.6, c=[color], label=f"{fam} (AUC={family_auc:.3f})", zorder=3)

    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5, label="Random")
    ax.set_xlabel("Per-Substrate Promiscuity (# enzymes catalyzing this substrate)")
    ax.set_ylabel("Per-Substrate AUC-ROC")
    ax.set_title("Per-Substrate AUC vs Promiscuity (Cross-Family)")
    ax.legend(fontsize=8, loc="lower right")
    ax.set_xscale("log")

    plt.tight_layout()
    plt.savefig(output_dir / "e6_promiscuity_vs_auc_scatter.png", dpi=150)
    plt.close()


def _plot_within_family_promiscuity(per_sub_df, output_dir):
    """Box plots of per-substrate AUC grouped by promiscuity bins, one subplot per family."""
    families = sorted(per_sub_df["family"].unique())
    n_fam = len(families)
    cols = min(3, n_fam)
    rows = (n_fam + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows), squeeze=False)

    for i, fam in enumerate(families):
        ax = axes[i // cols][i % cols]
        fam_data = per_sub_df[per_sub_df["family"] == fam].copy()

        if len(fam_data) < 3:
            ax.text(0.5, 0.5, f"{fam}\n(n={len(fam_data)}, too few)", ha="center", va="center", transform=ax.transAxes)
            ax.set_title(fam)
            continue

        # Create promiscuity bins
        prom_vals = fam_data["promiscuity"].dropna()
        if prom_vals.nunique() <= 3:
            fam_data["prom_bin"] = fam_data["promiscuity"].astype(str)
        else:
            try:
                fam_data["prom_bin"] = pd.qcut(fam_data["promiscuity"], q=min(4, prom_vals.nunique()), duplicates="drop")
            except ValueError:
                fam_data["prom_bin"] = pd.cut(fam_data["promiscuity"], bins=3)

        bins = sorted(fam_data["prom_bin"].unique(), key=lambda x: float(str(x).split(",")[0].strip("([")) if "," in str(x) else float(str(x)))
        box_data = [fam_data[fam_data["prom_bin"] == b]["auc"].values for b in bins]
        box_labels = [f"{b}\n(n={len(d)})" for b, d in zip(bins, box_data)]

        bp = ax.boxplot(box_data, labels=box_labels, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("steelblue")
            patch.set_alpha(0.5)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_ylabel("Per-Substrate AUC")
        ax.set_xlabel("Promiscuity Bin")
        ax.set_title(f"{fam} (n={len(fam_data)})")
        ax.tick_params(axis="x", labelsize=7)

    # Hide empty subplots
    for j in range(i + 1, rows * cols):
        axes[j // cols][j % cols].set_visible(False)

    plt.suptitle("Within-Family: Per-Substrate AUC by Promiscuity Bin", fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "e6_within_family_promiscuity.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    main()
