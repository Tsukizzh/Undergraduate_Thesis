"""
E6 Analysis: Two-Level Promiscuity-AUC Correlation
Level 1 (Inter-family): 7 families — Spearman(mean_promiscuity, AUC-ROC)
Level 2 (Intra-family): Within each family, bin by promiscuity (log-scale),
                         compute real AUC per bin, Spearman across bins.

IMPORTANT: The `positive_reactions` column in prediction CSVs is WRONG
(mirrors Substrate Index). True promiscuity must be computed from full
family data in tmp_lmdb/{Family}_data.csv.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- Paths ---
PROJECT_ROOT = Path(r"D:\EZSpecificity_Project")
PATHB = PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目" / "PathB_2026-02-12_P450数据集构建与结构优化"
E6_DIR = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E6_expansion_多家族推理"
P450_PRED = PATHB / "results" / "05_Step5_重构评估" / "predictions.csv"
LOCAL_TMP = PROJECT_ROOT / "tmp_lmdb"
OUTPUT_DIR = E6_DIR / "promiscuity_auc_correlation"

FAMILIES_6 = ["Duf", "Esterase", "Gt_acceptor", "Nitrilase", "Phosphatase", "Thiolase"]

# Log-scale bins (powers-of-2 style, half-open: (left, right])
BIN_EDGES = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
BIN_LABELS = ["1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256"]

CONFIDENCE_THRESHOLD = 10  # pos and neg both >= this for high-confidence AUC


def _compute_promiscuity_from_full_data(family: str) -> pd.Series:
    """Compute true promiscuity from full family data (label-based)."""
    full_csv = LOCAL_TMP / f"{family}_data.csv"
    full_df = pd.read_csv(full_csv)
    # columns: reaction, enzyme, label, structure_index, positive_reactions
    prom = (
        full_df.loc[full_df["label"] == 1]
        .groupby("reaction")["enzyme"]
        .nunique()
    )
    return prom  # index = reaction (= Substrate Index), value = count


def _compute_promiscuity_from_labels(df: pd.DataFrame) -> pd.Series:
    """Compute promiscuity from test-set labels (for P450)."""
    prom = (
        df.loc[df["Label"] == 1]
        .groupby("Substrate Index")["Enzyme Index"]
        .nunique()
    )
    return prom  # index = Substrate Index, value = count


def load_family_predictions(family: str) -> pd.DataFrame:
    """Load predictions CSV, add true promiscuity from full family data."""
    path = E6_DIR / f"{family}_predictions.csv"
    df = pd.read_csv(path)

    # Compute true promiscuity
    prom = _compute_promiscuity_from_full_data(family)
    df["promiscuity"] = df["Substrate Index"].map(prom).fillna(0).astype(int)

    # Drop the wrong column to avoid confusion
    if "positive_reactions" in df.columns:
        df = df.drop(columns=["positive_reactions"])

    return df


def load_p450_predictions() -> pd.DataFrame:
    """Load P450 predictions, compute promiscuity from test-set labels."""
    df = pd.read_csv(P450_PRED)

    prom = _compute_promiscuity_from_labels(df)
    df["promiscuity"] = df["Substrate Index"].map(prom).fillna(0).astype(int)

    if "positive_reactions" in df.columns:
        df = df.drop(columns=["positive_reactions"])

    return df


def compute_overall_auc(df: pd.DataFrame) -> float:
    return roc_auc_score(df["Label"], df["logit"])


def assign_promiscuity_bin(promiscuity: int) -> str | None:
    """Assign a promiscuity value to its log-scale bin. Returns None for 0."""
    if promiscuity <= 0:
        return None
    for i in range(len(BIN_EDGES) - 1):
        if BIN_EDGES[i] < promiscuity <= BIN_EDGES[i + 1]:
            return BIN_LABELS[i]
    return None


def compute_binned_auc(df: pd.DataFrame) -> list[dict]:
    """Compute AUC per promiscuity bin. All bins with both classes get AUC."""
    df = df.copy()
    df["prom_bin"] = df["promiscuity"].apply(assign_promiscuity_bin)
    df = df.dropna(subset=["prom_bin"])

    results = []
    for label in BIN_LABELS:
        sub = df[df["prom_bin"] == label]
        if len(sub) == 0:
            continue
        n_pos = int(sub["Label"].sum())
        n_neg = int((sub["Label"] == 0).sum())
        n_substrates = sub["Substrate Index"].nunique()
        prom_median = float(sub.drop_duplicates("Substrate Index")["promiscuity"].median())
        low_conf = n_pos < CONFIDENCE_THRESHOLD or n_neg < CONFIDENCE_THRESHOLD

        if sub["Label"].nunique() < 2:
            results.append({
                "bin": label, "n_pos": n_pos, "n_neg": n_neg,
                "n_substrates": n_substrates, "prom_median": prom_median,
                "auc": None, "low_confidence": low_conf, "note": "single class"
            })
            continue

        auc = float(roc_auc_score(sub["Label"], sub["logit"]))
        results.append({
            "bin": label, "n_pos": n_pos, "n_neg": n_neg,
            "n_substrates": n_substrates, "prom_median": prom_median,
            "auc": auc, "low_confidence": low_conf, "note": None
        })
    return results


def compute_intra_spearman(binned: list[dict]) -> dict:
    """Spearman on HIGH-CONFIDENCE bins only (pos>=10 and neg>=10)."""
    eligible = [b for b in binned
                if b["auc"] is not None and not b["low_confidence"]]
    if len(eligible) < 4:
        return {"rho": None, "p": None, "n_bins": len(eligible),
                "note": f"fewer than 4 high-confidence bins (have {len(eligible)})"}
    x = [b["prom_median"] for b in eligible]
    y = [b["auc"] for b in eligible]
    if len(set(x)) < 2 or len(set(y)) < 2:
        return {"rho": None, "p": None, "n_bins": len(eligible),
                "note": "constant values in bins"}
    rho, p = stats.spearmanr(x, y)
    return {"rho": float(rho), "p": float(p), "n_bins": len(eligible), "note": None}


def analyze_peak_pattern(binned: list[dict]) -> dict:
    """Descriptive peak-bin analysis. NO split Spearman (avoids circular analysis)."""
    valid = [b for b in binned if b["auc"] is not None]
    if len(valid) < 2:
        return {"peak_bin": None, "peak_auc": None, "n_valid_bins": len(valid),
                "pattern": "insufficient_data", "note": "fewer than 2 valid bins"}

    aucs = [b["auc"] for b in valid]
    peak_idx = max(range(len(valid)), key=lambda i: aucs[i])
    peak = valid[peak_idx]

    # Check actual monotonicity by counting direction changes
    increases = sum(1 for i in range(len(aucs) - 1) if aucs[i + 1] > aucs[i])
    decreases = sum(1 for i in range(len(aucs) - 1) if aucs[i + 1] < aucs[i])
    n_steps = len(aucs) - 1

    if n_steps == 0:
        pattern = "single_bin"
    elif increases == n_steps:
        pattern = "monotone_increasing"
    elif decreases == n_steps:
        pattern = "monotone_decreasing"
    elif peak_idx > 0 and peak_idx < len(valid) - 1:
        pattern = "interior_peak"
    else:
        pattern = "mixed"

    return {
        "peak_bin": peak["bin"],
        "peak_auc": float(peak["auc"]),
        "peak_position": f"{peak_idx + 1}/{len(valid)}",
        "n_valid_bins": len(valid),
        "n_increases": increases,
        "n_decreases": decreases,
        "pattern": pattern,
        "note": None,
    }


def plot_inter_family(summary: list[dict], out_path: Path):
    """Scatter plot: mean promiscuity vs AUC across 7 families."""
    fig, ax = plt.subplots(figsize=(8, 6))
    x = [s["promiscuity_mean"] for s in summary]
    y = [s["auc"] for s in summary]
    names = [s["family"] for s in summary]

    ax.scatter(x, y, s=100, zorder=5, edgecolors="black", linewidth=0.8)
    for i, name in enumerate(names):
        ax.annotate(name, (x[i], y[i]), textcoords="offset points",
                    xytext=(8, 5), fontsize=9)

    rho, p = stats.spearmanr(x, y)
    tau, p_tau = stats.kendalltau(x, y)
    ax.text(0.05, 0.95,
            f"Spearman \u03c1 = {rho:.3f} (p = {p:.3f})\n"
            f"Kendall \u03c4 = {tau:.3f} (p = {p_tau:.3f})",
            transform=ax.transAxes, fontsize=10, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    ax.set_xlabel("Mean Promiscuity (# positive enzymes per substrate)", fontsize=11)
    ax.set_ylabel("AUC-ROC", fontsize=11)
    ax.set_title("Inter-Family: Promiscuity vs Model Performance", fontsize=13)
    ax.set_ylim(0.45, 1.0)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_intra_family(all_binned: dict, out_path: Path):
    """Multi-panel: binned AUC vs promiscuity per family."""
    families = [f for f in all_binned if f != "P450" and all_binned[f]["binned"]]
    n = len(families)
    if n == 0:
        print("  No families with valid bins to plot.")
        return
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), squeeze=False)

    for idx, family in enumerate(families):
        ax = axes[idx // cols][idx % cols]
        binned = all_binned[family]["binned"]
        plottable = [b for b in binned if b["auc"] is not None]

        if not plottable:
            ax.text(0.5, 0.5, "No valid bins", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(family)
            continue

        labels_text = [b["bin"] for b in plottable]
        y = [b["auc"] for b in plottable]
        colors = ["lightsteelblue" if b["low_confidence"] else "steelblue"
                  for b in plottable]

        bars = ax.bar(range(len(plottable)), y, color=colors,
                      edgecolor="black", linewidth=0.5)
        for bar, b in zip(bars, plottable):
            bar.set_alpha(0.50 if b["low_confidence"] else 0.75)
            if b["low_confidence"]:
                bar.set_hatch("//")
        ax.set_xticks(range(len(plottable)))
        ax.set_xticklabels(labels_text, fontsize=8, rotation=30)
        ax.set_ylim(0.4, 1.05)
        ax.set_ylabel("AUC-ROC", fontsize=9)
        ax.set_xlabel("Promiscuity bin", fontsize=9)

        sp = all_binned[family]["spearman"]
        title = family
        if sp["rho"] is not None:
            title += f"\n\u03c1={sp['rho']:.3f}, p={sp['p']:.3f}"
        else:
            title += f"\n({sp['note']})"
        ax.set_title(title, fontsize=10)
        ax.grid(True, axis="y", alpha=0.3)

    for idx in range(n, rows * cols):
        axes[idx // cols][idx % cols].set_visible(False)

    fig.suptitle("Intra-Family: Binned AUC vs Promiscuity", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def leave_one_out_sensitivity(summary: list[dict]) -> list[dict]:
    """Leave-one-family-out Spearman sensitivity analysis."""
    results = []
    x_all = [s["promiscuity_mean"] for s in summary]
    y_all = [s["auc"] for s in summary]
    for i, s in enumerate(summary):
        x_loo = x_all[:i] + x_all[i + 1:]
        y_loo = y_all[:i] + y_all[i + 1:]
        rho, p = stats.spearmanr(x_loo, y_loo)
        results.append({"excluded": s["family"], "rho": float(rho), "p": float(p)})
    return results


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("E6 Analysis: Two-Level Promiscuity-AUC Correlation")
    print("=" * 60)

    # --- Load all predictions with CORRECT promiscuity ---
    print("\n[1] Loading predictions (with true promiscuity)...")
    all_preds = {}
    for family in FAMILIES_6:
        all_preds[family] = load_family_predictions(family)
        prom = all_preds[family].drop_duplicates("Substrate Index")["promiscuity"]
        print(f"  {family}: {len(all_preds[family])} pairs, "
              f"promiscuity [{prom.min()}-{prom.max()}] mean={prom.mean():.1f}")
    all_preds["P450"] = load_p450_predictions()
    prom_p450 = all_preds["P450"].drop_duplicates("Substrate Index")["promiscuity"]
    print(f"  P450: {len(all_preds['P450'])} pairs, "
          f"promiscuity [{prom_p450.min()}-{prom_p450.max()}] mean={prom_p450.mean():.2f}")

    # --- Level 1: Inter-family ---
    print("\n[2] Level 1: Inter-family correlation")
    summary = []
    for family, df in all_preds.items():
        auc = compute_overall_auc(df)
        prom_per_sub = df.drop_duplicates("Substrate Index")["promiscuity"]
        # Use only substrates with promiscuity > 0 for mean/median
        # (0-promiscuity substrates are negative-only fillers, not real promiscuity data)
        prom_positive = prom_per_sub[prom_per_sub > 0]
        summary.append({
            "family": family,
            "auc": float(auc),
            "promiscuity_mean": float(prom_positive.mean()) if len(prom_positive) > 0 else 0.0,
            "promiscuity_median": float(prom_positive.median()) if len(prom_positive) > 0 else 0.0,
            "n_pairs": len(df),
            "n_substrates": df["Substrate Index"].nunique(),
            "n_enzymes": df["Enzyme Index"].nunique(),
        })

    summary.sort(key=lambda s: s["auc"])
    print("\n  Family Summary (sorted by AUC):")
    print(f"  {'Family':<15} {'AUC':>7} {'Prom_mean':>10} {'Prom_med':>9} "
          f"{'#pairs':>7} {'#sub':>5} {'#enz':>5}")
    for s in summary:
        print(f"  {s['family']:<15} {s['auc']:>7.4f} {s['promiscuity_mean']:>10.2f} "
              f"{s['promiscuity_median']:>9.1f} {s['n_pairs']:>7} "
              f"{s['n_substrates']:>5} {s['n_enzymes']:>5}")

    x_inter = [s["promiscuity_mean"] for s in summary]
    y_inter = [s["auc"] for s in summary]
    rho_sp, p_sp = stats.spearmanr(x_inter, y_inter)
    tau_k, p_k = stats.kendalltau(x_inter, y_inter)
    x_log = np.log1p(x_inter)
    r_pear, p_pear = stats.pearsonr(x_log, y_inter)

    inter_results = {
        "spearman": {"rho": float(rho_sp), "p": float(p_sp)},
        "kendall": {"tau": float(tau_k), "p": float(p_k)},
        "pearson_log": {"r": float(r_pear), "p": float(p_pear)},
        "summary": summary,
    }

    print(f"\n  Spearman \u03c1 = {rho_sp:.4f} (p = {p_sp:.4f})")
    print(f"  Kendall  \u03c4 = {tau_k:.4f} (p = {p_k:.4f})")
    print(f"  Pearson  r = {r_pear:.4f} (p = {p_pear:.4f}) [on log(1+promiscuity)]")

    loo = leave_one_out_sensitivity(summary)
    print("\n  Leave-one-out sensitivity:")
    for entry in loo:
        print(f"    Exclude {entry['excluded']:<15} -> "
              f"\u03c1 = {entry['rho']:.4f} (p = {entry['p']:.4f})")
    inter_results["leave_one_out"] = loo

    plot_inter_family(summary, OUTPUT_DIR / "level1_inter_family_scatter.png")

    # --- Level 2: Intra-family ---
    print("\n[3] Level 2: Intra-family binned analysis")
    all_binned = {}
    for family, df in all_preds.items():
        print(f"\n  --- {family} ---")
        prom_per_sub = df.drop_duplicates("Substrate Index")["promiscuity"]
        print(f"  Promiscuity range: [{prom_per_sub.min()}, {prom_per_sub.max()}], "
              f"median={prom_per_sub.median()}, mean={prom_per_sub.mean():.1f}")

        if prom_per_sub.nunique() <= 1:
            print(f"  SKIP: no promiscuity variance (all = {prom_per_sub.iloc[0]})")
            all_binned[family] = {
                "binned": [],
                "spearman": {"rho": None, "p": None, "n_bins": 0,
                             "note": "no promiscuity variance"},
                "peak_pattern": analyze_peak_pattern([]),
                "promiscuity_range": [int(prom_per_sub.min()),
                                      int(prom_per_sub.max())]
            }
            continue

        binned = compute_binned_auc(df)
        sp = compute_intra_spearman(binned)
        peak = analyze_peak_pattern(binned)
        all_binned[family] = {
            "binned": binned,
            "spearman": sp,
            "peak_pattern": peak,
            "promiscuity_range": [int(prom_per_sub.min()), int(prom_per_sub.max())]
        }

        for b in binned:
            if b["auc"] is not None:
                flag = " *" if b["low_confidence"] else ""
                status = f"AUC={b['auc']:.4f}{flag}"
            else:
                status = b["note"]
            print(f"  Bin [{b['bin']:>7}]: pos={b['n_pos']:>4}, neg={b['n_neg']:>4}, "
                  f"subs={b['n_substrates']:>3}, prom_med={b['prom_median']:>5.1f}, "
                  f"{status}")
        if sp["rho"] is not None:
            print(f"  Spearman (high-conf bins): \u03c1 = {sp['rho']:.4f} "
                  f"(p = {sp['p']:.4f}), {sp['n_bins']} bins")
        else:
            print(f"  Spearman: N/A ({sp['note']})")
        pb = peak.get("peak_bin") or "N/A"
        pa = peak.get("peak_auc")
        pa_str = f"{pa:.4f}" if pa is not None else "N/A"
        print(f"  Peak pattern: {peak['pattern']} "
              f"(peak at bin '{pb}', AUC={pa_str})")

    plot_intra_family(all_binned, OUTPUT_DIR / "level2_intra_family_binned.png")

    # --- Save results ---
    output = {
        "level1_inter_family": inter_results,
        "level2_intra_family": all_binned,
    }
    out_json = OUTPUT_DIR / "promiscuity_auc_results.json"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False, default=str)
    print(f"\n  Saved: {out_json}")

    rows = []
    for family in all_binned:
        for b in all_binned[family].get("binned", []):
            rows.append({"family": family, **b})
    if rows:
        csv_path = OUTPUT_DIR / "intra_family_binned_auc.csv"
        pd.DataFrame(rows).to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")

    print("\n" + "=" * 60)
    print("Done.")


if __name__ == "__main__":
    main()
