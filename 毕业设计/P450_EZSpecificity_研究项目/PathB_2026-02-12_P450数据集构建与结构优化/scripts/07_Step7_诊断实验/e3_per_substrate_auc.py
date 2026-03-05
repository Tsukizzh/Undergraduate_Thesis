"""
E3: Per-Substrate AUC Analysis
Tests: D1 (diagnostic) — Is 0.52 uniform across all substrates, or do some substrates have high AUC?
Method: Group Step 5 predictions by substrate. Compute AUC for each substrate with >=1 pos + >=5 neg.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
PREDICTIONS_PATH = os.path.join(PATHB, "results", "05_Step5_重构评估", "predictions.csv")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E3_Per_Substrate_AUC")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(PREDICTIONS_PATH)
    print(f"Loaded {len(df)} predictions")

    # Group by substrate
    results_list = []
    for sub_idx, grp in df.groupby("Substrate Index"):
        n_pos = grp["Label"].sum()
        n_neg = len(grp) - n_pos
        if n_pos < 1 or n_neg < 5:
            continue
        try:
            auc = roc_auc_score(grp["Label"], grp["logit"])
        except ValueError:
            continue
        pos_mean = grp[grp["Label"] == 1]["logit"].mean()
        neg_mean = grp[grp["Label"] == 0]["logit"].mean()
        results_list.append({
            "substrate_index": int(sub_idx),
            "n_pos": int(n_pos),
            "n_neg": int(n_neg),
            "auc": float(auc),
            "pos_logit_mean": float(pos_mean),
            "neg_logit_mean": float(neg_mean),
            "score_gap": float(pos_mean - neg_mean),
        })

    res_df = pd.DataFrame(results_list)
    print(f"\nSubstrates with >= 1 pos + >= 5 neg: {len(res_df)}")
    print(f"  AUC median:  {res_df['auc'].median():.4f}")
    print(f"  AUC mean:    {res_df['auc'].mean():.4f}")
    print(f"  AUC std:     {res_df['auc'].std():.4f}")
    print(f"  AUC Q25:     {res_df['auc'].quantile(0.25):.4f}")
    print(f"  AUC Q75:     {res_df['auc'].quantile(0.75):.4f}")
    print(f"  AUC min:     {res_df['auc'].min():.4f}")
    print(f"  AUC max:     {res_df['auc'].max():.4f}")
    print(f"  AUC > 0.65:  {(res_df['auc'] > 0.65).sum()} substrates")
    print(f"  AUC > 0.75:  {(res_df['auc'] > 0.75).sum()} substrates")

    # Verdict
    median_auc = res_df["auc"].median()
    if median_auc > 0.65:
        verdict = "PARTIAL ABILITY: Model can distinguish for some substrates. Global AUC dragged down by failures."
    elif median_auc > 0.55:
        verdict = "WEAK SIGNAL: Slight above-random performance for most substrates."
    else:
        verdict = "UNIFORM FAILURE: Nearly all substrates show random-level AUC. Model uniformly unable to distinguish."

    print(f"\nVerdict: {verdict}")

    # Also compute per-enzyme AUC for comparison
    enzyme_results = []
    for enz_idx, grp in df.groupby("Enzyme Index"):
        n_pos = grp["Label"].sum()
        n_neg = len(grp) - n_pos
        if n_pos < 1 or n_neg < 5:
            continue
        try:
            auc = roc_auc_score(grp["Label"], grp["logit"])
        except ValueError:
            continue
        enzyme_results.append({
            "enzyme_index": int(enz_idx),
            "n_pos": int(n_pos),
            "n_neg": int(n_neg),
            "auc": float(auc),
        })
    enz_df = pd.DataFrame(enzyme_results)
    if len(enz_df) > 0:
        print(f"\nPer-enzyme AUC (>= 1 pos + >= 5 neg): {len(enz_df)} enzymes")
        print(f"  Median: {enz_df['auc'].median():.4f}, Mean: {enz_df['auc'].mean():.4f}")

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Histogram of per-substrate AUC
    ax = axes[0]
    ax.hist(res_df["auc"], bins=20, color="steelblue", edgecolor="white", alpha=0.8)
    ax.axvline(x=0.5, color="red", linestyle="--", label="Random (0.5)")
    ax.axvline(x=median_auc, color="green", linestyle="-", label=f"Median ({median_auc:.3f})")
    ax.set_xlabel("Per-Substrate AUC")
    ax.set_ylabel("Count")
    ax.set_title(f"Per-Substrate AUC Distribution (n={len(res_df)})")
    ax.legend()

    # Sorted AUC with CI-like bars (n_pos shown by dot size)
    ax = axes[1]
    sorted_df = res_df.sort_values("auc")
    ax.scatter(range(len(sorted_df)), sorted_df["auc"],
               s=sorted_df["n_pos"] * 20 + 10, alpha=0.6, c="steelblue", edgecolors="white")
    ax.axhline(y=0.5, color="red", linestyle="--", alpha=0.7)
    ax.set_xlabel("Substrate (sorted by AUC)")
    ax.set_ylabel("AUC")
    ax.set_title("Per-Substrate AUC (dot size = n_pos)")

    # Score gap distribution
    ax = axes[2]
    ax.hist(res_df["score_gap"], bins=20, color="coral", edgecolor="white", alpha=0.8)
    ax.axvline(x=0, color="red", linestyle="--")
    ax.set_xlabel("Score Gap (pos_mean - neg_mean)")
    ax.set_ylabel("Count")
    ax.set_title("Per-Substrate Score Separation")

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e3_per_substrate_auc.png"), dpi=150)
    plt.close()

    # Save
    res_df.to_csv(os.path.join(OUTPUT_DIR, "per_substrate_auc.csv"), index=False)
    if len(enz_df) > 0:
        enz_df.to_csv(os.path.join(OUTPUT_DIR, "per_enzyme_auc.csv"), index=False)

    results = {
        "experiment": "E3",
        "name": "Per-Substrate AUC",
        "tests": "D1 — uniform failure vs partial ability",
        "n_substrates_analyzed": len(res_df),
        "auc_median": float(median_auc),
        "auc_mean": float(res_df["auc"].mean()),
        "auc_std": float(res_df["auc"].std()),
        "auc_q25": float(res_df["auc"].quantile(0.25)),
        "auc_q75": float(res_df["auc"].quantile(0.75)),
        "n_above_065": int((res_df["auc"] > 0.65).sum()),
        "n_above_075": int((res_df["auc"] > 0.75).sum()),
        "verdict": verdict,
        "per_enzyme_median": float(enz_df["auc"].median()) if len(enz_df) > 0 else None,
    }
    with open(os.path.join(OUTPUT_DIR, "e3_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
