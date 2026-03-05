"""
E5: Substrate Base Frequency Check
Tests: M2 confirmation — substrate identity carries no label signal under 100% Direction A.
Method: For each substrate s, compute p(y=1|s) = n_pos / n_total. Use as prediction score. Compute AUC.
Expected: AUC ≈ 0.5 (confirming substrate identity is uninformative).
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
import json
import os

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
PREDICTIONS_PATH = os.path.join(PATHB, "results", "05_Step5_重构评估", "predictions.csv")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E5_底物基频检验")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(PREDICTIONS_PATH)
    print(f"Loaded {len(df)} predictions ({df['Label'].sum()} pos, {(1-df['Label']).sum():.0f} neg)")

    # Compute p(y=1|s) for each substrate
    substrate_stats = df.groupby("Substrate Index").agg(
        n_total=("Label", "count"),
        n_pos=("Label", "sum"),
    ).reset_index()
    substrate_stats["p_pos"] = substrate_stats["n_pos"] / substrate_stats["n_total"]

    print(f"\nSubstrate statistics:")
    print(f"  Total substrates: {len(substrate_stats)}")
    print(f"  p(y=1|s) mean:   {substrate_stats['p_pos'].mean():.4f}")
    print(f"  p(y=1|s) std:    {substrate_stats['p_pos'].std():.4f}")
    print(f"  p(y=1|s) min:    {substrate_stats['p_pos'].min():.4f}")
    print(f"  p(y=1|s) max:    {substrate_stats['p_pos'].max():.4f}")

    # Assign p(y=1|s) as "prediction score" for each sample
    df = df.merge(substrate_stats[["Substrate Index", "p_pos"]], on="Substrate Index", how="left")

    # Compute AUC using p(y=1|s) as prediction
    auc = roc_auc_score(df["Label"], df["p_pos"])
    print(f"\nAUC using p(y=1|s) as score: {auc:.4f}")

    # Bootstrap CI
    rng = np.random.RandomState(42)
    n_bootstrap = 1000
    aucs = []
    for _ in range(n_bootstrap):
        idx = rng.choice(len(df), len(df), replace=True)
        boot_df = df.iloc[idx]
        if boot_df["Label"].nunique() < 2:
            continue
        aucs.append(roc_auc_score(boot_df["Label"], boot_df["p_pos"]))
    ci_lower, ci_upper = np.percentile(aucs, [2.5, 97.5])
    print(f"  95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]")

    # Interpretation
    if abs(auc - 0.5) < 0.02:
        verdict = "CONFIRMED: substrate identity carries no label signal (AUC ≈ 0.5)"
    elif auc > 0.55:
        verdict = "UNEXPECTED: substrate identity carries signal (AUC > 0.55)"
    else:
        verdict = f"MARGINAL: AUC = {auc:.4f}, weak but possibly non-zero signal"

    print(f"\nVerdict: {verdict}")

    # Save results
    results = {
        "experiment": "E5",
        "name": "Substrate Base Frequency Check",
        "tests": "M2 — substrate identity uninformative under 100% Direction A",
        "n_samples": len(df),
        "n_substrates": len(substrate_stats),
        "p_pos_mean": float(substrate_stats["p_pos"].mean()),
        "p_pos_std": float(substrate_stats["p_pos"].std()),
        "p_pos_min": float(substrate_stats["p_pos"].min()),
        "p_pos_max": float(substrate_stats["p_pos"].max()),
        "auc": float(auc),
        "ci_lower": float(ci_lower),
        "ci_upper": float(ci_upper),
        "verdict": verdict,
    }
    with open(os.path.join(OUTPUT_DIR, "e5_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    substrate_stats.to_csv(os.path.join(OUTPUT_DIR, "substrate_frequencies.csv"), index=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
