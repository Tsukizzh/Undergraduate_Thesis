"""
E4: Score Decomposition
Tests: D2 (diagnostic) — Is the prediction just enzyme_bias + substrate_bias (additive model)?
Method: Fit logit ~ mean_enzyme + mean_substrate, compute R². Compute residual AUC.
Uses leave-one-out means to avoid tautological inflation (per Codex suggestion).
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, r2_score
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
PREDICTIONS_PATH = os.path.join(PATHB, "results", "05_Step5_重构评估", "predictions.csv")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E4_分数分解")


def leave_one_out_mean(df, group_col, value_col):
    """Compute leave-one-out group mean for each row."""
    group_sum = df.groupby(group_col)[value_col].transform("sum")
    group_count = df.groupby(group_col)[value_col].transform("count")
    return (group_sum - df[value_col]) / (group_count - 1)


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = pd.read_csv(PREDICTIONS_PATH)
    print(f"Loaded {len(df)} predictions")

    # Standard marginal means
    enzyme_mean = df.groupby("Enzyme Index")["logit"].transform("mean")
    substrate_mean = df.groupby("Substrate Index")["logit"].transform("mean")
    df["enzyme_mean"] = enzyme_mean
    df["substrate_mean"] = substrate_mean

    # Leave-one-out marginal means (avoid self-inclusion bias)
    enzyme_counts = df.groupby("Enzyme Index")["logit"].transform("count")
    substrate_counts = df.groupby("Substrate Index")["logit"].transform("count")
    # Only compute LOO where count > 1
    df["enzyme_mean_loo"] = np.where(
        enzyme_counts > 1,
        leave_one_out_mean(df, "Enzyme Index", "logit"),
        df["enzyme_mean"]
    )
    df["substrate_mean_loo"] = np.where(
        substrate_counts > 1,
        leave_one_out_mean(df, "Substrate Index", "logit"),
        df["substrate_mean"]
    )

    # Additive prediction = enzyme_mean + substrate_mean - grand_mean
    grand_mean = df["logit"].mean()
    df["additive_pred"] = df["enzyme_mean_loo"] + df["substrate_mean_loo"] - grand_mean

    # R² of additive model
    r2_additive = r2_score(df["logit"], df["additive_pred"])
    # Also standard R² (for comparison)
    df["additive_pred_std"] = df["enzyme_mean"] + df["substrate_mean"] - grand_mean
    r2_standard = r2_score(df["logit"], df["additive_pred_std"])

    # R² for enzyme-only and substrate-only models
    r2_enzyme = r2_score(df["logit"], df["enzyme_mean_loo"])
    r2_substrate = r2_score(df["logit"], df["substrate_mean_loo"])

    # Residuals and residual AUC
    df["residual"] = df["logit"] - df["additive_pred"]
    residual_auc = roc_auc_score(df["Label"], df["residual"])
    original_auc = roc_auc_score(df["Label"], df["logit"])

    # AUC using only enzyme_mean and substrate_mean
    enzyme_auc = roc_auc_score(df["Label"], df["enzyme_mean_loo"])
    substrate_auc = roc_auc_score(df["Label"], df["substrate_mean_loo"])
    additive_auc = roc_auc_score(df["Label"], df["additive_pred"])

    print(f"\n--- R-squared Analysis ---")
    print(f"R2 (enzyme-only LOO):       {r2_enzyme:.4f}")
    print(f"R2 (substrate-only LOO):    {r2_substrate:.4f}")
    print(f"R2 (additive LOO):          {r2_additive:.4f}")
    print(f"R2 (additive standard):     {r2_standard:.4f}")

    print(f"\n--- AUC Analysis ---")
    print(f"Original AUC (logit):       {original_auc:.4f}")
    print(f"Enzyme-mean AUC:            {enzyme_auc:.4f}")
    print(f"Substrate-mean AUC:         {substrate_auc:.4f}")
    print(f"Additive model AUC:         {additive_auc:.4f}")
    print(f"Residual AUC:               {residual_auc:.4f}")

    # Interpretation
    if r2_additive > 0.8:
        verdict = ("MARGINAL DOMINANCE: Predictions are almost entirely explained by enzyme and "
                   "substrate marginal effects. No meaningful enzyme-substrate interaction learned.")
    elif r2_additive > 0.5:
        verdict = ("PARTIAL INTERACTION: About half the variance comes from marginal effects. "
                   "Some enzyme-substrate interaction signal exists.")
    else:
        verdict = ("INTERACTION PRESENT: Marginal effects explain < 50% of variance. "
                   "Model may capture some enzyme-substrate specificity.")

    print(f"\nVerdict: {verdict}")

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Actual vs predicted (additive model)
    ax = axes[0]
    ax.scatter(df["additive_pred"], df["logit"], s=5, alpha=0.3, c=df["Label"].map({0: "blue", 1: "red"}))
    min_val = min(df["additive_pred"].min(), df["logit"].min())
    max_val = max(df["additive_pred"].max(), df["logit"].max())
    ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.5)
    ax.set_xlabel("Additive Prediction (enzyme + substrate mean)")
    ax.set_ylabel("Actual Logit")
    ax.set_title(f"Additive Model R²={r2_additive:.3f}")

    # R² bar chart
    ax = axes[1]
    categories = ["Enzyme\nonly", "Substrate\nonly", "Additive\n(both)"]
    r2_values = [r2_enzyme, r2_substrate, r2_additive]
    bars = ax.bar(categories, r2_values, color=["steelblue", "coral", "green"], alpha=0.7)
    for bar, val in zip(bars, r2_values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01, f"{val:.3f}",
                ha="center", fontsize=10)
    ax.set_ylabel("R²")
    ax.set_title("Variance Explained by Marginal Effects")
    ax.set_ylim(0, 1)

    # Residual distribution by label
    ax = axes[2]
    ax.hist(df[df["Label"] == 1]["residual"], bins=30, alpha=0.6, density=True,
            label=f"Positive (n={df['Label'].sum()})", color="green")
    ax.hist(df[df["Label"] == 0]["residual"], bins=30, alpha=0.6, density=True,
            label=f"Negative (n={(1-df['Label']).sum():.0f})", color="red")
    ax.axvline(x=0, color="gray", linestyle="--")
    ax.set_xlabel("Residual (logit - additive_pred)")
    ax.set_ylabel("Density")
    ax.set_title(f"Residual AUC = {residual_auc:.4f}")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e4_score_decomposition.png"), dpi=150)
    plt.close()

    # Save
    results = {
        "experiment": "E4",
        "name": "Score Decomposition",
        "tests": "D2 — additive bias vs interaction",
        "r2_enzyme_loo": float(r2_enzyme),
        "r2_substrate_loo": float(r2_substrate),
        "r2_additive_loo": float(r2_additive),
        "r2_additive_standard": float(r2_standard),
        "auc_original": float(original_auc),
        "auc_enzyme_mean": float(enzyme_auc),
        "auc_substrate_mean": float(substrate_auc),
        "auc_additive": float(additive_auc),
        "auc_residual": float(residual_auc),
        "verdict": verdict,
    }
    with open(os.path.join(OUTPUT_DIR, "e4_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
