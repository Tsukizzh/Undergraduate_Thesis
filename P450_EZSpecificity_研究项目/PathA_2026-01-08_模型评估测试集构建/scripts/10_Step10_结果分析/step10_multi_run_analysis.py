"""
Multi-Run Analysis Script
=========================
Analyze all inference runs and create comparison summary.

Usage:
    python step10_multi_run_analysis.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

# Path setup
SCRIPT_DIR = Path(__file__).parent.resolve()
PATHA_ROOT = SCRIPT_DIR.parent.parent
DATA_ROOT = PATHA_ROOT / "data"
STEP9_DIR = DATA_ROOT / "09_Step9_模型推理"
STEP10_DIR = DATA_ROOT / "10_Step10_结果分析"

# Find all runs
def find_all_runs():
    """Find all Run directories in Step9."""
    runs = []
    for d in sorted(STEP9_DIR.iterdir()):
        if d.is_dir() and d.name.startswith("Run"):
            pred_csv = d / "predictions.csv"
            if pred_csv.exists():
                runs.append({
                    "run_id": d.name,
                    "predictions_csv": pred_csv,
                    "checkpoint": d.name.split("_", 1)[1] if "_" in d.name else "unknown"
                })
    return runs


def compute_metrics(df):
    """Compute metrics for a single run."""
    from sklearn.metrics import roc_auc_score, average_precision_score

    y_true = df["Label"].values
    y_prob = df["prob"].values
    y_logit = df["logit"].values

    # Check for both classes
    unique_labels = np.unique(y_true)
    has_both = len(unique_labels) >= 2 and 0 in unique_labels and 1 in unique_labels

    if not has_both:
        return {
            "AUC-ROC": np.nan,
            "AUC-PR": np.nan,
            "Logit_Mean": np.mean(y_logit),
            "Logit_Std": np.std(y_logit),
            "Prob_Mean": np.mean(y_prob),
            "N_Samples": len(df),
        }

    return {
        "AUC-ROC": roc_auc_score(y_true, y_prob),
        "AUC-PR": average_precision_score(y_true, y_prob),
        "Logit_Mean": np.mean(y_logit),
        "Logit_Std": np.std(y_logit),
        "Prob_Mean": np.mean(y_prob),
        "N_Samples": len(df),
    }


def main():
    print("="*60)
    print("Multi-Run Analysis")
    print("="*60)

    runs = find_all_runs()
    print(f"\nFound {len(runs)} runs:")
    for r in runs:
        print(f"  - {r['run_id']} ({r['checkpoint']})")

    # Compute metrics for each run
    results = []
    for r in runs:
        print(f"\nAnalyzing {r['run_id']}...")
        df = pd.read_csv(r["predictions_csv"])
        metrics = compute_metrics(df)
        metrics["Run_ID"] = r["run_id"]
        metrics["Checkpoint"] = r["checkpoint"]
        results.append(metrics)

    # Create summary DataFrame
    summary_df = pd.DataFrame(results)
    cols = ["Run_ID", "Checkpoint", "AUC-ROC", "AUC-PR", "Logit_Mean", "Logit_Std", "Prob_Mean", "N_Samples"]
    summary_df = summary_df[cols]

    # Save summary
    output_csv = STEP10_DIR / "multi_run_comparison.csv"
    summary_df.to_csv(output_csv, index=False)
    print(f"\nSaved: {output_csv}")

    # Print summary
    print("\n" + "="*80)
    print("MULTI-RUN COMPARISON SUMMARY")
    print("="*80)
    print(summary_df.to_string(index=False))

    # Highlight best performer
    print("\n" + "-"*40)
    best_idx = summary_df["AUC-ROC"].idxmax()
    best_run = summary_df.loc[best_idx]
    print(f"Best AUC-ROC: {best_run['Run_ID']} ({best_run['Checkpoint']}) = {best_run['AUC-ROC']:.4f}")

    # Create markdown summary
    md_content = f"""# 多模型推理结果对比

## 汇总表

| 运行编号 | 检查点 | AUC-ROC | AUC-PR | Logit均值 | Logit标准差 | Prob均值 | 样本数 |
|---------|--------|---------|--------|-----------|-------------|----------|--------|
"""
    for _, row in summary_df.iterrows():
        md_content += f"| {row['Run_ID']} | {row['Checkpoint']} | {row['AUC-ROC']:.4f} | {row['AUC-PR']:.4f} | {row['Logit_Mean']:.4f} | {row['Logit_Std']:.4f} | {row['Prob_Mean']:.4f} | {row['N_Samples']} |\n"

    md_content += f"""
## 最佳模型

**{best_run['Run_ID']}** ({best_run['Checkpoint']}) 达到最高 AUC-ROC: **{best_run['AUC-ROC']:.4f}**

## 检查点说明

| 检查点 | 说明 |
|--------|------|
| best-checkpoint | 训练最终保存的最佳模型（默认使用） |
| best-checkpoint-v1 | 训练过程中第1次刷新验证集最佳分数时保存 |
| best-checkpoint-v2 | 第2次刷新时保存 |
| best-checkpoint-v3 | 第3次刷新时保存 |
| best-checkpoint-v4 | 第4次刷新时保存 |

v1-v4是训练过程中的中间检查点，可能代表不同训练阶段的模型状态。通常后保存的版本（v4 → best-checkpoint）性能更好，但在特定数据集上可能有差异。

---

*生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}*
"""

    md_path = STEP10_DIR / "多模型推理结果对比.md"
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)
    print(f"Saved: {md_path}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
