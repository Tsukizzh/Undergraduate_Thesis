"""
Batch Analysis for All Runs
===========================
Generate full Step10 analysis for each Run directory.

Usage:
    python step10_batch_full_analysis.py
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score,
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, matthews_corrcoef, brier_score_loss,
    roc_auc_score
)

# Path setup
SCRIPT_DIR = Path(__file__).parent.resolve()
PATHA_ROOT = SCRIPT_DIR.parent.parent
DATA_ROOT = PATHA_ROOT / "data"
STEP9_DIR = DATA_ROOT / "09_Step9_模型推理"
STEP10_DIR = DATA_ROOT / "10_Step10_结果分析"

# Colors
COLOR_POSITIVE = '#2ecc71'
COLOR_NEGATIVE = '#e74c3c'


def has_both_classes(y_true):
    unique = np.unique(y_true)
    return len(unique) >= 2 and 0 in unique and 1 in unique


def safe_roc_auc(y_true, y_prob):
    if not has_both_classes(y_true):
        return np.nan
    return roc_auc_score(y_true, y_prob)


def safe_average_precision(y_true, y_prob):
    if not has_both_classes(y_true):
        return np.nan
    return average_precision_score(y_true, y_prob)


def compute_global_metrics(y_true, y_prob, y_logit, threshold=0.5):
    y_pred = (y_prob >= threshold).astype(int)
    both_classes = has_both_classes(y_true)

    return {
        'AUC-ROC': safe_roc_auc(y_true, y_prob),
        'AUC-PR (AUPR)': safe_average_precision(y_true, y_prob),
        'Accuracy': accuracy_score(y_true, y_pred),
        'Precision': precision_score(y_true, y_pred, zero_division=0),
        'Recall (Sensitivity)': recall_score(y_true, y_pred, zero_division=0),
        'Specificity': recall_score(1 - y_true, 1 - y_pred, zero_division=0) if both_classes else np.nan,
        'F1 Score': f1_score(y_true, y_pred, zero_division=0),
        'MCC': matthews_corrcoef(y_true, y_pred) if both_classes else np.nan,
        'Brier Score': brier_score_loss(y_true, y_prob),
        'Logit Mean': np.mean(y_logit),
        'Logit Std': np.std(y_logit),
        'Logit Median': np.median(y_logit),
        'Prob Mean': np.mean(y_prob),
        'Prob Std': np.std(y_prob),
    }


def find_optimal_threshold(y_true, y_prob):
    if not has_both_classes(y_true):
        return (0.5, np.nan)
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    j_scores = tpr - fpr
    optimal_idx = np.argmax(j_scores)
    return thresholds[optimal_idx], j_scores[optimal_idx]


def create_main_figure(y_true, y_prob, y_logit, run_id):
    fig = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
    both_classes = has_both_classes(y_true)

    # Panel A: ROC Curve
    ax1 = fig.add_subplot(gs[0, 0])
    if both_classes:
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        roc_auc = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=COLOR_POSITIVE, lw=2, label=f'ROC (AUC = {roc_auc:.3f})')
        ax1.fill_between(fpr, tpr, alpha=0.2, color=COLOR_POSITIVE)
    ax1.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Random')
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1.02])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('A. ROC Curve', fontweight='bold', loc='left')
    ax1.legend(loc='lower right')
    ax1.set_aspect('equal')

    # Panel B: PR Curve
    ax2 = fig.add_subplot(gs[0, 1])
    baseline = y_true.sum() / len(y_true)
    if both_classes:
        precision, recall, _ = precision_recall_curve(y_true, y_prob)
        ap = average_precision_score(y_true, y_prob)
        ax2.plot(recall, precision, color=COLOR_NEGATIVE, lw=2, label=f'PR (AP = {ap:.3f})')
        ax2.fill_between(recall, precision, alpha=0.2, color=COLOR_NEGATIVE)
    ax2.axhline(y=baseline, color='k', linestyle='--', lw=1, alpha=0.5, label=f'Baseline ({baseline:.3f})')
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1.02])
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('B. Precision-Recall Curve', fontweight='bold', loc='left')
    ax2.legend(loc='lower left')
    ax2.set_aspect('equal')

    # Panel C: Score distribution histogram
    ax3 = fig.add_subplot(gs[0, 2])
    pos_logits = y_logit[y_true == 1]
    neg_logits = y_logit[y_true == 0]
    bins = np.linspace(min(y_logit.min(), -10), max(y_logit.max(), 5), 40)
    ax3.hist(pos_logits, bins=bins, alpha=0.6, color=COLOR_POSITIVE, label=f'Positive (n={len(pos_logits)})', density=True)
    ax3.hist(neg_logits, bins=bins, alpha=0.6, color=COLOR_NEGATIVE, label=f'Negative (n={len(neg_logits)})', density=True)
    ax3.axvline(x=0, color='k', linestyle='--', lw=1, alpha=0.5)
    ax3.set_xlabel('Logit Score')
    ax3.set_ylabel('Density')
    ax3.set_title('C. Score Distribution', fontweight='bold', loc='left')
    ax3.legend(loc='upper right')

    # Panel D: Prob distribution
    ax4 = fig.add_subplot(gs[1, 0])
    pos_probs = y_prob[y_true == 1]
    neg_probs = y_prob[y_true == 0]
    bins_prob = np.linspace(0, 1, 30)
    ax4.hist(pos_probs, bins=bins_prob, alpha=0.6, color=COLOR_POSITIVE, label='Positive', density=True)
    ax4.hist(neg_probs, bins=bins_prob, alpha=0.6, color=COLOR_NEGATIVE, label='Negative', density=True)
    ax4.axvline(x=0.5, color='k', linestyle='--', lw=1, alpha=0.5, label='Threshold=0.5')
    ax4.set_xlabel('Predicted Probability')
    ax4.set_ylabel('Density')
    ax4.set_title('D. Probability Distribution', fontweight='bold', loc='left')
    ax4.legend(loc='upper center')

    # Panel E: Confusion Matrix
    ax5 = fig.add_subplot(gs[1, 1])
    y_pred = (y_prob >= 0.5).astype(int)
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
    row_sums = cm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    cm_norm = cm.astype('float') / row_sums
    im = ax5.imshow(cm_norm, interpolation='nearest', cmap='Blues', vmin=0, vmax=1)
    for i in range(2):
        for j in range(2):
            text = f'{cm[i, j]}\n({cm_norm[i, j]:.1%})'
            ax5.text(j, i, text, ha='center', va='center', fontsize=11)
    ax5.set_xticks([0, 1])
    ax5.set_yticks([0, 1])
    ax5.set_xticklabels(['Pred Neg', 'Pred Pos'])
    ax5.set_yticklabels(['True Neg', 'True Pos'])
    ax5.set_title('E. Confusion Matrix', fontweight='bold', loc='left')

    # Panel F: Logit boxplot
    ax6 = fig.add_subplot(gs[1, 2])
    boxplot_data = [neg_logits, pos_logits]
    bp = ax6.boxplot(boxplot_data, tick_labels=['Negative', 'Positive'], patch_artist=True, widths=0.6)
    for box, color in zip(bp['boxes'], [COLOR_NEGATIVE, COLOR_POSITIVE]):
        box.set_facecolor(color)
        box.set_alpha(0.6)
    ax6.axhline(y=0, color='k', linestyle='--', lw=1, alpha=0.5)
    ax6.set_ylabel('Logit Score')
    ax6.set_title('F. Logit Distribution by Class', fontweight='bold', loc='left')

    fig.suptitle(f'EZSpecificity Model Evaluation - {run_id}', fontsize=14, fontweight='bold', y=0.98)
    return fig


def analyze_single_run(run_dir, output_dir):
    """Analyze a single run and save results to output_dir."""
    run_id = run_dir.name
    pred_csv = run_dir / "predictions.csv"

    print(f"  Loading {pred_csv}...")
    df = pd.read_csv(pred_csv)

    y_true = df["Label"].values
    y_prob = df["prob"].values
    y_logit = df["logit"].values

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)

    # Compute metrics
    print(f"  Computing metrics...")
    metrics = compute_global_metrics(y_true, y_prob, y_logit)
    optimal_threshold, youden_j = find_optimal_threshold(y_true, y_prob)

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(output_dir / "metrics_summary.csv", index=False)

    # Create main figure
    print(f"  Creating figures...")
    fig = create_main_figure(y_true, y_prob, y_logit, run_id)
    fig.savefig(figures_dir / "main_analysis_figure.png", dpi=300, bbox_inches='tight')
    fig.savefig(figures_dir / "main_analysis_figure.pdf", bbox_inches='tight')
    plt.close(fig)

    # Error analysis
    y_pred = (y_prob >= 0.5).astype(int)
    fp_mask = (y_pred == 1) & (y_true == 0)
    fn_mask = (y_pred == 0) & (y_true == 1)
    error_df = df.copy()
    error_df['Predicted'] = y_pred
    error_df['Error_Type'] = 'Correct'
    error_df.loc[fp_mask, 'Error_Type'] = 'False Positive'
    error_df.loc[fn_mask, 'Error_Type'] = 'False Negative'
    error_df[error_df['Error_Type'] != 'Correct'].to_csv(output_dir / "error_analysis.csv", index=False)

    # Analysis report
    report_lines = [
        "=" * 70,
        f"EZSpecificity Model Evaluation Report - {run_id}",
        "=" * 70,
        "",
        "1. GLOBAL PERFORMANCE METRICS",
        "-" * 40,
    ]
    for key, value in metrics.items():
        if isinstance(value, float):
            report_lines.append(f"  {key}: {value:.4f}" if not np.isnan(value) else f"  {key}: N/A")
        else:
            report_lines.append(f"  {key}: {value}")
    report_lines.extend([
        "",
        "2. THRESHOLD ANALYSIS",
        "-" * 40,
        f"  Optimal threshold (Youden's J): {optimal_threshold:.4f}",
        f"  Youden's J statistic: {youden_j:.4f}" if not np.isnan(youden_j) else "  Youden's J statistic: N/A",
        "",
        "3. ERROR ANALYSIS (threshold=0.5)",
        "-" * 40,
        f"  False Positives: {fp_mask.sum()}",
        f"  False Negatives: {fn_mask.sum()}",
        "",
        "=" * 70,
    ])

    with open(output_dir / "analysis_report.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    # README
    checkpoint = run_id.split("_", 1)[1] if "_" in run_id else "unknown"
    readme = f"""# {run_id} 分析结果

## 检查点
`{checkpoint}.ckpt`

## 核心指标

| 指标 | 值 |
|------|-----|
| AUC-ROC | {metrics['AUC-ROC']:.4f} |
| AUC-PR | {metrics['AUC-PR (AUPR)']:.4f} |
| Accuracy | {metrics['Accuracy']:.4f} |
| Recall | {metrics['Recall (Sensitivity)']:.4f} |
| F1 Score | {metrics['F1 Score']:.4f} |

## 输出文件

- `metrics_summary.csv`: 指标汇总
- `error_analysis.csv`: 错误样本详情
- `analysis_report.txt`: 文本报告
- `figures/main_analysis_figure.png`: 分析图
"""
    with open(output_dir / "README.md", "w", encoding="utf-8") as f:
        f.write(readme)

    return metrics


def main():
    print("=" * 60)
    print("Batch Full Analysis for All Runs")
    print("=" * 60)

    # Find all runs
    runs = sorted([d for d in STEP9_DIR.iterdir() if d.is_dir() and d.name.startswith("Run")])
    print(f"\nFound {len(runs)} runs to analyze")

    for run_dir in runs:
        run_id = run_dir.name
        output_dir = STEP10_DIR / run_id

        print(f"\nAnalyzing {run_id}...")
        try:
            analyze_single_run(run_dir, output_dir)
            print(f"  Done: {output_dir}")
        except Exception as e:
            print(f"  ERROR: {e}")

    print("\n" + "=" * 60)
    print("Batch analysis complete!")


if __name__ == "__main__":
    main()
