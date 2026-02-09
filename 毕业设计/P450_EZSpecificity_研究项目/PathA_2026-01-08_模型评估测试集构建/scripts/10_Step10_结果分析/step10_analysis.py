"""
Step 10: EZSpecificity P450 Model Evaluation Analysis
======================================================
Comprehensive analysis of model predictions with publication-quality visualizations.

Based on Codex analysis plan and Gemini visualization design.

Usage:
    D:/anaconda3/envs/torch/python.exe step10_analysis.py

Output:
    data/10_Step10_结果分析/
    ├── figures/
    │   ├── main_analysis_figure.png
    │   ├── main_analysis_figure.pdf
    │   ├── per_enzyme_analysis.png
    │   └── score_distribution_detail.png
    ├── metrics_summary.csv
    ├── per_enzyme_metrics.csv
    ├── per_substrate_metrics.csv
    ├── error_analysis.csv
    └── analysis_report.txt
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
import seaborn as sns
from scipy import stats

# Scikit-learn metrics
from sklearn.metrics import (
    roc_curve, auc, precision_recall_curve, average_precision_score,
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, matthews_corrcoef, brier_score_loss,
    roc_auc_score, classification_report
)
from sklearn.calibration import calibration_curve

# Only suppress specific non-critical warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")


# ============================================================================
# Safe Metric Functions (handle edge cases)
# ============================================================================
def has_both_classes(y_true):
    """Check if both classes (0 and 1) are present in y_true."""
    unique = np.unique(y_true)
    return len(unique) >= 2 and 0 in unique and 1 in unique


def safe_roc_auc(y_true, y_prob):
    """Compute ROC-AUC with single-class handling."""
    if not has_both_classes(y_true):
        return np.nan
    return roc_auc_score(y_true, y_prob)


def safe_average_precision(y_true, y_prob):
    """Compute Average Precision with single-class handling."""
    if not has_both_classes(y_true):
        return np.nan
    return average_precision_score(y_true, y_prob)

# ============================================================================
# Path Configuration
# ============================================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
PATHA_ROOT = SCRIPT_DIR.parent.parent  # PathA_2026-01-08_模型评估测试集构建

# Input files
DATA_ROOT = PATHA_ROOT / "data"
PREDICTIONS_CSV = DATA_ROOT / "09_Step9_模型推理" / "predictions.csv"
ENZYMES_CSV = DATA_ROOT / "04_Step4_格式修正后数据" / "Enzymes.csv"
SUBSTRATES_CSV = DATA_ROOT / "04_Step4_格式修正后数据" / "Substrates.csv"

# Output directories
OUTPUT_DIR = DATA_ROOT / "10_Step10_结果分析"
FIGURES_DIR = OUTPUT_DIR / "figures"

# ============================================================================
# Visualization Settings (Gemini's design)
# ============================================================================
# Okabe-Ito colorblind-friendly palette
COLOR_POSITIVE = "#D55E00"  # Vermilion for positive class
COLOR_NEGATIVE = "#0072B2"  # Blue for negative class
COLOR_NEUTRAL = "#009E73"   # Bluish green for neutral/threshold
COLOR_HIGHLIGHT = "#F0E442" # Yellow for highlights

# Matplotlib publication settings
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica'],
    'font.size': 10,
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 100,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 1.0,
    'axes.grid': True,
    'grid.alpha': 0.3,
})

# Set style
plt.style.use('seaborn-v0_8-whitegrid')


def load_data():
    """Load prediction results and metadata."""
    print("Loading data...")

    pred_df = pd.read_csv(PREDICTIONS_CSV)
    enzymes_df = pd.read_csv(ENZYMES_CSV)
    substrates_df = pd.read_csv(SUBSTRATES_CSV)

    print(f"  Predictions: {len(pred_df)} samples")
    print(f"  Enzymes: {len(enzymes_df)} unique")
    print(f"  Substrates: {len(substrates_df)} unique")

    return pred_df, enzymes_df, substrates_df


def compute_global_metrics(y_true, y_prob, y_logit, threshold=0.5):
    """Compute comprehensive performance metrics."""
    if len(y_true) == 0:
        return {'Error': 'Empty dataset'}

    y_pred = (y_prob >= threshold).astype(int)
    both_classes = has_both_classes(y_true)

    metrics = {
        # Threshold-independent
        'AUC-ROC': safe_roc_auc(y_true, y_prob),
        'AUC-PR (AUPR)': safe_average_precision(y_true, y_prob),

        # Threshold-dependent (default 0.5)
        'Accuracy': accuracy_score(y_true, y_pred),
        'Precision': precision_score(y_true, y_pred, zero_division=0),
        'Recall (Sensitivity)': recall_score(y_true, y_pred, zero_division=0),
        'Specificity': recall_score(1 - y_true, 1 - y_pred, zero_division=0) if both_classes else np.nan,
        'F1 Score': f1_score(y_true, y_pred, zero_division=0),
        'MCC': matthews_corrcoef(y_true, y_pred) if both_classes else np.nan,

        # Calibration
        'Brier Score': brier_score_loss(y_true, y_prob),

        # Score statistics
        'Logit Mean': np.mean(y_logit),
        'Logit Std': np.std(y_logit),
        'Logit Median': np.median(y_logit),
        'Prob Mean': np.mean(y_prob),
        'Prob Std': np.std(y_prob),
    }

    return metrics


def find_optimal_threshold(y_true, y_prob):
    """Find optimal threshold using Youden's J statistic."""
    if not has_both_classes(y_true):
        return (0.5, np.nan)  # Default threshold
    fpr, tpr, thresholds = roc_curve(y_true, y_prob)
    j_scores = tpr - fpr  # Youden's J
    optimal_idx = np.argmax(j_scores)
    return thresholds[optimal_idx], j_scores[optimal_idx]


def compute_per_enzyme_metrics(pred_df, enzymes_df):
    """Compute metrics per enzyme."""
    results = []

    for enz_idx in pred_df['Enzyme Index'].unique():
        subset = pred_df[pred_df['Enzyme Index'] == enz_idx]

        if len(subset) < 2:
            continue

        y_true = subset['Label'].values
        y_prob = subset['prob'].values

        # Skip if only one class present
        if len(np.unique(y_true)) < 2:
            auc_roc = np.nan
        else:
            auc_roc = roc_auc_score(y_true, y_prob)

        # Get enzyme metadata
        enz_info = enzymes_df[enzymes_df['Enzyme_Index'] == enz_idx]
        uniprot_id = enz_info['UniProt_ID'].values[0] if len(enz_info) > 0 else 'Unknown'
        organism = enz_info['Organism'].values[0] if len(enz_info) > 0 and 'Organism' in enz_info.columns else 'Unknown'

        results.append({
            'Enzyme_Index': enz_idx,
            'UniProt_ID': uniprot_id,
            'Organism': organism,
            'N_samples': len(subset),
            'N_positive': int(y_true.sum()),
            'N_negative': int(len(y_true) - y_true.sum()),
            'AUC-ROC': auc_roc,
            'Mean_prob_positive': subset[subset['Label'] == 1]['prob'].mean() if (subset['Label'] == 1).any() else np.nan,
            'Mean_prob_negative': subset[subset['Label'] == 0]['prob'].mean() if (subset['Label'] == 0).any() else np.nan,
        })

    return pd.DataFrame(results)


def compute_per_substrate_metrics(pred_df, substrates_df):
    """Compute metrics per substrate."""
    results = []

    for sub_idx in pred_df['Substrate Index'].unique():
        subset = pred_df[pred_df['Substrate Index'] == sub_idx]

        if len(subset) < 2:
            continue

        y_true = subset['Label'].values
        y_prob = subset['prob'].values

        # Skip if only one class present
        if len(np.unique(y_true)) < 2:
            auc_roc = np.nan
        else:
            auc_roc = roc_auc_score(y_true, y_prob)

        # Get substrate SMILES
        sub_info = substrates_df[substrates_df['Substrate_Index'] == sub_idx]
        smiles = sub_info['Substrate_SMILES'].values[0] if len(sub_info) > 0 else 'Unknown'

        results.append({
            'Substrate_Index': sub_idx,
            'SMILES': smiles[:50] + '...' if len(smiles) > 50 else smiles,
            'N_samples': len(subset),
            'N_positive': int(y_true.sum()),
            'N_negative': int(len(y_true) - y_true.sum()),
            'AUC-ROC': auc_roc,
            'Mean_prob_positive': subset[subset['Label'] == 1]['prob'].mean() if (subset['Label'] == 1).any() else np.nan,
            'Mean_prob_negative': subset[subset['Label'] == 0]['prob'].mean() if (subset['Label'] == 0).any() else np.nan,
        })

    return pd.DataFrame(results)


def analyze_errors(pred_df, threshold=0.5):
    """Analyze false positives and false negatives."""
    if len(pred_df) == 0:
        print("\nError Analysis: Empty dataset, skipping.")
        return pd.DataFrame()

    y_true = pred_df['Label'].values
    y_prob = pred_df['prob'].values
    y_pred = (y_prob >= threshold).astype(int)

    # Identify errors
    fp_mask = (y_pred == 1) & (y_true == 0)  # False positive
    fn_mask = (y_pred == 0) & (y_true == 1)  # False negative

    error_df = pred_df.copy()
    error_df['Predicted'] = y_pred
    error_df['Error_Type'] = 'Correct'
    error_df.loc[fp_mask, 'Error_Type'] = 'False Positive'
    error_df.loc[fn_mask, 'Error_Type'] = 'False Negative'

    # Summary
    fp_count = fp_mask.sum()
    fn_count = fn_mask.sum()
    n_total = len(pred_df)

    print(f"\nError Analysis (threshold={threshold}):")
    print(f"  False Positives: {fp_count} ({100*fp_count/n_total:.1f}%)")
    print(f"  False Negatives: {fn_count} ({100*fn_count/n_total:.1f}%)")

    return error_df[error_df['Error_Type'] != 'Correct']


def create_main_figure(y_true, y_prob, y_logit):
    """Create the main 2x3 analysis figure (Gemini's design)."""
    fig = plt.figure(figsize=(14, 9))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

    both_classes = has_both_classes(y_true)

    # =========== Panel A: ROC Curve ===========
    ax1 = fig.add_subplot(gs[0, 0])
    if both_classes:
        fpr, tpr, _ = roc_curve(y_true, y_prob)
        roc_auc = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=COLOR_POSITIVE, lw=2, label=f'ROC (AUC = {roc_auc:.3f})')
        ax1.fill_between(fpr, tpr, alpha=0.2, color=COLOR_POSITIVE)
    else:
        ax1.text(0.5, 0.5, 'Single class\n(ROC undefined)', ha='center', va='center', fontsize=12)
    ax1.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Random')
    ax1.set_xlim([0, 1])
    ax1.set_ylim([0, 1.02])
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('A. ROC Curve', fontweight='bold', loc='left')
    ax1.legend(loc='lower right')
    ax1.set_aspect('equal')

    # =========== Panel B: Precision-Recall Curve ===========
    ax2 = fig.add_subplot(gs[0, 1])
    baseline = y_true.sum() / len(y_true) if len(y_true) > 0 else 0.5
    if both_classes:
        precision, recall, _ = precision_recall_curve(y_true, y_prob)
        ap = average_precision_score(y_true, y_prob)
        ax2.plot(recall, precision, color=COLOR_NEGATIVE, lw=2, label=f'PR (AP = {ap:.3f})')
        ax2.fill_between(recall, precision, alpha=0.2, color=COLOR_NEGATIVE)
    else:
        ax2.text(0.5, 0.5, 'Single class\n(PR undefined)', ha='center', va='center', fontsize=12)
    ax2.axhline(y=baseline, color='k', linestyle='--', lw=1, alpha=0.5, label=f'Baseline ({baseline:.3f})')
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1.02])
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('B. Precision-Recall Curve', fontweight='bold', loc='left')
    ax2.legend(loc='lower left')
    ax2.set_aspect('equal')

    # =========== Panel C: Calibration Curve ===========
    ax3 = fig.add_subplot(gs[0, 2])
    prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=10, strategy='uniform')
    brier = brier_score_loss(y_true, y_prob)

    ax3.plot(prob_pred, prob_true, 's-', color=COLOR_NEUTRAL, lw=2, markersize=6,
             label=f'Model (Brier = {brier:.3f})')
    ax3.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Perfectly calibrated')
    ax3.set_xlim([0, 1])
    ax3.set_ylim([0, 1.02])
    ax3.set_xlabel('Mean Predicted Probability')
    ax3.set_ylabel('Fraction of Positives')
    ax3.set_title('C. Calibration Curve', fontweight='bold', loc='left')
    ax3.legend(loc='lower right')
    ax3.set_aspect('equal')

    # =========== Panel D: Score Distribution ===========
    ax4 = fig.add_subplot(gs[1, 0])

    # Separate by label
    pos_probs = y_prob[y_true == 1]
    neg_probs = y_prob[y_true == 0]

    # KDE plots
    if len(pos_probs) > 1:
        sns.kdeplot(pos_probs, ax=ax4, color=COLOR_POSITIVE, fill=True, alpha=0.4,
                    label=f'Positive (n={len(pos_probs)})', linewidth=2)
    if len(neg_probs) > 1:
        sns.kdeplot(neg_probs, ax=ax4, color=COLOR_NEGATIVE, fill=True, alpha=0.4,
                    label=f'Negative (n={len(neg_probs)})', linewidth=2)

    ax4.axvline(x=0.5, color='k', linestyle='--', lw=1, alpha=0.7, label='Threshold (0.5)')
    ax4.set_xlim([0, 1])
    ax4.set_xlabel('Predicted Probability')
    ax4.set_ylabel('Density')
    ax4.set_title('D. Score Distribution', fontweight='bold', loc='left')
    ax4.legend(loc='upper center')

    # =========== Panel E: Confusion Matrix ===========
    ax5 = fig.add_subplot(gs[1, 1])
    y_pred = (y_prob >= 0.5).astype(int)
    # Use labels=[0,1] to ensure 2x2 matrix even with single class
    cm = confusion_matrix(y_true, y_pred, labels=[0, 1])

    # Normalize for display (handle division by zero)
    row_sums = cm.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    cm_norm = cm.astype('float') / row_sums

    # Plot heatmap
    im = ax5.imshow(cm_norm, interpolation='nearest', cmap='Blues', vmin=0, vmax=1)

    # Add text annotations
    for i in range(2):
        for j in range(2):
            text = f'{cm[i, j]}\n({cm_norm[i, j]:.1%})'
            color = 'white' if cm_norm[i, j] > 0.5 else 'black'
            ax5.text(j, i, text, ha='center', va='center', fontsize=11, color=color)

    ax5.set_xticks([0, 1])
    ax5.set_yticks([0, 1])
    ax5.set_xticklabels(['Negative', 'Positive'])
    ax5.set_yticklabels(['Negative', 'Positive'])
    ax5.set_xlabel('Predicted')
    ax5.set_ylabel('Actual')
    ax5.set_title('E. Confusion Matrix', fontweight='bold', loc='left')

    # =========== Panel F: Logit Distribution ===========
    ax6 = fig.add_subplot(gs[1, 2])

    pos_logits = y_logit[y_true == 1]
    neg_logits = y_logit[y_true == 0]

    # Box plot - handle empty arrays
    boxplot_data = []
    boxplot_labels = []
    box_colors = []

    if len(neg_logits) > 0:
        boxplot_data.append(neg_logits)
        boxplot_labels.append('Negative')
        box_colors.append(COLOR_NEGATIVE)
    if len(pos_logits) > 0:
        boxplot_data.append(pos_logits)
        boxplot_labels.append('Positive')
        box_colors.append(COLOR_POSITIVE)

    if len(boxplot_data) > 0:
        bp = ax6.boxplot(boxplot_data, tick_labels=boxplot_labels, patch_artist=True, widths=0.6)
        for box, color in zip(bp['boxes'], box_colors):
            box.set_facecolor(color)
            box.set_alpha(0.6)
    else:
        ax6.text(0.5, 0.5, 'No data', ha='center', va='center', fontsize=12)

    ax6.axhline(y=0, color='k', linestyle='--', lw=1, alpha=0.5)
    ax6.set_ylabel('Logit Score')
    ax6.set_title('F. Logit Distribution by Class', fontweight='bold', loc='left')

    # Add statistical annotation only if both classes have data
    if len(pos_logits) > 0 and len(neg_logits) > 0:
        stat, pval = stats.mannwhitneyu(pos_logits, neg_logits, alternative='greater')
        ax6.annotate(f'Mann-Whitney U\np < 0.001' if pval < 0.001 else f'p = {pval:.3f}',
                     xy=(0.5, 0.95), xycoords='axes fraction',
                     ha='center', va='top', fontsize=9,
                     bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Add overall title
    fig.suptitle('EZSpecificity Model Evaluation on P450 Test Set',
                 fontsize=14, fontweight='bold', y=0.98)

    return fig


def create_per_enzyme_figure(per_enzyme_df):
    """Create per-enzyme analysis figure."""
    # Handle empty DataFrame
    if len(per_enzyme_df) == 0 or 'AUC-ROC' not in per_enzyme_df.columns:
        print("No per-enzyme data available for figure.")
        return None

    # Filter enzymes with valid AUC
    valid_df = per_enzyme_df.dropna(subset=['AUC-ROC']).sort_values('AUC-ROC', ascending=True)

    if len(valid_df) == 0:
        print("No enzymes with valid AUC-ROC for per-enzyme figure.")
        return None

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Panel A: AUC-ROC distribution
    ax1 = axes[0]
    ax1.barh(range(len(valid_df)), valid_df['AUC-ROC'], color=COLOR_POSITIVE, alpha=0.7)
    ax1.axvline(x=0.5, color='k', linestyle='--', lw=1, alpha=0.5, label='Random (0.5)')
    ax1.axvline(x=valid_df['AUC-ROC'].median(), color=COLOR_NEUTRAL, linestyle='-', lw=2,
                label=f'Median ({valid_df["AUC-ROC"].median():.3f})')
    ax1.set_xlabel('AUC-ROC')
    ax1.set_ylabel('Enzyme')
    ax1.set_title('A. Per-Enzyme AUC-ROC', fontweight='bold', loc='left')
    ax1.set_yticks([])
    ax1.set_xlim([0, 1])
    ax1.legend(loc='lower right')

    # Panel B: Sample count vs AUC-ROC
    ax2 = axes[1]
    scatter = ax2.scatter(valid_df['N_samples'], valid_df['AUC-ROC'],
                          c=valid_df['N_positive'] / valid_df['N_samples'],
                          cmap='RdYlBu_r', alpha=0.7, s=60, edgecolors='k', linewidth=0.5)
    ax2.axhline(y=0.5, color='k', linestyle='--', lw=1, alpha=0.5)
    ax2.set_xlabel('Number of Samples')
    ax2.set_ylabel('AUC-ROC')
    ax2.set_title('B. AUC-ROC vs Sample Size', fontweight='bold', loc='left')
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label('Positive Fraction')

    fig.suptitle('Per-Enzyme Performance Analysis', fontsize=12, fontweight='bold')
    plt.tight_layout()

    return fig


def generate_report(metrics, per_enzyme_df, per_substrate_df, error_df, optimal_threshold):
    """Generate text analysis report."""
    report_lines = [
        "=" * 70,
        "EZSpecificity P450 Model Evaluation Report",
        "=" * 70,
        "",
        "1. GLOBAL PERFORMANCE METRICS",
        "-" * 40,
    ]

    for key, value in metrics.items():
        if isinstance(value, float):
            if np.isnan(value):
                report_lines.append(f"  {key}: N/A (undefined)")
            else:
                report_lines.append(f"  {key}: {value:.4f}")
        else:
            report_lines.append(f"  {key}: {value}")

    report_lines.extend([
        "",
        "2. THRESHOLD ANALYSIS",
        "-" * 40,
    ])
    if np.isnan(optimal_threshold[1]):
        report_lines.append("  Optimal threshold: N/A (single class)")
    else:
        report_lines.extend([
            f"  Optimal threshold (Youden's J): {optimal_threshold[0]:.4f}",
            f"  Youden's J statistic: {optimal_threshold[1]:.4f}",
        ])

    report_lines.extend([
        "",
        "3. PER-ENZYME ANALYSIS",
        "-" * 40,
        f"  Total enzymes analyzed: {len(per_enzyme_df)}",
    ])

    if len(per_enzyme_df) > 0 and 'AUC-ROC' in per_enzyme_df.columns:
        valid_auc_count = per_enzyme_df['AUC-ROC'].notna().sum()
        report_lines.append(f"  Enzymes with valid AUC-ROC: {valid_auc_count}")

        if valid_auc_count > 0:
            report_lines.extend([
                f"  Median AUC-ROC: {per_enzyme_df['AUC-ROC'].median():.4f}",
                f"  AUC-ROC range: [{per_enzyme_df['AUC-ROC'].min():.4f}, {per_enzyme_df['AUC-ROC'].max():.4f}]",
            ])

    report_lines.extend([
        "",
        "4. PER-SUBSTRATE ANALYSIS",
        "-" * 40,
        f"  Total substrates analyzed: {len(per_substrate_df)}",
    ])

    if len(per_substrate_df) > 0 and 'AUC-ROC' in per_substrate_df.columns:
        valid_auc_count = per_substrate_df['AUC-ROC'].notna().sum()
        report_lines.append(f"  Substrates with valid AUC-ROC: {valid_auc_count}")

        if valid_auc_count > 0:
            report_lines.extend([
                f"  Median AUC-ROC: {per_substrate_df['AUC-ROC'].median():.4f}",
                f"  AUC-ROC range: [{per_substrate_df['AUC-ROC'].min():.4f}, {per_substrate_df['AUC-ROC'].max():.4f}]",
            ])

    report_lines.extend([
        "",
        "5. ERROR ANALYSIS (threshold=0.5)",
        "-" * 40,
        f"  Total errors: {len(error_df)}",
    ])
    if len(error_df) > 0 and 'Error_Type' in error_df.columns:
        report_lines.extend([
            f"  False Positives: {(error_df['Error_Type'] == 'False Positive').sum()}",
            f"  False Negatives: {(error_df['Error_Type'] == 'False Negative').sum()}",
        ])

    report_lines.extend([
        "",
        "=" * 70,
        "Report generated by step10_analysis.py",
        "=" * 70,
    ])

    return "\n".join(report_lines)


def main():
    print("=" * 60)
    print("Step 10: EZSpecificity P450 Model Evaluation Analysis")
    print("=" * 60)

    # Create output directories
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    pred_df, enzymes_df, substrates_df = load_data()

    # Extract arrays
    y_true = pred_df['Label'].values
    y_prob = pred_df['prob'].values
    y_logit = pred_df['logit'].values

    # Compute global metrics
    print("\n" + "=" * 60)
    print("Computing global metrics...")
    metrics = compute_global_metrics(y_true, y_prob, y_logit)

    print("\nGlobal Performance Metrics:")
    for key, value in metrics.items():
        print(f"  {key}: {value:.4f}")

    # Find optimal threshold
    optimal_threshold = find_optimal_threshold(y_true, y_prob)
    print(f"\nOptimal threshold (Youden's J): {optimal_threshold[0]:.4f}")

    # Per-enzyme analysis
    print("\n" + "=" * 60)
    print("Computing per-enzyme metrics...")
    per_enzyme_df = compute_per_enzyme_metrics(pred_df, enzymes_df)
    print(f"  Analyzed {len(per_enzyme_df)} enzymes")

    # Per-substrate analysis
    print("\n" + "=" * 60)
    print("Computing per-substrate metrics...")
    per_substrate_df = compute_per_substrate_metrics(pred_df, substrates_df)
    print(f"  Analyzed {len(per_substrate_df)} substrates")

    # Error analysis
    print("\n" + "=" * 60)
    error_df = analyze_errors(pred_df)

    # Create visualizations
    print("\n" + "=" * 60)
    print("Creating visualizations...")

    # Main analysis figure
    fig_main = create_main_figure(y_true, y_prob, y_logit)
    fig_main.savefig(FIGURES_DIR / "main_analysis_figure.png", dpi=300)
    fig_main.savefig(FIGURES_DIR / "main_analysis_figure.pdf")
    print(f"  Saved: main_analysis_figure.png/pdf")
    plt.close(fig_main)

    # Per-enzyme figure
    fig_enzyme = create_per_enzyme_figure(per_enzyme_df)
    if fig_enzyme:
        fig_enzyme.savefig(FIGURES_DIR / "per_enzyme_analysis.png", dpi=300)
        print(f"  Saved: per_enzyme_analysis.png")
        plt.close(fig_enzyme)

    # Save data tables
    print("\n" + "=" * 60)
    print("Saving data tables...")

    # Metrics summary
    metrics_df = pd.DataFrame([metrics])
    metrics_df.to_csv(OUTPUT_DIR / "metrics_summary.csv", index=False)
    print(f"  Saved: metrics_summary.csv")

    # Per-enzyme metrics
    per_enzyme_df.to_csv(OUTPUT_DIR / "per_enzyme_metrics.csv", index=False)
    print(f"  Saved: per_enzyme_metrics.csv")

    # Per-substrate metrics
    per_substrate_df.to_csv(OUTPUT_DIR / "per_substrate_metrics.csv", index=False)
    print(f"  Saved: per_substrate_metrics.csv")

    # Error analysis
    error_df.to_csv(OUTPUT_DIR / "error_analysis.csv", index=False)
    print(f"  Saved: error_analysis.csv")

    # Generate and save report
    report = generate_report(metrics, per_enzyme_df, per_substrate_df, error_df, optimal_threshold)
    with open(OUTPUT_DIR / "analysis_report.txt", "w", encoding="utf-8") as f:
        f.write(report)
    print(f"  Saved: analysis_report.txt")

    # Print final summary
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nKey Results:")
    print(f"  AUC-ROC: {metrics['AUC-ROC']:.4f}")
    print(f"  AUC-PR:  {metrics['AUC-PR (AUPR)']:.4f}")
    print(f"  F1:      {metrics['F1 Score']:.4f}")
    print(f"  MCC:     {metrics['MCC']:.4f}")
    print(f"\nOutput directory: {OUTPUT_DIR}")

    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
