# Phase 1: Ratio Ablation Report

**Date**: 2026-03-04 03:20
**Branch**: pathb-ablation

## Experimental Design

Isolate the effect of positive:negative ratio on AUC-ROC,
holding structure source (Vina) and negative type (random) constant.

| Experiment | Ratio | Pos | Neg | Structure | Neg Type |
|-----------|-------|-----|-----|-----------|----------|
| ABL-01 | 1:1.0 | 265 | 265 | Vina | Random |
| ABL-02 | 1:3.0 | 265 | 795 | Vina | Random |
| ABL-03 | 1:9.4 | 265 | 2501 | Vina | Random |

## Results

| Experiment | AUC-ROC | 95% CI | AUC-PR | Score Sep |
|-----------|---------|--------|--------|----------|
| ABL-01 | **0.5154** | [0.4655, 0.5644] | 0.5148 | 0.3320 |
| ABL-02 | **0.5132** | [0.4731, 0.5533] | 0.2752 | 0.2735 |
| ABL-03 | **0.5170** | [0.4804, 0.5521] | 0.1116 | 0.2648 |
| EXP01 | **0.7115** | [0.6657, 0.7569] | 0.7356 | 2.5126 |

### ABL-01 Multi-Seed Results (5 seeds)

- Mean AUC-ROC: **0.5176** ± 0.0160
- Range: [0.4944, 0.5447]
- Seeds: [42, 123, 456, 789, 2026]
- AUCs: [0.5151, 0.5447, 0.5182, 0.4944, 0.5154]

### ABL-02 Multi-Seed Results (5 seeds)

- Mean AUC-ROC: **0.5146** ± 0.0051
- Range: [0.5093, 0.5234]
- Seeds: [42, 123, 456, 789, 2026]
- AUCs: [0.5234, 0.5167, 0.5093, 0.5102, 0.5132]

## Interpretation

**Ratio Effect**: ABL-01 (1:1) - ABL-03 (1:9) = -0.0016

Ratio has **minimal effect** on AUC-ROC. The performance drop from EXP01 to Step 5 is NOT primarily due to class imbalance.

**Next step**: Phase 2 — isolate structure source vs negative type effect.

## Files Generated

- `phase1_metrics.csv` — All experiment metrics
- `phase1_roc.png` — ROC curve comparison
- `phase1_auc_vs_ratio.png` — AUC vs ratio plot
- `phase1_score_dist.png` — Score distributions
- `ABL-*/predictions*.csv` — Per-experiment predictions
