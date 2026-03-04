# Phase 2: Structure-Source Ablation Report

**Date**: 2026-03-04 10:11
**Branch**: pathb-ablation

## 2×2 Matrix Design

| | Inhibitor Neg | Random Neg |
|---|---|---|
| **Crystal** | EXP01 = **0.7115** | N/A (impossible) |
| **Vina** | ABL-04 = **0.6885** | ABL-03 = **0.5170** |

## All Experiment Results

| Experiment | AUC-ROC | 95% CI | AUC-PR | Pos | Neg | Score Sep |
|-----------|---------|--------|--------|-----|-----|----------|
| EXP01 | **0.7115** | [0.6657, 0.7569] | 0.7356 | 267 | 228 | 2.5126 |
| ABL-04 | **0.6885** | [0.6392, 0.7348] | 0.7180 | 265 | 222 | 2.3107 |

## Gap Decomposition

**Total gap** = EXP01 - ABL-03 = 0.1945

**Structure effect** = EXP01 - ABL-04 = 0.0230
  (Crystal → Vina, holding neg type = inhibitor)

**Neg-type effect** = ABL-04 - ABL-03 = 0.1715
  (Inhibitor → Random, holding structure = Vina)

**Consistency check**: 0.0230 + 0.1715 = 0.1945 (total gap = 0.1945)

### Relative Contributions

- Structure source: 11.8% of total gap
- Negative type: 88.2% of total gap

## Interpretation

**Negative type is the dominant factor.** Switching from inhibitor to random negatives causes most of the performance drop. This confirms that the task difficulty (inhibitor vs random) is the primary driver, not structure quality.

## Matched Cohort Analysis (EXP01 vs ABL-04)

- Matched pairs: 470 (261 pos + 209 neg)
- EXP01 AUC (matched): 0.7129
- ABL-04 AUC (matched): 0.6805
- AUC difference: 0.0324 [0.0138, 0.0521]

The matched cohort controls for enzyme-substrate pair composition, isolating the pure structure-source effect.

## Files Generated

- `phase2_metrics.csv` — Experiment metrics
- `phase2_roc.png` — ROC curve comparison
- `phase2_heatmap.png` — 2×2 ablation heatmap
- `phase2_score_dist.png` — Score distributions
- `ABL-04_*/predictions.csv` — ABL-04 predictions
