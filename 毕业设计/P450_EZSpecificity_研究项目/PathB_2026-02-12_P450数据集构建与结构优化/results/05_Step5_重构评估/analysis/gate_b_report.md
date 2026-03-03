# Gate B Decision Report

**Date**: 2026-03-03
**Decision**: INFORMATIVE FAIL

## Summary

| Dataset | Neg Type | Structure | n_total | Prevalence | AUC-ROC (95% CI) | AUC-PR (95% CI) |
|---------|----------|-----------|---------|------------|-------------------|-----------------|
| Step 5 (Random neg, Vina) | Random | Vina-docked | 2766 | 0.096 | **0.5170** [0.4804, 0.5521] | 0.1116 [0.0926, 0.1389] |
| EXP01 (Inhibitor neg, Crystal) | Inhibitor | PDB crystal | 495 | 0.539 | **0.7115** [0.6657, 0.7569] | 0.7356 [0.6777, 0.7928] |
| Paper (Unknown E+S) | Random | AutoDock | ~65K | ~0.096 | **0.7198** | — |

> **Note on AUC-PR comparability**: Step 5 prevalence (0.096) differs substantially from EXP01 (0.539)
> AUC-PR is prevalence-sensitive; cross-dataset comparison should focus on AUC-ROC.

## Key Finding

**AUC-ROC = 0.5170 [0.4804, 0.5521]** / 0.7115 (EXP01), near random baseline (0.5).

The model shows negligible discriminative power between real substrates
and randomly paired molecules when both are Vina-docked into unseen P450 enzymes.

## Root Cause Analysis

### Score Distribution

| Dataset | Positive mean | Negative mean | Separation |
|---|---|---|---|
| Step 5 (Random neg, Vina) | -2.9369 | -3.2017 | 0.2648 |
| EXP01 (Inhibitor neg, Crystal) | -3.0017 | -5.5144 | 2.5126 |

Positive scores are consistent across settings (~-3.0).
EXP01 negatives (inhibitors, crystal) score much lower (-5.51), providing clear separation.
Step 5 negatives (random, Vina-docked) score close to positives (-3.20), collapsing the separation to 0.26.

### Hypothesized Root Causes (ranked by plausibility)

1. **Dockability != Catalysis (likely primary factor)**: Vina optimizes binding poses
   for all pairs, producing physically plausible complexes regardless of true catalytic
   relationship. The model may interpret these plausible poses as genuine binding.

2. **OOD Enzymes (0% overlap)**: None of our 152 P450 enzymes appear in the training
   set. The model lacks learned P450-specific binding preferences, compounding (1).

### Design Limitation: Confounded Variables

Step 5 changed **two variables simultaneously** relative to EXP01:

| | Crystal structure | Vina-docked |
|---|---|---|
| **Inhibitor negatives** | EXP01 (AUC=0.7115) | — |
| **Random negatives** | — (not tested) | Step 5 (AUC=0.5170) |

Without the two missing cells (crystal+random, Vina+inhibitor), we cannot
definitively attribute the AUC drop to negative type vs. structure source.
The observed result is **consistent with** the dockability hypothesis but does
not constitute causal proof.

### Per-Enzyme Analysis (min_support>=3)

- Enzymes analyzed: 21
- Mean per-enzyme AUC: 0.5543
- Median per-enzyme AUC: 0.5469
- AUC > 0.7: 8 (38.1%)
- AUC > 0.5: 11 (52.4%)
- AUC < 0.5: 9 (42.9%)

## Gate B Decision

**INFORMATIVE FAIL**: The random negative + Vina docking strategy yields
AUC-ROC = 0.5170, indistinguishable from random.
This is scientifically informative:

1. It suggests the model relies heavily on structural feature quality
2. It indicates that Vina docking poses alone do not encode catalytic specificity
3. It is consistent with the 0% enzyme overlap being a fundamental limitation

## Implications for Path C (Model Training)

- **Fine-tuning appears necessary**: The pretrained model does not generalize to P450s
- **Training data strategy**: P450-specific training data with known substrates would help
- **Negative sampling**: Random negatives with Vina docking may create ambiguous signal;
  alternative strategies (non-docked negatives, harder negative mining) should be explored
- **Structure quality**: Crystal structures likely provide more reliable evaluation signal
  than Vina-docked structures (but this needs controlled testing)

## Recommended Next Steps

1. **Control experiment (optional)**: Crystal positives + Vina random negatives
   (isolate structure-source effect from negative-definition effect)
2. **Proceed to Path C**: Use accumulated data for P450-specific fine-tuning
3. **Alternative evaluation**: Consider per-enzyme ranking metrics alongside global AUC

## Methodology Notes

- Bootstrap CI: 2000 iterations, seed=42, percentile method
- Threshold-based metrics (accuracy, precision, recall, F1) use in-sample Youden J
  optimization and should NOT be interpreted as out-of-sample performance
- Per-enzyme AUC minimum support: n_pos >= 3 AND n_neg >= 3
- Model checkpoint: EZSpecificity pretrained (Nature 2025), no P450 fine-tuning
