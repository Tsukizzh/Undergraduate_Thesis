# EXP004_paper_baseline_unified

**Purpose**: Run the paper's pre-trained checkpoint on our P450 test set as an external baseline. Test samples whose enzyme overlaps with the 389 ESIBank P450 training enzymes are removed to prevent leakage.

## Checkpoints

| File | Source | Notes |
|---|---|---|
| `results/checkpoints/paper_best-checkpoint.ckpt` | `saved_model/model/run_0/models/best-checkpoint.ckpt` (Nature 2025) | PL 1.9.0, random_split, edge-mode=legacy_bug implicit, md5=f4d87ea08fc64b62700aadef8bf151cf |
| `results/checkpoints/ours_EXP001_ep43_auc0.9250.ckpt` | symlink → EXP001_allfix_unified ep43 best | Our allfix_unified bare baseline, edge-mode=fixed |

## Cache

```
--cache-dir /root/autodl-tmp/.../data/pt_cache_allfix_unified_paperfilter/random
```

Non-destructive overlay of `pt_cache_allfix_unified/random/` — everything except `test/index.pt` is a symlink. The new `test/index.pt` drops 3036 samples (27.6%) whose `enzyme_ids` is in the 356-enzyme blacklist, leaving 7963 samples (1125 unique enzymes) for evaluation.

## Scripts

- `scripts/run_train.sh` — single run, main result: paper ckpt + legacy_bug + filtered test
- `scripts/run_test_grid.sh` — 4-run diagnostic grid:
  1. `paper_legacy_filtered` (main external baseline)
  2. `paper_fixed_filtered` (edge sensitivity for paper weights)
  3. `ours_legacy_filtered` (symmetric control — our weights should prefer fixed)
  4. `ours_fixed_filtered` (our native mode on filtered subset, compare to full 0.9320)

All runs write to `results/test_eval_<tag>.json`.

## Preflight results (all PASS)

- **Blacklist**: 356/389 ESIBank P450 are in our dataset (91.5%). 33 not in our data → no residual leakage from those.
- **Filter**: 3036/10999 test samples dropped (27.6%). Kept 7963 samples, 1125 unique enzymes. Below 30% threshold → acceptable external baseline.
- **Cache verify** (A/B/C/D): enzyme_id→seq_len mapping intact, index arrays self-consistent, 17 spread samples load cleanly, original cache untouched.
- **Ckpt preflight** (local): 76/76 state_dict keys match between paper ckpt and `SS(config)`, 0 missing / 0 unexpected / 0 shape mismatches, `strict=True` PASSED. torch 2.3.0 + PL 1.9.0 (matches paper training environment).

## Known limitations

- Sequence-hash fallback is a no-op because we only have the 389 ESIBank UniProt IDs, not their sequences. Any aliasing (same sequence, different accession) is not caught. Worst case: 33 hidden matches hidden in our 165 synthetic-ID (`ENZ_G*`) enzymes.
- Filter is by enzyme only, not substrate. The paper's random_split does have substrate-level overlap with our test; this is not controlled here.
- Paper ckpt trained with edge-mode=legacy_bug (because that was the buggy default). We match by running with `--edge-mode legacy_bug`. `paper_fixed_filtered` in the grid shows the counterfactual.

## Expected runtime on 4×RTX4090

- Single eval: ~3-5 minutes (7963 samples, batch_size=88)
- Full 4-run grid: ~15-20 minutes total
- User must spin up GPU manually before running.
