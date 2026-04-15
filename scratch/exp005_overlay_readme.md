# pt_cache_dualgraph_allfix_unified (overlay cache for EXP005)

**Status**: HALF-BUILT (2026-04-16)
**Base cache**: `pt_cache_allfix_unified/random` (unchanged)
**Purpose**: dualgraph (2+) experiment — atom-level + residue-level GVP

## What's here

```
random/
├── enzymes     -> symlink base (unchanged)
├── substrates  -> symlink base (unchanged)
├── manifest.pt -> symlink base (unchanged)
├── train/
│   ├── samples        -> symlink base     (sample .pt files, untouched)
│   ├── index.pt       -> symlink base     (base row-aligned index)
│   └── dock_sidecar.pt  ← NEW (Step 2)
├── val/    (same layout)
└── test/   (same layout)
```

## What's NOT here yet

- `gvp_cache/` — residue-level GVP features keyed by dock_index. Will be populated by Step 3 (batch GVP feature extraction on all pocket PDBs).

## dock_sidecar.pt format

Each split has its own sidecar:

```python
{
    "sample_ids":   int32 [N],   # copy from base index.pt for integrity
    "dock_indices": int32 [N],   # dock_index per sample row, row-aligned with base
    "source_csv":   str,         # which split CSV was used
    "source_index": str,         # which base index.pt was used
}
```

Usage: at runtime, `dock_sidecar["dock_indices"][k]` gives the dock_index for
the k-th sample in the base split index. Use this to look up per-sample GVP
features (once Step 3 populates gvp_cache/).

## Build provenance (Step 2, 2026-04-16)

- **Script**: `scratch/exp005_build_dock_sidecar.py`
- **Method**: for each row k in base index.pt, read split CSV row sample_ids[k],
  extract Dock Index / Enzyme Index / Substrate Index, verify full-match
  against base.enzyme_ids and base.substrate_ids, then save.
- **Integrity check passed**:
  - train: 22,083/22,083 (enzyme, substrate) match; 0 missing pocket PDB
  - val:   11,008/11,008 match; 0 missing PDB
  - test:  10,999/10,999 match; 0 missing PDB
- **Pre-condition**: Step 1 verified sample_id IS the original CSV row index
  (60/60 spot check across 3 splits).

## Do NOT modify

This overlay is non-destructive. If you need to revert, just delete this
`pt_cache_dualgraph_allfix_unified/` directory entirely — the base cache is
untouched. The symlinks are safe to delete; they do not affect base files.
