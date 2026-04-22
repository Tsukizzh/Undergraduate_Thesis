"""Phase 3: Generate 6 new index.pt files (3 natural + 3 unified).

For each base _fixed cache and each split:
1. Load index.pt
2. Filter rows where substrate_id is NOT in shared substrates index (drops sub_id=8)
3. Output "natural" filtered index.pt (keeping original sample_ids/graph_shards/graph_rows)

Then compute unified sample set:
- unified_sample_ids = intersection of the 3 natural sets (by sample_id)
- For each base, produce a "unified" index.pt by projecting unified_sample_ids back
  to its own rows (keeping base-specific graph_shards/graph_rows)

Output paths (written to temp staging, copied to _allfix dirs in Phase 4):
    /tmp/allfix_indices/{name}/{split}/index.pt
where name in {bare_natural, heme_natural, geom_natural, bare_unified, heme_unified, geom_unified}

Codex-reviewed round 6.
"""
import torch
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
STAGING = Path("/tmp/allfix_indices")
SPLITS = ['train', 'val', 'test']

BASE_CACHES = {
    'bare': BASE / "data/pt_cache_fixed/random",
    'heme': BASE / "data/pt_cache_heme_fixed/random",
    'geom': BASE / "data/pt_cache_geom_fixed",
}

# Load shared substrates index
shared_sub_idx = torch.load(
    BASE / "data/pt_cache_allfix_shared/substrates/substrates_index.pt",
    map_location="cpu", weights_only=False
)
shared_sub_set = set(shared_sub_idx['index'].keys())
print(f"Shared substrates: {len(shared_sub_set)} (has 8? {8 in shared_sub_set})")
assert 8 not in shared_sub_set


def filter_natural(base_idx: dict, allowed_subs: set) -> dict:
    """Filter rows where substrate_id is in allowed_subs. Keep all fields as-is (no renumber)."""
    sids = base_idx['sample_ids']
    ezids = base_idx['enzyme_ids']
    sbids = base_idx['substrate_ids']
    shards = base_idx['graph_shards']
    rows = base_idx['graph_rows']

    keep_mask = torch.zeros(len(sids), dtype=torch.bool)
    for i, sb in enumerate(sbids.tolist()):
        if int(sb) in allowed_subs:
            keep_mask[i] = True

    return {
        'sample_ids':    sids[keep_mask].clone(),
        'enzyme_ids':    ezids[keep_mask].clone(),
        'substrate_ids': sbids[keep_mask].clone(),
        'graph_shards':  shards[keep_mask].clone(),
        'graph_rows':    rows[keep_mask].clone(),
    }


def project_to_unified(base_idx: dict, unified_sample_ids: set) -> dict:
    """Keep only rows whose sample_id is in unified_sample_ids."""
    sids = base_idx['sample_ids']
    mask = torch.tensor([int(s) in unified_sample_ids for s in sids.tolist()], dtype=torch.bool)
    return {
        'sample_ids':    sids[mask].clone(),
        'enzyme_ids':    base_idx['enzyme_ids'][mask].clone(),
        'substrate_ids': base_idx['substrate_ids'][mask].clone(),
        'graph_shards':  base_idx['graph_shards'][mask].clone(),
        'graph_rows':    base_idx['graph_rows'][mask].clone(),
    }


STAGING.mkdir(parents=True, exist_ok=True)

print("\n=== Phase 3a: Natural filtered indices ===")
natural = {}  # {cache_name: {split: filtered_idx}}
for name, cache_path in BASE_CACHES.items():
    natural[name] = {}
    for split in SPLITS:
        idx_path = cache_path / split / "index.pt"
        base_idx = torch.load(idx_path, map_location="cpu", weights_only=False)
        filtered = filter_natural(base_idx, shared_sub_set)
        natural[name][split] = filtered
        n_before = len(base_idx['sample_ids'])
        n_after = len(filtered['sample_ids'])
        dropped = n_before - n_after
        print(f"  {name}/{split}: {n_before} -> {n_after}  (-{dropped} dropped)")

        # Sanity: no sub_id=8 remains
        assert 8 not in filtered['substrate_ids'].tolist()

        # Save to staging
        out_dir = STAGING / f"{name}_natural" / split
        out_dir.mkdir(parents=True, exist_ok=True)
        torch.save(filtered, out_dir / "index.pt")

print("\n=== Phase 3b: Compute unified sample_id sets ===")
unified_sample_ids = {}
for split in SPLITS:
    bare_sids = set(natural['bare'][split]['sample_ids'].tolist())
    heme_sids = set(natural['heme'][split]['sample_ids'].tolist())
    geom_sids = set(natural['geom'][split]['sample_ids'].tolist())
    union = bare_sids | heme_sids | geom_sids
    inter = bare_sids & heme_sids & geom_sids
    unified_sample_ids[split] = inter
    print(f"  {split}: union={len(union)}, intersection={len(inter)}  (bare={len(bare_sids)}, heme={len(heme_sids)}, geom={len(geom_sids)})")

print("\n=== Phase 3c: Project each cache to unified set ===")
for name in BASE_CACHES:
    for split in SPLITS:
        projected = project_to_unified(natural[name][split], unified_sample_ids[split])
        n_before = len(natural[name][split]['sample_ids'])
        n_after = len(projected['sample_ids'])
        print(f"  {name}_unified/{split}: {n_before} -> {n_after}")

        # Sanity: no sub_id=8
        assert 8 not in projected['substrate_ids'].tolist()

        # Sanity: sample_ids matches the unified set
        actual_sids = set(projected['sample_ids'].tolist())
        assert actual_sids == unified_sample_ids[split], (
            f"unified projection sample_id mismatch for {name}/{split}"
        )

        out_dir = STAGING / f"{name}_unified" / split
        out_dir.mkdir(parents=True, exist_ok=True)
        torch.save(projected, out_dir / "index.pt")

print("\n=== Phase 3d: Cross-cache consistency check on unified indices ===")
# Verify all 3 unified indices for the same split have the same sorted sample_ids
for split in SPLITS:
    b = torch.load(STAGING / "bare_unified" / split / "index.pt", weights_only=False)
    h = torch.load(STAGING / "heme_unified" / split / "index.pt", weights_only=False)
    g = torch.load(STAGING / "geom_unified" / split / "index.pt", weights_only=False)

    bs = sorted(b['sample_ids'].tolist())
    hs = sorted(h['sample_ids'].tolist())
    gs = sorted(g['sample_ids'].tolist())

    assert bs == hs == gs, f"unified sample_id mismatch in {split}"
    print(f"  {split}: 3 unified indices have identical sample_ids ({len(bs)} each)")

print("\n=== Phase 3 DONE ===")
print(f"Staging directory: {STAGING}")
for sub in sorted(STAGING.iterdir()):
    print(f"  {sub.name}/")
    for split_dir in sorted(sub.iterdir()):
        idx_file = split_dir / "index.pt"
        if idx_file.exists():
            d = torch.load(idx_file, weights_only=False)
            print(f"    {split_dir.name}: {len(d['sample_ids'])} rows")
