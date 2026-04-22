"""Verify sample_id semantics are consistent across bare/heme/geom _fixed caches.

If sample_id=N in all 3 caches refers to the SAME (enzyme_id, substrate_id) pair,
then unified sample set can be computed by intersecting sample_id sets.

Otherwise, need to use (enzyme_id, substrate_id) tuple as cross-cache key.
"""
import torch

BASE = "/root/autodl-tmp/EZSpecificity/PathC/P450"

def load(path):
    return torch.load(path, map_location="cpu", weights_only=False)

def report(name, idx):
    n = len(idx['sample_ids'])
    print(f"  {name}: {n} rows")
    return n

print("=== Loading _fixed index.pt for all three caches ===")
caches = {}
for split in ['train', 'val', 'test']:
    print(f"\n--- {split} ---")
    bare = load(f"{BASE}/data/pt_cache_fixed/random/{split}/index.pt")
    heme = load(f"{BASE}/data/pt_cache_heme_fixed/random/{split}/index.pt")
    geom = load(f"{BASE}/data/pt_cache_geom_fixed/{split}/index.pt")

    report('bare', bare)
    report('heme', heme)
    report('geom', geom)

    caches[split] = {'bare': bare, 'heme': heme, 'geom': geom}

print("\n\n=== Sample ID consistency check ===")
for split, c in caches.items():
    print(f"\n--- {split} ---")
    # Build dicts: sample_id -> (enzyme_id, substrate_id)
    dicts = {}
    for name in ['bare', 'heme', 'geom']:
        idx = c[name]
        sids = idx['sample_ids'].tolist()
        ezids = idx['enzyme_ids'].tolist()
        sbids = idx['substrate_ids'].tolist()
        d = {int(s): (int(e), int(b)) for s, e, b in zip(sids, ezids, sbids)}
        dicts[name] = d
        # Sanity: no duplicate sample_ids within a cache
        assert len(d) == len(sids), f"{name} has duplicate sample_ids in {split}"
        print(f"  {name}: {len(d)} unique sample_ids, range {min(d.keys())}..{max(d.keys())}")

    # For each sample_id that appears in multiple caches, check agreement
    all_sids = set(dicts['bare'].keys()) | set(dicts['heme'].keys()) | set(dicts['geom'].keys())
    common_all = set(dicts['bare'].keys()) & set(dicts['heme'].keys()) & set(dicts['geom'].keys())
    print(f"  any: {len(all_sids)}, common to all 3: {len(common_all)}")

    # Pairwise intersections
    bh = set(dicts['bare'].keys()) & set(dicts['heme'].keys())
    bg = set(dicts['bare'].keys()) & set(dicts['geom'].keys())
    hg = set(dicts['heme'].keys()) & set(dicts['geom'].keys())
    print(f"  bare∩heme: {len(bh)}, bare∩geom: {len(bg)}, heme∩geom: {len(hg)}")

    # Check agreement on (enzyme_id, substrate_id) for common sample_ids
    disagreements = {'bh': 0, 'bg': 0, 'hg': 0, 'all': 0}
    for sid in common_all:
        b = dicts['bare'][sid]
        h = dicts['heme'][sid]
        g = dicts['geom'][sid]
        if b != h: disagreements['bh'] += 1
        if b != g: disagreements['bg'] += 1
        if h != g: disagreements['hg'] += 1
        if not (b == h == g): disagreements['all'] += 1

    print(f"  disagreements (among {len(common_all)} common sample_ids):")
    print(f"    bare vs heme: {disagreements['bh']}")
    print(f"    bare vs geom: {disagreements['bg']}")
    print(f"    heme vs geom: {disagreements['hg']}")
    print(f"    any disagreement: {disagreements['all']}")

    if disagreements['all'] == 0:
        print(f"  [OK] sample_id consistent across all 3 caches in {split}")
    else:
        print(f"  [WARNING] sample_id semantics differ across caches in {split}")
        # Show first few disagreements
        shown = 0
        for sid in sorted(common_all):
            b = dicts['bare'][sid]
            h = dicts['heme'][sid]
            g = dicts['geom'][sid]
            if not (b == h == g):
                print(f"    sid={sid}: bare={b}, heme={h}, geom={g}")
                shown += 1
                if shown >= 5:
                    break
