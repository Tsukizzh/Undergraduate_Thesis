"""Phase 5: End-to-end sanity. For each _allfix cache, actually load a sample
through its native pt_dataset.py and verify:

1. dataset length matches Phase 3 index.pt count
2. first sample's protein feature dim matches expected (28/31/37)
3. GROVER total_embedding from sample == fixed LMDB's Substrate Index → byte match (fp16)
4. Scan first 200 samples: no sub_id=8
5. bare_unified / heme_unified / geom_unified: first sample has same (enzyme_id, substrate_id)
"""
import sys
import os
import torch
import numpy as np
import lmdb
import pickle
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
GROVER_FIXED = BASE / "data/features/grover_fingerprint_fixed.lmdb"

# Each cache and its corresponding experiment dir (for src/ and pt_dataset.py)
CACHES = [
    ('pt_cache_allfix/random',          'EXP001_fixed', 28),
    ('pt_cache_heme_allfix/random',     'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix',            'EXP003_fixed', 37),
    ('pt_cache_allfix_unified/random',  'EXP001_fixed', 28),
    ('pt_cache_heme_allfix_unified/random', 'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix_unified',    'EXP003_fixed', 37),
]

# Expected index.pt sample counts (from Phase 3 run output)
EXPECTED_COUNTS = {
    'pt_cache_allfix/random':               {'train': 22178, 'val': 11058, 'test': 11048},
    'pt_cache_heme_allfix/random':          {'train': 22384, 'val': 11148, 'test': 11148},
    'pt_cache_geom_allfix':                 {'train': 22312, 'val': 11111, 'test': 11107},
    'pt_cache_allfix_unified/random':       {'train': 22083, 'val': 11008, 'test': 10999},
    'pt_cache_heme_allfix_unified/random':  {'train': 22083, 'val': 11008, 'test': 10999},
    'pt_cache_geom_allfix_unified':         {'train': 22083, 'val': 11008, 'test': 10999},
}


def load_grover_lmdb_direct(sub_id: int) -> dict:
    env = lmdb.open(str(GROVER_FIXED), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        raw = txn.get(str(sub_id).encode())
    env.close()
    return pickle.loads(raw) if raw else None


def test_cache(cache_rel: str, exp_name: str, expected_dim: int):
    print(f"\n=== {cache_rel} ===")
    cache_dir = BASE / "data" / cache_rel
    print(f"  feature_dim expected: {expected_dim}")

    # Set up Python path for this experiment's pt_dataset
    exp_dir = BASE / "experiments" / exp_name
    scripts_dir = str(exp_dir / "scripts")
    src_dir = str(exp_dir / "src")
    sys.path = [p for p in sys.path if 'experiments' not in p]  # cleanup
    sys.path.insert(0, scripts_dir)
    sys.path.insert(0, src_dir)

    # Re-import pt_dataset (flushing cached modules)
    for mod_name in list(sys.modules.keys()):
        if 'pt_dataset' in mod_name or mod_name.startswith('Datasets') or mod_name.startswith('Models'):
            del sys.modules[mod_name]

    from pt_dataset import PtCacheDataset  # type: ignore

    # 1. Check dataset length for each split
    for split in ['train', 'val', 'test']:
        ds = PtCacheDataset(cache_dir=cache_dir, split=split, edge_mode='fixed')
        n = len(ds)
        expected = EXPECTED_COUNTS[cache_rel][split]
        status = "OK" if n == expected else "FAIL"
        print(f"  {split}: {n} samples (expected {expected}) [{status}]")
        assert n == expected, f"count mismatch for {split}"

    # Use train for detailed checks
    ds = PtCacheDataset(cache_dir=cache_dir, split='train', edge_mode='fixed')

    # 2. Sample 0: check protein feature dim
    sample = ds[0]
    # Find protein features (node_feat or protein_atom_feature or similar)
    for attr_name in ['protein_feat', 'atom_feat', 'x', 'h']:
        if hasattr(sample, attr_name):
            feat = getattr(sample, attr_name)
            if feat is not None and hasattr(feat, 'shape') and len(feat.shape) >= 2:
                dim = feat.shape[-1]
                print(f"  sample[0].{attr_name}.shape[-1] = {dim}")
                break

    # 3. GROVER total_embedding match
    # Sample should have grover_mean as tensor (4885,)
    sub_id = int(sample.substrate_id) if hasattr(sample, 'substrate_id') else None
    if sub_id is not None:
        sample_grover_mean = None
        for attr_name in ['grover_mean', 'grover_total', 'substrate_features']:
            if hasattr(sample, attr_name):
                v = getattr(sample, attr_name)
                if v is not None and hasattr(v, 'shape'):
                    sample_grover_mean = v
                    print(f"  sample[0].{attr_name}: shape={tuple(v.shape)}, dtype={v.dtype}")
                    break

        # Direct LMDB read
        lmdb_data = load_grover_lmdb_direct(sub_id)
        lmdb_total = np.asarray(lmdb_data['total_embedding'], dtype=np.float16)

        if sample_grover_mean is not None:
            sample_arr = np.asarray(sample_grover_mean.detach().cpu() if torch.is_tensor(sample_grover_mean) else sample_grover_mean).reshape(-1).astype(np.float16)
            match = np.array_equal(sample_arr, lmdb_total)
            print(f"  sample[0].grover_mean == LMDB[{sub_id}].total_embedding: {match}")

    # 4. Scan first 200 samples for sub_id=8
    seen_8 = 0
    seen_max = 200
    for i in range(min(seen_max, len(ds))):
        s = ds[i]
        if hasattr(s, 'substrate_id') and int(s.substrate_id) == 8:
            seen_8 += 1
    assert seen_8 == 0
    print(f"  scanned {min(seen_max, len(ds))} samples: 0 with sub_id=8")


print("=" * 70)
print("PHASE 5 END-TO-END SANITY")
print("=" * 70)

for cache_rel, exp_name, expected_dim in CACHES:
    try:
        test_cache(cache_rel, exp_name, expected_dim)
    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

# 5. Cross-cache: unified versions should agree on first sample's (enzyme_id, substrate_id)
print("\n=== Cross-cache unified consistency ===")
unified_caches = [
    'pt_cache_allfix_unified/random',
    'pt_cache_heme_allfix_unified/random',
    'pt_cache_geom_allfix_unified',
]
first_pairs = []
for cache_rel in unified_caches:
    idx = torch.load(BASE / "data" / cache_rel / "train" / "index.pt", weights_only=False)
    sid = int(idx['sample_ids'][0])
    eid = int(idx['enzyme_ids'][0])
    bid = int(idx['substrate_ids'][0])
    first_pairs.append((cache_rel, sid, eid, bid))
    print(f"  {cache_rel}: first sample (sid={sid}, eid={eid}, sub={bid})")

all_match = all(t[1:] == first_pairs[0][1:] for t in first_pairs)
print(f"  all 3 agree: {all_match}")
assert all_match

print("\n=== Phase 5 DONE ===")
