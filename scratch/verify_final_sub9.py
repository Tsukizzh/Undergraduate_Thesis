"""Final verification: find samples with substrate_id=9 (shift start) in each cache,
load them through PtCacheDataset, verify grover_mean matches fixed LMDB byte-for-byte.

This is the critical test: sub_id=9 is where the old LMDB bug started shifting data.
If this works, the fix is proven at the exact bug boundary.

Also verify a sample with sub_id in the 'critical region' (e.g., 10, 50, 2100, 2124).
"""
import sys
import torch
import numpy as np
import lmdb
import pickle
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
GROVER_FIXED = BASE / "data/features/grover_fingerprint_fixed.lmdb"
GROVER_OLD = BASE / "data/features/grover_fingerprint.lmdb"


def setup_exp_path(exp_name: str) -> None:
    exp_dir = BASE / "experiments" / exp_name
    sys.path = [p for p in sys.path if 'experiments' not in p]
    sys.path.insert(0, str(exp_dir / "scripts"))
    sys.path.insert(0, str(exp_dir / "src"))
    for mod_name in list(sys.modules.keys()):
        if 'pt_dataset' in mod_name or mod_name.startswith('Datasets') or mod_name.startswith('Models') or mod_name == 'build_pt_cache':
            del sys.modules[mod_name]


def load_lmdb(path, sub_id):
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        raw = txn.get(str(sub_id).encode())
    env.close()
    return pickle.loads(raw) if raw else None


CACHES = [
    ('pt_cache_allfix/random',              'EXP001_fixed'),
    ('pt_cache_heme_allfix/random',         'EXP002a_fixed'),
    ('pt_cache_geom_allfix',                'EXP003_fixed'),
    ('pt_cache_allfix_unified/random',      'EXP001_fixed'),
    ('pt_cache_heme_allfix_unified/random', 'EXP002a_fixed'),
    ('pt_cache_geom_allfix_unified',        'EXP003_fixed'),
]

# Critical substrate ids to check (shift boundary + spot checks)
CRITICAL_SUB_IDS = [7, 9, 10, 100, 1000, 2123, 2124]


for cache_rel, exp_name in CACHES:
    print(f"\n=== {cache_rel} ===")
    setup_exp_path(exp_name)
    from pt_dataset import PtCacheDataset

    cache_dir = BASE / "data" / cache_rel
    ds = PtCacheDataset(cache_dir=cache_dir, split='train', edge_mode='fixed')
    idx = torch.load(cache_dir / "train" / "index.pt", weights_only=False)
    sub_ids_tensor = idx['substrate_ids']

    # For each critical sub_id, find FIRST sample where substrate_ids == target
    for target_sub in CRITICAL_SUB_IDS:
        matches = (sub_ids_tensor == target_sub).nonzero(as_tuple=True)[0]
        if len(matches) == 0:
            print(f"  sub_id={target_sub}: no sample in this cache's train split")
            continue

        row = int(matches[0])
        sample = ds[row]
        gm = sample['grover_mean']
        sample_arr = gm.detach().cpu().numpy().reshape(-1).astype(np.float16)

        # Compare to fixed LMDB
        fixed_data = load_lmdb(GROVER_FIXED, target_sub)
        fixed_total = np.asarray(fixed_data['total_embedding'], dtype=np.float16)

        match_fixed = np.array_equal(sample_arr, fixed_total)

        # Also load OLD LMDB to show it WOULD NOT match
        old_data = load_lmdb(GROVER_OLD, target_sub)
        old_total = np.asarray(old_data['total_embedding'], dtype=np.float16) if old_data else None

        if old_total is not None:
            match_old = np.array_equal(sample_arr, old_total)
        else:
            match_old = None

        status = "FIX OK" if match_fixed and not match_old else ("BUG STILL" if match_old else "UNCLEAR")
        print(f"  sub_id={target_sub:4d} (row={row}): "
              f"sample==fixed_lmdb[{target_sub}]={match_fixed}, "
              f"sample==old_lmdb[{target_sub}]={match_old}  [{status}]")
        assert match_fixed, f"sub_id={target_sub}: sample doesn't match FIXED LMDB"

        # For k >= 9, we expect sample != old_lmdb (since old had shifted data)
        # For k <= 7, sample should equal both (old and fixed agree on these)
        if target_sub <= 7:
            assert match_old, f"sub_id={target_sub}: old LMDB was already correct for this range, should match"
        elif target_sub == 2124:
            # old LMDB has no key 2124 at all
            assert old_data is None, f"sub_id=2124: old LMDB should have no entry"
        else:
            # The critical check: sample MUST NOT match old LMDB's wrong data
            assert not match_old, f"sub_id={target_sub}: sample matches OLD LMDB — bug not fixed!"

print("\n=== FINAL VERIFICATION PASS ===")
print("All 6 _allfix caches correctly load the FIXED GROVER data at every critical substrate index.")
print("None of them accidentally serve the OLD (shifted) data.")
