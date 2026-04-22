"""Phase 5 v3: use known sample keys (avoid dir() on PyG Data)."""
import sys
import os
import torch
import numpy as np
import lmdb
import pickle
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
GROVER_FIXED = BASE / "data/features/grover_fingerprint_fixed.lmdb"

SAMPLE_KEYS = ['ligand_index', 'protein_x', 'complex_edge_index', 'label',
               'sample_weight', 'morgan', 'ligand_x', 'ligand_mask',
               'protein_mask', 'tag', 'reaction_padding_mask', 'embedding',
               'complex_edge_attr', 'str_tag', 'grover_mean']


def setup_exp_path(exp_name: str) -> None:
    exp_dir = BASE / "experiments" / exp_name
    sys.path = [p for p in sys.path if 'experiments' not in p]
    sys.path.insert(0, str(exp_dir / "scripts"))
    sys.path.insert(0, str(exp_dir / "src"))
    for mod_name in list(sys.modules.keys()):
        if 'pt_dataset' in mod_name or mod_name.startswith('Datasets') or mod_name.startswith('Models') or mod_name == 'build_pt_cache':
            del sys.modules[mod_name]


def load_grover_lmdb(sub_id: int) -> dict:
    env = lmdb.open(str(GROVER_FIXED), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        raw = txn.get(str(sub_id).encode())
    env.close()
    return pickle.loads(raw) if raw else None


def deep_test(cache_rel: str, exp_name: str, expected_feat_dim: int):
    print(f"\n=== {cache_rel} (expects protein_x.shape[-1]={expected_feat_dim}) ===")
    setup_exp_path(exp_name)
    from pt_dataset import PtCacheDataset

    cache_dir = BASE / "data" / cache_rel

    # Load index.pt to get substrate_ids
    idx_path = cache_dir / "train" / "index.pt"
    idx = torch.load(idx_path, weights_only=False)

    ds = PtCacheDataset(cache_dir=cache_dir, split='train', edge_mode='fixed')
    print(f"  len(ds) = {len(ds)}")

    # Test sample 0 and a few more
    for test_i in [0, 100, 500]:
        if test_i >= len(ds):
            continue
        sample = ds[test_i]
        expected_sub_id = int(idx['substrate_ids'][test_i])
        expected_enz_id = int(idx['enzyme_ids'][test_i])
        expected_sid = int(idx['sample_ids'][test_i])

        # Access by getitem since dir() fails on this Data object
        protein_x = sample['protein_x'] if 'protein_x' in sample else None
        grover_mean = sample['grover_mean'] if 'grover_mean' in sample else None
        ligand_x = sample['ligand_x'] if 'ligand_x' in sample else None

        print(f"  sample[{test_i}] index says: sid={expected_sid}, enz={expected_enz_id}, sub={expected_sub_id}")

        if protein_x is not None:
            px = protein_x
            print(f"    protein_x: shape={tuple(px.shape)} dtype={px.dtype}")
            feat_dim = px.shape[-1]
            match = feat_dim == expected_feat_dim
            print(f"    feat_dim: {feat_dim} (expected {expected_feat_dim}) {'OK' if match else 'MISMATCH'}")
            assert match

        if grover_mean is not None and expected_sub_id is not None:
            gm = grover_mean
            print(f"    grover_mean: shape={tuple(gm.shape)} dtype={gm.dtype}")
            sample_arr = gm.detach().cpu().numpy().reshape(-1).astype(np.float16)

            lmdb_data = load_grover_lmdb(expected_sub_id)
            lmdb_total = np.asarray(lmdb_data['total_embedding'], dtype=np.float16)
            min_len = min(len(sample_arr), len(lmdb_total))
            match = np.array_equal(sample_arr[:min_len], lmdb_total[:min_len])
            print(f"    grover_mean vs lmdb[{expected_sub_id}].total_embedding: match={match} "
                  f"(sample_len={len(sample_arr)}, lmdb_len={len(lmdb_total)})")
            if not match:
                # Show diff
                diff = (sample_arr[:10].astype(np.float32) - lmdb_total[:10].astype(np.float32))
                print(f"    first-10 diff: {diff}")
                print(f"    sample[:5]: {sample_arr[:5]}")
                print(f"    lmdb[:5]:   {lmdb_total[:5]}")
            assert match, f"grover_mean mismatch for sub_id={expected_sub_id}"


CACHES = [
    ('pt_cache_allfix/random',              'EXP001_fixed', 28),
    ('pt_cache_heme_allfix/random',         'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix',                'EXP003_fixed', 37),
    ('pt_cache_allfix_unified/random',      'EXP001_fixed', 28),
    ('pt_cache_heme_allfix_unified/random', 'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix_unified',        'EXP003_fixed', 37),
]

for cache_rel, exp, dim in CACHES:
    try:
        deep_test(cache_rel, exp, dim)
    except Exception as e:
        print(f"  ERROR in {cache_rel}: {type(e).__name__}: {e}")
        import traceback; traceback.print_exc()
