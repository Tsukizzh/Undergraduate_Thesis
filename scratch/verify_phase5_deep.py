"""Phase 5 deep: inspect sample objects and verify feature_dim + GROVER byte-equality."""
import sys
import os
import torch
import numpy as np
import lmdb
import pickle
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
GROVER_FIXED = BASE / "data/features/grover_fingerprint_fixed.lmdb"


def setup_exp_path(exp_name: str) -> None:
    exp_dir = BASE / "experiments" / exp_name
    sys.path = [p for p in sys.path if 'experiments' not in p and 'EZSpecificity' not in p]
    sys.path.insert(0, str(exp_dir / "scripts"))
    sys.path.insert(0, str(exp_dir / "src"))
    for mod_name in list(sys.modules.keys()):
        if 'pt_dataset' in mod_name or mod_name.startswith('Datasets') or mod_name.startswith('Models') or mod_name == 'build_pt_cache':
            del sys.modules[mod_name]


def load_grover(sub_id: int) -> dict:
    env = lmdb.open(str(GROVER_FIXED), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        raw = txn.get(str(sub_id).encode())
    env.close()
    return pickle.loads(raw) if raw else None


def deep_test(cache_rel: str, exp_name: str, expected_feat_dim: int):
    print(f"\n=== {cache_rel} (expects feat_dim={expected_feat_dim}) ===")
    setup_exp_path(exp_name)
    from pt_dataset import PtCacheDataset

    ds = PtCacheDataset(cache_dir=BASE / "data" / cache_rel, split='train', edge_mode='fixed')
    sample = ds[0]

    # Inspect sample structure
    print(f"  sample type: {type(sample).__name__}")
    if hasattr(sample, 'keys'):
        keys = list(sample.keys()) if callable(sample.keys) else sample.keys
        print(f"  sample keys: {list(keys)[:15]}")

    # PyG Data object: has attributes like x, pos, edge_index, etc.
    all_attrs = [a for a in dir(sample) if not a.startswith('_') and not callable(getattr(sample, a, None))]
    tensor_attrs = {}
    for a in all_attrs:
        try:
            v = getattr(sample, a)
            if torch.is_tensor(v):
                tensor_attrs[a] = tuple(v.shape)
        except Exception:
            pass
    print(f"  tensor attributes (name -> shape):")
    for name, shape in sorted(tensor_attrs.items()):
        print(f"    {name}: {shape}")

    # Find substrate_id / sub_id
    sub_id = None
    for k in ['substrate_id', 'sub_id', 'reaction_id', 'substrate_ids', 'substrate_global_id']:
        if hasattr(sample, k):
            v = getattr(sample, k)
            if torch.is_tensor(v):
                sub_id = int(v.item()) if v.numel() == 1 else int(v.flatten()[0])
            else:
                sub_id = int(v)
            print(f"  sub_id via .{k}: {sub_id}")
            break

    # Find the feature dim from likely attributes
    for k in ['protein_feat', 'protein_atom_feature', 'atom_feature', 'x_p', 'node_feat']:
        if hasattr(sample, k):
            v = getattr(sample, k)
            if torch.is_tensor(v) and len(v.shape) >= 2:
                dim = v.shape[-1]
                match = dim == expected_feat_dim
                print(f"  feat_dim via .{k}: {dim} (expected {expected_feat_dim}) {'OK' if match else 'MISMATCH'}")

    # GROVER mean comparison: look for likely fields
    for k in ['grover_mean', 'grover_total', 'substrate_grover_mean']:
        if hasattr(sample, k):
            v = getattr(sample, k)
            if torch.is_tensor(v):
                sample_arr = v.detach().cpu().numpy().reshape(-1).astype(np.float16)
                if sub_id is not None:
                    lmdb_data = load_grover(sub_id)
                    lmdb_total = np.asarray(lmdb_data['total_embedding'], dtype=np.float16)
                    min_len = min(len(sample_arr), len(lmdb_total))
                    match = np.array_equal(sample_arr[:min_len], lmdb_total[:min_len])
                    print(f"  .{k}[:50] vs lmdb[{sub_id}].total_embedding[:50]: match={match} "
                          f"(sample.shape={v.shape}, lmdb.shape={lmdb_total.shape})")
                    if not match:
                        diff = (sample_arr[:10] - lmdb_total[:10])
                        print(f"    first-10 diff: {diff}")
            break


CACHES = [
    ('pt_cache_allfix/random',          'EXP001_fixed', 28),
    ('pt_cache_heme_allfix/random',     'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix',            'EXP003_fixed', 37),
    ('pt_cache_allfix_unified/random',  'EXP001_fixed', 28),
    ('pt_cache_heme_allfix_unified/random', 'EXP002a_fixed', 31),
    ('pt_cache_geom_allfix_unified',    'EXP003_fixed', 37),
]

for cache_rel, exp, dim in CACHES:
    try:
        deep_test(cache_rel, exp, dim)
    except Exception as e:
        print(f"  ERROR: {type(e).__name__}: {e}")
        import traceback; traceback.print_exc()
