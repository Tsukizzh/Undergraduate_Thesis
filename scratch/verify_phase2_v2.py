"""Phase 2 verification v2: understand flatbin format then do byte comparison."""
import torch
import numpy as np
import lmdb
import pickle

BASE = "/root/autodl-tmp/EZSpecificity/PathC/P450"
GROVER_FIXED = f"{BASE}/data/features/grover_fingerprint_fixed.lmdb"
SUB_DIR = f"{BASE}/data/pt_cache_allfix_shared/substrates"

idx = torch.load(f"{SUB_DIR}/substrates_index.pt", map_location="cpu", weights_only=False)
meta = torch.load(f"{SUB_DIR}/substrates_meta.pt", map_location="cpu", weights_only=False)

# Inspect idx['index']
inner_idx = idx['index']
print(f"idx['index'] type: {type(inner_idx).__name__}")
if isinstance(inner_idx, dict):
    sample_keys = list(inner_idx.keys())[:5]
    print(f"  first 5 keys: {sample_keys}")
    first_val = inner_idx[sample_keys[0]]
    print(f"  first value type: {type(first_val).__name__}, value: {first_val}")
    # Try to see how it maps substrate_id -> row/offset
    all_keys = sorted(int(k) if isinstance(k, (int, np.integer)) else int(k) for k in inner_idx.keys())
elif isinstance(inner_idx, (list, tuple)):
    print(f"  len: {len(inner_idx)}")
    print(f"  first 5: {inner_idx[:5]}")

# All keys
print(f"\nidx['index'] covered keys: {len(all_keys)}")
print(f"  min: {all_keys[0]}, max: {all_keys[-1]}")
print(f"  has 0: {0 in all_keys}")
print(f"  has 7: {7 in all_keys}")
print(f"  has 8: {8 in all_keys}  <- should be False")
print(f"  has 9: {9 in all_keys}")
print(f"  has 2124: {2124 in all_keys}")

expected = set(range(8)) | set(range(9, 2125))
missing = expected - set(all_keys)
extra = set(all_keys) - expected
print(f"  missing from expected: {len(missing)}")
print(f"  extra not in expected: {len(extra)}")

assert len(missing) == 0, f"missing substrate_ids: {sorted(missing)[:10]}"
assert len(extra) == 0, f"unexpected substrate_ids: {sorted(extra)[:10]}"
print("\n[CHECK] substrates_index.pt covers exactly {0..7, 9..2124}: PASS")

# Byte-level comparison: read flatbin for a few substrates and compare to LMDB
env = lmdb.open(GROVER_FIXED, subdir=False, readonly=True, lock=False)

# Read grover atom flatbin
grover_atom_dim = idx['grover_atom_dim']
print(f"\ngrover_atom_dim: {grover_atom_dim}")
print(f"format: {idx.get('format')}, dtype: {idx.get('dtype')}, endianness: {idx.get('endianness')}")

bin_path = f"{SUB_DIR}/substrates_grover.bin"
bin_size = open(bin_path, "rb").seek(0, 2)
print(f"bin file size: {bin_size}")

# Figure out how to read.
# format says idx['index'][sub_id] -> (offset_bytes, n_atoms) or (row, n_atoms) ?
test_keys = [0, 7, 9, 100, 500, 1000, 2000, 2124]
print("\nInspecting idx['index'] for test keys:")
for k in test_keys:
    v = inner_idx.get(k, None)
    print(f"  [{k}]: {v}")
