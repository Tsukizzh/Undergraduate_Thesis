"""Phase 2 verification: compare substrates flatbin content to fixed LMDB.

For each substrate_id in [0..7, 9..2124]:
1. Read from fixed LMDB directly (via pickle.loads)
2. Read from flatbin via substrates_index.pt offsets
3. Byte-level compare grover_atom embedding and grover_mean / morgan
4. Verify sub_idx=8 is absent from the flatbin index
"""
import torch
import numpy as np
import lmdb
import pickle

BASE = "/root/autodl-tmp/EZSpecificity/PathC/P450"
GROVER_FIXED = f"{BASE}/data/features/grover_fingerprint_fixed.lmdb"
MORGAN = f"{BASE}/data/features/morgan_fingerprint.npy"
SUB_DIR = f"{BASE}/data/pt_cache_allfix_shared/substrates"

# Load flatbin artifacts
idx = torch.load(f"{SUB_DIR}/substrates_index.pt", map_location="cpu", weights_only=False)
meta = torch.load(f"{SUB_DIR}/substrates_meta.pt", map_location="cpu", weights_only=False)

print("substrates_index.pt keys:", list(idx.keys()) if isinstance(idx, dict) else type(idx).__name__)
print("substrates_meta.pt keys:", list(meta.keys()) if isinstance(meta, dict) else type(meta).__name__)

# Inspect format of the index
if isinstance(idx, dict):
    for k, v in idx.items():
        if hasattr(v, 'shape'):
            print(f"  idx[{k}]: shape={tuple(v.shape)} dtype={v.dtype}")
        elif hasattr(v, '__len__'):
            print(f"  idx[{k}]: len={len(v)}")

if isinstance(meta, dict):
    for k, v in meta.items():
        if hasattr(v, 'shape'):
            print(f"  meta[{k}]: shape={tuple(v.shape)} dtype={v.dtype}")
        elif hasattr(v, '__len__'):
            print(f"  meta[{k}]: len={len(v)}")

# Extract substrate_id list from the index
# Format is likely {substrate_global_id: (offset, n_atoms, meta_row)} or similar
if 'substrate_ids' in idx:
    sub_ids_in_flatbin = set(int(x) for x in idx['substrate_ids'].tolist())
    sub_id_to_row = {int(sid): row for row, sid in enumerate(idx['substrate_ids'].tolist())}
elif 'keys' in idx:
    sub_ids_in_flatbin = set(int(x) for x in idx['keys'])
    sub_id_to_row = None
else:
    # Assume idx itself is a dict {sub_id: ...}
    sub_ids_in_flatbin = set(int(k) for k in idx.keys())
    sub_id_to_row = None

print(f"\nSubstrate count in flatbin: {len(sub_ids_in_flatbin)}")
print(f"  min: {min(sub_ids_in_flatbin)}, max: {max(sub_ids_in_flatbin)}")
print(f"  has 0: {0 in sub_ids_in_flatbin}")
print(f"  has 7: {7 in sub_ids_in_flatbin}")
print(f"  has 8: {8 in sub_ids_in_flatbin}  <- should be False")
print(f"  has 9: {9 in sub_ids_in_flatbin}")
print(f"  has 2124: {2124 in sub_ids_in_flatbin}")

expected = set(range(8)) | set(range(9, 2125))
missing = expected - sub_ids_in_flatbin
extra = sub_ids_in_flatbin - expected
print(f"  missing from expected: {len(missing)} (first 5: {sorted(missing)[:5]})")
print(f"  extra not in expected: {len(extra)} (first 5: {sorted(extra)[:5]})")

# Byte-level comparison for a few substrates
print("\n--- Byte comparison: flatbin vs fixed LMDB ---")
env = lmdb.open(GROVER_FIXED, subdir=False, readonly=True, lock=False)

# Raw flatbin read
bin_path = f"{SUB_DIR}/substrates_grover.bin"
bin_fp = open(bin_path, "rb")
try:
    bin_size = bin_fp.seek(0, 2)
    bin_fp.seek(0)
    print(f"substrates_grover.bin size: {bin_size/1024/1024:.1f} MB")

    # We need to decode the flatbin. Its format depends on substrates_index.pt.
    # If idx has {'offsets', 'lengths'} or similar, we can slice.
    print()
    print("Index has fields:", list(idx.keys()) if isinstance(idx, dict) else '-')
except Exception as e:
    print(f"failed to seek bin: {e}")
finally:
    bin_fp.close()
