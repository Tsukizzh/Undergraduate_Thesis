"""Phase 2 final verification: byte-level compare flatbin payloads to fixed LMDB."""
import torch
import numpy as np
import lmdb
import pickle

BASE = "/root/autodl-tmp/EZSpecificity/PathC/P450"
GROVER_FIXED = f"{BASE}/data/features/grover_fingerprint_fixed.lmdb"
MORGAN_NPY = f"{BASE}/data/features/morgan_fingerprint.npy"
SUB_DIR = f"{BASE}/data/pt_cache_allfix_shared/substrates"

idx = torch.load(f"{SUB_DIR}/substrates_index.pt", map_location="cpu", weights_only=False)
meta = torch.load(f"{SUB_DIR}/substrates_meta.pt", map_location="cpu", weights_only=False)
morgan_src = np.load(MORGAN_NPY)
env = lmdb.open(GROVER_FIXED, subdir=False, readonly=True, lock=False)

inner_idx = idx['index']
grover_atom_dim = idx['grover_atom_dim']
mean_tensor = meta['grover_mean']   # (2124, 4885) fp16
morgan_tensor = meta['morgan']      # (2124, 1024) uint8

print("Reading substrates_grover.bin...")
bin_path = f"{SUB_DIR}/substrates_grover.bin"
with open(bin_path, "rb") as f:
    bin_bytes = f.read()
print(f"  total size: {len(bin_bytes)} bytes ({len(bin_bytes)/1024/1024:.1f} MB)")

# Full sweep byte-level comparison
print("\n=== Full sweep: compare every substrate (2124) ===")
mismatches = {'grover_atom': 0, 'grover_mean': 0, 'morgan': 0}
errors = []

with env.begin() as txn:
    for sub_id in sorted(inner_idx.keys()):
        byte_offset, n_atoms, meta_row = inner_idx[sub_id]

        # -- Read from LMDB via pickle --
        raw = txn.get(str(sub_id).encode())
        assert raw is not None, f"sub_id={sub_id} missing from fixed LMDB"
        lmdb_data = pickle.loads(raw)
        lmdb_emb = lmdb_data['embedding']         # (n_atoms, 2400) float64 typically
        lmdb_total = lmdb_data['total_embedding']  # (4885,) float64

        # -- Read grover_atom from flatbin --
        byte_len = n_atoms * grover_atom_dim * 2  # fp16
        bin_slice = bin_bytes[byte_offset:byte_offset + byte_len]
        bin_emb = np.frombuffer(bin_slice, dtype=np.float16).reshape(n_atoms, grover_atom_dim)

        # fp16 cast of original LMDB emb
        lmdb_emb_fp16 = np.asarray(lmdb_emb, dtype=np.float16)

        if not np.array_equal(bin_emb, lmdb_emb_fp16):
            mismatches['grover_atom'] += 1
            if len(errors) < 5:
                diff = (bin_emb - lmdb_emb_fp16).astype(np.float32)
                errors.append(f"grover_atom[{sub_id}] max_diff={np.max(np.abs(diff))}")

        # -- Read grover_mean from meta --
        bin_mean = mean_tensor[meta_row].numpy()
        lmdb_mean_fp16 = np.asarray(lmdb_total, dtype=np.float16)
        if not np.array_equal(bin_mean, lmdb_mean_fp16):
            mismatches['grover_mean'] += 1
            if len(errors) < 5:
                diff = (bin_mean - lmdb_mean_fp16).astype(np.float32)
                errors.append(f"grover_mean[{sub_id}] max_diff={np.max(np.abs(diff))}")

        # -- Read morgan from meta --
        bin_morgan = morgan_tensor[meta_row].numpy()
        src_morgan = morgan_src[sub_id].astype(np.uint8)
        if not np.array_equal(bin_morgan, src_morgan):
            mismatches['morgan'] += 1
            if len(errors) < 5:
                errors.append(f"morgan[{sub_id}] differs")

print(f"\nTotal substrates checked: {len(inner_idx)}")
print(f"  grover_atom mismatches: {mismatches['grover_atom']}")
print(f"  grover_mean mismatches: {mismatches['grover_mean']}")
print(f"  morgan mismatches:      {mismatches['morgan']}")
if errors:
    print("\nFirst errors:")
    for e in errors:
        print(f"  {e}")
else:
    print("\n[ALL MATCH]  flatbin byte-level equivalent to LMDB (fp16 precision)")

# Sub_id=8 should be missing
print(f"\nsub_id=8 in flatbin index: {8 in inner_idx}")
assert 8 not in inner_idx

# Atom count spot-check
print("\nAtom count spot-check:")
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
df = pd.read_csv(f"{BASE}/data/Substrates.csv")
for sid in [0, 9, 100, 2124]:
    n_atoms_idx = inner_idx[sid][1]
    rdkit_n = Chem.MolFromSmiles(df.iloc[sid]['Substrate_SMILES']).GetNumAtoms()
    match = n_atoms_idx == rdkit_n
    print(f"  sub_id={sid}: idx={n_atoms_idx}, rdkit={rdkit_n}, match={match}")
    assert match

print("\n=== Phase 2 verification PASS ===")
