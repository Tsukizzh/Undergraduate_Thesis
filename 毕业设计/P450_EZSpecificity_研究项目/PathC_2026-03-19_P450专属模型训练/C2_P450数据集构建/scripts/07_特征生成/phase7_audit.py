"""Phase 7 Audit: count exactly how many pairs are usable for training."""
import lmdb, pandas as pd, numpy as np

P450 = "/root/rivermind-data/EZSpecificity/PathC/P450"

data = pd.read_csv(f"{P450}/data.csv")
enzymes = pd.read_csv(f"{P450}/Enzymes.csv")
substrates = pd.read_csv(f"{P450}/Substrates.csv")

pos_count = int((data["Label"] == 1).sum())
neg_count = int((data["Label"] == 0).sum())
print("=== INPUT ===")
print(f"Enzymes.csv: {len(enzymes)}")
print(f"Substrates.csv: {len(substrates)}")
print(f"data.csv: {len(data)} pairs (pos={pos_count}, neg={neg_count})")

with open(f"{P450}/high_quality_id.txt") as f:
    hq_ids = set(int(l.strip()) for l in f if l.strip())
print(f"\n=== STRUCTURE ===")
print(f"high_quality_id (pocket+ligand): {len(hq_ids)}")

env = lmdb.open(f"{P450}/structure/structure_features.lmdb", readonly=True, subdir=False, lock=False)
with env.begin() as txn:
    str_keys = set()
    for key, _ in txn.cursor(): str_keys.add(int(key.decode()))
env.close()
print(f"structure_features.lmdb entries: {len(str_keys)}")

env = lmdb.open(f"{P450}/enzyme_features.lmdb", readonly=True, subdir=False, lock=False)
with env.begin() as txn:
    esm_keys = set()
    for key, _ in txn.cursor(): esm_keys.add(int(key.decode()))
env.close()
missing_esm_reasons = []
for i, seq in enumerate(enzymes["Protein sequence"]):
    if i not in esm_keys:
        if len(seq) > 1000:
            missing_esm_reasons.append("too_long")
        else:
            missing_esm_reasons.append("non_standard_aa")

print(f"\n=== ESM ===")
print(f"enzyme_features.lmdb: {len(esm_keys)} / {len(enzymes)}")
from collections import Counter
rc = Counter(missing_esm_reasons)
for r, c in rc.most_common():
    print(f"  Missing: {c} ({r})")

env = lmdb.open(f"{P450}/reaction_features.lmdb", readonly=True, subdir=False, lock=False)
with env.begin() as txn: rxn_count = txn.stat()["entries"]
env.close()

env = lmdb.open(f"{P450}/grover_fingerprint.lmdb", readonly=True, subdir=False, lock=False)
with env.begin() as txn: grov_count = txn.stat()["entries"]
env.close()

morgan = np.load(f"{P450}/morgan_fingerprint.npy")

print(f"\n=== SUBSTRATE FEATURES ===")
print(f"reaction_features.lmdb: {rxn_count} / {len(substrates)}")
print(f"grover_fingerprint.lmdb: {grov_count} / {len(substrates)}")
print(f"morgan_fingerprint.npy: {morgan.shape[0]} / {len(substrates)}")

# Pair-level analysis
no_hq = no_str = no_esm = usable = 0
usable_pos = usable_neg = 0
enz_used = set()
sub_used = set()
loss_reasons = Counter()

for _, row in data.iterrows():
    dk = int(row["Dock Index"])
    ei = int(row["Enzyme Index"])
    si = int(row["Substrate Index"])
    label = int(row["Label"])

    if dk not in hq_ids:
        no_hq += 1
        loss_reasons["no_pocket_or_ligand"] += 1
        continue
    if dk not in str_keys:
        no_str += 1
        loss_reasons["no_structure_feature"] += 1
        continue
    if ei not in esm_keys:
        no_esm += 1
        loss_reasons["no_esm"] += 1
        continue
    usable += 1
    enz_used.add(ei)
    sub_used.add(si)
    if label == 1: usable_pos += 1
    else: usable_neg += 1

print(f"\n{'='*50}")
print(f"=== FINAL USABLE PAIRS ===")
print(f"{'='*50}")
print(f"Total in data.csv:        {len(data):>8}")
print(f"  - No pocket/ligand:     {no_hq:>8}")
print(f"  - No structure feature:  {no_str:>8}")
print(f"  - No ESM embedding:     {no_esm:>8}")
print(f"  ─────────────────────────────")
print(f"  = USABLE:               {usable:>8} ({usable/len(data)*100:.1f}%)")
print(f"    Positive:             {usable_pos:>8}")
print(f"    Negative:             {usable_neg:>8}")
ratio = usable_neg / max(usable_pos, 1)
print(f"    Pos:Neg ratio:          1:{ratio:.1f}")
print(f"    Unique enzymes:       {len(enz_used):>8} / {len(enzymes)}")
print(f"    Unique substrates:    {len(sub_used):>8} / {len(substrates)}")
