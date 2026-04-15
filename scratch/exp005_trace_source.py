"""
Trace which structure source was used for failing vs passing docks.

The receptor_manifest.csv is keyed by global_enzyme_id (ENZ_Gxxxxxx).
Enzymes.csv rows have "uniprots" column with either the real uniprot OR ENZ_G*.
But the global_enzyme_id in receptor_manifest is ENZ_G*, so we need to join
via uniprot_id column in both.

Look at:
  receptor_manifest.csv has columns: global_enzyme_id, canonical_uniprot_id, structure_source, ...

Actually global_enzymes.csv probably has the full mapping. Let's find both.
"""
from __future__ import annotations
import csv
from pathlib import Path

BASE_C2 = Path("d:/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C2_P450数据集构建/data")
COMBINED = BASE_C2 / "combined"
MANIFEST = BASE_C2 / "structures/manifests/receptor_manifest.csv"

# Load receptor manifest
with open(MANIFEST, encoding="utf-8-sig") as f:
    manifest = list(csv.DictReader(f))
print(f"receptor_manifest rows: {len(manifest)}")
print(f"columns: {list(manifest[0].keys())[:10]}")
print()

# Load global_enzymes.csv
global_csv = COMBINED / "global_enzymes.csv"
if global_csv.exists():
    with open(global_csv, encoding="utf-8-sig") as f:
        g_rows = list(csv.DictReader(f))
    print(f"global_enzymes.csv rows: {len(g_rows)}")
    print(f"columns: {list(g_rows[0].keys())[:10]}")
else:
    print(f"global_enzymes.csv NOT FOUND at {global_csv}")
    import os
    print(f"files in {COMBINED}:")
    for f in os.listdir(COMBINED):
        print(f"  {f}")
print()

# Build mapping uniprot -> source
uni_to_info = {}
for r in manifest:
    u = r.get("canonical_uniprot_id", "").strip()
    g = r.get("global_enzyme_id", "").strip()
    src = r.get("structure_source", "").strip()
    if u:
        uni_to_info[u] = (g, src)

# Failing docks → find enzyme uniprot via Enzymes.csv
# Enzymes.csv is the cached local version
ENZ_CSV_LOCAL = Path("d:/EZSpecificity_Project/scratch/P450_Enzymes_server.csv")
with open(ENZ_CSV_LOCAL, encoding="utf-8-sig") as f:
    enz_rows = list(csv.DictReader(f))

# For dock_indices we need the split CSV. Let me pick the 3 failing enz_ids we know:
# dock=30789, 985, 24 → enzyme 93(NO, that was passing), 838, and we don't know rand_0 enz_id yet
# Let's load dock→enzyme from any split csv (not downloaded locally, so fetch)
# Actually, from previous output:
#   dock=3 → enzyme_id=93 (PASSING)
#   dock=8444 → enzyme_id=659 (PASSING)
#   dock=24 → enzyme_id=838 (FAILING)
# For dock=30789 and dock=985 we'd need the server

# Check the enzymes we know
KNOWN = {
    3: ("PASSING", 93),
    8444: ("PASSING", 659),
    24: ("FAILING", 838),
}

print("=" * 70)
print("Known dock → enzyme → structure source")
print("=" * 70)
for dock, (status, enz_id) in KNOWN.items():
    row = enz_rows[enz_id]
    uniprot = row["uniprots"].strip()
    seq = row["Protein sequence"].strip()
    src_info = uni_to_info.get(uniprot, ("(not in manifest)", "(none)"))
    print(f"dock={dock} [{status}]")
    print(f"  enzyme_id:  {enz_id}")
    print(f"  uniprot:    {uniprot}")
    print(f"  seq[:3]:    {seq[:3]}")
    print(f"  seq[:50]:   {seq[:50]}")
    print(f"  global_id:  {src_info[0]}")
    print(f"  source:     {src_info[1]}")
    print()

# Now: aggregate all seq[:1] (the first residue) by structure_source
# to understand which sources have structures that start with non-M residues
print("=" * 70)
print("First residue in Enzymes.csv seq, grouped by structure_source")
print("=" * 70)

# Build enz_id by uniprot lookup (both Enzymes.csv rows and manifest rows)
uni_to_enz_row = {}
for i, r in enumerate(enz_rows):
    u = r["uniprots"].strip()
    uni_to_enz_row[u] = (i, r["Protein sequence"].strip())

from collections import Counter
src_first_aa = {}  # structure_source -> Counter(first_aa)
for r in manifest:
    u = r.get("canonical_uniprot_id", "").strip()
    src = r.get("structure_source", "").strip() or "(empty)"
    if u not in uni_to_enz_row:
        continue
    _, seq = uni_to_enz_row[u]
    first_aa = seq[:1] if seq else ""
    src_first_aa.setdefault(src, Counter())[first_aa] += 1

for src, cnt in sorted(src_first_aa.items()):
    total = sum(cnt.values())
    m_count = cnt.get("M", 0)
    non_m_count = total - m_count
    print(f"\n  {src}: {total} enzymes")
    print(f"    start with M: {m_count}")
    print(f"    start with non-M: {non_m_count}")
    if non_m_count > 0:
        print(f"    top first-AA: {cnt.most_common(5)}")
