"""
Probe receptor_manifest.csv locally:
- What's the column format?
- For experimental_pdb rows: how is the PDB ID stored?
- How does global_enzyme_id map to the P450/data/Enzymes.csv enzyme_id (integer row index)?

We need enzyme_id (int row in Enzymes.csv) -> (structure_source, pdb_id if experimental).
"""
from __future__ import annotations
import csv
from collections import Counter
from pathlib import Path

BASE_C2 = Path("d:/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C2_P450数据集构建/data")
MANIFEST = BASE_C2 / "structures/manifests/receptor_manifest.csv"
GLOBAL_ENZ = BASE_C2 / "combined/global_enzymes.csv"

with open(MANIFEST, encoding="utf-8-sig") as f:
    man = list(csv.DictReader(f))
print(f"manifest rows: {len(man)}")
print(f"columns (first 30): {list(man[0].keys())[:30]}")
print(f"columns (all): {list(man[0].keys())}")
print()

# Count by status / structure_source
status_c = Counter(r.get("status", "") for r in man)
src_c = Counter(r.get("structure_source", "") for r in man)

print("status distribution:")
for s, c in status_c.most_common():
    print(f"  {s!r:45s}: {c}")
print()
print("structure_source distribution:")
for s, c in src_c.most_common():
    print(f"  {s!r:45s}: {c}")
print()

# Find the PDB-ID field for experimental_pdb rows
exp = [r for r in man if "experimental" in (r.get("status", "") + r.get("structure_source", "")).lower()]
print(f"experimental_pdb-ish rows: {len(exp)}")
if exp:
    r0 = exp[0]
    print("first experimental row non-empty fields:")
    for k, v in r0.items():
        if v and str(v).strip() and str(v) != "nan":
            print(f"  {k:40s} = {str(v)[:120]}")
    print()
    # Print 10 examples of existing_pdb_path and any *pdb_id* / PDB-looking field
    print("first 10 experimental rows with key id fields:")
    for r in exp[:10]:
        gid = r.get("global_enzyme_id", "")
        uid = r.get("canonical_uniprot_id", "")
        ep = r.get("existing_pdb_path", "")
        src = r.get("structure_source", "")
        st = r.get("status", "")
        print(f"  gid={gid} uid={uid} src={src!r} st={st!r} ep={ep!r}")

# How does global_enzyme_id relate to Enzymes.csv row index?
# Load global_enzymes.csv and see
print()
print("=" * 60)
with open(GLOBAL_ENZ, encoding="utf-8-sig") as f:
    genz = list(csv.DictReader(f))
print(f"global_enzymes rows: {len(genz)}")
print(f"columns: {list(genz[0].keys())}")
print()
print("first 3 rows:")
for r in genz[:3]:
    for k, v in r.items():
        if v and str(v).strip() != "":
            print(f"  {k:35s} = {str(v)[:100]}")
    print()
