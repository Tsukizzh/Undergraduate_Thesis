"""
Local preprocessing: build enzyme_meta.json for all 1622 enzymes.

For each enzyme:
  enzyme_id -> {
    "uniprot": str,
    "source": str (structure_source from manifest),
    "pdb_id": str or None (only for experimental_pdb),
  }

Saves to scratch/enzyme_meta.json, ready to scp to server.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

BASE_C2 = Path("d:/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C2_P450数据集构建/data")
MANIFEST = BASE_C2 / "structures/manifests/receptor_manifest.csv"

with open(MANIFEST, encoding="utf-8-sig") as f:
    rows = list(csv.DictReader(f))

meta = {}
for r in rows:
    gid = r.get("global_enzyme_id", "").strip()
    if not gid.startswith("ENZ_G"):
        continue
    enz_id = int(gid[5:]) - 1  # ENZ_G000001 -> row 0 in Enzymes.csv
    source = r.get("structure_source", "").strip() or r.get("status", "").strip() or "unknown"
    pdb_id = None
    if source == "experimental_pdb":
        ep = r.get("existing_pdb_path", "").strip()
        # "s1_rcsb/1JPZ" -> "1JPZ"
        if ep.startswith("s1_rcsb/"):
            pdb_id = ep.split("/", 1)[1].lower()
    meta[enz_id] = {
        "uniprot": r.get("canonical_uniprot_id", "").strip(),
        "source": source,
        "pdb_id": pdb_id,
    }

out_path = Path("d:/EZSpecificity_Project/scratch/enzyme_meta.json")
out_path.write_text(json.dumps(meta, indent=None), encoding="utf-8")
print(f"wrote {len(meta)} enzymes to {out_path}  ({out_path.stat().st_size} bytes)")

# Quick stats
from collections import Counter
c = Counter(m["source"] for m in meta.values())
for s, n in c.most_common():
    print(f"  {s:40s}: {n}")
exp_with_pdb = sum(1 for m in meta.values() if m["source"] == "experimental_pdb" and m["pdb_id"])
print(f"experimental_pdb with pdb_id: {exp_with_pdb}")
