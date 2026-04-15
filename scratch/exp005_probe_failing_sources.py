"""Look up source type + PDB ID for the 3 failing-case enzymes (54, 838, 952)."""
import csv
from pathlib import Path

BASE_C2 = Path("d:/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C2_P450数据集构建/data")

with open(BASE_C2 / "structures/manifests/receptor_manifest.csv", encoding="utf-8-sig") as f:
    man = list(csv.DictReader(f))
# key by gid
by_gid = {r["global_enzyme_id"]: r for r in man}

# enzyme_id -> ENZ_G{enzyme_id+1:06d}
for enz_id in [54, 838, 952, 93, 659]:
    gid = f"ENZ_G{enz_id+1:06d}"
    r = by_gid.get(gid, {})
    print(f"enzyme_id={enz_id} ({gid})")
    print(f"  uniprot:     {r.get('canonical_uniprot_id','?')}")
    print(f"  source:      {r.get('structure_source','?')}")
    print(f"  status:      {r.get('status','?')}")
    print(f"  existing_pdb_path: {r.get('existing_pdb_path','?')}")
    print(f"  alphafill_id: {r.get('alphafill_id','?')}")
    print(f"  heme_transplant_pdb_path: {r.get('heme_transplant_pdb_path','?')}")
    print(f"  alphafold_pdb_path: {r.get('alphafold_pdb_path','?')}")
    print()
