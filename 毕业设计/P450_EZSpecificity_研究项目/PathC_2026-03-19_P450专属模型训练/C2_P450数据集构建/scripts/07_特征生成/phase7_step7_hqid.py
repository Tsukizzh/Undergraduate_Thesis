"""Phase 7 Step 7: Generate high_quality_id.txt."""
from pathlib import Path

pocket_dir = Path("/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/pocket")
ligand_dir = Path("/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/ligand")
out_txt = Path("/root/rivermind-data/EZSpecificity/PathC/P450/high_quality_id.txt")

pocket_ids = {p.stem for p in pocket_dir.glob("*.pdb")}
ligand_ids = {p.stem for p in ligand_dir.glob("*.sdf")}
ids = sorted(pocket_ids & ligand_ids, key=lambda x: int(x))

with open(out_txt, "w") as f:
    for x in ids:
        f.write(f"{x}\n")

print(f"high_quality_ids={len(ids)}")
