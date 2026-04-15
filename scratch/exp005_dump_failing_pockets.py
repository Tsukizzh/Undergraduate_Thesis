"""
Pure data dump: for each selected dock, print the real content.
No scoring, no inference. Just facts.

For each dock: enzyme_id, uniprot, seq_len, full sequence (head+tail),
all pocket residues in order with seq[resid-1] comparison.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
COMPLEX_DIR = BASE / "data/structure/complex"
ENZ_CSV = BASE / "data/Enzymes.csv"
SPLIT_CSVS = [
    BASE / "data/splits/random/training_datas_0_pt.csv",
    BASE / "data/splits/random/val_datas_0_pt.csv",
    BASE / "data/splits/random/testing_datas_0_pt.csv",
]

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

with open(ENZ_CSV, encoding="utf-8-sig") as f:
    ENZ_ROWS = list(csv.DictReader(f))

# Build DOCK_TO_ENZYME and also DOCK_TO_SUBSTRATE
DOCK_INFO = {}  # dock -> (enzyme_id, substrate_id, label, which_split)
for split_file in SPLIT_CSVS:
    split_name = split_file.stem
    with open(split_file, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            d = int(row["Dock Index"])
            if d not in DOCK_INFO:
                DOCK_INFO[d] = (
                    int(row["Enzyme Index"]),
                    int(row["Substrate Index"]),
                    int(row["Label"]),
                    split_name,
                )

BACKBONE = {"N", "CA", "C", "O"}


def parse_pocket_residues(pdb_path):
    """Return ordered list of (resid, resname, chain). All kept standard AA with full backbone."""
    residue_atoms = {}  # (chain, resseq, resname) -> set of atom names
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26]
        if icode != " ":
            continue
        if resname not in THREE_TO_ONE:
            continue
        residue_atoms.setdefault((chain, resseq, resname), set()).add(atom_name)

    kept = []
    for (chain, resseq, resname), atoms in residue_atoms.items():
        if BACKBONE.issubset(atoms):
            kept.append((resseq, resname, chain))
    kept.sort()
    return kept


def get_complex_first_atom_line(dock):
    """Get first ATOM line from complex PDB for source inspection."""
    cpath = COMPLEX_DIR / f"{dock}.pdb"
    if not cpath.exists():
        return f"(complex/{dock}.pdb not found)"
    with open(cpath) as f:
        for line in f:
            if line.startswith("ATOM"):
                return line.rstrip()
    return "(no ATOM line)"


def dump_dock(dock, label):
    print()
    print("=" * 78)
    print(f"DUMP: {label}  (dock={dock})")
    print("=" * 78)

    if dock not in DOCK_INFO:
        print(f"  dock {dock} NOT in any split CSV")
        return
    enz_id, sub_id, lbl, split = DOCK_INFO[dock]
    print(f"  split:       {split}")
    print(f"  enzyme_id:   {enz_id}")
    print(f"  substrate_id:{sub_id}")
    print(f"  label:       {lbl}")

    row = ENZ_ROWS[enz_id]
    uniprot = row["uniprots"]
    seq = row["Protein sequence"].strip()
    print(f"  uniprot:     {uniprot}")
    print(f"  seq_len:     {len(seq)}")
    print(f"  seq[0:60]:   {seq[:60]}")
    print(f"  seq[-60:]:   {seq[-60:]}")

    # Also check complex PDB first ATOM line
    print(f"  complex first ATOM: {get_complex_first_atom_line(dock)[:100]}")

    pocket_path = POCKET_DIR / f"{dock}.pdb"
    if not pocket_path.exists():
        print(f"  pocket/{dock}.pdb not found!")
        return
    kept = parse_pocket_residues(pocket_path)
    print(f"  pocket residues kept: {len(kept)}")
    if not kept:
        return

    print()
    print(f"  {'pocket_resid':>12}  {'pocket_aa':>10}  {'seq[resid-1]':>14}  {'match':>6}")
    print(f"  {'-'*12}  {'-'*10}  {'-'*14}  {'-'*6}")
    for resid, resname, chain in kept:
        expected = THREE_TO_ONE[resname]
        if 0 <= resid - 1 < len(seq):
            actual = seq[resid - 1]
        else:
            actual = "(out)"
        mark = "OK" if expected == actual else "DIFF"
        print(f"  {resid:>12}  {resname} ({expected}){'':>3}  seq[{resid-1:>4}]={actual}{'':>4}  {mark}")


# 3 failing + 2 passing + 1 edge for reference
FAILING = [30789, 985, 24]
PASSING = [3, 8444]
EDGE = [19432]  # min_N case

for d in FAILING:
    dump_dock(d, "FAILING")
for d in PASSING:
    dump_dock(d, "PASSING (control)")
for d in EDGE:
    dump_dock(d, "EDGE N=1")
'''


def main():
    import subprocess, sys
    print("[local driver] dumping failing pocket data from server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=120,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
