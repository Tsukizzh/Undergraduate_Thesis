"""
For each test dock, read its complex PDB, compute best delta offset
that aligns structure's residues with Enzymes.csv sequence.

Small test first: 5 known cases (3 passing + 3 failing + 1 more PCPD for safety).
Report per-dock delta and final match rate after fix.

If delta approach works universally on these 6, proceed to all enzymes.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
from pathlib import Path
from collections import defaultdict

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
COMPLEX_DIR = BASE / "data/structure/complex"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
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
ENZ_SEQ = [r["Protein sequence"].strip() for r in ENZ_ROWS]

DOCK_TO_ENZYME = {}
for p in SPLIT_CSVS:
    with open(p, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            DOCK_TO_ENZYME[int(row["Dock Index"])] = int(row["Enzyme Index"])


def parse_pdb_residues(pdb_path, is_complex=False):
    """Return sorted list of (resid, one_letter_aa) from protein ATOM lines."""
    residue_names = {}  # (chain, resid) -> resname (first seen)
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
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
        key = (chain, resseq)
        if key not in residue_names:
            residue_names[key] = resname

    # Check chains
    chains = set(c for c, _ in residue_names.keys())
    if len(chains) > 1:
        # Multi-chain: take the first chain
        first_chain = sorted(chains)[0]
        residue_names = {k: v for k, v in residue_names.items() if k[0] == first_chain}

    result = []
    for (c, rseq), rname in residue_names.items():
        result.append((rseq, THREE_TO_ONE[rname]))
    result.sort()
    return result


def find_best_delta(residues, seq, delta_range=200):
    """Return (best_delta, best_match, best_total) where best_match / best_total is
    the fraction of residues that align when using seq[resid + delta - 1]."""
    seq_len = len(seq)
    best_delta = 0
    best_match = -1
    best_total = 0
    for delta in range(-delta_range, delta_range + 1):
        matches = 0
        total = 0
        for resid, aa in residues:
            idx = resid + delta - 1
            if 0 <= idx < seq_len:
                total += 1
                if seq[idx] == aa:
                    matches += 1
        if total == 0:
            continue
        if matches > best_match or (matches == best_match and total > best_total):
            best_match = matches
            best_total = total
            best_delta = delta
    return best_delta, best_match, best_total


# Test cases
TEST_DOCKS = {
    3:     "PASSING (alphafill)",
    8444:  "PASSING (pcpd_predicted)",
    34:    "PASSING (short_seq)",
    24:    "FAILING (experimental_pdb)",
    985:   "FAILING",
    30789: "FAILING",
}

print("=" * 78)
print("Delta detection on test docks")
print("=" * 78)

for dock, label in TEST_DOCKS.items():
    print(f"\n--- dock={dock}  [{label}] ---")

    if dock not in DOCK_TO_ENZYME:
        print(f"  NOT in any split")
        continue
    enz_id = DOCK_TO_ENZYME[dock]
    seq = ENZ_SEQ[enz_id]
    print(f"  enzyme_id={enz_id}, seq_len={len(seq)}")
    print(f"  seq[0:5]={seq[:5]}")

    # 1. Parse complex PDB (has more residues)
    cpath = COMPLEX_DIR / f"{dock}.pdb"
    if cpath.exists():
        complex_res = parse_pdb_residues(cpath, is_complex=True)
        print(f"  complex residues: {len(complex_res)} (first: resid={complex_res[0][0]} {complex_res[0][1]}, last: resid={complex_res[-1][0]} {complex_res[-1][1]})")
        delta_c, match_c, total_c = find_best_delta(complex_res, seq)
        print(f"  COMPLEX delta detection: delta={delta_c:+d}, match={match_c}/{total_c} ({100*match_c/total_c:.1f}%)")
    else:
        print(f"  complex/{dock}.pdb NOT FOUND")
        complex_res = []

    # 2. Parse pocket PDB
    ppath = POCKET_DIR / f"{dock}.pdb"
    if ppath.exists():
        pocket_res = parse_pdb_residues(ppath)
        print(f"  pocket residues: {len(pocket_res)}")
        delta_p, match_p, total_p = find_best_delta(pocket_res, seq)
        print(f"  POCKET  delta detection: delta={delta_p:+d}, match={match_p}/{total_p} ({100*match_p/total_p:.1f}%)")
    else:
        print(f"  pocket/{dock}.pdb NOT FOUND")
'''


def main():
    import subprocess, sys
    print("[local driver] running delta detection on server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=300,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
