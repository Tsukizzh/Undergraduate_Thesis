#!/usr/bin/env python3
"""Convert clean PDBs to PDBQT using local MGLTools. Split protein/HEM, convert, merge."""

import os, subprocess, tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

CLEAN_DIR = Path(r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathC_2026-03-19_P450专属模型训练\C2_P450数据集构建\data\structures\receptors_clean")
PDBQT_DIR = Path(r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathC_2026-03-19_P450专属模型训练\C2_P450数据集构建\data\structures\receptors_pdbqt")
MGLTOOLS_PYTHON = r"D:\autodock\MGLTools-1.5.7\python.exe"
PREPARE_RECEPTOR = r"D:\autodock\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py"

HEM_ATOM_TYPES = {
    "FE": ("FE", 0.400),
    "NA": ("NA", -0.180), "NB": ("NA", -0.180), "NC": ("NA", -0.180), "ND": ("NA", -0.180),
    "CHA": ("A", 0.100), "CHB": ("A", 0.100), "CHC": ("A", 0.100), "CHD": ("A", 0.100),
    "C1A": ("A", 0.060), "C2A": ("A", -0.010), "C3A": ("A", 0.060), "C4A": ("A", 0.060),
    "C1B": ("A", 0.060), "C2B": ("A", -0.010), "C3B": ("A", 0.060), "C4B": ("A", 0.060),
    "C1C": ("A", 0.060), "C2C": ("A", -0.010), "C3C": ("A", 0.060), "C4C": ("A", 0.060),
    "C1D": ("A", 0.060), "C2D": ("A", -0.010), "C3D": ("A", 0.060), "C4D": ("A", 0.060),
    "CMA": ("C", -0.110), "CMB": ("C", -0.110), "CMC": ("C", -0.110), "CMD": ("C", -0.110),
    "CAA": ("C", -0.010), "CBA": ("C", -0.010), "CGA": ("C", 0.270),
    "CAD": ("C", -0.010), "CBD": ("C", -0.010), "CGD": ("C", 0.270),
    "O1A": ("OA", -0.490), "O2A": ("OA", -0.490),
    "O1D": ("OA", -0.490), "O2D": ("OA", -0.490),
    "CAB": ("C", 0.000), "CBB": ("C", 0.000),
    "CAC": ("C", 0.000), "CBC": ("C", 0.000),
}


def convert_one(pdb_path):
    stem = pdb_path.stem
    out_pdbqt = PDBQT_DIR / f"{stem}.pdbqt"
    if out_pdbqt.exists():
        return stem, "skip", ""

    protein_lines = []
    hem_lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip().upper()
                if res_name in ("HEM", "HEC", "HEA", "HEB", "HEO"):
                    hem_lines.append(line)
                elif line.startswith("ATOM"):
                    if res_name == "MSE":
                        line = line[:17] + "MET" + line[20:]
                    protein_lines.append(line)
            elif line.startswith(("TER", "END")):
                protein_lines.append(line)

    if not protein_lines:
        return stem, "fail", "no protein atoms"

    # Use a temp dir without Chinese chars
    tmp_dir = Path(tempfile.gettempdir()) / "pdbqt_conv"
    tmp_dir.mkdir(exist_ok=True)
    tmp_pdb = tmp_dir / f"{stem}.pdb"
    tmp_pdbqt = tmp_dir / f"{stem}.pdbqt"

    with open(tmp_pdb, "w") as f:
        f.writelines(protein_lines)

    try:
        proc = subprocess.run(
            [MGLTOOLS_PYTHON, PREPARE_RECEPTOR, "-r", str(tmp_pdb), "-o", str(tmp_pdbqt),
             "-A", "hydrogens", "-U", "nphs_lps_waters"],
            capture_output=True, text=True, timeout=120
        )
        if proc.returncode != 0 or not tmp_pdbqt.exists():
            return stem, "fail", (proc.stderr or proc.stdout or "unknown")[:200]
    except Exception as e:
        return stem, "fail", str(e)[:200]
    finally:
        tmp_pdb.unlink(missing_ok=True)

    with open(tmp_pdbqt) as f:
        pdbqt_lines = f.readlines()
    tmp_pdbqt.unlink(missing_ok=True)

    hem_pdbqt = []
    serial = 90000
    for line in hem_lines:
        atom_name = line[12:16].strip().upper()
        x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        atype, charge = HEM_ATOM_TYPES.get(atom_name, ("C", 0.000))
        hem_pdbqt.append(
            f"HETATM{serial:>5d}  {atom_name:<3s} HEM A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00    {charge:>+6.3f} {atype}\n"
        )
        serial += 1

    with open(out_pdbqt, "w") as f:
        for line in pdbqt_lines:
            if line.startswith("END"):
                continue
            f.write(line)
        for line in hem_pdbqt:
            f.write(line)
        f.write("END\n")

    return stem, "ok", ""


if __name__ == "__main__":
    PDBQT_DIR.mkdir(parents=True, exist_ok=True)
    todo = [p for p in sorted(CLEAN_DIR.glob("*.pdb")) if not (PDBQT_DIR / f"{p.stem}.pdbqt").exists()]
    print(f"Need to convert: {len(todo)}")

    ok = fail = 0
    with ProcessPoolExecutor(max_workers=6) as executor:
        futures = {executor.submit(convert_one, p): p for p in todo}
        for fut in as_completed(futures):
            stem, status, err = fut.result()
            if status == "ok":
                ok += 1
            elif status == "fail":
                fail += 1
                if fail <= 5:
                    print(f"  FAIL {stem}: {err}")
            if (ok + fail) % 50 == 0 and (ok + fail) > 0:
                print(f"  Progress: {ok+fail}/{len(todo)}: {ok} ok, {fail} fail")

    print(f"\nDone: {ok} ok, {fail} fail")
    print(f"Total PDBQTs: {len(list(PDBQT_DIR.glob('*.pdbqt')))}")
