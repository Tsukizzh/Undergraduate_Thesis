"""Step3A: Validate docking environment (Meeko/RDKit/BioPython/Vina/MGLTools)."""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import tempfile
from pathlib import Path

LOG = logging.getLogger(__name__)

DEFAULT_VINA_EXE = Path(r"D:\autodock\vina.exe")
DEFAULT_MGL_PYTHON = Path(r"D:\autodock\MGLTools-1.5.7\python.exe")


def check_meeko_roundtrip() -> tuple[bool, str]:
    try:
        import meeko
        from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
    except Exception as exc:
        return False, f"Meeko import failed: {exc}"

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            return False, "RDKit embed failed"
        AllChem.MMFFOptimizeMolecule(mol)

        setups = MoleculePreparation().prepare(mol)
        pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            return False, f"PDBQT write failed: {err}"

        with tempfile.NamedTemporaryFile("w", suffix=".pdbqt", delete=False) as tmp:
            tmp.write(pdbqt_str)
            tmp_path = Path(tmp.name)
        try:
            pdbqt_mol = PDBQTMolecule.from_file(str(tmp_path))
            mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
            if not mols or mols[0] is None:
                return False, "PDBQT roundtrip produced no molecule"
        finally:
            tmp_path.unlink(missing_ok=True)

        return True, f"Meeko {meeko.__version__} roundtrip OK"
    except Exception as exc:
        return False, f"Meeko roundtrip error: {exc}"


def check_rdkit_3d() -> tuple[bool, str]:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            return False, "ETKDGv3 embedding failed"
        return True, "RDKit 3D embedding OK"
    except Exception as exc:
        return False, f"RDKit error: {exc}"


def check_biopython() -> tuple[bool, str]:
    try:
        from Bio.PDB import PDBParser
        PDBParser(QUIET=True)
        return True, "BioPython PDBParser OK"
    except Exception as exc:
        return False, f"BioPython error: {exc}"


def check_vina(exe: Path, timeout: int) -> tuple[bool, str]:
    if not exe.exists():
        return False, f"Not found: {exe}"
    try:
        proc = subprocess.run([str(exe), "--version"], capture_output=True, text=True, timeout=timeout)
        out = (proc.stdout or proc.stderr or "").strip()
        return (proc.returncode == 0, out or "Vina OK")
    except Exception as exc:
        return False, f"Vina error: {exc}"


def check_mgl(exe: Path, timeout: int) -> tuple[bool, str]:
    if not exe.exists():
        return False, f"Not found: {exe}"
    try:
        proc = subprocess.run([str(exe), "-c", "print('MGLTools OK')"],
                              capture_output=True, text=True, timeout=timeout)
        out = (proc.stdout or proc.stderr or "").strip()
        return (proc.returncode == 0, out or "MGLTools OK")
    except Exception as exc:
        return False, f"MGLTools error: {exc}"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vina_exe", type=Path, default=DEFAULT_VINA_EXE)
    p.add_argument("--mgl_python", type=Path, default=DEFAULT_MGL_PYTHON)
    p.add_argument("--timeout", type=int, default=15)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    checks = [
        ("Meeko PDBQT roundtrip", check_meeko_roundtrip),
        ("RDKit 3D embedding", check_rdkit_3d),
        ("BioPython PDBParser", check_biopython),
        ("Vina CLI", lambda: check_vina(args.vina_exe, args.timeout)),
        ("MGLTools Python", lambda: check_mgl(args.mgl_python, args.timeout)),
    ]

    results = []
    for name, fn in checks:
        ok, msg = fn()
        results.append((name, ok, msg))
        LOG.info("[%s] %s: %s", "PASS" if ok else "FAIL", name, msg)

    passed = sum(ok for _, ok, _ in results)
    print(f"\n=== Environment Validation: {passed}/{len(results)} passed ===")
    for name, ok, msg in results:
        print(f"  {'PASS' if ok else 'FAIL'} | {name} | {msg}")

    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
