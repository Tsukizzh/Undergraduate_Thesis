"""Prepare protein receptor PDB → PDBQT using MGLTools (with HEM preserved)."""
from __future__ import annotations

import logging
import re
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)

WATER_RESNAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})
HEME_RESNAMES = frozenset({"HEM", "HEC", "HEA", "HEO", "DHE", "HEB"})
ION_RESNAMES = frozenset({
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO",
    "NI", "CD", "HG", "AL", "CS", "SR", "BA", "PB", "LI", "RB",
    "IOD", "BR", "F", "SO4", "PO4", "NO3",
})

DEFAULT_TMP_DIR = Path(r"D:\autodock\tmp")
DEFAULT_MGL_PYTHON = Path(r"D:\autodock\MGLTools-1.5.7\python.exe")
DEFAULT_PREPARE_SCRIPT = Path(
    r"D:\autodock\MGLTools-1.5.7\Lib\site-packages"
    r"\AutoDockTools\Utilities24\prepare_receptor4.py"
)


@dataclass
class ReceptorPrepResult:
    pdbqt_path: Optional[Path]
    success: bool
    error: str = ""
    atom_count: int = 0


def _safe_stem(name: str) -> str:
    """Sanitize filename to ASCII-safe characters."""
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._")
    return cleaned or "receptor"


class _CleanupSelect(Select):
    """Keep standard protein residues + HEM only.

    Removes water, ions, co-crystal ligands, and all other HETATM residues.
    Co-crystal ligands must be removed so Vina can dock into the active site.
    """

    def accept_residue(self, residue: Residue) -> int:
        resname = residue.get_resname().strip().upper()
        if resname in HEME_RESNAMES:
            return 1
        hetflag = residue.get_id()[0].strip()
        if hetflag:
            return 0
        return 1

    def accept_atom(self, atom) -> int:
        altloc = atom.get_altloc()
        if altloc and altloc.strip() and altloc.strip() != "A":
            return 0
        return 1


def _strip_altloc(pdb_path: Path) -> None:
    """Blank altloc indicator (column 17) before running MGLTools."""
    lines = pdb_path.read_text(encoding="utf-8").splitlines(keepends=True)
    with pdb_path.open("w", encoding="utf-8") as fh:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
                line = line[:16] + " " + line[17:]
            fh.write(line)


def _count_pdbqt_atoms(pdbqt_path: Path) -> int:
    count = 0
    with pdbqt_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                count += 1
    return count


def _has_fe_atom(pdbqt_path: Path) -> bool:
    with pdbqt_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom_name = line[12:16].strip().upper()
                if atom_name == "FE":
                    return True
    return False


def prepare_receptor_pdbqt(
    pdb_path: Path,
    output_pdbqt_path: Path,
    tmp_dir: Path = DEFAULT_TMP_DIR,
    mgl_python: Path = DEFAULT_MGL_PYTHON,
    prepare_script: Path = DEFAULT_PREPARE_SCRIPT,
    timeout_sec: int = 300,
) -> ReceptorPrepResult:
    """Convert raw PDB → cleaned PDBQT via MGLTools.

    Steps: parse PDB → remove water/ions → keep HEM → copy to ASCII tmp dir
    → run prepare_receptor4.py → verify FE presence → copy result back.
    """
    result = ReceptorPrepResult(pdbqt_path=output_pdbqt_path, success=False)

    if not pdb_path.exists():
        result.error = f"Input PDB not found: {pdb_path}"
        return result
    if not mgl_python.exists():
        result.error = f"MGLTools python.exe not found: {mgl_python}"
        return result
    if not prepare_script.exists():
        result.error = f"prepare_receptor4.py not found: {prepare_script}"
        return result

    try:
        structure = PDBParser(QUIET=True).get_structure("receptor", str(pdb_path))
    except Exception as exc:
        result.error = f"PDB parse failed: {exc}"
        return result

    token = uuid.uuid4().hex[:10]
    stem = _safe_stem(pdb_path.stem)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cleaned_path = tmp_dir / f"{stem}_{token}_clean.pdb"
    pdbqt_tmp_path = tmp_dir / f"{stem}_{token}.pdbqt"
    temp_files = [cleaned_path, pdbqt_tmp_path]

    try:
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(cleaned_path), select=_CleanupSelect())
        _strip_altloc(cleaned_path)

        cmd = [
            str(mgl_python), str(prepare_script),
            "-r", str(cleaned_path),
            "-o", str(pdbqt_tmp_path),
            "-A", "hydrogens",
            "-U", "nphs_lps_waters",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)

        if proc.returncode != 0:
            result.error = (proc.stderr or proc.stdout or "Unknown MGLTools error").strip()
            return result

        if not pdbqt_tmp_path.exists():
            result.error = "MGLTools finished but no PDBQT generated"
            return result

        if not _has_fe_atom(pdbqt_tmp_path):
            result.error = "Receptor PDBQT missing FE atom — HEM may have been removed"
            return result

        output_pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pdbqt_tmp_path, output_pdbqt_path)
        result.atom_count = _count_pdbqt_atoms(output_pdbqt_path)
        result.success = True
        return result

    except subprocess.TimeoutExpired:
        result.error = f"MGLTools timed out after {timeout_sec}s"
        return result
    except Exception as exc:
        result.error = f"Unexpected receptor prep error: {exc}"
        return result
    finally:
        for tmp_file in temp_files:
            try:
                if tmp_file.exists():
                    tmp_file.unlink()
            except OSError:
                pass
