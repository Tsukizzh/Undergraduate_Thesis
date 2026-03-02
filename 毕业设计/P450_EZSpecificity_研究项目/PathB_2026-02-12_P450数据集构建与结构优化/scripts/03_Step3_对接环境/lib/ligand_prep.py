"""Prepare ligand SMILES → 3D conformer → PDBQT using RDKit + Meeko."""
from __future__ import annotations

import logging
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

ASCII_TMP_ROOT = Path(r"D:\autodock\tmp")


@dataclass
class LigandPrepResult:
    pdbqt_path: Optional[Path]
    sdf_path: Optional[Path]
    success: bool
    error: str = ""
    atom_count: int = 0


def _embed_3d(mol: Chem.Mol, seed: int) -> bool:
    """Attempt 3D embedding with ETKDGv3 → ETKDG → random coords fallback."""
    params_v3 = AllChem.ETKDGv3()
    params_v3.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params_v3) == 0:
        return True

    params_v2 = AllChem.ETKDG()
    params_v2.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params_v2) == 0:
        return True

    return AllChem.EmbedMolecule(
        mol, randomSeed=seed, useRandomCoords=True, maxAttempts=200,
    ) == 0


def _needs_ascii_workaround(path: Path) -> bool:
    """Check if path contains non-ASCII characters (RDKit SDWriter fails on these)."""
    try:
        str(path).encode("ascii")
        return False
    except UnicodeEncodeError:
        return True


def _minimize(mol: Chem.Mol) -> None:
    """Energy-minimize with MMFF94, fallback to UFF."""
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            if AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94") == 0:
                return
    except Exception:
        pass
    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception:
        logger.debug("Both MMFF94 and UFF failed; using embedded geometry directly")


def prepare_ligand_from_smiles(
    smiles: str,
    ligand_id: str,
    pdbqt_dir: Path,
    sdf_dir: Path,
    random_seed: int = 2026,
) -> LigandPrepResult:
    """Convert SMILES → 3D SDF + PDBQT.

    Args:
        smiles: Input SMILES string.
        ligand_id: Identifier used for output filenames (e.g. "substrate_42").
        pdbqt_dir: Directory for .pdbqt output.
        sdf_dir: Directory for .sdf output (3D reference).
        random_seed: Seed for conformer generation reproducibility.
    """
    safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", ligand_id).strip("._") or "ligand"
    pdbqt_path = pdbqt_dir / f"{safe_id}.pdbqt"
    sdf_path = sdf_dir / f"{safe_id}.sdf"
    result = LigandPrepResult(pdbqt_path=pdbqt_path, sdf_path=sdf_path, success=False)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result.error = f"Invalid SMILES: {smiles}"
        return result

    mol = Chem.AddHs(mol)
    if not _embed_3d(mol, seed=random_seed):
        result.error = "3D embedding failed (ETKDGv3/ETKDG/random all failed)"
        return result

    _minimize(mol)

    pdbqt_dir.mkdir(parents=True, exist_ok=True)
    sdf_dir.mkdir(parents=True, exist_ok=True)

    use_tmp = _needs_ascii_workaround(sdf_path) or _needs_ascii_workaround(pdbqt_path)
    if use_tmp:
        ASCII_TMP_ROOT.mkdir(parents=True, exist_ok=True)
        tmp_dir = Path(tempfile.mkdtemp(dir=ASCII_TMP_ROOT, prefix="lig_"))
    else:
        tmp_dir = None

    try:
        write_sdf = tmp_dir / f"{safe_id}.sdf" if tmp_dir else sdf_path
        writer = Chem.SDWriter(str(write_sdf))
        if writer is None:
            result.error = f"Failed to create SDF writer: {write_sdf}"
            return result
        writer.write(mol)
        writer.close()

        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            result.error = "Meeko returned no setup objects"
            return result

        pdbqt_text, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            result.error = f"Meeko PDBQT conversion failed: {err}"
            return result

        write_pdbqt = tmp_dir / f"{safe_id}.pdbqt" if tmp_dir else pdbqt_path
        write_pdbqt.write_text(pdbqt_text, encoding="utf-8")

        if tmp_dir:
            shutil.copy2(str(write_sdf), str(sdf_path))
            shutil.copy2(str(write_pdbqt), str(pdbqt_path))

    except Exception as exc:
        result.error = f"Ligand conversion error: {exc}"
        return result
    finally:
        if tmp_dir and tmp_dir.exists():
            shutil.rmtree(str(tmp_dir), ignore_errors=True)

    result.atom_count = mol.GetNumAtoms()
    result.success = True
    return result
