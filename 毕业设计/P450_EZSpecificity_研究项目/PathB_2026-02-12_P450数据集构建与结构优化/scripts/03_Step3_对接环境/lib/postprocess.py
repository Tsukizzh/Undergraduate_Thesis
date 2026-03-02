"""Post-process Vina docked PDBQT → pocket PDB + ligand SDF for Step 8-10 pipeline."""
from __future__ import annotations

import logging
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Set

import numpy as np
from Bio.PDB import NeighborSearch, PDBIO, PDBParser, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem

logger = logging.getLogger(__name__)
ASCII_TMP_ROOT = Path(r"D:\autodock\tmp")

WATER_RESNAMES = frozenset({"HOH", "WAT", "DOD", "H2O"})
HEME_RESNAMES = frozenset({"HEM", "HEC", "HEA", "HEO", "DHE", "HEB"})
DEFAULT_POCKET_RADIUS = 10.0
FE_DISTANCE_MIN = 2.0
FE_DISTANCE_MAX = 15.0


@dataclass
class PostprocessResult:
    dock_index: int
    pocket_path: Optional[Path]
    ligand_path: Optional[Path]
    success: bool
    error: str = ""
    pocket_atoms: int = 0
    ligand_atoms: int = 0
    fe_ligand_distance: Optional[float] = None


class _PocketSelect(Select):
    """Select only atoms belonging to the extracted pocket."""

    def __init__(self, pocket_full_ids: Set[tuple], target_model_id: int):
        self._ids = pocket_full_ids
        self._model_id = target_model_id

    def accept_model(self, model: Model) -> int:
        return 1 if model.get_id() == self._model_id else 0

    def accept_atom(self, atom: Atom) -> int:
        return 1 if atom.get_full_id() in self._ids else 0


def _resolve_altloc(atom: Atom) -> Atom:
    """Pick the best alt-loc conformer (prefer A or highest occupancy)."""
    if not atom.is_disordered():
        return atom
    alts = atom.disordered_get_list()

    def _key(a: Atom):
        alt = (a.get_altloc() or " ").strip() or " "
        occ = a.get_occupancy()
        occ = -1.0 if occ is None else float(occ)
        return (0 if alt == " " else 1 if alt == "A" else 2, -occ)

    return sorted(alts, key=_key)[0]


def _collect_protein_atoms(model: Model) -> list[Atom]:
    """Collect standard protein heavy atoms with altloc resolution.

    Excludes water, HEM (Gate A: noHeme), all other HETATM, and hydrogen.
    """
    atoms: list[Atom] = []
    for residue in model.get_residues():
        resname = residue.get_resname().strip().upper()
        hetflag = residue.get_id()[0].strip()
        if resname in WATER_RESNAMES:
            continue
        if resname in HEME_RESNAMES:
            continue
        if hetflag:
            continue
        for atom in residue.get_atoms():
            resolved = _resolve_altloc(atom)
            element = (getattr(resolved, "element", "") or "").strip().upper()
            if element == "H":
                continue
            atoms.append(resolved)
    return atoms


def _find_heme_fe(model: Model) -> Optional[np.ndarray]:
    """Return FE atom coordinates from HEM residue, or None."""
    for residue in model.get_residues():
        if residue.get_resname().strip().upper() not in HEME_RESNAMES:
            continue
        for atom in residue.get_atoms():
            if atom.get_name().strip().upper() == "FE":
                return np.asarray(atom.get_coord(), dtype=float)
    return None


def _extract_pocket(
    protein_atoms: list[Atom],
    ligand_coords: np.ndarray,
    radius: float,
) -> Set[Atom]:
    """Find protein atoms within radius of any ligand atom."""
    if not protein_atoms:
        return set()
    ns = NeighborSearch(protein_atoms)
    pocket: Set[Atom] = set()
    for coord in ligand_coords:
        pocket.update(ns.search(coord, radius, level="A"))
    return pocket


def _ligand_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    coords = np.zeros((mol.GetNumAtoms(), 3), dtype=float)
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        coords[i] = [p.x, p.y, p.z]
    return coords


def _needs_ascii_workaround(path: Path) -> bool:
    try:
        str(path).encode("ascii")
        return False
    except UnicodeEncodeError:
        return True


def _cleanup_on_failure(pocket_path: Path, ligand_path: Path) -> None:
    """Remove partial output files to prevent downstream confusion."""
    for p in (pocket_path, ligand_path):
        try:
            if p.exists():
                p.unlink()
        except OSError:
            pass


def postprocess_docked_pose(
    docked_pdbqt: Path,
    original_pdb: Path,
    dock_index: int,
    pocket_dir: Path,
    ligand_dir: Path,
    pocket_radius: float = DEFAULT_POCKET_RADIUS,
    expected_heavy_atoms: Optional[int] = None,
    fe_dist_min: float = FE_DISTANCE_MIN,
    fe_dist_max: float = FE_DISTANCE_MAX,
) -> PostprocessResult:
    """Convert docked PDBQT → pocket PDB + ligand SDF.

    Steps:
    1. Meeko reads docked PDBQT → RDKit mol (best pose)
    2. QC: ligand atom count, FE-ligand distance
    3. Extract pocket from original PDB (protein heavy atoms within radius)
    4. Save ligand SDF + pocket PDB only after all QC passes
    """
    pocket_path = pocket_dir / f"{dock_index}.pdb"
    ligand_path = ligand_dir / f"{dock_index}.sdf"
    result = PostprocessResult(
        dock_index=dock_index,
        pocket_path=pocket_path,
        ligand_path=ligand_path,
        success=False,
    )

    if not docked_pdbqt.exists():
        result.error = f"Docked PDBQT not found: {docked_pdbqt}"
        return result
    if not original_pdb.exists():
        result.error = f"Original PDB not found: {original_pdb}"
        return result

    try:
        pocket_dir.mkdir(parents=True, exist_ok=True)
        ligand_dir.mkdir(parents=True, exist_ok=True)

        # 1. Meeko: docked PDBQT → RDKit mol
        pdbqt_mol = PDBQTMolecule.from_file(str(docked_pdbqt))
        poses = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        if not poses or poses[0] is None:
            result.error = "Meeko failed to recover ligand from docked PDBQT"
            return result

        best_pose = poses[0]
        result.ligand_atoms = best_pose.GetNumAtoms()

        # 2a. QC: ligand heavy atom count (optional)
        if expected_heavy_atoms is not None:
            actual_heavy = sum(1 for a in best_pose.GetAtoms() if a.GetAtomicNum() > 1)
            if actual_heavy != expected_heavy_atoms:
                result.error = (
                    f"Heavy atom mismatch: expected {expected_heavy_atoms}, "
                    f"got {actual_heavy}"
                )
                return result

        # 3. Parse original PDB and extract pocket
        coords = _ligand_coords(best_pose)
        structure = PDBParser(QUIET=True).get_structure("receptor", str(original_pdb))
        model = next(structure.get_models(), None)
        if model is None:
            result.error = "No model in receptor structure"
            return result

        # 2b. QC: FE-ligand distance
        fe_coord = _find_heme_fe(model)
        if fe_coord is None:
            result.error = "HEM FE not found in original PDB"
            return result

        dists = np.linalg.norm(coords - fe_coord[None, :], axis=1)
        result.fe_ligand_distance = float(np.min(dists))
        if not (fe_dist_min <= result.fe_ligand_distance <= fe_dist_max):
            result.error = (
                f"FE-ligand distance {result.fe_ligand_distance:.2f}Å "
                f"outside [{fe_dist_min}, {fe_dist_max}]Å"
            )
            return result

        protein_atoms = _collect_protein_atoms(model)
        if not protein_atoms:
            result.error = "No protein heavy atoms after filtering"
            return result

        pocket_atoms = _extract_pocket(protein_atoms, coords, pocket_radius)
        result.pocket_atoms = len(pocket_atoms)
        if result.pocket_atoms == 0:
            result.error = "No pocket atoms within radius"
            return result

        # 4. All QC passed — write outputs (via ASCII temp dir if needed)
        use_tmp = _needs_ascii_workaround(pocket_path) or _needs_ascii_workaround(ligand_path)
        if use_tmp:
            ASCII_TMP_ROOT.mkdir(parents=True, exist_ok=True)
            tmp_dir = Path(tempfile.mkdtemp(dir=ASCII_TMP_ROOT, prefix="post_"))
            tmp_lig = tmp_dir / ligand_path.name
            tmp_pkt = tmp_dir / pocket_path.name
        else:
            tmp_dir = None
            tmp_lig = ligand_path
            tmp_pkt = pocket_path

        try:
            writer = Chem.SDWriter(str(tmp_lig))
            writer.write(best_pose)
            writer.close()

            ids = {a.get_full_id() for a in pocket_atoms}
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(tmp_pkt), select=_PocketSelect(ids, model.get_id()))

            if tmp_dir:
                shutil.copy2(str(tmp_lig), str(ligand_path))
                shutil.copy2(str(tmp_pkt), str(pocket_path))
        finally:
            if tmp_dir and tmp_dir.exists():
                shutil.rmtree(str(tmp_dir), ignore_errors=True)

        result.success = True
        return result

    except Exception as exc:
        result.error = f"Postprocess error: {exc}"
        _cleanup_on_failure(pocket_path, ligand_path)
        return result
