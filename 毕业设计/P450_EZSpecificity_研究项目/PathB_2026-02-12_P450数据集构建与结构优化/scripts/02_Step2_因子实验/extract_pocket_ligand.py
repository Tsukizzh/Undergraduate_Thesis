"""
PathB Step 2: Extract pocket and ligand from PDB/CIF files.

Extends PathA's extraction with two factorial controls:
  --include_heme  : include heme-family residues in pocket output
  --pocket_radius : configurable neighbor-search radius (default 10 A)

Design: heme atoms are always added to pocket when enabled,
independent of radius filter, for clean factorial interpretation.
"""

import argparse
import csv
import logging
from collections import Counter
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch, PDBIO, PDBParser, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from rdkit import Chem

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_POCKET_RADIUS = 10.0
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
COFACTOR_RESNAMES = {"HEM", "FAD", "NAD", "NAP", "FMN"}
HEME_RESNAMES = {"HEM", "HEC", "HEA", "HEO", "DHE", "HEB"}


@dataclass
class ExtractionResult:
    dock_index: int
    pdb_id: str
    ligand_ccd: str
    success: bool
    pocket_atoms: int = 0
    ligand_atoms: int = 0
    heme_atoms_added: int = 0
    error: str = ""


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def load_dock_mapping(csv_path: Path) -> List[Dict]:
    with open(csv_path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def find_structure_file(pdb_dir: Path, pdb_id: str) -> Optional[Path]:
    pdb_id_upper = pdb_id.upper()
    for ext in (".pdb", ".cif"):
        p = pdb_dir / f"{pdb_id_upper}{ext}"
        if p.exists():
            return p
    return None


def parse_structure(path: Path) -> Optional[Structure]:
    try:
        parser = MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
        return parser.get_structure("s", str(path))
    except Exception as e:
        logger.error("Failed to parse %s: %s", path, e)
        return None


# ---------------------------------------------------------------------------
# AltLoc handling
# ---------------------------------------------------------------------------
def resolve_altloc(atom: Atom) -> Atom:
    if not atom.is_disordered():
        return atom
    alts = atom.disordered_get_list()

    def _key(a: Atom) -> Tuple[int, float]:
        alt = (a.get_altloc() or " ").strip() or " "
        occ = a.get_occupancy()
        occ = -1.0 if occ is None else float(occ)
        return (0 if alt == " " else 1 if alt == "A" else 2, -occ)

    return sorted(alts, key=_key)[0]


def _normalize_altloc(alt: Optional[str]) -> str:
    alt = (alt or "").strip()
    return "" if alt in {"", ".", "?"} else alt


def _preferred_residue_altloc(residue: Residue) -> str:
    counts: Counter = Counter()
    for atom in residue.get_atoms():
        if atom.is_disordered():
            for a in atom.disordered_get_list():
                alt = _normalize_altloc(a.get_altloc())
                if alt:
                    counts[alt] += 1
        else:
            alt = _normalize_altloc(atom.get_altloc())
            if alt:
                counts[alt] += 1
    if not counts:
        return ""
    best_n = max(counts.values())
    best_alts = sorted([a for a, n in counts.items() if n == best_n],
                       key=lambda a: (0 if a == "A" else 1, a))
    return best_alts[0]


# ---------------------------------------------------------------------------
# Ligand selection
# ---------------------------------------------------------------------------
def find_ligand_residues(model: Model, ligand_ccd: str) -> List[Residue]:
    ligand_ccd = ligand_ccd.strip().upper()
    candidates = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            if resname == ligand_ccd and resname not in COFACTOR_RESNAMES:
                candidates.append(residue)
    return candidates


def select_best_ligand(candidates: List[Residue], model: Model) -> Optional[Residue]:
    if not candidates:
        return None
    if len(candidates) == 1:
        return candidates[0]

    # P450-specific: prefer ligand closest to HEM FE atom
    fe_coord = None
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip().upper() == "HEM":
                for atom in residue:
                    if atom.get_name().strip().upper() == "FE":
                        fe_coord = atom.get_coord()
                        break
                if fe_coord is not None:
                    break
        if fe_coord is not None:
            break

    if fe_coord is not None:
        def _dist(res: Residue) -> float:
            dists = [np.linalg.norm(resolve_altloc(a).get_coord() - fe_coord)
                     for a in res]
            return min(dists) if dists else float("inf")
        return min(candidates, key=_dist)

    return max(candidates, key=lambda r: sum(1 for _ in r.get_atoms()))


# ---------------------------------------------------------------------------
# Atom collection (core modification vs PathA)
# ---------------------------------------------------------------------------
def collect_protein_atoms(
    model: Model,
    ligand_residue: Residue,
    include_heme: bool = False,
) -> Tuple[List[Atom], List[Atom]]:
    """Collect standard protein atoms and optionally heme atoms.

    Returns (regular_atoms, heme_atoms). Heme atoms are collected
    separately so they can be added to pocket independent of radius.
    """
    regular: List[Atom] = []
    heme: List[Atom] = []
    lig_id = ligand_residue.get_id()
    lig_chain = ligand_residue.get_parent().get_id()

    for chain in model:
        for residue in chain:
            hetflag, _, _ = residue.get_id()
            resname = residue.get_resname().strip().upper()

            if resname in WATER_RESNAMES:
                continue
            if chain.get_id() == lig_chain and residue.get_id() == lig_id:
                continue

            is_standard = (hetflag == " ")
            is_heme = include_heme and (hetflag != " ") and (resname in HEME_RESNAMES)

            if not is_standard and not is_heme:
                continue

            target = heme if is_heme else regular
            for atom in residue:
                resolved = resolve_altloc(atom)
                if (getattr(resolved, "element", "") or "").strip().upper() == "H":
                    continue
                target.append(resolved)

    return regular, heme


# ---------------------------------------------------------------------------
# Pocket extraction
# ---------------------------------------------------------------------------
def extract_pocket_atoms(
    protein_atoms: List[Atom],
    ligand_coords: List[np.ndarray],
    radius: float = DEFAULT_POCKET_RADIUS,
) -> Set[Atom]:
    if not protein_atoms:
        return set()
    ns = NeighborSearch(protein_atoms)
    pocket: Set[Atom] = set()
    for coord in ligand_coords:
        pocket.update(ns.search(coord, radius, level="A"))
    return pocket


# ---------------------------------------------------------------------------
# PDB / SDF writers
# ---------------------------------------------------------------------------
class PocketSelect(Select):
    def __init__(self, pocket_full_ids: Set[tuple], target_model_id: int):
        self._ids = pocket_full_ids
        self._model_id = target_model_id

    def accept_model(self, model):
        return 1 if model.get_id() == self._model_id else 0

    def accept_atom(self, atom):
        return 1 if atom.get_full_id() in self._ids else 0


def write_pocket_pdb(model: Model, pocket_atoms: Set[Atom], path: Path) -> bool:
    try:
        ids = {a.get_full_id() for a in pocket_atoms}
        io = PDBIO()
        io.set_structure(model.get_parent())
        io.save(str(path), select=PocketSelect(ids, model.get_id()))
        return True
    except Exception as e:
        logger.error("Failed to write pocket PDB: %s", e)
        return False


def ligand_to_pdb_block(model: Model, ligand_residue: Residue) -> str:
    class LigandSelect(Select):
        def __init__(self, target: Residue, model_id: int):
            self.target_id = target.get_id()
            self.target_chain = target.get_parent().get_id()
            self._model_id = model_id
            self._best_alt = _preferred_residue_altloc(target)

        def accept_model(self, m):
            return 1 if m.get_id() == self._model_id else 0

        def accept_chain(self, chain):
            return 1 if chain.get_id() == self.target_chain else 0

        def accept_residue(self, residue):
            return 1 if residue.get_id() == self.target_id else 0

        def accept_atom(self, atom):
            alt = _normalize_altloc(atom.get_altloc())
            if alt == "" or self._best_alt == "":
                return 1
            return 1 if alt == self._best_alt else 0

    io = PDBIO()
    io.set_structure(model.get_parent())
    buf = StringIO()
    io.save(buf, select=LigandSelect(ligand_residue, model.get_id()))
    return buf.getvalue()


def pdb_block_to_mol(pdb_block: str) -> Optional[Chem.Mol]:
    try:
        return Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, flavor=1)
    except Exception as e:
        logger.error("RDKit failed to parse PDB block: %s", e)
        return None


def write_ligand_sdf(mol: Chem.Mol, path: Path) -> bool:
    try:
        block = Chem.MolToMolBlock(mol)
        with open(path, "w", encoding="utf-8") as f:
            f.write(block)
            f.write("$$$$\n")
        return True
    except Exception as e:
        logger.error("Failed to write ligand SDF: %s", e)
        return False


# ---------------------------------------------------------------------------
# Single-record extraction pipeline
# ---------------------------------------------------------------------------
def extract_single(
    dock_index: int,
    pdb_id: str,
    ligand_ccd: str,
    pdb_dir: Path,
    pocket_dir: Path,
    ligand_dir: Path,
    radius: float = DEFAULT_POCKET_RADIUS,
    include_heme: bool = False,
) -> ExtractionResult:
    result = ExtractionResult(dock_index=dock_index, pdb_id=pdb_id,
                              ligand_ccd=ligand_ccd, success=False)

    structure_path = find_structure_file(pdb_dir, pdb_id)
    if structure_path is None:
        result.error = f"Structure file not found for {pdb_id}"
        return result

    structure = parse_structure(structure_path)
    if structure is None:
        result.error = f"Failed to parse {structure_path}"
        return result

    try:
        model = next(structure.get_models())
    except StopIteration:
        result.error = "Structure has no models"
        return result

    candidates = find_ligand_residues(model, ligand_ccd)
    if not candidates:
        result.error = f"Ligand {ligand_ccd} not found in {pdb_id}"
        return result

    ligand_residue = select_best_ligand(candidates, model)
    if ligand_residue is None:
        result.error = f"Failed to select ligand from {len(candidates)} candidates"
        return result

    ligand_coords = [resolve_altloc(a).get_coord() for a in ligand_residue]
    if not ligand_coords:
        result.error = "Ligand has no atoms"
        return result

    # Collect atoms: regular protein + optional heme
    regular_atoms, heme_atoms = collect_protein_atoms(model, ligand_residue, include_heme)
    result.heme_atoms_added = len(heme_atoms)

    if not regular_atoms and not heme_atoms:
        result.error = "No protein/heme atoms found"
        return result

    # Radius-based pocket from regular protein atoms only
    pocket_atoms = extract_pocket_atoms(regular_atoms, ligand_coords, radius)

    # Heme atoms always added when enabled (independent of radius)
    if include_heme and heme_atoms:
        pocket_atoms.update(heme_atoms)

    if not pocket_atoms:
        result.error = "No pocket atoms within radius"
        return result

    # Write outputs
    if not write_pocket_pdb(model, pocket_atoms, pocket_dir / f"{dock_index}.pdb"):
        result.error = "Failed to write pocket PDB"
        return result

    pdb_block = ligand_to_pdb_block(model, ligand_residue)
    mol = pdb_block_to_mol(pdb_block)
    if mol is None:
        result.error = "RDKit failed to parse ligand"
        return result

    if not write_ligand_sdf(mol, ligand_dir / f"{dock_index}.sdf"):
        result.error = "Failed to write ligand SDF"
        return result

    result.success = True
    result.pocket_atoms = len(pocket_atoms)
    result.ligand_atoms = mol.GetNumAtoms()
    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Extract pocket/ligand for PathB factorial experiments.")
    p.add_argument("--mapping_csv", required=True, help="Path to dock_index_mapping.csv")
    p.add_argument("--pdb_dir", required=True, help="Directory with source PDB/CIF files")
    p.add_argument("--output_dir", required=True, help="Output base dir (creates pocket/ and raw_ligand/ subdirs)")
    p.add_argument("--pocket_radius", type=float, default=DEFAULT_POCKET_RADIUS,
                    help=f"Pocket radius in angstroms (default: {DEFAULT_POCKET_RADIUS})")
    p.add_argument("--include_heme", action="store_true",
                    help="Include heme-family residues (HEM/HEC/HEA/HEO/DHE/HEB) in pocket")
    p.add_argument("--experiment_name", default="", help="Experiment name for logging/summary")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    mapping_csv = Path(args.mapping_csv)
    pdb_dir = Path(args.pdb_dir)
    output_base = Path(args.output_dir)
    pocket_dir = output_base / "pocket"
    ligand_dir = output_base / "raw_ligand"

    if not mapping_csv.exists():
        logger.error("Mapping CSV not found: %s", mapping_csv)
        return 2
    if not pdb_dir.exists():
        logger.error("PDB directory not found: %s", pdb_dir)
        return 2
    if args.pocket_radius <= 0:
        logger.error("--pocket_radius must be > 0, got %s", args.pocket_radius)
        return 2

    pocket_dir.mkdir(parents=True, exist_ok=True)
    ligand_dir.mkdir(parents=True, exist_ok=True)

    exp_name = args.experiment_name.strip() or "unnamed"
    logger.info("=" * 60)
    logger.info("PathB extract_pocket_ligand  experiment=%s", exp_name)
    logger.info("  pocket_radius=%.1f  include_heme=%s", args.pocket_radius, args.include_heme)
    logger.info("  mapping_csv=%s", mapping_csv)
    logger.info("  pdb_dir=%s", pdb_dir)
    logger.info("  output_dir=%s", output_base)
    logger.info("=" * 60)

    mapping = load_dock_mapping(mapping_csv)
    logger.info("Loaded %d records", len(mapping))

    results: List[ExtractionResult] = []
    ok = 0

    for i, row in enumerate(mapping):
        di = int(row["Dock_Index"])
        pid = row["pdb_id"]
        lig = row["ligand_ccd"]

        if (i + 1) % 50 == 0:
            logger.info("Processing %d/%d: %s/%s", i + 1, len(mapping), pid, lig)

        try:
            r = extract_single(di, pid, lig, pdb_dir, pocket_dir, ligand_dir,
                               radius=args.pocket_radius, include_heme=args.include_heme)
        except Exception as e:
            r = ExtractionResult(dock_index=di, pdb_id=pid, ligand_ccd=lig,
                                 success=False, error=f"Unexpected: {e}")
            logger.error("Dock_Index %d: %s", di, e)

        results.append(r)
        if r.success:
            ok += 1
        else:
            logger.warning("Dock_Index %d (%s/%s): %s", di, pid, lig, r.error)

    logger.info("Extraction complete: %d/%d successful", ok, len(mapping))

    # Failures report
    failures = [r for r in results if not r.success]
    if failures:
        fp = output_base / "extraction_failures.csv"
        with open(fp, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["Dock_Index", "pdb_id", "ligand_ccd", "error"])
            for r in failures:
                w.writerow([r.dock_index, r.pdb_id, r.ligand_ccd, r.error])
        logger.info("Failures: %s", fp)

    # Full summary
    sp = output_base / "extraction_summary.csv"
    with open(sp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Dock_Index", "pdb_id", "ligand_ccd", "success",
                     "pocket_atoms", "ligand_atoms", "heme_atoms_added",
                     "include_heme", "pocket_radius", "experiment_name"])
        for r in results:
            w.writerow([r.dock_index, r.pdb_id, r.ligand_ccd, r.success,
                         r.pocket_atoms, r.ligand_atoms, r.heme_atoms_added,
                         args.include_heme, args.pocket_radius, exp_name])
    logger.info("Summary: %s", sp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
