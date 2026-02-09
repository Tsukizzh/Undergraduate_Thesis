"""
Step 8.1: Extract pocket and ligand from PDB/CIF files.

This script extracts:
1. Pocket: protein atoms within 10A of ligand (ATOM records only, for PDBProtein compatibility)
2. Ligand: raw SDF file with experimental coordinates

Input:
- dock_index_mapping.csv: mapping from Dock_Index to (pdb_id, ligand_ccd)
- PDB/CIF source files

Output:
- str_tmp_data/pocket/{Dock_Index}.pdb
- str_tmp_data/raw_ligand/{Dock_Index}.sdf
"""

import os
import csv
import logging
from collections import Counter
from pathlib import Path
from io import StringIO
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, NeighborSearch
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from rdkit import Chem

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
POCKET_RADIUS = 10.0  # Angstroms
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O"}
COFACTOR_RESNAMES = {"HEM", "FAD", "NAD", "NAP", "FMN"}  # Common cofactors to exclude from ligand


@dataclass
class ExtractionResult:
    """Result of pocket/ligand extraction for one Dock_Index."""
    dock_index: int
    pdb_id: str
    ligand_ccd: str
    success: bool
    pocket_atoms: int = 0
    ligand_atoms: int = 0
    error: str = ""


def load_dock_mapping(csv_path: Path) -> List[Dict]:
    """Load dock_index_mapping.csv."""
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        return list(reader)


def find_structure_file(pdb_dir: Path, pdb_id: str) -> Optional[Path]:
    """Find PDB or CIF file for given PDB ID."""
    pdb_id_upper = pdb_id.upper()

    # Try PDB first
    pdb_file = pdb_dir / f"{pdb_id_upper}.pdb"
    if pdb_file.exists():
        return pdb_file

    # Try CIF
    cif_file = pdb_dir / f"{pdb_id_upper}.cif"
    if cif_file.exists():
        return cif_file

    return None


def parse_structure(structure_path: Path) -> Optional[Structure]:
    """Parse PDB or CIF file."""
    try:
        if structure_path.suffix.lower() == '.cif':
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        return parser.get_structure("structure", str(structure_path))
    except Exception as e:
        logger.error(f"Failed to parse {structure_path}: {e}")
        return None


def resolve_altloc(atom: Atom) -> Atom:
    """Resolve disordered atom to best altLoc."""
    if not atom.is_disordered():
        return atom

    # Get all altlocs
    alts = atom.disordered_get_list()

    def sort_key(a):
        alt = (a.get_altloc() or " ").strip() or " "
        occ = a.get_occupancy()
        occ = -1.0 if occ is None else float(occ)
        return (
            0 if alt == " " else 1 if alt == "A" else 2,
            -occ,
        )

    return sorted(alts, key=sort_key)[0]


def _normalize_altloc(alt: Optional[str]) -> str:
    """Normalize altLoc across PDB/mmCIF conventions.

    mmCIF commonly uses '.' or '?' for missing values.
    Treat those as blank altLoc.
    """
    alt = (alt or "").strip()
    if alt in {"", ".", "?"}:
        return ""
    return alt


def _preferred_residue_altloc(residue: Residue) -> str:
    """Pick a preferred altLoc for a residue (most common; tie-break prefer 'A')."""
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
    best_alts = [a for a, n in counts.items() if n == best_n]
    best_alts.sort(key=lambda a: (0 if a == "A" else 1, a))
    return best_alts[0]


def find_ligand_residues(model: Model, ligand_ccd: str) -> List[Residue]:
    """Find all residues matching the ligand CCD code."""
    ligand_ccd = ligand_ccd.strip().upper()
    candidates = []

    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            if resname == ligand_ccd:
                # Skip if it's actually a cofactor we want to exclude
                if resname in COFACTOR_RESNAMES:
                    continue
                candidates.append(residue)

    return candidates


def select_best_ligand(candidates: List[Residue], model: Model) -> Optional[Residue]:
    """Select the best ligand instance if multiple exist.

    Strategy:
    1. If only one candidate, return it
    2. For P450, prefer candidate closest to HEM iron (FE atom)
    3. Fallback: select candidate with most atoms
    """
    if not candidates:
        return None

    if len(candidates) == 1:
        return candidates[0]

    # Try to find HEM FE atom for P450-specific selection
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
        # Select candidate closest to FE
        def min_dist_to_fe(residue):
            min_dist = float('inf')
            for atom in residue:
                dist = np.linalg.norm(resolve_altloc(atom).get_coord() - fe_coord)
                min_dist = min(min_dist, dist)
            return min_dist

        return min(candidates, key=min_dist_to_fe)

    # Fallback: select by atom count
    return max(candidates, key=lambda r: sum(1 for _ in r.get_atoms()))


def get_ligand_coords(residue: Residue) -> List[np.ndarray]:
    """Get coordinates of all atoms in ligand residue."""
    coords = []
    for atom in residue:
        resolved = resolve_altloc(atom)
        coords.append(resolved.get_coord())
    return coords


def collect_protein_atoms(model: Model, ligand_residue: Residue) -> List[Atom]:
    """Collect all protein ATOM records (excluding water, ligand, and HETATM)."""
    protein_atoms = []
    ligand_id = ligand_residue.get_id()
    ligand_chain = ligand_residue.get_parent().get_id()

    for chain in model:
        for residue in chain:
            # Get residue identifier
            hetflag, resseq, icode = residue.get_id()
            resname = residue.get_resname().strip().upper()

            # Skip water
            if resname in WATER_RESNAMES:
                continue

            # Skip the target ligand itself
            if chain.get_id() == ligand_chain and residue.get_id() == ligand_id:
                continue

            # Only include standard amino acids (hetflag == ' ')
            # This ensures only ATOM records, not HETATM
            if hetflag != ' ':
                continue

            for atom in residue:
                resolved = resolve_altloc(atom)
                # Skip hydrogens
                if resolved.element == 'H':
                    continue
                protein_atoms.append(resolved)

    return protein_atoms


def extract_pocket_atoms(protein_atoms: List[Atom], ligand_coords: List[np.ndarray],
                         radius: float = POCKET_RADIUS) -> Set[Atom]:
    """Find protein atoms within radius of any ligand atom."""
    if not protein_atoms:
        return set()

    ns = NeighborSearch(protein_atoms)
    pocket_atoms = set()

    for coord in ligand_coords:
        nearby = ns.search(coord, radius, level="A")
        pocket_atoms.update(nearby)

    return pocket_atoms


class PocketSelect(Select):
    """PDBIO selector for pocket atoms using full_id for robust matching."""

    def __init__(self, pocket_full_ids: Set[tuple], target_model_id: int):
        self._pocket_ids = pocket_full_ids
        self._model_id = target_model_id

    def accept_model(self, model):
        return 1 if model.get_id() == self._model_id else 0

    def accept_atom(self, atom):
        return 1 if atom.get_full_id() in self._pocket_ids else 0


def write_pocket_pdb(model: Model, pocket_atoms: Set[Atom], output_path: Path) -> bool:
    """Write pocket atoms to PDB file."""
    try:
        pocket_ids = {a.get_full_id() for a in pocket_atoms}
        io = PDBIO()
        io.set_structure(model.get_parent())
        io.save(str(output_path), select=PocketSelect(pocket_ids, model.get_id()))
        return True
    except Exception as e:
        logger.error(f"Failed to write pocket PDB: {e}")
        return False


def ligand_to_pdb_block(model: Model, ligand_residue: Residue) -> str:
    """Extract ligand residue as PDB block with consistent altLoc selection.

    Fixes:
    - altLoc selection: previously only ""/'A' were accepted, dropping B/C/D-only ligands
    - mmCIF compatibility: normalize '.'/'?' altLoc markers to avoid dropping all atoms
    """
    class LigandSelect(Select):
        def __init__(self, target_residue, model_id):
            self.target_id = target_residue.get_id()
            self.target_chain = target_residue.get_parent().get_id()
            self._model_id = model_id
            self._best_alt = _preferred_residue_altloc(target_residue)

        def accept_model(self, m):
            return 1 if m.get_id() == self._model_id else 0

        def accept_chain(self, chain):
            return 1 if chain.get_id() == self.target_chain else 0

        def accept_residue(self, residue):
            return 1 if residue.get_id() == self.target_id else 0

        def accept_atom(self, atom):
            alt = _normalize_altloc(atom.get_altloc())
            if alt == "":
                return 1
            if self._best_alt == "":
                return 1
            return 1 if alt == self._best_alt else 0

    io = PDBIO()
    io.set_structure(model.get_parent())
    buf = StringIO()
    io.save(buf, select=LigandSelect(ligand_residue, model.get_id()))
    return buf.getvalue()


def pdb_block_to_mol(pdb_block: str) -> Optional[Chem.Mol]:
    """Convert PDB block to RDKit Mol.

    v3.1 fix: Use flavor=1 to read element symbols from columns 77-78.
    This fixes parsing issues when altLoc (column 17) is non-blank,
    causing RDKit to misread the residue name (e.g., "B5LW" instead of "5LW").
    """
    try:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, flavor=1)
        return mol
    except Exception as e:
        logger.error(f"RDKit failed to parse PDB block: {e}")
        return None


def write_ligand_sdf(mol: Chem.Mol, output_path: Path) -> bool:
    """Write ligand to SDF file.

    Uses MolToMolBlock instead of SDWriter to avoid RDKit's
    inability to handle non-ASCII characters in file paths.
    """
    try:
        mol_block = Chem.MolToMolBlock(mol)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(mol_block)
            f.write("$$$$\n")
        return True
    except Exception as e:
        logger.error(f"Failed to write ligand SDF: {e}")
        return False


def extract_single(
    dock_index: int,
    pdb_id: str,
    ligand_ccd: str,
    pdb_dir: Path,
    pocket_dir: Path,
    ligand_dir: Path
) -> ExtractionResult:
    """Extract pocket and ligand for a single Dock_Index."""
    result = ExtractionResult(
        dock_index=dock_index,
        pdb_id=pdb_id,
        ligand_ccd=ligand_ccd,
        success=False
    )

    # Find structure file
    structure_path = find_structure_file(pdb_dir, pdb_id)
    if structure_path is None:
        result.error = f"Structure file not found for {pdb_id}"
        return result

    # Parse structure
    structure = parse_structure(structure_path)
    if structure is None:
        result.error = f"Failed to parse {structure_path}"
        return result

    # Get first model (use iterator for safety - model id may not be 0)
    try:
        model = next(structure.get_models())
    except StopIteration:
        result.error = "Structure has no models"
        return result

    # Find ligand residues
    candidates = find_ligand_residues(model, ligand_ccd)
    if not candidates:
        result.error = f"Ligand {ligand_ccd} not found in {pdb_id}"
        return result

    # Select best ligand
    ligand_residue = select_best_ligand(candidates, model)
    if ligand_residue is None:
        result.error = f"Failed to select ligand from {len(candidates)} candidates"
        return result

    # Get ligand coordinates
    ligand_coords = get_ligand_coords(ligand_residue)
    if not ligand_coords:
        result.error = "Ligand has no atoms"
        return result

    # Collect protein atoms
    protein_atoms = collect_protein_atoms(model, ligand_residue)
    if not protein_atoms:
        result.error = "No protein atoms found"
        return result

    # Extract pocket
    pocket_atoms = extract_pocket_atoms(protein_atoms, ligand_coords)
    if not pocket_atoms:
        result.error = "No pocket atoms within radius"
        return result

    # Write pocket PDB
    pocket_path = pocket_dir / f"{dock_index}.pdb"
    if not write_pocket_pdb(model, pocket_atoms, pocket_path):
        result.error = "Failed to write pocket PDB"
        return result

    # Extract ligand as PDB block
    pdb_block = ligand_to_pdb_block(model, ligand_residue)

    # Convert to RDKit mol
    mol = pdb_block_to_mol(pdb_block)
    if mol is None:
        result.error = "RDKit failed to parse ligand"
        return result

    # Write ligand SDF
    ligand_path = ligand_dir / f"{dock_index}.sdf"
    if not write_ligand_sdf(mol, ligand_path):
        result.error = "Failed to write ligand SDF"
        return result

    # Success
    result.success = True
    result.pocket_atoms = len(pocket_atoms)
    result.ligand_atoms = mol.GetNumAtoms()

    return result


def main():
    """Main entry point."""
    # Paths (relative to this script)
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    data_dir = project_dir / "data"

    # Input paths
    mapping_csv = data_dir / "04_Step4_格式修正后数据" / "dock_index_mapping.csv"
    pdb_dir = data_dir / "01_Step1_PDB文件"

    # Output paths
    output_base = data_dir / "08_Step8_结构特征" / "str_tmp_data"
    pocket_dir = output_base / "pocket"
    raw_ligand_dir = output_base / "raw_ligand"

    # Ensure output directories exist
    pocket_dir.mkdir(parents=True, exist_ok=True)
    raw_ligand_dir.mkdir(parents=True, exist_ok=True)

    # Load mapping
    logger.info(f"Loading mapping from {mapping_csv}")
    mapping = load_dock_mapping(mapping_csv)
    logger.info(f"Loaded {len(mapping)} records")

    # Process each record
    results = []
    success_count = 0

    for i, row in enumerate(mapping):
        dock_index = int(row['Dock_Index'])
        pdb_id = row['pdb_id']
        ligand_ccd = row['ligand_ccd']

        if (i + 1) % 50 == 0:
            logger.info(f"Processing {i + 1}/{len(mapping)}: {pdb_id}/{ligand_ccd}")

        try:
            result = extract_single(
                dock_index=dock_index,
                pdb_id=pdb_id,
                ligand_ccd=ligand_ccd,
                pdb_dir=pdb_dir,
                pocket_dir=pocket_dir,
                ligand_dir=raw_ligand_dir
            )
        except Exception as e:
            result = ExtractionResult(
                dock_index=dock_index, pdb_id=pdb_id,
                ligand_ccd=ligand_ccd, success=False,
                error=f"Unexpected: {e}"
            )
            logger.error(f"Dock_Index {dock_index}: unexpected error: {e}")

        results.append(result)

        if result.success:
            success_count += 1
        else:
            logger.warning(f"Dock_Index {dock_index} ({pdb_id}/{ligand_ccd}): {result.error}")

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"Extraction complete: {success_count}/{len(mapping)} successful")

    # Write failures report
    failures = [r for r in results if not r.success]
    if failures:
        failures_path = output_base / "extraction_failures.csv"
        with open(failures_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Dock_Index', 'pdb_id', 'ligand_ccd', 'error'])
            for r in failures:
                writer.writerow([r.dock_index, r.pdb_id, r.ligand_ccd, r.error])
        logger.info(f"Failures written to {failures_path}")

    # Write success summary
    summary_path = output_base / "extraction_summary.csv"
    with open(summary_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Dock_Index', 'pdb_id', 'ligand_ccd', 'success', 'pocket_atoms', 'ligand_atoms'])
        for r in results:
            writer.writerow([r.dock_index, r.pdb_id, r.ligand_ccd, r.success, r.pocket_atoms, r.ligand_atoms])
    logger.info(f"Summary written to {summary_path}")


if __name__ == "__main__":
    main()
