"""
PathB Step 8.3: Generate Structure Features LMDB.

Fork of PathA's step8_generate_structure_lmdb.py with key modifications:
  1. PDBProtein supports both ATOM and HETATM records (for Heme factor experiment)
  2. All paths configurable via CLI arguments
  3. HETATM statistics logging for gate verification
  4. Optional dataset filtering (only process Dock Indexes present in B6 data.csv)

Input:
  - pocket/{Dock_Index}.pdb  (from Step 8.1, may contain HETATM if --include_heme)
  - ligand/{Dock_Index}.sdf  (from Step 8.2, aligned with AtomMapNum)
  - alignment_summary.csv    (from Step 8.2, for success filtering)
  - [optional] data.csv      (B6 dataset, for Dock Index filtering)

Output:
  - structure_features.lmdb
  - high_quality_id.txt
  - structure_build_summary.csv (with HETATM statistics)
"""

import argparse
import csv
import logging
import pickle
import sys
import types
from pathlib import Path

import lmdb
import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, HybridizationType
from torch_geometric.data import Data

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


# ============================================================================
# StructureComplexData (pickle-compatible with downstream training code)
# ============================================================================
class StructureComplexData(Data):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_protein_ligand_dicts(protein_dict=None, ligand_dict=None, **kwargs):
        instance = StructureComplexData(**kwargs)
        if protein_dict is not None:
            for key, item in protein_dict.items():
                instance["protein_" + key] = item
        if ligand_dict is not None:
            for key, item in ligand_dict.items():
                instance["ligand_" + key] = item
        return instance

    def __inc__(self, key, value, *args, **kwargs):
        if key == "complex_edge_index":
            return self["mask_ligand"].size(0)
        return super().__inc__(key, value, *args, **kwargs)


StructureComplexData.__module__ = "Datasets.Structure.utils"

if "Datasets" not in sys.modules:
    sys.modules["Datasets"] = types.ModuleType("Datasets")
if "Datasets.Structure" not in sys.modules:
    sys.modules["Datasets.Structure"] = types.ModuleType("Datasets.Structure")
if "Datasets.Structure.utils" not in sys.modules:
    fake_utils = types.ModuleType("Datasets.Structure.utils")
    sys.modules["Datasets.Structure.utils"] = fake_utils
sys.modules["Datasets.Structure.utils"].StructureComplexData = StructureComplexData


def torchify_dict(data):
    output = {}
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            output[k] = torch.from_numpy(v)
        else:
            output[k] = v
    return output


# ============================================================================
# PDBProtein — supports ATOM + HETATM (core PathB modification)
# ============================================================================
class PDBProtein:
    """PDB pocket parser that reads both ATOM and HETATM records.

    Differences from PathA version:
      - Parses HETATM lines (same column layout as ATOM per wwPDB standard)
      - Tracks record_type per atom for statistics
      - Residue grouping key includes record_type to prevent collisions
    """

    AA_NAME_SYM = {
        "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
        "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
        "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
        "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
        "UNK": "X",
    }
    AA_NAME_NUMBER = {k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())}
    BACKBONE_NAMES = ["CA", "C", "N", "O"]

    HEME_RESNAMES = {"HEM", "HEC", "HEA", "HEO", "DHE", "HEB"}

    def __init__(self, data, mode="auto"):
        self.fn = data
        if (data[-4:].lower() == ".pdb" and mode == "auto") or mode == "path":
            with open(data, "r") as f:
                self.block = f.read()
        else:
            self.block = data

        self.ptable = Chem.GetPeriodicTable()
        self.title = None
        self.atoms = []
        self.element = []
        self.atomic_weight = []
        self.pos = []
        self.atom_name = []
        self.is_backbone = []
        self.atom_to_aa_type = []

        # HETATM statistics
        self.n_atom = 0
        self.n_hetatm = 0
        self.n_heme_atoms = 0
        self.n_fe_atoms = 0

        self._parse()

    def _parse(self):
        for line in self.block.splitlines():
            record_type = line[0:6].strip()

            if record_type in ("ATOM", "HETATM"):
                element_symb = line[76:78].strip().capitalize()
                if len(element_symb) == 0:
                    # Fallback: atom name cols 13-14 (handles 2-letter elements like FE)
                    element_symb = line[12:14].strip().capitalize()

                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()

                atom = {
                    "record_type": record_type,
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element_symb": element_symb,
                }
                self.atoms.append(atom)

                atomic_number = self.ptable.GetAtomicNumber(atom["element_symb"])
                self.element.append(atomic_number)
                self.atomic_weight.append(self.ptable.GetAtomicWeight(atomic_number))
                self.pos.append(np.array([atom["x"], atom["y"], atom["z"]], dtype=np.float32))
                self.atom_name.append(atom_name)
                self.is_backbone.append(atom_name in self.BACKBONE_NAMES and record_type == "ATOM")

                if res_name not in self.AA_NAME_NUMBER:
                    self.atom_to_aa_type.append(self.AA_NAME_NUMBER["UNK"])
                else:
                    self.atom_to_aa_type.append(self.AA_NAME_NUMBER[res_name])

                # Statistics
                if record_type == "ATOM":
                    self.n_atom += 1
                else:
                    self.n_hetatm += 1
                    if res_name.upper() in self.HEME_RESNAMES:
                        self.n_heme_atoms += 1
                    if element_symb.upper() == "FE":
                        self.n_fe_atoms += 1

            elif record_type == "HEADER":
                self.title = line[10:].strip().lower()
            elif record_type == "ENDMDL":
                break

    def to_dict_atom(self):
        return {
            "element": np.array(self.element, dtype=int),
            "molecule_name": self.title,
            "pos": np.array(self.pos, dtype=np.float32) if self.pos else np.zeros((0, 3), dtype=np.float32),
            "is_backbone": np.array(self.is_backbone, dtype=bool),
            "atom_name": self.atom_name,
            "atom_to_aa_type": np.array(self.atom_to_aa_type, dtype=int),
        }


# ============================================================================
# Ligand feature extraction (standalone, no torch_scatter)
# ============================================================================
def get_ligand_atom_features(rdmol):
    num_atoms = rdmol.GetNumAtoms()
    atomic_number, aromatic, hybrid, degree = [], [], [], []
    HYBRID_TYPES = {t: i for i, t in enumerate(HybridizationType.names.values())}

    for atom_idx in range(num_atoms):
        atom = rdmol.GetAtomWithIdx(atom_idx)
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybrid.append(HYBRID_TYPES[atom.GetHybridization()])
        degree.append(atom.GetDegree())

    row, col = [], []
    for bond in rdmol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]

    node_type = np.array(atomic_number)
    hs = (node_type == 1).astype(float)
    num_hs = np.zeros(num_atoms, dtype=float)
    for src, dst in zip(row, col):
        num_hs[dst] += hs[src]

    feat_mat = np.array(
        [atomic_number, aromatic, degree, num_hs.astype(int), hybrid], dtype=int
    ).transpose()
    return feat_mat


def parse_sdf_file_mol(path, mol=None, heavy_only=True):
    if mol is None:
        mol = next(iter(Chem.SDMolSupplier(path, removeHs=heavy_only, sanitize=False)))
    feat_mat = get_ligand_atom_features(mol)
    ptable = Chem.GetPeriodicTable()

    num_atoms = mol.GetNumAtoms()
    pos = mol.GetConformer().GetPositions()
    element, indexs = [], []
    accum_pos = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    accum_mass = 0.0

    for atom_idx in range(num_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        atomic_number = atom.GetAtomicNum()
        element.append(atomic_number)
        indexs.append(atom.GetAtomMapNum())
        x, y, z = pos[atom_idx]
        atomic_weight = ptable.GetAtomicWeight(atomic_number)
        accum_pos += np.array([x, y, z]) * atomic_weight
        accum_mass += atomic_weight

    center_of_mass = np.array(accum_pos / accum_mass, dtype=np.float32)
    element = np.array(element, dtype=int)
    pos = np.array(pos, dtype=np.float32)
    indexs = np.array(indexs, dtype=int)

    row, col, edge_type = [], [], []
    BOND_TYPES = {
        BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3,
        BondType.AROMATIC: 4, BondType.UNSPECIFIED: 5,
    }
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [BOND_TYPES[bond.GetBondType()]]

    edge_index = np.array([row, col], dtype=int)
    edge_type = np.array(edge_type, dtype=int)
    perm = (edge_index[0] * num_atoms + edge_index[1]).argsort()
    edge_index = edge_index[:, perm]
    edge_type = edge_type[perm]

    return {
        "element": element,
        "pos": pos,
        "bond_index": edge_index,
        "bond_type": edge_type,
        "center_of_mass": center_of_mass,
        "atom_feature": feat_mat,
        "index": indexs,
    }


# ============================================================================
# SDF loader (Chinese path compatible)
# ============================================================================
def load_sdf_as_mol(sdf_path: Path, heavy_only: bool = True):
    try:
        with open(sdf_path, "r") as f:
            content = f.read()
        mol_block = content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        if mol is None:
            return None
        if heavy_only:
            mol = Chem.RemoveHs(mol)
        return mol
    except Exception as e:
        logger.error("Failed to load %s: %s", sdf_path, e)
        return None


# ============================================================================
# Single record processing
# ============================================================================
def process_single(dock_index: int, pocket_path: Path, ligand_path: Path):
    """Process one pocket-ligand pair.

    Returns:
        (success, data_or_None, error_str, stats_dict)
    """
    stats = {"n_atom": 0, "n_hetatm": 0, "n_heme_atoms": 0, "n_fe_atoms": 0, "protein_atoms": 0}

    if not pocket_path.exists():
        return False, None, f"Pocket not found: {pocket_path.name}", stats
    if not ligand_path.exists():
        return False, None, f"Ligand not found: {ligand_path.name}", stats

    # Parse pocket PDB (with HETATM support)
    try:
        pdb = PDBProtein(str(pocket_path))
        pocket_dict = pdb.to_dict_atom()
        stats["n_atom"] = pdb.n_atom
        stats["n_hetatm"] = pdb.n_hetatm
        stats["n_heme_atoms"] = pdb.n_heme_atoms
        stats["n_fe_atoms"] = pdb.n_fe_atoms
        stats["protein_atoms"] = len(pocket_dict["element"])
        if len(pocket_dict["element"]) == 0:
            return False, None, "Pocket has no atoms", stats
    except Exception as e:
        return False, None, f"PDBProtein parse error: {e}", stats

    # Load ligand
    mol = load_sdf_as_mol(ligand_path, heavy_only=True)
    if mol is None:
        return False, None, "RDKit failed to load ligand SDF", stats
    if mol.GetNumConformers() == 0:
        return False, None, "Ligand has no conformer", stats

    # Verify AtomMapNum integrity (skip for single-atom ligands: 0-based map is valid)
    map_nums = [mol.GetAtomWithIdx(i).GetAtomMapNum() for i in range(mol.GetNumAtoms())]
    if mol.GetNumAtoms() > 1 and all(n == 0 for n in map_nums):
        return False, None, "Ligand AtomMapNum all zero (alignment may have failed)", stats

    # Parse ligand features
    try:
        ligand_dict = parse_sdf_file_mol(None, mol=mol, heavy_only=True)
    except Exception as e:
        return False, None, f"parse_sdf_file_mol error: {e}", stats

    # Build StructureComplexData
    try:
        data = StructureComplexData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        if data.protein_pos.size(0) == 0:
            return False, None, "protein_pos is empty", stats
    except Exception as e:
        return False, None, f"StructureComplexData error: {e}", stats

    return True, data, "", stats


# ============================================================================
# CLI & main
# ============================================================================
def parse_args():
    p = argparse.ArgumentParser(description="PathB Step 8.3: Generate structure LMDB with HETATM support")
    p.add_argument("--pocket_dir", required=True, help="Directory with pocket PDB files")
    p.add_argument("--ligand_dir", required=True, help="Directory with aligned ligand SDF files")
    p.add_argument("--alignment_summary", required=True, help="alignment_summary.csv from Step 8.2")
    p.add_argument("--output_dir", required=True, help="Output directory for LMDB and high_quality_id.txt")
    p.add_argument("--dataset_csv", default="", help="Optional: B6 data.csv to filter Dock Indexes")
    p.add_argument("--experiment_name", default="", help="Experiment name for logging")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    pocket_dir = Path(args.pocket_dir)
    ligand_dir = Path(args.ligand_dir)
    alignment_csv = Path(args.alignment_summary)
    output_dir = Path(args.output_dir)
    exp_name = args.experiment_name.strip() or "unnamed"

    output_dir.mkdir(parents=True, exist_ok=True)
    lmdb_path = output_dir / "structure_features.lmdb"
    hq_path = output_dir / "high_quality_id.txt"
    summary_path = output_dir / "structure_build_summary.csv"

    # --- Load alignment success list ---
    if not alignment_csv.exists():
        logger.error("alignment_summary.csv not found: %s", alignment_csv)
        return 2
    with open(alignment_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        aligned_docks = {int(r["Dock_Index"]) for r in reader if r["success"].strip().lower() == "true"}
    logger.info("Aligned Dock_Indexes: %d", len(aligned_docks))

    # --- Optional: filter by B6 dataset ---
    if args.dataset_csv:
        ds_path = Path(args.dataset_csv)
        if not ds_path.exists():
            logger.error("dataset_csv not found: %s", ds_path)
            return 2
        with open(ds_path, "r", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            b6_docks = {int(r["Dock Index"]) for r in reader}
        logger.info("B6 Dock_Indexes: %d", len(b6_docks))
        target_docks = sorted(aligned_docks & b6_docks)
        logger.info("Intersection (aligned ∩ B6): %d", len(target_docks))
    else:
        target_docks = sorted(aligned_docks)

    if not target_docks:
        logger.error("No target Dock_Indexes after filtering — cannot build LMDB")
        return 2

    logger.info("=" * 60)
    logger.info("PathB step8_generate_structure_lmdb  experiment=%s", exp_name)
    logger.info("  pocket_dir=%s", pocket_dir)
    logger.info("  ligand_dir=%s", ligand_dir)
    logger.info("  target records=%d", len(target_docks))
    logger.info("=" * 60)

    # --- Remove stale LMDB if exists ---
    if lmdb_path.exists():
        lmdb_path.unlink()
    lock_path = Path(str(lmdb_path) + "-lock")
    if lock_path.exists():
        lock_path.unlink()

    # --- Open LMDB ---
    db = lmdb.open(
        str(lmdb_path),
        map_size=10 * (1024 * 1024 * 1024),
        create=True,
        subdir=False,
        readonly=False,
    )

    # --- Process ---
    success_count = 0
    high_quality_ids = []
    summary_rows = []

    total_n_atom = 0
    total_n_hetatm = 0
    total_n_heme = 0
    total_n_fe = 0

    for i, dock_index in enumerate(target_docks):
        if (i + 1) % 50 == 0:
            logger.info("Processing %d/%d (Dock_Index=%d)", i + 1, len(target_docks), dock_index)

        pocket_path = pocket_dir / f"{dock_index}.pdb"
        ligand_path = ligand_dir / f"{dock_index}.sdf"

        success, data, error, stats = process_single(dock_index, pocket_path, ligand_path)

        summary_rows.append({
            "Dock_Index": dock_index,
            "success": success,
            "error": error,
            "protein_atoms": stats["protein_atoms"],
            "n_atom": stats["n_atom"],
            "n_hetatm": stats["n_hetatm"],
            "n_heme_atoms": stats["n_heme_atoms"],
            "n_fe_atoms": stats["n_fe_atoms"],
        })

        if success:
            with db.begin(write=True, buffers=True) as txn:
                txn.put(key=str(dock_index).encode(), value=pickle.dumps(data))
            high_quality_ids.append(dock_index)
            success_count += 1
            total_n_atom += stats["n_atom"]
            total_n_hetatm += stats["n_hetatm"]
            total_n_heme += stats["n_heme_atoms"]
            total_n_fe += stats["n_fe_atoms"]
        else:
            logger.warning("Dock_Index %d: %s", dock_index, error)

    db.close()

    # --- Write high quality IDs ---
    with open(hq_path, "w") as f:
        for idx in sorted(high_quality_ids):
            f.write(f"{idx}\n")

    # --- Write summary CSV ---
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Dock_Index", "success", "error", "protein_atoms",
            "n_atom", "n_hetatm", "n_heme_atoms", "n_fe_atoms",
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    # --- Final report ---
    logger.info("=" * 60)
    logger.info("LMDB generation complete: %d/%d successful", success_count, len(target_docks))
    logger.info("HETATM statistics (successful records):")
    logger.info("  total ATOM lines:   %d", total_n_atom)
    logger.info("  total HETATM lines: %d", total_n_hetatm)
    logger.info("  total Heme atoms:   %d", total_n_heme)
    logger.info("  total Fe atoms:     %d", total_n_fe)
    logger.info("Output: %s", lmdb_path)
    logger.info("High quality IDs: %s", hq_path)
    logger.info("Summary: %s", summary_path)
    logger.info("=" * 60)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
