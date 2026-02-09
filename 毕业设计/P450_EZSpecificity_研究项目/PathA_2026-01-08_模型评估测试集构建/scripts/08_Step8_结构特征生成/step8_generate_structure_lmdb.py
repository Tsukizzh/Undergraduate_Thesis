"""
Step 8.2: Generate Structure Features LMDB.

This script:
1. Reads pocket PDB files using PDBProtein
2. Reads aligned ligand SDF files using parse_sdf_file_mol (with Chinese path workaround)
3. Creates StructureComplexData objects (pickle-compatible with Datasets.Structure.utils)
4. Stores results in LMDB keyed by Dock_Index

REQUIREMENTS:
- Full EZSpecificity environment with torch_geometric

Input:
- str_tmp_data/pocket/{Dock_Index}.pdb
- str_tmp_data/ligand/{Dock_Index}.sdf (aligned, with AtomMapNum set)
- alignment_summary.csv (for success list)

Output:
- structure_features.lmdb
- high_quality_id.txt (Dock_Index list where both pocket and ligand exist)
"""

import os
import sys
import csv
import pickle
import logging
from pathlib import Path

import lmdb
import torch
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, HybridizationType
from torch_geometric.data import Data


class StructureComplexData(Data):
    """StructureComplexData compatible with Datasets.Structure.utils for pickle.

    This class is defined locally but with __module__ set to 'Datasets.Structure.utils'
    so that pickle stores the correct module path, ensuring compatibility with
    downstream training code that unpickles using the original module.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @staticmethod
    def from_protein_ligand_dicts(protein_dict=None, ligand_dict=None, **kwargs):
        instance = StructureComplexData(**kwargs)
        if protein_dict is not None:
            for key, item in protein_dict.items():
                instance['protein_' + key] = item
        if ligand_dict is not None:
            for key, item in ligand_dict.items():
                instance['ligand_' + key] = item
        return instance

    def __inc__(self, key, value, *args, **kwargs):
        if key == 'complex_edge_index':
            return self['mask_ligand'].size(0)
        else:
            return super().__inc__(key, value, *args, **kwargs)


# Set __module__ for pickle compatibility with downstream training code
StructureComplexData.__module__ = 'Datasets.Structure.utils'

# Create fake module hierarchy so pickle can serialize without importing real Datasets package
import types
if 'Datasets' not in sys.modules:
    sys.modules['Datasets'] = types.ModuleType('Datasets')
if 'Datasets.Structure' not in sys.modules:
    sys.modules['Datasets.Structure'] = types.ModuleType('Datasets.Structure')
if 'Datasets.Structure.utils' not in sys.modules:
    fake_utils = types.ModuleType('Datasets.Structure.utils')
    sys.modules['Datasets.Structure.utils'] = fake_utils
sys.modules['Datasets.Structure.utils'].StructureComplexData = StructureComplexData


def torchify_dict(data):
    """Convert numpy arrays to torch tensors."""
    output = {}
    for k, v in data.items():
        if isinstance(v, np.ndarray):
            output[k] = torch.from_numpy(v)
        else:
            output[k] = v
    return output


# ========== Standalone implementations (no torch_scatter dependency) ==========

class PDBProtein:
    """Simplified PDBProtein parser - only to_dict_atom() needed."""

    AA_NAME_SYM = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
        'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
        'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y', 'UNK': 'X'
    }
    AA_NAME_NUMBER = {k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())}
    BACKBONE_NAMES = ["CA", "C", "N", "O"]

    def __init__(self, data, mode='auto'):
        self.fn = data
        if (data[-4:].lower() == '.pdb' and mode == 'auto') or mode == 'path':
            with open(data, 'r') as f:
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
        self._parse()

    def _parse(self):
        for line in self.block.splitlines():
            if line[0:6].strip() == 'ATOM':
                element_symb = line[76:78].strip().capitalize()
                if len(element_symb) == 0:
                    element_symb = line[13:14]
                atom = {
                    'atom_name': line[12:16].strip(),
                    'res_name': line[17:20].strip(),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element_symb': element_symb,
                }
                self.atoms.append(atom)
                atomic_number = self.ptable.GetAtomicNumber(atom['element_symb'])
                self.element.append(atomic_number)
                self.atomic_weight.append(self.ptable.GetAtomicWeight(atomic_number))
                self.pos.append(np.array([atom['x'], atom['y'], atom['z']], dtype=np.float32))
                self.atom_name.append(atom['atom_name'])
                self.is_backbone.append(atom['atom_name'] in self.BACKBONE_NAMES)
                if atom['res_name'] not in self.AA_NAME_NUMBER:
                    self.atom_to_aa_type.append(self.AA_NAME_NUMBER['UNK'])
                else:
                    self.atom_to_aa_type.append(self.AA_NAME_NUMBER[atom['res_name']])
            elif line[0:6].strip() == 'HEADER':
                self.title = line[10:].strip().lower()
            elif line[0:6].strip() == 'ENDMDL':
                break

    def to_dict_atom(self):
        return {
            'element': np.array(self.element, dtype=int),
            'molecule_name': self.title,
            'pos': np.array(self.pos, dtype=np.float32) if self.pos else np.zeros((0, 3), dtype=np.float32),
            'is_backbone': np.array(self.is_backbone, dtype=bool),
            'atom_name': self.atom_name,
            'atom_to_aa_type': np.array(self.atom_to_aa_type, dtype=int)
        }


def get_ligand_atom_features(rdmol):
    """Compute atom features without torch_scatter."""
    num_atoms = rdmol.GetNumAtoms()
    atomic_number = []
    aromatic = []
    hybrid = []
    degree = []

    HYBRID_TYPES = {t: i for i, t in enumerate(HybridizationType.names.values())}

    for atom_idx in range(num_atoms):
        atom = rdmol.GetAtomWithIdx(atom_idx)
        atomic_number.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        hybrid.append(HYBRID_TYPES[atom.GetHybridization()])
        degree.append(atom.GetDegree())

    # Build edge list
    row, col = [], []
    for bond in rdmol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]

    # Count hydrogen neighbors (replacement for torch_scatter)
    node_type = np.array(atomic_number)
    hs = (node_type == 1).astype(float)
    num_hs = np.zeros(num_atoms, dtype=float)
    for src, dst in zip(row, col):
        num_hs[dst] += hs[src]

    feat_mat = np.array([atomic_number, aromatic, degree, num_hs.astype(int), hybrid], dtype=int).transpose()
    return feat_mat


def parse_sdf_file_mol(path, mol=None, heavy_only=True):
    """Parse SDF file or use pre-loaded mol."""
    if mol is None:
        mol = next(iter(Chem.SDMolSupplier(path, removeHs=heavy_only, sanitize=False)))

    feat_mat = get_ligand_atom_features(mol)
    ptable = Chem.GetPeriodicTable()

    num_atoms = mol.GetNumAtoms()
    pos = mol.GetConformer().GetPositions()

    element = []
    indexs = []
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
        BondType.AROMATIC: 4, BondType.UNSPECIFIED: 5
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
        'element': element,
        'pos': pos,
        'bond_index': edge_index,
        'bond_type': edge_type,
        'center_of_mass': center_of_mass,
        'atom_feature': feat_mat,
        'index': indexs
    }


# ========== End standalone implementations ==========

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_sdf_as_mol(sdf_path: Path, heavy_only: bool = True):
    """Load SDF file using Python open() to avoid RDKit's non-ASCII path issue.

    Args:
        sdf_path: Path to SDF file
        heavy_only: If True, remove hydrogens

    Returns:
        RDKit Mol object or None if parsing fails
    """
    try:
        with open(sdf_path, 'r') as f:
            content = f.read()
        # Split on $$$$ and take first mol block (do NOT strip!)
        mol_block = content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)

        if mol is None:
            return None

        # Enforce heavy_only since we bypassed SDMolSupplier
        if heavy_only:
            mol = Chem.RemoveHs(mol)

        return mol
    except Exception as e:
        logger.error(f"Failed to load {sdf_path}: {e}")
        return None


def process_single(
    dock_index: int,
    pocket_path: Path,
    ligand_path: Path
) -> tuple:
    """Process a single pocket-ligand pair into StructureComplexData.

    Returns:
        (success: bool, data: StructureComplexData or None, error: str)
    """
    # Check files exist
    if not pocket_path.exists():
        return False, None, f"Pocket file not found: {pocket_path}"

    if not ligand_path.exists():
        return False, None, f"Ligand file not found: {ligand_path}"

    # Parse pocket PDB
    try:
        pocket_dict = PDBProtein(str(pocket_path)).to_dict_atom()
        if len(pocket_dict['element']) == 0:
            return False, None, "Pocket has no atoms"
    except Exception as e:
        return False, None, f"PDBProtein parse error: {e}"

    # Load ligand mol
    mol = load_sdf_as_mol(ligand_path, heavy_only=True)
    if mol is None:
        return False, None, "RDKit failed to load ligand SDF"

    # Check conformer exists
    if mol.GetNumConformers() == 0:
        return False, None, "Ligand has no conformer (no 3D coords)"

    # Parse ligand using parse_sdf_file_mol with pre-loaded mol
    try:
        ligand_dict = parse_sdf_file_mol(None, mol=mol, heavy_only=True)
    except Exception as e:
        return False, None, f"parse_sdf_file_mol error: {e}"

    # Create StructureComplexData
    try:
        data = StructureComplexData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        # Verify protein_pos exists and has atoms
        if data.protein_pos.size(0) == 0:
            return False, None, "protein_pos is empty"
    except Exception as e:
        return False, None, f"StructureComplexData creation error: {e}"

    return True, data, ""


def main():
    """Main entry point."""
    # Paths
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    data_dir = project_dir / "data"

    # Input paths
    str_tmp_dir = data_dir / "08_Step8_结构特征" / "str_tmp_data"
    pocket_dir = str_tmp_dir / "pocket"
    ligand_dir = str_tmp_dir / "ligand"
    summary_csv = str_tmp_dir / "alignment_summary.csv"

    # Output paths
    output_dir = data_dir / "08_Step8_结构特征"
    lmdb_path = output_dir / "structure_features.lmdb"
    high_quality_path = output_dir / "high_quality_id.txt"

    # Ensure output dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load alignment summary to get success list
    logger.info(f"Loading alignment summary from {summary_csv}")
    with open(summary_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        records = [r for r in reader if r['success'] == 'True']
    logger.info(f"Found {len(records)} successfully aligned records")

    # Open LMDB
    logger.info(f"Creating LMDB at {lmdb_path}")
    db = lmdb.open(
        str(lmdb_path),
        map_size=10 * (1024 * 1024 * 1024),  # 10GB
        create=True,
        subdir=False,
        readonly=False,
    )

    # Process each record
    success_count = 0
    high_quality_ids = []
    failures = []

    for i, row in enumerate(records):
        dock_index = int(row['Dock_Index'])
        pocket_path = pocket_dir / f"{dock_index}.pdb"
        ligand_path = ligand_dir / f"{dock_index}.sdf"

        if (i + 1) % 50 == 0:
            logger.info(f"Processing {i + 1}/{len(records)}")

        success, data, error = process_single(dock_index, pocket_path, ligand_path)

        if success:
            # Write to LMDB
            with db.begin(write=True, buffers=True) as txn:
                txn.put(
                    key=str(dock_index).encode(),
                    value=pickle.dumps(data)
                )
            high_quality_ids.append(dock_index)
            success_count += 1
        else:
            failures.append((dock_index, error))
            logger.warning(f"Dock_Index {dock_index}: {error}")

    db.close()

    # Write high quality IDs
    logger.info(f"Writing {len(high_quality_ids)} high quality IDs to {high_quality_path}")
    with open(high_quality_path, 'w') as f:
        for idx in sorted(high_quality_ids):
            f.write(f"{idx}\n")

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"LMDB generation complete: {success_count}/{len(records)} successful")
    logger.info(f"LMDB path: {lmdb_path}")
    logger.info(f"High quality IDs: {high_quality_path}")

    if failures:
        failures_path = str_tmp_dir / "lmdb_failures.csv"
        with open(failures_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow(['Dock_Index', 'error'])
            for idx, err in failures:
                writer.writerow([idx, err])
        logger.info(f"Failures written to {failures_path}")


if __name__ == "__main__":
    main()
