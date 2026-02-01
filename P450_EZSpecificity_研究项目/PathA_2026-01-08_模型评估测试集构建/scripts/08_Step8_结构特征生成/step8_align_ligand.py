"""
Step 8.1b: Align PDB-extracted ligands to SMILES atom ordering.

This script:
1. Loads SMILES from Substrates.csv
2. Loads PDB-extracted ligand SDF from raw_ligand/
3. Uses RDKit AssignBondOrdersFromTemplate to fix bond perception
4. Uses GetSubstructMatches to map SMILES atoms to PDB atoms
5. Sets AtomMapNum for each atom (= SMILES atom index)
6. Saves aligned SDF to ligand/

Input:
- str_tmp_data/raw_ligand/{Dock_Index}.sdf
- dock_index_mapping.csv (for Dock_Index -> Substrate_Index mapping)
- Substrates.csv (for Substrate_Index -> SMILES)

Output:
- str_tmp_data/ligand/{Dock_Index}.sdf (with AtomMapNum set)
- alignment_summary.csv
"""

import csv
import logging
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    dock_index: int
    substrate_index: int
    success: bool
    num_atoms_smiles: int = 0
    num_atoms_pdb: int = 0
    num_matched: int = 0
    error: str = ""


def load_smiles_mapping(substrates_csv: Path) -> dict:
    """Load Substrate_Index -> SMILES mapping."""
    mapping = {}
    with open(substrates_csv, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx = int(row.get('Substrate_Index', -1))
            smiles = row.get('Substrate_SMILES', '')
            if idx >= 0 and smiles:
                mapping[idx] = smiles
    return mapping


def load_dock_mapping(dock_csv: Path) -> list:
    """Load dock_index_mapping.csv."""
    with open(dock_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        return list(reader)


def load_sdf_mol(sdf_path: Path) -> Optional[Chem.Mol]:
    """Load SDF file using Python open() for Chinese path compatibility."""
    try:
        with open(sdf_path, 'r') as f:
            content = f.read()
        mol_block = content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        return mol
    except Exception as e:
        logger.error(f"Failed to load {sdf_path}: {e}")
        return None


def write_sdf_mol(mol: Chem.Mol, output_path: Path) -> bool:
    """Write mol to SDF using Python open() for Chinese path compatibility."""
    try:
        mol_block = Chem.MolToMolBlock(mol)
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(mol_block)
            f.write("$$$$\n")
        return True
    except Exception as e:
        logger.error(f"Failed to write {output_path}: {e}")
        return False


def align_pdb_to_smiles(
    smiles: str,
    pdb_mol: Chem.Mol
) -> Tuple[Optional[Chem.Mol], str]:
    """Align PDB-extracted mol to SMILES atom ordering.

    Returns:
        (aligned_mol, error_message)
        aligned_mol is None on failure
    """
    # Create reference mol from SMILES
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return None, "Invalid SMILES"
    ref = Chem.RemoveHs(ref)

    # Remove hydrogens from PDB mol (may fail on invalid valence)
    try:
        pdb = Chem.RemoveHs(pdb_mol)
    except Exception as e:
        return None, f"RemoveHs failed: {e}"

    # Fix bond orders using SMILES as template
    try:
        pdb_fixed = AllChem.AssignBondOrdersFromTemplate(ref, pdb)
    except Exception as e:
        return None, f"AssignBondOrdersFromTemplate failed: {e}"

    if pdb_fixed is None:
        return None, "AssignBondOrdersFromTemplate returned None"

    # Find substructure match (ref -> pdb_fixed)
    # This returns tuple of pdb atom indices for each ref atom
    matches = pdb_fixed.GetSubstructMatches(ref, useChirality=True)
    if not matches:
        matches = pdb_fixed.GetSubstructMatches(ref, useChirality=False)

    if not matches:
        return None, "No substructure match found"

    # Use first match (deterministic)
    match = matches[0]

    # Verify full match
    if len(match) != ref.GetNumAtoms():
        return None, f"Partial match: {len(match)}/{ref.GetNumAtoms()} atoms"

    # Set AtomMapNum on pdb_fixed atoms
    # match[smiles_idx] = pdb_idx
    # We want pdb atom to carry the smiles index as its map number
    for smiles_idx, pdb_idx in enumerate(match):
        pdb_fixed.GetAtomWithIdx(pdb_idx).SetAtomMapNum(smiles_idx)

    return pdb_fixed, ""


def process_single(
    dock_index: int,
    substrate_index: int,
    smiles: str,
    raw_ligand_dir: Path,
    aligned_ligand_dir: Path
) -> AlignmentResult:
    """Process a single ligand alignment."""
    result = AlignmentResult(
        dock_index=dock_index,
        substrate_index=substrate_index,
        success=False
    )

    # Load raw ligand
    raw_path = raw_ligand_dir / f"{dock_index}.sdf"
    if not raw_path.exists():
        result.error = "Raw ligand file not found"
        return result

    pdb_mol = load_sdf_mol(raw_path)
    if pdb_mol is None:
        result.error = "Failed to load raw ligand SDF"
        return result

    result.num_atoms_pdb = pdb_mol.GetNumAtoms()

    # Get SMILES atom count
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        result.error = "Invalid SMILES"
        return result
    ref = Chem.RemoveHs(ref)
    result.num_atoms_smiles = ref.GetNumAtoms()

    # Align
    aligned_mol, error = align_pdb_to_smiles(smiles, pdb_mol)
    if aligned_mol is None:
        result.error = error
        return result

    result.num_matched = aligned_mol.GetNumAtoms()

    # Write aligned ligand
    output_path = aligned_ligand_dir / f"{dock_index}.sdf"
    if not write_sdf_mol(aligned_mol, output_path):
        result.error = "Failed to write aligned SDF"
        return result

    result.success = True
    return result


def main():
    """Main entry point."""
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    data_dir = project_dir / "data"

    # Input paths
    dock_mapping_csv = data_dir / "04_Step4_格式修正后数据" / "dock_index_mapping.csv"
    substrates_csv = data_dir / "04_Step4_格式修正后数据" / "Substrates.csv"
    raw_ligand_dir = data_dir / "08_Step8_结构特征" / "str_tmp_data" / "raw_ligand"

    # Output paths
    aligned_ligand_dir = data_dir / "08_Step8_结构特征" / "str_tmp_data" / "ligand"
    summary_csv = data_dir / "08_Step8_结构特征" / "str_tmp_data" / "alignment_summary.csv"

    # Ensure output dir exists
    aligned_ligand_dir.mkdir(parents=True, exist_ok=True)

    # Load mappings
    logger.info(f"Loading SMILES from {substrates_csv}")
    smiles_map = load_smiles_mapping(substrates_csv)
    logger.info(f"Loaded {len(smiles_map)} SMILES")

    logger.info(f"Loading dock mapping from {dock_mapping_csv}")
    dock_records = load_dock_mapping(dock_mapping_csv)
    logger.info(f"Loaded {len(dock_records)} dock records")

    # Process each dock record
    results = []
    success_count = 0

    for i, row in enumerate(dock_records):
        dock_index = int(row['Dock_Index'])
        substrate_index = int(row['Substrate_Index'])

        if substrate_index not in smiles_map:
            result = AlignmentResult(
                dock_index=dock_index,
                substrate_index=substrate_index,
                success=False,
                error=f"Substrate_Index {substrate_index} not in SMILES map"
            )
            results.append(result)
            continue

        smiles = smiles_map[substrate_index]

        if (i + 1) % 50 == 0:
            logger.info(f"Processing {i + 1}/{len(dock_records)}")

        try:
            result = process_single(
                dock_index=dock_index,
                substrate_index=substrate_index,
                smiles=smiles,
                raw_ligand_dir=raw_ligand_dir,
                aligned_ligand_dir=aligned_ligand_dir
            )
        except Exception as e:
            result = AlignmentResult(
                dock_index=dock_index,
                substrate_index=substrate_index,
                success=False,
                error=f"Unexpected: {e}"
            )
            logger.error(f"Dock_Index {dock_index}: unexpected error: {e}")

        results.append(result)

        if result.success:
            success_count += 1
        else:
            logger.warning(f"Dock_Index {dock_index}: {result.error}")

    # Write summary
    logger.info(f"Writing summary to {summary_csv}")
    with open(summary_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([
            'Dock_Index', 'Substrate_Index', 'success',
            'num_atoms_smiles', 'num_atoms_pdb', 'num_matched', 'error'
        ])
        for r in results:
            writer.writerow([
                r.dock_index, r.substrate_index, r.success,
                r.num_atoms_smiles, r.num_atoms_pdb, r.num_matched, r.error
            ])

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"Alignment complete: {success_count}/{len(dock_records)} successful")
    logger.info(f"Aligned ligands: {aligned_ligand_dir}")
    logger.info(f"Summary: {summary_csv}")


if __name__ == "__main__":
    main()
