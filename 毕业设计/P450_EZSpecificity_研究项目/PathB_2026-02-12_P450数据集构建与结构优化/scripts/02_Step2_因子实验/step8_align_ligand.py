"""
PathB Step 8.2: Align PDB-extracted ligands to SMILES atom ordering.

Fork of PathA's step8_align_ligand.py with configurable paths via CLI.
Core alignment logic unchanged — only I/O paths are parameterized.

Input:
  - raw_ligand/{Dock_Index}.sdf  (from Step 8.1)
  - dock_index_mapping.csv       (Dock_Index -> Substrate_Index mapping)
  - Substrates.csv               (Substrate_Index -> SMILES)

Output:
  - ligand/{Dock_Index}.sdf      (aligned, with AtomMapNum set)
  - alignment_summary.csv
"""

import argparse
import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
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
    mapping = {}
    with open(substrates_csv, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            idx = int(row.get("Substrate_Index", -1))
            smiles = row.get("Substrate_SMILES", "")
            if idx >= 0 and smiles:
                mapping[idx] = smiles
    return mapping


def load_dock_mapping(dock_csv: Path) -> list:
    with open(dock_csv, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def load_sdf_mol(sdf_path: Path) -> Optional[Chem.Mol]:
    try:
        with open(sdf_path, "r") as f:
            content = f.read()
        mol_block = content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        return mol
    except Exception as e:
        logger.error("Failed to load %s: %s", sdf_path, e)
        return None


def write_sdf_mol(mol: Chem.Mol, output_path: Path) -> bool:
    try:
        mol_block = Chem.MolToMolBlock(mol)
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(mol_block)
            f.write("$$$$\n")
        return True
    except Exception as e:
        logger.error("Failed to write %s: %s", output_path, e)
        return False


def align_pdb_to_smiles(smiles: str, pdb_mol: Chem.Mol) -> Tuple[Optional[Chem.Mol], str]:
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return None, "Invalid SMILES"
    ref = Chem.RemoveHs(ref)

    try:
        pdb = Chem.RemoveHs(pdb_mol)
    except Exception as e:
        return None, f"RemoveHs failed: {e}"

    try:
        pdb_fixed = AllChem.AssignBondOrdersFromTemplate(ref, pdb)
    except Exception as e:
        return None, f"AssignBondOrdersFromTemplate failed: {e}"

    if pdb_fixed is None:
        return None, "AssignBondOrdersFromTemplate returned None"

    matches = pdb_fixed.GetSubstructMatches(ref, useChirality=True)
    if not matches:
        matches = pdb_fixed.GetSubstructMatches(ref, useChirality=False)
    if not matches:
        return None, "No substructure match found"

    match = matches[0]
    if len(match) != ref.GetNumAtoms():
        return None, f"Partial match: {len(match)}/{ref.GetNumAtoms()} atoms"

    for smiles_idx, pdb_idx in enumerate(match):
        pdb_fixed.GetAtomWithIdx(pdb_idx).SetAtomMapNum(smiles_idx)

    return pdb_fixed, ""


def process_single(
    dock_index: int,
    substrate_index: int,
    smiles: str,
    raw_ligand_dir: Path,
    aligned_ligand_dir: Path,
) -> AlignmentResult:
    result = AlignmentResult(dock_index=dock_index, substrate_index=substrate_index, success=False)

    raw_path = raw_ligand_dir / f"{dock_index}.sdf"
    if not raw_path.exists():
        result.error = "Raw ligand file not found"
        return result

    pdb_mol = load_sdf_mol(raw_path)
    if pdb_mol is None:
        result.error = "Failed to load raw ligand SDF"
        return result
    result.num_atoms_pdb = pdb_mol.GetNumAtoms()

    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        result.error = "Invalid SMILES"
        return result
    ref = Chem.RemoveHs(ref)
    result.num_atoms_smiles = ref.GetNumAtoms()

    aligned_mol, error = align_pdb_to_smiles(smiles, pdb_mol)
    if aligned_mol is None:
        result.error = error
        return result
    result.num_matched = aligned_mol.GetNumAtoms()

    output_path = aligned_ligand_dir / f"{dock_index}.sdf"
    if not write_sdf_mol(aligned_mol, output_path):
        result.error = "Failed to write aligned SDF"
        return result

    result.success = True
    return result


def parse_args():
    p = argparse.ArgumentParser(description="PathB Step 8.2: Align ligands to SMILES ordering")
    p.add_argument("--raw_ligand_dir", required=True, help="Input: raw_ligand/ from Step 8.1")
    p.add_argument("--aligned_ligand_dir", required=True, help="Output: aligned ligand/ directory")
    p.add_argument("--mapping_csv", required=True, help="dock_index_mapping.csv")
    p.add_argument("--substrates_csv", required=True, help="Substrates.csv with SMILES")
    p.add_argument("--summary_csv", default="", help="Output: alignment_summary.csv path (default: aligned_ligand_dir/../alignment_summary.csv)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    raw_ligand_dir = Path(args.raw_ligand_dir)
    aligned_ligand_dir = Path(args.aligned_ligand_dir)
    mapping_csv = Path(args.mapping_csv)
    substrates_csv = Path(args.substrates_csv)
    summary_csv = Path(args.summary_csv) if args.summary_csv else aligned_ligand_dir.parent / "alignment_summary.csv"

    if not raw_ligand_dir.exists():
        logger.error("raw_ligand_dir not found: %s", raw_ligand_dir)
        return 2
    if not mapping_csv.exists():
        logger.error("mapping_csv not found: %s", mapping_csv)
        return 2
    if not substrates_csv.exists():
        logger.error("substrates_csv not found: %s", substrates_csv)
        return 2

    aligned_ligand_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading SMILES from %s", substrates_csv)
    smiles_map = load_smiles_mapping(substrates_csv)
    logger.info("Loaded %d SMILES", len(smiles_map))

    logger.info("Loading dock mapping from %s", mapping_csv)
    dock_records = load_dock_mapping(mapping_csv)
    logger.info("Loaded %d dock records", len(dock_records))

    results = []
    success_count = 0

    for i, row in enumerate(dock_records):
        dock_index = int(row["Dock_Index"])
        substrate_index = int(row["Substrate_Index"])

        if substrate_index not in smiles_map:
            result = AlignmentResult(
                dock_index=dock_index, substrate_index=substrate_index,
                success=False, error=f"Substrate_Index {substrate_index} not in SMILES map",
            )
            results.append(result)
            continue

        smiles = smiles_map[substrate_index]

        if (i + 1) % 50 == 0:
            logger.info("Processing %d/%d", i + 1, len(dock_records))

        try:
            result = process_single(dock_index, substrate_index, smiles, raw_ligand_dir, aligned_ligand_dir)
        except Exception as e:
            result = AlignmentResult(
                dock_index=dock_index, substrate_index=substrate_index,
                success=False, error=f"Unexpected: {e}",
            )
            logger.error("Dock_Index %d: unexpected error: %s", dock_index, e)

        results.append(result)
        if result.success:
            success_count += 1
        else:
            logger.warning("Dock_Index %d: %s", dock_index, result.error)

    # Write summary
    logger.info("Writing summary to %s", summary_csv)
    with open(summary_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Dock_Index", "Substrate_Index", "success",
                          "num_atoms_smiles", "num_atoms_pdb", "num_matched", "error"])
        for r in results:
            writer.writerow([r.dock_index, r.substrate_index, r.success,
                              r.num_atoms_smiles, r.num_atoms_pdb, r.num_matched, r.error])

    logger.info("=" * 60)
    logger.info("Alignment complete: %d/%d successful", success_count, len(dock_records))
    logger.info("Aligned ligands: %s", aligned_ligand_dir)
    logger.info("Summary: %s", summary_csv)

    # Fail if success rate is too low (likely upstream data issue)
    if len(dock_records) > 0 and success_count / len(dock_records) < 0.5:
        logger.error("Alignment success rate %.1f%% is below 50%% threshold — likely upstream issue",
                      100 * success_count / len(dock_records))
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
