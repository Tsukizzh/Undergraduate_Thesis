"""Test script for step8_extract_pocket_ligand.py - small scale test."""

import sys
from pathlib import Path

# Add script directory to path
sys.path.insert(0, str(Path(__file__).parent))
from step8_extract_pocket_ligand import *

def main():
    # Paths
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    data_dir = project_dir / "data"

    mapping_csv = data_dir / "04_Step4_格式修正后数据" / "dock_index_mapping.csv"
    pdb_dir = data_dir / "01_Step1_PDB文件"
    output_base = data_dir / "08_Step8_结构特征" / "str_tmp_data"
    pocket_dir = output_base / "pocket"
    raw_ligand_dir = output_base / "raw_ligand"

    # Ensure output directories exist
    pocket_dir.mkdir(parents=True, exist_ok=True)
    raw_ligand_dir.mkdir(parents=True, exist_ok=True)

    # Load mapping
    mapping = load_dock_mapping(mapping_csv)
    print(f"Total records: {len(mapping)}")

    # Test first 5
    results = []
    for i, row in enumerate(mapping[:5]):
        dock_index = int(row['Dock_Index'])
        pdb_id = row['pdb_id']
        ligand_ccd = row['ligand_ccd']
        print(f"Processing {i+1}/5: Dock_Index={dock_index}, PDB={pdb_id}, Ligand={ligand_ccd}")

        result = extract_single(
            dock_index=dock_index,
            pdb_id=pdb_id,
            ligand_ccd=ligand_ccd,
            pdb_dir=pdb_dir,
            pocket_dir=pocket_dir,
            ligand_dir=raw_ligand_dir
        )
        results.append(result)

        if result.success:
            print(f"  SUCCESS: pocket={result.pocket_atoms} atoms, ligand={result.ligand_atoms} atoms")
        else:
            print(f"  FAILED: {result.error}")

    print(f"\nSummary: {sum(1 for r in results if r.success)}/5 successful")

    # Verify output files
    print("\n--- Verification ---")
    for r in results:
        if r.success:
            pocket_file = pocket_dir / f"{r.dock_index}.pdb"
            ligand_file = raw_ligand_dir / f"{r.dock_index}.sdf"
            print(f"Dock_Index {r.dock_index}:")
            print(f"  Pocket: {pocket_file.exists()} ({pocket_file.stat().st_size if pocket_file.exists() else 0} bytes)")
            print(f"  Ligand: {ligand_file.exists()} ({ligand_file.stat().st_size if ligand_file.exists() else 0} bytes)")

if __name__ == "__main__":
    main()
