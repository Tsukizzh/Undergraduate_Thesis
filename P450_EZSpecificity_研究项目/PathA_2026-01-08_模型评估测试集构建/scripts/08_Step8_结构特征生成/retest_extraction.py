"""Re-extract specific failed ligands to verify v3.0 fix works."""

import sys
from pathlib import Path

# Add script directory to path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from step8_extract_pocket_ligand import extract_single, ExtractionResult

# Directories
BASE_DIR = script_dir.parent.parent
PDB_DIR = BASE_DIR / "data" / "01_Step1_PDB文件"
OUTPUT_DIR = BASE_DIR / "data" / "08_Step8_结构特征" / "str_tmp_data"
POCKET_DIR = OUTPUT_DIR / "pocket_retest"
LIGAND_DIR = OUTPUT_DIR / "raw_ligand_retest"

# Create output dirs
POCKET_DIR.mkdir(parents=True, exist_ok=True)
LIGAND_DIR.mkdir(parents=True, exist_ok=True)

# Test cases that previously failed with 0 atoms
TEST_CASES = [
    (75, "5EAH", "5LW"),
    (76, "5EAH", "5LX"),
    (77, "5EAH", "5LY"),
    (161, "6CPP", "CAE"),
    (402, "2P85", "IND"),
    (409, "6WGW", "CAM"),
]

print("=" * 60)
print("Re-extracting failed ligands with v3.0 fix")
print("=" * 60)

for dock_index, pdb_id, ligand_ccd in TEST_CASES:
    print(f"\nDock {dock_index}: {pdb_id}/{ligand_ccd}")

    result = extract_single(
        dock_index=dock_index,
        pdb_id=pdb_id,
        ligand_ccd=ligand_ccd,
        pdb_dir=PDB_DIR,
        pocket_dir=POCKET_DIR,
        ligand_dir=LIGAND_DIR
    )

    if result.success:
        print(f"  [OK] Success! pocket_atoms={result.pocket_atoms}, ligand_atoms={result.ligand_atoms}")

        # Check the output file
        sdf_path = LIGAND_DIR / f"{dock_index}.sdf"
        if sdf_path.exists():
            with open(sdf_path, 'r') as f:
                content = f.read()
            # Count atoms in V2000 format (line after "V2000" contains atom count)
            lines = content.split('\n')
            for i, line in enumerate(lines):
                if 'V2000' in line or 'V3000' in line:
                    parts = line.split()
                    if parts:
                        print(f"  SDF file: {sdf_path.stat().st_size} bytes, MOL block line: {line.strip()[:40]}")
                    break
    else:
        print(f"  [FAIL] Failed: {result.error}")

print("\n" + "=" * 60)
print("Done!")
