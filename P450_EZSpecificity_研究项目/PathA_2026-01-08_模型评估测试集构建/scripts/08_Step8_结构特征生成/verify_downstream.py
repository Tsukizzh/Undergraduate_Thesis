"""Verify generated pocket/ligand files are compatible with downstream code.

Standalone verification - does not depend on the full EZSpecificity environment.
Tests:
1. Pocket PDB: must contain ATOM lines parseable by PDBProtein logic
2. Ligand SDF: must be readable by RDKit SDMolSupplier
"""

from pathlib import Path
from rdkit import Chem
import numpy as np


def verify_pocket_pdb(pdb_path: Path) -> dict:
    """Simulate PDBProtein parsing to verify pocket PDB format."""
    with open(pdb_path, 'r') as f:
        block = f.read()

    atoms = []
    for line in block.splitlines():
        if line[0:6].strip() == 'ATOM':
            try:
                element_symb = line[76:78].strip().capitalize()
                if len(element_symb) == 0:
                    element_symb = line[13:14]
                atom_info = {
                    'atom_name': line[12:16].strip(),
                    'res_name': line[17:20].strip(),
                    'chain': line[21:22].strip(),
                    'res_id': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element_symb': element_symb,
                }
                atoms.append(atom_info)
            except (ValueError, IndexError) as e:
                return {'ok': False, 'error': f'Parse error at line: {line[:30]}... : {e}'}

    if len(atoms) == 0:
        return {'ok': False, 'error': 'No ATOM lines found'}

    # Check for backbone atoms
    backbone_names = {"CA", "C", "N", "O"}
    has_backbone = any(a['atom_name'] in backbone_names for a in atoms)

    # Count unique residues
    residues = set()
    for a in atoms:
        residues.add((a['chain'], a['res_id'], a['res_name']))

    return {
        'ok': True,
        'num_atoms': len(atoms),
        'num_residues': len(residues),
        'has_backbone': has_backbone,
    }


def verify_ligand_sdf(sdf_path: Path) -> dict:
    """Verify ligand SDF file is readable by RDKit.

    Uses MolFromMolBlock with Python file read to avoid RDKit's
    inability to handle non-ASCII paths in SDMolSupplier.
    NOTE: Do NOT strip() the mol_block - it removes header lines!
    """
    try:
        with open(sdf_path, 'r') as f:
            content = f.read()
        # Split on $$$$ delimiter, keep the mol block as-is (no strip!)
        mol_block = content.split("$$$$")[0]
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        if mol is None:
            return {'ok': False, 'error': 'RDKit returned None from MolBlock'}

        num_atoms = mol.GetNumAtoms()
        has_conformer = mol.GetNumConformers() > 0

        if has_conformer:
            conf = mol.GetConformer()
            positions = conf.GetPositions()
            all_zero = np.allclose(positions, 0.0)
        else:
            all_zero = True

        return {
            'ok': True,
            'num_atoms': num_atoms,
            'has_3d': has_conformer and not all_zero,
        }
    except Exception as e:
        return {'ok': False, 'error': str(e)}


def main():
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent.parent
    data_dir = project_dir / "data" / "08_Step8_结构特征" / "str_tmp_data"
    pocket_dir = data_dir / "pocket"
    ligand_dir = data_dir / "raw_ligand"

    print("=== Pocket PDB Verification (PDBProtein compatible) ===")
    pocket_ok = 0
    for i in range(5):
        pocket_file = pocket_dir / f"{i}.pdb"
        if not pocket_file.exists():
            print(f"Dock_Index {i}: FILE NOT FOUND")
            continue
        result = verify_pocket_pdb(pocket_file)
        if result['ok']:
            print(f"Dock_Index {i}: {result['num_atoms']} atoms, "
                  f"{result['num_residues']} residues, "
                  f"backbone={'YES' if result['has_backbone'] else 'NO'}")
            pocket_ok += 1
        else:
            print(f"Dock_Index {i}: FAILED - {result['error']}")

    print()
    print("=== Ligand SDF Verification (RDKit compatible) ===")
    ligand_ok = 0
    for i in range(5):
        ligand_file = ligand_dir / f"{i}.sdf"
        if not ligand_file.exists():
            print(f"Dock_Index {i}: FILE NOT FOUND")
            continue
        result = verify_ligand_sdf(ligand_file)
        if result['ok']:
            print(f"Dock_Index {i}: {result['num_atoms']} atoms, "
                  f"3D coords={'YES' if result['has_3d'] else 'NO'}")
            ligand_ok += 1
        else:
            print(f"Dock_Index {i}: FAILED - {result['error']}")

    print()
    print(f"=== Summary ===")
    print(f"Pocket PDB: {pocket_ok}/5 OK")
    print(f"Ligand SDF: {ligand_ok}/5 OK")

    if pocket_ok == 5 and ligand_ok == 5:
        print("ALL TESTS PASSED - ready for full run")
    else:
        print("SOME TESTS FAILED - investigate before full run")


if __name__ == "__main__":
    main()
