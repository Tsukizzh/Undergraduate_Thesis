"""Quick validation of Meeko installation and PDBQT conversion pipeline."""
import sys

def main():
    # 1. Meeko imports
    import meeko
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    print(f"Meeko version: {meeko.__version__}")

    from meeko import PDBQTMolecule, RDKitMolCreate
    print("All Meeko imports: OK")

    # 2. RDKit + Meeko integration: SMILES -> PDBQT
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles("c1ccccc1")  # benzene
    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    assert status == 0, f"EmbedMolecule failed with status {status}"
    AllChem.MMFFOptimizeMolecule(mol)

    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)
    pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
    assert is_ok, f"PDBQT write failed: {err}"
    assert len(pdbqt_str) > 100, "PDBQT string too short"
    print(f"SMILES->PDBQT: OK (length={len(pdbqt_str)} chars)")

    # 3. Write PDBQT to temp file, then read back with PDBQTMolecule
    import tempfile, os
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdbqt", delete=False) as f:
        f.write(pdbqt_str)
        tmp_path = f.name

    try:
        pdbqt_mol = PDBQTMolecule.from_file(tmp_path)
        rdkit_mols = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        assert len(rdkit_mols) > 0, "No molecules recovered from PDBQT"
        recovered = rdkit_mols[0]
        n_atoms = recovered.GetNumAtoms()
        print(f"PDBQT->RDKit Mol: OK (recovered {n_atoms} atoms)")
    finally:
        os.unlink(tmp_path)

    # 4. BioPython check
    from Bio.PDB import PDBParser
    print("BioPython PDBParser: OK")

    # 5. Vina CLI check
    import subprocess
    result = subprocess.run(
        ["D:\\autodock\\vina.exe", "--version"],
        capture_output=True, text=True, timeout=10
    )
    print(f"Vina CLI: {result.stdout.strip()}")

    # 6. MGLTools Python check
    result = subprocess.run(
        ["D:\\autodock\\MGLTools-1.5.7\\python.exe", "-c", "print('MGLTools Python: OK')"],
        capture_output=True, text=True, timeout=10
    )
    print(result.stdout.strip())

    print("\n=== ALL ENVIRONMENT CHECKS PASSED ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
