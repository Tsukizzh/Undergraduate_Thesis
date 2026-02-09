"""Debug the PDB block to RDKit mol conversion."""

import sys
from pathlib import Path
from io import StringIO

# Add script directory to path
script_dir = Path(__file__).parent
sys.path.insert(0, str(script_dir))

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Residue import Residue
from collections import Counter
from rdkit import Chem

# Directories
BASE_DIR = script_dir.parent.parent
PDB_DIR = BASE_DIR / "data" / "01_Step1_PDB文件"


def _normalize_altloc(alt):
    alt = (alt or "").strip()
    if alt in {"", ".", "?"}:
        return ""
    return alt


def _preferred_residue_altloc(residue):
    counts = Counter()
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


def debug_conversion(pdb_id, ligand_ccd):
    print(f"\n{'='*60}")
    print(f"Debug: {pdb_id} / {ligand_ccd}")
    print(f"{'='*60}")

    # Parse structure
    pdb_path = PDB_DIR / f"{pdb_id}.pdb"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("test", str(pdb_path))
    model = next(structure.get_models())

    # Find ligand
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip().upper() == ligand_ccd.upper():
                print(f"Found ligand in chain {chain.id}")

                best_alt = _preferred_residue_altloc(residue)
                print(f"Preferred altLoc: '{best_alt}'")

                # Generate PDB block
                class LigandSelect(Select):
                    def __init__(self, target_residue, model_id, best_alt):
                        self.target_id = target_residue.get_id()
                        self.target_chain = target_residue.get_parent().get_id()
                        self._model_id = model_id
                        self._best_alt = best_alt

                    def accept_model(self, m):
                        return 1 if m.get_id() == self._model_id else 0

                    def accept_chain(self, c):
                        return 1 if c.get_id() == self.target_chain else 0

                    def accept_residue(self, r):
                        return 1 if r.get_id() == self.target_id else 0

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
                io.save(buf, select=LigandSelect(residue, model.get_id(), best_alt))
                pdb_block = buf.getvalue()

                # Show PDB block
                hetatm_lines = [l for l in pdb_block.split('\n') if l.startswith('HETATM')]
                print(f"PDB block has {len(hetatm_lines)} HETATM lines")
                print("First 3 lines:")
                for line in hetatm_lines[:3]:
                    print(f"  {line}")

                # Try RDKit parsing
                print("\nTrying RDKit parsing...")
                mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
                if mol is None:
                    print("  RDKit returned None!")
                else:
                    print(f"  RDKit parsed: {mol.GetNumAtoms()} atoms")

                    # Try to get conformer
                    if mol.GetNumConformers() > 0:
                        conf = mol.GetConformer()
                        print(f"  Conformer has {conf.GetNumAtoms()} atoms")

                # Try with different options
                print("\nTrying with different RDKit options...")

                # Option 1: proximity bonding
                mol2 = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, proximityBonding=True)
                if mol2:
                    print(f"  proximityBonding=True: {mol2.GetNumAtoms()} atoms")
                else:
                    print("  proximityBonding=True: None")

                # Option 2: flavor=1 (use element from residue name)
                mol3 = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, flavor=1)
                if mol3:
                    print(f"  flavor=1: {mol3.GetNumAtoms()} atoms")
                else:
                    print("  flavor=1: None")

                # Save PDB block to file for manual inspection
                debug_pdb = script_dir / f"debug_{pdb_id}_{ligand_ccd}.pdb"
                with open(debug_pdb, 'w') as f:
                    f.write(pdb_block)
                print(f"\nSaved PDB block to: {debug_pdb}")

                return  # Only process first match


# Test cases
TEST_CASES = [
    ("5EAH", "5LW"),
    ("6CPP", "CAE"),
]

for pdb_id, ligand_ccd in TEST_CASES:
    debug_conversion(pdb_id, ligand_ccd)

print("\n" + "=" * 60)
print("Done!")
