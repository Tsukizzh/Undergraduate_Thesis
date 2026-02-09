"""Debug script to diagnose ligand extraction issues."""

import sys
sys.path.insert(0, str(__file__).replace('\\', '/').rsplit('/', 1)[0])

from pathlib import Path
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.Residue import Residue
from io import StringIO
from collections import Counter
# from rdkit import Chem  # Skip RDKit for now

# Test cases: PDB files where ligand extraction produced 0 atoms
TEST_CASES = [
    ("5EAH", "5LW"),   # Dock 75 - altLoc B
    ("6CPP", "CAE"),   # Dock 161 - altLoc B
    ("2P85", "IND"),   # Dock 402 - altLoc A
    ("6WGW", "CAM"),   # Dock 409 - this PDB has 2 ligands (CAH and CAM)
]

PDB_DIR = Path(__file__).parent.parent.parent / "data" / "01_Step1_PDB文件"


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


def analyze_pdb(pdb_id, ligand_ccd):
    print(f"\n{'='*60}")
    print(f"Analyzing {pdb_id} / {ligand_ccd}")
    print(f"{'='*60}")

    # Find file
    pdb_path = PDB_DIR / f"{pdb_id}.pdb"
    cif_path = PDB_DIR / f"{pdb_id}.cif"

    if pdb_path.exists():
        structure_path = pdb_path
        parser = PDBParser(QUIET=True)
    elif cif_path.exists():
        structure_path = cif_path
        parser = MMCIFParser(QUIET=True)
    else:
        print(f"ERROR: Neither {pdb_path} nor {cif_path} exists")
        return

    print(f"Using: {structure_path}")
    structure = parser.get_structure("test", str(structure_path))
    model = next(structure.get_models())

    # Find all HETATM residues
    print(f"\n--- All HETATM residues ---")
    hetatm_count = 0
    for chain in model:
        for residue in chain:
            het_flag = residue.get_id()[0]
            if het_flag.startswith("H_") or het_flag == "W":
                resname = residue.get_resname()
                if resname not in ["HOH", "WAT"]:
                    atom_count = sum(1 for _ in residue.get_atoms())
                    hetatm_count += 1
                    print(f"  Chain {chain.id}: {resname} (id={residue.get_id()}, atoms={atom_count})")
    print(f"Total HETATM residues: {hetatm_count}")

    # Find target ligand
    print(f"\n--- Searching for ligand '{ligand_ccd}' ---")
    candidates = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            if resname == ligand_ccd.upper():
                candidates.append(residue)

    print(f"Found {len(candidates)} candidate(s)")

    if not candidates:
        # Debug: show all residue names
        all_resnames = set()
        for chain in model:
            for residue in chain:
                all_resnames.add(residue.get_resname().strip().upper())
        print(f"All unique residue names: {sorted(all_resnames)}")
        return

    for i, residue in enumerate(candidates):
        print(f"\n--- Candidate {i+1}: {residue.get_resname()} ---")
        print(f"ID: {residue.get_id()}")
        print(f"Chain: {residue.get_parent().get_id()}")

        # Analyze atoms
        print(f"\nAtom analysis:")
        altloc_counts = Counter()
        atom_list = []
        for atom in residue.get_atoms():
            atom_name = atom.get_name()
            is_disordered = atom.is_disordered()

            if is_disordered:
                for a in atom.disordered_get_list():
                    alt = a.get_altloc()
                    norm_alt = _normalize_altloc(alt)
                    altloc_counts[norm_alt or "(blank)"] += 1
                    atom_list.append((atom_name, alt, a.get_coord()[:2]))
            else:
                alt = atom.get_altloc()
                norm_alt = _normalize_altloc(alt)
                altloc_counts[norm_alt or "(blank)"] += 1
                atom_list.append((atom_name, alt, atom.get_coord()[:2]))

        print(f"  Total atoms (via get_atoms): {len(list(residue.get_atoms()))}")
        print(f"  AltLoc distribution: {dict(altloc_counts)}")

        best_alt = _preferred_residue_altloc(residue)
        print(f"  Preferred altLoc: '{best_alt}'")

        # Show first few atoms
        print(f"  First 5 atoms: {atom_list[:5]}")

        # Try to generate PDB block
        print(f"\n  Generating PDB block...")

        class LigandSelect(Select):
            def __init__(self, target_residue, model_id, best_alt):
                self.target_id = target_residue.get_id()
                self.target_chain = target_residue.get_parent().get_id()
                self._model_id = model_id
                self._best_alt = best_alt

            def accept_model(self, m):
                return 1 if m.get_id() == self._model_id else 0

            def accept_chain(self, chain):
                return 1 if chain.get_id() == self.target_chain else 0

            def accept_residue(self, res):
                return 1 if res.get_id() == self.target_id else 0

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

        # Count HETATM lines in PDB block
        hetatm_lines = [l for l in pdb_block.split('\n') if l.startswith('HETATM')]
        print(f"  PDB block: {len(pdb_block)} bytes, {len(hetatm_lines)} HETATM lines")
        if hetatm_lines:
            print(f"  First HETATM: {hetatm_lines[0][:60]}...")

        # Try RDKit parsing (skip if not available)
        # mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
        # if mol is None:
        #     print(f"  RDKit: Failed to parse PDB block!")
        # else:
        #     print(f"  RDKit: Parsed OK, {mol.GetNumAtoms()} atoms")


if __name__ == "__main__":
    for pdb_id, ligand_ccd in TEST_CASES:
        analyze_pdb(pdb_id, ligand_ccd)
