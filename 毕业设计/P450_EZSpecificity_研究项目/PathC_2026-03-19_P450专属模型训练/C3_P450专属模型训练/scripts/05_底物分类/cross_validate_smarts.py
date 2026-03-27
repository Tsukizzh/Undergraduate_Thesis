#!/usr/bin/env python3
"""
Cross-validate NPClassifier results with SMARTS structural rules.

For each compound, runs structural checks and compares with NPClassifier label.
Reports: confirmed / contradicted / no_opinion per class.

Usage:
    python cross_validate_smarts.py
"""

import csv
import sys
from pathlib import Path
from collections import Counter, defaultdict

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
INPUT_FILE = PROJECT_DIR / "data" / "classifications" / "substrate_classifications_final.csv"

# ── SMARTS patterns for structural verification ───────────────────────
# Each entry: class_name -> list of (SMARTS, description)
# A compound is "confirmed" if ANY pattern matches.

SMARTS_RULES = {
    "Coumarin": [
        ("O=c1ccc2ccccc2o1", "coumarin core (2H-chromen-2-one)"),
        ("O=c1oc2ccccc2cc1", "coumarin core variant"),
        ("o1c(=O)c2ccccc2c1", "coumarin aromatic form"),
    ],
    "Flavonoid": [
        ("O=C1CC(c2ccccc2)Oc2ccccc21", "flavanone core"),
        ("O=c1cc(oc2ccccc12)-c1ccccc1", "flavone core"),
        ("O=c1c(oc2ccccc12)-c1ccccc1", "isoflavone core"),
        ("OC1Cc2ccccc2OC1c1ccccc1", "flavan-3-ol core"),
        ("O=C(/C=C/c1ccccc1)c1ccccc1", "chalcone core"),
        # More permissive: chromanone + phenyl
        ("O=C1CCOc2ccccc21", "chromanone (broad flavanone)"),
    ],
    "Amino_acid": [
        ("[NX3,NX4+;!$([N]C=O)][CX4H]([*])[CX3](=O)[OX2H1,OX1-]", "alpha-amino acid"),
        ("[NX3,NX4+;!$([N]C=O)][CH2][CX3](=O)[OX2H1,OX1-]", "glycine-like"),
        ("[NX3][CX4H](Cc1c[nH]c2ccccc12)[CX3](=O)[O]", "tryptophan-like"),
    ],
    "Fatty_acid": [
        ("[CX4,CX3;!a][CX4,CX3;!a][CX4,CX3;!a][CX4,CX3;!a][CX4,CX3;!a][CX3](=O)[OX2H1,OX1-]",
         "C5+ chain + carboxylate"),
    ],
    "Macrolide": [
        # Detected via ring analysis, not SMARTS alone
    ],
    "Steroid": [
        # Detected via ring topology analysis, not SMARTS alone
    ],
}

# ── Exclusion patterns (things that should NOT be in a class) ─────────
EXCLUSION_RULES = {
    "Fatty_acid": [
        # Exclude if contains steroid-like fused ring system
        ("[R3]", "polycyclic (unlikely pure fatty acid)"),
    ],
    "Amino_acid": [
        # Exclude if has multiple amide bonds (peptide)
        ("[NX3][CX3](=O)[NX3][CX3](=O)", "multiple amide bonds (peptide)"),
    ],
}


def check_smarts(mol, label):
    """Check if mol matches SMARTS patterns for the given label.
    Returns: 'confirmed', 'contradicted', or 'no_opinion'
    """
    patterns = SMARTS_RULES.get(label, [])
    if not patterns:
        return "no_opinion"

    for smarts_str, desc in patterns:
        pat = Chem.MolFromSmarts(smarts_str)
        if pat is None:
            continue
        if mol.HasSubstructMatch(pat):
            # Check exclusion rules
            excl = EXCLUSION_RULES.get(label, [])
            excluded = False
            for excl_smarts, excl_desc in excl:
                excl_pat = Chem.MolFromSmarts(excl_smarts)
                if excl_pat and mol.HasSubstructMatch(excl_pat):
                    excluded = True
                    break
            if not excluded:
                return "confirmed"
    return "not_matched"


def check_steroid_topology(mol):
    """Check for steroid-like fused ring system (6-6-6-5 or similar)."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    if len(rings) < 4:
        return False

    # Find fused ring groups
    ring_sizes = sorted([len(r) for r in rings])

    # Count 6-membered and 5-membered rings
    six_rings = sum(1 for r in rings if len(r) == 6)
    five_rings = sum(1 for r in rings if len(r) == 5)

    # Steroid: at least 3 six-membered + 1 five-membered, all fused
    if six_rings >= 3 and five_rings >= 1:
        # Check that they're fused (share atoms)
        fused_count = 0
        for i in range(len(rings)):
            for j in range(i + 1, len(rings)):
                if rings[i] & rings[j]:
                    fused_count += 1
        if fused_count >= 3:
            return True
    return False


def check_macrolide(mol):
    """Check for macrocyclic lactone (ring size >= 12 with ester in ring)."""
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) < 12:
            continue
        # Check for ester bond in ring
        ring_set = set(ring)
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "O":
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring_set:
                        if neighbor.GetSymbol() == "C":
                            # Check for C=O
                            for n2 in neighbor.GetNeighbors():
                                if n2.GetSymbol() == "O" and n2.GetIdx() not in ring_set:
                                    bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2.GetIdx())
                                    if bond and bond.GetBondTypeAsDouble() == 2.0:
                                        return True
    return False


def check_alkaloid_heuristic(mol):
    """Heuristic: contains ring nitrogen, not amino acid / nucleoside."""
    has_ring_n = False
    ri = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "N" and ri.NumAtomRings(atom.GetIdx()) > 0:
            has_ring_n = True
            break
    return has_ring_n


def validate_compound(smiles, label):
    """Validate a single compound's classification.
    Returns: (verdict, method, detail)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "error", "parse_failed", "RDKit cannot parse SMILES"

    # Class-specific validation
    if label == "Steroid":
        if check_steroid_topology(mol):
            return "confirmed", "ring_topology", "6-6-6-5 fused ring system"
        else:
            return "not_matched", "ring_topology", "no steroid ring system found"

    elif label == "Macrolide":
        if check_macrolide(mol):
            return "confirmed", "ring_analysis", "macrocyclic lactone >= 12"
        else:
            return "not_matched", "ring_analysis", "no macrocyclic lactone found"

    elif label == "Alkaloid":
        if check_alkaloid_heuristic(mol):
            return "confirmed", "heuristic", "ring nitrogen present"
        else:
            return "not_matched", "heuristic", "no ring nitrogen"

    elif label in SMARTS_RULES and SMARTS_RULES[label]:
        result = check_smarts(mol, label)
        if result == "confirmed":
            return "confirmed", "smarts", f"matched SMARTS for {label}"
        else:
            return "not_matched", "smarts", f"no SMARTS match for {label}"

    else:
        return "no_opinion", "no_rule", f"no validation rule for {label}"


def main():
    sys.stdout.reconfigure(encoding="utf-8")

    rows = []
    with open(INPUT_FILE, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)

    print(f"Loaded {len(rows)} compounds from {INPUT_FILE.name}")
    print()

    # Run validation
    results = []
    class_stats = defaultdict(lambda: Counter())

    for row in rows:
        label = row["controlled_label_final"]
        smiles = row["canonical_smiles"]
        verdict, method, detail = validate_compound(smiles, label)
        results.append({
            "id": row["global_compound_id"],
            "name": row["canonical_name"],
            "smiles": smiles,
            "npc_label": label,
            "verdict": verdict,
            "method": method,
            "detail": detail,
        })
        class_stats[label][verdict] += 1

    # Summary
    print(f"{'Class':<20} {'Total':>6} {'Confirmed':>10} {'NotMatch':>10} "
          f"{'NoOpinion':>10} {'Error':>6} {'ConfirmRate':>12}")
    print("-" * 80)

    total_confirmed = 0
    total_with_opinion = 0
    for label in sorted(class_stats.keys(),
                        key=lambda x: -(class_stats[x]["confirmed"] + class_stats[x]["not_matched"])):
        s = class_stats[label]
        total = sum(s.values())
        confirmed = s["confirmed"]
        not_matched = s["not_matched"]
        no_opinion = s["no_opinion"]
        error = s["error"]
        has_opinion = confirmed + not_matched
        rate = f"{100 * confirmed / has_opinion:.1f}%" if has_opinion > 0 else "N/A"

        total_confirmed += confirmed
        total_with_opinion += has_opinion

        print(f"{label:<20} {total:>6} {confirmed:>10} {not_matched:>10} "
              f"{no_opinion:>10} {error:>6} {rate:>12}")

    print("-" * 80)
    overall = f"{100 * total_confirmed / total_with_opinion:.1f}%" if total_with_opinion > 0 else "N/A"
    print(f"{'OVERALL':<20} {len(rows):>6} {total_confirmed:>10} "
          f"{total_with_opinion - total_confirmed:>10} "
          f"{len(rows) - total_with_opinion:>10} "
          f"{'':>6} {overall:>12}")

    # Show disagreements (not_matched) for classes with rules
    print()
    print("=== Disagreements (NPClassifier label vs SMARTS/topology) ===")
    for label in sorted(class_stats.keys()):
        disagreements = [r for r in results
                         if r["npc_label"] == label and r["verdict"] == "not_matched"]
        if disagreements:
            print(f"\n[{label}] {len(disagreements)} not matched:")
            for d in disagreements[:5]:
                name = d["name"][:45] or d["smiles"][:45]
                print(f"  {d['id']}: {name} ({d['method']}: {d['detail']})")
            if len(disagreements) > 5:
                print(f"  ... and {len(disagreements) - 5} more")


if __name__ == "__main__":
    main()
