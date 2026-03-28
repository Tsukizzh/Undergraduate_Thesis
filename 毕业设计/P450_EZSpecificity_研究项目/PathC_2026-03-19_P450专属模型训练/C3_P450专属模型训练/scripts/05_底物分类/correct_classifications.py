#!/usr/bin/env python3
"""
Deterministic correction layer on top of NPClassifier-derived substrate labels.

Strategy (Codex Round 1-2 consensus):
1. Manual overrides for known data/SMILES errors
2. Re-derive from raw NPC cache with EXPANDED superclass mapping (fixes parenthetical name mismatch)
3. Structural vetoes:
   - Alkaloid: must have nitrogen
   - Fatty_acid: must have free carboxylate + aliphatic chain >=4C
   - Amino_acid: reject multi-peptide/macrocycle; DKP → Alkaloid only if NPC agrees
   - Phenylpropanoid: pathway-only without superclass → veto then rescue with family SMARTS
4. Unclassified rescue: re-read raw NPC class/superclass
5. Phenylpropanoid family rescue: cinnamate, coumarin, stilbene, lignan SMARTS (only for truly unresolved)

Usage:
    python correct_classifications.py [--dry-run]
"""

import csv
import json
import re
import argparse
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Tuple

from rdkit import Chem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "05_底物分类"
INPUT_CSV = DATA_DIR / "substrate_classifications_final.csv"
CACHE_JSON = DATA_DIR / "npclassifier_cache.json"
OUTPUT_CSV = DATA_DIR / "substrate_classifications_corrected.csv"

VALID_LABELS = frozenset({
    "Alkaloid", "Amino_acid", "Fatty_acid", "Phenylpropanoid",
    "Steroid", "Diterpenoid", "Triterpenoid", "Polyketide",
    "Sesquiterpenoid", "Monoterpenoid", "Terpenoid_other",
    "Flavonoid", "Macrolide", "Coumarin", "Unclassified",
})

# ── SMARTS patterns ─────────────────────────────────────────────────────
PAT_FREE_CARBOXYLATE = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
PAT_AMIDE_BOND = Chem.MolFromSmarts("[NX3][CX3](=O)[#6]")
PAT_DKP_1 = Chem.MolFromSmarts("O=C1CNC(=O)CN1")
PAT_DKP_2 = Chem.MolFromSmarts("O=C1[CH]NC(=O)[CH]N1")
# Natural product aromatic indicators (hydroxyl or methoxy on aromatic ring)
PAT_AROMATIC_OH = Chem.MolFromSmarts("[OH]c")
PAT_AROMATIC_OME = Chem.MolFromSmarts("COc")

# Phenylpropanoid family rescue SMARTS (for truly unresolved pathway-only cases)
PHENYLPROPANOID_RESCUE = {
    "cinnamate": [
        Chem.MolFromSmarts("c1ccccc1/C=C/C(=O)"),
        Chem.MolFromSmarts("c1ccccc1C=CC(=O)"),
        Chem.MolFromSmarts("c1ccccc1/C=C/C=O"),
        Chem.MolFromSmarts("c1ccccc1/C=C/CO"),
        Chem.MolFromSmarts("c1cc(O)ccc1/C=C/C"),
        Chem.MolFromSmarts("c1cc(OC)ccc1/C=C/C"),
    ],
    "coumarin": [
        Chem.MolFromSmarts("O=c1ccc2ccccc2o1"),
        Chem.MolFromSmarts("c1cc2OC(=O)C=Cc2cc1"),
    ],
    "stilbene": [
        Chem.MolFromSmarts("c1ccccc1/C=C/c2ccccc2"),
        Chem.MolFromSmarts("c1ccc(CCc2ccccc2)cc1"),
    ],
    "lignan": [
        Chem.MolFromSmarts("c1cc(CC(CO)Cc2ccccc2)ccc1"),
    ],
    "phenol_simple": [
        Chem.MolFromSmarts("c1cc(O)cc(OC)c1"),
        Chem.MolFromSmarts("c1cc(O)cc(CC=C)c1"),
        Chem.MolFromSmarts("c1cc(O)c(OC)c(OC)c1"),
    ],
    "phenylpropanoid_broad": [
        Chem.MolFromSmarts("c1ccccc1CCC(=O)O"),
        Chem.MolFromSmarts("c1ccccc1CCC=O"),
        Chem.MolFromSmarts("c1ccccc1CCCO"),
        Chem.MolFromSmarts("c1cc(O)ccc1CC"),
    ],
}

# ── Manual overrides (data errors only) ─────────────────────────────────
MANUAL_OVERRIDES: Dict[str, Tuple[str, str]] = {
    "CMP_G001809": ("Unclassified", "data_error: SMILES encodes glucose but name is 'all tested hydroperoxides'"),
    "CMP_G001812": ("Unclassified", "data_error: SMILES encodes I3- triiodide ion but name is glucosinolate"),
}

# ── Expanded NPC superclass mapping ─────────────────────────────────────
# Handles NPClassifier's actual naming with parenthetical annotations.
# Key: normalized superclass name → our 15-class label.
# Normalization: strip parenthetical suffix, strip trailing whitespace.
EXPANDED_SUPERCLASS_MAP = {
    # Terpenoids
    "Monoterpenoids": "Monoterpenoid",
    "Sesquiterpenoids": "Sesquiterpenoid",
    "Diterpenoids": "Diterpenoid",
    "Sesterterpenoids": "Terpenoid_other",
    "Triterpenoids": "Triterpenoid",
    "Steroids": "Steroid",
    "Jasmonoids": "Terpenoid_other",
    "Prenol lipids": "Terpenoid_other",
    # Phenylpropanoids (the key fix)
    "Phenylpropanoids": "Phenylpropanoid",          # was "Phenylpropanoids (C6-C3)"
    "Phenolic acids": "Phenylpropanoid",             # was "Phenolic acids (C6-C1)" — C6-C1, shikimate-derived
    "Stilbenoids": "Phenylpropanoid",                # was missed (dict had "Stilbenes")
    "Stilbenes": "Phenylpropanoid",
    "Phenylethanoids": "Phenylpropanoid",            # was "Phenylethanoids (C6-C2)"
    "Styrylpyrones": "Phenylpropanoid",
    "Lignans": "Phenylpropanoid",
    "Coumarins": "Coumarin",
    "Chromones": "Phenylpropanoid",
    "Cinnamic acids and derivatives": "Phenylpropanoid",
    "Phenylpropanoids": "Phenylpropanoid",
    # Flavonoids
    "Flavonoids": "Flavonoid",
    "Isoflavonoids": "Flavonoid",
    "Chalcones and dihydrochalcones": "Flavonoid",
    # Alkaloids
    "Tyrosine alkaloids": "Alkaloid",
    "Tryptophan alkaloids": "Alkaloid",
    "Pseudoalkaloids": "Alkaloid",
    "Histidine alkaloids": "Alkaloid",
    "Nicotinic acid alkaloids": "Alkaloid",
    "Ornithine alkaloids": "Alkaloid",
    "Lysine alkaloids": "Alkaloid",
    "Purine alkaloids": "Alkaloid",
    "Steroidal alkaloids": "Alkaloid",
    "Terpenoid alkaloids": "Alkaloid",
    "Anthranilic acid alkaloids": "Alkaloid",
    # Fatty acids
    "Fatty acyls": "Fatty_acid",
    "Fatty esters": "Fatty_acid",
    "Fatty amides": "Fatty_acid",
    "Fatty alcohols": "Fatty_acid",
    # Amino acids
    "Amino acids and peptides": "Amino_acid",
    "Small peptides": "Amino_acid",
    # Others
    "Polyketides": "Polyketide",
    "Cyclic polyketides": "Polyketide",
    "Macrolides": "Macrolide",
    # Unsupported → stay Unclassified
    "Nucleosides": None,
    "Carbohydrates": None,
    "Porphyrins": None,
    "Glycerolipids": None,
    "Glycerophospholipids": None,
    "Sphingolipids": None,
    "Saccharolipids": None,
    "Polyphenols": None,
    # Denylist: aromatic buckets that should NOT auto-map to Phenylpropanoid
    "Phenanthrenoids": None,     # PAH-like, needs individual review
    "Terphenyls": None,          # aromatic, not phenylpropanoid
    "Xanthones": None,           # closer to polyketide but ambiguous
    "Diphenyl ethers": None,     # not phenylpropanoid
}

# NPC class-level rescue for Unclassified compounds
NPC_CLASS_RESCUE = {
    "Anthracyclines": "Polyketide",
    "Anthraquinones": "Polyketide",
    "Naphthoquinones": "Polyketide",
    "Tryptophan and derivatives": "Amino_acid",
    "Tyrosine and derivatives": "Amino_acid",
    "Phenylalanine and derivatives": "Amino_acid",
    "Camphane monoterpenoids": "Monoterpenoid",
    "Menthane monoterpenoids": "Monoterpenoid",
    "Acyclic monoterpenoids": "Monoterpenoid",
    "Cinnamic acids and derivatives": "Phenylpropanoid",
    "Simple phenolic acids": "Phenylpropanoid",
    "Phenols": "Phenylpropanoid",
    "Monomeric stilbenes": "Phenylpropanoid",
}


# ── Helper functions ────────────────────────────────────────────────────

def normalize_npc_name(name: str) -> str:
    """Strip parenthetical suffixes: 'Phenylpropanoids (C6-C3)' → 'Phenylpropanoids'"""
    return re.sub(r"\s*\([^)]*\)\s*$", "", name).strip()


def safe_mol(smiles: str) -> Optional[Chem.Mol]:
    if not smiles or "*" in smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def has_pattern(mol, pat) -> bool:
    return pat is not None and mol.HasSubstructMatch(pat)


def count_pattern(mol, pat) -> int:
    return len(mol.GetSubstructMatches(pat)) if pat is not None else 0


def has_nitrogen(mol) -> bool:
    return any(a.GetAtomicNum() == 7 for a in mol.GetAtoms())


def has_free_carboxylate(mol) -> bool:
    return has_pattern(mol, PAT_FREE_CARBOXYLATE)


def longest_aliphatic_chain_from_carboxyl(mol) -> int:
    """Max acyclic non-aromatic carbon chain from any carboxylate C."""
    matches = mol.GetSubstructMatches(PAT_FREE_CARBOXYLATE)
    if not matches:
        return 0
    best = 0
    for match in matches:
        c_idx = match[0]
        visited = set(match)
        c_atom = mol.GetAtomWithIdx(c_idx)
        stack = []
        for nbr in c_atom.GetNeighbors():
            nid = nbr.GetIdx()
            if nid not in visited and nbr.GetAtomicNum() == 6 \
               and not nbr.GetIsAromatic() and not nbr.IsInRing():
                stack.append((nid, 1))
                visited.add(nid)
        while stack:
            idx, depth = stack.pop()
            best = max(best, depth)
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                nid = nbr.GetIdx()
                if nid not in visited and nbr.GetAtomicNum() == 6 \
                   and not nbr.GetIsAromatic() and not nbr.IsInRing():
                    visited.add(nid)
                    stack.append((nid, depth + 1))
    return best


def is_diketopiperazine(mol) -> bool:
    return has_pattern(mol, PAT_DKP_1) or has_pattern(mol, PAT_DKP_2)


def has_macrocycle(mol, min_size: int = 8) -> bool:
    return any(len(r) >= min_size for r in mol.GetRingInfo().AtomRings())


def rescue_phenylpropanoid(mol) -> Optional[str]:
    for family, patterns in PHENYLPROPANOID_RESCUE.items():
        for pat in patterns:
            if has_pattern(mol, pat):
                return family
    return None


# ── Re-derive label from raw NPC cache ──────────────────────────────────

def rederive_from_cache(cache_entry: Optional[dict]) -> Optional[Tuple[str, str]]:
    """Try to derive a label from raw NPC superclass/class with expanded mapping."""
    if not cache_entry:
        return None

    # Try superclass first (most reliable)
    scs = cache_entry.get("superclass", [])
    if isinstance(scs, list):
        for sc in scs:
            norm = normalize_npc_name(sc)
            if norm in EXPANDED_SUPERCLASS_MAP:
                label = EXPANDED_SUPERCLASS_MAP[norm]
                if label and label in VALID_LABELS:
                    return label, f"rederive_superclass: '{sc}' -> {label}"

    # Try class-level
    classes = cache_entry.get("class", [])
    if isinstance(classes, list):
        for c in classes:
            if c in NPC_CLASS_RESCUE:
                label = NPC_CLASS_RESCUE[c]
                if label in VALID_LABELS:
                    return label, f"rederive_class: '{c}' -> {label}"

    return None


# ── Core correction logic ───────────────────────────────────────────────

def correct_one(row: dict, cache: dict) -> Tuple[str, str, str]:
    """Returns (corrected_label, reason, confidence)."""
    cid = row["global_compound_id"]
    smiles = row["canonical_smiles"]
    label = row["controlled_label_final"]
    source = row.get("label_source", "")
    mol = safe_mol(smiles)
    cache_entry = cache.get(smiles)

    # Step 1: Manual overrides
    if cid in MANUAL_OVERRIDES:
        new_label, reason = MANUAL_OVERRIDES[cid]
        return new_label, f"manual: {reason}", "high"

    if mol is None:
        return label, "unchanged", "high" if source in ("name_exact", "name_keyword") else "low"

    # Step 2: Re-derive from raw NPC with expanded mapping
    # This fixes the parenthetical name mismatch (e.g. "Phenylpropanoids (C6-C3)")
    if source == "pathway":
        rederived = rederive_from_cache(cache_entry)
        if rederived:
            new_label, reason = rederived
            # Still apply structural vetoes on rederived labels
            if new_label == "Alkaloid" and not has_nitrogen(mol):
                return "Unclassified", f"{reason} BUT veto: no nitrogen", "high"
            if new_label == "Fatty_acid":
                if not has_free_carboxylate(mol):
                    return "Unclassified", f"{reason} BUT veto: no carboxylate", "high"
                if longest_aliphatic_chain_from_carboxyl(mol) < 4:
                    return "Unclassified", f"{reason} BUT veto: chain < 4C", "high"
            # Phenylpropanoid from Phenolic acids: require natural product indicators
            # (hydroxyl or methoxy on aromatic ring). Plain/halogenated benzoic acids are xenobiotics.
            if new_label == "Phenylpropanoid" and "Phenolic acids" in reason:
                has_ar_oh = has_pattern(mol, PAT_AROMATIC_OH)
                has_ar_ome = has_pattern(mol, PAT_AROMATIC_OME)
                if not (has_ar_oh or has_ar_ome):
                    return "Unclassified", f"{reason} BUT veto: no ArOH/ArOMe (likely xenobiotic)", "high"
            return new_label, reason, "high"

    # Step 3: Structural vetoes on original labels

    # 3a. Alkaloid: must have nitrogen
    if label == "Alkaloid" and not has_nitrogen(mol):
        return "Unclassified", "veto: Alkaloid requires nitrogen", "high"

    # 3b. Fatty_acid: must have free carboxylate + chain >=4C
    if label == "Fatty_acid" and source not in ("name_exact", "name_keyword"):
        if not has_free_carboxylate(mol):
            return "Unclassified", "veto: Fatty_acid requires free carboxylate", "high"
        if longest_aliphatic_chain_from_carboxyl(mol) < 4:
            return "Unclassified", "veto: Fatty_acid requires aliphatic chain >= 4C", "high"

    # 3c. Amino_acid: reject multi-peptide/macrocycle
    if label == "Amino_acid" and source not in ("name_exact", "name_keyword"):
        if is_diketopiperazine(mol):
            # Only reclassify to Alkaloid if NPC also suggests alkaloid
            npc_pws = cache_entry.get("pathway", []) if cache_entry else []
            if "Alkaloids" in npc_pws:
                return "Alkaloid", "reclass: DKP + NPC pathway=Alkaloids", "medium"
            else:
                return "Unclassified", "veto: DKP from Amino_acid, NPC does not confirm Alkaloid", "medium"
        if count_pattern(mol, PAT_AMIDE_BOND) >= 3:
            return "Unclassified", "veto: Amino_acid excludes multi-peptide compounds", "high"
        if has_macrocycle(mol, 8):
            return "Unclassified", "veto: Amino_acid excludes macrocyclic peptides", "medium"

    # 3d. Phenylpropanoid: pathway-only without superclass → veto then rescue
    if label == "Phenylpropanoid" and source == "pathway":
        # If we got here, rederive_from_cache didn't find anything → truly unresolved
        family = rescue_phenylpropanoid(mol)
        if family:
            return "Phenylpropanoid", f"pathway_rescued: {family} SMARTS match", "medium"
        return "Unclassified", "veto: pathway-only Phenylpropanoid not rescued", "high"

    # Step 4: Rescue Unclassified
    if label == "Unclassified":
        rescued = rederive_from_cache(cache_entry)
        if rescued:
            new_label, reason = rescued
            # Validate rescued labels structurally
            if new_label == "Fatty_acid" and not (has_free_carboxylate(mol) and
                                                   longest_aliphatic_chain_from_carboxyl(mol) >= 4):
                return "Unclassified", f"rescue_rejected: {reason} but failed structural check", "medium"
            if new_label == "Alkaloid" and not has_nitrogen(mol):
                return "Unclassified", f"rescue_rejected: {reason} but no nitrogen", "medium"
            return new_label, reason, "medium"

    return label, "unchanged", "high"


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")

    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    rows = list(csv.DictReader(open(INPUT_CSV, "r", encoding="utf-8-sig")))
    cache = json.load(open(CACHE_JSON, "r", encoding="utf-8"))
    print(f"Loaded {len(rows)} compounds, {len(cache)} cache entries")

    out_rows = []
    changes = []
    reason_counts = Counter()
    conf_counts = Counter()

    for row in rows:
        corrected, reason, conf = correct_one(row, cache)
        row_out = dict(row)
        row_out["corrected_label"] = corrected
        row_out["correction_reason"] = reason
        row_out["confidence"] = conf
        out_rows.append(row_out)
        reason_counts[reason] += 1
        conf_counts[conf] += 1
        if corrected != row["controlled_label_final"]:
            changes.append((row["global_compound_id"], row["canonical_name"][:40],
                           row["controlled_label_final"], corrected, reason))

    print(f"\n{'='*70}")
    print(f"Total changes: {len(changes)}")

    print(f"\n--- Confidence ---")
    for c in ("high", "medium", "low"):
        print(f"  {c:<8} {conf_counts.get(c, 0):>5}")

    print(f"\n--- Correction reasons (non-unchanged, top 25) ---")
    for reason, cnt in reason_counts.most_common(30):
        if reason != "unchanged":
            print(f"  {cnt:>4}  {reason}")

    print(f"\n--- Transitions ---")
    trans = Counter()
    for _, _, old, new, _ in changes:
        trans[f"{old} -> {new}"] += 1
    for t, cnt in trans.most_common():
        print(f"  {cnt:>4}  {t}")

    print(f"\n--- New label distribution ---")
    final = Counter(r["corrected_label"] for r in out_rows)
    for label, cnt in final.most_common():
        print(f"  {label:<20} {cnt:>5}")

    print(f"\n--- Sample changes (first 3 per transition) ---")
    by_trans = {}
    for cid, name, old, new, reason in changes:
        key = f"{old} -> {new}"
        by_trans.setdefault(key, []).append((cid, name, reason))
    for key, items in sorted(by_trans.items()):
        print(f"\n  [{key}] ({len(items)} total)")
        for cid, name, reason in items[:3]:
            print(f"    {cid}: {name} | {reason}")

    if not args.dry_run:
        fieldnames = list(rows[0].keys()) + ["corrected_label", "correction_reason", "confidence"]
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out_rows)
        print(f"\nWrote: {OUTPUT_CSV}")
    else:
        print("\n[DRY RUN] No file written.")


if __name__ == "__main__":
    main()
