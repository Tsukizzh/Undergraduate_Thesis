#!/usr/bin/env python3
"""
Multi-source verification engine for substrate classifications.

Combines:
1. Structural SMARTS validators (15 classes)
2. Ring topology / carbon count analysis (terpenoids, steroids, macrolides)
3. ClassyFire crosswalk (from classyfire_full_results.csv)
4. NPClassifier raw cache
5. Consensus scoring → evidence tiers

Output: verified_classifications.csv with evidence tier and all sources.

Usage:
    python verify_classifications.py
"""

import csv
import json
import sys
import re
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "05_底物分类"
INPUT_CSV = DATA_DIR / "substrate_classifications_corrected.csv"
NPC_CACHE = DATA_DIR / "npclassifier_cache.json"
CF_RESULTS = DATA_DIR / "classyfire_full_results.csv"
OUTPUT_CSV = DATA_DIR / "verified_classifications.csv"

VALID_LABELS = frozenset({
    "Alkaloid", "Amino_acid", "Fatty_acid", "Phenylpropanoid",
    "Steroid", "Diterpenoid", "Triterpenoid", "Polyketide",
    "Sesquiterpenoid", "Monoterpenoid", "Terpenoid_other",
    "Flavonoid", "Macrolide", "Coumarin", "Unclassified",
})

# ═══════════════════════════════════════════════════════════════════════
# 1. SMARTS VALIDATOR LIBRARY
# ═══════════════════════════════════════════════════════════════════════

def _compile(smarts_list):
    return [Chem.MolFromSmarts(s) for s in smarts_list if s]

# --- Steroid: use ring topology (fused 6-6-6-5) ---
def check_steroid_topology(mol):
    """Check for steroid-like fused 6-6-6-5 tetracyclic ring system."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    ring_sizes = [len(r) for r in rings]
    six_rings = [r for r, s in zip(rings, ring_sizes) if s == 6]
    five_rings = [r for r, s in zip(rings, ring_sizes) if s == 5]
    if len(six_rings) < 3 or len(five_rings) < 1:
        return False
    # Check fusion: at least 3 six-rings sharing edges, plus 1 five-ring fused
    for i in range(len(six_rings)):
        fused_6 = [six_rings[j] for j in range(len(six_rings)) if j != i and len(six_rings[i] & six_rings[j]) >= 2]
        if len(fused_6) >= 2:
            all_6 = six_rings[i] | fused_6[0] | fused_6[1]
            for fr in five_rings:
                if len(fr & all_6) >= 2:
                    return True
    return False

SMARTS = {
    "Coumarin": {
        "definitive": _compile([
            "O=c1ccc2ccccc2o1",
            "O=c1oc2ccccc2cc1",
            "c1cc2OC(=O)C=Cc2cc1",
        ]),
    },
    "Flavonoid": {
        "definitive": _compile([
            "O=C1CC(c2ccccc2)Oc2ccccc21",       # flavanone
            "O=c1cc(-c2ccccc2)oc2ccccc12",       # flavone
            "O=C(/C=C/c1ccccc1)c1ccc(O)cc1",     # chalcone (hydroxylated)
            "O=C(/C=C/c1ccccc1)c1ccccc1",         # chalcone
        ]),
        "supporting": _compile([
            "OC1Cc2ccccc2OC1c1ccccc1",            # flavan-3-ol
        ]),
    },
    "Macrolide": {
        # Use topology: ring size >= 12 with ester in ring
        "definitive": [],
    },
    "Fatty_acid": {
        "definitive": _compile([
            "[CX3](=O)[OX1H0-,OX2H1]",           # free carboxylate
        ]),
    },
    "Amino_acid": {
        "definitive": _compile([
            "[NX3,NX4+][CH]([*])[CX3](=O)[OX1H0-,OX2H1]",  # alpha-AA
            "[NX3,NX4+]CC(=O)[OX1H0-,OX2H1]",               # glycine-like
        ]),
        "supporting": _compile([
            "[NX3][CH](Cc1c[nH]c2ccccc12)C(=O)",  # tryptophan-like
        ]),
    },
    "Alkaloid": {
        "supporting": _compile([
            "c1ccc2[nH]ccc2c1",                    # indole
            "c1ccc2c(c1)[nH]cc2",                  # indole variant
            "c1ccc2c(c1)CCN2",                     # indoline
            "n1ccccc1",                            # pyridine
            "N1CCCCC1",                            # piperidine
            "c1ncnc2[nH]cnc12",                    # purine
            "c1ncc2ccccc2c1",                      # quinoline
            "c1cnc2ccccc2c1",                      # isoquinoline
            "N1CCCC1",                             # pyrrolidine
            "C1CC2CCC1N2",                         # tropane-like
        ]),
        "exclusions": _compile([
            "[N+](=O)[O-]",                        # nitro (not alkaloid indicator)
            "NC(=O)Nc1ccccc1",                     # phenylurea
        ]),
    },
    "Phenylpropanoid": {
        "supporting": _compile([
            "c1ccccc1/C=C/C(=O)",                  # cinnamate
            "c1ccccc1C=CC(=O)",                    # cinnamate (unspecified)
            "c1ccccc1/C=C/c2ccccc2",               # stilbene
            "c1ccc(CCc2ccccc2)cc1",                # dihydrostilbene
            "O=c1ccc2ccccc2o1",                    # coumarin
        ]),
    },
    "Polyketide": {
        "supporting": _compile([
            "O=C1c2ccccc2C(=O)c2ccccc21",          # anthraquinone
            "O=c1c(O)c2ccccc2oc2ccccc12",          # xanthone
        ]),
    },
}

# ═══════════════════════════════════════════════════════════════════════
# 2. CLASSYFIRE CROSSWALK
# ═══════════════════════════════════════════════════════════════════════

CF_CROSSWALK = {
    "Steroids and steroid derivatives": "Steroid",
    "Monoterpenoids": "Monoterpenoid",
    "Sesquiterpenoids": "Sesquiterpenoid",
    "Diterpenoids": "Diterpenoid",
    "Triterpenoids": "Triterpenoid",
    "Prenol lipids": "Terpenoid_other",
    "Terpene lactones": "Terpenoid_other",
    "Phenylpropanoids and polyketides": "Phenylpropanoid",  # broad, needs context
    "Cinnamic acids and derivatives": "Phenylpropanoid",
    "Coumarins and derivatives": "Coumarin",
    "Stilbenes": "Phenylpropanoid",
    "Lignans, neolignans and related compounds": "Phenylpropanoid",
    "Flavonoids": "Flavonoid",
    "Flavones": "Flavonoid",
    "Flavanones": "Flavonoid",
    "Isoflavonoids": "Flavonoid",
    "Chalcones and dihydrochalcones": "Flavonoid",
    "Alkaloids and derivatives": "Alkaloid",
    "Indoles and derivatives": "Alkaloid",
    "Amino acids, peptides, and analogues": "Amino_acid",
    "Alpha amino acids and derivatives": "Amino_acid",
    "Fatty acids and conjugates": "Fatty_acid",
    "Long-chain fatty acids": "Fatty_acid",
    "Very long-chain fatty acids": "Fatty_acid",
    "Medium-chain fatty acids": "Fatty_acid",
    "Eicosanoids": "Fatty_acid",
    "Lineolic acids and derivatives": "Fatty_acid",
    "Fatty acid esters": "Fatty_acid",
    "Fatty alcohols": "Fatty_acid",
    "Polyketides": "Polyketide",
    "Macrolides and analogues": "Macrolide",
    "Macrolactams": "Macrolide",
    "Anthraquinones": "Polyketide",
    "Carbohydrates and carbohydrate conjugates": "Unclassified",
    "Nucleosides, nucleotides, and analogues": "Unclassified",
    "Tetrapyrroles and derivatives": "Unclassified",
    "Organohalogen compounds": "Unclassified",
    "Benzene and substituted derivatives": "Unclassified",
    "Organic oxygen compounds": "Unclassified",
    "Organic nitrogen compounds": "Unclassified",
    "Organosulfur compounds": "Unclassified",
    "Hydrocarbons": "Unclassified",
    "Organopnictogen compounds": "Unclassified",
    "Homogeneous non-metal compounds": "Unclassified",
    "Phenol ethers": "Unclassified",
    "Anisoles": "Unclassified",
    "Benzoic acids and derivatives": "Phenylpropanoid",  # needs ArOH check
    "Phenols": "Phenylpropanoid",
    "Naphthalenes": "Unclassified",
    "Dicarboxylic acids and derivatives": "Unclassified",
    "Keto acids and derivatives": "Unclassified",
    "Toluenes": "Unclassified",
    "Pyridines and derivatives": "Alkaloid",
    "Imidazopyrimidines": "Alkaloid",
    "Diazines": "Alkaloid",
    "Lactones": "Polyketide",
    "2-arylbenzofuran flavonoids": "Flavonoid",
    "Benzofurans": "Phenylpropanoid",
}


def map_classyfire(cf_row):
    """Map ClassyFire result to our 15 classes. Try direct_parent > class > superclass."""
    if not cf_row or cf_row.get("cf_status") != "found":
        return None
    for field in ["cf_direct_parent", "cf_class", "cf_subclass", "cf_superclass"]:
        val = cf_row.get(field, "")
        if val and val in CF_CROSSWALK:
            return CF_CROSSWALK[val], field, val
    return None


# ═══════════════════════════════════════════════════════════════════════
# 3. STRUCTURAL VALIDATORS
# ═══════════════════════════════════════════════════════════════════════

def safe_mol(smiles):
    if not smiles or "*" in smiles:
        return None
    return Chem.MolFromSmiles(smiles)


def has_any(mol, patterns):
    return any(mol.HasSubstructMatch(p) for p in patterns if p is not None)


def has_nitrogen(mol):
    return any(a.GetAtomicNum() == 7 for a in mol.GetAtoms())


def has_ring_nitrogen(mol):
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 7 and a.IsInRing():
            return True
    return False


def count_carbons(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)


def count_heteroatoms(mol):
    return sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))


def has_macrolactone(mol, min_ring=12):
    """Check for macrocyclic lactone (ester bond in ring >= min_ring)."""
    ri = mol.GetRingInfo()
    ester_pat = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if ester_pat is None:
        return False
    ester_matches = mol.GetSubstructMatches(ester_pat)
    for ring_atoms in ri.AtomRings():
        if len(ring_atoms) >= min_ring:
            ring_set = set(ring_atoms)
            for match in ester_matches:
                if match[0] in ring_set and match[2] in ring_set:
                    return True
    return False


def structural_evidence(mol, label):
    """Return (evidence_type, detail) for structural validation."""
    if mol is None:
        return None, ""

    if label == "Steroid":
        if check_steroid_topology(mol):
            return "S+", "steroid 6-6-6-5 topology"
        c = count_carbons(mol)
        if c >= 17 and c <= 35:
            return "S~", f"steroid-like carbon count ({c})"
        return "S-", f"no steroid topology (C={c})"

    if label == "Coumarin":
        pats = SMARTS.get("Coumarin", {}).get("definitive", [])
        if has_any(mol, pats):
            return "S+", "coumarin core"
        return "S-", "no coumarin core"

    if label == "Flavonoid":
        pats_d = SMARTS.get("Flavonoid", {}).get("definitive", [])
        pats_s = SMARTS.get("Flavonoid", {}).get("supporting", [])
        if has_any(mol, pats_d):
            return "S+", "flavonoid skeleton"
        if has_any(mol, pats_s):
            return "S~", "flavonoid-like supporting"
        return "S-", "no flavonoid skeleton"

    if label == "Macrolide":
        if has_macrolactone(mol, 12):
            return "S+", "macrocyclic lactone (>=12)"
        return "S-", "no macrolactone"

    if label == "Fatty_acid":
        pats = SMARTS.get("Fatty_acid", {}).get("definitive", [])
        if has_any(mol, pats):
            return "S+", "has free carboxylate"
        return "S-", "no free carboxylate"

    if label == "Amino_acid":
        pats_d = SMARTS.get("Amino_acid", {}).get("definitive", [])
        pats_s = SMARTS.get("Amino_acid", {}).get("supporting", [])
        if has_any(mol, pats_d):
            return "S+", "alpha-amino acid motif"
        if has_any(mol, pats_s):
            return "S~", "amino acid derivative"
        return "S-", "no amino acid motif"

    if label == "Alkaloid":
        if not has_nitrogen(mol):
            return "S-", "no nitrogen"
        pats_s = SMARTS.get("Alkaloid", {}).get("supporting", [])
        pats_e = SMARTS.get("Alkaloid", {}).get("exclusions", [])
        if has_any(mol, pats_e) and not has_any(mol, pats_s):
            return "S-", "exclusion pattern without alkaloid scaffold"
        if has_any(mol, pats_s):
            return "S~", "alkaloid subtype scaffold"
        if has_ring_nitrogen(mol):
            return "S~", "ring nitrogen"
        return "S-", "nitrogen but no alkaloid scaffold"

    if label == "Phenylpropanoid":
        pats_s = SMARTS.get("Phenylpropanoid", {}).get("supporting", [])
        if has_any(mol, pats_s):
            return "S~", "phenylpropanoid family motif"
        return "S-", "no phenylpropanoid motif"

    if label == "Polyketide":
        pats_s = SMARTS.get("Polyketide", {}).get("supporting", [])
        if has_any(mol, pats_s):
            return "S~", "polyketide subtype motif"
        return None, "no polyketide proof (class lacks universal skeleton)"

    # Terpenoids: use carbon count heuristic
    if label in ("Monoterpenoid", "Sesquiterpenoid", "Diterpenoid", "Triterpenoid", "Terpenoid_other"):
        c = count_carbons(mol)
        h = count_heteroatoms(mol)
        bins = {
            "Monoterpenoid": (8, 13),
            "Sesquiterpenoid": (13, 18),
            "Diterpenoid": (18, 23),
            "Triterpenoid": (25, 36),
            "Terpenoid_other": (5, 50),
        }
        lo, hi = bins.get(label, (0, 100))
        if lo <= c <= hi and h <= c * 0.5:
            return "S~", f"terpene-like (C={c}, hetero={h})"
        if label == "Terpenoid_other":
            return None, f"terpenoid_other catch-all (C={c})"
        return "S-", f"carbon count {c} outside {lo}-{hi} for {label}"

    if label == "Unclassified":
        return None, "unclassified sink"

    return None, "no validator"


# ═══════════════════════════════════════════════════════════════════════
# 4. CONSENSUS ENGINE
# ═══════════════════════════════════════════════════════════════════════

def determine_tier(label, struct_ev, cf_label, npc_source):
    """
    Assign evidence tier based on consensus.
    Returns: tier (1-4), reason
    """
    has_struct_proof = struct_ev in ("S+",)
    has_struct_support = struct_ev in ("S+", "S~")
    has_struct_contra = struct_ev == "S-"
    has_cf = cf_label is not None
    cf_agrees = cf_label == label if cf_label else False
    npc_strong = npc_source in ("superclass",)

    # Tier 1: definitive structure + at least one external agrees
    if has_struct_proof and (cf_agrees or npc_strong):
        return 1, "structural_proof + external_agree"
    if has_struct_proof and label in ("Coumarin", "Macrolide", "Steroid"):
        return 1, "structural_proof (definitive class)"

    # Tier 2: two independent sources agree
    if cf_agrees and npc_strong and not has_struct_contra:
        return 2, "CF + NPC superclass agree"
    if has_struct_support and cf_agrees:
        return 2, "structural_support + CF agree"
    if has_struct_support and npc_strong:
        return 2, "structural_support + NPC superclass"

    # Tier 2b: strong structural + one source for reliable classes
    if has_struct_proof:
        return 2, "structural_proof alone"
    if has_struct_support and label in ("Fatty_acid", "Amino_acid") and npc_strong:
        return 2, f"structural_support + NPC for {label}"

    # Tier 3: one source or weak agreement
    if cf_agrees and not has_struct_contra:
        return 3, "CF agrees only"
    if npc_strong and not has_struct_contra:
        return 3, "NPC superclass only"
    if has_struct_support:
        return 3, "structural_support only"

    # Tier 4: no support or contradiction
    if has_struct_contra:
        return 4, f"structural contradiction ({struct_ev})"
    if label == "Unclassified":
        return 3, "unclassified sink (no verification needed)"
    return 4, "insufficient evidence"


def main():
    rows = list(csv.DictReader(open(INPUT_CSV, "r", encoding="utf-8-sig")))
    npc_cache = json.load(open(NPC_CACHE, "r", encoding="utf-8"))

    # Load ClassyFire results if available
    cf_data = {}
    if CF_RESULTS.exists():
        for r in csv.DictReader(open(CF_RESULTS, "r", encoding="utf-8")):
            cf_data[r["global_compound_id"]] = r
        print(f"Loaded ClassyFire results: {len(cf_data)}")
    else:
        print("WARNING: ClassyFire results not found, running without CF data")

    print(f"Loaded {len(rows)} compounds")

    out_rows = []
    tier_counts = Counter()
    class_tier = {}  # label -> Counter of tiers

    for row in rows:
        cid = row["global_compound_id"]
        smi = row["canonical_smiles"]
        label = row["corrected_label"]
        npc_source = row.get("label_source", "")
        mol = safe_mol(smi)

        # Structural evidence
        s_type, s_detail = structural_evidence(mol, label)

        # ClassyFire evidence
        cf_row = cf_data.get(cid, {})
        cf_mapped = map_classyfire(cf_row)
        cf_label = cf_mapped[0] if cf_mapped else None
        cf_field = cf_mapped[1] if cf_mapped else ""
        cf_raw = cf_mapped[2] if cf_mapped else ""

        # NPC raw
        npc_entry = npc_cache.get(smi, {})
        npc_sc = npc_entry.get("superclass", [])
        npc_cls = npc_entry.get("class", [])

        # Consensus
        tier, tier_reason = determine_tier(label, s_type, cf_label, npc_source)

        tier_counts[tier] += 1
        class_tier.setdefault(label, Counter())[tier] += 1

        out_row = dict(row)
        out_row["structural_evidence"] = s_type or ""
        out_row["structural_detail"] = s_detail
        out_row["cf_label"] = cf_label or ""
        out_row["cf_raw"] = cf_raw
        out_row["cf_field"] = cf_field
        out_row["cf_status"] = cf_row.get("cf_status", "")
        out_row["evidence_tier"] = tier
        out_row["tier_reason"] = tier_reason
        out_row["needs_review"] = "yes" if tier >= 4 else "no"
        out_rows.append(out_row)

    # Write output
    fieldnames = list(rows[0].keys()) + [
        "structural_evidence", "structural_detail",
        "cf_label", "cf_raw", "cf_field", "cf_status",
        "evidence_tier", "tier_reason", "needs_review",
    ]
    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(out_rows)

    # Summary
    print(f"\n{'='*60}")
    print(f"Wrote: {OUTPUT_CSV}")
    print(f"\n--- Evidence Tier Distribution ---")
    for t in sorted(tier_counts.keys()):
        desc = {1: "definitive", 2: "strong", 3: "moderate", 4: "NEEDS REVIEW"}
        print(f"  Tier {t} ({desc.get(t,'?'):<14}): {tier_counts[t]:>5} ({100*tier_counts[t]/len(rows):.1f}%)")

    review_count = sum(1 for r in out_rows if r["needs_review"] == "yes")
    auto_count = len(rows) - review_count
    print(f"\n  Auto-verified: {auto_count} ({100*auto_count/len(rows):.1f}%)")
    print(f"  Needs review:  {review_count} ({100*review_count/len(rows):.1f}%)")

    print(f"\n--- Per-class Tier Distribution ---")
    dist = Counter(r["corrected_label"] for r in rows)
    for label, total in dist.most_common():
        tiers = class_tier.get(label, Counter())
        t1 = tiers.get(1, 0)
        t2 = tiers.get(2, 0)
        t3 = tiers.get(3, 0)
        t4 = tiers.get(4, 0)
        pct_ok = 100 * (t1 + t2) / total if total else 0
        print(f"  {label:<20} T1={t1:>3} T2={t2:>3} T3={t3:>3} T4={t4:>3}  verified={pct_ok:.0f}%")

    # Show CF coverage
    cf_found = sum(1 for r in out_rows if r["cf_status"] == "found")
    cf_agrees = sum(1 for r in out_rows if r["cf_label"] == r["corrected_label"])
    cf_disagrees = sum(1 for r in out_rows if r["cf_label"] and r["cf_label"] != r["corrected_label"])
    print(f"\n--- ClassyFire Coverage ---")
    print(f"  Found:    {cf_found}")
    print(f"  Agrees:   {cf_agrees}")
    print(f"  Disagrees: {cf_disagrees}")


if __name__ == "__main__":
    main()
