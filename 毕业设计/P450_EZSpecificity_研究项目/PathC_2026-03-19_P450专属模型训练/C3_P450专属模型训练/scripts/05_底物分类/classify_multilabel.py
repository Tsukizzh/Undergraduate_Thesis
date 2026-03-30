#!/usr/bin/env python3
"""
Assign strict multilabel natural-product classes to P450 substrates.

Tiers:
  - gold:   compound found in NPClassifier training data with matching label
  - auto:   high-precision automatic labels (SMARTS + NPC + ClassyFire)
  - review: evidence exists but not strong enough for auto
  - other:  no evidence for any of the 7 classes

Each class is detected independently -- no cross-class suppression.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "05_底物分类"
DEFAULT_INPUT = DATA_DIR / "substrate_classifications_final.csv"
DEFAULT_CACHE = DATA_DIR / "npclassifier_cache.json"
DEFAULT_CLASSYFIRE = DATA_DIR / "classyfire_full_results.csv"
DEFAULT_GOLD = None  # NPClassifier_dataset.tsv
DEFAULT_OUTPUT = DATA_DIR / "substrate_multilabel_7class.csv"

LABELS: Sequence[str] = (
    "Steroid",
    "Terpenoid",
    "Alkaloid",
    "Amino_acid",
    "Fatty_acid",
    "Phenylpropanoid",
    "Polyketide",
)

# ---------------------------------------------------------------------------
# Priority 2: NPC Superclass → label mapping (PRIMARY classification method)
# Keys are lowercase; lookup normalizes by stripping parenthetical suffixes.
# ---------------------------------------------------------------------------
SUPERCLASS_TO_LABELS: Dict[str, List[str]] = {
    # Steroid
    "steroids": ["Steroid"],
    # Terpenoid
    "monoterpenoids": ["Terpenoid"],
    "sesquiterpenoids": ["Terpenoid"],
    "diterpenoids": ["Terpenoid"],
    "sesterterpenoids": ["Terpenoid"],
    "triterpenoids": ["Terpenoid"],
    "carotenoids": ["Terpenoid"],       # raw: "Carotenoids (C40)"
    "apocarotenoids": ["Terpenoid"],
    "meroterpenoids": ["Terpenoid"],
    # Alkaloid
    "tryptophan alkaloids": ["Alkaloid"],
    "tyrosine alkaloids": ["Alkaloid"],
    "nicotinic acid alkaloids": ["Alkaloid"],
    "histidine alkaloids": ["Alkaloid"],
    "ornithine alkaloids": ["Alkaloid"],
    "lysine alkaloids": ["Alkaloid"],
    "anthranilic acid alkaloids": ["Alkaloid"],
    "proline alkaloids": ["Alkaloid"],
    "tetramate alkaloids": ["Alkaloid"],
    "pseudoalkaloids": ["Alkaloid"],
    "peptide alkaloids": ["Alkaloid", "Amino_acid"],
    # Amino_acid
    "small peptides": ["Amino_acid"],
    "oligopeptides": ["Amino_acid"],
    "amino acid glycosides": ["Amino_acid"],
    "\u03b2-lactams": ["Amino_acid"],
    "beta-lactams": ["Amino_acid"],     # alternate spelling
    # Fatty_acid  (Fatty acyls handled specially: carbon check)
    "fatty acids and conjugates": ["Fatty_acid"],
    "fatty acyls": ["Fatty_acid"],      # carbon check applied later
    "fatty esters": ["Fatty_acid"],
    "fatty amides": ["Fatty_acid"],
    "glycerophospholipids": ["Fatty_acid"],
    "eicosanoids": ["Fatty_acid"],
    "octadecanoids": ["Fatty_acid"],
    "docosanoids": ["Fatty_acid"],
    # Phenylpropanoid
    "phenylpropanoids": ["Phenylpropanoid"],  # raw: "Phenylpropanoids (C6-C3)"
    # "phenolic acids" removed from auto — C6-C1 compounds are often NOT
    # phenylpropanoids (spot-check: isovanillic acid, prenylhydroquinone,
    # methylthiobenzoic acid all wrong). Moved to SUPERCLASS_REVIEW.
    "phenylethanoids": ["Phenylpropanoid"],   # raw: "Phenylethanoids (C6-C2)"
    "coumarins": ["Phenylpropanoid"],
    "flavonoids": ["Phenylpropanoid"],
    "isoflavonoids": ["Phenylpropanoid"],
    "lignans": ["Phenylpropanoid"],
    "stilbenoids": ["Phenylpropanoid"],
    "styrylpyrones": ["Phenylpropanoid"],
    "phenanthrenoids": [],  # PAH, not phenylpropanoid (spot-check: phenanthrene/pyrene are wrong)
    # Polyketide
    "macrolides": ["Polyketide"],
    "cyclic polyketides": ["Polyketide"],
    "linear polyketides": ["Polyketide"],
    "aromatic polyketides": ["Polyketide"],
    "polycyclic aromatic polyketides": ["Polyketide"],
    "phloroglucinols": ["Polyketide"],
    "xanthones": ["Polyketide"],
    "naphthalenes": ["Polyketide"],
    "chromanes": ["Polyketide"],
    "alkylresorsinols": ["Polyketide"],
    "tropolones": ["Polyketide"],
    # Fatty_acid (additional lipid superclasses)
    "glycerolipids": ["Fatty_acid"],
    "spingolipids": ["Fatty_acid"],
    "fatty acyl glycosides": ["Fatty_acid"],
    # Amino_acid (additional)
    "_-lactams": ["Amino_acid"],
    "_-lactam-_-lactones": ["Amino_acid"],
    "mycosporine derivatives": ["Amino_acid"],
    # Alkaloid (additional)
    "guanidine alkaloids": ["Alkaloid"],
    # Phenylpropanoid (additional shikimate-derived)
    "diarylheptanoids": ["Phenylpropanoid"],
    "terphenyls": ["Polyketide"],  # biphenyl synthase is type III PKS (spot-check verified)
    "fluorenes": ["Phenylpropanoid"],
    "diazotetronic acids and derivatives": ["Phenylpropanoid"],
    # Terpenoid (additional)
    "polyethers": ["Terpenoid"],
    "polyprenols": ["Terpenoid"],
    # Other (empty label list → falls through to Other)
    "saccharides": [],
    "aminosugars and aminoglycosides": [],
    "nucleosides": [],
    "polyols": [],
}

# Ambiguous superclasses → send to review with hint
SUPERCLASS_REVIEW: Set[str] = {
    "diphenyl ethers",     # raw may have "(DPEs)"
    "phenolic acids",      # C6-C1; many are NOT phenylpropanoid-derived
}


def _smarts(pattern: str) -> Chem.Mol:
    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        raise ValueError(f"Invalid SMARTS: {pattern}")
    return mol


# --- SMARTS patterns ---

PAT_ALPHA_AA = _smarts("[NX3,NX4+;!$(NC=O)][CH1,CH2][CX3](=O)[OX2H1,OX1-,OX2]")
PAT_ALPHA_AA_AR = _smarts("[NX3,NX4+;!$(NC=O)][CH]([#6])[CX3](=O)[OX2H1,OX1-,OX2]")
PAT_PEPTIDE_BOND = _smarts("[NX3,NX4+][CX3](=O)[#6]")
PAT_OXIME = _smarts("[CX3]=[NX2][OX2H,OX1-]")
PAT_DKP = _smarts("O=C1NCC(=O)NC1")
PAT_ADENOSINE = _smarts("n1cnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O")

PAT_FREE_ACID = _smarts("[CX3](=O)[OX2H1,OX1-]")
PAT_ESTER = _smarts("[CX3](=O)O[#6]")
PAT_AMIDE = _smarts("[CX3](=O)N")
PAT_THIOESTER = _smarts("[CX3](=O)S[#6]")
PAT_COA_HINT = _smarts("SCCNC(=O)CCNC(=O)")
PAT_CARBONYL_OXYGENATED = _smarts("[CX3](=O)[O,N,S]")
PAT_KETONE = _smarts("[#6][CX3](=O)[#6]")

PAT_CINNAMATE = _smarts("c1ccccc1/C=C/[CX3](=O)")
PAT_CINNAMATE2 = _smarts("c1ccccc1C=CC(=O)")
PAT_HYDROCINNAMATE = _smarts("c1ccccc1CCC(=O)[O,N,S]")
PAT_STILBENE = _smarts("c1ccccc1/C=C/c2ccccc2")
PAT_COUMARIN = _smarts("O=c1ccc2ccccc2o1")
PAT_COUMARIN2 = _smarts("O=C1C=Cc2ccccc2O1")
PAT_ISOCOUMARIN = _smarts("O=c1occc2ccccc12")
PAT_FLAVONE_CORE = _smarts("O=c1cc(-c2ccccc2)oc2ccccc12")
PAT_FLAVANONE = _smarts("O=C1CC(c2ccccc2)Oc2ccccc21")
PAT_CHALCONE = _smarts("O=C(/C=C/c1ccccc1)c1ccccc1")
PAT_ISOFLAVONE = _smarts("O=c1oc2ccccc2c(-c2ccccc2)c1")
PAT_BENZOIC_ACID = _smarts("c1ccccc1C(=O)[OX2H1,OX1-]")
PAT_HYDROXYBENZOIC = _smarts("[OH]c1ccccc1C(=O)[OX2H1,OX1-]")

TRUE_PHENYLPROPANOID_PATTERNS = (
    PAT_CINNAMATE, PAT_CINNAMATE2, PAT_STILBENE,
    PAT_COUMARIN, PAT_COUMARIN2,
    PAT_FLAVONE_CORE, PAT_FLAVANONE, PAT_CHALCONE, PAT_ISOFLAVONE,
)
# C6-C1 phenolic acids count as phenylpropanoid-derived
C6C1_PATTERNS = (PAT_BENZOIC_ACID, PAT_HYDROXYBENZOIC)

PAT_ANTHRAQUINONE = _smarts("O=C1c2ccccc2C(=O)c2ccccc21")
PAT_NAPHTHOQUINONE = _smarts("O=C1C=CC(=O)c2ccccc21")
PAT_XANTHONE = _smarts("O=c1oc2ccccc2c2ccccc12")
POLYKETIDE_STRUCTURES = (PAT_ANTHRAQUINONE, PAT_NAPHTHOQUINONE, PAT_XANTHONE)

PAT_INDOLE = _smarts("c1ccc2[nH]ccc2c1")
PAT_INDOLE_GENERIC = _smarts("c1cc2cc[nH]c2cc1")  # looser indole match
PAT_PHENYLETHYLAMINE = _smarts("c1ccccc1CCN")
PAT_TETRAHYDROISOQUINOLINE = _smarts("c1ccc2CCNCC2c1")
PAT_ISOQUINOLINE = _smarts("c1cnc2ccccc2c1")
PAT_TROPANE = _smarts("C1CC2CCC1N2")
PAT_PURINE = _smarts("c1ncnc2[nH]cnc12")
PAT_PYRIDINE = _smarts("n1ccccc1")
PAT_IMIDAZOLE = _smarts("[nH]1ccnc1")
PAT_PIPERIDINE = _smarts("N1CCCCC1")
PAT_QUINOLINE = _smarts("c1ccc2ncccc2c1")

HIGH_PRECISION_ALKALOID_SCAFFOLDS = (
    PAT_TETRAHYDROISOQUINOLINE, PAT_ISOQUINOLINE,
    PAT_TROPANE, PAT_PURINE,
)
# Broadened: simple indole and phenylethylamine now count as alkaloid scaffolds
BROAD_ALKALOID_SCAFFOLDS = (
    PAT_INDOLE, PAT_INDOLE_GENERIC, PAT_PHENYLETHYLAMINE,
)
GENERIC_ALKALOID_SCAFFOLDS = (
    PAT_PYRIDINE, PAT_IMIDAZOLE, PAT_PIPERIDINE, PAT_QUINOLINE,
)

PAT_CF3 = _smarts("[CX4](F)(F)F")
PAT_SULFONE = _smarts("[#16X4](=[OX1])(=[OX1])([#6])[#6]")
PAT_SULFONAMIDE = _smarts("[#16X4](=[OX1])(=[OX1])[NX3]")
PAT_AZO = _smarts("[NX2]=[NX2]")
PAT_TRIAZOLE_1 = _smarts("n1ccnn1")
PAT_TRIAZOLE_2 = _smarts("n1nccn1")
PAT_FLUORINE = _smarts("[F]")


# --- NPClassifier vocab ---

NPC_STEROID = {"steroids", "steroid"}
NPC_TERPENOID_AUTO = {
    "monoterpenoids", "sesquiterpenoids", "diterpenoids", "sesterterpenoids",
    "triterpenoids", "meroterpenoids", "apocarotenoids", "iridoids",
    "retinoids", "carotenoids",
}
NPC_TERPENOID_REVIEW = {"terpenoids"}
NPC_ALKALOID = {
    "alkaloids", "tyrosine alkaloids", "tryptophan alkaloids", "pseudoalkaloids",
    "histidine alkaloids", "nicotinic acid alkaloids", "ornithine alkaloids",
    "lysine alkaloids", "purine alkaloids", "steroidal alkaloids",
    "terpenoid alkaloids", "anthranilic acid alkaloids", "peptide alkaloids",
    "indole diketopiperazine alkaloids",
}
NPC_AMINO = {"amino acids and peptides", "small peptides", "oligopeptides"}
NPC_FATTY_SPECIFIC = {
    "fatty acyls", "fatty esters", "fatty amides", "fatty acids and conjugates",
    "octadecanoids", "eicosanoids", "docosanoids", "amino fatty acids",
}
NPC_FATTY_REVIEW = NPC_FATTY_SPECIFIC | {"fatty acids"}
NPC_JASMONOID = {"jasmonoids", "jasmonoid"}
NPC_PHENYLPROPANOID_AUTO = {
    "cinnamic acids and derivatives", "coumarins", "flavonoids", "isoflavonoids",
    "stilbenes", "stilbenoids", "lignans",
}
NPC_PHENYLPROPANOID_REVIEW = NPC_PHENYLPROPANOID_AUTO | {
    "phenylpropanoids", "chalcones and dihydrochalcones", "chromones",
    "flavones", "flavanones", "dihydroflavonols", "flavan-3-ols", "anthocyanins",
    "shikimates and phenylpropanoids",
}
NPC_POLYKETIDE_AUTO = {
    "macrolides", "ansamycins", "tetracyclines", "cyclic polyketides",
    "linear polyketides", "anthraquinones", "naphthoquinones", "polycyclic polyketides",
}
NPC_POLYKETIDE_REVIEW = NPC_POLYKETIDE_AUTO | {"polyketides"}

CF_TERPENE_TOKENS = ("terpenoid", "terpene")
CF_POLYKETIDE_CLASSES = {
    "macrolides and analogues", "anthraquinones", "naphthoquinones",
    "xanthones", "tetracyclines",
}
CF_PHENYLPROPANOID_CLASSES = {
    "cinnamic acids and derivatives", "flavonoids", "coumarins and derivatives",
    "stilbenes", "lignans, neolignans and related compounds",
}

# Gold standard uses SUPERCLASS_TO_LABELS (same mapping as P2).
# No separate gold maps needed — load_gold_standard() normalizes the
# TSV's Super_class column and looks it up in SUPERCLASS_TO_LABELS.


@dataclass
class Decision:
    auto: bool = False
    review: bool = False
    reason: str = ""


# --- Gold standard ---

def load_gold_standard(tsv_path: Path) -> Dict[str, Set[str]]:
    """Load NPClassifier training data: SMILES -> set of our 7 labels.

    Uses the Superclass column (not Pathway) for mapping, via the same
    SUPERCLASS_TO_LABELS dictionary used by P2.  Superclass is more precise
    than Pathway (70 categories vs 7) and each Class belongs to exactly one
    Superclass, so Superclass is the optimal granularity.
    """
    gold: Dict[str, Set[str]] = {}
    if not tsv_path.exists():
        return gold
    with tsv_path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            smi_raw = (row.get("index") or "").strip()
            if not smi_raw:
                continue
            mol = Chem.MolFromSmiles(smi_raw)
            if mol is None:
                continue
            csmi = Chem.MolToSmiles(mol)
            raw_sc = (row.get("Super_class") or "").strip()
            key = superclass_lookup_key(raw_sc)
            label_list = SUPERCLASS_TO_LABELS.get(key, [])
            if label_list:
                gold.setdefault(csmi, set()).update(label_list)
    return gold


# --- Helpers ---

def normalize_term(value: str) -> str:
    return re.sub(r"\s*\([^)]*\)\s*$", "", value or "").strip().lower()

def normalize_cf(value: str) -> str:
    return (value or "").strip().lower()

def normalize_npc(entry: Optional[dict]) -> Set[str]:
    if not entry: return set()
    terms: Set[str] = set()
    for key in ("pathway", "superclass", "class"):
        for value in entry.get(key, []) or []:
            terms.add(normalize_term(value))
    return terms

def split_npc_field(value: str) -> List[str]:
    if not value: return []
    return [item.strip() for item in value.split("|") if item.strip()]

def get_npc_entry(row: dict, cache: Dict[str, dict]) -> Optional[dict]:
    smiles = row.get("canonical_smiles", "").strip()
    if smiles in cache: return cache[smiles]
    pw = split_npc_field(row.get("npc_pathway", ""))
    sc = split_npc_field(row.get("npc_superclass", ""))
    cl = split_npc_field(row.get("npc_class", ""))
    if pw or sc or cl: return {"pathway": pw, "superclass": sc, "class": cl}
    return None

def get_raw_superclasses(npc_entry: Optional[dict]) -> List[str]:
    """Return raw superclass strings from NPC entry (not lowercased yet)."""
    if not npc_entry:
        return []
    return list(npc_entry.get("superclass", []) or [])

def get_raw_pathways(npc_entry: Optional[dict]) -> List[str]:
    if not npc_entry:
        return []
    return list(npc_entry.get("pathway", []) or [])

def superclass_lookup_key(raw: str) -> str:
    """Normalize superclass name for SUPERCLASS_TO_LABELS lookup:
    strip parenthetical suffix, lowercase, strip whitespace."""
    return re.sub(r"\s*\([^)]*\)\s*$", "", raw).strip().lower()

def map_superclass_labels(npc_entry: Optional[dict],
                          mol: Optional[Chem.Mol] = None,
                          ) -> Tuple[Dict[str, Decision], List[str]]:
    """Apply SUPERCLASS_TO_LABELS mapping.
    Returns (label->Decision dict, list of review hints).
    Special rule: 'Fatty acyls' requires C>=6."""
    results: Dict[str, Decision] = {}
    review_hints: List[str] = []
    raw_scs = get_raw_superclasses(npc_entry)
    matched_any = False
    for raw_sc in raw_scs:
        key = superclass_lookup_key(raw_sc)
        if key in SUPERCLASS_TO_LABELS:
            label_list = SUPERCLASS_TO_LABELS[key]
            # Special rule: fatty acyls needs C>=6
            if key == "fatty acyls" and mol is not None:
                s = summarize_mol(mol)
                if s["carbons"] < 6:
                    review_hints.append(f"fatty acyls but C={s['carbons']}<6")
                    continue
            for lbl in label_list:
                results[lbl] = Decision(auto=True, reason=f"NPC superclass: {raw_sc}")
            if not label_list:
                # Explicitly mapped to empty → Other (no action needed)
                pass
            matched_any = True
        elif key in SUPERCLASS_REVIEW:
            review_hints.append(f"ambiguous superclass: {raw_sc}")
            matched_any = True
    return results, review_hints


def load_classyfire_index(csv_path: Path) -> Dict[str, dict]:
    index: Dict[str, dict] = {}
    if not csv_path.exists(): return index
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            smiles = (row.get("canonical_smiles") or "").strip()
            if smiles: index[smiles] = row
    return index

def map_classyfire_labels(cf_row: Optional[dict]) -> Set[str]:
    if not cf_row: return set()
    sup = normalize_cf(cf_row.get("cf_superclass", ""))
    cl = normalize_cf(cf_row.get("cf_class", ""))
    sub = normalize_cf(cf_row.get("cf_subclass", ""))
    dp = normalize_cf(cf_row.get("cf_direct_parent", ""))
    af = (sup, cl, sub, dp)
    labels: Set[str] = set()
    if any("steroid" in x for x in af): labels.add("Steroid")
    if any(tok in x for x in af for tok in CF_TERPENE_TOKENS):
        labels.add("Terpenoid")
    if "alkaloid" in sup or "alkaloid" in cl or "alkaloid" in dp: labels.add("Alkaloid")
    if "amino acid" in sup or "amino acid" in cl or "peptide" in cl or "peptide" in sup: labels.add("Amino_acid")
    if "fatty" in sup or "fatty" in cl: labels.add("Fatty_acid")
    if cl in CF_PHENYLPROPANOID_CLASSES: labels.add("Phenylpropanoid")
    if cl in CF_POLYKETIDE_CLASSES: labels.add("Polyketide")
    if "phenylpropan" in sup: labels.add("Phenylpropanoid")
    if "polyketide" in sup: labels.add("Polyketide")
    return labels

def safe_mol(smiles: str) -> Optional[Chem.Mol]:
    if not smiles or "*" in smiles: return None
    return Chem.MolFromSmiles(smiles)

def has_any(mol: Chem.Mol, patterns: Iterable[Chem.Mol]) -> bool:
    return any(mol.HasSubstructMatch(p) for p in patterns)

def count_matches(mol: Chem.Mol, pattern: Chem.Mol) -> int:
    return len(mol.GetSubstructMatches(pattern))

def has_nitrogen(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() == 7 for a in mol.GetAtoms())

def has_ring_nitrogen(mol: Chem.Mol) -> bool:
    return any(a.GetAtomicNum() == 7 and a.IsInRing() for a in mol.GetAtoms())

def has_fluorine(mol):
    return mol is not None and mol.HasSubstructMatch(PAT_FLUORINE)

def has_valid_alpha_amino_acid_backbone(mol):
    """Reject oximes and require sp3 alpha carbon with single bonds."""
    if mol.HasSubstructMatch(PAT_OXIME):
        return False
    for pat in (PAT_ALPHA_AA, PAT_ALPHA_AA_AR):
        for match in mol.GetSubstructMatches(pat):
            if len(match) < 3: continue
            n_idx, alpha_idx, carbonyl_idx = match[:3]
            n_atom = mol.GetAtomWithIdx(n_idx)
            alpha_atom = mol.GetAtomWithIdx(alpha_idx)
            bond_na = mol.GetBondBetweenAtoms(n_idx, alpha_idx)
            bond_ac = mol.GetBondBetweenAtoms(alpha_idx, carbonyl_idx)
            if bond_na is None or bond_ac is None: continue
            if bond_na.GetBondType() != Chem.BondType.SINGLE: continue
            if bond_ac.GetBondType() != Chem.BondType.SINGLE: continue
            if alpha_atom.GetHybridization() != Chem.HybridizationType.SP3: continue
            if n_atom.GetIsAromatic(): continue
            return True
    return False

def summarize_mol(mol: Optional[Chem.Mol]) -> Dict[str, int]:
    if mol is None:
        return {"carbons":0,"nitrogens":0,"oxygens":0,"halogens":0,"heavy_atoms":0,"rings":0,"aromatic_rings":0}
    ri = mol.GetRingInfo()
    return {
        "carbons": sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6),
        "nitrogens": sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7),
        "oxygens": sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 8),
        "halogens": sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() in (9,17,35,53)),
        "heavy_atoms": mol.GetNumHeavyAtoms(),
        "rings": ri.NumRings(),
        "aromatic_rings": sum(1 for ring in ri.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)),
    }

def get_carbon_flags(mol: Optional[Chem.Mol]) -> Dict[str, bool]:
    s = summarize_mol(mol)
    return {"c18_c29": 18 <= s["carbons"] <= 29, "contains_oxygen": s["oxygens"] > 0,
            "contains_nitrogen": s["nitrogens"] > 0, "aromatic_present": s["aromatic_rings"] > 0,
            "terpenoid_plausible": s["carbons"] >= 10}

def is_probably_synthetic_xenobiotic(mol: Optional[Chem.Mol]) -> bool:
    if mol is None: return False
    score = 0
    if mol.HasSubstructMatch(PAT_CF3): score += 2
    if mol.HasSubstructMatch(PAT_SULFONE) or mol.HasSubstructMatch(PAT_SULFONAMIDE): score += 2
    if mol.HasSubstructMatch(PAT_AZO): score += 2
    if mol.HasSubstructMatch(PAT_TRIAZOLE_1) or mol.HasSubstructMatch(PAT_TRIAZOLE_2): score += 2
    s = summarize_mol(mol)
    if s["halogens"] >= 2: score += 1
    if any(a.GetAtomicNum() == 9 for a in mol.GetAtoms()): score += 1
    return score >= 2

def longest_aliphatic_chain(mol: Chem.Mol) -> int:
    """Find longest contiguous chain of sp3 carbons (not in ring)."""
    best = 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 6 or a.IsInRing() or a.GetIsAromatic():
            continue
        stack = [(a.GetIdx(), 1, {a.GetIdx()})]
        while stack:
            idx, depth, seen = stack.pop()
            best = max(best, depth)
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                nid = nbr.GetIdx()
                if nid in seen or nbr.GetAtomicNum() != 6 or nbr.IsInRing() or nbr.GetIsAromatic():
                    continue
                stack.append((nid, depth + 1, seen | {nid}))
    return best

def longest_acyl_chain(mol: Chem.Mol, carbonyl_pattern: Chem.Mol) -> int:
    best = 0
    for match in mol.GetSubstructMatches(carbonyl_pattern):
        ci = match[0]; ca = mol.GetAtomWithIdx(ci); blocked = set(match)
        for nbr in ca.GetNeighbors():
            if nbr.GetIdx() in blocked or nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic(): continue
            stack = [(nbr.GetIdx(), 1, blocked | {nbr.GetIdx()})]
            while stack:
                idx, depth, seen = stack.pop(); best = max(best, depth)
                for nxt in mol.GetAtomWithIdx(idx).GetNeighbors():
                    nid = nxt.GetIdx()
                    if nid in seen or nxt.GetAtomicNum() != 6 or nxt.GetIsAromatic(): continue
                    if nxt.IsInRing() and depth > 1: continue
                    stack.append((nid, depth+1, seen|{nid}))
    return best

def has_macrolactone(mol: Chem.Mol, min_ring_size: int = 12) -> bool:
    em = mol.GetSubstructMatches(PAT_ESTER)
    if not em: return False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) < min_ring_size: continue
        rs = set(ring)
        for m in em:
            if m[0] in rs and m[2] in rs: return True
    return False

def _aryl_side_chain_length(mol: Chem.Mol) -> int:
    """Max carbon chain length extending from any aromatic ring.
    Returns -1 if no aromatic ring found."""
    rings = mol.GetRingInfo().AtomRings()
    ar_atoms: Set[int] = set()
    for r in rings:
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r):
            ar_atoms.update(r)
    if not ar_atoms:
        return -1
    best = 0
    for idx in ar_atoms:
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            nid = nbr.GetIdx()
            if nid in ar_atoms or nbr.GetAtomicNum() in (8, 7, 16):
                continue
            # BFS counting carbon chain
            visited = ar_atoms | {nid}
            stack = [(nid, 1)]
            while stack:
                curr, depth = stack.pop()
                if mol.GetAtomWithIdx(curr).GetAtomicNum() == 6:
                    best = max(best, depth)
                for n2 in mol.GetAtomWithIdx(curr).GetNeighbors():
                    if n2.GetIdx() not in visited and n2.GetAtomicNum() == 6:
                        visited.add(n2.GetIdx())
                        stack.append((n2.GetIdx(), depth + 1))
    return best


def has_steroid_topology(mol: Chem.Mol) -> bool:
    """Check for 6-6-6-5 fused ring system with ALL-CARBON core."""
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    six = [r for r in rings if len(r)==6]; five = [r for r in rings if len(r)==5]
    if len(six)<3 or len(five)<1: return False
    for i in range(len(six)):
        fused = [six[j] for j in range(len(six)) if j!=i and len(six[i]&six[j])>=2]
        if len(fused)<2: continue
        core = six[i]|fused[0]|fused[1]
        for fr in five:
            if len(fr&core)>=2:
                full_core = core | fr
                has_N_in_core = any(
                    mol.GetAtomWithIdx(idx).GetAtomicNum() == 7
                    for idx in full_core
                )
                if has_N_in_core:
                    continue
                aromatic_in_core = sum(
                    1 for idx in full_core
                    if mol.GetAtomWithIdx(idx).GetIsAromatic()
                )
                if aromatic_in_core >= 4:
                    continue
                return True
    return False


# --- Detection functions (each class independent, no cross-class suppression) ---

def detect_steroid(mol, npc, cf_labels, cflags):
    if npc & NPC_STEROID:
        return Decision(auto=True, reason="NPC explicit steroid")
    if mol is not None and has_steroid_topology(mol):
        if cflags["c18_c29"]:
            return Decision(auto=True, reason="6-6-6-5 topology + C18-C29")
        return Decision(review=True, reason="steroid-like topology but outside C18-C29")
    if "Steroid" in cf_labels:
        return Decision(review=True, reason="ClassyFire steroid only")
    return Decision()

def detect_terpenoid(mol, npc, cf_labels):
    ta = sorted(npc & NPC_TERPENOID_AUTO)
    if ta:
        s = summarize_mol(mol)
        if s["carbons"] < 10:
            return Decision(review=True, reason=f"NPC terpenoid {ta[0]} but C<10")
        return Decision(auto=True, reason=f"NPC terpenoid superclass/class: {ta[0]}")
    if npc & NPC_TERPENOID_REVIEW:
        return Decision(review=True, reason="NPC pathway-level terpenoid only")
    if "Terpenoid" in cf_labels:
        return Decision(review=True, reason="ClassyFire terpenoid only")
    return Decision()

def detect_alkaloid(mol, npc, cf_labels):
    if mol is None:
        if npc & NPC_ALKALOID or "Alkaloid" in cf_labels:
            return Decision(review=True, reason="alkaloid call without structure")
        return Decision()
    if not has_nitrogen(mol): return Decision()
    syn = is_probably_synthetic_xenobiotic(mol)
    hp = has_any(mol, HIGH_PRECISION_ALKALOID_SCAFFOLDS)
    # Broadened: bare indole or phenylethylamine is enough for auto
    broad = has_any(mol, BROAD_ALKALOID_SCAFFOLDS)
    gs = has_any(mol, GENERIC_ALKALOID_SCAFFOLDS) or has_ring_nitrogen(mol)

    if hp and not syn:
        return Decision(auto=True, reason="high-precision natural-product alkaloid scaffold")
    if hp and syn:
        return Decision(review=True, reason="high-precision alkaloid scaffold but synthetic-like")
    if broad and not syn:
        return Decision(auto=True, reason="broad alkaloid scaffold (indole/phenylethylamine)")
    if broad and syn:
        return Decision(review=True, reason="broad alkaloid scaffold but synthetic-like")
    if gs:
        if npc & NPC_ALKALOID:
            return Decision(auto=True, reason="generic N-heterocycle + NPC alkaloid")
        if "Alkaloid" in cf_labels:
            return Decision(review=True, reason="generic N-heterocycle + ClassyFire alkaloid")
        return Decision(review=True, reason="generic N-heterocycle only")
    if npc & NPC_ALKALOID:
        return Decision(review=True, reason="NPC alkaloid without scaffold match")
    if "Alkaloid" in cf_labels:
        return Decision(review=True, reason="ClassyFire alkaloid without scaffold")
    return Decision()

def detect_amino_acid(mol, npc, cf_labels):
    if mol is None:
        if npc & NPC_AMINO or "Amino_acid" in cf_labels:
            return Decision(review=True, reason="amino-acid evidence without structure")
        return Decision()
    ha = has_valid_alpha_amino_acid_backbone(mol)
    hd = mol.HasSubstructMatch(PAT_DKP)
    pb = count_matches(mol, PAT_PEPTIDE_BOND)
    ad = mol.HasSubstructMatch(PAT_ADENOSINE)
    if ha:
        if ad: return Decision(review=True, reason="alpha-AA backbone but adenosine/cofactor-like")
        return Decision(auto=True, reason="intact alpha-amino-acid backbone")
    if hd:
        if has_fluorine(mol) or is_probably_synthetic_xenobiotic(mol):
            return Decision(review=True, reason="diketopiperazine but fluorinated/synthetic-like")
        return Decision(auto=True, reason="diketopiperazine")
    if pb >= 2 and (npc & NPC_AMINO): return Decision(auto=True, reason=f"peptide ({pb} amide bonds) + NPC")
    if pb >= 2 and "Amino_acid" in cf_labels: return Decision(review=True, reason=f"peptide ({pb} amide bonds) + CF")
    if pb >= 1 and (npc & NPC_AMINO or "Amino_acid" in cf_labels):
        return Decision(review=True, reason="partial peptide/amino-acid evidence")
    return Decision()

def detect_phenylpropanoid(mol, npc, cf_labels):
    if mol is None:
        if npc & NPC_PHENYLPROPANOID_REVIEW or "Phenylpropanoid" in cf_labels:
            return Decision(review=True, reason="phenylpropanoid evidence without structure")
        return Decision()
    syn = is_probably_synthetic_xenobiotic(mol)
    if has_any(mol, TRUE_PHENYLPROPANOID_PATTERNS):
        if syn: return Decision(review=True, reason="true PP scaffold but synthetic-like")
        return Decision(auto=True, reason="true phenylpropanoid scaffold")
    # C6-C1 phenolic acids (benzoic acid derivatives) are phenylpropanoid-derived
    if has_any(mol, C6C1_PATTERNS):
        if syn: return Decision(review=True, reason="C6-C1 phenolic acid but synthetic-like")
        return Decision(auto=True, reason="C6-C1 phenolic acid (benzoic acid derivative)")
    if mol.HasSubstructMatch(PAT_HYDROCINNAMATE):
        return Decision(auto=True, reason="hydrocinnamate pattern")
    na = sorted(npc & NPC_PHENYLPROPANOID_AUTO)
    if na:
        if syn: return Decision(review=True, reason=f"NPC PP class {na[0]} but synthetic-like")
        return Decision(auto=True, reason=f"NPC specific phenylpropanoid class: {na[0]}")
    if npc & NPC_PHENYLPROPANOID_REVIEW:
        return Decision(review=True, reason="NPC broad phenylpropanoid/pathway")
    if "Phenylpropanoid" in cf_labels:
        return Decision(review=True, reason="ClassyFire phenylpropanoid only")
    return Decision()

def detect_polyketide(mol, npc, cf_labels):
    if mol is None:
        if npc & NPC_POLYKETIDE_REVIEW or "Polyketide" in cf_labels:
            return Decision(review=True, reason="polyketide evidence without structure")
        return Decision()
    tri = mol.HasSubstructMatch(PAT_TRIAZOLE_1) or mol.HasSubstructMatch(PAT_TRIAZOLE_2)
    if has_any(mol, POLYKETIDE_STRUCTURES):
        return Decision(auto=True, reason="polyketide quinone/xanthone scaffold")
    na = sorted(npc & NPC_POLYKETIDE_AUTO)
    if na:
        if tri and has_macrolactone(mol, 12):
            return Decision(review=True, reason="NPC PK class but synthetic triazole macrolactone risk")
        return Decision(auto=True, reason=f"NPC specific polyketide class: {na[0]}")
    if has_macrolactone(mol, 12):
        syn = is_probably_synthetic_xenobiotic(mol)
        if tri or syn: return Decision(review=True, reason="macrolactone but synthetic-like")
        return Decision(review=True, reason="macrolactone only")
    if npc & NPC_POLYKETIDE_REVIEW:
        return Decision(review=True, reason="NPC pathway-level polyketide")
    if "Polyketide" in cf_labels:
        return Decision(review=True, reason="ClassyFire polyketide only")
    return Decision()

def detect_fatty_acid(mol, npc, cf_labels):
    # Broadened: include alkanes/alcohols/ketones/aldehydes with chain>=6
    if npc & NPC_JASMONOID:
        return Decision(auto=True, reason="jasmonoid (fatty-acid derived)")
    if mol is None:
        if npc & NPC_FATTY_REVIEW or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence without structure")
        return Decision()
    s = summarize_mol(mol)
    fs = bool(npc & NPC_FATTY_SPECIFIC)

    # Long aliphatic chain (>=6C) with no rings -> fatty even without carbonyl
    chain_len = longest_aliphatic_chain(mol)
    if chain_len >= 6 and s["rings"] == 0 and s["aromatic_rings"] == 0:
        return Decision(auto=True, reason=f"aliphatic chain={chain_len} no rings (broadened fatty)")

    # Acyl-group based detection
    ag = mol.HasSubstructMatch(PAT_CARBONYL_OXYGENATED)
    if not ag:
        if fs or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence but no acyl group")
        return Decision()
    bc, bn = 0, ""
    for pat, name in ((PAT_FREE_ACID,"free acid"),(PAT_ESTER,"ester"),(PAT_AMIDE,"amide"),(PAT_THIOESTER,"thioester")):
        c = longest_acyl_chain(mol, pat)
        if c > bc: bc, bn = c, name
    # Also check ketone/aldehyde chains
    kc = longest_acyl_chain(mol, PAT_KETONE)
    if kc > bc: bc, bn = kc, "ketone"

    if bc < 4:
        if npc & NPC_FATTY_REVIEW or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence but chain <4C")
        return Decision()
    if bn == "free acid" and bc >= 8 and s["rings"] == 0:
        return Decision(auto=True, reason=f"strong free fatty acid chain={bc} no rings")
    if s["aromatic_rings"] > 0 and not fs:
        return Decision(review=True, reason=f"fatty {bn} chain={bc} but aromatic")
    return Decision(auto=True, reason=f"fatty {bn} chain={bc}")


# --- Main ---

def classify_compound(smiles, npc_entry, cf_row=None, gold_labels=None):
    """Classify a compound using a prioritised pipeline:
    P1: Gold standard match → P2: NPC superclass mapping →
    P3: NPC pathway fallback (review) → P4: SMARTS structural detection →
    P5: Other.
    Each label is determined independently (no cross-class suppression).
    """
    mol = safe_mol(smiles)
    npc = normalize_npc(npc_entry)
    cfl = map_classyfire_labels(cf_row)
    cf = get_carbon_flags(mol)

    # Initialise all labels as empty decisions
    d: Dict[str, Decision] = {l: Decision() for l in LABELS}

    # --- Priority 1: Gold standard override ---
    has_gold = False
    if gold_labels:
        for label in gold_labels:
            if label in d:
                d[label] = Decision(auto=True, reason=f"gold: NPClassifier training data ({label})")
                has_gold = True

    # --- Priority 2: NPC Superclass mapping (PRIMARY) ---
    sc_decisions, sc_review_hints = map_superclass_labels(npc_entry, mol)
    has_superclass = bool(sc_decisions) or bool(sc_review_hints)
    for lbl, dec in sc_decisions.items():
        if lbl in d and not d[lbl].auto:
            d[lbl] = dec
    # Record review hints from ambiguous superclasses
    sc_review_labels: List[str] = []
    for hint in sc_review_hints:
        for lbl in LABELS:
            if not d[lbl].auto and not d[lbl].review:
                d[lbl] = Decision(review=True, reason=hint)
                sc_review_labels.append(lbl)
                break  # one hint → one review slot

    # --- Priority 3: NPC Pathway fallback (review only) ---
    if not has_gold and not has_superclass:
        raw_pws = get_raw_pathways(npc_entry)
        if raw_pws and not any(d[l].auto for l in LABELS):
            pw_hint = "|".join(raw_pws)
            # Map known pathways to review labels
            PW_HINT_MAP = {
                "terpenoids": "Terpenoid",
                "alkaloids": "Alkaloid",
                "amino acids and peptides": "Amino_acid",
                "fatty acids": "Fatty_acid",
                "shikimates and phenylpropanoids": "Phenylpropanoid",
                "polyketides": "Polyketide",
            }
            pw_mapped = False
            for pw_raw in raw_pws:
                pw_key = superclass_lookup_key(pw_raw)
                if pw_key in PW_HINT_MAP:
                    lbl = PW_HINT_MAP[pw_key]
                    if not d[lbl].auto and not d[lbl].review:
                        d[lbl] = Decision(review=True, reason=f"NPC pathway hint: {pw_raw}")
                        pw_mapped = True
            if not pw_mapped:
                # Unknown pathway → general review
                for lbl in LABELS:
                    if not d[lbl].auto and not d[lbl].review:
                        d[lbl] = Decision(review=True, reason=f"NPC pathway (unmapped): {pw_hint}")
                        break

    # --- Priority 4: Amino_acid SMARTS only (other P4 detectors removed) ---
    # Only Amino_acid structural detection is reliable (100% in spot-check).
    # Steroid/Alkaloid/Fatty_acid/Phenylpropanoid/Polyketide SMARTS had
    # 75-100% error rates and are removed.
    if not d["Amino_acid"].auto:
        det = detect_amino_acid(mol, npc, cfl)
        if det.auto or (det.review and not d["Amino_acid"].review):
            d["Amino_acid"] = det

    # --- Priority 4b: NPC alkaloid corrections ---
    # NPC labels substituted tryptophans as "Tryptophan alkaloid", but if the
    # alpha-amino acid backbone (NH2-CH-COOH) is still intact, the compound
    # is at the amino acid stage of biosynthesis, not yet an alkaloid.
    # Principle: classify by the actual skeleton at the time P450 acts on it.
    if d["Alkaloid"].auto and d["Amino_acid"].auto:
        if mol is not None and has_valid_alpha_amino_acid_backbone(mol):
            d["Alkaloid"] = Decision(
                review=True,
                reason="NPC alkaloid but amino acid backbone intact → Amino_acid primary"
            )

    # 4b-ii: Oxime compounds (CYP79 products like IAOx) are amino acid-derived
    # intermediates, not yet alkaloids.  Downgrade Alkaloid to review.
    if d["Alkaloid"].auto and mol is not None and mol.HasSubstructMatch(PAT_OXIME):
        d["Alkaloid"] = Decision(
            review=True,
            reason="NPC alkaloid but has oxime group → CYP79 amino acid intermediate"
        )

    # 4b-iii: Very simple indole derivatives (e.g., indole-3-carboxylic acid,
    # indole-3-acetic acid) are too small to be true alkaloids.  If the molecule
    # has ≤12 heavy atoms and NPC says alkaloid, downgrade to review.
    if d["Alkaloid"].auto and mol is not None:
        s = summarize_mol(mol)
        if s["heavy_atoms"] <= 12 and s["rings"] <= 2:
            d["Alkaloid"] = Decision(
                review=True,
                reason=f"NPC alkaloid but molecule too simple (HA={s['heavy_atoms']}, rings={s['rings']})"
            )

    # --- Priority 4c: Naphthalenes polyketide correction (Codex-reviewed) ---
    # NPC "Naphthalenes" is a mixed superclass: some are true polyketides
    # (rubrofusarin, biflaviolin, fonsecin), some are terpenoid quinones
    # (tanshinones). Only downgrade Polyketide when Terpenoid co-label exists
    # (indicating a terpenoid quinone, not a true polyketide).
    if d["Polyketide"].auto and d["Terpenoid"].auto:
        reason_pk = d["Polyketide"].reason
        if "Naphthalenes" in reason_pk:
            d["Polyketide"] = Decision(
                review=True,
                reason="Naphthalenes + Terpenoid co-label → likely terpenoid quinone, not polyketide"
            )

    # --- Priority 4d: Fatty_acid carbon count gate ---
    # NPC "Fatty acyls/esters" includes very small molecules (C<6) that are
    # not true fatty acids. IUPAC definition covers C4-C28, but C<6 compounds
    # in our dataset are mostly non-fatty (propanol, itaconate, etc.).
    if d["Fatty_acid"].auto and mol is not None:
        s_fa = summarize_mol(mol)
        if s_fa["carbons"] < 6:
            d["Fatty_acid"] = Decision(
                review=True,
                reason=f"NPC fatty but total carbons={s_fa['carbons']}<6 → too short"
            )

    # --- Priority 4e: Phenylpropanoid C6-C3 verification ---
    # NPC "Phenylpropanoids (C6-C3)" sometimes includes C6-C2 or C6-C1
    # compounds. Verify aryl side chain length; if <3, downgrade to review.
    if d["Phenylpropanoid"].auto and mol is not None:
        if "NPC superclass: Phenylpropanoids" in d["Phenylpropanoid"].reason:
            chain_len = _aryl_side_chain_length(mol)
            if 0 <= chain_len < 3:
                d["Phenylpropanoid"] = Decision(
                    review=True,
                    reason=f"NPC says C6-C3 but measured aryl chain={chain_len} → review"
                )

    # --- Build output ---
    flags = {l: int(d[l].auto) for l in LABELS}
    rc = [l for l in LABELS if d[l].review and not d[l].auto]
    if has_gold:
        tier = "gold"
    elif any(flags.values()):
        tier = "auto"
    elif rc:
        tier = "review"
    else:
        tier = "other"
    other = 1 if tier == "other" else 0
    reasons = {l: d[l].reason for l in LABELS if d[l].reason}
    flags["Other"] = other
    return flags, rc, reasons, tier

def run(input_csv, cache_json, classyfire_csv, output_csv, gold_tsv=None):
    sys.stdout.reconfigure(encoding="utf-8")
    with cache_json.open("r", encoding="utf-8") as f: cache = json.load(f)
    print(f"NPC cache: {len(cache)} entries")
    ci = load_classyfire_index(classyfire_csv) if classyfire_csv.exists() else {}
    if ci: print(f"ClassyFire: {len(ci)} entries")
    gold = {}
    if gold_tsv and gold_tsv.exists():
        gold = load_gold_standard(gold_tsv)
        print(f"Gold standard: {len(gold)} unique SMILES")
    rows = []
    with input_csv.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f); of = list(reader.fieldnames or [])
        for row in reader: rows.append(row)
    print(f"Compounds: {len(rows)}")
    outf = of + list(LABELS) + ["Other","tier","review_candidates","multilabel_vector","multilabel_reason"]
    tc, lc, rlc = Counter(), Counter(), Counter()
    mc = 0
    gold_hits = 0
    with output_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=outf); w.writeheader()
        for row in rows:
            smi = row.get("canonical_smiles","").strip()
            ne = get_npc_entry(row, cache); cr = ci.get(smi)
            # Look up gold standard via canonical SMILES
            gl = None
            if gold and smi:
                mol_tmp = Chem.MolFromSmiles(smi)
                if mol_tmp:
                    csmi = Chem.MolToSmiles(mol_tmp)
                    gl = gold.get(csmi)
            flags, rc, reasons, tier = classify_compound(smi, ne, cr, gl)
            out = dict(row)
            for l in LABELS: out[l] = flags[l]
            out["Other"] = flags["Other"]; out["tier"] = tier
            out["review_candidates"] = json.dumps(rc, ensure_ascii=False)
            out["multilabel_vector"] = json.dumps([flags[l] for l in LABELS])
            parts = []
            for l in LABELS:
                if l in reasons:
                    pfx = "AUTO" if flags[l] else "REVIEW"
                    parts.append(f"{l}={pfx}: {reasons[l]}")
            out["multilabel_reason"] = "; ".join(parts)
            w.writerow(out)
            tc[tier] += 1
            if tier == "gold": gold_hits += 1
            pos = [l for l in LABELS if flags[l]]
            for l in pos: lc[l] += 1
            if len(pos) > 1: mc += 1
            for r in rc: rlc[r] += 1
    print(f"\n=== Results ({len(rows)} compounds) ===")
    print("\n--- Tiers ---")
    for t in ("gold","auto","review","other"):
        print(f"  {t:10s} {tc[t]:5d} ({tc[t]/len(rows)*100:.1f}%)")
    print("\n--- Auto Labels (includes gold) ---")
    for l in LABELS: print(f"  {l:20s} {lc[l]:5d}")
    print(f"  {'Multi-label':20s} {mc:5d}")
    if gold_hits:
        print(f"\n  Gold standard hits: {gold_hits}")
    print("\n--- Review Candidates ---")
    for l in LABELS:
        if rlc[l]: print(f"  {l:20s} {rlc[l]:5d}")
    print(f"\nOutput: {output_csv}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--classyfire", type=Path, default=DEFAULT_CLASSYFIRE)
    p.add_argument("--gold", type=Path, default=DEFAULT_GOLD)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    a = p.parse_args()
    run(a.input, a.cache, a.classyfire, a.output, a.gold)

if __name__ == "__main__":
    main()
