#!/usr/bin/env python3
"""
Assign strict multilabel natural-product classes to P450 substrates.

This version uses a 3-tier output:
  - auto:   high-precision automatic labels
  - review: evidence exists, but not strong enough for auto-labeling
  - other:  no evidence for any of the 7 classes

Output labels are AUTO labels only.

Literature sources: 7 independent agent searches of IUPAC, Dewick, LIPID MAPS, Pelletier, Vogt, Hertweck.
See sessions/05_底物分类与验证/literature_definitions/ for full reports.
Codex multi-round review (Rounds 1-7). 400-sample audit informed final criteria.
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


def _smarts(pattern: str) -> Chem.Mol:
    mol = Chem.MolFromSmarts(pattern)
    if mol is None:
        raise ValueError(f"Invalid SMARTS: {pattern}")
    return mol


# --- SMARTS patterns ---

PAT_ALPHA_AA = _smarts("[NX3,NX4+;!$(NC=O)][CH1,CH2][CX3](=O)[OX2H1,OX1-,OX2]")
PAT_ALPHA_AA_AR = _smarts("[NX3,NX4+;!$(NC=O)][CH]([#6])[CX3](=O)[OX2H1,OX1-,OX2]")
PAT_PEPTIDE_BOND = _smarts("[NX3,NX4+][CX3](=O)[#6]")
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

TRUE_PHENYLPROPANOID_PATTERNS = (
    PAT_CINNAMATE, PAT_CINNAMATE2, PAT_STILBENE,
    PAT_COUMARIN, PAT_COUMARIN2,
    PAT_FLAVONE_CORE, PAT_FLAVANONE, PAT_CHALCONE, PAT_ISOFLAVONE,
)

PAT_ANTHRAQUINONE = _smarts("O=C1c2ccccc2C(=O)c2ccccc21")
PAT_NAPHTHOQUINONE = _smarts("O=C1C=CC(=O)c2ccccc21")
PAT_XANTHONE = _smarts("O=c1oc2ccccc2c2ccccc12")
POLYKETIDE_STRUCTURES = (PAT_ANTHRAQUINONE, PAT_NAPHTHOQUINONE, PAT_XANTHONE)

PAT_INDOLE = _smarts("c1ccc2[nH]ccc2c1")
PAT_TETRAHYDROISOQUINOLINE = _smarts("c1ccc2CCNCC2c1")
PAT_ISOQUINOLINE = _smarts("c1cnc2ccccc2c1")
PAT_TROPANE = _smarts("C1CC2CCC1N2")
PAT_PURINE = _smarts("c1ncnc2[nH]cnc12")
PAT_PYRIDINE = _smarts("n1ccccc1")
PAT_IMIDAZOLE = _smarts("[nH]1ccnc1")
PAT_PIPERIDINE = _smarts("N1CCCCC1")
PAT_QUINOLINE = _smarts("c1ccc2ncccc2c1")

HIGH_PRECISION_ALKALOID_SCAFFOLDS = (
    PAT_INDOLE, PAT_TETRAHYDROISOQUINOLINE, PAT_ISOQUINOLINE,
    PAT_TROPANE, PAT_PURINE,
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


@dataclass
class Decision:
    auto: bool = False
    review: bool = False
    reason: str = ""


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
        if "Steroid" not in labels: labels.add("Terpenoid")
    if "alkaloid" in sup or "alkaloid" in cl or "alkaloid" in dp: labels.add("Alkaloid")
    if "amino acid" in sup or "amino acid" in cl or "peptide" in cl or "peptide" in sup: labels.add("Amino_acid")
    if "fatty" in sup or "fatty" in cl: labels.add("Fatty_acid")
    if cl in CF_PHENYLPROPANOID_CLASSES: labels.add("Phenylpropanoid")
    if cl in CF_POLYKETIDE_CLASSES: labels.add("Polyketide")
    if sup != "phenylpropanoids and polyketides":
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
            "contains_nitrogen": s["nitrogens"] > 0, "aromatic_present": s["aromatic_rings"] > 0}

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

def has_steroid_topology(mol: Chem.Mol) -> bool:
    """Check for 6-6-6-5 fused ring system with ALL-CARBON core (no N in ring system).
    Real steroid gonane skeleton is pure C/H. Alkaloids and polyketides with similar
    ring topology are excluded by checking for heteroatoms in the fused core."""
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
                # Steroid core must be all-carbon: reject if N found in ring system
                has_N_in_core = any(
                    mol.GetAtomWithIdx(idx).GetAtomicNum() == 7
                    for idx in full_core
                )
                if has_N_in_core:
                    continue  # Not a steroid — likely alkaloid
                # Reject if core has aromatic rings (real steroids are saturated/partially saturated)
                aromatic_in_core = sum(
                    1 for idx in full_core
                    if mol.GetAtomWithIdx(idx).GetIsAromatic()
                )
                if aromatic_in_core >= 4:
                    continue  # Likely anthraquinone/xanthone polyketide, not steroid
                return True
    return False


# --- Detection functions ---

def detect_steroid(mol, npc, cf_labels, cflags):
    if npc & NPC_STEROID:
        return Decision(auto=True, reason="NPC explicit steroid")
    if mol is not None and has_steroid_topology(mol):
        if cflags["c18_c29"] and not (npc & NPC_TERPENOID_AUTO):
            return Decision(auto=True, reason="6-6-6-5 topology + C18-C29 + no NPC terpenoid")
        return Decision(review=True, reason="steroid-like topology but outside strict auto rule")
    if "Steroid" in cf_labels:
        return Decision(review=True, reason="ClassyFire steroid only")
    return Decision()

def detect_terpenoid(mol, npc, cf_labels, steroid_auto):
    if steroid_auto: return Decision()
    ta = sorted(npc & NPC_TERPENOID_AUTO)
    if ta: return Decision(auto=True, reason=f"NPC terpenoid superclass/class: {ta[0]}")
    if npc & NPC_TERPENOID_REVIEW:
        return Decision(review=True, reason="NPC pathway-level terpenoid only")
    if "Terpenoid" in cf_labels:
        return Decision(review=True, reason="ClassyFire terpenoid only")
    return Decision()

def detect_alkaloid(mol, npc, cf_labels, amino_auto):
    if mol is None:
        if npc & NPC_ALKALOID or "Alkaloid" in cf_labels:
            return Decision(review=True, reason="alkaloid call without structure")
        return Decision()
    if not has_nitrogen(mol): return Decision()
    syn = is_probably_synthetic_xenobiotic(mol)
    hp = has_any(mol, HIGH_PRECISION_ALKALOID_SCAFFOLDS)
    gs = has_any(mol, GENERIC_ALKALOID_SCAFFOLDS) or has_ring_nitrogen(mol)
    if amino_auto:
        if hp or npc & NPC_ALKALOID or "Alkaloid" in cf_labels:
            return Decision(review=True, reason="alkaloid evidence but intact amino-acid backbone")
        return Decision()
    if hp and not syn: return Decision(auto=True, reason="high-precision natural-product alkaloid scaffold")
    if hp and syn: return Decision(review=True, reason="high-precision alkaloid scaffold but synthetic-like")
    if gs:
        if npc & NPC_ALKALOID: return Decision(review=True, reason="generic N-heterocycle + NPC alkaloid")
        if "Alkaloid" in cf_labels: return Decision(review=True, reason="generic N-heterocycle + ClassyFire alkaloid")
        return Decision(review=True, reason="generic N-heterocycle only")
    if npc & NPC_ALKALOID: return Decision(review=True, reason="NPC alkaloid without high-precision scaffold")
    if "Alkaloid" in cf_labels: return Decision(review=True, reason="ClassyFire alkaloid without scaffold")
    return Decision()

def detect_amino_acid(mol, npc, cf_labels):
    if mol is None:
        if npc & NPC_AMINO or "Amino_acid" in cf_labels:
            return Decision(review=True, reason="amino-acid evidence without structure")
        return Decision()
    ha = mol.HasSubstructMatch(PAT_ALPHA_AA) or mol.HasSubstructMatch(PAT_ALPHA_AA_AR)
    hd = mol.HasSubstructMatch(PAT_DKP)
    pb = count_matches(mol, PAT_PEPTIDE_BOND)
    ad = mol.HasSubstructMatch(PAT_ADENOSINE)
    if ha:
        if ad: return Decision(review=True, reason="alpha-AA backbone but adenosine/cofactor-like")
        return Decision(auto=True, reason="intact alpha-amino-acid backbone")
    if hd: return Decision(auto=True, reason="diketopiperazine")
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
    iso = mol.HasSubstructMatch(PAT_ISOCOUMARIN)
    if has_any(mol, TRUE_PHENYLPROPANOID_PATTERNS):
        if syn: return Decision(review=True, reason="true PP scaffold but synthetic-like")
        return Decision(auto=True, reason="true phenylpropanoid scaffold")
    if mol.HasSubstructMatch(PAT_HYDROCINNAMATE):
        return Decision(review=True, reason="hydrocinnamate pattern only")
    na = sorted(npc & NPC_PHENYLPROPANOID_AUTO)
    if na:
        if syn: return Decision(review=True, reason=f"NPC PP class {na[0]} but synthetic-like")
        if iso: return Decision(review=True, reason=f"NPC coumarin-like {na[0]} but isocoumarin-like")
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

def detect_fatty_acid(mol, npc, cf_labels, steroid_auto, terpenoid_auto, polyketide_auto):
    if npc & NPC_JASMONOID:
        return Decision(review=True, reason="jasmonoid/fatty-derived evidence")
    if mol is None:
        if npc & NPC_FATTY_REVIEW or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence without structure")
        return Decision()
    cf = get_carbon_flags(mol)
    if not cf["contains_oxygen"]: return Decision()
    fs = bool(npc & NPC_FATTY_SPECIFIC)
    ag = mol.HasSubstructMatch(PAT_CARBONYL_OXYGENATED)
    if not ag:
        if npc & NPC_FATTY_REVIEW or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence but no acyl group")
        return Decision()
    if mol.HasSubstructMatch(PAT_KETONE) and not (
        mol.HasSubstructMatch(PAT_FREE_ACID) or mol.HasSubstructMatch(PAT_ESTER)
        or mol.HasSubstructMatch(PAT_AMIDE) or mol.HasSubstructMatch(PAT_THIOESTER)):
        return Decision()
    bc, bn = 0, ""
    for pat, name in ((PAT_FREE_ACID,"free acid"),(PAT_ESTER,"ester"),(PAT_AMIDE,"amide"),(PAT_THIOESTER,"thioester")):
        c = longest_acyl_chain(mol, pat)
        if c > bc: bc, bn = c, name
    if bc < 4:
        if npc & NPC_FATTY_REVIEW or "Fatty_acid" in cf_labels:
            return Decision(review=True, reason="fatty evidence but chain <4C")
        return Decision()
    occ = steroid_auto or terpenoid_auto or polyketide_auto
    if occ and not fs:
        return Decision(review=True, reason=f"fatty {bn} chain={bc} but occupied by another class")
    if cf["aromatic_present"] and not fs:
        return Decision(review=True, reason=f"fatty {bn} chain={bc} but aromatic")
    return Decision(auto=True, reason=f"fatty {bn} chain={bc}")


# --- Main ---

def classify_compound(smiles, npc_entry, cf_row=None):
    mol = safe_mol(smiles)
    npc = normalize_npc(npc_entry)
    cfl = map_classyfire_labels(cf_row)
    cf = get_carbon_flags(mol)
    d = {}
    d["Steroid"] = detect_steroid(mol, npc, cfl, cf)
    d["Amino_acid"] = detect_amino_acid(mol, npc, cfl)
    d["Alkaloid"] = detect_alkaloid(mol, npc, cfl, amino_auto=d["Amino_acid"].auto)
    d["Phenylpropanoid"] = detect_phenylpropanoid(mol, npc, cfl)
    d["Polyketide"] = detect_polyketide(mol, npc, cfl)
    d["Terpenoid"] = detect_terpenoid(mol, npc, cfl, steroid_auto=d["Steroid"].auto)
    d["Fatty_acid"] = detect_fatty_acid(mol, npc, cfl,
        steroid_auto=d["Steroid"].auto, terpenoid_auto=d["Terpenoid"].auto, polyketide_auto=d["Polyketide"].auto)
    flags = {l: int(d[l].auto) for l in LABELS}
    rc = [l for l in LABELS if d[l].review and not d[l].auto]
    if any(flags.values()): tier, other = "auto", 0
    elif rc: tier, other = "review", 0
    else: tier, other = "other", 1
    reasons = {l: d[l].reason for l in LABELS if d[l].reason}
    flags["Other"] = other
    return flags, rc, reasons, tier

def run(input_csv, cache_json, classyfire_csv, output_csv):
    sys.stdout.reconfigure(encoding="utf-8")
    with cache_json.open("r", encoding="utf-8") as f: cache = json.load(f)
    print(f"NPC cache: {len(cache)} entries")
    ci = load_classyfire_index(classyfire_csv) if classyfire_csv.exists() else {}
    if ci: print(f"ClassyFire: {len(ci)} entries")
    rows = []
    with input_csv.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f); of = list(reader.fieldnames or [])
        for row in reader: rows.append(row)
    print(f"Compounds: {len(rows)}")
    outf = of + list(LABELS) + ["Other","tier","review_candidates","multilabel_vector","multilabel_reason"]
    tc, lc, rlc = Counter(), Counter(), Counter()
    mc = 0
    with output_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=outf); w.writeheader()
        for row in rows:
            smi = row.get("canonical_smiles","").strip()
            ne = get_npc_entry(row, cache); cr = ci.get(smi)
            flags, rc, reasons, tier = classify_compound(smi, ne, cr)
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
            pos = [l for l in LABELS if flags[l]]
            for l in pos: lc[l] += 1
            if len(pos) > 1: mc += 1
            for r in rc: rlc[r] += 1
    print(f"\n=== Results ({len(rows)} compounds) ===")
    print("\n--- Tiers ---")
    for t in ("auto","review","other"):
        print(f"  {t:10s} {tc[t]:5d} ({tc[t]/len(rows)*100:.1f}%)")
    print("\n--- Auto Labels ---")
    for l in LABELS: print(f"  {l:20s} {lc[l]:5d}")
    print(f"  {'Multi-label':20s} {mc:5d}")
    print("\n--- Review Candidates ---")
    for l in LABELS:
        if rlc[l]: print(f"  {l:20s} {rlc[l]:5d}")
    print(f"\nOutput: {output_csv}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    p.add_argument("--classyfire", type=Path, default=DEFAULT_CLASSYFIRE)
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    a = p.parse_args()
    run(a.input, a.cache, a.classyfire, a.output)

if __name__ == "__main__":
    main()
