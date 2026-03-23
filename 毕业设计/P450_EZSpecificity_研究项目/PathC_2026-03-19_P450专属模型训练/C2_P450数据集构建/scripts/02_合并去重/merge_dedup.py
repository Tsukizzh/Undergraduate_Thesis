#!/usr/bin/env python3
"""
Merge and deduplicate the five primary P450 enzyme-substrate sources.

Outputs (to data/combined/):
    - global_enzymes.csv
    - global_compounds.csv
    - global_interactions.csv
    - enzyme_xref.csv        (optional)
    - compound_xref.csv      (optional)
    - interaction_evidence.csv (optional)
    - merge_audit.csv        (optional)

Design:
    Enzymes:      primary key = UniProt ID; fallback = exact sequence SHA256
    Compounds:    primary key = Standard InChIKey; fallback = canonical SMILES
    Interactions: primary key = (global_enzyme_id, global_compound_id)

Deterministic: groups sorted before ID assignment; outputs in sorted order.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

try:
    from rdkit.Chem.SaltRemover import SaltRemover as _SaltRemover
    _SALT_REMOVER = _SaltRemover()
    HAS_SALT_REMOVER = True
except ImportError:
    HAS_SALT_REMOVER = False

SCRIPT_PATH = Path(__file__).resolve()
PROJECT_ROOT = SCRIPT_PATH.parents[2]
DEFAULT_INPUT = PROJECT_ROOT / "data" / "sources"
DEFAULT_OUTPUT = PROJECT_ROOT / "data" / "combined"

SOURCE_CONFIG = [
    ("Source_RCSB",        "S1_RCSB"),
    ("Source_ESIBank",     "S2_ESIBank"),
    ("Source_P450Rdb",     "S3_P450Rdb"),
    ("Source_PlantP450DB", "S8_PlantP450DB"),
    ("Source_PCPD",        "S9_PCPD"),
]
SOURCE_ORDER = {tag: i for i, (_, tag) in enumerate(SOURCE_CONFIG)}
QUALITY_ORDER = {"A": 0, "B": 1, "C": 2, "": 99}


# ---------------------------------------------------------------------------
# Data records
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class EnzymeRec:
    source_tag: str
    source_id: str
    uniprot_id: str
    p450_symbol: str
    species: str
    species_class: str
    ec_number: str
    sequence: str
    seq_hash: str
    raw: dict


@dataclass(frozen=True)
class CompoundRec:
    source_tag: str
    source_id: str
    raw_smiles: str
    name: str
    canonical_smiles: str
    inchi: str
    inchikey: str
    parse_status: str
    raw: dict


@dataclass(frozen=True)
class InteractionRec:
    source_tag: str
    source_id: str
    enzyme_id: str
    compound_id: str
    label: str
    quality_tier: str
    num_reactions: int
    has_pmid: int
    has_products: int
    has_pdb: int
    pdb_id: str
    extra: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def norm(val: Any) -> str:
    if val is None:
        return ""
    t = str(val).strip()
    return "" if t.lower() in {"", "na", "n/a", "none", "null"} else t


def norm_seq(seq: Any) -> str:
    t = norm(seq).upper()
    return "".join(t.split()) if t else ""


def sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest() if seq else ""


def to_int(val: Any, default: int = 0) -> int:
    t = norm(val)
    if not t:
        return default
    try:
        return int(float(t))
    except ValueError:
        return default


def to_bool(val: Any) -> int:
    return 1 if norm(val).lower() in {"1", "true", "t", "yes", "y"} else 0


def read_csv(path: Path, required: set[str]) -> list[dict]:
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path}")
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        return [dict(r) for r in reader]


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log(f"  wrote {len(rows)} rows -> {path.name}")


def best_val(pairs: list[tuple[str, str]]) -> str:
    """Pick most-frequent non-empty value, breaking ties by source priority."""
    filt = [(v, s) for v, s in pairs if norm(v)]
    if not filt:
        return ""
    counts = Counter(v for v, _ in filt)
    ranked = sorted(counts, key=lambda v: (
        -counts[v],
        min(SOURCE_ORDER.get(s, 99) for val, s in filt if val == v),
        v,
    ))
    return ranked[0]


def aliases(vals: Iterable[str]) -> str:
    return "|".join(sorted({norm(v) for v in vals if norm(v)}))


def standardize_smiles(smi: str, strip_salts: bool) -> tuple[str, str, str, str]:
    """Returns (canonical_smiles, inchi, inchikey, status)."""
    raw = norm(smi)
    if not raw:
        return "", "", "", "missing"
    if not HAS_RDKIT:
        return raw, "", "", "no_rdkit"
    try:
        mol = Chem.MolFromSmiles(raw)
    except Exception:
        mol = None
    if mol is None:
        return raw, "", "", "invalid"
    if strip_salts and HAS_SALT_REMOVER:
        mol = _SALT_REMOVER.StripMol(mol, dontRemoveEverything=True) or mol
    try:
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return raw, "", "", "canon_fail"
    try:
        inchi = MolToInchi(mol) or ""
        inchikey = InchiToInchiKey(inchi) if inchi else ""
    except Exception:
        inchi, inchikey = "", ""
    return canonical, inchi, inchikey or "", "ok"


# ---------------------------------------------------------------------------
# Phase 1: Load
# ---------------------------------------------------------------------------
ENZ_COLS = {"enzyme_id", "uniprot_id", "p450_symbol", "species",
            "species_class", "ec_number", "sequence", "sequence_length"}
CMP_COLS = {"compound_id", "smiles", "name"}
INT_COLS = {"interaction_id", "enzyme_id", "compound_id", "label",
            "source", "quality_tier", "num_reactions", "has_pmid", "has_products"}


def load_enzymes(root: Path) -> list[EnzymeRec]:
    out = []
    for d, tag in SOURCE_CONFIG:
        for r in read_csv(root / d / "enzymes.csv", ENZ_COLS):
            seq = norm_seq(r.get("sequence", ""))
            out.append(EnzymeRec(
                source_tag=tag, source_id=norm(r["enzyme_id"]),
                uniprot_id=norm(r["uniprot_id"]).upper(),
                p450_symbol=norm(r.get("p450_symbol", "")),
                species=norm(r.get("species", "")),
                species_class=norm(r.get("species_class", "")),
                ec_number=norm(r.get("ec_number", "")),
                sequence=seq, seq_hash=sha256(seq), raw=r,
            ))
    return out


def load_compounds(root: Path, strip_salts: bool) -> list[CompoundRec]:
    out = []
    for d, tag in SOURCE_CONFIG:
        for r in read_csv(root / d / "compounds.csv", CMP_COLS):
            can, inchi, ikey, st = standardize_smiles(r.get("smiles", ""), strip_salts)
            out.append(CompoundRec(
                source_tag=tag, source_id=norm(r["compound_id"]),
                raw_smiles=norm(r.get("smiles", "")), name=norm(r.get("name", "")),
                canonical_smiles=can, inchi=inchi, inchikey=ikey,
                parse_status=st, raw=r,
            ))
    return out


def load_interactions(root: Path) -> tuple[list[InteractionRec], Counter]:
    out = []
    labels = Counter()
    for d, tag in SOURCE_CONFIG:
        for r in read_csv(root / d / "interactions.csv", INT_COLS):
            lab = norm(r.get("label", ""))
            labels[lab] += 1
            pdb = norm(r.get("pdb_id", ""))
            extra_keys = sorted(k for k in r if k not in INT_COLS and k != "pdb_id")
            extra = json.dumps({k: r.get(k, "") for k in extra_keys},
                               ensure_ascii=False, sort_keys=True)
            out.append(InteractionRec(
                source_tag=tag, source_id=norm(r["interaction_id"]),
                enzyme_id=norm(r["enzyme_id"]), compound_id=norm(r["compound_id"]),
                label=lab, quality_tier=norm(r.get("quality_tier", "")),
                num_reactions=to_int(r.get("num_reactions", 0)),
                has_pmid=to_bool(r.get("has_pmid", 0)),
                has_products=to_bool(r.get("has_products", 0)),
                has_pdb=int(bool(pdb)), pdb_id=pdb, extra=extra,
            ))
    return out, labels


# ---------------------------------------------------------------------------
# Phase 2: Enzyme resolution
# ---------------------------------------------------------------------------
def resolve_enzymes(recs: list[EnzymeRec]):
    log(f"[enzymes] {len(recs)} raw rows")
    by_uid: dict[str, list[EnzymeRec]] = defaultdict(list)
    no_uid: list[EnzymeRec] = []
    seq2uid: dict[str, set[str]] = defaultdict(set)
    audit = []

    for r in recs:
        if r.uniprot_id:
            by_uid[r.uniprot_id].append(r)
            if r.seq_hash:
                seq2uid[r.seq_hash].add(r.uniprot_id)
        else:
            no_uid.append(r)

    # Detect same-sequence → multiple UniProts
    for sh, uids in sorted(seq2uid.items()):
        if len(uids) > 1:
            audit.append(dict(entity_type="enzyme", global_id="", source="",
                              source_local_id="", merge_basis="seq_hash_conflict",
                              conflict_flags="same_seq_multi_uniprot",
                              notes=f"{sh[:16]}... -> {sorted(uids)}"))

    groups = []
    uid_index = {}

    # UniProt-based groups
    for uid in sorted(by_uid):
        members = sorted(by_uid[uid], key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id))
        g = dict(members=members, merge_basis="uniprot_id", uid=uid)
        groups.append(g)
        uid_index[uid] = g

    # Bridge no-UniProt records via exact sequence match
    remaining_by_seq: dict[str, list[EnzymeRec]] = defaultdict(list)
    orphans = []

    for r in no_uid:
        if r.seq_hash and r.seq_hash in seq2uid:
            candidates = sorted(seq2uid[r.seq_hash])
            if len(candidates) == 1:
                uid_index[candidates[0]]["members"].append(r)
                uid_index[candidates[0]]["merge_basis"] = "uniprot_id|seq_bridge"
                audit.append(dict(entity_type="enzyme", global_id="", source=r.source_tag,
                                  source_local_id=r.source_id, merge_basis="seq_bridge",
                                  conflict_flags="", notes=f"bridged to {candidates[0]}"))
            else:
                orphans.append(r)
                audit.append(dict(entity_type="enzyme", global_id="", source=r.source_tag,
                                  source_local_id=r.source_id, merge_basis="ambiguous_bridge",
                                  conflict_flags="multi_uniprot_match",
                                  notes=f"candidates={candidates}"))
        elif r.seq_hash:
            remaining_by_seq[r.seq_hash].append(r)
        else:
            orphans.append(r)

    # Sequence-hash groups (no UniProt anywhere)
    for sh in sorted(remaining_by_seq):
        members = sorted(remaining_by_seq[sh], key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id))
        groups.append(dict(members=members, merge_basis="sequence_hash", uid=""))

    # Orphans (no sequence, no UniProt)
    for r in sorted(orphans, key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id)):
        groups.append(dict(members=[r], merge_basis="singleton", uid=""))

    # Sort groups deterministically for ID assignment
    def sort_key(g):
        rank = 0 if g["uid"] else (1 if g["merge_basis"] == "sequence_hash" else 2)
        return (rank, g.get("uid", ""), min(SOURCE_ORDER[m.source_tag] for m in g["members"]))
    groups.sort(key=sort_key)

    # Assign global IDs and build outputs
    g_rows, xref_rows = [], []
    src2gid: dict[tuple[str, str], str] = {}

    for i, g in enumerate(groups, 1):
        gid = f"ENZ_G{i:06d}"
        members = sorted(g["members"], key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id))
        seqs = sorted({r.sequence for r in members if r.sequence})
        canon_seq = max(seqs, key=len) if seqs else ""
        canon_uid = g["uid"] or best_val([(r.uniprot_id, r.source_tag) for r in members])

        g_rows.append({
            "global_enzyme_id": gid,
            "canonical_uniprot_id": canon_uid,
            "canonical_p450_symbol": best_val([(r.p450_symbol, r.source_tag) for r in members]),
            "canonical_species": best_val([(r.species, r.source_tag) for r in members]),
            "canonical_species_class": best_val([(r.species_class, r.source_tag) for r in members]),
            "canonical_ec_number": best_val([(r.ec_number, r.source_tag) for r in members]),
            "canonical_sequence": canon_seq,
            "sequence_hash": sha256(canon_seq),
            "sequence_length": len(canon_seq) if canon_seq else "",
            "sources": "|".join(sorted({r.source_tag for r in members}, key=lambda t: SOURCE_ORDER[t])),
            "source_count": len({r.source_tag for r in members}),
            "member_count": len(members),
            "merge_basis": g["merge_basis"],
            "sequence_conflict": int(len(seqs) > 1),
            "species_conflict": int(len({r.species for r in members if r.species}) > 1),
            "ec_conflict": int(len({r.ec_number for r in members if r.ec_number}) > 1),
            "symbol_aliases": aliases(r.p450_symbol for r in members),
            "species_aliases": aliases(r.species for r in members),
            "ec_number_all": aliases(r.ec_number for r in members),
        })

        if len(seqs) > 1:
            audit.append(dict(entity_type="enzyme", global_id=gid, source="",
                              source_local_id="", merge_basis=g["merge_basis"],
                              conflict_flags="seq_conflict",
                              notes=f"n={len(seqs)} lens={[len(s) for s in seqs]}"))

        for r in members:
            src2gid[(r.source_tag, r.source_id)] = gid
            xref_rows.append({
                "global_enzyme_id": gid, "source": r.source_tag,
                "source_enzyme_id": r.source_id, "source_uniprot_id": r.uniprot_id,
                "source_p450_symbol": r.p450_symbol, "source_species": r.species,
                "source_ec_number": r.ec_number, "source_sequence_hash": r.seq_hash,
            })

    log(f"[enzymes] {len(recs)} -> {len(g_rows)} global enzymes "
        f"(uid_groups={len(by_uid)}, no_uid={len(no_uid)}, bridged={len(no_uid)-len(remaining_by_seq)-len(orphans)})")
    return g_rows, xref_rows, audit, src2gid


# ---------------------------------------------------------------------------
# Phase 3: Compound resolution
# ---------------------------------------------------------------------------
def resolve_compounds(recs: list[CompoundRec]):
    log(f"[compounds] {len(recs)} raw rows")
    grouped: dict[tuple[str, str], list[CompoundRec]] = defaultdict(list)
    audit = []

    for r in recs:
        if r.inchikey:
            key = ("inchikey", r.inchikey)
        elif r.canonical_smiles:
            key = ("canonical_smiles", r.canonical_smiles)
        elif r.raw_smiles:
            key = ("raw_smiles", r.raw_smiles)
        else:
            key = ("singleton", f"{r.source_tag}:{r.source_id}")
        grouped[key].append(r)

        if r.parse_status not in {"ok", "missing", "no_rdkit"}:
            audit.append(dict(entity_type="compound", global_id="", source=r.source_tag,
                              source_local_id=r.source_id, merge_basis="parse_issue",
                              conflict_flags=r.parse_status, notes=r.raw_smiles[:100]))

    g_rows, xref_rows = [], []
    src2gid: dict[tuple[str, str], str] = {}

    for i, gk in enumerate(sorted(grouped), 1):
        gid = f"CMP_G{i:06d}"
        members = sorted(grouped[gk], key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id))

        g_rows.append({
            "global_compound_id": gid,
            "canonical_name": best_val([(r.name, r.source_tag) for r in members]),
            "canonical_smiles": best_val([(r.canonical_smiles, r.source_tag) for r in members]),
            "standard_inchi": best_val([(r.inchi, r.source_tag) for r in members]),
            "standard_inchikey": best_val([(r.inchikey, r.source_tag) for r in members]),
            "sources": "|".join(sorted({r.source_tag for r in members}, key=lambda t: SOURCE_ORDER[t])),
            "source_count": len({r.source_tag for r in members}),
            "member_count": len(members),
            "merge_basis": gk[0],
            "structure_parse_status": "|".join(sorted({r.parse_status for r in members})),
            "name_aliases": aliases(r.name for r in members),
            "original_smiles_count": len({r.raw_smiles for r in members if r.raw_smiles}),
        })

        for r in members:
            src2gid[(r.source_tag, r.source_id)] = gid
            xref_rows.append({
                "global_compound_id": gid, "source": r.source_tag,
                "source_compound_id": r.source_id, "source_smiles": r.raw_smiles,
                "source_name": r.name, "canonical_smiles": r.canonical_smiles,
                "standard_inchikey": r.inchikey, "structure_parse_status": r.parse_status,
            })

    log(f"[compounds] {len(recs)} -> {len(g_rows)} global compounds")
    return g_rows, xref_rows, audit, src2gid


# ---------------------------------------------------------------------------
# Phase 4: Interaction resolution
# ---------------------------------------------------------------------------
def resolve_interactions(recs: list[InteractionRec],
                         enz_map: dict[tuple[str, str], str],
                         cmp_map: dict[tuple[str, str], str]):
    log(f"[interactions] {len(recs)} raw rows")
    pairs: dict[tuple[str, str], list[InteractionRec]] = defaultdict(list)
    audit, evidence = [], []
    skipped = 0

    for r in recs:
        ge = enz_map.get((r.source_tag, r.enzyme_id))
        gc = cmp_map.get((r.source_tag, r.compound_id))
        if not ge or not gc:
            skipped += 1
            audit.append(dict(entity_type="interaction", global_id="", source=r.source_tag,
                              source_local_id=r.source_id, merge_basis="skipped",
                              conflict_flags=f"missing_{'enzyme' if not ge else ''}{'_compound' if not gc else ''}",
                              notes=f"enz={r.enzyme_id} cmp={r.compound_id}"))
            continue
        pairs[(ge, gc)].append(r)

    g_rows = []
    for i, pk in enumerate(sorted(pairs), 1):
        gid = f"INT_G{i:06d}"
        members = sorted(pairs[pk], key=lambda r: (SOURCE_ORDER[r.source_tag], r.source_id))
        labs = sorted({r.label for r in members if r.label})

        if len(labs) > 1:
            audit.append(dict(entity_type="interaction", global_id=gid, source="",
                              source_local_id="", merge_basis="pair_merge",
                              conflict_flags="label_conflict", notes=f"labels={labs}"))

        tiers = [r.quality_tier for r in members if r.quality_tier]
        best_tier = sorted(tiers, key=lambda t: (QUALITY_ORDER.get(t, 50), t))[0] if tiers else ""

        g_rows.append({
            "global_interaction_id": gid,
            "global_enzyme_id": pk[0], "global_compound_id": pk[1],
            "label": labs[0] if len(labs) == 1 else "",
            "sources": "|".join(sorted({r.source_tag for r in members}, key=lambda t: SOURCE_ORDER[t])),
            "supporting_source_count": len({r.source_tag for r in members}),
            "evidence_count": len(members),
            "best_quality_tier": best_tier,
            "has_pmid_any": int(any(r.has_pmid for r in members)),
            "has_products_any": int(any(r.has_products for r in members)),
            "has_pdb_any": int(any(r.has_pdb for r in members)),
            "max_num_reactions": max((r.num_reactions for r in members), default=0),
            "conflicting_labels": int(len(labs) > 1),
        })

        for r in members:
            evidence.append({
                "global_interaction_id": gid, "source": r.source_tag,
                "source_interaction_id": r.source_id,
                "source_enzyme_id": r.enzyme_id, "source_compound_id": r.compound_id,
                "label": r.label, "quality_tier": r.quality_tier,
                "num_reactions": r.num_reactions, "has_pmid": r.has_pmid,
                "has_products": r.has_products, "pdb_id": r.pdb_id, "extra_json": r.extra,
            })

    if skipped:
        log(f"[interactions] WARNING: {skipped} rows skipped (missing enzyme/compound mapping)")
    log(f"[interactions] {len(recs)} -> {len(g_rows)} global interactions")
    return g_rows, evidence, audit


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Merge and deduplicate P450 sources.")
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--skip-optional", action="store_true")
    parser.add_argument("--strip-salts", action="store_true")
    args = parser.parse_args()

    inp = args.input_root.resolve()
    out = args.output_root.resolve()
    log(f"[setup] input={inp}")
    log(f"[setup] output={out}")
    log(f"[setup] rdkit={'yes' if HAS_RDKIT else 'NO'} strip_salts={args.strip_salts}")

    for d, _ in SOURCE_CONFIG:
        if not (inp / d).exists():
            raise FileNotFoundError(f"Missing: {inp / d}")

    # Resolve entities
    enz_recs = load_enzymes(inp)
    g_enz, enz_xref, enz_audit, enz_map = resolve_enzymes(enz_recs)

    cmp_recs = load_compounds(inp, args.strip_salts)
    g_cmp, cmp_xref, cmp_audit, cmp_map = resolve_compounds(cmp_recs)

    int_recs, lab_dist = load_interactions(inp)
    log(f"[interactions] label distribution: {dict(sorted(lab_dist.items()))}")
    g_int, int_evidence, int_audit = resolve_interactions(int_recs, enz_map, cmp_map)

    all_audit = sorted(enz_audit + cmp_audit + int_audit,
                       key=lambda r: (r.get("entity_type", ""), r.get("global_id", ""),
                                      SOURCE_ORDER.get(r.get("source", ""), 99)))

    # Write main outputs
    write_csv(out / "global_enzymes.csv", g_enz, [
        "global_enzyme_id", "canonical_uniprot_id", "canonical_p450_symbol",
        "canonical_species", "canonical_species_class", "canonical_ec_number",
        "canonical_sequence", "sequence_hash", "sequence_length",
        "sources", "source_count", "member_count", "merge_basis",
        "sequence_conflict", "species_conflict", "ec_conflict",
        "symbol_aliases", "species_aliases", "ec_number_all",
    ])
    write_csv(out / "global_compounds.csv", g_cmp, [
        "global_compound_id", "canonical_name", "canonical_smiles",
        "standard_inchi", "standard_inchikey",
        "sources", "source_count", "member_count", "merge_basis",
        "structure_parse_status", "name_aliases", "original_smiles_count",
    ])
    write_csv(out / "global_interactions.csv", g_int, [
        "global_interaction_id", "global_enzyme_id", "global_compound_id",
        "label", "sources", "supporting_source_count", "evidence_count",
        "best_quality_tier", "has_pmid_any", "has_products_any",
        "has_pdb_any", "max_num_reactions", "conflicting_labels",
    ])

    if not args.skip_optional:
        write_csv(out / "enzyme_xref.csv", enz_xref, [
            "global_enzyme_id", "source", "source_enzyme_id", "source_uniprot_id",
            "source_p450_symbol", "source_species", "source_ec_number", "source_sequence_hash",
        ])
        write_csv(out / "compound_xref.csv", cmp_xref, [
            "global_compound_id", "source", "source_compound_id", "source_smiles",
            "source_name", "canonical_smiles", "standard_inchikey", "structure_parse_status",
        ])
        write_csv(out / "interaction_evidence.csv", int_evidence, [
            "global_interaction_id", "source", "source_interaction_id",
            "source_enzyme_id", "source_compound_id", "label", "quality_tier",
            "num_reactions", "has_pmid", "has_products", "pdb_id", "extra_json",
        ])
        write_csv(out / "merge_audit.csv", all_audit, [
            "entity_type", "global_id", "source", "source_local_id",
            "merge_basis", "conflict_flags", "notes",
        ])

    # Summary
    log("=" * 60)
    log(f"[DONE] global enzymes:      {len(g_enz)}")
    log(f"[DONE] global compounds:    {len(g_cmp)}")
    log(f"[DONE] global interactions: {len(g_int)}")
    log(f"[DONE] audit rows:          {len(all_audit)}")
    log("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted")
        raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
