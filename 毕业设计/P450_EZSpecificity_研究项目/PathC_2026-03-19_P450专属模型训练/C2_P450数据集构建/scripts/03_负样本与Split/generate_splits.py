#!/usr/bin/env python3
"""
Phase 5 v3: Paper-aligned one-global-data.csv + 4 split families.

Design (matches EZSpecificity ESIBank):
  1. Generate ONE global negative set (single-direction: fix substrate, swap enzyme)
  2. Combine into ONE shared data.csv
  3. Split the SAME data.csv 4 ways:
     - random_split:   row-level random (strict)
     - reaction_split: substrate-group disjoint (strict)
     - enzyme_split:   enzyme sequence-hash-group disjoint (strict)
     - all_split:      soft joint-OOD (minimize enzyme+substrate overlap, report metrics)
"""
from __future__ import annotations

import argparse, csv, hashlib, json, random, sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

MASTER_SEED = 20260324
FOLDS = 4
DEFAULT_NEG_RATIO = 10

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
OUTPUT = PROJECT / "data" / "P450_Family"

SPLITS = ("random_split", "reaction_split", "enzyme_split", "all_split")
FIELDS = ["enzyme", "reaction", "label", "ecnumber", "difficulty", "fake_ecnumber", "structure_index"]

# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Enzyme:
    gid: str; seq: str; uniprot: str; seq_hash: str; ec: str

@dataclass(frozen=True)
class Compound:
    gid: str; smiles: str

@dataclass(frozen=True)
class PosInt:
    gid: str; enz_idx: int; sub_idx: int

@dataclass(frozen=True)
class Row:
    enzyme: int; reaction: int; label: int
    ecnumber: str; difficulty: int; fake_ecnumber: str; structure_index: int
    def to_dict(self):
        return {f: getattr(self, f) for f in FIELDS}

# ---------------------------------------------------------------------------
def log(m): print(m, file=sys.stderr, flush=True)
def norm(v):
    t = str(v).strip() if v is not None else ""
    return "" if t.lower() in {"","na","n/a","none","null"} else t

def make_seed(*p):
    return int(hashlib.sha256("|".join(str(x) for x in (MASTER_SEED,*p)).encode()).hexdigest()[:16], 16)
def make_rng(*p): return random.Random(make_seed(*p))

def read_csv(path, req):
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        missing = sorted(req - set(reader.fieldnames or []))
        if missing: raise ValueError(f"{path}: missing {missing}")
        return [dict(r) for r in reader]

def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for r in rows: w.writerow(r)

def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2, sort_keys=True)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load_enzymes(path):
    rows = read_csv(path, {"global_enzyme_id","canonical_sequence","canonical_uniprot_id","sequence_hash","canonical_ec_number"})
    enzymes, idx, sg = [], {}, defaultdict(list)
    for i, r in enumerate(rows):
        e = Enzyme(gid=norm(r["global_enzyme_id"]), seq=norm(r["canonical_sequence"]),
                   uniprot=norm(r["canonical_uniprot_id"]) or norm(r["global_enzyme_id"]),
                   seq_hash=norm(r["sequence_hash"]) or f"missing_{i}",
                   ec=norm(r["canonical_ec_number"]))
        enzymes.append(e); idx[e.gid] = i; sg[e.seq_hash].append(i)
    return enzymes, idx, dict(sg)

def load_compounds(path):
    rows = read_csv(path, {"global_compound_id","canonical_smiles"})
    compounds, idx = [], {}
    for i, r in enumerate(rows):
        c = Compound(gid=norm(r["global_compound_id"]), smiles=norm(r["canonical_smiles"]))
        compounds.append(c); idx[c.gid] = i
    return compounds, idx

def load_interactions(path, enz_idx, cmp_idx):
    rows = read_csv(path, {"global_interaction_id","global_enzyme_id","global_compound_id","label"})
    out = []
    for r in rows:
        if norm(r["label"]) != "1": raise ValueError(f"Non-positive: {r}")
        out.append(PosInt(gid=norm(r["global_interaction_id"]),
                          enz_idx=enz_idx[norm(r["global_enzyme_id"])],
                          sub_idx=cmp_idx[norm(r["global_compound_id"])]))
    return out

# ---------------------------------------------------------------------------
# Enzyme grouping
# ---------------------------------------------------------------------------
def build_enzyme_groups(seq_groups, positives):
    e2g = {}
    for sh, idxs in seq_groups.items():
        for i in idxs: e2g[i] = sh
    return e2g

# ---------------------------------------------------------------------------
# EC difficulty
# ---------------------------------------------------------------------------
def ec_difficulty(real_ec, fake_ec):
    if not real_ec or not fake_ec: return -1
    rp, fp = real_ec.split("."), fake_ec.split(".")
    if len(rp) < 4 or len(fp) < 4: return -1
    m = 0
    for a, b in zip(rp[:4], fp[:4]):
        if not a or not b or a == "-" or b == "-": break
        if a != b: break
        m += 1
    return m

# ---------------------------------------------------------------------------
# Global negative generation (single-direction, NO component restriction)
# ---------------------------------------------------------------------------
def generate_negatives(positives, enzymes, n_per_pos, max_retry=50):
    pos_pairs = {(p.enz_idx, p.sub_idx) for p in positives}
    pos_by_sub = defaultdict(set)
    for p in positives: pos_by_sub[p.sub_idx].add(p.enz_idx)

    all_enzyme_indices = list(range(len(enzymes)))
    used, negs = set(), []
    shortfall = 0

    for p in sorted(positives, key=lambda x: x.gid):
        real_ec = enzymes[p.enz_idx].ec
        rng = make_rng("neg", p.gid)
        got = 0
        for _ in range(n_per_pos * max_retry):
            if got >= n_per_pos: break
            fe = rng.choice(all_enzyme_indices)
            if fe == p.enz_idx: continue
            if fe in pos_by_sub[p.sub_idx]: continue
            pair = (fe, p.sub_idx)
            if pair in used: continue
            used.add(pair)
            fake_ec = enzymes[fe].ec
            negs.append(Row(enzyme=fe, reaction=p.sub_idx, label=0,
                            ecnumber=real_ec, difficulty=ec_difficulty(real_ec, fake_ec),
                            fake_ecnumber=fake_ec, structure_index=-1))
            got += 1
        if got < n_per_pos: shortfall += 1

    stats = {"target": len(positives)*n_per_pos, "generated": len(negs),
             "ratio": round(len(negs)/len(positives), 4) if positives else 0,
             "shortfall_positives": shortfall}
    return negs, stats

# ---------------------------------------------------------------------------
# Greedy fold assignment
# ---------------------------------------------------------------------------
def greedy_assign(weights, n=FOLDS):
    asgn, loads = {}, [0]*n
    for k, w in sorted(weights.items(), key=lambda x: (-x[1], str(x[0]))):
        f = min(range(n), key=lambda i: (loads[i], i))
        asgn[k] = f; loads[f] += w
    return asgn, loads

# ---------------------------------------------------------------------------
# Split assignment functions
# ---------------------------------------------------------------------------
def assign_random(rows):
    idxs = list(range(len(rows)))
    make_rng("random_split").shuffle(idxs)
    return {i: pos % FOLDS for pos, i in enumerate(idxs)}

def assign_by_substrate(rows):
    sw, sr = defaultdict(int), defaultdict(list)
    for i, r in enumerate(rows): sw[r.reaction] += 1; sr[r.reaction].append(i)
    s2f, _ = greedy_assign(dict(sw))
    return {i: s2f[r.reaction] for i, r in enumerate(rows)}, s2f

def assign_by_enzyme_group(rows, e2g):
    gw, gr = defaultdict(int), defaultdict(list)
    for i, r in enumerate(rows):
        g = e2g[r.enzyme]; gw[g] += 1; gr[g].append(i)
    g2f, _ = greedy_assign(dict(gw))
    return {i: g2f[e2g[r.enzyme]] for i, r in enumerate(rows)}, g2f

def assign_all_split_soft(rows, e2g):
    """
    Soft all_split: minimize enzyme+substrate overlap across folds.
    1. Independently assign enzyme groups and substrates to folds (greedy balanced)
    2. For each row, prefer fold matching both; if conflict, pick by penalty + balance
    All rows retained, overlap minimized but NOT guaranteed zero.
    """
    # Independent entity fold preferences
    egw = defaultdict(int)
    for r in rows: egw[e2g[r.enzyme]] += 1
    eg2f, _ = greedy_assign(dict(egw))

    sw = defaultdict(int)
    for r in rows: sw[r.reaction] += 1
    s2f, _ = greedy_assign(dict(sw))

    # Assign each row: prefer fold matching both, then one, then lightest
    loads = [0]*FOLDS
    r2f = {}
    # Sort by descending row weight to assign heavy items first (greedy)
    row_order = sorted(range(len(rows)), key=lambda i: (
        -(1 if eg2f[e2g[rows[i].enzyme]] == s2f[rows[i].reaction] else 0),
        i
    ))

    for i in row_order:
        r = rows[i]
        ef = eg2f[e2g[r.enzyme]]
        sf = s2f[r.reaction]
        if ef == sf:
            best = ef
        else:
            # Prefer the fold with fewer rows among {ef, sf}
            candidates = sorted([ef, sf], key=lambda f: (loads[f], f))
            best = candidates[0]
        r2f[i] = best
        loads[best] += 1

    return r2f, eg2f, s2f

# ---------------------------------------------------------------------------
# Partition + validation + overlap metrics
# ---------------------------------------------------------------------------
def make_partitions(rows, r2f):
    by_fold = defaultdict(list)
    for i, f in r2f.items(): by_fold[f].append(rows[i])
    parts = {}
    for k in range(FOLDS):
        te, va = k, (k+1)%FOLDS
        tr_folds = [f for f in range(FOLDS) if f not in {te, va}]
        train = []
        for f in tr_folds: train.extend(by_fold.get(f, []))
        srt = lambda rs: sorted(rs, key=lambda r: (r.reaction, r.enzyme, -r.label))
        parts[k] = {"train": srt(train), "val": srt(by_fold.get(va,[])), "test": srt(by_fold.get(te,[]))}
    return parts

def validate_strict(split_type, parts, e2g):
    issues = []
    for k, p in parts.items():
        tr, va, te = p["train"], p["val"], p["test"]
        if split_type in ("enzyme_split","all_split"):
            trg = {e2g[r.enzyme] for r in tr}; vag = {e2g[r.enzyme] for r in va}; teg = {e2g[r.enzyme] for r in te}
            if trg & teg: issues.append(f"fold{k}: train/test enzyme leak ({len(trg&teg)} groups)")
            if trg & vag: issues.append(f"fold{k}: train/val enzyme leak ({len(trg&vag)} groups)")
        if split_type in ("reaction_split","all_split"):
            trs = {r.reaction for r in tr}; vas = {r.reaction for r in va}; tes = {r.reaction for r in te}
            if trs & tes: issues.append(f"fold{k}: train/test substrate leak ({len(trs&tes)} substrates)")
            if trs & vas: issues.append(f"fold{k}: train/val substrate leak ({len(trs&vas)} substrates)")
    return issues

def compute_overlap_metrics(parts, e2g):
    """Compute enzyme/substrate overlap for all_split (soft)."""
    metrics = {}
    for k, p in parts.items():
        tr, va, te = p["train"], p["val"], p["test"]
        tr_eg = {e2g[r.enzyme] for r in tr}; va_eg = {e2g[r.enzyme] for r in va}; te_eg = {e2g[r.enzyme] for r in te}
        tr_s = {r.reaction for r in tr}; va_s = {r.reaction for r in va}; te_s = {r.reaction for r in te}
        metrics[str(k)] = {
            "train_test_enzyme_overlap": len(tr_eg & te_eg),
            "train_test_enzyme_total": len(tr_eg | te_eg),
            "train_test_enzyme_overlap_frac": round(len(tr_eg & te_eg) / max(len(tr_eg | te_eg), 1), 4),
            "train_test_substrate_overlap": len(tr_s & te_s),
            "train_test_substrate_total": len(tr_s | te_s),
            "train_test_substrate_overlap_frac": round(len(tr_s & te_s) / max(len(tr_s | te_s), 1), 4),
            "train_val_enzyme_overlap": len(tr_eg & va_eg),
            "train_val_substrate_overlap": len(tr_s & va_s),
        }
    return metrics

# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------
def write_split(split_type, rows, parts, out_dir):
    sd = out_dir / split_type; sd.mkdir(parents=True, exist_ok=True)
    write_csv(sd / "data.csv", [r.to_dict() for r in rows], FIELDS)
    fold_stats = {}
    for k in range(FOLDS):
        fold_stats[k] = {}
        for pn in ("train","val","test"):
            pr = parts[k][pn]
            pos = sum(1 for r in pr if r.label==1)
            neg = sum(1 for r in pr if r.label==0)
            fn = {"train":f"training_datas_{k}.csv","val":f"val_datas_{k}.csv","test":f"testing_datas_{k}.csv"}[pn]
            write_csv(sd / fn, [r.to_dict() for r in pr], FIELDS)
            fold_stats[k][pn] = {"pos": pos, "neg": neg, "total": pos+neg}
            log(f"  [{split_type}] fold={k} {pn}: {pos}+ {neg}-")
    return fold_stats

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    global MASTER_SEED
    parser = argparse.ArgumentParser()
    parser.add_argument("--combined-dir", type=Path, default=COMBINED)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--seed", type=int, default=MASTER_SEED)
    parser.add_argument("--neg-ratio", type=int, default=DEFAULT_NEG_RATIO)
    args = parser.parse_args()
    MASTER_SEED = args.seed

    cdir, odir = args.combined_dir.resolve(), args.output_dir.resolve()
    log(f"[setup] seed={MASTER_SEED} neg_ratio={args.neg_ratio}")

    enzymes, ei, sg = load_enzymes(cdir / "global_enzymes.csv")
    compounds, ci = load_compounds(cdir / "global_compounds.csv")
    positives = load_interactions(cdir / "global_interactions.csv", ei, ci)
    e2g = build_enzyme_groups(sg, positives)
    log(f"[data] {len(enzymes)} enzymes, {len(compounds)} compounds, {len(positives)} positives")

    # Positive rows
    pos_rows = [Row(enzyme=p.enz_idx, reaction=p.sub_idx, label=1,
                    ecnumber=enzymes[p.enz_idx].ec, difficulty=0,
                    fake_ecnumber="", structure_index=-1)
                for p in sorted(positives, key=lambda x: x.gid)]

    # Global negatives (single-direction, no component restriction)
    neg_rows, neg_stats = generate_negatives(positives, enzymes, args.neg_ratio)
    log(f"[negatives] {neg_stats['generated']}/{neg_stats['target']} (ratio={neg_stats['ratio']})")

    # ONE master data.csv
    master = sorted(pos_rows + neg_rows, key=lambda r: (r.reaction, r.enzyme, -r.label))
    log(f"[master] {len(master)} rows ({len(pos_rows)} pos + {len(neg_rows)} neg)")

    # Shared entity files
    write_csv(odir / "Enzymes.csv",
              [{"Protein sequence": e.seq, "uniprots": e.uniprot} for e in enzymes],
              ["Protein sequence", "uniprots"])
    write_csv(odir / "Substrates.csv",
              [{"Substrate_SMILES": c.smiles} for c in compounds],
              ["Substrate_SMILES"])

    # Assignments
    a_random = assign_random(master)
    a_reaction, s2f_reaction = assign_by_substrate(master)
    a_enzyme, g2f_enzyme = assign_by_enzyme_group(master, e2g)
    a_all, eg2f_all, s2f_all = assign_all_split_soft(master, e2g)

    assigns = {"random_split": a_random, "reaction_split": a_reaction,
               "enzyme_split": a_enzyme, "all_split": a_all}

    stats = {"seed": MASTER_SEED, "neg_ratio": args.neg_ratio,
             "input": {"enzymes": len(enzymes), "compounds": len(compounds),
                       "positives": len(pos_rows), "negatives": len(neg_rows),
                       "master_rows": len(master)},
             "negative_generation": neg_stats, "splits": {}}

    for stype in SPLITS:
        log(f"\n[{stype}]")
        parts = make_partitions(master, assigns[stype])

        if stype == "all_split":
            # Soft: report overlap metrics instead of strict validation
            overlap = compute_overlap_metrics(parts, e2g)
            log(f"  all_split overlap metrics:")
            for k, m in overlap.items():
                log(f"    fold{k}: enzyme overlap={m['train_test_enzyme_overlap_frac']*100:.1f}%, substrate overlap={m['train_test_substrate_overlap_frac']*100:.1f}%")
            issues = []  # soft split, no strict leakage check
            extra = {"overlap_metrics": overlap}
        else:
            issues = validate_strict(stype, parts, e2g)
            log(f"  leakage: {'PASS' if not issues else 'FAIL '+str(issues)}")
            extra = {}

        fs = write_split(stype, master, parts, odir)
        stats["splits"][stype] = {"leakage": issues, "folds": {str(k): v for k, v in fs.items()}, **extra}

    write_json(cdir / "negatives" / "split_stats.json", stats)
    log(f"\n{'='*60}")
    log(f"[DONE] {len(master)} rows x 4 splits (same data.csv)")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
