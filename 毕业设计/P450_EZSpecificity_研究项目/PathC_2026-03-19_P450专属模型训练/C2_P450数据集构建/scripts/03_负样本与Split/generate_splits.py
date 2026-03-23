#!/usr/bin/env python3
"""
Generate EZSpecificity-formatted 4-way splits with bidirectional negatives.

Inputs:  data/combined/global_{enzymes,compounds,interactions}.csv
Outputs: data/P450_Family/{Enzymes,Substrates}.csv + 4 split dirs
         data/combined/negatives/split_stats.json
"""
from __future__ import annotations

import argparse, csv, hashlib, json, random, sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

MASTER_SEED = 20260324
FOLDS = 4
NEG_PER_DIR = 5          # negatives per direction per positive
MAX_RETRY = 50           # retries per single negative sample
DOCK_PLACEHOLDER = -1

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
OUTPUT = PROJECT / "data" / "P450_Family"

SPLITS = ("random_split", "reaction_split", "enzyme_split", "all_split")
PARTS = ("train", "val", "test")

# ---------------------------------------------------------------------------
# Data types
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Enzyme:
    gid: str; seq: str; uniprot: str; seq_hash: str

@dataclass(frozen=True)
class Compound:
    gid: str; smiles: str

@dataclass(frozen=True)
class PosInt:
    gid: str; enz_id: str; cmp_id: str; enz_idx: int; sub_idx: int

@dataclass(frozen=True)
class Row:
    sub: int; enz: int; label: int; dock: int; pos_rxn: int
    def to_dict(self):
        return {"Substrate Index": self.sub, "Enzyme Index": self.enz,
                "Label": self.label, "Dock Index": self.dock,
                "positive_reactions": self.pos_rxn}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def log(msg): print(msg, file=sys.stderr, flush=True)
def norm(v):
    t = str(v).strip() if v is not None else ""
    return "" if t.lower() in {"","na","n/a","none","null"} else t

def make_seed(*parts):
    h = hashlib.sha256("|".join(str(p) for p in (MASTER_SEED, *parts)).encode()).hexdigest()
    return int(h[:16], 16)

def make_rng(*parts): return random.Random(make_seed(*parts))

def read_csv(path, req):
    with path.open("r", encoding="utf-8-sig", newline="", errors="replace") as f:
        reader = csv.DictReader(f)
        missing = sorted(req - set(reader.fieldnames or []))
        if missing: raise ValueError(f"{path}: missing {missing}")
        return [dict(r) for r in reader]

def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows: w.writerow(r)

def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2, sort_keys=True)

# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
def load_enzymes(path):
    rows = read_csv(path, {"global_enzyme_id","canonical_sequence","canonical_uniprot_id","sequence_hash"})
    enzymes, idx_map, seq_groups = [], {}, defaultdict(list)
    for i, r in enumerate(rows):
        e = Enzyme(gid=norm(r["global_enzyme_id"]), seq=norm(r["canonical_sequence"]),
                   uniprot=norm(r["canonical_uniprot_id"]) or norm(r["global_enzyme_id"]),
                   seq_hash=norm(r["sequence_hash"]))
        enzymes.append(e); idx_map[e.gid] = i; seq_groups[e.seq_hash].append(i)
    return enzymes, idx_map, dict(seq_groups)

def load_compounds(path):
    rows = read_csv(path, {"global_compound_id","canonical_smiles"})
    compounds, idx_map = [], {}
    for i, r in enumerate(rows):
        c = Compound(gid=norm(r["global_compound_id"]), smiles=norm(r["canonical_smiles"]))
        compounds.append(c); idx_map[c.gid] = i
    return compounds, idx_map

def load_interactions(path, enz_idx, cmp_idx):
    rows = read_csv(path, {"global_interaction_id","global_enzyme_id","global_compound_id","label"})
    out = []
    for r in rows:
        if norm(r["label"]) != "1":
            raise ValueError(f"Non-positive label: {r}")
        eid, cid = norm(r["global_enzyme_id"]), norm(r["global_compound_id"])
        out.append(PosInt(gid=norm(r["global_interaction_id"]),
                          enz_id=eid, cmp_id=cid,
                          enz_idx=enz_idx[eid], sub_idx=cmp_idx[cid]))
    return out

# ---------------------------------------------------------------------------
# Fold assignment
# ---------------------------------------------------------------------------
def greedy_assign(group_weights: dict, n_folds=FOLDS):
    """Greedy bin-packing: heaviest-first → lightest fold."""
    assignments = {}
    loads = [0]*n_folds
    for gk, w in sorted(group_weights.items(), key=lambda x: (-x[1], str(x[0]))):
        fold = min(range(n_folds), key=lambda f: (loads[f], f))
        assignments[gk] = fold
        loads[fold] += w
    return assignments, loads

def assign_random(interactions):
    ids = sorted(p.gid for p in interactions)
    rng = make_rng("random_split")
    rng.shuffle(ids)
    mapping = {}
    for i, gid in enumerate(ids):
        mapping[gid] = i % FOLDS
    return mapping

def assign_enzyme(interactions, enz_to_group, group_weights):
    g2f, loads = greedy_assign(group_weights)
    return {p.gid: g2f[enz_to_group[p.enz_idx]] for p in interactions}, loads

def assign_reaction(interactions):
    sub_w = defaultdict(int)
    for p in interactions: sub_w[p.sub_idx] += 1
    s2f, loads = greedy_assign(dict(sub_w))
    return {p.gid: s2f[p.sub_idx] for p in interactions}, loads

def assign_all_split(interactions, enz_to_group, group_weights):
    """Optimized all_split: greedy enzyme assignment + substrate follows enzymes."""
    g2f, e_loads = greedy_assign(group_weights)

    # For each substrate, count interactions per enzyme-fold
    sub_to_efold = defaultdict(Counter)
    sub_total = defaultdict(int)
    for p in interactions:
        ef = g2f[enz_to_group[p.enz_idx]]
        sub_to_efold[p.sub_idx][ef] += 1
        sub_total[p.sub_idx] += 1

    # Assign substrates: maximize same-fold interactions, respect balance
    target = len(interactions) / FOLDS
    cap = max(int(target * 1.15), max(sub_total.values(), default=0))
    s_loads = [0]*FOLDS
    s2f = {}
    for sidx, total in sorted(sub_total.items(), key=lambda x: (-x[1], x[0])):
        counts = sub_to_efold[sidx]
        candidates = [f for f in range(FOLDS) if s_loads[f] + total <= cap] or list(range(FOLDS))
        best = max(candidates, key=lambda f: (counts.get(f, 0), -s_loads[f], -f))
        s2f[sidx] = best
        s_loads[best] += total

    # Keep only interactions where enzyme_fold == substrate_fold
    mapping = {}
    retained, dropped = 0, 0
    per_fold = [0]*FOLDS
    for p in interactions:
        ef = g2f[enz_to_group[p.enz_idx]]
        sf = s2f[p.sub_idx]
        if ef == sf:
            mapping[p.gid] = ef
            retained += 1
            per_fold[ef] += 1
        else:
            dropped += 1

    stats = {"total": len(interactions), "retained": retained, "dropped": dropped,
             "fraction": round(retained/len(interactions), 4),
             "per_fold": per_fold, "e_loads": e_loads, "s_loads": s_loads}
    log(f"[all_split] retained={retained}/{len(interactions)} ({stats['fraction']*100:.1f}%), dropped={dropped}")
    return mapping, stats

# ---------------------------------------------------------------------------
# Split partitions
# ---------------------------------------------------------------------------
def make_partitions(interactions, fold_map):
    by_fold = defaultdict(list)
    for p in interactions:
        if p.gid in fold_map:
            by_fold[fold_map[p.gid]].append(p)
    parts = {}
    for k in range(FOLDS):
        test_f, val_f = k, (k+1) % FOLDS
        train_fs = [f for f in range(FOLDS) if f not in {test_f, val_f}]
        train = []
        for f in train_fs: train.extend(by_fold.get(f, []))
        parts[k] = {"train": sorted(train, key=lambda p: p.gid),
                     "val": sorted(by_fold.get(val_f, []), key=lambda p: p.gid),
                     "test": sorted(by_fold.get(test_f, []), key=lambda p: p.gid)}
    return parts

# ---------------------------------------------------------------------------
# Negative generation
# ---------------------------------------------------------------------------
def gen_negatives(positives, split_type, fold_k, part_name, all_pos_pairs):
    if not positives:
        return [], {"pos": 0, "neg": 0, "a_gen": 0, "b_gen": 0}
    enz_set = sorted({p.enz_idx for p in positives})
    sub_set = sorted({p.sub_idx for p in positives})
    rng = make_rng(split_type, fold_k, part_name)
    used = set()
    negs = []
    a_gen = b_gen = 0

    for p in positives:
        # Direction A: fix substrate, swap enzyme
        if len(enz_set) > 1:
            got = 0
            for _ in range(NEG_PER_DIR * MAX_RETRY):
                if got >= NEG_PER_DIR: break
                e = rng.choice(enz_set)
                if e == p.enz_idx: continue
                pair = (e, p.sub_idx)
                if pair in all_pos_pairs or pair in used: continue
                used.add(pair)
                negs.append(Row(sub=p.sub_idx, enz=e, label=0, dock=DOCK_PLACEHOLDER, pos_rxn=-1))
                got += 1; a_gen += 1

        # Direction B: fix enzyme, swap substrate
        if len(sub_set) > 1:
            got = 0
            for _ in range(NEG_PER_DIR * MAX_RETRY):
                if got >= NEG_PER_DIR: break
                s = rng.choice(sub_set)
                if s == p.sub_idx: continue
                pair = (p.enz_idx, s)
                if pair in all_pos_pairs or pair in used: continue
                used.add(pair)
                negs.append(Row(sub=s, enz=p.enz_idx, label=0, dock=DOCK_PLACEHOLDER, pos_rxn=-1))
                got += 1; b_gen += 1

    return negs, {"pos": len(positives), "neg": len(negs),
                  "a_gen": a_gen, "b_gen": b_gen,
                  "a_target": len(positives)*NEG_PER_DIR,
                  "b_target": len(positives)*NEG_PER_DIR,
                  "enzymes": len(enz_set), "substrates": len(sub_set)}

# ---------------------------------------------------------------------------
# Leakage validation
# ---------------------------------------------------------------------------
def validate_leakage(split_type, partitions, seq_groups):
    idx2group = {}
    for sh, idxs in seq_groups.items():
        for i in idxs: idx2group[i] = sh
    issues = []
    for k, parts in partitions.items():
        train, val, test = parts["train"], parts["val"], parts["test"]
        if split_type in ("enzyme_split", "all_split"):
            tr_g = {idx2group[p.enz_idx] for p in train}
            va_g = {idx2group[p.enz_idx] for p in val}
            te_g = {idx2group[p.enz_idx] for p in test}
            if tr_g & te_g: issues.append(f"fold{k}: train/test enzyme leak")
            if tr_g & va_g: issues.append(f"fold{k}: train/val enzyme leak")
        if split_type in ("reaction_split", "all_split"):
            tr_s = {p.sub_idx for p in train}
            va_s = {p.sub_idx for p in val}
            te_s = {p.sub_idx for p in test}
            if tr_s & te_s: issues.append(f"fold{k}: train/test substrate leak")
            if tr_s & va_s: issues.append(f"fold{k}: train/val substrate leak")
    return issues

# ---------------------------------------------------------------------------
# Write split files
# ---------------------------------------------------------------------------
FIELDS = ["Substrate Index", "Enzyme Index", "Label", "Dock Index", "positive_reactions"]

def write_split(split_type, partitions, out_dir, pos_pairs):
    split_dir = out_dir / split_type
    all_rows = []
    fold_stats = {}
    for k in range(FOLDS):
        fold_stats[k] = {}
        for part_name in PARTS:
            positives = partitions[k][part_name]
            pos_rows = [Row(sub=p.sub_idx, enz=p.enz_idx, label=1,
                            dock=DOCK_PLACEHOLDER, pos_rxn=p.sub_idx) for p in positives]
            neg_rows, ns = gen_negatives(positives, split_type, k, part_name, pos_pairs)
            combined = sorted(pos_rows + neg_rows, key=lambda r: (r.sub, r.enz, -r.label))
            all_rows.extend(combined)
            fname = {"train": f"training_datas_{k}.csv",
                     "val": f"val_datas_{k}.csv",
                     "test": f"testing_datas_{k}.csv"}[part_name]
            write_csv(split_dir / fname, [r.to_dict() for r in combined], FIELDS)
            fold_stats[k][part_name] = ns
            log(f"  [{split_type}] fold={k} {part_name}: {ns['pos']}+ {ns['neg']}-")

    # data.csv = deduplicated union
    unique = {}
    for r in all_rows:
        unique[(r.sub, r.enz, r.label)] = r
    data_rows = sorted(unique.values(), key=lambda r: (r.sub, r.enz, -r.label))
    write_csv(split_dir / "data.csv", [r.to_dict() for r in data_rows], FIELDS)
    log(f"  [{split_type}] data.csv: {len(data_rows)} rows")
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
    args = parser.parse_args()
    MASTER_SEED = args.seed

    cdir, odir = args.combined_dir.resolve(), args.output_dir.resolve()
    log(f"[setup] combined={cdir}")
    log(f"[setup] output={odir}")

    enzymes, enz_idx, seq_groups = load_enzymes(cdir / "global_enzymes.csv")
    compounds, cmp_idx = load_compounds(cdir / "global_compounds.csv")
    positives = load_interactions(cdir / "global_interactions.csv", enz_idx, cmp_idx)
    pos_pairs = {(p.enz_idx, p.sub_idx) for p in positives}
    log(f"[data] {len(enzymes)} enzymes, {len(compounds)} compounds, {len(positives)} positives")

    # Build enzyme sequence-hash groups
    enz2group = {}
    group_w = defaultdict(int)
    for sh, idxs in seq_groups.items():
        for i in idxs: enz2group[i] = sh
    for p in positives: group_w[enz2group[p.enz_idx]] += 1

    # Write shared entity files
    write_csv(odir / "Enzymes.csv",
              [{"Protein sequence": e.seq, "uniprots": e.uniprot} for e in enzymes],
              ["Protein sequence", "uniprots"])
    write_csv(odir / "Substrates.csv",
              [{"Substrate_SMILES": c.smiles} for c in compounds],
              ["Substrate_SMILES"])
    log(f"[write] Enzymes.csv ({len(enzymes)}), Substrates.csv ({len(compounds)})")

    # Assign folds
    log("[splits] assigning folds...")
    fold_maps = {}
    all_stats = {}

    fold_maps["random_split"] = assign_random(positives)
    fold_maps["enzyme_split"], e_loads = assign_enzyme(positives, enz2group, dict(group_w))
    fold_maps["reaction_split"], r_loads = assign_reaction(positives)
    fold_maps["all_split"], all_stats["all_split"] = assign_all_split(positives, enz2group, dict(group_w))

    # Build partitions and validate
    stats = {"seed": MASTER_SEED, "input": {"enzymes": len(enzymes), "compounds": len(compounds),
             "positives": len(positives), "seq_groups_multi": sum(1 for v in seq_groups.values() if len(v)>1)},
             "splits": {}}

    for stype in SPLITS:
        log(f"\n[{stype}] building partitions...")
        parts = make_partitions(positives, fold_maps[stype])
        issues = validate_leakage(stype, parts, seq_groups)
        if issues:
            log(f"  WARNING: leakage detected: {issues}")
        else:
            log(f"  leakage check: PASS")

        fold_stats = write_split(stype, parts, odir, pos_pairs)
        stats["splits"][stype] = {
            "leakage_issues": issues,
            "fold_stats": {str(k): v for k, v in fold_stats.items()},
        }
        if stype == "all_split":
            stats["splits"]["all_split"]["assignment"] = all_stats.get("all_split", {})

    write_json(cdir / "negatives" / "split_stats.json", stats)
    log(f"\n{'='*60}")
    log(f"[DONE] Output: {odir}")
    log(f"[DONE] Stats: {cdir / 'negatives' / 'split_stats.json'}")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}")
        import traceback; traceback.print_exc()
        raise SystemExit(1)
