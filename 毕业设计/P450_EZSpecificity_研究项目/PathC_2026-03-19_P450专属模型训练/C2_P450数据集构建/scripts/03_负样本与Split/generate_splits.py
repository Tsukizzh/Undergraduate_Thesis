#!/usr/bin/env python3
"""
Phase 5 v3: Paper-aligned one-global-data.csv + 4 split families.

Design (matches EZSpecificity ESIBank):
  1. Generate ONE global negative set (bidirectional: 5A swap enzyme + 5B swap substrate)
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

SHARED_SPLITS = ("random_split", "reaction_split", "enzyme_split")
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
# Global negative generation (BIDIRECTIONAL, NO component restriction)
#   Direction A: fix substrate, swap enzyme (5 per positive)
#   Direction B: fix enzyme, swap substrate (5 per positive)
# ---------------------------------------------------------------------------
def generate_negatives(positives, enzymes, compounds, n_per_dir=5, max_retry=50):
    pos_pairs = {(p.enz_idx, p.sub_idx) for p in positives}
    pos_by_sub = defaultdict(set)  # substrate → set of positive enzymes
    pos_by_enz = defaultdict(set)  # enzyme → set of positive substrates
    for p in positives:
        pos_by_sub[p.sub_idx].add(p.enz_idx)
        pos_by_enz[p.enz_idx].add(p.sub_idx)

    all_enz = list(range(len(enzymes)))
    all_sub = list(range(len(compounds)))
    used, negs = set(), []
    a_gen = b_gen = a_short = b_short = 0

    for p in sorted(positives, key=lambda x: x.gid):
        real_ec = enzymes[p.enz_idx].ec
        rng = make_rng("neg", p.gid)

        # Direction A: fix substrate, swap enzyme
        got = 0
        for _ in range(n_per_dir * max_retry):
            if got >= n_per_dir: break
            fe = rng.choice(all_enz)
            if fe == p.enz_idx: continue
            if fe in pos_by_sub[p.sub_idx]: continue
            pair = (fe, p.sub_idx)
            if pair in used: continue
            used.add(pair)
            fake_ec = enzymes[fe].ec
            negs.append(Row(enzyme=fe, reaction=p.sub_idx, label=0,
                            ecnumber=real_ec, difficulty=ec_difficulty(real_ec, fake_ec),
                            fake_ecnumber=fake_ec, structure_index=-1))
            got += 1; a_gen += 1
        if got < n_per_dir: a_short += 1

        # Direction B: fix enzyme, swap substrate
        got = 0
        for _ in range(n_per_dir * max_retry):
            if got >= n_per_dir: break
            fs = rng.choice(all_sub)
            if fs == p.sub_idx: continue
            if fs in pos_by_enz[p.enz_idx]: continue
            pair = (p.enz_idx, fs)
            if pair in used: continue
            used.add(pair)
            negs.append(Row(enzyme=p.enz_idx, reaction=fs, label=0,
                            ecnumber=real_ec, difficulty=-1,
                            fake_ecnumber="", structure_index=-1))
            got += 1; b_gen += 1
        if got < n_per_dir: b_short += 1

    total = a_gen + b_gen
    stats = {"target_per_dir": len(positives)*n_per_dir,
             "a_generated": a_gen, "b_generated": b_gen, "total": total,
             "ratio": round(total/len(positives), 4) if positives else 0,
             "a_shortfall": a_short, "b_shortfall": b_short}
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

def build_strict_all_split(positives, enzymes, compounds, e2g, n_per_dir=5):
    """
    Strict all_split: 0% enzyme AND substrate overlap across folds.
    Uses optimized assignment: greedy enzyme groups → substrates follow enzymes.
    Rows where enzyme_fold != substrate_fold are DROPPED.
    Returns its own separate data.csv (positives + regenerated negatives within subset).
    """
    # Step 1: Assign enzyme groups to folds (greedy balanced on positive count)
    group_pos_w = defaultdict(int)
    for p in positives: group_pos_w[e2g[p.enz_idx]] += 1
    eg2f, e_loads = greedy_assign(dict(group_pos_w))

    # Step 2: Assign substrates — maximize same-fold as their enzymes
    sub_to_efold = defaultdict(Counter)
    sub_total = defaultdict(int)
    for p in positives:
        ef = eg2f[e2g[p.enz_idx]]
        sub_to_efold[p.sub_idx][ef] += 1
        sub_total[p.sub_idx] += 1

    target = len(positives) / FOLDS
    cap = max(int(target * 1.15), max(sub_total.values(), default=0))
    s_loads = [0]*FOLDS
    s2f = {}
    for sidx, total in sorted(sub_total.items(), key=lambda x: (-x[1], x[0])):
        counts = sub_to_efold[sidx]
        candidates = [f for f in range(FOLDS) if s_loads[f] + total <= cap] or list(range(FOLDS))
        best = max(candidates, key=lambda f: (counts.get(f, 0), -s_loads[f], -f))
        s2f[sidx] = best
        s_loads[best] += total

    # Step 3: Keep only positives where enzyme_fold == substrate_fold
    kept_pos, dropped = [], 0
    per_fold = [0]*FOLDS
    for p in positives:
        ef = eg2f[e2g[p.enz_idx]]
        sf = s2f[p.sub_idx]
        if ef == sf:
            kept_pos.append((p, ef))
            per_fold[ef] += 1
        else:
            dropped += 1

    # Step 4: Generate negatives WITHIN the strict subset
    # For each kept positive, sample negatives only from same-fold enzymes/substrates
    fold_enzymes = defaultdict(list)  # fold → enzyme indices in that fold
    fold_substrates = defaultdict(list)  # fold → substrate indices in that fold
    for sh, idxs in defaultdict(list).items(): pass  # placeholder
    # Build fold membership
    for sh, fold in eg2f.items():
        # sh is a sequence hash group, get all enzyme indices in it
        pass
    # Actually build from kept positives
    for p, fold in kept_pos:
        fold_enzymes[fold].append(p.enz_idx)
        fold_substrates[fold].append(p.sub_idx)
    for f in range(FOLDS):
        fold_enzymes[f] = sorted(set(fold_enzymes[f]))
        fold_substrates[f] = sorted(set(fold_substrates[f]))

    pos_pairs = {(p.enz_idx, p.sub_idx) for p in positives}  # global positives for exclusion
    pos_by_sub = defaultdict(set)
    pos_by_enz = defaultdict(set)
    for p in positives:
        pos_by_sub[p.sub_idx].add(p.enz_idx)
        pos_by_enz[p.enz_idx].add(p.sub_idx)

    all_rows = []
    used = set()
    neg_a = neg_b = 0

    for p, fold in sorted(kept_pos, key=lambda x: x[0].gid):
        # Positive row
        all_rows.append(Row(enzyme=p.enz_idx, reaction=p.sub_idx, label=1,
                            ecnumber=enzymes[p.enz_idx].ec, difficulty=0,
                            fake_ecnumber="", structure_index=-1))
        rng = make_rng("allsplit_neg", p.gid)

        # Direction A: fix substrate, swap enzyme (from same fold)
        fe_pool = fold_enzymes[fold]
        got = 0
        for _ in range(n_per_dir * 50):
            if got >= n_per_dir: break
            fe = rng.choice(fe_pool)
            if fe == p.enz_idx or fe in pos_by_sub[p.sub_idx]: continue
            pair = (fe, p.sub_idx)
            if pair in used: continue
            used.add(pair)
            fake_ec = enzymes[fe].ec
            all_rows.append(Row(enzyme=fe, reaction=p.sub_idx, label=0,
                                ecnumber=enzymes[p.enz_idx].ec,
                                difficulty=ec_difficulty(enzymes[p.enz_idx].ec, fake_ec),
                                fake_ecnumber=fake_ec, structure_index=-1))
            got += 1; neg_a += 1

        # Direction B: fix enzyme, swap substrate (from same fold)
        fs_pool = fold_substrates[fold]
        got = 0
        for _ in range(n_per_dir * 50):
            if got >= n_per_dir: break
            fs = rng.choice(fs_pool)
            if fs == p.sub_idx or fs in pos_by_enz[p.enz_idx]: continue
            pair = (p.enz_idx, fs)
            if pair in used: continue
            used.add(pair)
            all_rows.append(Row(enzyme=p.enz_idx, reaction=fs, label=0,
                                ecnumber=enzymes[p.enz_idx].ec, difficulty=-1,
                                fake_ecnumber="", structure_index=-1))
            got += 1; neg_b += 1

    all_rows.sort(key=lambda r: (r.reaction, r.enzyme, -r.label))
    n_pos = sum(1 for r in all_rows if r.label == 1)
    n_neg = sum(1 for r in all_rows if r.label == 0)

    stats = {
        "original_positives": len(positives),
        "retained_positives": len(kept_pos), "dropped_positives": dropped,
        "retention_rate": round(len(kept_pos)/len(positives), 4),
        "per_fold_positives": per_fold,
        "negatives_a": neg_a, "negatives_b": neg_b, "negatives_total": neg_a+neg_b,
        "total_rows": len(all_rows), "pos": n_pos, "neg": n_neg,
        "ratio": round(n_neg/max(n_pos,1), 4),
    }
    return all_rows, stats, eg2f, s2f

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

    # Global negatives (bidirectional: 5A swap enzyme + 5B swap substrate)
    neg_rows, neg_stats = generate_negatives(positives, enzymes, compounds, n_per_dir=args.neg_ratio // 2)
    log(f"[negatives] A={neg_stats['a_generated']} B={neg_stats['b_generated']} total={neg_stats['total']} (ratio={neg_stats['ratio']})")

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

    # Assignments for shared splits
    a_random = assign_random(master)
    a_reaction, s2f_reaction = assign_by_substrate(master)
    a_enzyme, g2f_enzyme = assign_by_enzyme_group(master, e2g)

    assigns = {"random_split": a_random, "reaction_split": a_reaction,
               "enzyme_split": a_enzyme}

    stats = {"seed": MASTER_SEED, "neg_ratio": args.neg_ratio,
             "input": {"enzymes": len(enzymes), "compounds": len(compounds),
                       "positives": len(pos_rows), "negatives": len(neg_rows),
                       "master_rows": len(master)},
             "negative_generation": neg_stats, "splits": {}}

    # --- random/reaction/enzyme: shared data.csv ---
    for stype in ("random_split", "reaction_split", "enzyme_split"):
        log(f"\n[{stype}]")
        parts = make_partitions(master, assigns[stype])
        issues = validate_strict(stype, parts, e2g)
        log(f"  leakage: {'PASS' if not issues else 'FAIL '+str(issues)}")
        fs = write_split(stype, master, parts, odir)
        stats["splits"][stype] = {"leakage": issues, "folds": {str(k): v for k, v in fs.items()}}

    # --- all_split: STRICT, separate data.csv ---
    log(f"\n[all_split] (strict, separate dataset)")
    all_rows, all_stats, all_eg2f, all_s2f = build_strict_all_split(
        positives, enzymes, compounds, e2g, n_per_dir=args.neg_ratio // 2)
    log(f"  retained {all_stats['retained_positives']}/{all_stats['original_positives']} positives ({all_stats['retention_rate']*100:.1f}%)")
    log(f"  negatives: A={all_stats['negatives_a']} B={all_stats['negatives_b']} total={all_stats['negatives_total']}")
    log(f"  total rows: {all_stats['total_rows']} (ratio={all_stats['ratio']})")
    log(f"  per-fold positives: {all_stats['per_fold_positives']}")

    # Assign all_split rows to folds (by enzyme group)
    all_r2f = {}
    for i, r in enumerate(all_rows):
        all_r2f[i] = all_eg2f[e2g[r.enzyme]]
    all_parts = make_partitions(all_rows, all_r2f)
    all_issues = validate_strict("all_split", all_parts, e2g)
    log(f"  leakage: {'PASS' if not all_issues else 'FAIL '+str(all_issues)}")
    all_fs = write_split("all_split", all_rows, all_parts, odir)
    stats["splits"]["all_split"] = {
        "leakage": all_issues, "strict": True,
        "assignment": all_stats,
        "folds": {str(k): v for k, v in all_fs.items()},
    }

    write_json(cdir / "negatives" / "split_stats.json", stats)
    log(f"\n{'='*60}")
    log(f"[DONE] shared data.csv: {len(master)} rows (random/reaction/enzyme)")
    log(f"[DONE] all_split data.csv: {all_stats['total_rows']} rows (strict, {all_stats['retention_rate']*100:.1f}% retained)")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
