"""
EXP005 offset diagnostic: run global offset detection on all 50177 pockets.

For each pocket:
  1. Parse kept residues -> list of (resid, one_letter_aa)
  2. Get associated enzyme sequence from Enzymes.csv via DOCK_TO_ENZYME
  3. Score delta=0
  4. Search best constant offset delta in [-max_range, +max_range]
  5. If offset score still poor, run substring match (ignore resid, slide pocket
     amino-acid string over enzyme seq)

Classify each pocket into buckets:
  - direct_ok       : score(delta=0) >= 0.95
  - offset_ok       : best_offset_score >= 0.95 (and different from delta=0)
  - substring_ok    : substring match score >= 0.95 (numbering issue not simple offset)
  - suspect         : all scores < 0.6 (probably wrong join / different sequence)
  - weird           : empty pocket, no valid residues, etc.

Output global stats + per-bucket examples + per-enzyme offset modal analysis.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
import sys
import time
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
ENZ_CSV = BASE / "data/Enzymes.csv"
SPLIT_CSVS = [
    BASE / "data/splits/random/training_datas_0_pt.csv",
    BASE / "data/splits/random/val_datas_0_pt.csv",
    BASE / "data/splits/random/testing_datas_0_pt.csv",
]

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

# Load enzyme sequences
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    ENZ_ROWS = list(csv.DictReader(f))
ENZ_SEQ = [row["Protein sequence"].strip() for row in ENZ_ROWS]
ENZ_UNI = [row["uniprots"] for row in ENZ_ROWS]

# Build DOCK_TO_ENZYME
DOCK_TO_ENZYME = {}
for p in SPLIT_CSVS:
    with open(p, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            DOCK_TO_ENZYME[int(row["Dock Index"])] = int(row["Enzyme Index"])

BACKBONE = {"N", "CA", "C", "O"}


def parse_pocket_residues(pdb_path):
    """Fast text parse. Return list of (resid, one_letter_aa) for kept residues."""
    residue_atoms = {}  # (chain, resseq, resname) -> set of atom names
    chains = set()
    for line in open(pdb_path):
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21]
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        icode = line[26]
        if icode != " ":
            continue  # skip insertion codes
        if resname not in THREE_TO_ONE:
            continue
        if resseq <= 0:
            continue
        chains.add(chain)
        key = (chain, resseq, resname)
        residue_atoms.setdefault(key, set()).add(atom_name)

    if len(chains) > 1:
        return None, "multi_chain"

    # Keep residues with full backbone
    kept = []
    for (cid, rseq, rname), atoms in residue_atoms.items():
        if BACKBONE.issubset(atoms):
            kept.append((rseq, THREE_TO_ONE[rname]))
    if not kept:
        return None, "empty"

    # Sort by resid for deterministic substring match
    kept.sort()
    return kept, None


def score_offset(kept, seq, delta):
    """Return (matches, valid_mapped) for given offset."""
    seq_len = len(seq)
    matches = 0
    valid = 0
    for resid, aa in kept:
        idx = resid - 1 + delta
        if 0 <= idx < seq_len:
            valid += 1
            if seq[idx] == aa:
                matches += 1
    if valid == 0:
        return 0.0, 0
    return matches / valid, valid


def best_constant_offset(kept, seq, delta_range=200):
    """Scan offset delta in [-delta_range, +delta_range], return best_delta, best_score, best_valid."""
    best_score = -1.0
    best_delta = 0
    best_valid = 0
    for delta in range(-delta_range, delta_range + 1):
        score, valid = score_offset(kept, seq, delta)
        if valid == 0:
            continue
        if score > best_score or (score == best_score and valid > best_valid):
            best_score = score
            best_delta = delta
            best_valid = valid
    return best_delta, best_score, best_valid


def substring_match(kept, seq):
    """
    Ignore resid. Build pocket AA string in residue order, slide over enzyme
    seq, return best fraction of matches.

    Note: pocket residues may not be contiguous in original sequence. This
    check assumes they might form a contiguous block in sequence after
    offset. For non-contiguous pockets, this will underestimate.

    Better approach: use pocket AAs as a set/histogram match. But for offset
    cases the residues ARE contiguous in original sequence.
    """
    pocket_aa = "".join(aa for _, aa in kept)
    seq_len = len(seq)
    n = len(pocket_aa)
    if n == 0 or n > seq_len:
        return 0.0
    best = 0
    for start in range(seq_len - n + 1):
        match = sum(1 for i in range(n) if seq[start + i] == pocket_aa[i])
        if match > best:
            best = match
    return best / n


def classify(direct_score, best_off_score, substring_score):
    if direct_score >= 0.95:
        return "direct_ok"
    if best_off_score >= 0.95:
        return "offset_ok"
    if substring_score >= 0.95:
        return "substring_ok"
    if max(direct_score, best_off_score, substring_score) < 0.6:
        return "suspect"
    return "partial"  # 0.6-0.95 zone


def main():
    pdb_files = sorted(POCKET_DIR.glob("*.pdb"), key=lambda p: int(p.stem))
    print(f"Diagnosing {len(pdb_files)} pockets...", flush=True)

    bucket_counts = Counter()
    per_enzyme_deltas = defaultdict(list)  # enzyme_id -> list of best_deltas
    failure_examples = defaultdict(list)
    global_deltas = Counter()
    direct_scores = []
    offset_scores = []
    orphan_no_enzyme = 0
    parse_errors = 0

    t0 = time.time()
    for i, pdb in enumerate(pdb_files):
        dock = int(pdb.stem)
        enzyme_id = DOCK_TO_ENZYME.get(dock, -1)
        if enzyme_id < 0:
            orphan_no_enzyme += 1
            continue

        kept, err = parse_pocket_residues(pdb)
        if kept is None:
            parse_errors += 1
            bucket_counts[f"parse_{err}"] += 1
            continue

        seq = ENZ_SEQ[enzyme_id]
        direct_score, direct_valid = score_offset(kept, seq, 0)
        best_delta, best_score, _ = best_constant_offset(kept, seq, delta_range=200)

        # Substring match only if offset is poor
        if best_score >= 0.95:
            substring_score = 0.0  # not needed
        else:
            substring_score = substring_match(kept, seq)

        bucket = classify(direct_score, best_score, substring_score)
        bucket_counts[bucket] += 1
        direct_scores.append(direct_score)
        offset_scores.append(best_score)
        global_deltas[best_delta] += 1
        per_enzyme_deltas[enzyme_id].append(best_delta)

        # Save examples
        if len(failure_examples[bucket]) < 5:
            failure_examples[bucket].append({
                "dock": dock,
                "enz_id": enzyme_id,
                "uniprot": ENZ_UNI[enzyme_id],
                "seq_len": len(seq),
                "N": len(kept),
                "direct_score": round(direct_score, 3),
                "best_delta": best_delta,
                "best_score": round(best_score, 3),
                "substring_score": round(substring_score, 3),
            })

        if (i + 1) % 10000 == 0:
            print(f"  ... {i+1}/{len(pdb_files)}  ({time.time()-t0:.0f}s)", flush=True)

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")
    print()

    # Global stats
    print("=" * 60)
    print("BUCKET COUNTS")
    print("=" * 60)
    total_analyzed = sum(bucket_counts.values())
    for b, c in bucket_counts.most_common():
        print(f"  {b:20s}: {c:7d}  ({100*c/total_analyzed:.1f}%)")
    print(f"\n  orphan_no_enzyme (not in any split): {orphan_no_enzyme}")
    print(f"  parse_errors:                        {parse_errors}")

    print()
    print("=" * 60)
    print("DIRECT (delta=0) SCORE DISTRIBUTION")
    print("=" * 60)
    d_arr = np.array(direct_scores)
    print(f"  mean:   {d_arr.mean():.3f}")
    print(f"  median: {np.median(d_arr):.3f}")
    print(f"  >=0.95: {(d_arr >= 0.95).sum()} / {len(d_arr)} ({100*(d_arr>=0.95).mean():.1f}%)")
    print(f"  <0.5:   {(d_arr < 0.5).sum()} / {len(d_arr)} ({100*(d_arr<0.5).mean():.1f}%)")

    print()
    print("=" * 60)
    print("BEST OFFSET SCORE DISTRIBUTION")
    print("=" * 60)
    o_arr = np.array(offset_scores)
    print(f"  mean:   {o_arr.mean():.3f}")
    print(f"  median: {np.median(o_arr):.3f}")
    print(f"  >=0.95: {(o_arr >= 0.95).sum()} / {len(o_arr)} ({100*(o_arr>=0.95).mean():.1f}%)")

    print()
    print("=" * 60)
    print("DELTA DISTRIBUTION (top 20)")
    print("=" * 60)
    for d, c in global_deltas.most_common(20):
        print(f"  delta={d:+5d}: {c:6d}  ({100*c/total_analyzed:.1f}%)")

    print()
    print("=" * 60)
    print("PER-ENZYME OFFSET STABILITY")
    print("=" * 60)
    per_enzyme_stability = Counter()
    unstable_enzymes = []
    for enz_id, deltas in per_enzyme_deltas.items():
        if len(deltas) == 0:
            continue
        unique = set(deltas)
        if len(unique) == 1:
            per_enzyme_stability["stable"] += 1
        else:
            per_enzyme_stability["unstable"] += 1
            if len(unstable_enzymes) < 10:
                unstable_enzymes.append({
                    "enz_id": enz_id,
                    "n_docks": len(deltas),
                    "unique_deltas": sorted(unique)[:10],
                })
    print(f"  stable per-enzyme offset:  {per_enzyme_stability['stable']}")
    print(f"  unstable (varies by dock): {per_enzyme_stability['unstable']}")
    if unstable_enzymes:
        print(f"  first 10 unstable enzymes:")
        for e in unstable_enzymes:
            print(f"    {e}")

    print()
    print("=" * 60)
    print("BUCKET EXAMPLES (up to 5 each)")
    print("=" * 60)
    for bucket, examples in failure_examples.items():
        print(f"\n  {bucket}:")
        for e in examples:
            print(f"    {e}")

    sys.stdout.flush()


if __name__ == "__main__":
    main()
'''


def main():
    import subprocess, sys
    print("[local driver] running offset diagnostic on server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=3600,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
