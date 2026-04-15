"""
Sanity check the built enzyme_resid_map.pt on pocket residues.

For each of 100 randomly sampled enzymes (from gold/trusted tiers):
  - Pick one of the enzyme's docks
  - Read its pocket PDB
  - For each pocket residue, use resid_map to get uniprot_pos
  - Check seq[uniprot_pos] == pocket_aa (tolerate mutations in trusted tier)

Also report:
  - Pocket coverage: how many pocket residues got a uniprot_pos
  - aa_match rate per enzyme
  - Any enzyme with < 100% pocket coverage OR < 95% aa_match -> flag
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, random, sys
from pathlib import Path
from collections import defaultdict, Counter
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
COMPLEX_DIR = BASE / "data/structure/complex"
OVERLAY = BASE / "data/pt_cache_dualgraph_allfix_unified"
RESID_MAP_PT = OVERLAY / "enzyme_resid_map.pt"
SIDECAR_PATHS = [
    OVERLAY / "random/train/dock_sidecar.pt",
    OVERLAY / "random/val/dock_sidecar.pt",
    OVERLAY / "random/test/dock_sidecar.pt",
]

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "MSE": "M",
}

# Load map
data = torch.load(str(RESID_MAP_PT), weights_only=False)
print(f"loaded resid_map for {len(data)} enzymes")

# Report the special tiers
tier_counts = Counter()
for e, d in data.items():
    tier_counts[d["tier"]] += 1
print(f"tiers: {dict(tier_counts)}")

# Show partial + align_failed
for e, d in data.items():
    if d["tier"] in ("partial", "suspect"):
        print(f"\n[{d['tier']}] enzyme {e}:")
        print(f"  source={d['source']}, uniprot={d['uniprot']}")
        print(f"  identity={d['identity']:.4f}, coverage={d['coverage']:.4f}")
        print(f"  per_chain={d['per_chain']}")
        print(f"  rejected={d['rejected_chains']}")

# Load enzyme sequences
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    enz_rows = list(csv.DictReader(f))

# Build enzyme -> list of docks from sidecars
enz_to_docks = defaultdict(list)
for sp in SIDECAR_PATHS:
    sc = torch.load(str(sp), weights_only=False)
    base_idx = torch.load(str(sp).replace("dock_sidecar.pt", "index.pt"), weights_only=False)
    dids = sc["dock_indices"].tolist()
    eids = base_idx["enzyme_ids"].tolist()
    for e, d in zip(eids, dids):
        enz_to_docks[int(e)].append(int(d))


def parse_pocket_residues(pdb_path):
    """Return list of (chain, resid, icode, resname, aa1). Same parser as before."""
    seen = {}
    order = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "COMPND":
                break
            if rec != "ATOM":
                continue
            resname = line[17:20].strip()
            if resname not in THREE_TO_ONE:
                continue
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            key = (chain, resseq, icode)
            if key not in seen:
                seen[key] = resname
                order.append(key)
    result = []
    for key in order:
        c, r, ic = key
        rn = seen[key]
        result.append((c, r, ic, rn, THREE_TO_ONE[rn]))
    return result


# ---- Random 100-enzyme sanity ----
random.seed(42)
candidates = [e for e, d in data.items() if d["tier"] in ("gold", "trusted") and e in enz_to_docks]
print(f"\neligible enzymes (gold/trusted with docks): {len(candidates)}")
sample = random.sample(candidates, min(100, len(candidates)))

print(f"\n{'='*70}")
print(f"Pocket-level sanity: {len(sample)} random enzymes, 3 random pockets each")
print(f"{'='*70}")

total_pockets = 0
total_pocket_res = 0
total_mapped = 0
total_aa_match = 0
enzyme_aa_match_rates = []
per_tier_rates = defaultdict(list)
mutant_count = 0
bad_coverage = []

for enz_id in sample:
    d = data[enz_id]
    seq = enz_rows[enz_id]["Protein sequence"].strip()
    resid_map = {}
    for k, v in d["resid_map"].items():
        c, r_, ic = k.split("|")
        resid_map[(c, int(r_), ic)] = v

    docks = enz_to_docks[enz_id][:3]  # up to 3 pockets
    for dock in docks:
        ppath = POCKET_DIR / f"{dock}.pdb"
        if not ppath.exists():
            continue
        residues = parse_pocket_residues(ppath)
        if not residues:
            continue
        total_pockets += 1
        total_pocket_res += len(residues)
        n_mapped = 0
        n_match = 0
        for c, r, ic, rn, a1 in residues:
            up = resid_map.get((c, r, ic))
            if up is None:
                continue
            n_mapped += 1
            if 0 <= up < len(seq) and seq[up] == a1:
                n_match += 1
        total_mapped += n_mapped
        total_aa_match += n_match
        coverage = n_mapped / len(residues) if residues else 0
        aa_rate = n_match / n_mapped if n_mapped else 0
        enzyme_aa_match_rates.append((enz_id, d["tier"], dock, len(residues), n_mapped, n_match, aa_rate))
        per_tier_rates[d["tier"]].append(aa_rate)
        if d["tier"] == "trusted" and aa_rate < 1.0:
            mutant_count += 1
        if coverage < 1.0:
            bad_coverage.append((enz_id, d["tier"], dock, len(residues), n_mapped, coverage))

# Aggregate
print(f"\n[AGGREGATE]")
print(f"  sampled enzymes:    {len(sample)}")
print(f"  pockets parsed:     {total_pockets}")
print(f"  pocket residues:    {total_pocket_res}")
print(f"  mapped (coverage):  {total_mapped}/{total_pocket_res} ({100*total_mapped/max(total_pocket_res,1):.2f}%)")
print(f"  aa_match:           {total_aa_match}/{max(total_mapped,1)} ({100*total_aa_match/max(total_mapped,1):.2f}%)")

print(f"\n[PER TIER]")
for tier in ["gold", "trusted"]:
    rates = per_tier_rates[tier]
    if rates:
        avg = sum(rates) / len(rates)
        below_100 = sum(1 for r in rates if r < 1.0)
        below_95 = sum(1 for r in rates if r < 0.95)
        print(f"  {tier}: n={len(rates)}, avg_aa_rate={avg:.4f}, <100%={below_100}, <95%={below_95}")

if bad_coverage:
    print(f"\n[LOW COVERAGE] (pocket coverage < 100%):")
    for row in bad_coverage[:10]:
        print(f"  {row}")
else:
    print(f"\n[COVERAGE] all 100% pocket residues mapped")

# Show any enzyme with aa_rate < 95% (potentially bad)
bad_ones = [r for r in enzyme_aa_match_rates if r[6] < 0.95]
if bad_ones:
    print(f"\n[AA MATCH < 95%]:")
    for row in bad_ones[:10]:
        print(f"  enz={row[0]} tier={row[1]} dock={row[2]} N={row[3]} mapped={row[4]} match={row[5]} rate={row[6]:.3f}")
else:
    print(f"\n[AA MATCH] all >= 95%")

sys.stdout.flush()
'''


def main():
    import subprocess, sys
    print("[local driver] running sanity on server...")
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=600,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
