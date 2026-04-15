"""Sanity check all trusted/partial/align_failed enzymes."""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, sys
from pathlib import Path
from collections import defaultdict
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
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

data = torch.load(str(RESID_MAP_PT), weights_only=False)
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    enz_rows = list(csv.DictReader(f))

# enzyme -> docks from sidecars
enz_to_docks = defaultdict(list)
for sp in SIDECAR_PATHS:
    sc = torch.load(str(sp), weights_only=False)
    base_idx = torch.load(str(sp).replace("dock_sidecar.pt", "index.pt"), weights_only=False)
    dids = sc["dock_indices"].tolist()
    eids = base_idx["enzyme_ids"].tolist()
    for e, d in zip(eids, dids):
        enz_to_docks[int(e)].append(int(d))


def parse_pocket(pdb):
    seen = {}; order = []
    with open(pdb) as f:
        for line in f:
            if line[:6].strip() == "COMPND": break
            if line[:6].strip() != "ATOM": continue
            rn = line[17:20].strip()
            if rn not in THREE_TO_ONE: continue
            c = line[21]
            try: r = int(line[22:26])
            except: continue
            ic = line[26]
            k = (c, r, ic)
            if k not in seen:
                seen[k] = rn; order.append(k)
    return [(c, r, ic, seen[(c,r,ic)], THREE_TO_ONE[seen[(c,r,ic)]]) for (c,r,ic) in order]


# Find edge cases
edge_cases = []
# align_failed ones: enzymes present in sidecar but NOT in data (i.e., built failed)
for e in enz_to_docks:
    if e not in data:
        edge_cases.append((e, "align_failed"))
for e, d in data.items():
    if d["tier"] in ("trusted", "partial", "suspect"):
        edge_cases.append((e, d["tier"]))

print(f"edge cases to check: {len(edge_cases)}")

total_res = 0
total_mapped = 0
total_match = 0
for enz_id, label in edge_cases:
    print(f"\n--- enzyme {enz_id} [{label}] ---")
    if enz_id not in data:
        print(f"  NOT in resid_map (align_failed)")
        # Peek the source
        continue
    d = data[enz_id]
    print(f"  source={d['source']}, uniprot={d['uniprot']}")
    print(f"  identity={d['identity']:.4f}, coverage={d['coverage']:.4f}")
    print(f"  per_chain={d['per_chain']}")
    seq = enz_rows[enz_id]["Protein sequence"].strip()
    rmap = {}
    for k, v in d["resid_map"].items():
        c, r_, ic = k.split("|")
        rmap[(c, int(r_), ic)] = v

    # All this enzyme's docks, sum up stats
    enz_total_res = enz_total_mapped = enz_total_match = 0
    n_pockets = 0
    for dock in enz_to_docks.get(enz_id, [])[:5]:
        ppath = POCKET_DIR / f"{dock}.pdb"
        if not ppath.exists(): continue
        residues = parse_pocket(ppath)
        if not residues: continue
        n_pockets += 1
        for c, r, ic, rn, a1 in residues:
            enz_total_res += 1
            total_res += 1
            up = rmap.get((c, r, ic))
            if up is None:
                continue
            enz_total_mapped += 1
            total_mapped += 1
            if 0 <= up < len(seq) and seq[up] == a1:
                enz_total_match += 1
                total_match += 1
    cov = enz_total_mapped / max(enz_total_res, 1)
    rate = enz_total_match / max(enz_total_mapped, 1)
    print(f"  pockets={n_pockets}, residues={enz_total_res}, "
          f"mapped={enz_total_mapped} ({100*cov:.1f}%), "
          f"aa_match={enz_total_match}/{enz_total_mapped} ({100*rate:.1f}%)")

print(f"\n{'='*60}")
print(f"[ALL EDGE CASES] residues={total_res}, "
      f"mapped={total_mapped} ({100*total_mapped/max(total_res,1):.1f}%), "
      f"aa_match={total_match} ({100*total_match/max(total_mapped,1):.1f}%)")

# Find align_failed enzymes (from the enz_to_docks set minus data keys)
failed_enz = set(enz_to_docks.keys()) - set(data.keys())
print(f"\nalign_failed enzymes: {len(failed_enz)}")
for e in failed_enz:
    print(f"  enz {e}: docks={enz_to_docks[e][:3]}")

sys.stdout.flush()
'''


def main():
    import subprocess, sys
    print("[local driver] checking edge cases...")
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=300,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
