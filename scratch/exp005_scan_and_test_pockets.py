"""
EXP005 Step 3ab: scan all 50177 pocket PDBs for edge cases + run GVP
feature extraction on 10 selected pockets (min-N, max-N, random, edge cases).

Phase 1 (fast, ~1-2 min even on 0.5 vCPU): pure Python text parse of all
pocket PDBs. Counts:
  - chain_count (multi-chain?)
  - insertion code count
  - hetero count
  - missing backbone count
  - altloc-backbone count
  - kept residue count (N)
  - raw residue count

Output: global stats + list of problematic dock_ids + N distribution.

Phase 2: Select 10 pockets = [min N, max N, 5 random, 1 short seq, 1 long seq,
1 edge-case pocket if any]. Run full GVP feature extraction + penetration
test (verify pocket_residue_idx matches enzyme sequence amino acids).

Entirely on server via ssh stdin. No persistent files, no local src/ touches.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
import math
import random
import sys
import time
from pathlib import Path
from collections import Counter

import numpy as np
import torch
import torch.nn.functional as F
from Bio.PDB import PDBParser

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
ENZ_CSV = BASE / "data/Enzymes.csv"

# Load enzyme sequences for penetration test
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    ENZ_ROWS = list(csv.DictReader(f))

# Load any split CSV to get (dock -> enzyme_id) map (random train has largest)
SPLIT_CSVS = [
    BASE / "data/splits/random/training_datas_0_pt.csv",
    BASE / "data/splits/random/val_datas_0_pt.csv",
    BASE / "data/splits/random/testing_datas_0_pt.csv",
]
DOCK_TO_ENZYME = {}
for p in SPLIT_CSVS:
    with open(p, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            d = int(row["Dock Index"])
            e = int(row["Enzyme Index"])
            DOCK_TO_ENZYME[d] = e

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

# ---------------------------------------------------------------------------
# Phase 1: fast metadata scan using raw text parsing (no BioPython)
# ---------------------------------------------------------------------------
BACKBONE = {"N", "CA", "C", "O"}

def scan_pocket_metadata(pdb_path):
    """
    Fast text-based scan. Returns dict with:
      chains: set of chain IDs
      kept_residues: count of standard AA residues with full backbone
      total_residues: count of all unique (chain, resid, icode) residues
      insertion_code_count: residues with non-blank icode
      hetero_count: records with hetero flag
      missing_backbone_count: residues missing any of N/CA/C/O
      altloc_backbone_count: backbone atoms with non-blank altloc
      non_standard_resname_count: residues with non-standard resname
    """
    chains = set()
    residue_atoms = {}  # (chain, resseq, icode) -> {atom_name: altloc}
    icode_count = 0
    hetero_count = 0
    altloc_bb = 0

    with open(pdb_path) as f:
        for line in f:
            if not line.startswith("ATOM") and not line.startswith("HETATM"):
                continue
            if line.startswith("HETATM"):
                hetero_count += 1
                continue
            atom_name = line[12:16].strip()
            altloc = line[16]
            resname = line[17:20].strip()
            chain_id = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            chains.add(chain_id)
            if icode != " ":
                icode_count += 1
            key = (chain_id, resseq, icode, resname)
            if key not in residue_atoms:
                residue_atoms[key] = {}
            if atom_name in BACKBONE:
                residue_atoms[key][atom_name] = altloc
                if altloc != " ":
                    altloc_bb += 1

    total_residues = len(residue_atoms)
    kept_residues = 0
    missing_bb = 0
    non_std = 0
    for (cid, rseq, ico, rname), atoms in residue_atoms.items():
        if rname not in THREE_TO_ONE:
            non_std += 1
            continue
        if ico != " ":
            continue  # already counted in icode_count
        if not BACKBONE.issubset(atoms.keys()):
            missing_bb += 1
            continue
        kept_residues += 1

    return {
        "chains": chains,
        "total_residues": total_residues,
        "kept_residues": kept_residues,
        "insertion_code_count": icode_count,
        "hetero_count": hetero_count,
        "missing_backbone_count": missing_bb,
        "altloc_backbone_count": altloc_bb,
        "non_standard_resname_count": non_std,
    }


def phase1_scan_all():
    pdb_files = sorted(POCKET_DIR.glob("*.pdb"), key=lambda p: int(p.stem))
    total = len(pdb_files)
    print(f"Phase 1: scanning {total} pocket PDBs...", flush=True)

    stats_per_dock = {}  # dock_index -> stats
    agg = Counter()
    n_dist = []
    problematic_docks = {"multi_chain": [], "icode": [], "hetero": [], "missing_bb": [], "altloc_bb": [], "non_std": [], "empty_kept": []}

    t0 = time.time()
    for i, pdb in enumerate(pdb_files):
        dock = int(pdb.stem)
        try:
            st = scan_pocket_metadata(pdb)
        except Exception as e:
            print(f"  scan failed for {dock}: {e}", flush=True)
            continue
        stats_per_dock[dock] = st
        n_dist.append(st["kept_residues"])

        if len(st["chains"]) > 1:
            problematic_docks["multi_chain"].append(dock)
        if st["insertion_code_count"] > 0:
            problematic_docks["icode"].append(dock)
        if st["hetero_count"] > 0:
            problematic_docks["hetero"].append(dock)
        if st["missing_backbone_count"] > 0:
            problematic_docks["missing_bb"].append(dock)
        if st["altloc_backbone_count"] > 0:
            problematic_docks["altloc_bb"].append(dock)
        if st["non_standard_resname_count"] > 0:
            problematic_docks["non_std"].append(dock)
        if st["kept_residues"] == 0:
            problematic_docks["empty_kept"].append(dock)

        agg["multi_chain"] += (len(st["chains"]) > 1)
        agg["icode_any"] += (st["insertion_code_count"] > 0)
        agg["hetero_any"] += (st["hetero_count"] > 0)
        agg["missing_bb_any"] += (st["missing_backbone_count"] > 0)
        agg["altloc_bb_any"] += (st["altloc_backbone_count"] > 0)
        agg["non_std_any"] += (st["non_standard_resname_count"] > 0)
        agg["empty_kept"] += (st["kept_residues"] == 0)

        if (i + 1) % 10000 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1}/{total}  ({elapsed:.0f}s, {(i+1)/elapsed:.0f} pdbs/s)", flush=True)

    elapsed = time.time() - t0
    n_arr = np.array(n_dist)
    print(f"\nPhase 1 done in {elapsed:.0f}s.")
    print(f"Scanned: {len(stats_per_dock)}/{total}")
    print()
    print("=== aggregate edge-case counts (any occurrence in a pocket) ===")
    for k, v in sorted(agg.items()):
        print(f"  {k:20s}: {v}")
    print()
    print("=== N distribution (kept residues per pocket) ===")
    print(f"  min:    {n_arr.min()}")
    print(f"  max:    {n_arr.max()}")
    print(f"  mean:   {n_arr.mean():.1f}")
    print(f"  median: {int(np.median(n_arr))}")
    print(f"  p5:     {int(np.percentile(n_arr, 5))}")
    print(f"  p95:    {int(np.percentile(n_arr, 95))}")
    print()
    print("=== problematic dock counts (first 10 each) ===")
    for cat, lst in problematic_docks.items():
        print(f"  {cat:20s}: {len(lst):6d}  first10: {lst[:10]}")

    return stats_per_dock, problematic_docks, n_arr


# ---------------------------------------------------------------------------
# Phase 2: GVP featurization (reuse from single-pocket script)
# ---------------------------------------------------------------------------
def _normalize(t, dim=-1):
    return torch.nan_to_num(torch.div(t, torch.norm(t, dim=dim, keepdim=True)))

def _rbf(D, D_min=0.0, D_max=20.0, D_count=16):
    D_mu = torch.linspace(D_min, D_max, D_count).view([1, -1])
    D_sigma = (D_max - D_min) / D_count
    D_expand = torch.unsqueeze(D, -1)
    return torch.exp(-(((D_expand - D_mu) / D_sigma) ** 2))

def _dihedrals(X, eps=1e-7):
    X = torch.reshape(X[:, :3], [3 * X.shape[0], 3])
    dX = X[1:] - X[:-1]
    U = _normalize(dX, dim=-1)
    u_2, u_1, u_0 = U[:-2], U[1:-1], U[2:]
    n_2 = _normalize(torch.cross(u_2, u_1, dim=-1), dim=-1)
    n_1 = _normalize(torch.cross(u_1, u_0, dim=-1), dim=-1)
    cosD = torch.clamp(torch.sum(n_2 * n_1, -1), -1 + eps, 1 - eps)
    D = torch.sign(torch.sum(u_2 * n_1, -1)) * torch.acos(cosD)
    D = F.pad(D, [1, 2])
    D = torch.reshape(D, [-1, 3])
    return torch.cat([torch.cos(D), torch.sin(D)], 1)

def _orientations(X_ca):
    forward = _normalize(X_ca[1:] - X_ca[:-1])
    backward = _normalize(X_ca[:-1] - X_ca[1:])
    forward = F.pad(forward, [0, 0, 0, 1])
    backward = F.pad(backward, [0, 0, 1, 0])
    return torch.cat([forward.unsqueeze(-2), backward.unsqueeze(-2)], -2)

def _sidechains(X):
    n, origin, c = X[:, 0], X[:, 1], X[:, 2]
    c, n = _normalize(c - origin), _normalize(n - origin)
    bisector = _normalize(c + n)
    perp = _normalize(torch.cross(c, n, dim=-1))
    return -bisector * math.sqrt(1 / 3) - perp * math.sqrt(2 / 3)

def _positional_embeddings(edge_index, num_embeddings=16):
    d = edge_index[0] - edge_index[1]
    frequency = torch.exp(
        torch.arange(0, num_embeddings, 2, dtype=torch.float32)
        * -(np.log(10000.0) / num_embeddings))
    angles = d.unsqueeze(-1) * frequency
    return torch.cat((torch.cos(angles), torch.sin(angles)), -1)

def knn_graph_custom(X_ca, k):
    N = X_ca.shape[0]
    if N < 2:
        return torch.zeros(2, 0, dtype=torch.long)
    dist = torch.cdist(X_ca, X_ca, p=2)
    dist.fill_diagonal_(float("inf"))
    k_eff = min(k, N - 1)
    _, nbr_idx = torch.topk(dist, k_eff, dim=1, largest=False, sorted=True)
    src = torch.arange(N).unsqueeze(1).expand(-1, k_eff).reshape(-1)
    dst = nbr_idx.reshape(-1)
    return torch.stack([src, dst], dim=0)


class PocketError(Exception):
    pass


def extract_pocket_residues_strict(pocket_pdb_path):
    """Strict BioPython parse. Raises on edge cases."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", str(pocket_pdb_path))
    models = list(structure.get_models())
    if len(models) != 1:
        raise PocketError(f"expected 1 model, got {len(models)}")
    chains = list(models[0].get_chains())
    if len(chains) != 1:
        raise PocketError(f"expected 1 chain, got {len(chains)}")

    out = []
    for res in chains[0].get_residues():
        hetero, resseq, icode = res.id
        if hetero.strip() != "":
            continue
        if icode.strip() != "":
            raise PocketError(f"insertion code on resid {resseq}{icode}")
        if res.resname not in THREE_TO_ONE:
            continue
        if resseq <= 0:
            raise PocketError(f"non-positive resid {resseq}")
        if not all(a in res for a in ("N", "CA", "C", "O")):
            continue
        N_c = np.array(res["N"].coord, dtype=np.float32)
        CA_c = np.array(res["CA"].coord, dtype=np.float32)
        C_c = np.array(res["C"].coord, dtype=np.float32)
        O_c = np.array(res["O"].coord, dtype=np.float32)
        out.append((int(resseq), res.resname, N_c, CA_c, C_c, O_c))
    if not out:
        raise PocketError("no valid residues")
    return out


def build_gvp_features(pocket_pdb_path, k=30):
    residues = extract_pocket_residues_strict(pocket_pdb_path)
    N = len(residues)
    coords_list = [np.stack([r[2], r[3], r[4], r[5]]) for r in residues]
    pocket_residue_idx = [r[0] - 1 for r in residues]
    coords = torch.from_numpy(np.stack(coords_list)).float()
    X_ca = coords[:, 1]
    if not torch.isfinite(coords).all():
        raise PocketError("non-finite coords")
    node_s = _dihedrals(coords)
    orientations = _orientations(X_ca)
    sidechains = _sidechains(coords)
    node_v = torch.cat([orientations, sidechains.unsqueeze(-2)], dim=-2)
    edge_index = knn_graph_custom(X_ca, k=k)
    if edge_index.shape[1] == 0:
        raise PocketError(f"empty edges (N={N})")
    E_vectors = X_ca[edge_index[0]] - X_ca[edge_index[1]]
    rbf = _rbf(E_vectors.norm(dim=-1))
    pos_emb = _positional_embeddings(edge_index)
    edge_s = torch.cat([rbf, pos_emb], dim=-1)
    edge_v = _normalize(E_vectors).unsqueeze(-2)
    node_s, node_v, edge_s, edge_v = map(torch.nan_to_num, (node_s, node_v, edge_s, edge_v))
    return {
        "node_s": node_s.float(),
        "node_v": node_v.float(),
        "edge_index": edge_index.long(),
        "edge_s": edge_s.float(),
        "edge_v": edge_v.float(),
        "pocket_residue_idx": torch.tensor(pocket_residue_idx, dtype=torch.long),
        "residue_meta": residues,  # for verification
    }


def run_penetration_test(dock, feats):
    """Verify pocket_residue_idx matches enzyme sequence amino acids."""
    if dock not in DOCK_TO_ENZYME:
        return None  # dock not in any split
    enzyme_id = DOCK_TO_ENZYME[dock]
    seq = ENZ_ROWS[enzyme_id]["Protein sequence"].strip()
    seq_len = len(seq)
    N = feats["node_s"].shape[0]
    max_idx = feats["pocket_residue_idx"].max().item()
    if max_idx >= seq_len:
        return {"enzyme_id": enzyme_id, "seq_len": seq_len, "N": N,
                "max_idx": max_idx, "error": "max_idx >= seq_len"}
    n_match = 0
    for i, (resid, resname, *_) in enumerate(feats["residue_meta"]):
        expected = THREE_TO_ONE[resname]
        actual = seq[resid - 1]
        if expected == actual:
            n_match += 1
    return {"enzyme_id": enzyme_id, "seq_len": seq_len, "N": N,
            "max_idx": max_idx, "aa_match": f"{n_match}/{N}",
            "ok": n_match == N}


def phase2_test_10(stats, problematic, n_arr):
    """Select 10 pockets and run full GVP extraction + verification."""
    pdb_files = sorted(POCKET_DIR.glob("*.pdb"), key=lambda p: int(p.stem))
    dock_list = [int(p.stem) for p in pdb_files]

    # Find min N and max N pockets
    valid_docks = [d for d, s in stats.items() if s["kept_residues"] > 0 and len(s["chains"]) == 1]
    min_dock = min(valid_docks, key=lambda d: stats[d]["kept_residues"])
    max_dock = max(valid_docks, key=lambda d: stats[d]["kept_residues"])

    # Random 5 (seed fixed for reproducibility)
    random.seed(1234)
    rand_docks = random.sample(valid_docks, 5)

    # Short / long seq candidates
    short_seq_docks = []
    long_seq_docks = []
    for d in valid_docks:
        if d in DOCK_TO_ENZYME:
            enz = DOCK_TO_ENZYME[d]
            sl = len(ENZ_ROWS[enz]["Protein sequence"].strip())
            if sl < 300:
                short_seq_docks.append(d)
            elif sl > 800:
                long_seq_docks.append(d)
    short_dock = short_seq_docks[0] if short_seq_docks else None
    long_dock = long_seq_docks[0] if long_seq_docks else None

    # Edge case pocket (first insertion code if any; else first multi-chain etc.)
    edge_dock = None
    for cat in ("icode", "multi_chain", "hetero", "missing_bb", "altloc_bb"):
        if problematic[cat]:
            edge_dock = problematic[cat][0]
            break

    selected = {"min_N": min_dock, "max_N": max_dock}
    for i, d in enumerate(rand_docks):
        selected[f"rand_{i}"] = d
    if short_dock is not None:
        selected["short_seq"] = short_dock
    if long_dock is not None:
        selected["long_seq"] = long_dock
    if edge_dock is not None:
        selected["edge_case"] = edge_dock

    print()
    print("=" * 60)
    print(f"Phase 2: running full GVP extraction on {len(selected)} selected pockets")
    print("=" * 60)

    total_pass = 0
    total_fail = 0
    for label, dock in selected.items():
        pdb = POCKET_DIR / f"{dock}.pdb"
        print(f"\n--- {label}: dock={dock} ---")
        try:
            feats = build_gvp_features(pdb)
            N = feats["node_s"].shape[0]
            E = feats["edge_index"].shape[1]
            expected_E = N * min(30, N - 1)
            n_match = sum(
                THREE_TO_ONE[r[1]] == ENZ_ROWS[DOCK_TO_ENZYME.get(dock, -1)]["Protein sequence"].strip()[r[0]-1]
                for r in feats["residue_meta"]
                if dock in DOCK_TO_ENZYME and 0 <= r[0]-1 < len(ENZ_ROWS[DOCK_TO_ENZYME[dock]]["Protein sequence"].strip())
            ) if dock in DOCK_TO_ENZYME else -1
            degree_min = feats["edge_index"][0].bincount(minlength=N).min().item()
            degree_max = feats["edge_index"][0].bincount(minlength=N).max().item()
            print(f"  N={N}, E={E} (expected {expected_E}), "
                  f"degree=[{degree_min},{degree_max}], "
                  f"aa_match={n_match}/{N}")
            ok = (E == expected_E) and (n_match == N if dock in DOCK_TO_ENZYME else True)
            if ok:
                print(f"  [OK]")
                total_pass += 1
            else:
                print(f"  [FAIL]")
                total_fail += 1
        except PocketError as e:
            print(f"  [PARSE_ERR]: {e}  (expected for edge cases with strict parser)")
            total_fail += 1
        except Exception as e:
            print(f"  [UNEXPECTED]: {e}")
            total_fail += 1

    print()
    print("=" * 60)
    print(f"Phase 2 result: {total_pass}/{total_pass+total_fail} passed")
    print("=" * 60)


def main():
    stats, problematic, n_arr = phase1_scan_all()
    phase2_test_10(stats, problematic, n_arr)

if __name__ == "__main__":
    main()
'''


def main():
    import subprocess, sys
    print("[local driver] running scan + 10-pocket GVP test on server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=1200,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
