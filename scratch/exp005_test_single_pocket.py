"""
EXP005 single-pocket GVP feature extraction — smoke test.

Runs entirely on server via ssh stdin. Local footprint = this file only.

Purpose:
  1. Read ONE pocket PDB file
  2. Parse residues, handle edge cases (insertion code / altloc / multi-chain)
  3. Build GVP features (node_s, node_v, edge_s, edge_v, edge_index, pocket_residue_idx)
  4. Use custom kNN (no torch_cluster dependency)
  5. Print shapes, dtypes, value ranges, and pocket_residue_idx vs csv sequence check

If this passes on one pocket, we can scale to batch processing 50177 pockets.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
import math
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from Bio.PDB import PDBParser

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
TEST_DOCK_INDEX = 3  # dock=3 → enzyme_id=93 (from Step 1 verification)
K_KNN = 30
NUM_RBF = 16
NUM_POSITIONAL = 16

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
ENZ_CSV = BASE / "data/Enzymes.csv"

THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}

# ---------------------------------------------------------------------------
# GVP math (ported from EnzymeCAGE gvp/data.py)
# ---------------------------------------------------------------------------
def _normalize(tensor, dim=-1):
    return torch.nan_to_num(
        torch.div(tensor, torch.norm(tensor, dim=dim, keepdim=True))
    )


def _rbf(D, D_min=0.0, D_max=20.0, D_count=16, device="cpu"):
    D_mu = torch.linspace(D_min, D_max, D_count, device=device).view([1, -1])
    D_sigma = (D_max - D_min) / D_count
    D_expand = torch.unsqueeze(D, -1)
    return torch.exp(-(((D_expand - D_mu) / D_sigma) ** 2))


def _dihedrals(X, eps=1e-7):
    # X: [N, 3, 3] — (N, CA, C) per residue
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
    return torch.cat([torch.cos(D), torch.sin(D)], 1)  # [N, 6]


def _orientations(X_ca):
    forward = _normalize(X_ca[1:] - X_ca[:-1])
    backward = _normalize(X_ca[:-1] - X_ca[1:])
    forward = F.pad(forward, [0, 0, 0, 1])
    backward = F.pad(backward, [0, 0, 1, 0])
    return torch.cat([forward.unsqueeze(-2), backward.unsqueeze(-2)], -2)  # [N, 2, 3]


def _sidechains(X):
    # X: [N, 4, 3]
    n, origin, c = X[:, 0], X[:, 1], X[:, 2]
    c, n = _normalize(c - origin), _normalize(n - origin)
    bisector = _normalize(c + n)
    perp = _normalize(torch.cross(c, n, dim=-1))
    return -bisector * math.sqrt(1 / 3) - perp * math.sqrt(2 / 3)  # [N, 3]


def _positional_embeddings(edge_index, num_embeddings=NUM_POSITIONAL, device="cpu"):
    d = edge_index[0] - edge_index[1]
    frequency = torch.exp(
        torch.arange(0, num_embeddings, 2, dtype=torch.float32, device=device)
        * -(np.log(10000.0) / num_embeddings)
    )
    angles = d.unsqueeze(-1) * frequency
    return torch.cat((torch.cos(angles), torch.sin(angles)), -1)  # [E, num_embeddings]


# ---------------------------------------------------------------------------
# Custom kNN (replaces torch_cluster.knn_graph)
# ---------------------------------------------------------------------------
def knn_graph_custom(X_ca, k):
    """
    Returns edge_index [2, E] with E = min(N-1, k) * N.
    For each node, connect to k nearest neighbors (excluding self).
    If N-1 < k, connect to all other nodes.
    Deterministic: stable sort by distance, then by index.
    """
    N = X_ca.shape[0]
    if N < 2:
        return torch.zeros(2, 0, dtype=torch.long, device=X_ca.device)

    # Pairwise squared distance
    dist = torch.cdist(X_ca, X_ca, p=2)  # [N, N]
    # Set self-distance to +inf so topk excludes self
    dist.fill_diagonal_(float("inf"))

    k_eff = min(k, N - 1)
    # topk with smallest distances
    _, nbr_idx = torch.topk(dist, k_eff, dim=1, largest=False, sorted=True)  # [N, k_eff]

    # Build edge_index: for each i, edges (i, j) where j is neighbor
    src = torch.arange(N, device=X_ca.device).unsqueeze(1).expand(-1, k_eff).reshape(-1)
    dst = nbr_idx.reshape(-1)
    edge_index = torch.stack([src, dst], dim=0)  # [2, N*k_eff]
    return edge_index


# ---------------------------------------------------------------------------
# Residue extraction from pocket PDB
# ---------------------------------------------------------------------------
class PocketParseError(Exception):
    pass


def extract_pocket_residues(pocket_pdb_path, strict=True):
    """
    Parse a pocket PDB, return list of (resid, resname, N_coord, CA_coord, C_coord, O_coord).

    Handles:
      - Standard AA filtering (via THREE_TO_ONE)
      - Single chain assertion (one chain only)
      - Insertion code rejection
      - Altloc: picks the first altloc per atom (BioPython default)
      - Full backbone (N/CA/C/O) required
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", str(pocket_pdb_path))

    # Check single model
    models = list(structure.get_models())
    if len(models) != 1:
        if strict:
            raise PocketParseError(f"expected 1 model, got {len(models)}")
    model = models[0]

    # Check chains
    chains = list(model.get_chains())
    if len(chains) > 1:
        if strict:
            raise PocketParseError(f"multi-chain pocket: {[c.id for c in chains]}")

    out = []
    chains_seen = set()
    for chain in chains:
        chains_seen.add(chain.id)
        for res in chain.get_residues():
            # res.id = (hetero_flag, resseq, icode)
            hetero, resseq, icode = res.id

            # Skip hetero (water, ligand, etc)
            if hetero.strip() != "":
                continue

            # Reject insertion codes
            if icode.strip() != "":
                if strict:
                    raise PocketParseError(
                        f"residue {resseq}{icode} has insertion code"
                    )
                continue

            # Standard AA only
            if res.resname not in THREE_TO_ONE:
                continue

            # Reject non-positive resid
            if resseq <= 0:
                if strict:
                    raise PocketParseError(f"non-positive resid {resseq}")
                continue

            # Full backbone required
            if not all(a in res for a in ("N", "CA", "C", "O")):
                continue

            # Extract backbone coords (BioPython picks default altloc)
            N = np.array(res["N"].coord, dtype=np.float32)
            CA = np.array(res["CA"].coord, dtype=np.float32)
            C = np.array(res["C"].coord, dtype=np.float32)
            O = np.array(res["O"].coord, dtype=np.float32)

            out.append((int(resseq), res.resname, N, CA, C, O))

    if not out:
        raise PocketParseError("no valid residues found")

    return out, list(chains_seen)


# ---------------------------------------------------------------------------
# build_gvp_features — the core single-pocket function
# ---------------------------------------------------------------------------
def build_gvp_features(pocket_pdb_path, k=K_KNN, strict=True):
    residues, chains = extract_pocket_residues(pocket_pdb_path, strict=strict)
    N = len(residues)

    # Build coords [N, 4, 3] in order (N, CA, C, O) per residue
    coords_list = []
    pocket_residue_idx_list = []
    for resid, _resname, N_c, CA_c, C_c, O_c in residues:
        coords_list.append(np.stack([N_c, CA_c, C_c, O_c], axis=0))  # [4, 3]
        pocket_residue_idx_list.append(resid - 1)

    coords = torch.from_numpy(np.stack(coords_list, axis=0)).float()  # [N, 4, 3]
    X_ca = coords[:, 1]  # [N, 3]

    # Validate all coords finite
    if not torch.isfinite(coords).all():
        if strict:
            raise PocketParseError("non-finite backbone coordinates")

    # Node features
    node_s = _dihedrals(coords)  # [N, 6]
    orientations = _orientations(X_ca)  # [N, 2, 3]
    sidechains = _sidechains(coords)  # [N, 3]
    node_v = torch.cat([orientations, sidechains.unsqueeze(-2)], dim=-2)  # [N, 3, 3]

    # Edges via custom kNN
    edge_index = knn_graph_custom(X_ca, k=k)  # [2, E]
    if edge_index.shape[1] == 0:
        raise PocketParseError("empty edge_index (pocket too small)")

    # Edge features
    E_vectors = X_ca[edge_index[0]] - X_ca[edge_index[1]]  # [E, 3]
    rbf = _rbf(E_vectors.norm(dim=-1), D_count=NUM_RBF)  # [E, NUM_RBF]
    pos_embeddings = _positional_embeddings(edge_index, num_embeddings=NUM_POSITIONAL)  # [E, NUM_POS]
    edge_s = torch.cat([rbf, pos_embeddings], dim=-1)  # [E, 32]
    edge_v = _normalize(E_vectors).unsqueeze(-2)  # [E, 1, 3]

    # nan-safe
    node_s, node_v, edge_s, edge_v = map(
        torch.nan_to_num, (node_s, node_v, edge_s, edge_v)
    )

    # Enforce dtypes
    out = {
        "node_s": node_s.float(),
        "node_v": node_v.float(),
        "edge_index": edge_index.long(),
        "edge_s": edge_s.float(),
        "edge_v": edge_v.float(),
        "pocket_residue_idx": torch.tensor(pocket_residue_idx_list, dtype=torch.long),
    }
    return out, chains


# ---------------------------------------------------------------------------
# Test on a single pocket + validate pocket_residue_idx vs sequence
# ---------------------------------------------------------------------------
def main():
    pocket_path = POCKET_DIR / f"{TEST_DOCK_INDEX}.pdb"
    if not pocket_path.exists():
        print(f"ERROR: {pocket_path} not found")
        sys.exit(1)

    print("=" * 60)
    print(f"Test on pocket: {pocket_path.name}")
    print("=" * 60)

    try:
        feats, chains = build_gvp_features(pocket_path)
    except PocketParseError as e:
        print(f"PARSE ERROR: {e}")
        sys.exit(1)

    N = feats["node_s"].shape[0]
    E = feats["edge_index"].shape[1]
    print(f"chains: {chains}")
    print(f"residues (N): {N}")
    print(f"edges (E):    {E}  (k=30 kNN → expected ~{min(N-1, 30)*N})")

    print()
    print("=== shapes and dtypes ===")
    for k, v in feats.items():
        print(f"  {k:22s} shape={tuple(v.shape)!s:18s} dtype={v.dtype}")

    print()
    print("=== value sanity ===")
    print(f"  node_s range:             [{feats['node_s'].min():.4f}, {feats['node_s'].max():.4f}]")
    print(f"  node_v norm range:        [{feats['node_v'].norm(dim=-1).min():.4f}, {feats['node_v'].norm(dim=-1).max():.4f}]")
    print(f"  edge_s range:             [{feats['edge_s'].min():.4f}, {feats['edge_s'].max():.4f}]")
    print(f"  edge_v norm range:        [{feats['edge_v'].norm(dim=-1).min():.4f}, {feats['edge_v'].norm(dim=-1).max():.4f}]")
    print(f"  pocket_residue_idx range: [{feats['pocket_residue_idx'].min().item()}, {feats['pocket_residue_idx'].max().item()}]")
    print(f"  edge_index range:         [{feats['edge_index'].min().item()}, {feats['edge_index'].max().item()}]")

    # Verify edges are local indices within [0, N-1]
    assert feats["edge_index"].min() >= 0 and feats["edge_index"].max() < N, \
        "edge_index out of range"
    print(f"  edge_index local: OK (all in [0, {N-1}])")

    print()
    print("=== pocket_residue_idx → enzyme sequence cross check ===")
    # dock=3 corresponds to enzyme_id=93 (from Step 1)
    with open(ENZ_CSV, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    enzyme_id = 93  # hard-coded for dock=3
    seq = rows[enzyme_id]["Protein sequence"].strip()
    uniprot = rows[enzyme_id]["uniprots"]
    print(f"  dock={TEST_DOCK_INDEX} → enzyme_id={enzyme_id}, uniprot={uniprot}, seq_len={len(seq)}")

    pidx = feats["pocket_residue_idx"].tolist()
    # Re-parse residues to get (resid, resname) pairs in same order
    residues, _ = extract_pocket_residues(pocket_path)
    n_match = 0
    for i in range(N):
        resid = residues[i][0]
        resname = residues[i][1]
        expected_aa = THREE_TO_ONE.get(resname, "?")
        actual_aa = seq[resid - 1] if 0 <= resid - 1 < len(seq) else "?"
        if expected_aa == actual_aa and pidx[i] == resid - 1:
            n_match += 1
    print(f"  residue-level aa match: {n_match}/{N}")
    assert n_match == N, "pocket_residue_idx alignment BROKEN"

    print()
    print("=" * 60)
    print(f"SMOKE TEST PASSED — single pocket feature extraction works")
    print(f"N={N}, E={E}, all features finite, all indices valid, aa match {n_match}/{N}")
    print("=" * 60)


if __name__ == "__main__":
    main()
'''


def main():
    import subprocess, sys
    print("[local driver] running single-pocket GVP test on server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=300,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
