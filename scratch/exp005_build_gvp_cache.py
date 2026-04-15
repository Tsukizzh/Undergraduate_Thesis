"""
EXP005 Phase 2: Build GVP feature cache for all 44090 training docks.

For each (enzyme, dock) pair in pt_cache_allfix_unified sidecars:
  1. Parse pocket PDB, extract residues with complete (N,CA,C,O) backbone
  2. Look up pocket_residue_idx from enzyme_resid_map.pt (built in Phase 1)
  3. Build GVP features:
     - node_s [N, 6]  (phi, psi, omega) as (cos, sin)
     - node_v [N, 3, 3]  (forward, backward, sidechain unit vectors)
     - edge_index [2, E]  (kNN on CA, k=30)
     - edge_s [E, 32]  (16 RBF + 16 positional)
     - edge_v [E, 1, 3]  (CA-CA unit vector)
     - pocket_residue_idx [N]  (0-indexed UniProt position, -1 if unmapped)
  4. Save to gvp_cache/samples/{dock//1000:03d}/gvp_{dock:06d}.pt

Multi-process with N_WORKERS=20. Each worker reads its own pocket and writes its
own .pt file. Failure log records any pocket that couldn't be built.

Mode:
  --mode smoke     5 known cases (dock 3/8444/24/985/30789)
  --mode full      all 44090 docks
  --mode verify    load all built files + sanity
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import argparse, csv, gc, json, math, os, sys, time, traceback
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import torch
import torch.nn.functional as F
import numpy as np

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"
OVERLAY = BASE / "data/pt_cache_dualgraph_allfix_unified"
RESID_MAP_PT = OVERLAY / "enzyme_resid_map.pt"
GVP_CACHE = OVERLAY / "gvp_cache"
SAMPLES_DIR = GVP_CACHE / "samples"
FAILURES_LOG = GVP_CACHE / "failures.log"
MANIFEST_PT = GVP_CACHE / "manifest.pt"
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
AA_TO_IDX = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")}

N_WORKERS = 20
KNN_K = 30
UNMAPPED_IDX = -1  # sentinel for residues without resid_map entry

BACKBONE_ATOMS = ("N", "CA", "C", "O")

TEST_DOCKS = [3, 8444, 24, 985, 30789]  # smoke test docks

# -------------------------------------------------------------------
# Pocket parser with backbone coordinates
# -------------------------------------------------------------------
def parse_pocket_backbone(pdb_path):
    """Return list of dicts per residue:
       {'chain', 'resid', 'icode', 'resname', 'aa1', 'coords': (N,CA,C,O) each (x,y,z)}
       Only residues with ALL 4 backbone atoms are kept. Sorted by (chain, resid, icode).
       First-occurrence wins per (chain, resid, icode, atom_name) (altloc A)."""
    atoms = {}
    resname_first = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line[:6].strip() == "COMPND":
                break
            if line[:6].strip() != "ATOM":
                continue
            rname = line[17:20].strip()
            if rname not in THREE_TO_ONE:
                continue
            atom_name = line[12:16].strip()
            if atom_name not in BACKBONE_ATOMS:
                continue
            chain = line[21]
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            key = (chain, resid, icode)
            if key not in resname_first:
                resname_first[key] = rname
            atom_map = atoms.setdefault(key, {})
            if atom_name not in atom_map:
                atom_map[atom_name] = (x, y, z)

    # Keep residues with full backbone
    kept = []
    for key in sorted(atoms.keys()):
        amap = atoms[key]
        if all(an in amap for an in BACKBONE_ATOMS):
            rname = resname_first[key]
            kept.append({
                "chain": key[0],
                "resid": key[1],
                "icode": key[2],
                "resname": rname,
                "aa1": THREE_TO_ONE[rname],
                "N":  amap["N"],
                "CA": amap["CA"],
                "C":  amap["C"],
                "O":  amap["O"],
            })
    return kept


# -------------------------------------------------------------------
# GVP helper math (ported from EnzymeCAGE gvp/data.py)
# -------------------------------------------------------------------
def _normalize(t, dim=-1, eps=1e-8):
    return torch.nan_to_num(t / (torch.norm(t, dim=dim, keepdim=True) + eps))


def _rbf(D, D_min=0.0, D_max=20.0, D_count=16):
    D_mu = torch.linspace(D_min, D_max, D_count).view(1, -1)
    D_sigma = (D_max - D_min) / D_count
    return torch.exp(-((D.unsqueeze(-1) - D_mu) / D_sigma) ** 2)


def _dihedrals(X, eps=1e-7):
    """X: [N, 3, 3] where the 3 atoms are (N, CA, C). Returns [N, 6]."""
    X = torch.reshape(X[:, :3], [3 * X.shape[0], 3])
    dX = X[1:] - X[:-1]
    U = _normalize(dX, dim=-1)
    u_2, u_1, u_0 = U[:-2], U[1:-1], U[2:]
    n_2 = _normalize(torch.cross(u_2, u_1, dim=-1), dim=-1)
    n_1 = _normalize(torch.cross(u_1, u_0, dim=-1), dim=-1)
    cosD = torch.clamp(torch.sum(n_2 * n_1, -1), -1 + eps, 1 - eps)
    D = torch.sign(torch.sum(u_2 * n_1, -1)) * torch.acos(cosD)
    D = F.pad(D, [1, 2])  # pad first phi, last psi, last omega
    D = torch.reshape(D, [-1, 3])  # [N, 3]
    return torch.cat([torch.cos(D), torch.sin(D)], 1)  # [N, 6]


def _orientations(X_ca):
    """X_ca: [N, 3]. Returns [N, 2, 3]."""
    forward = _normalize(X_ca[1:] - X_ca[:-1])
    backward = _normalize(X_ca[:-1] - X_ca[1:])
    forward = F.pad(forward, [0, 0, 0, 1])
    backward = F.pad(backward, [0, 0, 1, 0])
    return torch.cat([forward.unsqueeze(-2), backward.unsqueeze(-2)], -2)


def _sidechains(X):
    """X: [N, 4, 3] = (N, CA, C, O). Returns [N, 3] virtual sidechain direction."""
    n, origin, c = X[:, 0], X[:, 1], X[:, 2]
    c_n = _normalize(c - origin)
    n_n = _normalize(n - origin)
    bisector = _normalize(c_n + n_n)
    perp = _normalize(torch.cross(c_n, n_n, dim=-1))
    return -bisector * math.sqrt(1/3) - perp * math.sqrt(2/3)


def _positional_embeddings(edge_index, num_embeddings=16, period_range=(2, 1000)):
    """Positional encoding based on residue-index difference. Returns [E, 16]."""
    d = edge_index[0] - edge_index[1]
    frequency = torch.exp(
        torch.arange(0, num_embeddings, 2, dtype=torch.float32)
        * -(np.log(10000.0) / num_embeddings)
    )
    angles = d.unsqueeze(-1).float() * frequency
    return torch.cat((torch.cos(angles), torch.sin(angles)), -1)


def knn_graph_custom(X_ca, k=30):
    """kNN via torch.cdist + topk. Exclude self. Returns [2, E].

    Convention matches torch_cluster.knn_graph(flow='source_to_target'):
      edge_index[0] = source = neighbor
      edge_index[1] = target = self (center)
    So `ca[edge_index[0]] - ca[edge_index[1]]` = neighbor - self
    (vector FROM the receiving node TO the sending node), consistent
    with the standard GVP / EnzymeCAGE implementation.
    """
    N = X_ca.shape[0]
    if N < 2:
        return torch.zeros(2, 0, dtype=torch.long)
    dist = torch.cdist(X_ca, X_ca, p=2)
    dist.fill_diagonal_(float("inf"))
    k_eff = min(k, N - 1)
    _, nbr = torch.topk(dist, k_eff, dim=1, largest=False, sorted=True)
    # For each center i, add edges FROM each of i's k neighbors TO i.
    dst = torch.arange(N).unsqueeze(1).expand(-1, k_eff).reshape(-1)  # target = self
    src = nbr.reshape(-1)                                              # source = neighbor
    return torch.stack([src, dst], dim=0)


# -------------------------------------------------------------------
# GVP feature builder for a single pocket
# -------------------------------------------------------------------
def build_gvp_for_pocket(pdb_path, rmap):
    """
    Build GVP features for one pocket PDB using a pre-built resid_map.
    `rmap`: dict (chain, resid, icode) -> uniprot_pos_0idx.
    Returns dict or None on failure.
    """
    residues = parse_pocket_backbone(pdb_path)
    if not residues:
        return None, "empty_pocket"

    N = len(residues)
    if N < 2:
        return None, f"too_few_residues={N}"

    # Coordinates [N, 4, 3] (N, CA, C, O)
    coords = torch.tensor(
        [[r["N"], r["CA"], r["C"], r["O"]] for r in residues],
        dtype=torch.float32,
    )  # [N, 4, 3]

    # pocket_residue_idx
    prix = torch.full((N,), UNMAPPED_IDX, dtype=torch.long)
    n_mapped = 0
    for i, r in enumerate(residues):
        up = rmap.get((r["chain"], r["resid"], r["icode"]))
        if up is not None:
            prix[i] = int(up)
            n_mapped += 1

    # AA type [N] (for GVP seq_in embedding)
    aa_type = torch.tensor([AA_TO_IDX.get(r["aa1"], 0) for r in residues], dtype=torch.long)

    # Node scalar: dihedrals [N, 6]
    node_s = _dihedrals(coords)

    # Node vector: [N, 3, 3] = (forward_ca, backward_ca, sidechain)
    ca = coords[:, 1]  # [N, 3]
    ori = _orientations(ca)  # [N, 2, 3]
    side = _sidechains(coords).unsqueeze(-2)  # [N, 1, 3]
    node_v = torch.cat([ori, side], dim=-2)  # [N, 3, 3]

    # kNN edges
    edge_index = knn_graph_custom(ca, k=KNN_K)  # [2, E]
    if edge_index.shape[1] == 0:
        return None, "no_edges"

    # Edge scalar: RBF(distance) [E, 16] + positional [E, 16]
    src, dst = edge_index[0], edge_index[1]
    ca_diff = ca[src] - ca[dst]  # [E, 3]
    dist = torch.norm(ca_diff, dim=-1)  # [E]
    rbf_feat = _rbf(dist)  # [E, 16]
    pos_feat = _positional_embeddings(edge_index)  # [E, 16]
    edge_s = torch.cat([rbf_feat, pos_feat], dim=-1)  # [E, 32]

    # Edge vector: unit CA-CA diff [E, 1, 3]
    edge_v = _normalize(ca_diff).unsqueeze(-2)

    # Sanity: no NaN / Inf
    tensors = {"node_s": node_s, "node_v": node_v, "edge_s": edge_s, "edge_v": edge_v}
    for name, t in tensors.items():
        if torch.isnan(t).any() or torch.isinf(t).any():
            return None, f"nan_or_inf_in_{name}"

    return {
        "node_s": node_s.half(),               # [N, 6]
        "node_v": node_v.float(),              # [N, 3, 3]
        "edge_index": edge_index.int(),        # [2, E]
        "edge_s": edge_s.half(),               # [E, 32]
        "edge_v": edge_v.float(),              # [E, 1, 3]
        "pocket_residue_idx": prix,            # [N] long
        "aa_type": aa_type.to(torch.uint8),    # [N] uint8
        "n_mapped": n_mapped,
        "n_total": N,
    }, "ok"


# -------------------------------------------------------------------
# Worker
# -------------------------------------------------------------------
_rmap_by_enz = None

def _worker_init(rmap_path):
    # IMPORTANT: limit torch/OMP to a single thread per worker so that
    # 20 processes do not fight for cores (each torch op would otherwise
    # spawn as many threads as vCPUs, causing massive contention).
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    torch.set_num_threads(1)
    global _rmap_by_enz
    data = torch.load(str(rmap_path), weights_only=False)
    _rmap_by_enz = {}
    for enz_id, d in data.items():
        m = {}
        for k, v in d["resid_map"].items():
            c, r, ic = k.split("|")
            m[(c, int(r), ic)] = int(v)
        _rmap_by_enz[enz_id] = m


def _worker_process(job):
    enz_id, dock_index = job
    try:
        rmap = _rmap_by_enz.get(enz_id)
        if rmap is None:
            return {"dock": dock_index, "enz": enz_id, "status": "no_rmap"}
        ppath = POCKET_DIR / f"{dock_index}.pdb"
        if not ppath.exists():
            return {"dock": dock_index, "enz": enz_id, "status": "no_pocket_pdb"}
        feat, note = build_gvp_for_pocket(ppath, rmap)
        if feat is None:
            return {"dock": dock_index, "enz": enz_id, "status": note}
        # Save
        shard = dock_index // 1000
        out_dir = SAMPLES_DIR / f"{shard:03d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"gvp_{dock_index:06d}.pt"
        tmp_path = out_dir / f"tmp_{dock_index:06d}.pt"
        torch.save(feat, str(tmp_path))
        os.rename(str(tmp_path), str(out_path))
        return {
            "dock": dock_index,
            "enz": enz_id,
            "status": "ok",
            "N": feat["n_total"],
            "E": feat["edge_index"].shape[1],
            "n_mapped": feat["n_mapped"],
        }
    except Exception as e:
        return {
            "dock": dock_index,
            "enz": enz_id,
            "status": "exception",
            "error": f"{type(e).__name__}: {e}",
            "trace": traceback.format_exc()[:300],
        }


# -------------------------------------------------------------------
# Main
# -------------------------------------------------------------------
def load_jobs(mode):
    """Return list of (enz_id, dock_index) pairs, deduplicated by dock_index.

    Deduplication is a safety measure: for random_split each dock belongs to
    exactly one split, so there should be zero dupes, but if there ever are
    duplicates two workers would race to write the same file and the second
    writer would silently overwrite the first.
    """
    seen_docks = {}
    for sp in SIDECAR_PATHS:
        if not sp.exists():
            print(f"[WARN] sidecar missing: {sp}", flush=True)
            continue
        sc = torch.load(str(sp), weights_only=False)
        bi = torch.load(str(sp).replace("dock_sidecar.pt", "index.pt"), weights_only=False)
        eids = bi["enzyme_ids"].tolist()
        dids = sc["dock_indices"].tolist()
        for e, d in zip(eids, dids):
            d_int = int(d)
            if d_int in seen_docks:
                prev_e = seen_docks[d_int]
                if prev_e != int(e):
                    print(f"[WARN] dock {d_int} has two different enzymes: {prev_e} vs {int(e)}", flush=True)
                continue
            seen_docks[d_int] = int(e)

    jobs = [(enz, dock) for dock, enz in seen_docks.items()]

    if mode == "smoke":
        dock_set = set(TEST_DOCKS)
        jobs = [j for j in jobs if j[1] in dock_set]
    return jobs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["smoke", "full"], default="smoke")
    ap.add_argument("--workers", type=int, default=N_WORKERS)
    args = ap.parse_args()

    print("=" * 70)
    print(f"EXP005 GVP cache builder  mode={args.mode}  workers={args.workers}")
    print("=" * 70)
    t0 = time.time()

    SAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    jobs = load_jobs(args.mode)
    print(f"jobs: {len(jobs)}")
    if args.mode == "smoke":
        print(f"smoke docks: {sorted(set(j[1] for j in jobs))}")
    print(flush=True)

    # Run
    results = []
    stats = defaultdict(int)
    failures = []
    t_ali = time.time()
    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=_worker_init,
        initargs=(RESID_MAP_PT,),
    ) as ex:
        futs = [ex.submit(_worker_process, j) for j in jobs]
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            stats[r["status"]] += 1
            results.append(r)
            if r["status"] != "ok":
                failures.append(r)
            done += 1
            if done % 2000 == 0:
                print(f"  {done}/{len(jobs)}  ({time.time()-t_ali:.1f}s)  {dict(stats)}", flush=True)
    print(f"\ndone in {time.time()-t0:.1f}s", flush=True)
    print(f"stats: {dict(stats)}")

    # Manifest
    n_ok = stats.get("ok", 0)
    total_N = sum(r.get("N", 0) for r in results if r["status"] == "ok")
    total_E = sum(r.get("E", 0) for r in results if r["status"] == "ok")
    total_mapped = sum(r.get("n_mapped", 0) for r in results if r["status"] == "ok")
    total_unmapped = total_N - total_mapped
    manifest = {
        "n_jobs": len(jobs),
        "n_ok": n_ok,
        "n_failed": len(failures),
        "stats": dict(stats),
        "total_residues": total_N,
        "total_edges": total_E,
        "total_mapped": total_mapped,
        "total_unmapped": total_unmapped,
        "mode": args.mode,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    torch.save(manifest, str(MANIFEST_PT))
    print(f"manifest saved to {MANIFEST_PT}")

    # Failure log
    if failures:
        with open(FAILURES_LOG, "w") as f:
            for r in failures:
                f.write(json.dumps(r) + "\n")
        print(f"failures logged to {FAILURES_LOG}")
    else:
        # Remove stale log if any
        if FAILURES_LOG.exists():
            FAILURES_LOG.unlink()

    # Smoke mode: print detailed result per test dock
    if args.mode == "smoke":
        print()
        print("=" * 70)
        print("SMOKE RESULTS")
        print("=" * 70)
        for r in sorted(results, key=lambda x: x["dock"]):
            print(f"  dock={r['dock']:6d} enz={r['enz']:5d}: {r}")

    # Aggregate
    print()
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  total jobs:      {len(jobs)}")
    print(f"  ok:              {n_ok}")
    print(f"  failed:          {len(failures)}")
    print(f"  total residues:  {total_N}")
    print(f"  total edges:     {total_E}")
    print(f"  mapped/total:    {total_mapped}/{total_N} ({100*total_mapped/max(total_N,1):.2f}%)")
    if failures:
        by_status = defaultdict(int)
        for f_ in failures:
            by_status[f_["status"]] += 1
        print(f"  failure types:   {dict(by_status)}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
'''


def main():
    import argparse, subprocess, sys

    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], default="smoke")
    parser.add_argument("--workers", type=int, default=20)
    args = parser.parse_args()

    print(f"[local driver] running GVP cache build mode={args.mode} on server...")
    # Also export OMP_NUM_THREADS=1 at the parent level so it's inherited by
    # all child worker processes.
    remote_cmd = (
        f"export PATH=/root/miniconda3/bin:$PATH && "
        f"export OMP_NUM_THREADS=1 && export MKL_NUM_THREADS=1 && "
        f"python - --mode {args.mode} --workers {args.workers}"
    )
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj", remote_cmd],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=3600,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
