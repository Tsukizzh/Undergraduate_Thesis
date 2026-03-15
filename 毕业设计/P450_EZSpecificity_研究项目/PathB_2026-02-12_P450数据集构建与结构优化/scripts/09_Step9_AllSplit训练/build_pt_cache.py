# NOTE: This script only processes complex (protein+ligand) samples.
# Ligand-only samples (str_tag='ligand') are not supported because
# the current training config uses full_data=False.

"""
build_pt_cache.py — EZSpecificity .pt preprocessing pipeline (v1)

Converts LMDB / npy / csv data into a three-tier .pt cache that eliminates
runtime LMDB reads and pre-computes k-NN graphs, enabling 2-3× data-loading
throughput during training.

Three-tier storage
──────────────────
  ezspec_pt_v1/
  ├── manifest.pt                  # Schema + hyper-params
  ├── enzymes/
  │   └── esm_XXXX.pt              # fp16 ESM embeddings (unpadded, sharded)
  ├── substrates/
  │   └── grover_XXXX.pt           # fp16 GROVER embeddings + uint8 Morgan
  ├── train/
  │   ├── graph_XXXX.pt            # Pre-computed k-NN graphs (compact dtypes)
  │   └── index.pt                 # row-pointer index for the dataset loader
  ├── val/   (same structure)
  └── test/  (same structure)

Usage
─────
  python build_pt_cache.py --config train_allsplit_config.yml \\
                           --output-dir ./ezspec_pt_v1

  # Override shard size and worker count:
  python build_pt_cache.py --config train_allsplit_config.yml \\
                           --output-dir ./ezspec_pt_v1 \\
                           --shard-size 2048 --num-workers 0

Design decisions
────────────────
  * Enzyme/substrate/graph shards all support multi-worker (ProcessPoolExecutor)
    and resume (skip existing shard files). Each worker opens its own LMDB handles.
  * enzyme_idx and substrate_idx are LOCAL to each dataset (0-indexed).
    Global identity is (dataset_id, local_idx); enzyme_shard and
    substrate_shard are keyed by (dataset_id, local_idx) as well.
  * valid_idx logic exactly matches StructureSequence in data_representer.py:
      final = intersection(sequence_valid, structure_valid)
    where structure_valid = dock_idx present in complex LMDB keys
          sequence_valid  = enzyme_idx AND substrate_idx (AND grover) present
"""
from __future__ import annotations

import argparse
import math
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any

import lmdb
import numpy as np
import torch
import torch.nn.functional as F
import yaml
from torch_geometric.nn import knn_graph
from tqdm import tqdm

import pandas as pd

# Ensure src/ is on sys.path so pickle can deserialize LMDB objects
# (e.g. Datasets.Structure.utils.StructureComplexData)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_PROJECT_ROOT = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", "..", "..", "..", ".."))
_SRC_DIR = os.path.join(_PROJECT_ROOT, "src")
if _SRC_DIR not in sys.path:
    sys.path.insert(0, _SRC_DIR)

# ---------------------------------------------------------------------------
# Constants (mirrors transforms.py and cache_utils.py)
# ---------------------------------------------------------------------------
PROTEIN_ATOMIC_NUMBERS = [1, 6, 7, 8, 16, 34]      # H, C, N, O, S, Se
LIGAND_ATOMIC_NUMBERS  = [1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53]  # H,C,N,O,F,P,S,Cl,Fe,Br,I

# ATOM_FEATS from protein_ligand.py (order matters for ligand_atom_feature columns)
# {'AtomicNumber': 1, 'Aromatic': 1, 'Degree': 6, 'NumHs': 6, 'Hybridization': N}
# The hybridization count is rdkit version-dependent; we infer it at runtime.
ATOM_FEATS_STATIC = [
    ("AtomicNumber", 1),
    ("Aromatic",     1),
    ("Degree",       6),
    ("NumHs",        6),
]  # Hybridization appended after checking the LMDB

# ligand_atom_aux packs: aromatic(1), degree(1), num_hs(1), hybridization(1)
LIGAND_AUX_COLS = [1, 2, 3, 4]   # column indices in ligand_atom_feature (0 = AtomicNumber)

DATASET_NAMES = ["brenda", "Duf", "Esterase", "Gt_acceptor",
                 "Nitrilase", "Phosphatase", "Thiolase"]

STR_TAG_VOCAB = {"complex": 0, "ligand": 1}

# ---------------------------------------------------------------------------
# LMDB helpers
# ---------------------------------------------------------------------------

def _open_lmdb(path: str) -> lmdb.Environment:
    return lmdb.open(
        path,
        map_size=256 * 1024 ** 3,
        create=False,
        subdir=False,
        readonly=True,
        lock=False,
        readahead=False,
        meminit=False,
    )


def _lmdb_keys(path: str) -> set[bytes]:
    if not os.path.exists(path):
        return set()
    env = _open_lmdb(path)
    try:
        with env.begin(write=False) as txn:
            return set(txn.cursor().iternext(values=False))
    finally:
        env.close()


def _lmdb_get(env: lmdb.Environment, key: int | str) -> Any:
    with env.begin(write=False) as txn:
        raw = txn.get(str(key).encode())
    if raw is None:
        return None
    return pickle.loads(raw)


# ---------------------------------------------------------------------------
# Valid-idx computation (exact mirror of StructureSequence logic)
# ---------------------------------------------------------------------------

def compute_valid_idx(df: pd.DataFrame,
                      enzyme_lmdb_path: str,
                      reaction_lmdb_path: str,
                      grover_lmdb_path: str,
                      structure_lmdb_path: str) -> list[int]:
    """Return list of row-indices that are valid in both sequence and structure."""

    enzyme_keys   = _lmdb_keys(enzyme_lmdb_path)
    reaction_keys = _lmdb_keys(reaction_lmdb_path)
    grover_keys   = _lmdb_keys(grover_lmdb_path)
    complex_keys  = _lmdb_keys(structure_lmdb_path)

    seq_valid  = []
    str_valid  = set()

    for idx in range(len(df)):
        enz_idx = int(df.iloc[idx]["Enzyme Index"])
        sub_idx = int(df.iloc[idx]["Substrate Index"])
        dock_idx = int(df.iloc[idx]["Dock Index"])

        # Sequence validity: enzyme + reaction (+ grover if present)
        ok = (
            str(enz_idx).encode() in enzyme_keys
            and str(sub_idx).encode() in reaction_keys
        )
        if grover_keys:
            ok = ok and str(sub_idx).encode() in grover_keys
        if ok:
            seq_valid.append(idx)

        # Structure validity: dock_idx in complex LMDB
        if str(dock_idx).encode() in complex_keys:
            str_valid.add(idx)

    # Intersection preserving seq order (matches ordered_intersection in cache_utils.py)
    return [i for i in seq_valid if i in str_valid]


# ---------------------------------------------------------------------------
# k-NN graph computation (exact mirror of _build_knn_graph in transforms.py)
# ---------------------------------------------------------------------------

def build_knn_graph(ligand_pos: np.ndarray,
                    protein_pos: np.ndarray,
                    ligand_bond_index: np.ndarray,
                    k: int,
                    knn_device: str = "cpu") -> np.ndarray:
    """
    Return knn_edge_index (2, E) using torch_geometric knn_graph, with existing
    chemical bond edges removed — exactly mirroring transforms.py
    EdgeConnection._build_knn_graph.

    Args:
        knn_device: "cpu" or "cuda" / "cuda:0". GPU accelerates k-NN ~10x.
    """
    if protein_pos.shape[0] == 0:
        return np.zeros((2, 0), dtype=np.int64)

    all_pos = np.concatenate([ligand_pos, protein_pos], axis=0)  # (N, 3)
    pos_t = torch.from_numpy(all_pos.astype(np.float32)).to(knn_device)

    with torch.no_grad():
        knn_index = knn_graph(pos_t, k=k, flow='target_to_source')  # (2, N*k)

    # Move to CPU before Python loop (GPU .item() per-edge would be very slow)
    if knn_index.is_cuda:
        knn_index = knn_index.cpu()

    # Filter knn edges that already exist as bond edges (vectorized)
    if ligand_bond_index is not None and ligand_bond_index.shape[1] > 0:
        max_node = int(knn_index.max().item()) + 1
        knn_hash = knn_index[0] * max_node + knn_index[1]
        bond_t = torch.from_numpy(ligand_bond_index.astype(np.int64))
        bond_hash = bond_t[0] * max_node + bond_t[1]
        mask = ~torch.isin(knn_hash, bond_hash)
        knn_index = knn_index[:, mask]

    if knn_index.shape[1] == 0:
        return np.zeros((2, 0), dtype=np.int64)
    return knn_index.numpy().astype(np.int64)


# ---------------------------------------------------------------------------
# Enzyme shard builder
# ---------------------------------------------------------------------------

def _process_single_enzyme_shard(args_tuple):
    """Process a single enzyme shard in a worker process. Returns index metadata."""
    (shard_idx, batch_keys, enz_dir_str) = args_tuple

    S = len(batch_keys)
    enzyme_ids_arr = np.empty(S, dtype=np.int32)
    lengths_arr    = np.empty(S, dtype=np.int16)
    offsets_arr    = np.zeros(S + 1, dtype=np.int64)

    # Group keys by LMDB path to minimise open/close overhead
    path_to_indices: dict[str, list[tuple[int, int, int]]] = {}
    for row_in_shard, (ds_id, enz_idx, path) in enumerate(batch_keys):
        path_to_indices.setdefault(path, []).append((row_in_shard, ds_id, enz_idx))

    emb_by_row: list[np.ndarray | None] = [None] * S

    for path, entries in path_to_indices.items():
        env = _open_lmdb(path)
        try:
            with env.begin(write=False) as txn:
                for row_in_shard, ds_id, enz_idx in entries:
                    raw = txn.get(str(enz_idx).encode())
                    if raw is None:
                        raise RuntimeError(
                            f"Key {enz_idx} missing from enzyme LMDB {path}"
                        )
                    data = pickle.loads(raw)
                    emb = np.array(data["embedding"], dtype=np.float32)
                    enzyme_ids_arr[row_in_shard] = ds_id * 10_000_000 + enz_idx
                    lengths_arr[row_in_shard]    = emb.shape[0]
                    emb_by_row[row_in_shard]     = emb
        finally:
            env.close()

    # Build offsets and concatenate
    for i in range(S):
        offsets_arr[i + 1] = offsets_arr[i] + emb_by_row[i].shape[0]
    embedding = np.concatenate(emb_by_row, axis=0).astype(np.float16)

    shard = {
        "enzyme_ids": torch.from_numpy(enzyme_ids_arr),  # int32
        "lengths":    torch.from_numpy(lengths_arr),      # int16
        "offsets":    torch.from_numpy(offsets_arr),      # int64
        "embedding":  torch.from_numpy(embedding),        # float16
    }
    torch.save(shard, os.path.join(enz_dir_str, f"esm_{shard_idx:04d}.pt"))

    # Return index metadata: list of (global_id, shard_idx, row_in_shard)
    index_meta = []
    for row_in_shard in range(S):
        index_meta.append((int(enzyme_ids_arr[row_in_shard]), shard_idx, row_in_shard))

    del emb_by_row, embedding, shard
    return index_meta


def build_enzyme_shards(output_dir: Path,
                        config: dict,
                        shard_size: int,
                        num_workers: int = 0) -> None:
    """Iterate all enzyme LMDBs and write fp16 ESM shards.

    Streaming two-pass approach to avoid loading all embeddings into RAM:
      Pass 1: scan keys only (no values) to build (ds_id, enz_idx) list and assign shards.
      Pass 2: for each shard, load only its entries, write, then free.
    Supports resume (skips existing shards) and multi-worker via ProcessPoolExecutor.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    enz_dir = output_dir / "enzymes"
    enz_dir.mkdir(parents=True, exist_ok=True)
    enz_dir_str = str(enz_dir)

    enz_paths: list[str] = config["data"]["enzyme_lmdb_path"]

    print("\n[1/5] Building enzyme ESM shards ...")

    # Pass 1: collect (ds_id, enz_idx, path) tuples — keys only, no value loading
    all_keys: list[tuple[int, int, str]] = []   # (ds_id, enz_idx, lmdb_path)
    for ds_id, path in enumerate(tqdm(enz_paths, desc="Scanning enzyme keys")):
        if not os.path.exists(path):
            print(f"  WARNING: enzyme LMDB not found: {path}")
            continue
        env = _open_lmdb(path)
        try:
            with env.begin(write=False) as txn:
                for key_bytes in txn.cursor().iternext(keys=True, values=False):
                    all_keys.append((ds_id, int(key_bytes.decode()), path))
        finally:
            env.close()

    total = len(all_keys)
    total_shards = math.ceil(total / shard_size) if total > 0 else 0
    print(f"  {total} enzymes → {total_shards} shards")

    # Index arrays for enzymes/index.pt
    all_global_ids: list[int] = []
    all_shard_ids:  list[int] = []
    all_row_ids:    list[int] = []

    # Pass 2: build tasks, with resume check
    tasks = []
    skipped = 0
    for shard_idx in range(total_shards):
        shard_path = enz_dir / f"esm_{shard_idx:04d}.pt"
        batch_keys = all_keys[shard_idx * shard_size: (shard_idx + 1) * shard_size]

        if shard_path.exists():
            # Resume: recover index metadata from existing shard
            existing = torch.load(shard_path, map_location="cpu", weights_only=False)
            S_existing = len(existing["enzyme_ids"])
            for row_in_shard in range(S_existing):
                all_global_ids.append(int(existing["enzyme_ids"][row_in_shard].item()))
                all_shard_ids.append(shard_idx)
                all_row_ids.append(row_in_shard)
            del existing
            skipped += 1
            continue

        tasks.append((shard_idx, batch_keys, enz_dir_str))

    if skipped > 0:
        print(f"  Skipped {skipped} existing shards (resume).")

    if not tasks:
        print(f"  All enzyme shards already exist.")
    elif num_workers <= 1:
        # Single-process fallback
        for task in tqdm(tasks, desc="Writing enzyme shards"):
            index_meta = _process_single_enzyme_shard(task)
            for (gid, si, ri) in index_meta:
                all_global_ids.append(gid)
                all_shard_ids.append(si)
                all_row_ids.append(ri)
    else:
        # Multi-process
        print(f"  Using {num_workers} workers for {len(tasks)} enzyme shards ...")
        completed = 0
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(_process_single_enzyme_shard, t): t[0] for t in tasks}
            for future in as_completed(futures):
                shard_idx = futures[future]
                index_meta = future.result()
                for (gid, si, ri) in index_meta:
                    all_global_ids.append(gid)
                    all_shard_ids.append(si)
                    all_row_ids.append(ri)
                completed += 1
                print(f"    Enzyme shard {shard_idx:04d} done ({completed}/{len(tasks)})")

    # Write enzymes/index.pt for O(1) lookup at training time
    if all_global_ids:
        enz_index = {
            "enzyme_ids": torch.tensor(all_global_ids, dtype=torch.int32),   # int32[M]
            "shard_ids":  torch.tensor(all_shard_ids,  dtype=torch.int16),   # uint16[M]
            "row_ids":    torch.tensor(all_row_ids,    dtype=torch.int32),   # uint32[M]
        }
        torch.save(enz_index, enz_dir / "index.pt")

    print(f"  Enzyme shards written to {enz_dir}")


# ---------------------------------------------------------------------------
# Substrate shard builder
# ---------------------------------------------------------------------------

def _process_single_substrate_shard(args_tuple):
    """Process a single substrate shard in a worker process. Returns index metadata."""
    (shard_idx, batch_keys, morgan_paths, sub_dir_str) = args_tuple

    S = len(batch_keys)
    substrate_ids_arr = np.empty(S, dtype=np.int32)
    atom_lengths_arr  = np.empty(S, dtype=np.int16)
    atom_offsets_arr  = np.zeros(S + 1, dtype=np.int64)
    grover_atom_by_row: list[np.ndarray | None] = [None] * S
    grover_mean_by_row: list[np.ndarray | None] = [None] * S
    morgan_by_row:      list[np.ndarray | None] = [None] * S

    # Each worker loads its own Morgan mmaps
    morgan_mmap: list[np.ndarray | None] = []
    for mpath in morgan_paths:
        if os.path.exists(mpath):
            morgan_mmap.append(np.load(mpath, mmap_mode='r'))
        else:
            morgan_mmap.append(None)

    # Group by LMDB path to minimise open/close overhead
    path_to_indices: dict[str, list[tuple[int, int, int]]] = {}
    for row_in_shard, (ds_id, sub_idx, path) in enumerate(batch_keys):
        path_to_indices.setdefault(path, []).append((row_in_shard, ds_id, sub_idx))

    for path, entries in path_to_indices.items():
        env = _open_lmdb(path)
        try:
            with env.begin(write=False) as txn:
                for row_in_shard, ds_id, sub_idx in entries:
                    raw = txn.get(str(sub_idx).encode())
                    if raw is None:
                        raise RuntimeError(
                            f"Key {sub_idx} missing from grover LMDB {path}"
                        )
                    data = pickle.loads(raw)
                    g_atom = np.array(data["embedding"],       dtype=np.float32)
                    g_mean = np.array(data["total_embedding"], dtype=np.float32)

                    mmap = morgan_mmap[ds_id] if ds_id < len(morgan_mmap) else None
                    morgan = (mmap[sub_idx].copy()
                              if mmap is not None and sub_idx < len(mmap)
                              else np.zeros(1024, dtype=np.float32))

                    substrate_ids_arr[row_in_shard] = ds_id * 10_000_000 + sub_idx
                    atom_lengths_arr[row_in_shard]  = g_atom.shape[0]
                    grover_atom_by_row[row_in_shard] = g_atom
                    grover_mean_by_row[row_in_shard] = g_mean
                    morgan_by_row[row_in_shard]      = morgan
        finally:
            env.close()

    for i in range(S):
        atom_offsets_arr[i + 1] = atom_offsets_arr[i] + grover_atom_by_row[i].shape[0]

    grover_atom = np.concatenate(grover_atom_by_row, axis=0).astype(np.float16)
    grover_mean = np.stack(grover_mean_by_row, axis=0).astype(np.float16)
    morgan_data = np.stack(morgan_by_row, axis=0)

    # Store morgan as uint8 if all values are binary (0/1), else float16
    if np.all((morgan_data == 0) | (morgan_data == 1)):
        morgan_out = torch.from_numpy(morgan_data.astype(np.uint8))
    else:
        morgan_out = torch.from_numpy(morgan_data.astype(np.float16))

    shard = {
        "substrate_ids":  torch.from_numpy(substrate_ids_arr),
        "atom_lengths":   torch.from_numpy(atom_lengths_arr),
        "atom_offsets":   torch.from_numpy(atom_offsets_arr),
        "grover_atom":    torch.from_numpy(grover_atom),
        "grover_mean":    torch.from_numpy(grover_mean),
        "morgan":         morgan_out,
    }
    torch.save(shard, os.path.join(sub_dir_str, f"grover_{shard_idx:04d}.pt"))

    # Return index metadata: list of (global_id, shard_idx, row_in_shard)
    index_meta = []
    for row_in_shard in range(S):
        index_meta.append((int(substrate_ids_arr[row_in_shard]), shard_idx, row_in_shard))

    del grover_atom_by_row, grover_mean_by_row, morgan_by_row
    del grover_atom, grover_mean, morgan_data, shard
    return index_meta


def build_substrate_shards(output_dir: Path,
                            config: dict,
                            shard_size: int,
                            num_workers: int = 0) -> None:
    """Iterate all substrate LMDBs / npy and write fp16 GROVER + uint8 Morgan shards.

    Streaming two-pass approach to avoid loading all embeddings into RAM:
      Pass 1: scan GROVER LMDB keys only (no value loading) to assign shards.
      Pass 2: for each shard, load only its entries, write, then free.
    Supports resume (skips existing shards) and multi-worker via ProcessPoolExecutor.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    sub_dir = output_dir / "substrates"
    sub_dir.mkdir(parents=True, exist_ok=True)
    sub_dir_str = str(sub_dir)

    grover_paths: list[str] = config["data"]["grover_path"]
    morgan_paths: list[str] = config["data"]["morgan_path"]

    print("\n[2/5] Building substrate GROVER/Morgan shards ...")

    # Pass 1: collect (ds_id, sub_idx, grover_path) — keys only
    all_keys: list[tuple[int, int, str]] = []   # (ds_id, sub_idx, lmdb_path)
    for ds_id, (gpath, _) in enumerate(tqdm(
            zip(grover_paths, morgan_paths), desc="Scanning substrate keys",
            total=len(grover_paths))):
        if not os.path.exists(gpath):
            print(f"  WARNING: grover LMDB not found: {gpath}")
            continue
        env = _open_lmdb(gpath)
        try:
            with env.begin(write=False) as txn:
                for key_bytes in txn.cursor().iternext(keys=True, values=False):
                    all_keys.append((ds_id, int(key_bytes.decode()), gpath))
        finally:
            env.close()

    total = len(all_keys)
    total_shards = math.ceil(total / shard_size) if total > 0 else 0
    print(f"  {total} substrates → {total_shards} shards")

    # Index arrays for substrates/index.pt
    all_global_ids: list[int] = []
    all_shard_ids:  list[int] = []
    all_row_ids:    list[int] = []

    # Pass 2: build tasks, with resume check
    tasks = []
    skipped = 0
    for shard_idx in range(total_shards):
        shard_path = sub_dir / f"grover_{shard_idx:04d}.pt"
        batch_keys = all_keys[shard_idx * shard_size: (shard_idx + 1) * shard_size]

        if shard_path.exists():
            # Resume: recover index metadata from existing shard
            existing = torch.load(shard_path, map_location="cpu", weights_only=False)
            S_existing = len(existing["substrate_ids"])
            for row_in_shard in range(S_existing):
                all_global_ids.append(int(existing["substrate_ids"][row_in_shard].item()))
                all_shard_ids.append(shard_idx)
                all_row_ids.append(row_in_shard)
            del existing
            skipped += 1
            continue

        tasks.append((shard_idx, batch_keys, morgan_paths, sub_dir_str))

    if skipped > 0:
        print(f"  Skipped {skipped} existing shards (resume).")

    if not tasks:
        print(f"  All substrate shards already exist.")
    elif num_workers <= 1:
        # Single-process fallback
        for task in tqdm(tasks, desc="Writing substrate shards"):
            index_meta = _process_single_substrate_shard(task)
            for (gid, si, ri) in index_meta:
                all_global_ids.append(gid)
                all_shard_ids.append(si)
                all_row_ids.append(ri)
    else:
        # Multi-process
        print(f"  Using {num_workers} workers for {len(tasks)} substrate shards ...")
        completed = 0
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(_process_single_substrate_shard, t): t[0] for t in tasks}
            for future in as_completed(futures):
                shard_idx = futures[future]
                index_meta = future.result()
                for (gid, si, ri) in index_meta:
                    all_global_ids.append(gid)
                    all_shard_ids.append(si)
                    all_row_ids.append(ri)
                completed += 1
                print(f"    Substrate shard {shard_idx:04d} done ({completed}/{len(tasks)})")

    # Write substrates/index.pt for O(1) lookup at training time
    if all_global_ids:
        sub_index = {
            "substrate_ids": torch.tensor(all_global_ids, dtype=torch.int32),  # int32[M]
            "shard_ids":     torch.tensor(all_shard_ids,  dtype=torch.int16),  # uint16[M]
            "row_ids":       torch.tensor(all_row_ids,    dtype=torch.int32),  # uint32[M]
        }
        torch.save(sub_index, sub_dir / "index.pt")

    print(f"  Substrate shards written to {sub_dir}")


# ---------------------------------------------------------------------------
# Per-sample structure extractor
# ---------------------------------------------------------------------------

def extract_structure_sample(structure_data: Any,
                              k: int,
                              knn_device: str = "cpu") -> dict[str, np.ndarray]:
    """
    Extract compact arrays from one StructureComplexData (loaded from LMDB).

    Returns dict of numpy arrays matching the graph_XXXX.pt column spec.
    structure_data is a StructureComplexData (PyG Data-like) or a plain dict.
    """
    # Helper to get attribute from Data object or dict
    def _get(obj, attr):
        if isinstance(obj, dict):
            return obj.get(attr)
        return getattr(obj, attr, None)

    def _to_np(x):
        if x is None:
            return None
        if isinstance(x, torch.Tensor):
            return x.numpy()
        return np.asarray(x)

    ligand_pos      = _to_np(_get(structure_data, "ligand_pos"))      # (N_lig, 3)
    protein_pos     = _to_np(_get(structure_data, "protein_pos"))     # (N_prot, 3)
    ligand_element  = _to_np(_get(structure_data, "ligand_element"))  # (N_lig,)
    protein_element = _to_np(_get(structure_data, "protein_element")) # (N_prot,)
    protein_aa_type = _to_np(_get(structure_data, "protein_atom_to_aa_type"))  # (N_prot,)
    protein_is_bb   = _to_np(_get(structure_data, "protein_is_backbone"))      # (N_prot,)
    ligand_atom_feat = _to_np(_get(structure_data, "ligand_atom_feature"))     # (N_lig, >=5)
    ligand_bond_idx  = _to_np(_get(structure_data, "ligand_bond_index"))       # (2, E_bond)
    ligand_bond_type = _to_np(_get(structure_data, "ligand_bond_type"))        # (E_bond,)
    ligand_index     = _to_np(_get(structure_data, "ligand_index"))            # (N_lig,) atom-map nums

    if ligand_pos is None or protein_pos is None:
        raise ValueError("structure_data missing ligand_pos or protein_pos")

    N_lig  = ligand_pos.shape[0]
    N_prot = protein_pos.shape[0]

    # Ensure float32 positions
    ligand_pos  = ligand_pos.astype(np.float32)
    protein_pos = protein_pos.astype(np.float32)

    # Element as uint8
    ligand_element_u8  = ligand_element.astype(np.uint8)
    protein_element_u8 = protein_element.astype(np.uint8)

    # Protein atom-type and backbone flag
    if protein_aa_type is None:
        protein_aa_type = np.zeros(N_prot, dtype=np.uint8)
    if protein_is_bb is None:
        protein_is_bb = np.zeros(N_prot, dtype=np.uint8)
    protein_aa_type_u8 = protein_aa_type.astype(np.uint8)
    protein_is_bb_u8   = protein_is_bb.astype(np.uint8)

    # Ligand index (atom map numbers, 0..279)
    if ligand_index is None:
        ligand_index = np.zeros(N_lig, dtype=np.int32)
    # Guard against uint16 overflow before casting (values should be 0-279)
    if np.any(ligand_index > 65535):
        raise ValueError(
            f"ligand_index contains values > 65535: max={ligand_index.max()}"
        )
    ligand_index_u16 = ligand_index.astype(np.int32)

    # ligand_atom_aux: aromatic(col1), degree(col2), num_hs(col3), hybridization(col4)
    if ligand_atom_feat is not None and ligand_atom_feat.shape[1] >= 5:
        ligand_atom_aux = ligand_atom_feat[:, 1:5].astype(np.uint8)   # (N_lig, 4)
    elif ligand_atom_feat is not None:
        # Pad to 4 columns if fewer available
        ncols = ligand_atom_feat.shape[1] - 1
        aux = ligand_atom_feat[:, 1:].astype(np.uint8)
        ligand_atom_aux = np.pad(aux, ((0, 0), (0, max(0, 4 - ncols))))[:, :4]
    else:
        ligand_atom_aux = np.zeros((N_lig, 4), dtype=np.uint8)

    # Bonds
    if ligand_bond_idx is None or ligand_bond_idx.shape[1] == 0:
        bond_index = np.zeros((2, 0), dtype=np.int32)
        bond_type  = np.zeros(0, dtype=np.uint8)
    else:
        bond_index = ligand_bond_idx.astype(np.int32)        # (2, E_bond)
        bond_type  = (ligand_bond_type.astype(np.uint8)
                      if ligand_bond_type is not None
                      else np.zeros(bond_index.shape[1], dtype=np.uint8))

    # Pre-compute k-NN edges (mirrors transforms.py _build_knn_graph exactly).
    # Pass ligand_bond_idx (int64 view) so bond edges are filtered from k-NN result.
    knn_ei = build_knn_graph(ligand_pos, protein_pos,
                             ligand_bond_index=ligand_bond_idx.astype(np.int64)
                             if ligand_bond_idx is not None else np.zeros((2, 0), dtype=np.int64),
                             k=k, knn_device=knn_device)

    # Fix 5: assert indices fit in uint16 before casting
    if knn_ei.shape[1] > 0:
        assert knn_ei.max() < 65535, (
            f"Node index {knn_ei.max()} exceeds uint16 range"
        )
    if bond_index.shape[1] > 0:
        assert bond_index.max() < 65535, (
            f"Bond index {bond_index.max()} exceeds uint16 range"
        )

    knn_edge_index = knn_ei.astype(np.int32)   # node indices ≤ 65535

    return {
        "ligand_pos":         ligand_pos,           # float32 (N_lig, 3)
        "ligand_element":     ligand_element_u8,    # uint8   (N_lig,)
        "ligand_atom_aux":    ligand_atom_aux,       # uint8   (N_lig, 4)
        "ligand_index_raw":   ligand_index_u16,     # uint16  (N_lig,)
        "protein_pos":        protein_pos,           # float32 (N_prot, 3)
        "protein_element":    protein_element_u8,   # uint8   (N_prot,)
        "protein_aa_type":    protein_aa_type_u8,   # uint8   (N_prot,)
        "protein_is_backbone": protein_is_bb_u8,   # uint8   (N_prot,)
        "bond_index":         bond_index,            # uint16  (2, E_bond)
        "bond_type":          bond_type,             # uint8   (E_bond,)
        "knn_edge_index":     knn_edge_index,        # uint16  (2, E_knn)
    }


# ---------------------------------------------------------------------------
# Graph shard builder
# ---------------------------------------------------------------------------

def _process_single_shard(args_tuple):
    """Process a single shard in a worker process. Returns index metadata."""
    (shard_idx, batch_rows, df_rows_data, structure_paths, split_dir, k) = args_tuple

    import os as _os
    # Each worker opens its own LMDB handles
    worker_envs = []
    for path in structure_paths:
        if _os.path.exists(path):
            worker_envs.append(_open_lmdb(path))
        else:
            worker_envs.append(None)

    S = len(batch_rows)
    s_enzyme_ids    = np.zeros(S, dtype=np.int32)
    s_substrate_ids = np.zeros(S, dtype=np.int32)
    s_dataset_ids   = np.zeros(S, dtype=np.uint8)
    s_str_tag_codes = np.zeros(S, dtype=np.uint8)
    s_labels        = np.zeros(S, dtype=np.uint8)
    s_sample_weights = np.ones(S, dtype=np.float16)

    lig_pos_chunks, lig_elem_chunks, lig_aux_chunks, lig_idx_chunks = [], [], [], []
    prot_pos_chunks, prot_elem_chunks, prot_aa_chunks, prot_bb_chunks = [], [], [], []
    bond_idx_chunks, bond_type_chunks, knn_ei_chunks = [], [], []

    lig_ptr  = np.zeros(S + 1, dtype=np.int64)
    prot_ptr = np.zeros(S + 1, dtype=np.int64)
    bond_ptr = np.zeros(S + 1, dtype=np.int64)
    knn_ptr  = np.zeros(S + 1, dtype=np.int64)

    index_meta = []  # (sample_id, shard_idx, row, enz_id, sub_id)

    for row_in_shard, (idx, ds_id, enz_idx, sub_idx, dock_idx, label) in enumerate(
        zip(batch_rows, *zip(*df_rows_data))
    ):
        env = worker_envs[ds_id] if ds_id < len(worker_envs) else None
        if env is None:
            raise RuntimeError(f"Structure LMDB unavailable for dataset_id={ds_id}")
        raw = _lmdb_get(env, dock_idx)
        if raw is None:
            raise RuntimeError(f"dock_idx={dock_idx} not found (ds={ds_id})")
        struct_arrays = extract_structure_sample(raw, k=k, knn_device="cpu")

        s_enzyme_ids[row_in_shard]    = enz_idx
        s_substrate_ids[row_in_shard] = sub_idx
        s_dataset_ids[row_in_shard]   = ds_id
        s_labels[row_in_shard]        = label
        s_str_tag_codes[row_in_shard] = STR_TAG_VOCAB["complex"]

        N_lig  = struct_arrays["ligand_pos"].shape[0]
        N_prot = struct_arrays["protein_pos"].shape[0]
        E_bond = struct_arrays["bond_index"].shape[1]
        E_knn  = struct_arrays["knn_edge_index"].shape[1]

        lig_ptr[row_in_shard + 1]  = lig_ptr[row_in_shard]  + N_lig
        prot_ptr[row_in_shard + 1] = prot_ptr[row_in_shard] + N_prot
        bond_ptr[row_in_shard + 1] = bond_ptr[row_in_shard] + E_bond
        knn_ptr[row_in_shard + 1]  = knn_ptr[row_in_shard]  + E_knn

        lig_pos_chunks.append(struct_arrays["ligand_pos"])
        lig_elem_chunks.append(struct_arrays["ligand_element"])
        lig_aux_chunks.append(struct_arrays["ligand_atom_aux"])
        lig_idx_chunks.append(struct_arrays["ligand_index_raw"])
        prot_pos_chunks.append(struct_arrays["protein_pos"])
        prot_elem_chunks.append(struct_arrays["protein_element"])
        prot_aa_chunks.append(struct_arrays["protein_aa_type"])
        prot_bb_chunks.append(struct_arrays["protein_is_backbone"])
        bond_idx_chunks.append(struct_arrays["bond_index"])
        bond_type_chunks.append(struct_arrays["bond_type"])
        knn_ei_chunks.append(struct_arrays["knn_edge_index"])

        index_meta.append((idx, shard_idx, row_in_shard, enz_idx, sub_idx))

    # Concatenate
    def _cat(chunks, dtype):
        if not chunks: return np.zeros((0,), dtype=dtype)
        return np.concatenate([np.asarray(c, dtype=dtype) for c in chunks], axis=0)
    def _cat_2d_rows(chunks, dtype):
        if not chunks: return np.zeros((0, 4), dtype=dtype)
        return np.concatenate([np.asarray(c, dtype=dtype) for c in chunks], axis=0)
    def _cat_2d_cols(chunks, dtype):
        if not chunks or all(c.shape[1] == 0 for c in chunks):
            return np.zeros((2, 0), dtype=dtype)
        return np.concatenate([np.asarray(c, dtype=dtype) for c in chunks], axis=1)

    shard_data = {
        "num_samples": S,
        "enzyme_ids": torch.from_numpy(s_enzyme_ids),
        "substrate_ids": torch.from_numpy(s_substrate_ids),
        "dataset_ids": torch.from_numpy(s_dataset_ids),
        "str_tag_codes": torch.from_numpy(s_str_tag_codes),
        "labels": torch.from_numpy(s_labels),
        "sample_weights": torch.from_numpy(s_sample_weights),
        "ligand_ptr": torch.from_numpy(lig_ptr),
        "protein_ptr": torch.from_numpy(prot_ptr),
        "bond_ptr": torch.from_numpy(bond_ptr),
        "knn_ptr": torch.from_numpy(knn_ptr),
        "ligand_pos": torch.from_numpy(_cat(lig_pos_chunks, np.float32)).reshape(-1, 3),
        "ligand_element": torch.from_numpy(_cat(lig_elem_chunks, np.uint8)),
        "ligand_atom_aux": torch.from_numpy(_cat_2d_rows(lig_aux_chunks, np.uint8)),
        "ligand_index_raw": torch.from_numpy(_cat(lig_idx_chunks, np.int32)),
        "protein_pos": torch.from_numpy(_cat(prot_pos_chunks, np.float32)).reshape(-1, 3),
        "protein_element": torch.from_numpy(_cat(prot_elem_chunks, np.uint8)),
        "protein_aa_type": torch.from_numpy(_cat(prot_aa_chunks, np.uint8)),
        "protein_is_backbone": torch.from_numpy(_cat(prot_bb_chunks, np.uint8)),
        "bond_index": torch.from_numpy(_cat_2d_cols(bond_idx_chunks, np.int32)),
        "bond_type": torch.from_numpy(_cat(bond_type_chunks, np.uint8)),
        "knn_edge_index": torch.from_numpy(_cat_2d_cols(knn_ei_chunks, np.int32)),
    }
    shard_path = str(split_dir / f"graph_{shard_idx:04d}.pt")
    torch.save(shard_data, shard_path)

    # Close worker LMDB handles
    for env in worker_envs:
        if env is not None:
            env.close()

    return index_meta


def build_graph_shards(split_name: str,
                        df: pd.DataFrame,
                        valid_idx: list[int],
                        config: dict,
                        output_dir: Path,
                        shard_size: int,
                        k: int,
                        knn_device: str = "cpu",
                        num_workers: int = 0) -> None:
    """Build graph shards and index.pt for one data split."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    split_dir = output_dir / split_name
    split_dir.mkdir(parents=True, exist_ok=True)

    structure_paths: list[str] = config["data"]["structure_processed_path"]

    sample_ids_list:     list[int] = []
    graph_shards_list:   list[int] = []
    graph_rows_list:     list[int] = []
    enzyme_ids_list:     list[int] = []
    substrate_ids_list:  list[int] = []

    n_valid = len(valid_idx)
    n_shards = math.ceil(n_valid / shard_size) if n_valid > 0 else 0

    # Pre-extract df row data as plain lists (picklable for multiprocessing)
    all_ds_ids   = df["dataset_id"].values
    all_enz_ids  = df["Enzyme Index"].values
    all_sub_ids  = df["Substrate Index"].values
    all_dock_ids = df["Dock Index"].values
    all_labels   = df["Label"].values

    # Build shard tasks
    tasks = []
    skipped = 0
    for shard_idx in range(n_shards):
        shard_path = split_dir / f"graph_{shard_idx:04d}.pt"
        if shard_path.exists():
            # Resume: recover index metadata from existing shard
            existing = torch.load(shard_path, map_location="cpu", weights_only=False)
            S_existing = int(existing["num_samples"])
            batch_rows = valid_idx[shard_idx * shard_size: (shard_idx + 1) * shard_size]
            for r, idx in enumerate(batch_rows[:S_existing]):
                sample_ids_list.append(idx)
                graph_shards_list.append(shard_idx)
                graph_rows_list.append(r)
                enzyme_ids_list.append(int(existing["enzyme_ids"][r]))
                substrate_ids_list.append(int(existing["substrate_ids"][r]))
            del existing
            skipped += 1
            continue

        batch_rows = valid_idx[shard_idx * shard_size: (shard_idx + 1) * shard_size]
        df_rows_data = [
            (int(all_ds_ids[i]), int(all_enz_ids[i]), int(all_sub_ids[i]),
             int(all_dock_ids[i]), int(all_labels[i]))
            for i in batch_rows
        ]
        tasks.append((shard_idx, batch_rows, df_rows_data, structure_paths, split_dir, k))

    if skipped > 0:
        print(f"    Skipped {skipped} existing shards (resume).")

    if not tasks:
        print(f"    All shards already exist.")
    elif num_workers <= 1:
        # Single-process fallback
        for task in tqdm(tasks, desc=f"  [{split_name}] graph shards"):
            index_meta = _process_single_shard(task)
            for (idx, si, ri, ei, subi) in index_meta:
                sample_ids_list.append(idx)
                graph_shards_list.append(si)
                graph_rows_list.append(ri)
                enzyme_ids_list.append(ei)
                substrate_ids_list.append(subi)
    else:
        # Multi-process
        print(f"    Using {num_workers} workers for {len(tasks)} shards ...")
        completed = 0
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(_process_single_shard, t): t[0] for t in tasks}
            for future in as_completed(futures):
                shard_idx = futures[future]
                index_meta = future.result()
                for (idx, si, ri, ei, subi) in index_meta:
                    sample_ids_list.append(idx)
                    graph_shards_list.append(si)
                    graph_rows_list.append(ri)
                    enzyme_ids_list.append(ei)
                    substrate_ids_list.append(subi)
                completed += 1
                print(f"    [{split_name}] {completed}/{len(tasks)} shards done (shard {shard_idx:04d})")

    # Write index.pt
    if sample_ids_list:
        index_data = {
            "sample_ids":    torch.tensor(sample_ids_list,    dtype=torch.int32),
            "graph_shards":  torch.tensor(graph_shards_list,  dtype=torch.int16),
            "graph_rows":    torch.tensor(graph_rows_list,    dtype=torch.int32),
            "enzyme_ids":    torch.tensor(enzyme_ids_list,    dtype=torch.int32),
            "substrate_ids": torch.tensor(substrate_ids_list, dtype=torch.int32),
        }
    else:
        index_data = {
            "sample_ids":    torch.zeros(0, dtype=torch.int32),
            "graph_shards":  torch.zeros(0, dtype=torch.int16),
            "graph_rows":    torch.zeros(0, dtype=torch.int32),
            "enzyme_ids":    torch.zeros(0, dtype=torch.int32),
            "substrate_ids": torch.zeros(0, dtype=torch.int32),
        }
    torch.save(index_data, split_dir / "index.pt")
    print(f"  [{split_name}] {len(sample_ids_list)} samples → {n_shards} shards")


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------

def write_manifest(output_dir: Path, config: dict, k: int) -> None:
    cutoff      = config.get("transform", {}).get("cutoff", 10.0)
    n_gaussian  = config.get("transform", {}).get("num_r_gaussian", 32)
    max_enz_len = config.get("data", {}).get("max_enzyme_length", 1450)
    max_sub_len = config.get("data", {}).get("max_substrate_length", 280)

    manifest = {
        "version":            1,
        "k":                  k,
        "cutoff":             float(cutoff),
        "num_r_gaussian":     int(n_gaussian),
        "max_enzyme_len":     int(max_enz_len),
        "max_substrate_len":  int(max_sub_len),
        "protein_elements":   PROTEIN_ATOMIC_NUMBERS,
        "ligand_elements":    LIGAND_ATOMIC_NUMBERS,
        "dataset_tags":       [str(i) for i in range(len(DATASET_NAMES))],
        "str_tag_vocab":      ["complex", "ligand"],
    }
    torch.save(manifest, output_dir / "manifest.pt")
    print(f"  manifest.pt written to {output_dir}")


# ---------------------------------------------------------------------------
# Self-test: load one sample from train and verify shapes
# ---------------------------------------------------------------------------

def self_test(output_dir: Path) -> None:
    """
    Full reconstruction self-test for one sample from the train split.

    Steps:
      1. Load one sample from a graph shard.
      2. Load corresponding ESM embedding from an enzyme shard.
      3. Load corresponding GROVER/Morgan from a substrate shard.
      4. Reconstruct protein_x (one-hot element + aa_type + backbone flag).
      5. Reconstruct ligand_x (one-hot element + aux features).
      6. Reconstruct complex_edge_index (bond + knn, matching runtime order in transforms.py).
      7. Recompute edge distances from coords.
      8. Build ligand_mask and protein_mask.
      9. Print all tensor shapes and verify expected dimensions.
    """
    print("\n[5/5] Self-test ...")
    train_dir = output_dir / "train"
    index_path = train_dir / "index.pt"
    if not index_path.exists():
        print("  SKIP: train/index.pt not found")
        return

    index = torch.load(index_path, weights_only=False)
    n = len(index["sample_ids"])
    if n == 0:
        print("  SKIP: 0 train samples")
        return

    # ── 1. Load graph shard ──────────────────────────────────────────────────
    shard_id = int(index["graph_shards"][0].item())
    row_id   = int(index["graph_rows"][0].item())
    shard = torch.load(train_dir / f"graph_{shard_id:04d}.pt", weights_only=False)

    s_lig  = int(shard["ligand_ptr"][row_id].item())
    e_lig  = int(shard["ligand_ptr"][row_id + 1].item())
    s_prot = int(shard["protein_ptr"][row_id].item())
    e_prot = int(shard["protein_ptr"][row_id + 1].item())
    s_bond = int(shard["bond_ptr"][row_id].item())
    e_bond = int(shard["bond_ptr"][row_id + 1].item())
    s_knn  = int(shard["knn_ptr"][row_id].item())
    e_knn  = int(shard["knn_ptr"][row_id + 1].item())

    N_lig  = e_lig  - s_lig
    N_prot = e_prot - s_prot
    E_bond = e_bond - s_bond
    E_knn  = e_knn  - s_knn
    N_total = N_lig + N_prot

    ligand_pos      = shard["ligand_pos"][s_lig:e_lig].float()        # (N_lig, 3)
    protein_pos     = shard["protein_pos"][s_prot:e_prot].float()     # (N_prot, 3)
    ligand_elem     = shard["ligand_element"][s_lig:e_lig].long()     # (N_lig,)
    protein_elem    = shard["protein_element"][s_prot:e_prot].long()  # (N_prot,)
    protein_aa      = shard["protein_aa_type"][s_prot:e_prot].long()  # (N_prot,)
    protein_bb      = shard["protein_is_backbone"][s_prot:e_prot].long()  # (N_prot,)
    ligand_aux      = shard["ligand_atom_aux"][s_lig:e_lig].long()    # (N_lig, 4)
    bond_idx_raw    = shard["bond_index"][:, s_bond:e_bond].long()    # (2, E_bond)
    knn_idx_raw     = shard["knn_edge_index"][:, s_knn:e_knn].long()  # (2, E_knn)

    print(f"  Sample 0: N_lig={N_lig}, N_prot={N_prot}, E_bond={E_bond}, E_knn={E_knn}")
    print(f"  label={int(shard['labels'][row_id].item())}, "
          f"ds={int(shard['dataset_ids'][row_id].item())}, "
          f"str_tag={int(shard['str_tag_codes'][row_id].item())}")

    # ── 2. Load ESM embedding matching this sample's enzyme_id ───────────────
    enz_dir = output_dir / "enzymes"
    sample_enz_id = int(shard["enzyme_ids"][row_id].item())
    esm_emb: torch.Tensor | None = None
    enz_index_path = enz_dir / "index.pt"
    if enz_index_path.exists():
        enz_index = torch.load(enz_index_path, weights_only=False)
        mask = enz_index["enzyme_ids"] == sample_enz_id
        if mask.any():
            loc = mask.nonzero(as_tuple=True)[0][0].item()
            target_shard = int(enz_index["shard_ids"][loc].item())
            target_row = int(enz_index["row_ids"][loc].item())
            enz_shard = torch.load(enz_dir / f"esm_{target_shard:03d}.pt", weights_only=False)
            e_off_s = int(enz_shard["offsets"][target_row].item())
            e_off_e = int(enz_shard["offsets"][target_row + 1].item())
            esm_emb = enz_shard["embedding"][e_off_s:e_off_e].float()
            print(f"  ESM embedding shape: {esm_emb.shape}  (enzyme_id={sample_enz_id}, fp16→fp32)")
        else:
            print(f"  WARNING: enzyme_id={sample_enz_id} not found in index")
    else:
        print("  WARNING: enzymes/index.pt not found, skipping ESM check")

    # ── 3. Load GROVER/Morgan matching this sample's substrate_id ────────────
    sub_dir = output_dir / "substrates"
    sample_sub_id = int(shard["substrate_ids"][row_id].item())
    grover_mean: torch.Tensor | None = None
    grover_atom: torch.Tensor | None = None
    morgan: torch.Tensor | None = None
    sub_index_path = sub_dir / "index.pt"
    if sub_index_path.exists():
        sub_index = torch.load(sub_index_path, weights_only=False)
        mask = sub_index["substrate_ids"] == sample_sub_id
        if mask.any():
            loc = mask.nonzero(as_tuple=True)[0][0].item()
            target_shard = int(sub_index["shard_ids"][loc].item())
            target_row = int(sub_index["row_ids"][loc].item())
            sub_shard = torch.load(sub_dir / f"grover_{target_shard:03d}.pt", weights_only=False)
            a_off_s = int(sub_shard["atom_offsets"][target_row].item())
            a_off_e = int(sub_shard["atom_offsets"][target_row + 1].item())
            grover_atom = sub_shard["grover_atom"][a_off_s:a_off_e].float()
            grover_mean = sub_shard["grover_mean"][target_row].float()
            print(f"  GROVER atom shape:   {grover_atom.shape}  (substrate_id={sample_sub_id})")
        else:
            print(f"  WARNING: substrate_id={sample_sub_id} not found in index")
    elif sorted(sub_dir.glob("grover_*.pt")):
        sub_shard = torch.load(sorted(sub_dir.glob("grover_*.pt"))[0], weights_only=False)
        a_off_s = int(sub_shard["atom_offsets"][0].item())
        a_off_e = int(sub_shard["atom_offsets"][1].item())
        grover_atom = sub_shard["grover_atom"][a_off_s:a_off_e].float()
        grover_mean = sub_shard["grover_mean"][0].float()
        morgan      = sub_shard["morgan"][sub_idx_in_shard].float()       # (1024,)
        print(f"  GROVER atom shape:   {grover_atom.shape}")
        print(f"  GROVER mean shape:   {grover_mean.shape}")
        print(f"  Morgan shape:        {morgan.shape}")

    # ── 4. Reconstruct protein_x: one-hot(element, 6) + one-hot(aa_type, 21) + bb(1)
    # Mirrors FeaturizeProteinAtom in transforms.py:
    #   element  = protein_element.view(-1,1) == atomic_numbers.view(1,-1)  → (N_prot, 6)
    #   amino_acid = F.one_hot(protein_atom_to_aa_type, num_classes=21)     → (N_prot, 21)
    #   is_backbone = protein_is_backbone.view(-1,1)                        → (N_prot, 1)
    # Then in EdgeConnection.__call__:
    #   protein_x = cat([zeros(N_lig, F_prot), protein_atom_feature], dim=0)
    PROT_ELEM_T = torch.tensor(PROTEIN_ATOMIC_NUMBERS, dtype=torch.long)  # [6]
    prot_elem_oh = (protein_elem.unsqueeze(1) == PROT_ELEM_T.unsqueeze(0)).float()  # (N_prot, 6)
    prot_aa_oh   = F.one_hot(protein_aa.clamp(0, 20), num_classes=21).float()       # (N_prot, 21)
    prot_bb_f    = protein_bb.float().unsqueeze(1)                                  # (N_prot, 1)
    protein_atom_feature = torch.cat([prot_elem_oh, prot_aa_oh, prot_bb_f], dim=1) # (N_prot, 28)
    protein_x = torch.cat([
        torch.zeros(N_lig, protein_atom_feature.shape[1]),
        protein_atom_feature,
    ], dim=0)  # (N_total, 28)
    print(f"  protein_x shape:     {protein_x.shape}  (expected ({N_total}, 28))")

    # ── 5. Reconstruct ligand_x: one-hot(element, 11) + ATOM_FEATS encodings ─
    # Mirrors FeaturizeLigandAtom in transforms.py, iterating ATOM_FEATS:
    #   AtomicNumber (v=1): feat / 100.0               → (N_lig, 1)
    #   Aromatic     (v=1): raw 0/1                    → (N_lig, 1)
    #   Degree       (v=6): one-hot over range(6)      → (N_lig, 6)
    #   NumHs        (v=6): one-hot over range(6)      → (N_lig, 6)
    #   Hybridization(v=9): one-hot over range(9)      → (N_lig, 9)
    # ligand_atom_aux stores cols [aromatic, degree, num_hs, hybridization] (uint8)
    # ligand_element stores the raw atomic number (uint8)
    # Total ligand feature dim = 11 + 1 + 1 + 6 + 6 + 9 = 34
    LIG_ELEM_T = torch.tensor(LIGAND_ATOMIC_NUMBERS, dtype=torch.long)  # [11]
    lig_elem_oh = (ligand_elem.unsqueeze(1) == LIG_ELEM_T.unsqueeze(0)).float()  # (N_lig, 11)

    atomic_num_feat  = ligand_elem.float().unsqueeze(1) / 100.0            # (N_lig, 1)
    aromatic_feat    = ligand_aux[:, 0].float().unsqueeze(1)               # (N_lig, 1)
    degree_oh        = F.one_hot(ligand_aux[:, 1].clamp(0, 5),  num_classes=6).float()  # (N_lig, 6)
    num_hs_oh        = F.one_hot(ligand_aux[:, 2].clamp(0, 5),  num_classes=6).float()  # (N_lig, 6)
    hybrid_oh        = F.one_hot(ligand_aux[:, 3].clamp(0, 8),  num_classes=9).float()  # (N_lig, 9)

    ligand_atom_feature_full = torch.cat(
        [lig_elem_oh, atomic_num_feat, aromatic_feat, degree_oh, num_hs_oh, hybrid_oh], dim=1
    )  # (N_lig, 34)
    ligand_x = torch.cat([
        ligand_atom_feature_full,
        torch.zeros(N_prot, ligand_atom_feature_full.shape[1]),
    ], dim=0)  # (N_total, 34)
    LIG_FEAT_DIM = ligand_atom_feature_full.shape[1]
    print(f"  ligand_x shape:      {ligand_x.shape}  (expected ({N_total}, {LIG_FEAT_DIM}))")

    # ── 6. Reconstruct complex_edge_index (runtime order: bond then knn) ────
    # transforms.py line 147: complex_edge_index = cat([ligand_bond_index, knn_index], dim=-1)
    # Note: knn_index here is ALREADY filtered (no bond edges), matching our cache.
    complex_edge_index = torch.cat([bond_idx_raw, knn_idx_raw], dim=1)  # (2, E_bond+E_knn)
    E_total = complex_edge_index.shape[1]
    print(f"  complex_edge_index:  {complex_edge_index.shape}  "
          f"(bond={E_bond}, knn={E_knn}, total={E_total})")

    # ── 7. Recompute edge distances from coords ──────────────────────────────
    pos_all = torch.cat([ligand_pos, protein_pos], dim=0)  # (N_total, 3)
    if E_total > 0:
        dst, src = complex_edge_index
        dists = torch.norm(pos_all[dst] - pos_all[src], p=2, dim=-1)
        print(f"  edge dist range:     [{dists.min():.2f}, {dists.max():.2f}] Å")

    # ── 8. Build padding masks (matches EdgeConnection output) ───────────────
    ligand_mask  = torch.cat([torch.ones(N_lig),  torch.zeros(N_prot)])  # (N_total,)
    protein_mask = torch.cat([torch.zeros(N_lig), torch.ones(N_prot)])   # (N_total,)
    print(f"  ligand_mask  sum:    {ligand_mask.sum().int()}  (expected {N_lig})")
    print(f"  protein_mask sum:    {protein_mask.sum().int()}  (expected {N_prot})")

    # ── 9. Verify shapes match ss.py forward() expectations ─────────────────
    # ss.py: G.embedding → (B*1450, 1280),  G.grover → (B*280, 2400)
    # enzyme_padding_mask: (B, 1450),  reaction_padding_mask: (B, 280)
    if esm_emb is not None:
        L_enz = esm_emb.shape[0]
        assert esm_emb.shape[1] == 1280, f"ESM dim mismatch: {esm_emb.shape}"
        assert L_enz <= 1450, f"ESM length {L_enz} > max_enzyme_len 1450"
        # enzyme_padding_mask construction: True for padding positions
        enzyme_padding_mask = torch.zeros(1, 1450, dtype=torch.bool)
        enzyme_padding_mask[0, L_enz:] = True
        print(f"  enzyme_padding_mask: {enzyme_padding_mask.shape}, "
              f"padding={enzyme_padding_mask.sum().item()} positions")

    if grover_atom is not None:
        N_atoms = grover_atom.shape[0]
        assert grover_atom.shape[1] == 2400, f"GROVER atom dim mismatch: {grover_atom.shape}"
        assert N_atoms <= 280, f"GROVER atoms {N_atoms} > max_substrate_len 280"
        # reaction_padding_mask: True for padding positions
        reaction_padding_mask = torch.zeros(1, 280, dtype=torch.bool)
        reaction_padding_mask[0, N_atoms:] = True
        print(f"  reaction_padding_mask: {reaction_padding_mask.shape}, "
              f"padding={reaction_padding_mask.sum().item()} positions")

    # Core shape assertions
    assert protein_x.shape == (N_total, 28),        f"protein_x shape mismatch: {protein_x.shape}"
    assert ligand_x.shape  == (N_total, LIG_FEAT_DIM), f"ligand_x shape mismatch: {ligand_x.shape}"
    assert ligand_mask.shape  == (N_total,),         f"ligand_mask shape mismatch"
    assert protein_mask.shape == (N_total,),         f"protein_mask shape mismatch"
    assert complex_edge_index.shape[0] == 2,         f"edge_index row count wrong"
    if E_knn > 0:
        assert knn_idx_raw.max() < 65535, "knn index exceeds uint16 range"
    # Bond edges must NOT appear in knn edges (no duplicates)
    if E_bond > 0 and E_knn > 0:
        bond_set_test = set(zip(bond_idx_raw[0].tolist(), bond_idx_raw[1].tolist()))
        knn_set_test  = set(zip(knn_idx_raw[0].tolist(),  knn_idx_raw[1].tolist()))
        overlap = bond_set_test & knn_set_test
        assert len(overlap) == 0, (
            f"Bond/knn edge overlap: {len(overlap)} duplicate edges found"
        )

    print("  Self-test PASSED")


# ---------------------------------------------------------------------------
# Summary stats
# ---------------------------------------------------------------------------

def print_summary(output_dir: Path, t_start: float) -> None:
    print("\n========== Build Summary ==========")
    total_bytes = 0
    for root, _, files in os.walk(output_dir):
        for f in files:
            total_bytes += os.path.getsize(os.path.join(root, f))
    total_gb = total_bytes / (1024 ** 3)

    for split in ("train", "val", "test"):
        idx_path = output_dir / split / "index.pt"
        if idx_path.exists():
            idx = torch.load(idx_path, weights_only=False)
            n_graphs = len(sorted((output_dir / split).glob("graph_*.pt")))
            print(f"  {split:5s}: {len(idx['sample_ids']):7d} samples, {n_graphs} graph shards")

    enz_shards = len(sorted((output_dir / "enzymes").glob("esm_*.pt")))
    sub_shards = len(sorted((output_dir / "substrates").glob("grover_*.pt")))
    print(f"  enzymes: {enz_shards} shards | substrates: {sub_shards} shards")
    print(f"  Total size: {total_gb:.2f} GB")
    print(f"  Elapsed: {(time.time() - t_start):.0f}s")
    print("====================================")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build EZSpecificity .pt cache from LMDB/npy/csv sources."
    )
    p.add_argument(
        "--config",
        default=os.path.join(os.path.dirname(__file__), "train_allsplit_config.yml"),
        help="Path to train_allsplit_config.yml (or compatible YAML).",
    )
    p.add_argument(
        "--output-dir",
        default=os.path.join(os.path.dirname(__file__), "ezspec_pt_v1"),
        help="Root directory for output .pt files.",
    )
    p.add_argument(
        "--data-dir",
        default=None,
        help="If set, override the base data directory (replaces the common prefix "
             "in all paths inside the YAML). Useful for running on a different machine.",
    )
    p.add_argument(
        "--shard-size",
        type=int,
        default=2048,
        help="Number of samples per graph shard (default: 2048).",
    )
    p.add_argument(
        "--enzyme-shard-size",
        type=int,
        default=4096,
        help="Number of enzymes per ESM shard (default: 4096).",
    )
    p.add_argument(
        "--substrate-shard-size",
        type=int,
        default=4096,
        help="Number of substrates per GROVER shard (default: 4096).",
    )
    p.add_argument(
        "--num-workers",
        type=int,
        default=0,
        help="Number of parallel workers for shard generation (0 = sequential).",
    )
    p.add_argument(
        "--splits",
        nargs="+",
        default=["train", "val", "test"],
        choices=["train", "val", "test"],
        help="Which splits to build (default: all three).",
    )
    p.add_argument(
        "--skip-enzyme-shards",
        action="store_true",
        help="Skip (re-)building enzyme shards.",
    )
    p.add_argument(
        "--skip-substrate-shards",
        action="store_true",
        help="Skip (re-)building substrate shards.",
    )
    p.add_argument(
        "--knn-device",
        default="cpu",
        help="Device for knn_graph computation: cpu, cuda, cuda:0. GPU is ~10x faster.",
    )
    return p.parse_args()


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def maybe_remap_paths(config: dict, data_dir: str | None) -> None:
    """If --data-dir is given, rewrite all path strings under config['data']."""
    if data_dir is None:
        return
    data_dir = os.path.abspath(data_dir)

    def _remap(obj):
        if isinstance(obj, str) and (os.sep in obj or "/" in obj or "\\" in obj):
            # Detect the common project root and replace it
            # We look for the segment after which 'data/' appears
            parts = obj.replace("\\", "/").split("/")
            try:
                data_pos = next(i for i, p in enumerate(parts) if p == "data")
                return os.path.join(data_dir, *parts[data_pos:])
            except StopIteration:
                return obj
        if isinstance(obj, list):
            return [_remap(x) for x in obj]
        return obj

    for key, val in config["data"].items():
        config["data"][key] = _remap(val)


def main() -> None:
    args = parse_args()
    t_start = time.time()

    # Load config
    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")
    config = load_config(config_path)
    maybe_remap_paths(config, args.data_dir)

    output_dir = Path(os.path.abspath(args.output_dir))
    output_dir.mkdir(parents=True, exist_ok=True)

    k = config.get("transform", {}).get("k", 48)

    # Validate knn-device
    if args.knn_device.startswith("cuda") and not torch.cuda.is_available():
        print(f"WARNING: --knn-device={args.knn_device} but CUDA not available. Falling back to CPU.")
        args.knn_device = "cpu"

    print(f"Output directory : {output_dir}")
    print(f"Config           : {config_path}")
    print(f"k-NN k           : {k}")
    print(f"k-NN device      : {args.knn_device}")
    print(f"Graph shard size : {args.shard_size}")

    # ---------------------
    # Step 0: Manifest
    # ---------------------
    print("\n[0/5] Writing manifest.pt ...")
    write_manifest(output_dir, config, k)

    # ---------------------
    # Step 1: Enzyme shards
    # ---------------------
    if not args.skip_enzyme_shards:
        build_enzyme_shards(output_dir, config, shard_size=args.enzyme_shard_size,
                            num_workers=args.num_workers)
    else:
        print("\n[1/5] Skipping enzyme shards (--skip-enzyme-shards).")

    # ---------------------
    # Step 2: Substrate shards
    # ---------------------
    if not args.skip_substrate_shards:
        build_substrate_shards(output_dir, config, shard_size=args.substrate_shard_size,
                               num_workers=args.num_workers)
    else:
        print("\n[2/5] Skipping substrate shards (--skip-substrate-shards).")

    # ---------------------
    # Step 3: Graph shards per split
    # ---------------------
    split_cfg_keys = {
        "train": "train_data_path",
        "val":   "val_data_path",
        "test":  "test_data_path",
    }

    print("\n[3/5] Building graph shards ...")
    for split_name in args.splits:
        csv_paths: list[str] = config["data"][split_cfg_keys[split_name]]
        enz_paths: list[str] = config["data"]["enzyme_lmdb_path"]
        rxn_paths: list[str] = config["data"]["reaction_lmdb_path"]
        grv_paths: list[str] = config["data"]["grover_path"]
        str_paths: list[str] = config["data"]["structure_processed_path"]

        print(f"\n  Building split: {split_name}")

        # Load and concatenate CSVs (mirrors read_datasets in Datasets/utils.py)
        dfs = []
        dataset_ids_col = []
        for ds_id, csv_path in enumerate(csv_paths):
            if not os.path.exists(csv_path):
                print(f"  WARNING: CSV not found: {csv_path}")
                continue
            df_part = pd.read_csv(csv_path)
            df_part["dataset_id"] = ds_id
            df_part["tag"] = ds_id
            dfs.append(df_part)

        if not dfs:
            print(f"  SKIP: No CSVs found for split={split_name}")
            continue

        df = pd.concat(dfs, ignore_index=True).reset_index(drop=True)
        print(f"  Total rows: {len(df)}")

        # Compute valid_idx per dataset (intersection of sequence+structure)
        all_valid: list[int] = []
        offset = 0
        for ds_id in range(len(csv_paths)):
            ds_mask = df["dataset_id"] == ds_id
            ds_df   = df[ds_mask].reset_index(drop=True)
            ds_global_indices = list(df.index[ds_mask])

            if len(ds_df) == 0:
                offset += 0
                continue

            if ds_id >= len(enz_paths):
                print(f"  WARNING: no enzyme path for dataset {ds_id}")
                continue

            valid_local = compute_valid_idx(
                df=ds_df,
                enzyme_lmdb_path=enz_paths[ds_id],
                reaction_lmdb_path=rxn_paths[ds_id],
                grover_lmdb_path=grv_paths[ds_id],
                structure_lmdb_path=str_paths[ds_id],
            )
            # Map local indices back to global df indices
            for local_idx in valid_local:
                all_valid.append(ds_global_indices[local_idx])

        print(f"  Valid samples after filtering: {len(all_valid)}")

        build_graph_shards(
            split_name=split_name,
            df=df,
            valid_idx=all_valid,
            config=config,
            output_dir=output_dir,
            shard_size=args.shard_size,
            k=k,
            knn_device=args.knn_device,
            num_workers=args.num_workers,
        )

    # ---------------------
    # Step 4: Self-test
    # ---------------------
    print("\n[4/5] Self-test ...")
    self_test(output_dir)

    # ---------------------
    # Summary
    # ---------------------
    print_summary(output_dir, t_start)


if __name__ == "__main__":
    main()
