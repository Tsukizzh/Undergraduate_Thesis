#!/usr/bin/env python3
"""Build Q10 PT inference cache from Uni-Dock poses and ColabFold receptors.

This script creates the missing graph part of a model-readable PT cache. Dense
enzyme/substrate features are reused from q10_feature_base through symlinks.
The graph samples follow the EXP008 PtCacheDataset per-sample schema.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem.rdchem import BondType, HybridizationType


AA_NAME_SYM = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F", "GLY": "G", "HIS": "H",
    "ILE": "I", "LYS": "K", "LEU": "L", "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q",
    "ARG": "R", "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
    "UNK": "X",
}
AA_NAME_NUMBER = {k: i for i, (k, _) in enumerate(AA_NAME_SYM.items())}
BACKBONE_NAMES = {"CA", "C", "N", "O"}

PROTEIN_ATOMIC_NUMBERS = [1, 6, 7, 8, 16, 34]
LIGAND_ATOMIC_NUMBERS = [1, 6, 7, 8, 9, 15, 16, 17, 26, 35, 53]
STR_TAG_COMPLEX = 0
DATASET_ID_Q10 = 0
SUBSTRATE_ID_DIOSGENIN = 0


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    tmp_path.replace(path)


def ensure_symlink_or_copy(src: Path, dst: Path) -> str:
    if dst.exists() or dst.is_symlink():
        if dst.is_symlink() and Path(os.readlink(dst)) == src:
            return "existing_symlink"
        return "existing"
    try:
        os.symlink(src, dst, target_is_directory=src.is_dir())
        return "symlink"
    except OSError:
        if src.is_dir():
            shutil.copytree(src, dst)
            return "copied_dir"
        shutil.copy2(src, dst)
        return "copied_file"


def setup_output_cache(feature_base: Path, output_cache: Path) -> Dict[str, str]:
    output_cache.mkdir(parents=True, exist_ok=True)
    statuses = {}
    for name in ["enzymes", "substrates"]:
        statuses[name] = ensure_symlink_or_copy(feature_base / name, output_cache / name)
    manifest_src = feature_base / "manifest.pt"
    manifest_dst = output_cache / "manifest.pt"
    if not manifest_dst.exists():
        shutil.copy2(manifest_src, manifest_dst)
        statuses["manifest.pt"] = "copied_file"
    else:
        statuses["manifest.pt"] = "existing"
    return statuses


def parse_protein_atoms(pdb_path: Path) -> Dict[str, np.ndarray]:
    ptable = Chem.GetPeriodicTable()
    pos, element, aa_type, is_backbone, atom_name = [], [], [], [], []
    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[0:6].strip()
            if rec == "ENDMDL":
                break
            if rec != "ATOM":
                continue
            name = line[12:16].strip()
            res_name = line[17:20].strip().upper()
            elem_sym = line[76:78].strip().capitalize()
            if not elem_sym:
                elem_sym = line[13:14].strip().capitalize()
            try:
                atomic_num = ptable.GetAtomicNumber(elem_sym)
            except Exception:
                atomic_num = 0
            if atomic_num <= 0:
                continue
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                continue
            pos.append(xyz)
            element.append(atomic_num)
            aa_type.append(AA_NAME_NUMBER.get(res_name, AA_NAME_NUMBER["UNK"]))
            is_backbone.append(1 if name in BACKBONE_NAMES else 0)
            atom_name.append(name)
    if not pos:
        raise ValueError(f"No ATOM records parsed from {pdb_path}")
    return {
        "pos": np.asarray(pos, dtype=np.float32),
        "element": np.asarray(element, dtype=np.uint8),
        "aa_type": np.asarray(aa_type, dtype=np.uint8),
        "is_backbone": np.asarray(is_backbone, dtype=np.uint8),
        "atom_name": np.asarray(atom_name, dtype=object),
    }


def parse_ligand_sdf(sdf_path: Path) -> Dict[str, np.ndarray]:
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        raise ValueError(f"RDKit failed to parse ligand SDF: {sdf_path}")
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # The original code often used sanitize=False for docking outputs.
        pass
    try:
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
    except TypeError:
        mol_no_h = Chem.RemoveHs(mol)
    if mol_no_h is not None:
        mol = mol_no_h
    conf = mol.GetConformer()
    pos = np.asarray(conf.GetPositions(), dtype=np.float32)
    element, aromatic, degree, hybrid, atom_map = [], [], [], [], []
    hybrid_types = {t: i for i, t in enumerate(HybridizationType.names.values())}
    for atom in mol.GetAtoms():
        element.append(atom.GetAtomicNum())
        aromatic.append(1 if atom.GetIsAromatic() else 0)
        degree.append(atom.GetDegree())
        hybrid.append(hybrid_types.get(atom.GetHybridization(), 0))
        atom_map.append(atom.GetAtomMapNum())

    num_atoms = mol.GetNumAtoms()
    row, col, bond_type = [], [], []
    bond_map = {
        BondType.SINGLE: 1,
        BondType.DOUBLE: 2,
        BondType.TRIPLE: 3,
        BondType.AROMATIC: 4,
        BondType.UNSPECIFIED: 5,
    }
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row.extend([start, end])
        col.extend([end, start])
        bond_type.extend([bond_map.get(bond.GetBondType(), 5), bond_map.get(bond.GetBondType(), 5)])

    bond_index = np.asarray([row, col], dtype=np.int32)
    bond_type_arr = np.asarray(bond_type, dtype=np.uint8)
    if bond_index.size:
        perm = (bond_index[0] * num_atoms + bond_index[1]).argsort()
        bond_index = bond_index[:, perm]
        bond_type_arr = bond_type_arr[perm]

    explicit_h_counts = np.zeros(num_atoms, dtype=np.uint8)
    elem_arr = np.asarray(element, dtype=np.uint8)
    aux = np.stack([
        np.asarray(aromatic, dtype=np.uint8),
        np.asarray(degree, dtype=np.uint8),
        explicit_h_counts,
        np.asarray(hybrid, dtype=np.uint8),
    ], axis=1)
    return {
        "pos": pos,
        "element": elem_arr,
        "atom_aux": aux,
        "atom_map": np.asarray(atom_map, dtype=np.int32),
        "bond_index": bond_index,
        "bond_type": bond_type_arr,
    }


def select_pocket(protein: Dict[str, np.ndarray], ligand_pos: np.ndarray,
                  radius: float, min_atoms: int, fallback_atoms: int) -> Tuple[np.ndarray, str]:
    protein_pos = protein["pos"]
    dist = np.linalg.norm(protein_pos[:, None, :] - ligand_pos[None, :, :], axis=2).min(axis=1)
    mask = dist <= radius
    idx = np.where(mask)[0]
    method = f"radius_{radius:g}A"
    if len(idx) < min_atoms:
        order = np.argsort(dist)
        keep = min(max(fallback_atoms, min_atoms), len(order))
        idx = order[:keep]
        method = f"nearest_{keep}_fallback_from_radius_{radius:g}A"
    return np.sort(idx), method


def build_knn_graph(ligand_pos: np.ndarray, protein_pos: np.ndarray,
                    ligand_bond_index: np.ndarray, k: int) -> np.ndarray:
    pos = torch.from_numpy(np.concatenate([ligand_pos, protein_pos], axis=0).astype(np.float32))
    try:
        from torch_geometric.nn import knn_graph

        edge_index = knn_graph(pos, k=k, flow="target_to_source")
        if edge_index.is_cuda:
            edge_index = edge_index.cpu()
    except Exception:
        dist = torch.cdist(pos, pos)
        n = dist.shape[0]
        kk = min(k, max(1, n - 1))
        dist.fill_diagonal_(float("inf"))
        nbr = dist.topk(kk, largest=False).indices
        dst = torch.arange(n).repeat_interleave(kk)
        src = nbr.reshape(-1)
        edge_index = torch.stack([dst, src], dim=0)

    if ligand_bond_index is not None and ligand_bond_index.size:
        max_node = int(edge_index.max().item()) + 1 if edge_index.numel() else 1
        knn_hash = edge_index[0] * max_node + edge_index[1]
        bond_t = torch.from_numpy(ligand_bond_index.astype(np.int64))
        bond_hash = bond_t[0] * max_node + bond_t[1]
        edge_index = edge_index[:, ~torch.isin(knn_hash, bond_hash)]
    return edge_index.numpy().astype(np.int32)


def build_sample(row: Dict[str, str], enzyme_id: int, sample_id: int,
                 radius: float, min_atoms: int, fallback_atoms: int, k: int) -> Tuple[Dict[str, torch.Tensor], Dict[str, object]]:
    receptor_pdb = Path(row["receptor_pdb_path"])
    pose_sdf = Path(row["pose_sdf_path"])
    ligand = parse_ligand_sdf(pose_sdf)
    protein = parse_protein_atoms(receptor_pdb)
    pocket_idx, pocket_method = select_pocket(protein, ligand["pos"], radius, min_atoms, fallback_atoms)

    protein_pos = protein["pos"][pocket_idx]
    protein_element = protein["element"][pocket_idx]
    protein_aa_type = protein["aa_type"][pocket_idx]
    protein_is_backbone = protein["is_backbone"][pocket_idx]

    knn_edge_index = build_knn_graph(ligand["pos"], protein_pos, ligand["bond_index"], k)

    sample = {
        "ligand_pos": torch.from_numpy(ligand["pos"]).float(),
        "ligand_element": torch.from_numpy(ligand["element"]).long(),
        "ligand_atom_aux": torch.from_numpy(ligand["atom_aux"]).long(),
        "ligand_index_raw": torch.from_numpy(ligand["atom_map"]).long(),
        "protein_pos": torch.from_numpy(protein_pos).float(),
        "protein_element": torch.from_numpy(protein_element).long(),
        "protein_aa_type": torch.from_numpy(protein_aa_type).long(),
        "protein_is_backbone": torch.from_numpy(protein_is_backbone).long(),
        "bond_index": torch.from_numpy(ligand["bond_index"]).long(),
        "bond_type": torch.from_numpy(ligand["bond_type"]).long(),
        "knn_edge_index": torch.from_numpy(knn_edge_index).long(),
        "enzyme_id": torch.tensor(enzyme_id, dtype=torch.long),
        "substrate_id": torch.tensor(SUBSTRATE_ID_DIOSGENIN, dtype=torch.long),
        "dataset_id": torch.tensor(DATASET_ID_Q10, dtype=torch.long),
        "label": torch.tensor(0, dtype=torch.long),
        "str_tag_code": torch.tensor(STR_TAG_COMPLEX, dtype=torch.long),
        "sample_weight": torch.tensor(1.0, dtype=torch.float32),
    }
    meta = {
        "sample_id": sample_id,
        "candidate_id": row["candidate_id"],
        "input_list": row["input_list"],
        "enzyme_id": enzyme_id,
        "substrate_id": SUBSTRATE_ID_DIOSGENIN,
        "docking_score": row.get("docking_score", ""),
        "box_source": row.get("box_source", ""),
        "heme_fe_quality": row.get("heme_fe_quality", ""),
        "nearest_cys_sg_fe_distance": row.get("nearest_cys_sg_fe_distance", ""),
        "template_pdb_id": row.get("template_pdb_id", ""),
        "fit_rmsd": row.get("fit_rmsd", ""),
        "pocket_method": pocket_method,
        "n_ligand_atoms": int(ligand["pos"].shape[0]),
        "n_protein_atoms_full": int(protein["pos"].shape[0]),
        "n_protein_atoms_pocket": int(protein_pos.shape[0]),
        "n_bond_edges": int(ligand["bond_index"].shape[1]),
        "n_knn_edges": int(knn_edge_index.shape[1]),
        "receptor_pdb_path": str(receptor_pdb),
        "pose_sdf_path": str(pose_sdf),
    }
    return sample, meta


def load_candidate_to_enzyme_id(q10_dir: Path) -> Dict[str, int]:
    paths = [
        q10_dir / "features" / "enzyme_features" / "q10_candidate_enzymes_for_esm.csv",
        q10_dir / "structures" / "colabfold_inputs" / "q10_colabfold_input_manifest.csv",
    ]
    mapping: Dict[str, int] = {}
    for path in paths:
        if not path.exists():
            continue
        for row in read_csv(path):
            cid = row.get("candidate_id", "")
            gid = row.get("global_enzyme_id", "")
            if cid != "" and gid != "":
                mapping[cid] = int(gid)
        if mapping:
            break
    if not mapping:
        raise FileNotFoundError("Could not build candidate_id -> global_enzyme_id mapping")
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q10-dir", required=True)
    parser.add_argument("--output-cache-dir", default="")
    parser.add_argument("--split", default="test")
    parser.add_argument("--radius", type=float, default=10.0)
    parser.add_argument("--min-pocket-atoms", type=int, default=50)
    parser.add_argument("--nearest-fallback-atoms", type=int, default=300)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    q10_dir = Path(args.q10_dir)
    feature_base = q10_dir / "model_cache" / "q10_feature_base"
    output_cache = Path(args.output_cache_dir) if args.output_cache_dir else q10_dir / "model_cache" / "q10_unidock_pt_msa_m1r1_v1"
    split_dir = output_cache / args.split
    samples_dir = split_dir / "samples"
    samples_dir.mkdir(parents=True, exist_ok=True)

    setup_status = setup_output_cache(feature_base, output_cache)
    docking_results = q10_dir / "docking" / "manifests" / "q10_unidock_results_msa_m1r1_v1.csv"
    rows = [
        row for row in read_csv(docking_results)
        if row.get("docking_status", "").startswith("complete")
        and row.get("pose_sdf_path")
        and Path(row["pose_sdf_path"]).exists()
        and Path(row["receptor_pdb_path"]).exists()
    ]
    rows.sort(key=lambda r: r["candidate_id"])
    if args.limit:
        rows = rows[: args.limit]
    if not rows:
        raise RuntimeError(f"No complete docking rows found in {docking_results}")

    cid_to_gid = load_candidate_to_enzyme_id(q10_dir)
    index = {
        "sample_ids": [],
        "enzyme_ids": [],
        "substrate_ids": [],
        "graph_shards": [],
        "graph_rows": [],
    }
    manifest_rows: List[Dict[str, object]] = []
    failures: List[Dict[str, object]] = []

    for sample_id, row in enumerate(rows):
        candidate_id = row["candidate_id"]
        if candidate_id not in cid_to_gid:
            failures.append({"candidate_id": candidate_id, "reason": "missing_global_enzyme_id"})
            continue
        try:
            sample, meta = build_sample(
                row=row,
                enzyme_id=cid_to_gid[candidate_id],
                sample_id=sample_id,
                radius=args.radius,
                min_atoms=args.min_pocket_atoms,
                fallback_atoms=args.nearest_fallback_atoms,
                k=args.k,
            )
        except Exception as exc:
            failures.append({"candidate_id": candidate_id, "reason": repr(exc)})
            continue

        out_dir = samples_dir / f"{sample_id // 1000:03d}"
        out_dir.mkdir(parents=True, exist_ok=True)
        torch.save(sample, out_dir / f"sample_{sample_id:06d}.pt")
        index["sample_ids"].append(sample_id)
        index["enzyme_ids"].append(cid_to_gid[candidate_id])
        index["substrate_ids"].append(SUBSTRATE_ID_DIOSGENIN)
        index["graph_shards"].append(0)
        index["graph_rows"].append(sample_id)
        manifest_rows.append(meta)

    index_t = {
        "sample_ids": torch.tensor(index["sample_ids"], dtype=torch.int32),
        "enzyme_ids": torch.tensor(index["enzyme_ids"], dtype=torch.int32),
        "substrate_ids": torch.tensor(index["substrate_ids"], dtype=torch.int32),
        "graph_shards": torch.tensor(index["graph_shards"], dtype=torch.int16),
        "graph_rows": torch.tensor(index["graph_rows"], dtype=torch.int32),
    }
    torch.save(index_t, split_dir / "index.pt")

    manifest_csv = output_cache / "manifests" / f"q10_model_cache_manifest_{args.split}_msa_m1r1_v1.csv"
    fieldnames = [
        "sample_id", "candidate_id", "input_list", "enzyme_id", "substrate_id",
        "docking_score", "box_source", "heme_fe_quality", "nearest_cys_sg_fe_distance",
        "template_pdb_id", "fit_rmsd", "pocket_method", "n_ligand_atoms",
        "n_protein_atoms_full", "n_protein_atoms_pocket", "n_bond_edges",
        "n_knn_edges", "receptor_pdb_path", "pose_sdf_path",
    ]
    write_csv(manifest_csv, manifest_rows, fieldnames)
    failure_csv = output_cache / "manifests" / f"q10_model_cache_failures_{args.split}_msa_m1r1_v1.csv"
    write_csv(failure_csv, failures, ["candidate_id", "reason"]) if failures else None

    summary = {
        "q10_dir": str(q10_dir),
        "feature_base": str(feature_base),
        "output_cache": str(output_cache),
        "split": args.split,
        "docking_rows_considered": len(rows),
        "samples_written": len(manifest_rows),
        "failures": len(failures),
        "setup_status": setup_status,
        "manifest_csv": str(manifest_csv),
        "failure_csv": str(failure_csv) if failures else "",
        "radius": args.radius,
        "k": args.k,
    }
    summary_path = output_cache / "manifests" / f"q10_model_cache_summary_{args.split}_msa_m1r1_v1.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
