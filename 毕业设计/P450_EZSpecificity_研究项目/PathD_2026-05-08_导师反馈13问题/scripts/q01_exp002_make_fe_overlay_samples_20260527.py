#!/usr/bin/env python3
"""Build EXP002 Fe/HEM overlay sample files on top of EXP001 PT cache.

This script keeps the EXP001 cache as the base layout. The intended workflow is:

1. Create the EXP002 cache directory as a hard-linked copy of the EXP001 cache.
2. Run this script to replace only the samples whose brenda Dock Index can be
   aligned to an audited RCSB HEM/Fe structure.

The replacement is atomic: each enhanced sample is first written to a temporary
file and then moved over the hard link, so the original EXP001 cache is not
modified.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import pickle
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from Bio.Align import PairwiseAligner
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from torch_geometric.nn import knn_graph


HEME_COMPONENTS = {"HEM", "HEC"}
ELEMENT_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "S": 16,
    "SE": 34,
    "FE": 26,
}


@dataclass
class AlignmentResult:
    dock_index: int
    entry_id: str
    official_chain: str
    cif_label_chain: str
    cif_auth_chain: str
    n_pairs: int
    rmsd: float
    score: float
    rotation: np.ndarray
    translation: np.ndarray
    heme_atoms: list[dict]


def aa1(resname: str) -> str:
    return protein_letters_3to1_extended.get(str(resname).upper(), "X")


def as_list(value):
    return value if isinstance(value, list) else [value]


def read_official_pdb_ca(path: Path) -> dict[str, list[tuple[str, np.ndarray]]]:
    chains: dict[str, list[tuple[str, np.ndarray]]] = {}
    seen = set()
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21:22].strip() or "?"
            resid = (chain, line[22:26].strip(), line[26:27].strip())
            if resid in seen:
                continue
            seen.add(resid)
            resname = line[17:20].strip()
            xyz = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=np.float64,
            )
            chains.setdefault(chain, []).append((aa1(resname), xyz))
    return chains


def read_cif_ca_and_heme(path: Path) -> tuple[dict[str, list[tuple[str, np.ndarray, str]]], list[dict]]:
    mm = MMCIF2Dict(str(path))

    group = as_list(mm["_atom_site.group_PDB"])
    comp = as_list(mm["_atom_site.label_comp_id"])
    label_chain = as_list(mm["_atom_site.label_asym_id"])
    auth_chain = as_list(mm["_atom_site.auth_asym_id"])
    atom_name = as_list(mm["_atom_site.label_atom_id"])
    element = as_list(mm["_atom_site.type_symbol"])
    auth_seq = as_list(mm["_atom_site.auth_seq_id"])
    ins_code = as_list(mm.get("_atom_site.pdbx_PDB_ins_code", ["?"] * len(group)))
    x_vals = as_list(mm["_atom_site.Cartn_x"])
    y_vals = as_list(mm["_atom_site.Cartn_y"])
    z_vals = as_list(mm["_atom_site.Cartn_z"])

    ca: dict[str, list[tuple[str, np.ndarray, str]]] = {}
    heme_atoms: list[dict] = []
    seen_ca = set()

    for i in range(len(group)):
        xyz = np.array([float(x_vals[i]), float(y_vals[i]), float(z_vals[i])], dtype=np.float64)
        label = str(label_chain[i])
        auth = str(auth_chain[i])
        comp_id = str(comp[i]).upper()
        elem = str(element[i]).strip().upper()

        if group[i] == "ATOM" and str(atom_name[i]).strip() == "CA":
            key = (label, auth_seq[i], ins_code[i])
            if key not in seen_ca:
                seen_ca.add(key)
                ca.setdefault(label, []).append((aa1(comp_id), xyz, auth))

        if group[i] == "HETATM" and (comp_id in HEME_COMPONENTS or elem == "FE" or comp_id == "FE"):
            heme_atoms.append(
                {
                    "comp_id": comp_id,
                    "atom_name": str(atom_name[i]).strip(),
                    "element": elem,
                    "label_chain": label,
                    "auth_chain": auth,
                    "xyz": xyz,
                }
            )

    return ca, heme_atoms


def sequence_align_pairs(a: list[tuple[str, np.ndarray]], b: list[tuple[str, np.ndarray, str]]):
    seq_a = "".join(item[0] for item in a)
    seq_b = "".join(item[0] for item in b)
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    best = aligner.align(seq_a, seq_b)[0]
    pairs = []
    for (a0, a1), (b0, b1) in zip(*best.aligned):
        for ia, ib in zip(range(a0, a1), range(b0, b1)):
            if a[ia][0] != "X" and b[ib][0] != "X":
                pairs.append((ia, ib))
    return pairs, float(best.score)


def fit_one_structure(
    dock_index: int,
    pdb_path: Path,
    cif_path: Path,
    min_pairs: int,
    max_rmsd: float,
) -> AlignmentResult:
    official_ca = read_official_pdb_ca(pdb_path)
    cif_ca, all_heme_atoms = read_cif_ca_and_heme(cif_path)
    if not official_ca:
        raise RuntimeError("official PDB has no CA atoms")
    if not cif_ca:
        raise RuntimeError("RCSB mmCIF has no CA atoms")
    if not all_heme_atoms:
        raise RuntimeError("RCSB mmCIF has no HEM/HEC/Fe atoms")

    best = None
    for official_chain, official_records in official_ca.items():
        if len(official_records) < min_pairs:
            continue
        for cif_label_chain, cif_records in cif_ca.items():
            if len(cif_records) < min_pairs:
                continue
            pairs, score = sequence_align_pairs(official_records, cif_records)
            if len(pairs) < min_pairs:
                continue

            official_xyz = np.array([official_records[i][1] for i, _ in pairs], dtype=np.float64)
            cif_xyz = np.array([cif_records[j][1] for _, j in pairs], dtype=np.float64)
            sup = SVDSuperimposer()
            sup.set(official_xyz, cif_xyz)
            sup.run()
            rmsd = float(sup.get_rms())
            rot, tran = sup.get_rotran()
            auth_values = sorted({cif_records[j][2] for _, j in pairs})
            candidate = {
                "official_chain": official_chain,
                "cif_label_chain": cif_label_chain,
                "cif_auth_chain": auth_values[0] if auth_values else "",
                "n_pairs": len(pairs),
                "score": score,
                "rmsd": rmsd,
                "rotation": rot,
                "translation": tran,
                "auth_values": set(auth_values),
            }
            if best is None or (candidate["n_pairs"], -candidate["rmsd"]) > (best["n_pairs"], -best["rmsd"]):
                best = candidate

    if best is None:
        raise RuntimeError("no valid CA alignment")
    if best["rmsd"] > max_rmsd:
        raise RuntimeError(f"alignment RMSD too high: {best['rmsd']:.3f}")

    selected_atoms = [
        atom for atom in all_heme_atoms
        if atom["auth_chain"] in best["auth_values"]
    ]
    if not selected_atoms:
        selected_atoms = all_heme_atoms

    transformed = []
    for atom in selected_atoms:
        elem = atom["element"].upper()
        atomic_number = ELEMENT_NUMBERS.get(elem)
        if atomic_number is None:
            continue
        xyz = atom["xyz"] @ best["rotation"] + best["translation"]
        transformed.append({**atom, "atomic_number": atomic_number, "xyz_transformed": xyz.astype(np.float32)})

    if not transformed:
        raise RuntimeError("no transformable HEM/Fe atoms")

    return AlignmentResult(
        dock_index=dock_index,
        entry_id=cif_path.stem.upper(),
        official_chain=best["official_chain"],
        cif_label_chain=best["cif_label_chain"],
        cif_auth_chain=best["cif_auth_chain"],
        n_pairs=int(best["n_pairs"]),
        rmsd=float(best["rmsd"]),
        score=float(best["score"]),
        rotation=best["rotation"],
        translation=best["translation"],
        heme_atoms=transformed,
    )


def recompute_knn(ligand_pos: torch.Tensor, protein_pos: torch.Tensor, k: int) -> torch.Tensor:
    pos = torch.cat([ligand_pos.float(), protein_pos.float()], dim=0)
    with torch.no_grad():
        edge_index = knn_graph(pos, k=k, flow="target_to_source")
    return edge_index.to(torch.int32).cpu()


def sample_path(cache_dir: Path, split: str, sample_id: int) -> Path:
    return cache_dir / split / "samples" / f"{sample_id // 1000:03d}" / f"sample_{sample_id:06d}.pt"


def write_atomic_sample(path: Path, sample: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        torch.save(sample, tmp_path)
        os.replace(tmp_path, path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def update_manifest(cache_dir: Path) -> None:
    manifest_path = cache_dir / "manifest.pt"
    manifest = torch.load(manifest_path, map_location="cpu", weights_only=False)
    manifest["protein_elements"] = [1, 6, 7, 8, 16, 34, 26]
    manifest["protein_num_aa"] = 22
    manifest["protein_has_is_hetero"] = True
    manifest["exp002_fe_overlay"] = True
    fd, tmp_name = tempfile.mkstemp(prefix="manifest.", suffix=".tmp", dir=str(cache_dir))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        torch.save(manifest, tmp_path)
        os.replace(tmp_path, manifest_path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def iter_target_rows(audit_csv: Path, google_manifest_csv: Path) -> pd.DataFrame:
    audit = pd.read_csv(audit_csv)
    gd = pd.read_csv(google_manifest_csv)
    df = audit.merge(
        gd[["structure_index", "found_in_drive_metadata"]],
        on="structure_index",
        how="left",
    )
    mask = (
        df["exp002_atom_audit_usable_target_heme_fe"].fillna(False).astype(bool)
        & df["found_in_drive_metadata"].fillna(False).astype(bool)
    )
    return df.loc[mask].copy()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-cache", required=True)
    parser.add_argument("--out-cache", required=True)
    parser.add_argument("--normalized-splits", required=True)
    parser.add_argument("--audit-csv", required=True)
    parser.add_argument("--google-manifest-csv", required=True)
    parser.add_argument("--official-pdb-dir", required=True)
    parser.add_argument("--raw-mmcif-dir", required=True)
    parser.add_argument("--report-dir", required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--min-pairs", type=int, default=80)
    parser.add_argument("--max-rmsd", type=float, default=5.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--no-manifest-update", action="store_true")
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args()

    base_cache = Path(args.base_cache)
    out_cache = Path(args.out_cache)
    split_root = Path(args.normalized_splits)
    official_pdb_dir = Path(args.official_pdb_dir)
    raw_mmcif_dir = Path(args.raw_mmcif_dir)
    report_dir = Path(args.report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    if not out_cache.exists():
        raise SystemExit(f"out cache does not exist: {out_cache}")
    if not (out_cache / "train").exists():
        raise SystemExit("out cache must be initialized from base cache before running this script")

    target_rows = iter_target_rows(Path(args.audit_csv), Path(args.google_manifest_csv))
    target_by_dock = {
        int(row.structure_index): row
        for row in target_rows.itertuples(index=False)
    }

    plan = []
    for split in ["train", "val", "test"]:
        split_csv = split_root / "brenda" / f"{split}.csv"
        df = pd.read_csv(split_csv, usecols=["Dock Index"])
        index = torch.load(out_cache / split / "index.pt", map_location="cpu", weights_only=False)
        valid_sample_ids = set(int(x) for x in index["sample_ids"].tolist())
        for sample_id, dock_index in enumerate(df["Dock Index"].tolist()):
            if sample_id not in valid_sample_ids:
                continue
            dock_index = int(dock_index)
            if dock_index in target_by_dock:
                plan.append((split, sample_id, dock_index))

    if args.shard_count < 1:
        raise SystemExit("--shard-count must be >= 1")
    if not (0 <= args.shard_index < args.shard_count):
        raise SystemExit("--shard-index must satisfy 0 <= shard-index < shard-count")

    if args.shard_count > 1:
        plan = [
            item for i, item in enumerate(plan)
            if i % args.shard_count == args.shard_index
        ]

    if args.limit:
        plan = plan[: args.limit]

    if not args.no_manifest_update:
        update_manifest(out_cache)

    audit_path = report_dir / "fe_overlay_sample_audit.csv"
    fieldnames = [
        "split", "sample_id", "dock_index", "status", "entry_id",
        "n_pairs", "rmsd", "n_heme_atoms", "n_fe_atoms", "old_n_protein",
        "new_n_protein", "error",
    ]

    done_keys = set()
    if audit_path.exists():
        with audit_path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("status") == "ok":
                    done_keys.add((row["split"], int(row["sample_id"]), int(row["dock_index"])))

    start = time.time()
    counts = {"ok": 0, "skip_done": 0, "failed": 0, "missing_base_sample": 0}
    write_header = not audit_path.exists()
    with audit_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        for i, (split, sample_id, dock_index) in enumerate(plan, start=1):
            key = (split, sample_id, dock_index)
            if key in done_keys:
                counts["skip_done"] += 1
                continue

            base_sample_path = sample_path(base_cache, split, sample_id)
            out_sample_path = sample_path(out_cache, split, sample_id)
            row = target_by_dock[dock_index]

            record = {
                "split": split,
                "sample_id": sample_id,
                "dock_index": dock_index,
                "status": "",
                "entry_id": "",
                "n_pairs": "",
                "rmsd": "",
                "n_heme_atoms": "",
                "n_fe_atoms": "",
                "old_n_protein": "",
                "new_n_protein": "",
                "error": "",
            }

            if not base_sample_path.exists():
                record["status"] = "missing_base_sample"
                counts["missing_base_sample"] += 1
                writer.writerow(record)
                continue

            try:
                pdb_path = official_pdb_dir / f"{dock_index}.pdb"
                cif_path = raw_mmcif_dir / str(row.best_target_file)
                aln = fit_one_structure(dock_index, pdb_path, cif_path, args.min_pairs, args.max_rmsd)
                sample = torch.load(base_sample_path, map_location="cpu", weights_only=False)
                old_n = int(sample["protein_pos"].shape[0])

                heme_pos = torch.tensor(
                    np.stack([a["xyz_transformed"] for a in aln.heme_atoms], axis=0),
                    dtype=torch.float32,
                )
                heme_element = torch.tensor([a["atomic_number"] for a in aln.heme_atoms], dtype=torch.uint8)
                heme_aa = torch.full((len(aln.heme_atoms),), 21, dtype=torch.uint8)
                heme_backbone = torch.zeros(len(aln.heme_atoms), dtype=torch.uint8)
                old_hetero = sample.get(
                    "protein_is_hetero",
                    torch.zeros(old_n, dtype=torch.uint8),
                ).to(torch.uint8)
                heme_hetero = torch.ones(len(aln.heme_atoms), dtype=torch.uint8)

                sample["protein_pos"] = torch.cat([sample["protein_pos"].float(), heme_pos], dim=0)
                sample["protein_element"] = torch.cat([sample["protein_element"].to(torch.uint8), heme_element], dim=0)
                sample["protein_aa_type"] = torch.cat([sample["protein_aa_type"].to(torch.uint8), heme_aa], dim=0)
                sample["protein_is_backbone"] = torch.cat(
                    [sample["protein_is_backbone"].to(torch.uint8), heme_backbone], dim=0
                )
                sample["protein_is_hetero"] = torch.cat([old_hetero, heme_hetero], dim=0)
                sample["knn_edge_index"] = recompute_knn(
                    sample["ligand_pos"], sample["protein_pos"], args.k
                )

                write_atomic_sample(out_sample_path, sample)

                record.update(
                    {
                        "status": "ok",
                        "entry_id": aln.entry_id,
                        "n_pairs": aln.n_pairs,
                        "rmsd": f"{aln.rmsd:.4f}",
                        "n_heme_atoms": len(aln.heme_atoms),
                        "n_fe_atoms": sum(1 for a in aln.heme_atoms if a["atomic_number"] == 26),
                        "old_n_protein": old_n,
                        "new_n_protein": int(sample["protein_pos"].shape[0]),
                    }
                )
                counts["ok"] += 1
            except Exception as exc:
                record["status"] = "failed"
                record["error"] = repr(exc)
                counts["failed"] += 1

            writer.writerow(record)
            if i % args.progress_every == 0:
                elapsed = time.time() - start
                print(f"[{i}/{len(plan)}] {counts} elapsed={elapsed:.1f}s", flush=True)

    summary = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "base_cache": str(base_cache),
        "out_cache": str(out_cache),
        "planned_samples": len(plan),
        "shard_count": args.shard_count,
        "shard_index": args.shard_index,
        "counts": counts,
        "audit_csv": str(audit_path),
        "elapsed_sec": round(time.time() - start, 3),
    }
    summary_path = report_dir / "fe_overlay_sample_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if counts["failed"] == 0 and counts["missing_base_sample"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
