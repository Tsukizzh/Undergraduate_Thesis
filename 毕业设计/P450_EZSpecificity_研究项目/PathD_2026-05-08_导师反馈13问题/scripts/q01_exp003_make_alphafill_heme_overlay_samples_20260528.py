#!/usr/bin/env python3
"""Build EXP003 full-cache P450 HEM/Fe overlay samples using AlphaFill.

The script expects the output cache to be a hard-linked copy of EXP001.
It replaces only trainable 389-P450 samples whose official ESIBank PDB
exists and whose UniProt has an audited AlphaFill HEM transplant candidate.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer

from q01_exp002_make_fe_overlay_samples_20260527 import (
    ELEMENT_NUMBERS,
    aa1,
    as_list,
    read_official_pdb_ca,
    recompute_knn,
    sample_path,
    sequence_align_pairs,
    update_manifest,
    write_atomic_sample,
)


HEME_COMPONENTS = {"HEM", "HEC"}


@dataclass
class AlphaFillFit:
    uniprot: str
    dock_index: int
    official_chain: str
    alphafill_chain: str
    heme_asym_id: str
    n_pairs: int
    rmsd: float
    score: float
    heme_atoms: list[dict]


def read_alphafill_ca_and_heme(path: Path, heme_asym_id: str) -> tuple[dict[str, list[tuple[str, np.ndarray, str]]], list[dict]]:
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
        label = str(label_chain[i])
        auth = str(auth_chain[i])
        comp_id = str(comp[i]).upper()
        elem = str(element[i]).strip().upper()
        xyz = np.array([float(x_vals[i]), float(y_vals[i]), float(z_vals[i])], dtype=np.float64)

        if group[i] == "ATOM" and str(atom_name[i]).strip() == "CA":
            key = (label, auth_seq[i], ins_code[i])
            if key not in seen_ca:
                seen_ca.add(key)
                ca.setdefault(label, []).append((aa1(comp_id), xyz, auth))

        is_selected_heme = (
            group[i] == "HETATM"
            and (label == heme_asym_id or auth == heme_asym_id)
            and (comp_id in HEME_COMPONENTS or elem == "FE" or comp_id == "FE")
        )
        if is_selected_heme:
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


def fit_alphafill_heme(
    uniprot: str,
    dock_index: int,
    pdb_path: Path,
    cif_path: Path,
    heme_asym_id: str,
    min_pairs: int,
    max_rmsd: float,
) -> AlphaFillFit:
    official_ca = read_official_pdb_ca(pdb_path)
    alphafill_ca, heme_atoms = read_alphafill_ca_and_heme(cif_path, heme_asym_id)
    if not official_ca:
        raise RuntimeError("official PDB has no CA atoms")
    if not alphafill_ca:
        raise RuntimeError("AlphaFill CIF has no CA atoms")
    if not heme_atoms:
        raise RuntimeError(f"AlphaFill CIF has no selected HEM/Fe atoms for asym_id={heme_asym_id}")

    best = None
    for official_chain, official_records in official_ca.items():
        if len(official_records) < min_pairs:
            continue
        for af_chain, af_records in alphafill_ca.items():
            if len(af_records) < min_pairs:
                continue
            pairs, score = sequence_align_pairs(official_records, af_records)
            if len(pairs) < min_pairs:
                continue
            official_xyz = np.array([official_records[i][1] for i, _ in pairs], dtype=np.float64)
            af_xyz = np.array([af_records[j][1] for _, j in pairs], dtype=np.float64)
            sup = SVDSuperimposer()
            sup.set(official_xyz, af_xyz)
            sup.run()
            rmsd = float(sup.get_rms())
            rot, tran = sup.get_rotran()
            candidate = {
                "official_chain": official_chain,
                "af_chain": af_chain,
                "n_pairs": len(pairs),
                "score": float(score),
                "rmsd": rmsd,
                "rotation": rot,
                "translation": tran,
            }
            if best is None or (candidate["n_pairs"], -candidate["rmsd"]) > (best["n_pairs"], -best["rmsd"]):
                best = candidate

    if best is None:
        raise RuntimeError("no valid CA alignment")
    if best["rmsd"] > max_rmsd:
        raise RuntimeError(f"alignment RMSD too high: {best['rmsd']:.3f}")

    transformed = []
    for atom in heme_atoms:
        elem = atom["element"].upper()
        atomic_number = ELEMENT_NUMBERS.get(elem)
        if atomic_number is None:
            continue
        xyz = atom["xyz"] @ best["rotation"] + best["translation"]
        transformed.append({**atom, "atomic_number": atomic_number, "xyz_transformed": xyz.astype(np.float32)})
    if not transformed:
        raise RuntimeError("no transformable HEM/Fe atoms")
    if sum(1 for atom in transformed if atom["atomic_number"] == 26) < 1:
        raise RuntimeError("selected HEM transplant has no Fe atom")

    return AlphaFillFit(
        uniprot=uniprot,
        dock_index=dock_index,
        official_chain=best["official_chain"],
        alphafill_chain=best["af_chain"],
        heme_asym_id=heme_asym_id,
        n_pairs=int(best["n_pairs"]),
        rmsd=float(best["rmsd"]),
        score=float(best["score"]),
        heme_atoms=transformed,
    )


def load_plan(args) -> list[dict]:
    trainable = pd.read_csv(args.trainable_manifest)
    official = pd.read_csv(args.official_pdb_manifest)
    alpha = pd.read_csv(args.alphafill_summary)

    official = official.rename(columns={"structure_index": "dock_index"})
    official["dock_index"] = official["dock_index"].astype(int)
    alpha = alpha.rename(columns={"best_asym_id": "alphafill_heme_asym_id"})
    alpha["has_alphafill_heme"] = alpha["hem_candidates"].fillna(0).astype(int) > 0

    df = trainable.merge(
        official[["dock_index", "found_in_drive_metadata", "found_in_drivefs_metadata"]],
        on="dock_index",
        how="left",
    )
    df = df.merge(
        alpha[
            [
                "uniprot",
                "has_alphafill_heme",
                "alphafill_heme_asym_id",
                "best_pdb_id",
                "best_identity",
                "best_length",
                "best_global_rmsd",
                "best_local_rmsd",
                "best_clash_score",
                "best_atom_count",
            ]
        ],
        on="uniprot",
        how="left",
    )
    found = df["found_in_drive_metadata"].fillna(False).astype(bool) | df["found_in_drivefs_metadata"].fillna(False).astype(bool)
    mask = found & df["has_alphafill_heme"].fillna(False).astype(bool)
    plan = df.loc[mask].copy()
    plan["sample_id"] = plan["sample_id"].astype(int)
    plan["dock_index"] = plan["dock_index"].astype(int)
    return plan.to_dict("records")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-cache", required=True)
    parser.add_argument("--out-cache", required=True)
    parser.add_argument("--trainable-manifest", required=True)
    parser.add_argument("--official-pdb-manifest", required=True)
    parser.add_argument("--official-pdb-dir", required=True)
    parser.add_argument("--alphafill-summary", required=True)
    parser.add_argument("--alphafill-cif-dir", required=True)
    parser.add_argument("--report-dir", required=True)
    parser.add_argument("--k", type=int, default=48)
    parser.add_argument("--min-pairs", type=int, default=80)
    parser.add_argument("--max-rmsd", type=float, default=6.0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--no-manifest-update", action="store_true")
    parser.add_argument("--progress-every", type=int, default=100)
    args = parser.parse_args()

    base_cache = Path(args.base_cache)
    out_cache = Path(args.out_cache)
    official_pdb_dir = Path(args.official_pdb_dir)
    alphafill_cif_dir = Path(args.alphafill_cif_dir)
    report_dir = Path(args.report_dir)
    report_dir.mkdir(parents=True, exist_ok=True)

    if not out_cache.exists() or not (out_cache / "train").exists():
        raise SystemExit("out cache must be initialized from base cache before running this script")

    plan = load_plan(args)
    if args.shard_count < 1:
        raise SystemExit("--shard-count must be >= 1")
    if not (0 <= args.shard_index < args.shard_count):
        raise SystemExit("--shard-index must satisfy 0 <= shard-index < shard-count")
    if args.shard_count > 1:
        plan = [item for i, item in enumerate(plan) if i % args.shard_count == args.shard_index]
    if args.limit:
        plan = plan[: args.limit]

    if not args.no_manifest_update:
        update_manifest(out_cache)

    audit_path = report_dir / "alphafill_heme_overlay_sample_audit.csv"
    fieldnames = [
        "split", "sample_id", "dock_index", "uniprot", "status", "alphafill_template_pdb",
        "alphafill_heme_asym_id", "n_pairs", "rmsd", "n_heme_atoms", "n_fe_atoms",
        "old_n_protein", "new_n_protein", "best_identity", "best_local_rmsd", "best_clash_score", "error",
    ]
    done_keys = set()
    if audit_path.exists():
        with audit_path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("status") == "ok":
                    done_keys.add((row["split"], int(row["sample_id"]), int(row["dock_index"])))

    counts = {"ok": 0, "skip_done": 0, "failed": 0, "missing_base_sample": 0}
    start = time.time()
    write_header = not audit_path.exists()
    with audit_path.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        for i, row in enumerate(plan, start=1):
            split = row["split"]
            sample_id = int(row["sample_id"])
            dock_index = int(row["dock_index"])
            uniprot = str(row["uniprot"])
            key = (split, sample_id, dock_index)
            if key in done_keys:
                counts["skip_done"] += 1
                continue

            record = {
                "split": split,
                "sample_id": sample_id,
                "dock_index": dock_index,
                "uniprot": uniprot,
                "status": "",
                "alphafill_template_pdb": row.get("best_pdb_id", ""),
                "alphafill_heme_asym_id": row.get("alphafill_heme_asym_id", ""),
                "n_pairs": "",
                "rmsd": "",
                "n_heme_atoms": "",
                "n_fe_atoms": "",
                "old_n_protein": "",
                "new_n_protein": "",
                "best_identity": row.get("best_identity", ""),
                "best_local_rmsd": row.get("best_local_rmsd", ""),
                "best_clash_score": row.get("best_clash_score", ""),
                "error": "",
            }
            base_sample_path = sample_path(base_cache, split, sample_id)
            out_sample_path = sample_path(out_cache, split, sample_id)
            if not base_sample_path.exists():
                record["status"] = "missing_base_sample"
                counts["missing_base_sample"] += 1
                writer.writerow(record)
                continue

            try:
                fit = fit_alphafill_heme(
                    uniprot=uniprot,
                    dock_index=dock_index,
                    pdb_path=official_pdb_dir / f"{dock_index}.pdb",
                    cif_path=alphafill_cif_dir / f"{uniprot}.cif",
                    heme_asym_id=str(row["alphafill_heme_asym_id"]),
                    min_pairs=args.min_pairs,
                    max_rmsd=args.max_rmsd,
                )
                sample = torch.load(base_sample_path, map_location="cpu", weights_only=False)
                old_n = int(sample["protein_pos"].shape[0])
                heme_pos = torch.tensor(np.stack([a["xyz_transformed"] for a in fit.heme_atoms], axis=0), dtype=torch.float32)
                heme_element = torch.tensor([a["atomic_number"] for a in fit.heme_atoms], dtype=torch.uint8)
                heme_aa = torch.full((len(fit.heme_atoms),), 21, dtype=torch.uint8)
                heme_backbone = torch.zeros(len(fit.heme_atoms), dtype=torch.uint8)
                old_hetero = sample.get("protein_is_hetero", torch.zeros(old_n, dtype=torch.uint8)).to(torch.uint8)
                heme_hetero = torch.ones(len(fit.heme_atoms), dtype=torch.uint8)

                sample["protein_pos"] = torch.cat([sample["protein_pos"].float(), heme_pos], dim=0)
                sample["protein_element"] = torch.cat([sample["protein_element"].to(torch.uint8), heme_element], dim=0)
                sample["protein_aa_type"] = torch.cat([sample["protein_aa_type"].to(torch.uint8), heme_aa], dim=0)
                sample["protein_is_backbone"] = torch.cat([sample["protein_is_backbone"].to(torch.uint8), heme_backbone], dim=0)
                sample["protein_is_hetero"] = torch.cat([old_hetero, heme_hetero], dim=0)
                sample["knn_edge_index"] = recompute_knn(sample["ligand_pos"], sample["protein_pos"], args.k)
                write_atomic_sample(out_sample_path, sample)

                record.update(
                    {
                        "status": "ok",
                        "n_pairs": fit.n_pairs,
                        "rmsd": f"{fit.rmsd:.4f}",
                        "n_heme_atoms": len(fit.heme_atoms),
                        "n_fe_atoms": sum(1 for a in fit.heme_atoms if a["atomic_number"] == 26),
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
    summary_path = report_dir / "alphafill_heme_overlay_sample_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if counts["failed"] == 0 and counts["missing_base_sample"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
