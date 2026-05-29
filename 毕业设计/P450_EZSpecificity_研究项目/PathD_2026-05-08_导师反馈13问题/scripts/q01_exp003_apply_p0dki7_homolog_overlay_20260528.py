#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import torch

from q01_exp002_make_fe_overlay_samples_20260527 import recompute_knn, sample_path, write_atomic_sample
from q01_exp003_make_alphafill_heme_overlay_samples_20260528 import fit_alphafill_heme


BASE = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay")
BASE_CACHE = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1")
OUT_CACHE = BASE / "pt_cache/full_exp003_alphafill_heme_fe_v1"
TRAINABLE = BASE / "manifests/exp003_p450389_trainable_split_sample_manifest_20260527_214331.csv"
OFFICIAL_MANIFEST = BASE / "official_esibank_p450_pdb_20260528_package/exp003_p450_official_esibank_pdb_manifest_20260528.csv"
OFFICIAL_PDB = BASE / "official_esibank_p450_pdb_20260528_package/pdb"
ALPHAFILL_CIF = BASE / "alphafill_p450_216_20260528/cif/I3PLR1.cif"
REPORT_DIR = BASE / "overlay_reports/alphafill_full_v1_p0dki7_homolog"
AUDIT_CSV = REPORT_DIR / "p0dki7_homolog_heme_overlay_audit.csv"
SUMMARY_JSON = REPORT_DIR / "p0dki7_homolog_heme_overlay_summary.json"


def main() -> int:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    trainable_rows = list(csv.DictReader(TRAINABLE.open()))
    official_rows = list(csv.DictReader(OFFICIAL_MANIFEST.open()))
    found = {
        int(row["structure_index"])
        for row in official_rows
        if row["p450_uniprots"] == "P0DKI7"
        and (row["found_in_drive_metadata"].lower() == "true" or row["found_in_drivefs_metadata"].lower() == "true")
    }
    plan = [
        row
        for row in trainable_rows
        if row["uniprot"] == "P0DKI7" and int(row["dock_index"]) in found
    ]

    fieldnames = [
        "split", "sample_id", "dock_index", "uniprot", "status", "template_uniprot",
        "template_pdb", "heme_asym_id", "n_pairs", "rmsd", "n_heme_atoms", "n_fe_atoms",
        "old_n_protein", "new_n_protein", "error",
    ]
    counts = {"ok": 0, "failed": 0, "missing_base_sample": 0}
    rows_out = []
    for row in plan:
        split = row["split"]
        sample_id = int(row["sample_id"])
        dock_index = int(row["dock_index"])
        record = {
            "split": split,
            "sample_id": sample_id,
            "dock_index": dock_index,
            "uniprot": "P0DKI7",
            "status": "",
            "template_uniprot": "I3PLR1",
            "template_pdb": "8E83",
            "heme_asym_id": "E",
            "n_pairs": "",
            "rmsd": "",
            "n_heme_atoms": "",
            "n_fe_atoms": "",
            "old_n_protein": "",
            "new_n_protein": "",
            "error": "",
        }
        base_sample_path = sample_path(BASE_CACHE, split, sample_id)
        out_sample_path = sample_path(OUT_CACHE, split, sample_id)
        if not base_sample_path.exists():
            record["status"] = "missing_base_sample"
            counts["missing_base_sample"] += 1
            rows_out.append(record)
            continue
        try:
            fit = fit_alphafill_heme(
                uniprot="P0DKI7",
                dock_index=dock_index,
                pdb_path=OFFICIAL_PDB / f"{dock_index}.pdb",
                cif_path=ALPHAFILL_CIF,
                heme_asym_id="E",
                min_pairs=80,
                max_rmsd=6.0,
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
            sample["knn_edge_index"] = recompute_knn(sample["ligand_pos"], sample["protein_pos"], 48)
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
        rows_out.append(record)

    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)
    summary = {"planned_samples": len(plan), "counts": counts, "audit_csv": str(AUDIT_CSV)}
    SUMMARY_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if counts["ok"] == len(plan) and counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
