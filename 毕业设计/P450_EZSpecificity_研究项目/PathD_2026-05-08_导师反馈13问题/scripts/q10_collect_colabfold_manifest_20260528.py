#!/usr/bin/env python3
"""Collect Q10 ColabFold outputs into a structure manifest.

This script is intentionally read-mostly. It scans the Q10 ColabFold output
folders, records which candidate enzymes have usable PDBs, and optionally
creates selected-PDB symlinks under the Q10 workspace.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Optional


RUN_DIR_BY_LIST = {
    "MDP450_no_TM_removed": "CF_MDP450_msa_m1r1_v1",
    "ARATH_uniprot_p450": "CF_ARATH_msa_m1r1_v1",
}


STATUS_PRIORITY = {
    "complete": 5,
    "msa_or_prediction_in_progress": 4,
    "score_missing_pdb": 3,
    "done_missing_pdb": 2,
    "pending_or_not_started": 1,
}


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fasta_wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()[:16]


def first_existing(paths: Iterable[Path]) -> Optional[Path]:
    for path in paths:
        if path.exists():
            return path
    return None


def parse_scores(path: Optional[Path]) -> Dict[str, object]:
    if not path or not path.exists():
        return {
            "mean_plddt": "",
            "min_plddt": "",
            "max_plddt": "",
            "ptm": "",
            "score_parse_error": "",
        }
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        plddt = data.get("plddt") or []
        return {
            "mean_plddt": round(mean(plddt), 3) if plddt else "",
            "min_plddt": round(min(plddt), 3) if plddt else "",
            "max_plddt": round(max(plddt), 3) if plddt else "",
            "ptm": data.get("ptm", ""),
            "score_parse_error": "",
        }
    except Exception as exc:  # noqa: BLE001 - manifest should keep going.
        return {
            "mean_plddt": "",
            "min_plddt": "",
            "max_plddt": "",
            "ptm": "",
            "score_parse_error": repr(exc),
        }


def scan_run_dir(run_dir: Path, fasta_header: str) -> Dict[str, object]:
    done_path = first_existing(run_dir.glob(f"{fasta_header}.done.txt"))
    pdb_path = first_existing(sorted(run_dir.glob(f"{fasta_header}_unrelaxed_rank_*.pdb")))
    score_path = first_existing(sorted(run_dir.glob(f"{fasta_header}_scores_rank_*.json")))
    a3m_path = first_existing(run_dir.glob(f"{fasta_header}.a3m"))

    if pdb_path and score_path:
        status = "complete"
    elif done_path and not pdb_path:
        status = "done_missing_pdb"
    elif score_path and not pdb_path:
        status = "score_missing_pdb"
    elif a3m_path:
        status = "msa_or_prediction_in_progress"
    else:
        status = "pending_or_not_started"

    return {
        "run_dir_name": run_dir.name,
        "status": status,
        "done_path": done_path,
        "pdb_path": pdb_path,
        "score_path": score_path,
        "a3m_path": a3m_path,
    }


def candidate_run_dirs(q10_dir: Path, input_list: str) -> List[Path]:
    batch_dir = q10_dir / "structures" / "colabfold_batch"
    base_name = RUN_DIR_BY_LIST.get(input_list, f"CF_{input_list}_msa_m1r1_v1")
    run_dirs = [batch_dir / base_name]
    run_dirs.extend(sorted(batch_dir.glob("CF_ACCEL_*_msa_m1r1_v1")))
    seen = set()
    unique_dirs = []
    for run_dir in run_dirs:
        key = str(run_dir)
        if key in seen:
            continue
        seen.add(key)
        unique_dirs.append(run_dir)
    return unique_dirs


def create_symlink(target: Path, link: Path) -> str:
    link.parent.mkdir(parents=True, exist_ok=True)
    if link.exists() or link.is_symlink():
        try:
            if link.is_symlink() and link.resolve() == target.resolve():
                return "existing_symlink"
        except OSError:
            pass
        return "exists_not_changed"
    os.symlink(target, link)
    return "created_symlink"


def collect(args: argparse.Namespace) -> None:
    q10_dir = Path(args.q10_dir).resolve()
    input_manifest = q10_dir / "structures" / "colabfold_inputs" / "q10_colabfold_input_manifest.csv"
    sequence_source = q10_dir / "features" / "enzyme_features" / "q10_candidate_enzymes_for_esm.csv"
    out_manifest = q10_dir / "structures" / "manifests" / "q10_colabfold_structure_manifest_msa_m1r1_v1.csv"
    out_summary = q10_dir / "structures" / "manifests" / "q10_colabfold_structure_summary_msa_m1r1_v1.json"
    pending_fasta = q10_dir / "structures" / "manifests" / "q10_colabfold_pending_msa_m1r1_v1.fasta"
    selected_dir = q10_dir / "structures" / "selected_pdb_colabfold_msa_m1r1_v1"

    rows = read_csv(input_manifest)
    seq_by_candidate: Dict[str, str] = {}
    if sequence_source.exists():
        for seq_row in read_csv(sequence_source):
            seq_by_candidate[seq_row.get("candidate_id", "")] = (
                seq_row.get("sequences", "").replace(" ", "").replace("\n", "").strip().upper()
            )
    out_rows: List[Dict[str, object]] = []
    pending_records: List[Dict[str, str]] = []
    by_list: Dict[str, Dict[str, int]] = {}

    for row in rows:
        input_list = row["input_list"]
        fasta_header = row["fasta_header"]
        candidate_id = row["candidate_id"]
        sequence = seq_by_candidate.get(candidate_id, "")

        scans = [scan_run_dir(run_dir, fasta_header) for run_dir in candidate_run_dirs(q10_dir, input_list)]
        best = max(scans, key=lambda item: STATUS_PRIORITY[str(item["status"])])
        run_dir_name = str(best["run_dir_name"])
        done_path = best["done_path"]
        pdb_path = best["pdb_path"]
        score_path = best["score_path"]
        a3m_path = best["a3m_path"]
        status = str(best["status"])

        selected_path = ""
        selected_link_status = ""
        if status == "complete" and pdb_path and args.link_selected:
            link = selected_dir / f"{candidate_id}.pdb"
            selected_link_status = create_symlink(pdb_path, link)
            selected_path = str(link)
        elif status == "complete" and pdb_path:
            selected_path = str(pdb_path)

        scores = parse_scores(score_path)
        seq_len = int(row["sequence_length"])
        out_row: Dict[str, object] = {
            "global_enzyme_id": row["global_enzyme_id"],
            "input_list": input_list,
            "candidate_id": candidate_id,
            "fasta_header": fasta_header,
            "sequence_length": seq_len,
            "sequence_hash": sequence_hash(sequence or fasta_header),
            "colabfold_run_id": run_dir_name,
            "msa_mode": "mmseqs2_uniref_env",
            "num_models": 1,
            "num_recycle": 1,
            "model_type": "alphafold2_ptm",
            "prediction_status": status,
            "done_path": str(done_path) if done_path else "",
            "a3m_path": str(a3m_path) if a3m_path else "",
            "raw_pdb_path": str(pdb_path) if pdb_path else "",
            "score_json_path": str(score_path) if score_path else "",
            "selected_pdb_path": selected_path,
            "selected_link_status": selected_link_status,
            "mean_plddt": scores["mean_plddt"],
            "min_plddt": scores["min_plddt"],
            "max_plddt": scores["max_plddt"],
            "ptm": scores["ptm"],
            "score_parse_error": scores["score_parse_error"],
            "uniprots": row.get("uniprots", ""),
            "description": row.get("description", ""),
        }
        out_rows.append(out_row)

        bucket = by_list.setdefault(input_list, {"total": 0, "complete": 0, "pending": 0, "failed_or_partial": 0})
        bucket["total"] += 1
        if status == "complete":
            bucket["complete"] += 1
        elif status == "pending_or_not_started" or status == "msa_or_prediction_in_progress":
            bucket["pending"] += 1
        else:
            bucket["failed_or_partial"] += 1

        if status != "complete":
            pending_record = dict(row)
            pending_record["sequence"] = sequence
            pending_records.append(pending_record)

    fieldnames = [
        "global_enzyme_id",
        "input_list",
        "candidate_id",
        "fasta_header",
        "sequence_length",
        "sequence_hash",
        "colabfold_run_id",
        "msa_mode",
        "num_models",
        "num_recycle",
        "model_type",
        "prediction_status",
        "done_path",
        "a3m_path",
        "raw_pdb_path",
        "score_json_path",
        "selected_pdb_path",
        "selected_link_status",
        "mean_plddt",
        "min_plddt",
        "max_plddt",
        "ptm",
        "score_parse_error",
        "uniprots",
        "description",
    ]
    write_csv(out_manifest, out_rows, fieldnames)

    with pending_fasta.open("w", encoding="utf-8", newline="\n") as f:
        for row in pending_records:
            seq = row.get("sequence", "")
            if not seq:
                continue
            f.write(f">{row['fasta_header']}\n{fasta_wrap(seq)}\n")

    complete_plddt = [float(r["mean_plddt"]) for r in out_rows if r["prediction_status"] == "complete" and r["mean_plddt"] != ""]
    summary = {
        "q10_dir": str(q10_dir),
        "total": len(out_rows),
        "complete": sum(1 for r in out_rows if r["prediction_status"] == "complete"),
        "pending_or_partial": sum(1 for r in out_rows if r["prediction_status"] != "complete"),
        "by_input_list": by_list,
        "mean_plddt_complete": round(mean(complete_plddt), 3) if complete_plddt else "",
        "manifest_path": str(out_manifest),
        "summary_path": str(out_summary),
        "pending_fasta_path": str(pending_fasta),
        "selected_pdb_dir": str(selected_dir),
    }
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q10-dir", required=True)
    parser.add_argument("--link-selected", action="store_true")
    args = parser.parse_args()
    collect(args)


if __name__ == "__main__":
    main()
