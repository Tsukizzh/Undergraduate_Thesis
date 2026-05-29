#!/usr/bin/env python3
"""Audit candidate-id consistency for the full Q10 wet-lab scoring pipeline."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set


EXPECTED_TOTAL = 204
EXPECTED_BY_LIST = {
    "MDP450_no_TM_removed": 43,
    "ARATH_uniprot_p450": 161,
}


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def id_set(rows: Iterable[Dict[str, str]]) -> Set[str]:
    return {str(row.get("candidate_id", "")).strip() for row in rows if str(row.get("candidate_id", "")).strip()}


def count_by_list(rows: Iterable[Dict[str, str]]) -> Dict[str, int]:
    counts: Counter[str] = Counter()
    for row in rows:
        input_list = str(row.get("input_list", "")).strip()
        if input_list:
            counts[input_list] += 1
    return dict(sorted(counts.items()))


def add_set_errors(errors: List[str], name: str, ids: Set[str], expected: Set[str]) -> None:
    missing = sorted(expected - ids)
    extra = sorted(ids - expected)
    if missing:
        errors.append(f"{name} missing candidate_id count={len(missing)} first={missing[:10]}")
    if extra:
        errors.append(f"{name} extra candidate_id count={len(extra)} first={extra[:10]}")


def check_path_column(errors: List[str], rows: List[Dict[str, str]], column: str, table_name: str) -> None:
    missing = []
    for row in rows:
        value = str(row.get(column, "")).strip()
        if not value or not Path(value).exists():
            missing.append(str(row.get("candidate_id", "")))
    if missing:
        errors.append(f"{table_name}.{column} missing path count={len(missing)} first={missing[:10]}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q10-dir", required=True)
    parser.add_argument("--split", default="test")
    args = parser.parse_args()

    q10_dir = Path(args.q10_dir).resolve()
    input_manifest = q10_dir / "structures" / "colabfold_inputs" / "q10_colabfold_input_manifest.csv"
    structure_manifest = q10_dir / "structures" / "manifests" / "q10_colabfold_structure_manifest_msa_m1r1_v1.csv"
    receptor_manifest = q10_dir / "structures" / "manifests" / "q10_heme_fe_receptor_manifest_template_v1.csv"
    box_manifest = q10_dir / "docking" / "manifests" / "q10_unidock_box_manifest_msa_m1r1_v1.csv"
    docking_results = q10_dir / "docking" / "manifests" / "q10_unidock_results_msa_m1r1_v1.csv"
    pt_manifest = q10_dir / "model_cache" / "q10_unidock_pt_msa_m1r1_v1" / "manifests" / f"q10_model_cache_manifest_{args.split}_msa_m1r1_v1.csv"
    score_csv = q10_dir / "results" / "model_scores" / f"q10_model_scores_all_{args.split}_msa_m1r1_v1.csv"

    input_rows = read_csv(input_manifest)
    structure_rows = read_csv(structure_manifest)
    receptor_rows_all = read_csv(receptor_manifest)
    box_rows = read_csv(box_manifest)
    docking_rows_all = read_csv(docking_results)
    pt_rows = read_csv(pt_manifest)
    score_rows = read_csv(score_csv)

    structure_complete = [row for row in structure_rows if row.get("prediction_status") == "complete"]
    receptor_ok = [row for row in receptor_rows_all if row.get("status") == "ok"]
    docking_complete = [row for row in docking_rows_all if str(row.get("docking_status", "")).startswith("complete")]

    expected_ids = id_set(input_rows)
    tables = {
        "input_manifest": input_rows,
        "structure_complete": structure_complete,
        "receptor_ok": receptor_ok,
        "box_manifest": box_rows,
        "docking_complete": docking_complete,
        "pt_manifest": pt_rows,
        "score_csv": score_rows,
    }

    errors: List[str] = []
    if len(expected_ids) != EXPECTED_TOTAL:
        errors.append(f"input candidate_id unique count expected {EXPECTED_TOTAL}, got {len(expected_ids)}")
    for name, rows in tables.items():
        ids = id_set(rows)
        if len(ids) != EXPECTED_TOTAL:
            errors.append(f"{name} candidate_id unique count expected {EXPECTED_TOTAL}, got {len(ids)}")
        add_set_errors(errors, name, ids, expected_ids)
        by_list = count_by_list(rows)
        if by_list != EXPECTED_BY_LIST:
            errors.append(f"{name} input_list counts expected {EXPECTED_BY_LIST}, got {by_list}")

    check_path_column(errors, structure_complete, "selected_pdb_path", "structure_complete")
    check_path_column(errors, receptor_ok, "receptor_pdb", "receptor_ok")
    check_path_column(errors, box_rows, "receptor_pdb_path", "box_manifest")
    check_path_column(errors, docking_complete, "pose_sdf_path", "docking_complete")
    check_path_column(errors, docking_complete, "receptor_pdb_path", "docking_complete")
    check_path_column(errors, pt_rows, "pose_sdf_path", "pt_manifest")
    check_path_column(errors, pt_rows, "receptor_pdb_path", "pt_manifest")
    check_path_column(errors, score_rows, "pose_sdf_path", "score_csv")
    check_path_column(errors, score_rows, "receptor_pdb_path", "score_csv")

    ligand_atom_counts = Counter(str(row.get("n_ligand_atoms", "")).strip() for row in pt_rows)
    if ligand_atom_counts != Counter({"30": EXPECTED_TOTAL}):
        errors.append(f"pt_manifest n_ligand_atoms expected all 30, got {dict(ligand_atom_counts)}")

    receptor_fe_counts = Counter(str(row.get("n_fe_atoms", "")).strip() for row in receptor_ok)
    receptor_quality_counts = Counter(str(row.get("heme_fe_quality", "")).strip() for row in receptor_ok)
    if receptor_fe_counts != Counter({"1": EXPECTED_TOTAL}):
        errors.append(f"receptor_ok n_fe_atoms expected all 1, got {dict(receptor_fe_counts)}")
    if "" in receptor_quality_counts:
        errors.append(f"receptor_ok heme_fe_quality has blank values, got {dict(receptor_quality_counts)}")

    bad_scores = [row.get("candidate_id", "") for row in docking_complete if str(row.get("docking_score", "")).strip() == ""]
    if bad_scores:
        errors.append(f"docking_complete missing docking_score count={len(bad_scores)} first={bad_scores[:10]}")

    ranks = sorted(int(row["rank_overall"]) for row in score_rows if str(row.get("rank_overall", "")).isdigit())
    if ranks != list(range(1, len(score_rows) + 1)):
        errors.append("score_csv rank_overall is not continuous 1..n")

    within_by_list: Dict[str, List[int]] = defaultdict(list)
    for row in score_rows:
        input_list = row.get("input_list", "")
        rank = row.get("rank_within_input_list", "")
        if str(rank).isdigit():
            within_by_list[input_list].append(int(rank))
    for input_list, values in within_by_list.items():
        values_sorted = sorted(values)
        if values_sorted != list(range(1, len(values_sorted) + 1)):
            errors.append(f"score_csv rank_within_input_list not continuous for {input_list}")

    output = {
        "q10_dir": str(q10_dir),
        "expected_total": EXPECTED_TOTAL,
        "expected_by_input_list": EXPECTED_BY_LIST,
        "table_counts": {
            name: {
                "rows": len(rows),
                "unique_candidate_ids": len(id_set(rows)),
                "by_input_list": count_by_list(rows),
            }
            for name, rows in tables.items()
        },
        "docking_status_counts": dict(Counter(row.get("docking_status", "") for row in docking_rows_all)),
        "receptor_fe_counts": dict(receptor_fe_counts),
        "receptor_quality_counts": dict(receptor_quality_counts),
        "receptor_quality_note": "All receptor rows are scored if Fe is present; non-good Fe-Cys geometry is retained as a confidence flag, not silently dropped.",
        "pt_ligand_atom_counts": dict(ligand_atom_counts),
        "errors": errors,
        "ok": not errors,
    }
    out_dir = q10_dir / "results" / "audits"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"q10_full_output_audit_{args.split}_msa_m1r1_v1.json"
    out_path.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(output, ensure_ascii=False, indent=2))

    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
