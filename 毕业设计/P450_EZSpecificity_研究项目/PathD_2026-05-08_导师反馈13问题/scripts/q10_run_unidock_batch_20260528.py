#!/usr/bin/env python3
"""Run Uni-Dock for Q10 candidate P450 structures against diosgenin.

The script reads the Q10 box manifest, runs one enzyme receptor against the
single diosgenin ligand, and writes a row-level docking result table. It is
resume-friendly: existing result SDF files with a parsed docking score are
kept unless --rerun-existing is supplied.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional


RESULT_FIELDS = [
    "candidate_id",
    "input_list",
    "fasta_header",
    "receptor_pdb_path",
    "receptor_pdbqt_path",
    "ligand_sdf_path",
    "pose_sdf_path",
    "box_source",
    "heme_fe_quality",
    "nearest_cys_sg_fe_distance",
    "template_pdb_id",
    "fit_rmsd",
    "motif_text",
    "center_x",
    "center_y",
    "center_z",
    "size_x",
    "size_y",
    "size_z",
    "gpu_id",
    "docking_status",
    "docking_score",
    "runtime_sec",
    "stdout_log",
    "stderr_log",
    "failure_reason",
]


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_results(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=RESULT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp_path.replace(path)


def parse_docking_score(sdf_path: Path) -> Optional[float]:
    if not sdf_path.exists():
        return None
    lines = sdf_path.read_text(encoding="utf-8", errors="replace").splitlines()
    for i, line in enumerate(lines):
        if "<docking_score>" in line:
            for value_line in lines[i + 1 : i + 6]:
                value_line = value_line.strip()
                if not value_line:
                    continue
                try:
                    return float(value_line)
                except ValueError:
                    continue
        if "<Uni-Dock RESULT>" in line:
            for value_line in lines[i + 1 : i + 8]:
                match = re.search(r"ENERGY=\s*([-+]?\d+(?:\.\d+)?)", value_line)
                if match:
                    return float(match.group(1))
    return None


def load_existing_results(path: Path) -> Dict[str, Dict[str, object]]:
    if not path.exists():
        return {}
    rows = read_csv(path)
    return {row["candidate_id"]: row for row in rows if row.get("candidate_id")}


def strip_model_tags(input_pdb: Path, output_pdb: Path) -> None:
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with input_pdb.open(encoding="utf-8", errors="replace") as fin, output_pdb.open("w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(("MODEL", "ENDMDL")):
                continue
            fout.write(line)


def extract_heme_lines(input_pdb: Path) -> List[str]:
    lines: List[str] = []
    with input_pdb.open(encoding="utf-8", errors="replace") as fin:
        for line in fin:
            if line.startswith("HETATM") and line[17:20].strip() in {"HEM", "HEC"}:
                lines.append(line.rstrip("\n"))
    if not lines:
        raise RuntimeError(f"no HEM/HEC lines found in receptor PDB: {input_pdb}")
    return lines


def heme_atom_type(atom_name: str, element: str) -> str:
    raw = (element or "").strip()
    if raw.upper() == "FE" or atom_name.strip().upper() == "FE":
        return "Fe"
    elem = raw[:1].upper()
    if not elem:
        elem = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper()
    return {
        "C": "C",
        "N": "NA",
        "O": "OA",
        "S": "SA",
        "H": "HD",
        "P": "P",
    }.get(elem, "C")


def append_heme_to_pdbqt(protein_pdbqt: Path, heme_lines: List[str], receptor_pdbqt: Path) -> Dict[str, int]:
    protein_lines = protein_pdbqt.read_text(encoding="utf-8", errors="replace").splitlines()
    max_serial = 0
    for line in protein_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                pass

    out_lines = list(protein_lines)
    serial = max_serial + 1
    added = 0
    fe_count = 0
    for line in heme_lines:
        atom_name = line[12:16].strip()
        resname = line[17:20].strip() or "HEM"
        chain_id = (line[21].strip() or "H")[:1]
        try:
            residue_id = int(line[22:26])
        except ValueError:
            residue_id = 1
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        atom_type = heme_atom_type(atom_name, line[76:78].strip() if len(line) >= 78 else "")
        if atom_type == "Fe":
            fe_count += 1
        out_lines.append(
            f"ATOM  {serial:5d} {atom_name:^4s} {resname:3s} {chain_id:1s}{residue_id:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}    {0.0:6.3f} {atom_type:<2s}"
        )
        serial += 1
        added += 1
    out_lines.append("TER")
    receptor_pdbqt.parent.mkdir(parents=True, exist_ok=True)
    receptor_pdbqt.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return {"heme_atoms_added": added, "fe_atoms_added": fe_count}


def prepare_receptor_pdbqt(
    row: Dict[str, str],
    unidock_env: Path,
    amber_env: Path,
    rerun_existing: bool,
) -> Dict[str, object]:
    receptor_pdb = Path(row["receptor_pdb_path"])
    receptor_pdbqt = Path(row["receptor_pdbqt_path"])
    proteinprep_workdir = Path(row["proteinprep_workdir"])
    prepared_pdb = proteinprep_workdir / "receptor_single_model.pdb"
    protein_only_pdbqt = proteinprep_workdir / "protein_only.pdbqt"

    if receptor_pdbqt.exists() and not rerun_existing:
        text = receptor_pdbqt.read_text(encoding="utf-8", errors="replace")
        return {
            "status": "existing",
            "heme_atoms_added": sum(1 for line in text.splitlines() if line[17:20].strip() in {"HEM", "HEC"}),
            "fe_atoms_added": sum(1 for line in text.splitlines() if line[17:20].strip() in {"HEM", "HEC"} and line.rstrip().endswith("Fe")),
            "error": "",
        }

    proteinprep_workdir.mkdir(parents=True, exist_ok=True)
    strip_model_tags(receptor_pdb, prepared_pdb)
    heme_lines = extract_heme_lines(receptor_pdb)
    cmd = [
        "unidocktools",
        "proteinprep",
        "-r",
        str(prepared_pdb),
        "-ph",
        "-pr",
        "-wd",
        str(proteinprep_workdir),
        "-o",
        str(protein_only_pdbqt),
    ]
    env = os.environ.copy()
    env["PATH"] = f"{amber_env / 'bin'}:{unidock_env / 'bin'}:" + env.get("PATH", "")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=900, env=env)
    (proteinprep_workdir / "proteinprep.stdout.log").write_text(proc.stdout, encoding="utf-8")
    (proteinprep_workdir / "proteinprep.stderr.log").write_text(proc.stderr, encoding="utf-8")
    if proc.returncode != 0 or not protein_only_pdbqt.exists():
        return {
            "status": "failed",
            "heme_atoms_added": 0,
            "fe_atoms_added": 0,
            "error": f"proteinprep returncode={proc.returncode}; stderr_tail={(proc.stderr or '')[-500:]}",
        }

    stats = append_heme_to_pdbqt(protein_only_pdbqt, heme_lines, receptor_pdbqt)
    if stats["fe_atoms_added"] < 1:
        return {
            "status": "failed",
            **stats,
            "error": "receptor PDBQT has no appended Fe atom",
        }
    return {"status": "complete", **stats, "error": ""}


def build_command(row: Dict[str, str], ligand_sdf: Path, workdir: Path, savedir: Path, seed: int) -> List[str]:
    return [
        "unidock",
        "--receptor",
        row["receptor_pdbqt_path"],
        "--gpu_batch",
        str(ligand_sdf),
        "--center_x",
        row["center_x"],
        "--center_y",
        row["center_y"],
        "--center_z",
        row["center_z"],
        "--size_x",
        row["size_x"],
        "--size_y",
        row["size_y"],
        "--size_z",
        row["size_z"],
        "--dir",
        str(savedir),
        "--scoring",
        "vina",
        "--search_mode",
        "fast",
        "--num_modes",
        "1",
        "--seed",
        str(seed),
    ]


def run_one(
    row: Dict[str, str],
    q10_dir: Path,
    ligand_sdf: Path,
    unidock_env: Path,
    amber_env: Path,
    gpu_id: str,
    seed: int,
    rerun_existing: bool,
) -> Dict[str, object]:
    candidate_id = row["candidate_id"]
    savedir = q10_dir / "docking" / "results_unidock_msa_m1r1_v1" / candidate_id
    workdir = q10_dir / "docking" / "work_unidock_msa_m1r1_v1" / candidate_id
    logdir = q10_dir / "docking" / "logs_unidock_msa_m1r1_v1"
    logdir.mkdir(parents=True, exist_ok=True)
    stdout_log = logdir / f"{candidate_id}.stdout.log"
    stderr_log = logdir / f"{candidate_id}.stderr.log"
    pose_sdf = savedir / f"{ligand_sdf.stem}_out.sdf"

    if pose_sdf.exists() and not rerun_existing:
        score = parse_docking_score(pose_sdf)
        if score is not None:
            return result_row(row, ligand_sdf, pose_sdf, gpu_id, "complete_existing", score, 0.0, stdout_log, stderr_log, "")

    savedir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)
    prep = prepare_receptor_pdbqt(row, unidock_env, amber_env, rerun_existing)
    if prep["status"] == "failed":
        return result_row(
            row,
            ligand_sdf,
            pose_sdf,
            gpu_id,
            "failed_receptor_pdbqt",
            "",
            0.0,
            stdout_log,
            stderr_log,
            str(prep["error"]),
        )

    cmd = build_command(row, ligand_sdf, workdir, savedir, seed)
    env = os.environ.copy()
    env["PATH"] = f"{amber_env / 'bin'}:{unidock_env / 'bin'}:" + env.get("PATH", "")
    env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

    t0 = time.time()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=1800, env=env)
        runtime = time.time() - t0
        stdout_log.write_text(proc.stdout, encoding="utf-8")
        stderr_log.write_text(proc.stderr, encoding="utf-8")
        score = parse_docking_score(pose_sdf)
        if proc.returncode == 0 and score is not None:
            return result_row(row, ligand_sdf, pose_sdf, gpu_id, "complete", score, runtime, stdout_log, stderr_log, "")
        reason = f"returncode={proc.returncode}; score={score}; stderr_tail={(proc.stderr or '')[-300:]}"
        return result_row(row, ligand_sdf, pose_sdf, gpu_id, "failed", score if score is not None else "", runtime, stdout_log, stderr_log, reason)
    except Exception as exc:  # noqa: BLE001 - batch should keep going.
        runtime = time.time() - t0
        return result_row(row, ligand_sdf, pose_sdf, gpu_id, "failed_exception", "", runtime, stdout_log, stderr_log, repr(exc))


def result_row(
    row: Dict[str, str],
    ligand_sdf: Path,
    pose_sdf: Path,
    gpu_id: str,
    status: str,
    score: object,
    runtime: float,
    stdout_log: Path,
    stderr_log: Path,
    failure_reason: str,
) -> Dict[str, object]:
    return {
        "candidate_id": row["candidate_id"],
        "input_list": row["input_list"],
        "fasta_header": row["fasta_header"],
        "receptor_pdb_path": row["receptor_pdb_path"],
        "receptor_pdbqt_path": row["receptor_pdbqt_path"],
        "ligand_sdf_path": str(ligand_sdf),
        "pose_sdf_path": str(pose_sdf),
        "box_source": row["box_source"],
        "heme_fe_quality": row.get("heme_fe_quality", ""),
        "nearest_cys_sg_fe_distance": row.get("nearest_cys_sg_fe_distance", ""),
        "template_pdb_id": row.get("template_pdb_id", ""),
        "fit_rmsd": row.get("fit_rmsd", ""),
        "motif_text": row["motif_text"],
        "center_x": row["center_x"],
        "center_y": row["center_y"],
        "center_z": row["center_z"],
        "size_x": row["size_x"],
        "size_y": row["size_y"],
        "size_z": row["size_z"],
        "gpu_id": gpu_id,
        "docking_status": status,
        "docking_score": score,
        "runtime_sec": round(runtime, 3),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "failure_reason": failure_reason,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q10-dir", required=True)
    parser.add_argument("--unidock-env", default="/root/autodl-tmp/envs/unidock_env")
    parser.add_argument("--amber-env", default="/root/autodl-tmp/envs/ambertools_q10")
    parser.add_argument("--gpu-ids", default="0")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--input-list", default="")
    parser.add_argument("--rerun-existing", action="store_true")
    parser.add_argument(
        "--fresh-manifest",
        action="store_true",
        help="Start the result manifest from the current requested rows only.",
    )
    parser.add_argument("--seed", type=int, default=20260528)
    args = parser.parse_args()

    q10_dir = Path(args.q10_dir).resolve()
    unidock_env = Path(args.unidock_env).resolve()
    amber_env = Path(args.amber_env).resolve()
    ligand_sdf = q10_dir / "ligands" / "pubchem_3d" / "diosgenin_CID99474_pubchem_3d.sdf"
    box_manifest = q10_dir / "docking" / "manifests" / "q10_unidock_box_manifest_msa_m1r1_v1.csv"
    result_manifest = q10_dir / "docking" / "manifests" / "q10_unidock_results_msa_m1r1_v1.csv"
    summary_path = q10_dir / "docking" / "manifests" / "q10_unidock_results_summary_msa_m1r1_v1.json"

    rows = [row for row in read_csv(box_manifest) if Path(row["receptor_pdb_path"]).exists()]
    if args.input_list:
        rows = [row for row in rows if row["input_list"] == args.input_list]
    if args.limit:
        rows = rows[: args.limit]

    existing = {} if args.fresh_manifest else load_existing_results(result_manifest)
    results: List[Dict[str, object]] = list(existing.values())
    result_by_id: Dict[str, Dict[str, object]] = dict(existing)
    gpu_ids = [gpu.strip() for gpu in args.gpu_ids.split(",") if gpu.strip()]
    if not gpu_ids:
        gpu_ids = ["0"]

    lock = threading.Lock()
    pending = [row for row in rows if args.rerun_existing or row["candidate_id"] not in result_by_id or not str(result_by_id[row["candidate_id"]].get("docking_status", "")).startswith("complete")]

    def wrapped(index_row: tuple[int, Dict[str, str]]) -> Dict[str, object]:
        idx, row = index_row
        gpu_id = gpu_ids[idx % len(gpu_ids)]
        return run_one(row, q10_dir, ligand_sdf, unidock_env, amber_env, gpu_id, args.seed, args.rerun_existing)

    with ThreadPoolExecutor(max_workers=len(gpu_ids)) as executor:
        futures = [executor.submit(wrapped, item) for item in enumerate(pending)]
        for future in as_completed(futures):
            row = future.result()
            with lock:
                result_by_id[str(row["candidate_id"])] = row
                results = [result_by_id[k] for k in sorted(result_by_id)]
                write_results(result_manifest, results)

    final_rows = [result_by_id[k] for k in sorted(result_by_id)]
    write_results(result_manifest, final_rows)
    by_status: Dict[str, int] = {}
    by_input_list: Dict[str, int] = {}
    for row in final_rows:
        by_status[str(row["docking_status"])] = by_status.get(str(row["docking_status"]), 0) + 1
        by_input_list[str(row["input_list"])] = by_input_list.get(str(row["input_list"]), 0) + 1
    scores = [float(row["docking_score"]) for row in final_rows if str(row.get("docking_score", "")) not in {"", "None"}]
    summary = {
        "box_manifest": str(box_manifest),
        "result_manifest": str(result_manifest),
        "total_result_rows": len(final_rows),
        "requested_rows_this_run": len(rows),
        "pending_run_count": len(pending),
        "fresh_manifest": args.fresh_manifest,
        "by_status": by_status,
        "by_input_list": by_input_list,
        "best_score": min(scores) if scores else "",
        "worst_score": max(scores) if scores else "",
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
