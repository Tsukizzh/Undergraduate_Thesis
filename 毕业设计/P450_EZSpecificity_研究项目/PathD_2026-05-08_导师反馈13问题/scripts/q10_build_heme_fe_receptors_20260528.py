#!/usr/bin/env python3
"""Build Q10 HEM/Fe-aware receptor PDBs by template overlay.

This script is intentionally additive: it reads ColabFold apo structures and
audited AlphaFill P450 HEM templates, then writes new receptor PDB files under
the Q10 directory. It does not edit the original ColabFold outputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.cealign import CEAligner


PROJECT_SCRIPTS = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/scripts")
if PROJECT_SCRIPTS.exists():
    sys.path.insert(0, str(PROJECT_SCRIPTS))

from q01_exp003_make_alphafill_heme_overlay_samples_20260528 import (  # noqa: E402
    fit_alphafill_heme,
    read_alphafill_ca_and_heme,
)


DEFAULT_TEMPLATE_PDB_IDS = "6VBY,5YLW,7X2Q,6A15,8E83,3RUK,2VE3"


@dataclass
class Template:
    uniprot: str
    pdb_id: str
    heme_asym_id: str
    identity: float
    length: int
    local_rmsd: float
    clash_score: float
    cif_path: Path


@dataclass
class FitResult:
    n_pairs: int
    rmsd: float
    score: float
    heme_atoms: list[dict]


def parse_boolish(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "ok", "complete"}


def load_complete_targets(manifest_path: Path) -> list[dict]:
    rows = []
    with manifest_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("prediction_status") != "complete":
                continue
            raw_pdb = row.get("selected_pdb_path") or row.get("raw_pdb_path")
            if not raw_pdb:
                continue
            p = Path(raw_pdb)
            if not p.exists():
                continue
            rows.append(row)
    return rows


def numeric(value: str, default: float = 0.0) -> float:
    try:
        if value is None or value == "":
            return default
        return float(value)
    except Exception:
        return default


def load_templates(summary_csv: Path, cif_dir: Path, pdb_ids: set[str], max_templates_per_pdb: int) -> list[Template]:
    grouped: dict[str, list[Template]] = {}
    with summary_csv.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if numeric(row.get("hem_candidates"), 0) <= 0:
                continue
            pdb_id = str(row.get("best_pdb_id", "")).upper()
            if pdb_ids and pdb_id not in pdb_ids:
                continue
            uniprot = str(row["uniprot"])
            cif_path = cif_dir / f"{uniprot}.cif"
            if not cif_path.exists():
                continue
            item = Template(
                uniprot=uniprot,
                pdb_id=pdb_id,
                heme_asym_id=str(row["best_asym_id"]),
                identity=numeric(row.get("best_identity"), 0),
                length=int(numeric(row.get("best_length"), 0)),
                local_rmsd=numeric(row.get("best_local_rmsd"), 999),
                clash_score=numeric(row.get("best_clash_score"), 999),
                cif_path=cif_path,
            )
            grouped.setdefault(pdb_id, []).append(item)

    selected = []
    for pdb_id, items in sorted(grouped.items()):
        items.sort(key=lambda x: (x.identity, x.length, -x.local_rmsd, -x.clash_score), reverse=True)
        selected.extend(items[:max_templates_per_pdb])
    selected.sort(key=lambda x: (x.identity, x.length, -x.local_rmsd, -x.clash_score), reverse=True)
    return selected


def pdb_atom_serial_max(path: Path) -> int:
    max_serial = 0
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    max_serial = max(max_serial, int(line[6:11]))
                except Exception:
                    pass
    return max_serial


def read_cys_sg(path: Path) -> list[tuple[str, int, np.ndarray]]:
    atoms = []
    with path.open("r", errors="ignore") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[17:20].strip() != "CYS":
                continue
            if line[12:16].strip() != "SG":
                continue
            chain = line[21:22].strip() or "?"
            try:
                resseq = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
            except Exception:
                continue
            atoms.append((chain, resseq, xyz))
    return atoms


def nearest_cys_to_fe(path: Path, fe_xyz: np.ndarray | None) -> tuple[str, int | str, float | str]:
    if fe_xyz is None:
        return "", "", ""
    best = None
    for chain, resseq, xyz in read_cys_sg(path):
        dist = float(np.linalg.norm(xyz - fe_xyz))
        if best is None or dist < best[2]:
            best = (chain, resseq, dist)
    if best is None:
        return "", "", ""
    return best[0], best[1], round(best[2], 3)


def nearest_cys_to_xyz(path: Path, xyz: np.ndarray | None) -> tuple[str, int | str, float | str]:
    return nearest_cys_to_fe(path, xyz)


def first_fe_xyz(heme_atoms: list[dict]) -> np.ndarray | None:
    for atom in heme_atoms:
        if str(atom.get("element", "")).upper() == "FE":
            return np.asarray(atom["xyz_transformed"], dtype=float)
    return None


def fe_cys_quality(distance: float | str, ideal: float) -> tuple[int, float]:
    if distance == "" or distance is None:
        return 0, 999.0
    d = float(distance)
    delta = abs(d - ideal)
    if 1.8 <= d <= 3.5:
        return 3, delta
    if 1.2 <= d <= 5.0:
        return 2, delta
    if 0.8 <= d <= 6.0:
        return 1, delta
    return 0, delta


def pdb_hetatm_line(serial: int, atom_name: str, resname: str, chain: str, resseq: int, xyz: np.ndarray, element: str) -> str:
    atom_name = str(atom_name).strip()[:4] or element
    resname = str(resname).strip()[:3] or "HEM"
    element = str(element).strip().upper()[:2] or atom_name[0].upper()
    return (
        f"HETATM{serial:5d} {atom_name:<4} {resname:>3} {chain:1}{resseq:4d}    "
        f"{float(xyz[0]):8.3f}{float(xyz[1]):8.3f}{float(xyz[2]):8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          {element:>2}\n"
    )


def write_receptor_with_heme(target_pdb: Path, out_pdb: Path, heme_atoms: list[dict]) -> tuple[int, int, np.ndarray | None]:
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    source_lines = []
    with target_pdb.open("r", errors="ignore") as handle:
        for line in handle:
            if line.startswith("END"):
                continue
            source_lines.append(line if line.endswith("\n") else line + "\n")

    serial = pdb_atom_serial_max(target_pdb)
    n_fe = 0
    fe_xyz = None
    het_lines = []
    for atom in heme_atoms:
        elem = str(atom["element"]).upper()
        comp_id = str(atom.get("comp_id") or "HEM").upper()
        resname = "HEM" if comp_id in {"HEM", "HEC", "FE"} else comp_id[:3]
        xyz = np.asarray(atom["xyz_transformed"], dtype=float)
        serial += 1
        het_lines.append(pdb_hetatm_line(serial, atom.get("atom_name", elem), resname, "H", 1, xyz, elem))
        if elem == "FE":
            n_fe += 1
            if fe_xyz is None:
                fe_xyz = xyz

    out_pdb.write_text("".join(source_lines + het_lines + ["END\n"]), encoding="utf-8")
    return len(het_lines), n_fe, fe_xyz


def fit_template_heme_ce(target_pdb: Path, template: Template, window_size: int, max_gap: int) -> FitResult:
    """Align AlphaFill template structure onto target apo model with CEAligner."""
    pdb_parser = PDBParser(QUIET=True)
    cif_parser = MMCIFParser(QUIET=True)
    target_structure = pdb_parser.get_structure("target", str(target_pdb))
    mobile_structure = cif_parser.get_structure(f"template_{template.uniprot}", str(template.cif_path))
    _, heme_atoms = read_alphafill_ca_and_heme(template.cif_path, template.heme_asym_id)
    if not heme_atoms:
        raise RuntimeError(f"template has no selected HEM/Fe atoms for asym_id={template.heme_asym_id}")

    aligner = CEAligner(window_size=window_size, max_gap=max_gap)
    aligner.set_reference(target_structure)
    aligner.align(mobile_structure, transform=False)
    if not math.isfinite(float(aligner.rms)):
        raise RuntimeError("CE alignment produced non-finite RMSD")

    rot, tran = aligner._rigid_motion
    transformed = []
    for atom in heme_atoms:
        xyz = np.asarray(atom["xyz"], dtype=float) @ rot + tran
        transformed.append({**atom, "xyz_transformed": xyz.astype(np.float32)})
    if sum(1 for atom in transformed if str(atom["element"]).upper() == "FE") < 1:
        raise RuntimeError("selected HEM transplant has no Fe atom")

    n_pairs = min(len(aligner.refcoord), len(aligner._coord))
    return FitResult(n_pairs=n_pairs, rmsd=float(aligner.rms), score=float(n_pairs) - float(aligner.rms), heme_atoms=transformed)


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--q10-dir", required=True)
    parser.add_argument("--structure-manifest", default=None)
    parser.add_argument("--alphafill-summary", required=True)
    parser.add_argument("--alphafill-cif-dir", required=True)
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--template-pdb-ids", default=DEFAULT_TEMPLATE_PDB_IDS)
    parser.add_argument("--max-templates-per-pdb", type=int, default=2)
    parser.add_argument("--min-pairs", type=int, default=120)
    parser.add_argument("--max-rmsd", type=float, default=8.0)
    parser.add_argument("--alignment-mode", choices=["ce", "sequence"], default="ce")
    parser.add_argument("--ce-window-size", type=int, default=8)
    parser.add_argument("--ce-max-gap", type=int, default=30)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    q10_dir = Path(args.q10_dir)
    manifest_path = Path(args.structure_manifest) if args.structure_manifest else q10_dir / "structures/manifests/q10_colabfold_structure_manifest_msa_m1r1_v1.csv"
    out_dir = Path(args.out_dir) if args.out_dir else q10_dir / "structures/receptors_heme_fe_template_v1"
    receptor_dir = out_dir / "pdb"
    manifest_dir = q10_dir / "structures/manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    audit_path = manifest_dir / "q10_heme_fe_receptor_manifest_template_v1.csv"
    summary_path = manifest_dir / "q10_heme_fe_receptor_summary_template_v1.json"

    pdb_ids = {x.strip().upper() for x in args.template_pdb_ids.split(",") if x.strip()}
    templates = load_templates(Path(args.alphafill_summary), Path(args.alphafill_cif_dir), pdb_ids, args.max_templates_per_pdb)
    if not templates:
        raise SystemExit("no HEM templates loaded")

    targets = load_complete_targets(manifest_path)
    if args.limit:
        targets = targets[: args.limit]

    fields = [
        "global_enzyme_id", "input_list", "candidate_id", "status", "target_pdb", "receptor_pdb",
        "template_uniprot", "template_pdb_id", "template_heme_asym_id", "fit_pairs", "fit_rmsd",
        "fit_score", "n_heme_atoms", "n_fe_atoms", "fe_x", "fe_y", "fe_z",
        "nearest_cys_chain", "nearest_cys_resseq", "nearest_cys_sg_fe_distance",
        "heme_fe_quality", "mean_plddt", "ptm", "error",
    ]

    existing_ok = set()
    if audit_path.exists() and not args.overwrite:
        with audit_path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("status") == "ok":
                    existing_ok.add(row.get("candidate_id"))

    write_header = args.overwrite or not audit_path.exists()
    mode = "w" if args.overwrite else "a"
    counts = {"ok": 0, "skip_done": 0, "failed": 0, "target_complete": len(targets)}
    start = time.time()
    with audit_path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if write_header:
            writer.writeheader()

        for row in targets:
            candidate_id = row["candidate_id"]
            out_pdb = receptor_dir / f"{safe_name(candidate_id)}_heme_fe_receptor.pdb"
            if candidate_id in existing_ok and out_pdb.exists():
                counts["skip_done"] += 1
                continue

            target_pdb = Path(row.get("selected_pdb_path") or row.get("raw_pdb_path"))
            record = {
                "global_enzyme_id": row.get("global_enzyme_id", ""),
                "input_list": row.get("input_list", ""),
                "candidate_id": candidate_id,
                "status": "",
                "target_pdb": str(target_pdb),
                "receptor_pdb": str(out_pdb),
                "template_uniprot": "",
                "template_pdb_id": "",
                "template_heme_asym_id": "",
                "fit_pairs": "",
                "fit_rmsd": "",
                "fit_score": "",
                "n_heme_atoms": "",
                "n_fe_atoms": "",
                "fe_x": "",
                "fe_y": "",
                "fe_z": "",
                "nearest_cys_chain": "",
                "nearest_cys_resseq": "",
                "nearest_cys_sg_fe_distance": "",
                "heme_fe_quality": "",
                "mean_plddt": row.get("mean_plddt", ""),
                "ptm": row.get("ptm", ""),
                "error": "",
            }
            best = None
            errors = []
            for template in templates:
                try:
                    if args.alignment_mode == "ce":
                        fit = fit_template_heme_ce(target_pdb, template, args.ce_window_size, args.ce_max_gap)
                    else:
                        fit = fit_alphafill_heme(
                            uniprot=template.uniprot,
                            dock_index=int(row.get("global_enzyme_id") or 0),
                            pdb_path=target_pdb,
                            cif_path=template.cif_path,
                            heme_asym_id=template.heme_asym_id,
                            min_pairs=args.min_pairs,
                            max_rmsd=args.max_rmsd,
                        )
                    if fit.rmsd > args.max_rmsd:
                        raise RuntimeError(f"alignment RMSD too high: {fit.rmsd:.3f}")
                    _, _, cys_dist = nearest_cys_to_xyz(target_pdb, first_fe_xyz(fit.heme_atoms))
                    qlevel, qdelta = fe_cys_quality(cys_dist, 2.3)
                    candidate = (qlevel, -qdelta, -fit.rmsd, fit.n_pairs, fit.score, fit, template)
                    if best is None or candidate[:3] > best[:3]:
                        best = candidate
                except Exception as exc:
                    errors.append(f"{template.pdb_id}/{template.uniprot}:{exc}")

            if best is None:
                record["status"] = "failed"
                record["error"] = " | ".join(errors[:5])
                counts["failed"] += 1
                writer.writerow(record)
                continue

            _, _, _, _, _, fit, template = best
            n_heme, n_fe, fe_xyz = write_receptor_with_heme(target_pdb, out_pdb, fit.heme_atoms)
            cys_chain, cys_resseq, cys_dist = nearest_cys_to_fe(out_pdb, fe_xyz)
            qlevel, _ = fe_cys_quality(cys_dist, 2.3)
            quality = {3: "good_fe_cys", 2: "warning_fe_cys", 1: "poor_fe_cys", 0: "bad_or_missing_fe_cys"}[qlevel]
            record.update(
                {
                    "status": "ok",
                    "template_uniprot": template.uniprot,
                    "template_pdb_id": template.pdb_id,
                    "template_heme_asym_id": template.heme_asym_id,
                    "fit_pairs": fit.n_pairs,
                    "fit_rmsd": f"{fit.rmsd:.4f}",
                    "fit_score": f"{fit.score:.3f}",
                    "n_heme_atoms": n_heme,
                    "n_fe_atoms": n_fe,
                    "fe_x": "" if fe_xyz is None else f"{fe_xyz[0]:.3f}",
                    "fe_y": "" if fe_xyz is None else f"{fe_xyz[1]:.3f}",
                    "fe_z": "" if fe_xyz is None else f"{fe_xyz[2]:.3f}",
                    "nearest_cys_chain": cys_chain,
                    "nearest_cys_resseq": cys_resseq,
                    "nearest_cys_sg_fe_distance": cys_dist,
                    "heme_fe_quality": quality,
                }
            )
            counts["ok"] += 1
            writer.writerow(record)

    summary = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "q10_dir": str(q10_dir),
        "structure_manifest": str(manifest_path),
        "out_dir": str(out_dir),
        "audit_csv": str(audit_path),
        "templates_loaded": len(templates),
        "template_pdb_ids": sorted({t.pdb_id for t in templates}),
        "counts": counts,
        "elapsed_sec": round(time.time() - start, 3),
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if counts["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
