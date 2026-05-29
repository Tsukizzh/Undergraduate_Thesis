#!/usr/bin/env python3
"""Audit EXP003 P450 official ESIBank PDB downloads and HEM/Fe presence."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path


PROJECT = Path(r"D:\EZSpecificity_Project")
PATHD = PROJECT / r"毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题"

MANIFEST = PATHD / r"data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\official_pdb_manifest_20260528\exp003_p450_official_esibank_pdb_manifest_20260528.csv"
PDB_DIR = Path(r"D:\ESIBank\EXP003_p450_official_esibank_pdb_20260528\pdb")
OUT_DIR = PATHD / r"data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\official_pdb_manifest_20260528"
OUT_CSV = OUT_DIR / "exp003_p450_official_pdb_download_heme_fe_audit_20260528.csv"
OUT_JSON = OUT_DIR / "exp003_p450_official_pdb_download_heme_fe_summary_20260528.json"
LOCAL_JSON = Path(r"D:\ESIBank\EXP003_p450_official_esibank_pdb_20260528\local_audit_heme_fe_20260528.json")

HEME_NAMES = {
    "HEM",
    "HEA",
    "HEB",
    "HEC",
    "HDD",
    "HDM",
    "HEO",
    "HEV",
    "HAS",
}


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def as_bool(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def as_int(value: object, default: int = 0) -> int:
    try:
        text = str(value).strip()
        if not text:
            return default
        return int(float(text))
    except Exception:
        return default


def looks_like_pdb(path: Path) -> bool:
    try:
        head = path.read_bytes()[:8192]
    except OSError:
        return False
    return head.startswith(b"HEADER") or b"ATOM" in head or b"HETATM" in head or head.startswith(b"REMARK")


def scan_pdb(path: Path) -> dict[str, object]:
    result = {
        "atom_lines": 0,
        "hetatm_lines": 0,
        "heme_line_count": 0,
        "fe_line_count": 0,
        "heme_resnames": set(),
        "fe_resnames": set(),
    }
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("ATOM"):
                result["atom_lines"] += 1
            elif line.startswith("HETATM"):
                result["hetatm_lines"] += 1
            else:
                continue
            resname = line[17:20].strip().upper()
            atom_name = line[12:16].strip().upper()
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if resname in HEME_NAMES:
                result["heme_line_count"] += 1
                result["heme_resnames"].add(resname)
            if element == "FE" or atom_name == "FE":
                result["fe_line_count"] += 1
                result["fe_resnames"].add(resname or "UNKNOWN")
    result["heme_resnames"] = ";".join(sorted(result["heme_resnames"]))
    result["fe_resnames"] = ";".join(sorted(result["fe_resnames"]))
    result["has_heme"] = result["heme_line_count"] > 0
    result["has_fe"] = result["fe_line_count"] > 0
    result["has_heme_and_fe"] = bool(result["has_heme"] and result["has_fe"])
    return result


def main() -> int:
    rows = read_csv(MANIFEST)
    audit_rows: list[dict[str, object]] = []

    for row in rows:
        found = as_bool(row.get("found_in_drive_metadata")) or as_bool(row.get("found_in_drivefs_metadata"))
        expected_size = as_int(row.get("file_size"))
        pdb_name = str(row.get("pdb_name", "")).strip()
        path = PDB_DIR / pdb_name
        exists = path.exists()
        actual_size = path.stat().st_size if exists else 0
        size_match = bool(exists and expected_size > 0 and actual_size == expected_size)
        pdb_like = bool(exists and looks_like_pdb(path))
        scan = scan_pdb(path) if pdb_like else {
            "atom_lines": 0,
            "hetatm_lines": 0,
            "heme_line_count": 0,
            "fe_line_count": 0,
            "heme_resnames": "",
            "fe_resnames": "",
            "has_heme": False,
            "has_fe": False,
            "has_heme_and_fe": False,
        }
        audit_rows.append(
            {
                "structure_index": row.get("structure_index", ""),
                "pdb_name": pdb_name,
                "drive_file_id": row.get("drive_file_id", ""),
                "found_in_drive_metadata": found,
                "is_trainable_in_exp003_manifest": as_bool(row.get("is_trainable_in_exp003_manifest")),
                "trainable_rows": row.get("trainable_rows", "0"),
                "expected_size": expected_size,
                "exists": exists,
                "actual_size": actual_size,
                "size_match": size_match,
                "looks_like_pdb": pdb_like,
                **scan,
            }
        )

    fieldnames = [
        "structure_index",
        "pdb_name",
        "drive_file_id",
        "found_in_drive_metadata",
        "is_trainable_in_exp003_manifest",
        "trainable_rows",
        "expected_size",
        "exists",
        "actual_size",
        "size_match",
        "looks_like_pdb",
        "atom_lines",
        "hetatm_lines",
        "has_heme",
        "has_fe",
        "has_heme_and_fe",
        "heme_line_count",
        "fe_line_count",
        "heme_resnames",
        "fe_resnames",
    ]
    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in audit_rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})

    found_rows = [row for row in audit_rows if row["found_in_drive_metadata"]]
    trainable_rows = [row for row in audit_rows if row["is_trainable_in_exp003_manifest"]]
    found_trainable_rows = [
        row for row in trainable_rows
        if row["found_in_drive_metadata"]
    ]
    summary = {
        "created_at": now(),
        "manifest": str(MANIFEST),
        "pdb_dir": str(PDB_DIR),
        "audit_csv": str(OUT_CSV),
        "counts": {
            "manifest_rows": len(audit_rows),
            "drive_metadata_found_rows": len(found_rows),
            "download_exists_rows": sum(1 for row in found_rows if row["exists"]),
            "download_missing_rows": sum(1 for row in found_rows if not row["exists"]),
            "download_size_match_rows": sum(1 for row in found_rows if row["size_match"]),
            "download_bad_size_rows": sum(1 for row in found_rows if row["exists"] and not row["size_match"]),
            "download_pdb_like_rows": sum(1 for row in found_rows if row["looks_like_pdb"]),
            "all_p450_has_heme_rows": sum(1 for row in found_rows if row["has_heme"]),
            "all_p450_has_fe_rows": sum(1 for row in found_rows if row["has_fe"]),
            "all_p450_has_heme_and_fe_rows": sum(1 for row in found_rows if row["has_heme_and_fe"]),
            "trainable_unique_structure_index": len(trainable_rows),
            "trainable_drive_metadata_found_rows": len(found_trainable_rows),
            "trainable_download_exists_rows": sum(1 for row in found_trainable_rows if row["exists"]),
            "trainable_has_heme_rows": sum(1 for row in found_trainable_rows if row["has_heme"]),
            "trainable_has_fe_rows": sum(1 for row in found_trainable_rows if row["has_fe"]),
            "trainable_has_heme_and_fe_rows": sum(1 for row in found_trainable_rows if row["has_heme_and_fe"]),
        },
        "examples": {
            "missing_download": [row["pdb_name"] for row in found_rows if not row["exists"]][:20],
            "bad_size": [row["pdb_name"] for row in found_rows if row["exists"] and not row["size_match"]][:20],
            "has_heme_and_fe": [row["pdb_name"] for row in found_rows if row["has_heme_and_fe"]][:20],
        },
    }
    OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    try:
        LOCAL_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    except OSError as exc:
        summary["local_json_write_error"] = repr(exc)
        OUT_JSON.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
