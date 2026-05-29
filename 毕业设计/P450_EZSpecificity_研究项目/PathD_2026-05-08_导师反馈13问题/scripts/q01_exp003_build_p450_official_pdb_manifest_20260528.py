#!/usr/bin/env python3
"""Build the EXP003 P450 official ESIBank PDB download manifest.

This script does not read the mounted G: drive. It reads Google Drive Desktop's
local metadata SQLite cache and maps P450 ESIBank structure_index values to
the official brenda/structure/complex/{structure_index}.pdb file ids.
"""

from __future__ import annotations

import csv
import json
import re
import sqlite3
from collections import defaultdict
from datetime import datetime
from pathlib import Path


PROJECT = Path(r"D:\EZSpecificity_Project")
ESIBANK = Path(r"D:\ESIBank")
PATHD = PROJECT / r"毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题"

P450_CSV = PROJECT / r"毕业设计\提取P450过程日志\2026-01-02_01-46_P450精确验证\数据\P450酶列表_最终版389个.csv"
ENZYMES_CSV = ESIBANK / r"brenda\enzymes.csv"
DATA_CSV = ESIBANK / r"brenda\data.csv"
TRAINABLE_MANIFEST = PATHD / r"data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\manifests\exp003_p450389_trainable_split_sample_manifest_20260527_214331.csv"

OUT_DIR = PATHD / r"data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\official_pdb_manifest_20260528"
OUT_CSV = OUT_DIR / "exp003_p450_official_esibank_pdb_manifest_20260528.csv"
OUT_SUMMARY = OUT_DIR / "exp003_p450_official_esibank_pdb_manifest_summary_20260528.json"

DRIVEFS_METADATA = Path(
    r"C:\Users\Administrator\AppData\Local\Google\DriveFS\112184715553175865523\metadata_sqlite_db"
)
COMPLEX_FOLDER_ID = "133tzvNXYmGfelJve0VsUzfcsvQEw9YkC"
COMPLEX_DRIVE_PATH = "ESIBank/brenda/structure/complex"


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def as_int(value: object, default: int = -1) -> int:
    try:
        text = str(value).strip()
        if not text:
            return default
        return int(float(text))
    except Exception:
        return default


def load_p450_structure_rows() -> tuple[dict[int, dict[str, object]], dict[int, dict[str, object]]]:
    p450_rows = read_csv(P450_CSV)
    p450_ids = {row["uniprot_id"].strip() for row in p450_rows if row.get("uniprot_id", "").strip()}

    enzyme_rows = read_csv(ENZYMES_CSV)
    enzyme_to_uniprot = {
        idx: row.get("uniprots", "").strip()
        for idx, row in enumerate(enzyme_rows)
    }

    by_structure: dict[int, dict[str, object]] = {}
    split_counts = defaultdict(lambda: {"rows": 0, "positive": 0, "enzymes": set()})
    for row in read_csv(DATA_CSV):
        enzyme_index = as_int(row.get("enzyme"))
        uniprot = enzyme_to_uniprot.get(enzyme_index, "")
        if uniprot not in p450_ids:
            continue
        structure_index = as_int(row.get("structure_index"))
        if structure_index < 0:
            continue
        record = by_structure.setdefault(
            structure_index,
            {
                "structure_index": structure_index,
                "enzyme_indices": set(),
                "uniprots": set(),
                "data_rows": 0,
                "positive_rows": 0,
            },
        )
        record["enzyme_indices"].add(enzyme_index)
        record["uniprots"].add(uniprot)
        record["data_rows"] += 1
        if as_int(row.get("label")) == 1:
            record["positive_rows"] += 1

    trainable_by_structure: dict[int, dict[str, object]] = {}
    if TRAINABLE_MANIFEST.exists():
        for row in read_csv(TRAINABLE_MANIFEST):
            structure_index = as_int(row.get("dock_index"))
            if structure_index < 0:
                continue
            record = trainable_by_structure.setdefault(
                structure_index,
                {
                    "structure_index": structure_index,
                    "trainable_rows": 0,
                    "trainable_positive_rows": 0,
                    "trainable_splits": set(),
                },
            )
            record["trainable_rows"] += 1
            record["trainable_splits"].add(str(row.get("split", "")).strip())
            if as_int(row.get("label")) == 1:
                record["trainable_positive_rows"] += 1

    return by_structure, trainable_by_structure


def extract_pdb_name(proto: bytes) -> str:
    match = re.search(rb"([0-9]+\.pdb)\"", proto or b"")
    if not match:
        return ""
    return match.group(1).decode("ascii", "ignore")


def load_drive_complex_index() -> tuple[dict[str, dict[str, object]], dict[str, object]]:
    if not DRIVEFS_METADATA.exists():
        raise FileNotFoundError(DRIVEFS_METADATA)

    con = sqlite3.connect(f"file:{DRIVEFS_METADATA}?mode=ro", uri=True)
    con.text_factory = bytes
    try:
        row = con.execute("select stable_id from items where id=?", (COMPLEX_FOLDER_ID,)).fetchone()
        if not row:
            raise RuntimeError(f"complex folder id not found in DriveFS metadata: {COMPLEX_FOLDER_ID}")
        parent_stable_id = int(row[0])

        query = """
            select i.stable_id, i.id, i.file_size, i.proto
            from items i
            join stable_parents p on i.stable_id = p.item_stable_id
            where p.parent_stable_id = ? and i.is_folder = 0
        """
        by_name: dict[str, dict[str, object]] = {}
        duplicate_names: list[str] = []
        total_children = 0
        pdb_children = 0
        for stable_id, file_id, file_size, proto in con.execute(query, (parent_stable_id,)):
            total_children += 1
            name = extract_pdb_name(proto)
            if not name:
                continue
            pdb_children += 1
            if name in by_name:
                duplicate_names.append(name)
                continue
            by_name[name] = {
                "pdb_name": name,
                "drive_file_id": file_id.decode("ascii", "ignore"),
                "stable_id": int(stable_id),
                "file_size": int(file_size or 0),
                "drive_path": f"{COMPLEX_DRIVE_PATH}/{name}",
            }
    finally:
        con.close()

    summary = {
        "drivefs_metadata": str(DRIVEFS_METADATA),
        "complex_folder_id": COMPLEX_FOLDER_ID,
        "complex_parent_stable_id": parent_stable_id,
        "total_children_in_complex_parent": total_children,
        "pdb_children_with_parsed_name": pdb_children,
        "unique_pdb_names": len(by_name),
        "duplicate_pdb_name_count": len(duplicate_names),
        "duplicate_pdb_name_examples": duplicate_names[:20],
    }
    return by_name, summary


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    p450_structures, trainable_structures = load_p450_structure_rows()
    drive_index, drive_summary = load_drive_complex_index()

    rows: list[dict[str, object]] = []
    for structure_index in sorted(p450_structures):
        info = p450_structures[structure_index]
        trainable = trainable_structures.get(structure_index, {})
        pdb_name = f"{structure_index}.pdb"
        drive = drive_index.get(pdb_name, {})
        found = bool(drive)
        rows.append(
            {
                "structure_index": structure_index,
                "pdb_name": pdb_name,
                "drive_file_id": drive.get("drive_file_id", ""),
                "stable_id": drive.get("stable_id", ""),
                "file_size": drive.get("file_size", ""),
                "found_in_drivefs_metadata": found,
                "found_in_drive_metadata": found,
                "drive_path": drive.get("drive_path", f"{COMPLEX_DRIVE_PATH}/{pdb_name}"),
                "p450_data_rows": info["data_rows"],
                "p450_positive_rows": info["positive_rows"],
                "p450_unique_enzyme_count": len(info["enzyme_indices"]),
                "p450_uniprots": ";".join(sorted(info["uniprots"])),
                "is_trainable_in_exp003_manifest": structure_index in trainable_structures,
                "trainable_rows": trainable.get("trainable_rows", 0),
                "trainable_positive_rows": trainable.get("trainable_positive_rows", 0),
                "trainable_splits": ";".join(sorted(trainable.get("trainable_splits", []))),
            }
        )

    fieldnames = [
        "structure_index",
        "pdb_name",
        "drive_file_id",
        "stable_id",
        "file_size",
        "found_in_drivefs_metadata",
        "found_in_drive_metadata",
        "drive_path",
        "p450_data_rows",
        "p450_positive_rows",
        "p450_unique_enzyme_count",
        "p450_uniprots",
        "is_trainable_in_exp003_manifest",
        "trainable_rows",
        "trainable_positive_rows",
        "trainable_splits",
    ]
    write_csv(OUT_CSV, rows, fieldnames)

    found_rows = [row for row in rows if row["found_in_drivefs_metadata"]]
    trainable_rows = [row for row in rows if row["is_trainable_in_exp003_manifest"]]
    trainable_found = [row for row in trainable_rows if row["found_in_drivefs_metadata"]]
    summary = {
        "created_at": now(),
        "purpose": "EXP003 P450 official ESIBank brenda/structure/complex PDB manifest from DriveFS metadata",
        "inputs": {
            "p450_csv": str(P450_CSV),
            "enzymes_csv": str(ENZYMES_CSV),
            "data_csv": str(DATA_CSV),
            "trainable_manifest": str(TRAINABLE_MANIFEST),
        },
        "outputs": {
            "manifest_csv": str(OUT_CSV),
            "summary_json": str(OUT_SUMMARY),
        },
        "drivefs": drive_summary,
        "counts": {
            "p450_unique_structure_index": len(rows),
            "p450_found_in_drivefs_metadata": len(found_rows),
            "p450_missing_in_drivefs_metadata": len(rows) - len(found_rows),
            "p450_expected_download_gib": round(
                sum(int(row.get("file_size") or 0) for row in found_rows) / 1024**3,
                3,
            ),
            "trainable_unique_structure_index": len(trainable_rows),
            "trainable_found_in_drivefs_metadata": len(trainable_found),
            "trainable_missing_in_drivefs_metadata": len(trainable_rows) - len(trainable_found),
        },
        "missing_examples": [
            row["pdb_name"]
            for row in rows
            if not row["found_in_drivefs_metadata"]
        ][:20],
    }
    OUT_SUMMARY.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
