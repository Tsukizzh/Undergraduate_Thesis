#!/usr/bin/env python3
"""Download selected ESIBank PDB files from Google Drive by file id.

The script is intentionally append-only/resumable: existing files with the
expected size are skipped, partial files are stored as .part files, and a
progress JSON is rewritten after each completed item.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

import requests


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_manifest(path: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            found = str(row.get("found_in_drive_metadata", "")).strip().lower()
            file_id = str(row.get("drive_file_id", "")).strip()
            name = str(row.get("pdb_name", "")).strip()
            if found in {"true", "1", "yes"} and file_id and name:
                rows.append(row)
    return rows


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def looks_like_pdb(path: Path) -> bool:
    try:
        data = path.read_bytes()[:8192]
    except OSError:
        return False
    return data.startswith(b"HEADER") or b"ATOM" in data or b"HETATM" in data or data.startswith(b"REMARK")


def download_one(row: Dict[str, str], out_dir: Path, retries: int, timeout: int) -> Dict[str, object]:
    name = row["pdb_name"]
    file_id = row["drive_file_id"]
    expected_size = int(float(row.get("file_size") or 0))
    out = out_dir / name
    part = out_dir / f"{name}.part"
    url = f"https://drive.google.com/uc?export=download&id={file_id}"

    if out.exists() and expected_size > 0 and out.stat().st_size == expected_size and looks_like_pdb(out):
        return {
            "pdb_name": name,
            "drive_file_id": file_id,
            "status": "skipped_ok",
            "bytes": out.stat().st_size,
            "sha256": "",
            "message": "already exists with expected size",
        }

    if out.exists() and expected_size > 0 and out.stat().st_size != expected_size:
        out.rename(part)

    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            with requests.Session() as session:
                with session.get(url, stream=True, timeout=timeout, allow_redirects=True) as response:
                    response.raise_for_status()
                    content_type = response.headers.get("content-type", "")
                    tmp = part
                    with tmp.open("wb") as f:
                        for chunk in response.iter_content(chunk_size=1024 * 1024):
                            if chunk:
                                f.write(chunk)
                    size = tmp.stat().st_size
                    if expected_size > 0 and size != expected_size:
                        raise RuntimeError(f"size mismatch: got {size}, expected {expected_size}")
                    if not looks_like_pdb(tmp):
                        head = tmp.read_bytes()[:120].decode("utf-8", "ignore")
                        raise RuntimeError(f"downloaded file does not look like PDB, content_type={content_type}, head={head!r}")
                    tmp.replace(out)
                    return {
                        "pdb_name": name,
                        "drive_file_id": file_id,
                        "status": "downloaded",
                        "bytes": size,
                        "sha256": sha256_file(out),
                        "message": "",
                    }
        except Exception as exc:  # noqa: BLE001
            last_error = repr(exc)
            time.sleep(min(30, 2 ** attempt) + random.random())

    return {
        "pdb_name": name,
        "drive_file_id": file_id,
        "status": "failed",
        "bytes": part.stat().st_size if part.exists() else 0,
        "sha256": "",
        "message": last_error,
    }


def atomic_write_json(path: Path, payload: Dict[str, object]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--retries", type=int, default=5)
    parser.add_argument("--timeout", type=int, default=90)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--progress-json", type=Path, default=None)
    parser.add_argument("--audit-csv", type=Path, default=None)
    args = parser.parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(args.manifest)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_manifest(args.manifest)
    if args.limit:
        rows = rows[: args.limit]

    progress_json = args.progress_json or args.out_dir.parent / "download_progress.json"
    audit_csv = args.audit_csv or args.out_dir.parent / "download_audit.csv"

    start = time.time()
    results: Dict[str, Dict[str, object]] = {}
    completed = 0
    total = len(rows)
    print(f"[{now()}] manifest rows to download: {total}")
    print(f"[{now()}] output dir: {args.out_dir}")

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        future_map = {
            executor.submit(download_one, row, args.out_dir, args.retries, args.timeout): row["pdb_name"]
            for row in rows
        }
        for future in as_completed(future_map):
            name = future_map[future]
            try:
                result = future.result()
            except Exception as exc:  # noqa: BLE001
                result = {
                    "pdb_name": name,
                    "drive_file_id": "",
                    "status": "failed",
                    "bytes": 0,
                    "sha256": "",
                    "message": repr(exc),
                }
            results[str(result["pdb_name"])] = result
            completed += 1
            if completed % 50 == 0 or completed == total:
                counts: Dict[str, int] = {}
                for r in results.values():
                    counts[str(r["status"])] = counts.get(str(r["status"]), 0) + 1
                payload = {
                    "updated_at": now(),
                    "manifest": str(args.manifest),
                    "out_dir": str(args.out_dir),
                    "completed": completed,
                    "total": total,
                    "counts": counts,
                    "elapsed_sec": round(time.time() - start, 1),
                }
                atomic_write_json(progress_json, payload)
                print(f"[{now()}] {completed}/{total} {counts}", flush=True)

    with audit_csv.open("w", encoding="utf-8", newline="") as f:
        fieldnames = ["pdb_name", "drive_file_id", "status", "bytes", "sha256", "message"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(results):
            writer.writerow({k: results[key].get(k, "") for k in fieldnames})

    counts: Dict[str, int] = {}
    total_bytes = 0
    for r in results.values():
        counts[str(r["status"])] = counts.get(str(r["status"]), 0) + 1
        total_bytes += int(r.get("bytes") or 0)
    final = {
        "finished_at": now(),
        "manifest": str(args.manifest),
        "out_dir": str(args.out_dir),
        "audit_csv": str(audit_csv),
        "total": total,
        "counts": counts,
        "bytes": total_bytes,
        "gib": round(total_bytes / 1024**3, 3),
        "elapsed_sec": round(time.time() - start, 1),
        "free_space_out_drive_gib": round(shutil.disk_usage(args.out_dir).free / 1024**3, 3),
    }
    atomic_write_json(progress_json, final)
    print(json.dumps(final, ensure_ascii=False, indent=2))
    failed = counts.get("failed", 0)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
