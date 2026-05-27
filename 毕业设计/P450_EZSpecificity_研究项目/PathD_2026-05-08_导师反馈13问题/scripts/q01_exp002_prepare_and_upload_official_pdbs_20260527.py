#!/usr/bin/env python3
"""Prepare EXP002 official ESIBank PDB subset and upload it to the server."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


PROJECT = Path(r"D:\EZSpecificity_Project")
SCRIPTS = PROJECT / r"毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\scripts"
MANIFEST_DIR = PROJECT / r".q1_exp002_manifests\manifests_20260527"
MANIFEST_CSV = MANIFEST_DIR / "exp002_google_drive_pdb_manifest_for_rcsb_best1_20260527.csv"
MANIFEST_JSON = MANIFEST_DIR / "exp002_google_drive_pdb_manifest_for_rcsb_best1_20260527.json"
LOCAL_ROOT = Path(r"D:\ESIBank\EXP002_official_brenda_pdb_best1_20260527")
LOCAL_PDB_DIR = LOCAL_ROOT / "pdb"
LOCAL_MANIFEST_DIR = LOCAL_ROOT / "manifests"
ARCHIVE_DIR = LOCAL_ROOT.parent / "archives"
ARCHIVE = ARCHIVE_DIR / "EXP002_official_brenda_pdb_best1_20260527_package_v2.tar.gz"
DOWNLOAD_SCRIPT = SCRIPTS / "q01_exp002_download_drive_pdbs_20260527.py"

SSH_SKILL = Path(r"C:\Users\Administrator\.codex\skills\ssh-skill\scripts")
SSH_UPLOAD = SSH_SKILL / "ssh_upload.py"
SSH_EXEC = SSH_SKILL / "ssh_execute.py"
SERVER = "autodl-4x5090-bj"
REMOTE_BASE = "/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures"
REMOTE_ARCHIVE = f"{REMOTE_BASE}/archives/EXP002_official_brenda_pdb_best1_20260527_package_v2.tar.gz"
REMOTE_TARGET = f"{REMOTE_BASE}/official_esibank_brenda_pdb_best1_20260527"


def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def run(cmd: list[str], *, env: dict[str, str] | None = None, timeout: int | None = None) -> None:
    log("RUN " + " ".join(f'"{c}"' if " " in str(c) else str(c) for c in cmd))
    subprocess.run(cmd, check=True, env=env, timeout=timeout)


def read_manifest_counts() -> dict[str, int | float]:
    found = 0
    missing = 0
    bytes_found = 0
    with MANIFEST_CSV.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            is_found = str(row.get("found_in_drive_metadata", "")).strip().lower() in {"true", "1", "yes"}
            if is_found:
                found += 1
                bytes_found += int(float(row.get("file_size") or 0))
            else:
                missing += 1
    return {
        "found": found,
        "missing": missing,
        "bytes_found": bytes_found,
        "gib_found": round(bytes_found / 1024**3, 3),
    }


def local_audit(expected_found: int) -> dict[str, object]:
    files = sorted(LOCAL_PDB_DIR.glob("*.pdb"))
    total_bytes = sum(p.stat().st_size for p in files)
    bad = []
    for p in files[:]:
        head = p.read_bytes()[:8192]
        if not (head.startswith(b"HEADER") or b"ATOM" in head or b"HETATM" in head or head.startswith(b"REMARK")):
            bad.append(p.name)
            if len(bad) >= 20:
                break
    return {
        "pdb_files": len(files),
        "expected_found": expected_found,
        "bytes": total_bytes,
        "gib": round(total_bytes / 1024**3, 3),
        "bad_pdb_head_sample": bad,
        "complete": len(files) == expected_found and not bad,
    }


def write_readme(counts: dict[str, object], audit: dict[str, object]) -> None:
    text = f"""# EXP002 official ESIBank PDB subset

Generated at: {now()}

Purpose:

- Prepare official ESIBank `brenda/structure/complex/*.pdb` files that match the EXP002 RCSB best1 structure-index list.
- Keep these files separate from the original `D:\\ESIBank\\brenda` directory.

Source manifest:

- `{MANIFEST_CSV}`

Local directory:

- `{LOCAL_ROOT}`

Server target:

- `{REMOTE_TARGET}`

Counts:

| item | value |
|---|---:|
| requested structure_index rows | {counts['found'] + counts['missing']} |
| found in Google Drive brenda/structure/complex | {counts['found']} |
| missing from brenda/structure/complex | {counts['missing']} |
| expected downloaded size | {counts['gib_found']} GiB |
| downloaded PDB files | {audit['pdb_files']} |
| downloaded size | {audit['gib']} GiB |

Notes:

- The 588 missing items were not found under the official `brenda/structure/complex` folder by Google Drive metadata.
- Same-name PDB files in `small_family` were intentionally not used as substitutes.
- RCSB HEM/Fe mmCIF files are already stored separately on the server under `raw_mmcif_rcsb_heme_fe_v1`.
"""
    (LOCAL_ROOT / "README_EXP002_official_pdb_package.md").write_text(text, encoding="utf-8")


def make_archive() -> None:
    if ARCHIVE.exists():
        archive_mtime = ARCHIVE.stat().st_mtime
        newest = max((p.stat().st_mtime for p in LOCAL_ROOT.rglob("*") if p.is_file() and p != ARCHIVE), default=0)
        if archive_mtime >= newest:
            log(f"archive exists and is fresh: {ARCHIVE}")
            return
        ARCHIVE.unlink()
    tar = shutil.which("tar")
    if not tar:
        raise RuntimeError("tar executable not found on Windows PATH")
    members = [
        "pdb",
        "manifests",
        "README_EXP002_official_pdb_package.md",
        "local_audit_20260527.json",
        "download_progress.json",
        "download_audit.csv",
    ]
    run([
        tar,
        "-czf",
        str(ARCHIVE),
        "-C",
        str(LOCAL_ROOT),
        *members,
    ])
    log(f"archive created: {ARCHIVE} ({ARCHIVE.stat().st_size / 1024**3:.3f} GiB)")


def upload_and_extract() -> None:
    env = os.environ.copy()
    env["MSYS_NO_PATHCONV"] = "1"
    remote_setup = f"mkdir -p '{REMOTE_BASE}/archives' '{REMOTE_TARGET}' && df -h /root/autodl-tmp"
    run([sys.executable, str(SSH_EXEC), SERVER, remote_setup, "--timeout", "120", "--no-daemon"])
    run([sys.executable, str(SSH_UPLOAD), SERVER, str(ARCHIVE), REMOTE_ARCHIVE, "--no-progress"], env=env)
    remote_extract = f"""
set -e
BASE='{REMOTE_BASE}'
TARGET='{REMOTE_TARGET}'
ARCHIVE='{REMOTE_ARCHIVE}'
mkdir -p "$TARGET"
tar -xzf "$ARCHIVE" -C "$TARGET"
if [ ! -d "$TARGET/pdb" ]; then
  echo "ERROR: target pdb dir missing after extraction" >&2
  exit 2
fi
PDB_COUNT=$(find "$TARGET/pdb" -maxdepth 1 -type f -name '*.pdb' | wc -l)
PDB_BYTES=$(find "$TARGET/pdb" -maxdepth 1 -type f -name '*.pdb' -printf '%s\\n' | awk '{{s+=$1}} END{{print s+0}}')
cat > "$TARGET/server_audit_20260527.json" <<EOF
{{
  "created_at": "$(date '+%F %T')",
  "target": "$TARGET",
  "archive": "$ARCHIVE",
  "pdb_count": $PDB_COUNT,
  "pdb_bytes": $PDB_BYTES
}}
EOF
cat "$TARGET/server_audit_20260527.json"
df -h /root/autodl-tmp | tail -1
"""
    run([sys.executable, str(SSH_EXEC), SERVER, remote_extract, "--timeout", "3600", "--no-daemon"])


def main() -> int:
    if not MANIFEST_CSV.exists():
        raise FileNotFoundError(MANIFEST_CSV)
    if not DOWNLOAD_SCRIPT.exists():
        raise FileNotFoundError(DOWNLOAD_SCRIPT)
    LOCAL_ROOT.mkdir(parents=True, exist_ok=True)
    LOCAL_PDB_DIR.mkdir(parents=True, exist_ok=True)
    LOCAL_MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MANIFEST_CSV, LOCAL_MANIFEST_DIR / MANIFEST_CSV.name)
    if MANIFEST_JSON.exists():
        shutil.copy2(MANIFEST_JSON, LOCAL_MANIFEST_DIR / MANIFEST_JSON.name)

    counts = read_manifest_counts()
    log(f"manifest counts: {json.dumps(counts, ensure_ascii=False)}")
    progress_json = LOCAL_ROOT / "download_progress.json"
    audit_csv = LOCAL_ROOT / "download_audit.csv"
    run([
        sys.executable,
        str(DOWNLOAD_SCRIPT),
        "--manifest", str(MANIFEST_CSV),
        "--out-dir", str(LOCAL_PDB_DIR),
        "--workers", "12",
        "--retries", "6",
        "--timeout", "120",
        "--progress-json", str(progress_json),
        "--audit-csv", str(audit_csv),
    ])

    audit = local_audit(int(counts["found"]))
    (LOCAL_ROOT / "local_audit_20260527.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    write_readme(counts, audit)
    log(f"local audit: {json.dumps(audit, ensure_ascii=False)}")
    if not audit["complete"]:
        raise RuntimeError("local audit failed; refusing to upload incomplete package")

    make_archive()
    upload_and_extract()
    log("EXP002 official PDB subset prepared and uploaded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
