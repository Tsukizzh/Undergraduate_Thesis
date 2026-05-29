#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import re
import subprocess
import sys
import csv
from datetime import datetime
from pathlib import Path

import torch


P450_ROOT = Path("/root/autodl-tmp/EZSpecificity/PathD/P450")
EXP_DIR = P450_ROOT / "experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay"
CACHE_DIR = (
    P450_ROOT
    / "data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1"
)
REPORT_DIR = EXP_DIR / "preflight_reports"
CONFIG = EXP_DIR / "configs/config_exp003_p450_subset_fe_overlay.yml"
RUN_SCRIPT = EXP_DIR / "scripts/run_exp003_p450_subset_fe_overlay_train.sh"
TRAIN_SCRIPT = EXP_DIR / "scripts/main_training_pt_fe_overlay.py"
DATASET_SCRIPT = EXP_DIR / "scripts/pt_dataset_fe_overlay.py"


def sh(cmd: list[str], timeout: int = 20) -> dict:
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return {
            "ok": proc.returncode == 0,
            "returncode": proc.returncode,
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip(),
        }
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "error": repr(exc)}


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def add_check(checks: list[dict], name: str, ok: bool, detail=None) -> None:
    checks.append({"name": name, "ok": bool(ok), "detail": detail})


def load_torch(path: Path):
    return torch.load(path, map_location="cpu", weights_only=False)


def sample_path(cache_dir: Path, split: str, sample_id: int) -> Path:
    bucket = sample_id // 1000
    return cache_dir / split / "samples" / f"{bucket:03d}" / f"sample_{sample_id:06d}.pt"


def main() -> int:
    now = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    checks: list[dict] = []
    report: dict = {
        "timestamp": now,
        "exp_dir": str(EXP_DIR),
        "cache_dir": str(CACHE_DIR),
        "config": str(CONFIG),
        "run_script": str(RUN_SCRIPT),
    }

    for path_name, path in [
        ("EXP003 experiment directory", EXP_DIR),
        ("EXP003 cache directory", CACHE_DIR),
        ("EXP003 config", CONFIG),
        ("EXP003 run script", RUN_SCRIPT),
        ("EXP003 training script", TRAIN_SCRIPT),
        ("EXP003 dataset script", DATASET_SCRIPT),
    ]:
        add_check(checks, f"{path_name} exists", path.exists(), str(path))

    if CONFIG.exists():
        config_text = read_text(CONFIG)
        tag_match = re.search(r"^\s*tag:\s*(.+?)\s*$", config_text, flags=re.M)
        report["config_tag"] = tag_match.group(1) if tag_match else None
        add_check(
            checks,
            "config tag is EXP003",
            report["config_tag"] == "Q1-EXP003-p450389-subset-fe-heme-overlay",
            report["config_tag"],
        )
    else:
        config_text = ""

    if RUN_SCRIPT.exists():
        run_text = read_text(RUN_SCRIPT)
        forbidden = ["EXP002", "exp002_fe_heme_structures", "00_EXP002", "Q1_EXP002"]
        residuals = {token: token in run_text for token in forbidden}
        report["run_script_forbidden_tokens"] = residuals
        add_check(
            checks,
            "run script has no EXP002 residue",
            not any(residuals.values()),
            residuals,
        )
        add_check(
            checks,
            "run script points to EXP003 cache",
            str(CACHE_DIR) in run_text,
            str(CACHE_DIR),
        )
        add_check(
            checks,
            "run script writes EXP003 results",
            "00_EXP003_P450_SUBSET_FE_OVERLAY" in run_text,
            "00_EXP003_P450_SUBSET_FE_OVERLAY",
        )
    else:
        run_text = ""

    py_compile = sh([sys.executable, "-m", "py_compile", str(TRAIN_SCRIPT), str(DATASET_SCRIPT)], timeout=60)
    report["py_compile"] = py_compile
    add_check(checks, "training and dataset scripts compile", py_compile.get("ok", False), py_compile)

    bash_check = sh(["bash", "-n", str(RUN_SCRIPT)], timeout=30)
    report["bash_n"] = bash_check
    add_check(checks, "run script bash syntax is valid", bash_check.get("ok", False), bash_check)

    manifest_path = CACHE_DIR / "manifest.pt"
    add_check(checks, "cache manifest exists", manifest_path.exists(), str(manifest_path))
    if manifest_path.exists():
        manifest = load_torch(manifest_path)
        audit_path = manifest.get("exp003_fe_overlay_audit")
        report["manifest_keys_subset"] = {
            "protein_elements": manifest.get("protein_elements"),
            "protein_has_is_hetero": manifest.get("protein_has_is_hetero"),
            "exp003_fe_overlay_from_exp002_rules": manifest.get("exp003_fe_overlay_from_exp002_rules"),
            "exp003_fe_overlay_sample_count": manifest.get("exp003_fe_overlay_sample_count"),
            "exp003_fe_overlay_audit": audit_path,
            "has_exp002_fe_overlay_key": "exp002_fe_overlay" in manifest,
        }
        add_check(checks, "manifest includes Fe element 26", 26 in (manifest.get("protein_elements") or []), manifest.get("protein_elements"))
        add_check(checks, "manifest has hetero flag", manifest.get("protein_has_is_hetero") is True, manifest.get("protein_has_is_hetero"))
        add_check(
            checks,
            "manifest has EXP003 overlay flag",
            manifest.get("exp003_fe_overlay_from_exp002_rules") is True,
            manifest.get("exp003_fe_overlay_from_exp002_rules"),
        )
        add_check(
            checks,
            "manifest overlay count is 1140",
            manifest.get("exp003_fe_overlay_sample_count") == 1140,
            manifest.get("exp003_fe_overlay_sample_count"),
        )
        add_check(checks, "manifest has no misleading EXP002 overlay key", "exp002_fe_overlay" not in manifest, "exp002_fe_overlay" in manifest)

        audit = Path(audit_path) if audit_path else None
        add_check(checks, "manifest points to overlay audit CSV", bool(audit_path), audit_path)
        add_check(checks, "overlay audit CSV exists", bool(audit and audit.exists()), str(audit) if audit else None)
        if audit and audit.exists():
            audit_rows = []
            with audit.open(newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    audit_rows.append(row)
            ok_rows = [row for row in audit_rows if row.get("status") == "ok"]
            unique_keys = {(row.get("split"), row.get("sample_id")) for row in ok_rows}
            report["overlay_audit_summary"] = {
                "rows": len(audit_rows),
                "ok_rows": len(ok_rows),
                "unique_ok_keys": len(unique_keys),
            }
            add_check(checks, "overlay audit has 1140 unique ok sample keys", len(unique_keys) == 1140, report["overlay_audit_summary"])

            sample_checks = []
            for row in ok_rows[:5]:
                split = row["split"]
                sample_id = int(row["sample_id"])
                path = sample_path(CACHE_DIR, split, sample_id)
                item = {"split": split, "sample_id": sample_id, "path": str(path), "exists": path.exists()}
                if path.exists():
                    data = load_torch(path)
                    hetero = data.get("protein_is_hetero")
                    item["has_protein_is_hetero"] = hetero is not None
                    try:
                        item["hetero_sum"] = int(hetero.sum().item()) if hetero is not None else None
                    except Exception as exc:  # noqa: BLE001
                        item["hetero_sum_error"] = repr(exc)
                sample_checks.append(item)
            report["overlay_sample_checks"] = sample_checks
            add_check(
                checks,
                "sampled overlay PT files contain hetero atoms",
                all(item.get("exists") and item.get("has_protein_is_hetero") and (item.get("hetero_sum") or 0) > 0 for item in sample_checks),
                sample_checks,
            )

    split_summary = {}
    expected_counts = {"train": 6263, "val": 607, "test": 1359}
    for split, expected_n in expected_counts.items():
        index_path = CACHE_DIR / split / "index.pt"
        add_check(checks, f"{split} index exists", index_path.exists(), str(index_path))
        if index_path.exists():
            index = load_torch(index_path)
            n = len(index.get("sample_ids", []))
            split_summary[split] = {"n": n}
            add_check(checks, f"{split} index count is expected", n == expected_n, {"actual": n, "expected": expected_n})
    report["split_summary"] = split_summary

    nvidia = sh(
        ["nvidia-smi", "--query-gpu=index,utilization.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"],
        timeout=20,
    )
    report["nvidia_smi"] = nvidia

    proc = sh(
        [
            "bash",
            "-lc",
            "ps -eo pid,ppid,cmd | grep -E 'EXP002|EXP003|main_training_pt_fe_overlay|torchrun' "
            "| grep -v grep | grep -v q01_exp003_preflight_check_20260527.py || true",
        ],
        timeout=20,
    )
    report["training_processes"] = proc
    proc_lines = [line for line in proc.get("stdout", "").splitlines() if line.strip()]
    exp003_training_lines = [
        line
        for line in proc_lines
        if (
            ("Q1_EXP003" in line)
            or ("00_EXP003" in line)
            or ("EXP003_full_p450_fe_heme_overlay" in line and "main_training_pt_fe_overlay.py" in line)
        )
    ]
    report["exp003_training_process_lines"] = exp003_training_lines
    add_check(checks, "EXP003 training process is not already running", not exp003_training_lines, exp003_training_lines)

    report["checks"] = checks
    report["ok"] = all(item["ok"] for item in checks)

    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    report_path = REPORT_DIR / f"exp003_preflight_{now}.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print(json.dumps({"ok": report["ok"], "report_path": str(report_path), "failed": [c for c in checks if not c["ok"]]}, ensure_ascii=False, indent=2))
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
