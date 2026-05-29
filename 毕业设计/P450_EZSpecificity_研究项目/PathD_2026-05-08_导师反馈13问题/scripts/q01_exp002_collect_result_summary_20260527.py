#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

import torch


RESULT_DIR = Path(
    "/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/"
    "EXP002_fe_heme_overlay/results/00_EXP002_FE_OVERLAY_20260527_184505"
)
CKPT_DIR = RESULT_DIR / "checkpoints"
METRICS_CSV = RESULT_DIR / "metrics.csv"
TEST_EVAL = RESULT_DIR / "test_eval.json"


def run(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def read_metrics() -> list[dict]:
    if not METRICS_CSV.exists():
        return []
    with METRICS_CSV.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def to_float(value):
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def read_checkpoint(path: Path) -> dict:
    ckpt = torch.load(path, map_location="cpu", weights_only=False)
    current_score = None
    best_score = None
    best_model_path = None
    for value in ckpt.get("callbacks", {}).values():
        if not isinstance(value, dict):
            continue
        if "current_score" in value:
            raw = value["current_score"]
            try:
                current_score = float(raw.cpu())
            except Exception:  # noqa: BLE001
                current_score = str(raw)
        if "best_model_score" in value:
            raw = value["best_model_score"]
            try:
                best_score = float(raw.cpu())
            except Exception:  # noqa: BLE001
                best_score = str(raw)
        if "best_model_path" in value:
            best_model_path = value["best_model_path"]
    return {
        "file": path.name,
        "path": str(path),
        "epoch": ckpt.get("epoch"),
        "global_step": ckpt.get("global_step"),
        "current_auc_val": current_score,
        "best_auc_val_so_far": best_score,
        "best_model_path_recorded": best_model_path,
        "mtime": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
        "size_bytes": path.stat().st_size,
    }


def read_checkpoints() -> list[dict]:
    if not CKPT_DIR.exists():
        return []
    rows = []
    for path in sorted(CKPT_DIR.glob("*.ckpt")):
        rows.append(read_checkpoint(path))
    return rows


def choose_best_checkpoint(rows: list[dict]) -> dict | None:
    scored = [row for row in rows if isinstance(row.get("current_auc_val"), float)]
    if not scored:
        return None
    best = max(scored, key=lambda row: row["current_auc_val"])
    recorded_path = best.get("best_model_path_recorded")
    if recorded_path:
        recorded_name = Path(recorded_path).name
        for row in rows:
            if row.get("file") == recorded_name:
                row = dict(row)
                row["selected_by"] = "best_model_path_recorded"
                return row
    named = [row for row in scored if row.get("file") != "last.ckpt"]
    if named:
        best_named = max(named, key=lambda row: row["current_auc_val"])
        best_named = dict(best_named)
        best_named["selected_by"] = "max_current_auc_non_last"
        return best_named
    best = dict(best)
    best["selected_by"] = "max_current_auc"
    return best


def read_test_eval() -> dict | None:
    if not TEST_EVAL.exists():
        return None
    return json.loads(TEST_EVAL.read_text(encoding="utf-8"))


def write_markdown(summary: dict, path: Path) -> None:
    lines = []
    lines.append("# Q1 EXP002 FE/HEM Overlay Result Summary")
    lines.append("")
    lines.append(f"- Generated at: `{summary['generated_at']}`")
    lines.append(f"- Result dir: `{summary['result_dir']}`")
    lines.append(f"- Training process still running: `{summary['training_process_running']}`")
    lines.append(f"- Test eval exists: `{summary['test_eval_exists']}`")
    lines.append("")
    best = summary.get("best_checkpoint")
    if best:
        lines.append("## Best Checkpoint")
        lines.append("")
        lines.append(f"- File: `{best['file']}`")
        lines.append(f"- Epoch: `{best['epoch']}`")
        lines.append(f"- Global step: `{best['global_step']}`")
        lines.append(f"- Checkpoint `auc/val`: `{best['current_auc_val']}`")
        lines.append("")
    latest = summary.get("latest_metrics_row")
    if latest:
        lines.append("## Latest Metrics CSV Row")
        lines.append("")
        lines.append("| epoch | train_loss | val_loss | val_auc | val_aupr |")
        lines.append("|---:|---:|---:|---:|---:|")
        lines.append(
            f"| {latest.get('epoch')} | {latest.get('train_loss')} | {latest.get('val_loss')} | "
            f"{latest.get('val_auc')} | {latest.get('val_aupr')} |"
        )
        lines.append("")
    test_eval = summary.get("test_eval")
    if test_eval:
        lines.append("## Test Eval")
        lines.append("")
        lines.append(f"- checkpoint: `{test_eval.get('checkpoint')}`")
        lines.append(f"- test_auc_roc: `{test_eval.get('test_auc_roc')}`")
        lines.append(f"- test_aupr: `{test_eval.get('test_aupr')}`")
        lines.append("")
    else:
        lines.append("## Test Eval")
        lines.append("")
        lines.append("`test_eval.json` has not been generated yet.")
        lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    now = datetime.now().strftime("%Y%m%d_%H%M%S")
    metrics = read_metrics()
    ckpts = read_checkpoints()
    test_eval = read_test_eval()
    proc = run(["bash", "-lc", "pgrep -af 'Q1_EXP002|EXP002_fe_heme_overlay' | grep -v grep || true"])
    nvidia = run(["nvidia-smi", "--query-gpu=index,utilization.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"])
    best_ckpt = choose_best_checkpoint(ckpts)
    summary = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "result_dir": str(RESULT_DIR),
        "metrics_csv": str(METRICS_CSV),
        "metrics_rows": len(metrics),
        "latest_metrics_row": metrics[-1] if metrics else None,
        "checkpoint_count": len(ckpts),
        "checkpoints": ckpts,
        "best_checkpoint": best_ckpt,
        "test_eval_exists": test_eval is not None,
        "test_eval": test_eval,
        "training_process_running": bool(proc["stdout"].strip()),
        "training_processes": proc,
        "nvidia_smi": nvidia,
    }
    out_json = RESULT_DIR / f"result_summary_{now}.json"
    out_md = RESULT_DIR / f"result_summary_{now}.md"
    latest_json = RESULT_DIR / "result_summary_latest.json"
    latest_md = RESULT_DIR / "result_summary_latest.md"
    out_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    latest_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_markdown(summary, out_md)
    write_markdown(summary, latest_md)
    print(json.dumps({"summary_json": str(out_json), "summary_md": str(out_md), "best_checkpoint": best_ckpt, "test_eval_exists": test_eval is not None}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
