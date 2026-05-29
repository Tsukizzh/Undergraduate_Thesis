#!/usr/bin/env bash
set -u

ROOT=/root/autodl-tmp/EZSpecificity/PathD/P450
EXP=${ROOT}/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay
RESULT_DIR=${EXP}/results/00_EXP003_FULL_AF_HEME_FE_20260528_055549
RUN_NAME=Q1_EXP003_full_af_heme_fe_ddp2_b40_w4_20260528_055549
EVAL_SCRIPT=${EXP}/scripts/q01_exp003_eval_full_and_p450389_20260528.sh
LOG_DIR=${EXP}/logs
LOG=${LOG_DIR}/exp003_finish_watch_eval_shutdown_$(date +%Y%m%d_%H%M%S).log
SLEEP_SECONDS=${SLEEP_SECONDS:-300}
MAX_WAIT_SECONDS=${MAX_WAIT_SECONDS:-86400}

mkdir -p "${LOG_DIR}"
exec >>"${LOG}" 2>&1

echo "======================================================================"
echo "[watcher] started: $(date)"
echo "[watcher] run name: ${RUN_NAME}"
echo "[watcher] result dir: ${RESULT_DIR}"
echo "[watcher] eval script: ${EVAL_SCRIPT}"
echo "======================================================================"

start_ts=$(date +%s)

while true; do
  if pgrep -af "main_training_pt_fe_overlay.py.*${RUN_NAME}" >/dev/null; then
    now_ts=$(date +%s)
    elapsed=$((now_ts - start_ts))
    echo "[watcher] $(date) training still running, elapsed=${elapsed}s"
    if [ "${elapsed}" -gt "${MAX_WAIT_SECONDS}" ]; then
      echo "[watcher][ERROR] wait timeout; leaving server on for manual inspection"
      exit 2
    fi
    sleep "${SLEEP_SECONDS}"
  else
    echo "[watcher] $(date) training process not found; checking completion evidence"
    break
  fi
done

/root/miniconda3/bin/python - <<'PY'
from pathlib import Path
import csv
import re
import sys

result_dir = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/results/00_EXP003_FULL_AF_HEME_FE_20260528_055549")
metrics = result_dir / "metrics.csv"
ckpt_dir = result_dir / "checkpoints"

if not metrics.exists():
    raise SystemExit("[watcher][ERROR] metrics.csv missing")
rows = list(csv.DictReader(metrics.open()))
if not rows:
    raise SystemExit("[watcher][ERROR] metrics.csv empty")
epochs = sorted({int(r["epoch"]) for r in rows if r.get("epoch", "").isdigit()})
print("[watcher] metric epochs:", epochs[-10:])
if not epochs or max(epochs) < 34:
    raise SystemExit(f"[watcher][ERROR] training not complete; max metric epoch={max(epochs) if epochs else 'NA'}")

scored = []
for p in ckpt_dir.glob("*.ckpt"):
    if p.name == "last.ckpt":
        continue
    m_auc = re.search(r"auc([0-9.]+)\.ckpt$", p.name)
    m_ep = re.search(r"-ep(\d+)-", p.name)
    if m_auc and m_ep:
        scored.append((float(m_auc.group(1)), int(m_ep.group(1)), p))
if not scored:
    raise SystemExit("[watcher][ERROR] no scored checkpoint")
best = max(scored, key=lambda x: (x[0], x[1]))
print("[watcher] best checkpoint:", best[2])
print("[watcher] best val auc:", best[0], "epoch:", best[1])
PY

echo "[watcher] completion evidence passed; running final EXP003 evaluation"
bash "${EVAL_SCRIPT}"

/root/miniconda3/bin/python - <<'PY'
from pathlib import Path
import json
from datetime import datetime

root = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/results/q01_fe_embedding_patch")
session_log = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/sessions/q01_fe_embedding_patch/session_log.md")
dirs = sorted(root.glob("exp003_eval_*"), key=lambda p: p.stat().st_mtime, reverse=True)
if not dirs:
    raise SystemExit("[watcher][ERROR] no exp003 eval directory")
out = dirs[0]
summary_path = out / "EXP003_eval_summary.json"
full_path = out / "EXP003_full_test.json"
p450_path = out / "EXP003_p450389_test.json"
for p in [summary_path, full_path, p450_path]:
    if not p.exists():
        raise SystemExit(f"[watcher][ERROR] missing {p}")

summary = json.loads(summary_path.read_text())
full = summary["full_test"]
p450 = summary["p450389_test"]
checks = [
    ("full n_samples", full["n_samples"], 53588),
    ("p450389 n_samples", p450["n_samples"], 1359),
    ("p450389 n_positive", p450["n_positive"], 106),
    ("p450389 n_negative", p450["n_negative"], 1253),
]
for name, got, want in checks:
    if got != want:
        raise SystemExit(f"[watcher][ERROR] {name}: got {got}, want {want}")
print("[watcher] final eval summary path:", summary_path)
print(json.dumps(summary, indent=2))
(root / "EXP003_FINAL_EVAL_POINTER.txt").write_text(str(summary_path) + "\n", encoding="utf-8")

block = f"""

## {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}：EXP003 自动收尾完成

EXP003 已完成 35 epoch 训练后的最终测试。watcher 已核对 full test 与 389 P450 子 test 的 JSON 结果、样本数和正负样本数。

| 测试集 | 样本数 | 正样本 | 负样本 | 测试 AUC | 测试 AUPR |
|---|---:|---:|---:|---:|---:|
| 全量 test | {full['n_samples']} | {full['n_positive']} | {full['n_negative']} | {full['test_auc_roc']} | {full['test_aupr']} |
| 389 P450 子 test | {p450['n_samples']} | {p450['n_positive']} | {p450['n_negative']} | {p450['test_auc_roc']} | {p450['test_aupr']} |

最终 summary：

`{summary_path}`

最终 checkpoint：

`{summary['checkpoint']}`

watcher 在写入本记录后执行服务器关机。
"""
with session_log.open("a", encoding="utf-8") as f:
    f.write(block)
PY

echo "[watcher] final evaluation succeeded"
echo "[watcher] shutting down server in 2 minutes: $(date)"
sync
shutdown -h +2
