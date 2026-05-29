#!/usr/bin/env bash
set -euo pipefail

ROOT=/root/autodl-tmp/EZSpecificity/PathD/P450
EXP=${ROOT}/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay
RESULT_DIR=${EXP}/results/00_EXP003_FULL_AF_HEME_FE_20260528_055549
FULL_CACHE=${ROOT}/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/full_exp003_alphafill_heme_fe_v1
P450389_CACHE=${ROOT}/data/q01_fe_embedding_patch/eval_subsets/p450389_test_subset_20260528_033240/exp003_p450389_test_cache
CONFIG=${EXP}/configs/config_exp003_full_alphafill_heme_fe.yml
SCRIPT=${EXP}/scripts/main_training_pt_fe_overlay.py
OUT=${ROOT}/results/q01_fe_embedding_patch/exp003_eval_$(date +%Y%m%d_%H%M%S)

mkdir -p "${OUT}"

CKPT=$(/root/miniconda3/bin/python - <<'PY'
from pathlib import Path
import csv
import re

result_dir = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/results/00_EXP003_FULL_AF_HEME_FE_20260528_055549")
metrics = result_dir / "metrics.csv"
ckpt_dir = result_dir / "checkpoints"

if not metrics.exists():
    raise SystemExit("metrics.csv missing")
rows = list(csv.DictReader(metrics.open()))
epochs = sorted({int(r["epoch"]) for r in rows if r.get("epoch", "").isdigit()})
if not epochs or max(epochs) < 34:
    raise SystemExit(f"35-epoch training is not complete; max metric epoch={max(epochs) if epochs else 'NA'}")

items = []
for p in ckpt_dir.glob("*.ckpt"):
    if p.name == "last.ckpt":
        continue
    m_auc = re.search(r"auc([0-9.]+)\.ckpt$", p.name)
    m_ep = re.search(r"-ep(\d+)-", p.name)
    if m_auc and m_ep:
        items.append((float(m_auc.group(1)), int(m_ep.group(1)), p))
if not items:
    raise SystemExit("No scored checkpoint found")
auc, epoch, path = max(items, key=lambda x: (x[0], x[1]))
print(path)
PY
)

echo "[EXP003 eval] checkpoint: ${CKPT}"
echo "[EXP003 eval] output: ${OUT}"

CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0} /root/miniconda3/bin/python "${SCRIPT}" \
  --config "${CONFIG}" \
  --cache-dir "${FULL_CACHE}" \
  --edge-mode legacy_bug \
  --num-workers 8 \
  --batch-size 40 \
  --results-dir "${OUT}/tmp_exp003_full" \
  --test-only \
  --checkpoint "${CKPT}" \
  --output-json "${OUT}/EXP003_full_test.json"

CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-0} /root/miniconda3/bin/python "${SCRIPT}" \
  --config "${CONFIG}" \
  --cache-dir "${P450389_CACHE}" \
  --edge-mode legacy_bug \
  --num-workers 4 \
  --batch-size 40 \
  --results-dir "${OUT}/tmp_exp003_p450389" \
  --test-only \
  --checkpoint "${CKPT}" \
  --output-json "${OUT}/EXP003_p450389_test.json"

/root/miniconda3/bin/python - <<PY
import json
from pathlib import Path

out = Path("${OUT}")
full = json.loads((out / "EXP003_full_test.json").read_text())
p450 = json.loads((out / "EXP003_p450389_test.json").read_text())
summary = {
    "experiment": "Q1_EXP003_full_p450_fe_heme_overlay",
    "checkpoint": str(Path("${CKPT}")),
    "full_test": {
        "n_samples": full["n_samples"],
        "n_positive": full["n_positive"],
        "n_negative": full["n_negative"],
        "test_auc_roc": full["test_auc_roc"],
        "test_aupr": full["test_aupr"],
    },
    "p450389_test": {
        "n_samples": p450["n_samples"],
        "n_positive": p450["n_positive"],
        "n_negative": p450["n_negative"],
        "test_auc_roc": p450["test_auc_roc"],
        "test_aupr": p450["test_aupr"],
    },
}
(out / "EXP003_eval_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
PY

echo "[EXP003 eval] done: ${OUT}/EXP003_eval_summary.json"
