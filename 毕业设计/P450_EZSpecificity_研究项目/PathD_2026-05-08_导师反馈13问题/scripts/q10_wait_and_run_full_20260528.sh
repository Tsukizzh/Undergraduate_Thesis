#!/usr/bin/env bash
set -euo pipefail

Q10_DIR="/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
Q10_EXP="/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
PY="/root/miniconda3/bin/python"
CF_PY="/root/autodl-tmp/envs/colabfold_q10/bin/python"
COLLECT="${Q10_EXP}/scripts/q10_collect_colabfold_manifest_20260528.py"
FULL="${Q10_EXP}/scripts/q10_run_full_after_colabfold_20260528.sh"
SUMMARY="${Q10_DIR}/structures/manifests/q10_colabfold_structure_summary_msa_m1r1_v1.json"
LOG_DIR="${Q10_DIR}/logs/q10_full_pipeline"
LOCK_FILE="${LOG_DIR}/q10_full_pipeline.lock"
POLL_SECONDS="${POLL_SECONDS:-300}"

mkdir -p "${LOG_DIR}"
echo "[Q10 wait] started at $(date)"

while true; do
  "${CF_PY}" "${COLLECT}" --q10-dir "${Q10_DIR}" --link-selected >/tmp/q10_collect_latest.json
  cat /tmp/q10_collect_latest.json

  read -r total complete pending < <("${PY}" -c "import json; d=json.load(open('${SUMMARY}',encoding='utf-8')); print(d.get('total'), d.get('complete'), d.get('pending_or_partial'))")
  echo "[Q10 wait] $(date) total=${total} complete=${complete} pending_or_partial=${pending}"

  if [[ "${total}" == "204" && "${complete}" == "204" ]]; then
    echo "[Q10 wait] ColabFold complete. Starting full Q10 HEM/Fe-aware scoring pipeline."
    exec flock -n "${LOCK_FILE}" "${FULL}"
  fi

  sleep "${POLL_SECONDS}"
done
