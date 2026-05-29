#!/usr/bin/env bash
set -euo pipefail

Q10_DIR="/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
Q10_EXP="/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528"
CF_PY="/root/autodl-tmp/envs/colabfold_q10/bin/python"
CF_BIN="/root/autodl-tmp/envs/colabfold_q10/bin/colabfold_batch"
COLLECT="${Q10_EXP}/scripts/q10_collect_colabfold_manifest_20260528.py"
PENDING="${Q10_DIR}/structures/manifests/q10_colabfold_pending_msa_m1r1_v1.fasta"
ACCEL_INPUT_DIR="${Q10_DIR}/structures/colabfold_inputs_accel"
ACCEL_OUT_DIR="${Q10_DIR}/structures/colabfold_batch"
LOG_DIR="${Q10_DIR}/logs/colabfold_accel"
SHARDS="${SHARDS:-6}"

mkdir -p "${ACCEL_INPUT_DIR}" "${ACCEL_OUT_DIR}" "${LOG_DIR}"

"${CF_PY}" "${COLLECT}" --q10-dir "${Q10_DIR}" --link-selected

if [[ ! -s "${PENDING}" ]]; then
  echo "[Q10 accel] no pending FASTA records."
  exit 0
fi

"${CF_PY}" - "${PENDING}" "${ACCEL_INPUT_DIR}" "${SHARDS}" <<'PY'
import sys
from pathlib import Path

pending = Path(sys.argv[1])
out_dir = Path(sys.argv[2])
shards = int(sys.argv[3])

records = []
header = None
seq_lines = []
for line in pending.read_text(encoding="utf-8").splitlines():
    if line.startswith(">"):
        if header is not None:
            records.append((header, "".join(seq_lines)))
        header = line[1:].strip()
        seq_lines = []
    else:
        seq_lines.append(line.strip())
if header is not None:
    records.append((header, "".join(seq_lines)))

for idx in range(shards):
    shard_records = records[idx::shards]
    path = out_dir / f"q10_pending_accel_shard_{idx:02d}.fasta"
    with path.open("w", encoding="utf-8", newline="\n") as f:
        for h, s in shard_records:
            f.write(f">{h}\n")
            for start in range(0, len(s), 80):
                f.write(s[start:start+80] + "\n")
    print(f"{path}\t{len(shard_records)}")
PY

for shard in $(seq -f "%02g" 0 $((SHARDS - 1))); do
  fasta="${ACCEL_INPUT_DIR}/q10_pending_accel_shard_${shard}.fasta"
  if [[ ! -s "${fasta}" ]]; then
    continue
  fi
  out="${ACCEL_OUT_DIR}/CF_ACCEL_${shard}_msa_m1r1_v1"
  log="${LOG_DIR}/CF_ACCEL_${shard}_msa_m1r1_v1.log"
  if pgrep -f "CF_ACCEL_${shard}_msa_m1r1_v1" >/dev/null; then
    echo "[Q10 accel] shard ${shard} already running."
    continue
  fi
  gpu=$((10#${shard} % 2))
  echo "[Q10 accel] launch shard=${shard} gpu=${gpu} fasta=${fasta} out=${out}"
  CUDA_VISIBLE_DEVICES="${gpu}" nohup "${CF_BIN}" \
    --msa-mode mmseqs2_uniref_env \
    --num-models 1 \
    --num-recycle 1 \
    --model-type alphafold2_ptm \
    --sort-queries-by length \
    "${fasta}" "${out}" > "${log}" 2>&1 &
done

echo "[Q10 accel] launched pending acceleration shards."
