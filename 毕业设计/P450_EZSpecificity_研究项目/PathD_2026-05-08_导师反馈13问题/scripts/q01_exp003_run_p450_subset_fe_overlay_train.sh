#!/usr/bin/env bash
set -euo pipefail

cd /root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay

CACHE_DIR=/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1
CONFIG=configs/config_exp003_p450_subset_fe_overlay.yml
STAMP=$(date +%Y%m%d_%H%M%S)
DEVICES=${DEVICES:-1}
BATCH_SIZE=${BATCH_SIZE:-40}
NUM_WORKERS=${NUM_WORKERS:-8}
PREFETCH_FACTOR=${PREFETCH_FACTOR:-1}
MAX_EPOCHS=${MAX_EPOCHS:-200}
RUN_NAME=${RUN_NAME:-Q1_EXP003_p450_subset_fe_overlay_d${DEVICES}_b${BATCH_SIZE}_w${NUM_WORKERS}_${STAMP}}
RESULTS_DIR=${RESULTS_DIR:-results/00_EXP003_P450_SUBSET_FE_OVERLAY_${STAMP}}
LIMIT_TRAIN_B=${LIMIT_TRAIN_BATCHES:-}
LIMIT_VAL_B=${LIMIT_VAL_BATCHES:-}
EXTRA_ARGS=${EXTRA_ARGS:-}

args=(
  --config "$CONFIG"
  --cache-dir "$CACHE_DIR"
  --edge-mode legacy_bug
  --num-workers "$NUM_WORKERS"
  --prefetch-factor "$PREFETCH_FACTOR"
  --train-in-order false
  --batch-size "$BATCH_SIZE"
  --max-epochs "$MAX_EPOCHS"
  --devices "$DEVICES"
  --run-name "$RUN_NAME"
  --results-dir "$RESULTS_DIR"
)

if [[ -n "$LIMIT_TRAIN_B" ]]; then
  args+=(--limit-train-batches "$LIMIT_TRAIN_B")
fi
if [[ -n "$LIMIT_VAL_B" ]]; then
  args+=(--limit-val-batches "$LIMIT_VAL_B")
fi

# shellcheck disable=SC2206
extra=( $EXTRA_ARGS )

OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} \
MKL_NUM_THREADS=${MKL_NUM_THREADS:-1} \
OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} \
/root/miniconda3/bin/python scripts/main_training_pt_fe_overlay.py "${args[@]}" "${extra[@]}"
