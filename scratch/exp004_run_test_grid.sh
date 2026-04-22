#!/bin/bash
# EXP004 diagnostic grid: 4-way comparison to disambiguate edge-mode and weight effects.
# Run this AFTER the main run_train.sh to get the full sensitivity picture.

set -e

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${PATH}

EXP=~/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP004_paper_baseline_unified
CACHE=~/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_unified_paperfilter/random
PAPER_CKPT=${EXP}/results/checkpoints/paper_best-checkpoint.ckpt
OURS_CKPT=${EXP}/results/checkpoints/ours_EXP001_ep43_auc0.9250.ckpt

SCRIPT=${EXP}/scripts/main_training_pt.py
CFG=${EXP}/configs/config.yml

run_eval () {
  CKPT=$1; MODE=$2; TAG=$3
  OUT=${EXP}/results/test_eval_${TAG}.json
  echo ""
  echo "============================================================"
  echo "Run: ${TAG}"
  echo "  ckpt: $(basename ${CKPT})"
  echo "  edge: ${MODE}"
  echo "  out : ${OUT}"
  echo "============================================================"
  python ${SCRIPT} \
    --config ${CFG} \
    --cache-dir ${CACHE} \
    --test-only \
    --checkpoint ${CKPT} \
    --edge-mode ${MODE} \
    --batch-size 88 \
    --num-workers 6 \
    --run-name ${TAG} \
    --output-json ${OUT}
}

# 1. Paper ckpt + legacy_bug (main — matches paper training preprocessing)
run_eval ${PAPER_CKPT} legacy_bug paper_legacy_filtered

# 2. Paper ckpt + fixed (edge sensitivity for paper weights)
run_eval ${PAPER_CKPT} fixed paper_fixed_filtered

# 3. Ours ckpt + legacy_bug (symmetric control — our weights trained fixed, should prefer fixed)
run_eval ${OURS_CKPT}  legacy_bug ours_legacy_filtered

# 4. Ours ckpt + fixed (our native mode on filtered subset; compare to 0.9320 full)
run_eval ${OURS_CKPT}  fixed      ours_fixed_filtered

echo ""
echo "============================================================"
echo "ALL 4 DIAGNOSTIC RUNS COMPLETE"
echo "============================================================"
ls -la ${EXP}/results/test_eval_*.json
