#!/bin/bash
# EXP004 primary run: paper ckpt + legacy_bug edges + paperfilter test cache.
# Main leakage-controlled external P450 baseline.

set -e

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${PATH}

EXP=~/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP004_paper_baseline_unified
CACHE=~/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_unified_paperfilter/random
CKPT=${EXP}/results/checkpoints/paper_best-checkpoint.ckpt
OUT=${EXP}/results/test_eval_paper_legacy_filtered.json

python ${EXP}/scripts/main_training_pt.py \
  --config ${EXP}/configs/config.yml \
  --cache-dir ${CACHE} \
  --test-only \
  --checkpoint ${CKPT} \
  --edge-mode legacy_bug \
  --batch-size 88 \
  --num-workers 6 \
  --run-name EXP004_paper_legacy_filtered \
  --output-json ${OUT}
