#!/bin/bash
# EXP005: dualgraph 2+ (residue backfill + g_res bypass)
# Training recipe strictly matches EXP001_allfix_unified (Test AUC 0.9320).
# Only 3 variables differ: EXP path / cache path / run name.
# Model class differs via main_training_pt.py imports (SS -> SSDualgraph).
#
# NOTE: --preload was tested once and OOM'd due to Python refcount COW on
# dict-of-tensors preload + /dev/shm queue pressure with 4×6=24 workers.
# Baseline allfix series has NEVER used --preload, so we don't either.

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${PATH}

EXP=~/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP005_dualgraph_2plus_allfix_unified
CACHE=~/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_dualgraph_allfix_unified/random

python ${EXP}/scripts/main_training_pt.py \
  --config ${EXP}/configs/config.yml \
  --cache-dir ${CACHE} \
  --edge-mode fixed \
  --batch-size 88 \
  --max-epochs 200 \
  --devices 4 \
  --num-workers 6 \
  --run-name EXP005_dualgraph_2plus_allfix_unified \
  --shutdown
