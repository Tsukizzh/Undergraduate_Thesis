#!/bin/bash
# EXP003_fixed: EXP003 rerun on FIXED pt_cache (LMDB key alignment bug fixed)
#
# Purpose: establish the true baseline AUC after repairing the enzyme_features.lmdb
# key mapping bug. Training config is IDENTICAL to EXP003 (same hyperparams, same
# edge-mode, same max-epochs) so the ONLY variable is the training data:
#   - EXP003:       pt_cache_geom        (old flatbin with compressed keys 0..1576)
#   - EXP003_fixed: pt_cache_geom_fixed  (new flatbin with original CSV row keys)
#
# After training, compare test AUC with original EXP003 (0.7914) to measure the
# true impact of the data fix.

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${PATH}

EXP=~/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP003_fixed
CACHE=~/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_geom_fixed

python ${EXP}/scripts/main_training_pt.py \
  --config ${EXP}/configs/config.yml \
  --cache-dir ${CACHE} \
  --edge-mode fixed \
  --batch-size 88 \
  --max-epochs 200 \
  --devices 4 \
  --num-workers 6 \
  --run-name EXP003_fixed \
  --shutdown
