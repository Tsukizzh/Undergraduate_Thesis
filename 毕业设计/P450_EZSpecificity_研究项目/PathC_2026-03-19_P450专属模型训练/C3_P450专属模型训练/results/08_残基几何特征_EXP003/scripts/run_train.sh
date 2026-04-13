#!/bin/bash
# EXP003: Residue-level geometric features (phi/psi/chi1)
# Based on EXP002b hyperparams + EXP002a Fe/HEM + NEW angle features (dim=37)

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${PATH}

EXP=~/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP003_residue_geometry
CACHE=~/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_geom

python ${EXP}/scripts/main_training_pt.py \
  --config ${EXP}/configs/config.yml \
  --cache-dir ${CACHE} \
  --edge-mode fixed \
  --batch-size 88 \
  --max-epochs 200 \
  --devices 4 \
  --num-workers 6 \
   \
  --run-name EXP003_residue_geometry \
  --shutdown
