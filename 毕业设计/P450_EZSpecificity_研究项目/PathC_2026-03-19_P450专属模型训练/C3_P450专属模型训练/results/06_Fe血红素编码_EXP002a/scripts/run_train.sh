#!/bin/bash
# EXP002a training script — matches EXP001 setup exactly
# Only difference: --cache-dir points to pt_cache_heme

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PATH=/opt/conda/envs/ezspec/bin:$PATH

EXP=~/rivermind-data/EZSpecificity/PathC/P450/experiments/EXP002a_fe_heme_10A
CACHE=~/rivermind-data/EZSpecificity/PathC/P450/data/pt_cache_heme/random

python $EXP/scripts/main_training_pt.py \
  --config $EXP/configs/config.yml \
  --cache-dir $CACHE \
  --edge-mode fixed \
  --batch-size 56 \
  --max-epochs 200 \
  --devices 4 \
  --num-workers 6 \
  --preload \
  --run-name fe_heme_10A \
  --shutdown
