#!/bin/bash
# EXP002b: lr tuning (based on EXP002a, 2-GPU server)

ulimit -n 65536 2>/dev/null
export OMP_NUM_THREADS=1
export PATH=/opt/conda/envs/ezspec/bin:$PATH

EXP=~/rivermind-data/EZSpecificity/PathC/P450/experiments/EXP002b_lr_tuning
CACHE=~/rivermind-data/EZSpecificity/PathC/P450/data/pt_cache_heme/random

python $EXP/scripts/main_training_pt.py \
  --config $EXP/configs/config.yml \
  --cache-dir $CACHE \
  --edge-mode fixed \
  --batch-size 56 \
  --max-epochs 200 \
  --devices 2 \
  --num-workers 5 \
  --preload \
  --run-name lr_tuning \
  --shutdown
