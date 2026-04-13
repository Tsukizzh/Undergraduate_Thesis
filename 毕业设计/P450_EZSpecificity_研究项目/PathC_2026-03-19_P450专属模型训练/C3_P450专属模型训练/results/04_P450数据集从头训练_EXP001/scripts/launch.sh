#!/bin/bash
# EXP001: P450 random_split baseline - 4x RTX4090 DDP
set -e
cd /root/rivermind-data/EZSpecificity
ulimit -n 65536

EXP_DIR=PathC/P450/experiments/EXP001_random_baseline
CACHE_DIR=PathC/P450/data/pt_cache/random

export PYTHONPATH=/root/rivermind-data/EZSpecificity/src:

/opt/conda/envs/ezspec/bin/torchrun --nproc_per_node=4 /scripts/main_training_pt.py   --config /config.yml   --cache-dir    --edge-mode fixed   --devices 4   --batch-size 56   --num-workers 6   --max-epochs 50   --preload   --run-name EXP001_random_baseline   2>&1 | tee /train.log
