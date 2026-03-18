#!/bin/bash
cd /root/rivermind-data/EZSpecificity/src
export PYTHONPATH=/root/rivermind-data/EZSpecificity/src:$PYTHONPATH
ulimit -n 65536

echo '=== Starting GPU0: legacy_bug ==='
CUDA_VISIBLE_DEVICES=0 nohup /opt/conda/bin/python ../scripts/10_Step10_pt训练管线/main_training_pt.py   --config ../scripts/server_config.yml   --cache-dir ../data/10_Step10_pt训练/ezspec_pt_v1   --edge-mode legacy_bug   --devices 1 --num-workers 4 --batch-size 56   --max-epochs 50 > /root/train_legacy.log 2>&1 &
echo "GPU0 PID=$!"

echo '=== Starting GPU1: fixed ==='
CUDA_VISIBLE_DEVICES=1 nohup /opt/conda/bin/python ../scripts/10_Step10_pt训练管线/main_training_pt.py   --config ../scripts/server_config.yml   --cache-dir ../data/10_Step10_pt训练/ezspec_pt_v1   --edge-mode fixed   --devices 1 --num-workers 4 --batch-size 56   --max-epochs 50 > /root/train_fixed.log 2>&1 &
echo "GPU1 PID=$!"
