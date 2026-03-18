#!/bin/bash
cd /root/rivermind-data/EZSpecificity/src
export PYTHONPATH=/root/rivermind-data/EZSpecificity/src:$PYTHONPATH
ulimit -n 65536
nohup /opt/conda/bin/python ../scripts/10_Step10_pt训练管线/main_training_pt.py   --config ../scripts/server_config.yml   --cache-dir ../data/10_Step10_pt训练/ezspec_pt_v1   --edge-mode legacy_bug   --devices 2 --num-workers 6 --batch-size 56   --max-epochs 50 --resume last --shutdown > /root/train.log 2>&1 &
echo PID=$!
