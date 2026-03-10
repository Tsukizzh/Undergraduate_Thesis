#!/bin/bash
# Detached training launcher - survives parent process termination
cd "D:/EZSpecificity_Project/src"
LOG_DIR="D:/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/sessions/09_Step9_AllSplit训练"

echo "Starting training at $(date)" > "$LOG_DIR/training_detached.log"
echo "PID: $$" >> "$LOG_DIR/training_detached.log"

PYTHONUTF8=1 D:/anaconda3/envs/torch/python.exe -X utf8 -u \
  "../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/main_training.py" \
  >> "$LOG_DIR/training_detached.log" 2>&1

echo "Training finished with exit code $? at $(date)" >> "$LOG_DIR/training_detached.log"
