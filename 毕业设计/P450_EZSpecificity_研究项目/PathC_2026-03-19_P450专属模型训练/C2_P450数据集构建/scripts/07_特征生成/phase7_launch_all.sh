#!/bin/bash
# Phase 7: Launch all feature generation tasks in parallel
# GPU0: ESM | GPU1: GROVER | CPU: pocket extraction + graph + morgan
set -e

export ROOT=/root/rivermind-data/EZSpecificity
export PATHC=$ROOT/PathC
export P450=$PATHC/P450
export SRC=$ROOT/src
export GROVER_DIR=$SRC/other_softwares/grover_software
export PYTHONPATH=$SRC:$GROVER_DIR:$PYTHONPATH
export EZPY=/opt/conda/envs/ezspec/bin/python
export UNIPY=/opt/conda/envs/unidock/bin/python
LOGS=$PATHC/logs

mkdir -p $LOGS $P450/structure/str_tmp_data/{pocket,raw_ligand,ligand} $P450/grover_vocab

echo "[$(date)] ===== Phase 7 Launch ====="

# === 1. Pocket extraction (CPU, 10 workers, unidock env for rdkit) ===
echo "[$(date)] [1] Pocket extraction..."
$UNIPY $PATHC/scripts/phase7_step1_pocket.py \
  > $LOGS/phase7_pocket.log 2>&1 &
PID_POCKET=$!

# === 2. ESM embeddings (GPU 0, ezspec env) ===
echo "[$(date)] [2] ESM on GPU 0..."
CUDA_VISIBLE_DEVICES=0 $EZPY $PATHC/scripts/phase7_step2_esm.py \
  > $LOGS/phase7_esm.log 2>&1 &
PID_ESM=$!

# === 3+4. Graph features + Morgan (CPU, ezspec env) ===
echo "[$(date)] [3+4] Graph + Morgan..."
$EZPY $PATHC/scripts/phase7_step34_graph_morgan.py \
  > $LOGS/phase7_graph_morgan.log 2>&1 &
PID_GRAPH=$!

# === 5. GROVER (GPU 1, 3 steps) ===
echo "[$(date)] [5] GROVER on GPU 1..."
bash $PATHC/scripts/phase7_step5_grover.sh \
  > $LOGS/phase7_grover.log 2>&1 &
PID_GROVER=$!

echo "[$(date)] PIDs: pocket=$PID_POCKET esm=$PID_ESM graph=$PID_GRAPH grover=$PID_GROVER"
echo "[$(date)] Waiting for all tasks..."

wait $PID_POCKET; echo "[$(date)] [1] Pocket DONE (exit=$?)"
wait $PID_ESM; echo "[$(date)] [2] ESM DONE (exit=$?)"
wait $PID_GRAPH; echo "[$(date)] [3+4] Graph+Morgan DONE (exit=$?)"
wait $PID_GROVER; echo "[$(date)] [5] GROVER DONE (exit=$?)"

echo "[$(date)] ===== All parallel tasks complete ====="
echo "[$(date)] Starting Step 1b: Ligand alignment..."

# === 1b. Ligand alignment (CPU, depends on pocket being done) ===
$UNIPY $PATHC/scripts/phase7_step1b_align.py \
  > $LOGS/phase7_align.log 2>&1
echo "[$(date)] [1b] Alignment DONE"

# === 6. Structure features (CPU, serial, depends on 1b) ===
echo "[$(date)] [6] Structure features..."
$EZPY $PATHC/scripts/phase7_step6_structure.py \
  > $LOGS/phase7_structure.log 2>&1
echo "[$(date)] [6] Structure DONE"

# === 7. High quality ID ===
echo "[$(date)] [7] High quality ID..."
$EZPY $PATHC/scripts/phase7_step7_hqid.py \
  > $LOGS/phase7_hqid.log 2>&1
echo "[$(date)] [7] HQID DONE"

echo "[$(date)] ===== Phase 7 COMPLETE ====="
