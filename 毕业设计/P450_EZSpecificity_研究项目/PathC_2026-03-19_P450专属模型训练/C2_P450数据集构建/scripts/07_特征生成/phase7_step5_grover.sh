#!/bin/bash
# Phase 7 Step 5: GROVER embeddings on GPU 1 (3 sequential substeps)
set -e
export PATH=/opt/conda/envs/ezspec/bin:/usr/bin:/bin
export PYTHONPATH=/root/rivermind-data/EZSpecificity/src:/root/rivermind-data/EZSpecificity/src/other_softwares/grover_software
FEATURES=/root/rivermind-data/EZSpecificity/PathC/P450
GROVER_DIR=/root/rivermind-data/EZSpecificity/src/other_softwares/grover_software
CKPT=/root/rivermind-data/EZSpecificity/data/pretrain_model/grover_large.pt

# Prepare GROVER CSV
python -c "
import pandas as pd
df = pd.read_csv('/root/rivermind-data/EZSpecificity/PathC/dataset/csv/Substrates.csv').dropna(subset=['Substrate_SMILES'])
pd.DataFrame({'Substrate_SMILES': df['Substrate_SMILES'].values}).to_csv('$FEATURES/grover_substrates.csv', index=False)
print(f'GROVER CSV: {len(df)} rows')
"

cd $GROVER_DIR
echo "Step 5a: save_features"
python scripts/save_features.py --data_path $FEATURES/grover_substrates.csv --save_path $FEATURES/grover_substrates.npz --features_generator fgtasklabel --restart

echo "Step 5b: build_vocab"
mkdir -p $FEATURES/grover_vocab
python scripts/build_vocab.py --data_path $FEATURES/grover_substrates.csv --vocab_save_folder $FEATURES/grover_vocab --dataset_name test

echo "Step 5c: fingerprint"
CUDA_VISIBLE_DEVICES=1 python main.py fingerprint --data_path $FEATURES/grover_substrates.csv --features_path $FEATURES/grover_substrates.npz --checkpoint_path $CKPT --fingerprint_source both --output $FEATURES/grover_fingerprint.npz --save_lmdb_path $FEATURES/grover_fingerprint.lmdb

echo "GROVER DONE"
