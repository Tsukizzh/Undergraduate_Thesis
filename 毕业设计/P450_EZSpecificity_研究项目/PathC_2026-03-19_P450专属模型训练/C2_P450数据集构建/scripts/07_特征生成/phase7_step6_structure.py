"""Phase 7 Step 6: Structure features (serial, CPU)."""
import sys
sys.path.insert(0, '/root/rivermind-data/EZSpecificity/src')
from easydict import EasyDict as edict
import pandas as pd
from Datasets.Structure.structure import str_process, seq_process

cfg = edict({
    "data_df_path": "/root/rivermind-data/EZSpecificity/PathC/P450/data.csv",
    "ligand_df_path": "/root/rivermind-data/EZSpecificity/PathC/P450/Substrates.csv",
    "data": {
        "structure_processed_path": "/root/rivermind-data/EZSpecificity/PathC/P450/structure/structure_features.lmdb",
        "sequence_processed_path": "/root/rivermind-data/EZSpecificity/PathC/P450/structure/sequence_features.lmdb",
        "pdb_dir": "/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/pocket",
        "ligand_dir": "/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/ligand",
    }
})

print("str_process...", flush=True)
df = pd.read_csv(cfg.data_df_path)
str_process(cfg, df)
print("seq_process...", flush=True)
df = pd.read_csv(cfg.ligand_df_path)
seq_process(cfg, df)
print("DONE", flush=True)
