"""Phase 7 Step 2: ESM embeddings on GPU 0."""
import sys, lmdb, pandas as pd
sys.path.insert(0, '/root/rivermind-data/EZSpecificity/src')
from Datasets.create_features import generate_esm_embedding, convert_protein_sequence_to_number

enzyme_csv = '/root/rivermind-data/EZSpecificity/PathC/P450/Enzymes.csv'
out_lmdb = '/root/rivermind-data/EZSpecificity/PathC/P450/enzyme_features.lmdb'

df = pd.read_csv(enzyme_csv)
env = lmdb.open(out_lmdb, map_size=600*(1024**3), create=True, subdir=False, readonly=False)
data = []; uniprot_dict = {}
for idx, seq in enumerate(df['Protein sequence']):
    if len(seq) > 1000:
        print(f'skip {idx}: too long ({len(seq)})', flush=True); continue
    try:
        convert_protein_sequence_to_number(seq)
    except:
        print(f'skip {idx}: invalid aa', flush=True); continue
    uniprot_dict[str(idx)] = (len(uniprot_dict), 1)
    data.append((str(idx), seq))
print(f'ESM inputs: {len(data)}', flush=True)
generate_esm_embedding(env, data, uniprot_dict)
env.close()
print('ESM DONE', flush=True)
