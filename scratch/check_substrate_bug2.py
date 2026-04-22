import pandas as pd, lmdb, pickle, sys
from rdkit import Chem

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
df = pd.read_csv(f'{base}/Substrates.csv')

# Step 1: identify which Substrate Index rows were ACCEPTED into params (len<275)
# and iterate in df order. Since no dupes, each row either accepted or skipped (len>=275).
accepted_rows = []
skipped_len = []
for i, s in enumerate(df['Substrate_SMILES']):
    if pd.isna(s): continue
    if len(s) < 275:
        accepted_rows.append(i)
    else:
        skipped_len.append(i)
print('accepted (len<275):', len(accepted_rows))
print('skipped (len>=275):', skipped_len)
print('First accepted rows:', accepted_rows[:5])

# Step 2: Among accepted, which succeeded get_reaction_feature_single?
# reaction_features.lmdb has 2119 keys. Expected 2124 from len filter.
# 2124-2119 = 5 rows failed parsing.
# These 5 failure points cause all subsequent (Substrate Index -> LMDB key) offsets.

# Probe: grab LMDB key N, pickle-decode, reconstruct SMILES from the dict, find which accepted_rows[i] it corresponds to
env = lmdb.open(f'{base}/features/reaction_features.lmdb', subdir=False, readonly=True, lock=False, readahead=False)
def get(key):
    with env.begin() as txn:
        return pickle.loads(txn.get(str(key).encode()))

sample = get(0)
print('\nkey 0 dict keys:', list(sample.keys())[:15] if isinstance(sample, dict) else 'not dict')
# Check if there's a SMILES field
for k in ['smiles', 'reaction', 'mol', 'smile']:
    if k in sample:
        print(f'  {k}:', sample[k])

# If no SMILES, try comparing number of heavy atoms or similar
# The reaction data object typically has 'x' or 'num_atom' field
print('\nKey shapes for 0, 100, 1000, 2118:')
for k in [0, 100, 1000, 2118]:
    d = get(k)
    if isinstance(d, dict):
        print(f'  key {k}:', {kk: (vv.shape if hasattr(vv,'shape') else type(vv).__name__) for kk,vv in d.items() if kk in ['x','edge_index','num_atom']})
