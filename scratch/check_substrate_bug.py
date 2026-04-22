import pandas as pd, lmdb, pickle, sys
from rdkit import Chem

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
df = pd.read_csv(f'{base}/Substrates.csv')
print('Substrates.csv rows:', len(df))
print('Unique SMILES:', df['Substrate_SMILES'].nunique())
print('NaN SMILES:', df['Substrate_SMILES'].isna().sum())

# Replay script logic: unique + len<275
seen = set(); order = []
for i, s in enumerate(df['Substrate_SMILES']):
    if pd.isna(s): continue
    if s not in seen and len(s) < 275:
        seen.add(s); order.append((i, s))
print('Expected unique usable (replay):', len(order))
print('First 5 row_idx where unique SMILES first appears:', [r for r,_ in order[:5]])
print('First 5 original Substrate Index gap? any skip?')
first_rows = [r for r,_ in order[:20]]
print('First 20 row_idx:', first_rows)

# Check reaction LMDB keys
env = lmdb.open(f'{base}/features/reaction_features.lmdb', subdir=False, readonly=True, lock=False, readahead=False)
with env.begin() as txn:
    keys = [k.decode() for k,_ in txn.cursor()]
int_keys = sorted(int(k) for k in keys)
print('reaction_features.lmdb keys:', len(int_keys), 'range', int_keys[0], '-', int_keys[-1])
print('contiguous 0..N-1?', int_keys == list(range(len(int_keys))))

# Key meaning test: grab key "5" and see if its SMILES matches Substrate Index 5 OR the 6th unique SMILES
with env.begin() as txn:
    data5 = pickle.loads(txn.get(b'5'))
print('\nkey "5" type:', type(data5).__name__)
# data should be a graph Data obj with num_nodes attribute
try:
    print('key "5" num_nodes:', data5.num_nodes, 'num_edges:', data5.num_edges)
except Exception as e:
    print('inspect failed:', e)

# Does key "5" match Substrate Index 5's SMILES or the 6th unique SMILES?
sub5_smiles = df.iloc[5]['Substrate_SMILES']
uniq6_row, uniq6_smiles = order[5]  # 6th unique
print(f'\nSubstrate Index 5 SMILES: {sub5_smiles}')
print(f'  row_idx 5 atom count: {Chem.MolFromSmiles(sub5_smiles).GetNumAtoms() if Chem.MolFromSmiles(sub5_smiles) else "parse fail"}')
print(f'6th unique SMILES (row_idx {uniq6_row}): {uniq6_smiles}')
print(f'  6th unique atom count: {Chem.MolFromSmiles(uniq6_smiles).GetNumAtoms() if Chem.MolFromSmiles(uniq6_smiles) else "parse fail"}')
