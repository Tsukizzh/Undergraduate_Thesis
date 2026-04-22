import pandas as pd, lmdb, pickle, numpy as np

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
df = pd.read_csv(f'{base}/Substrates.csv')

# Find [H] or very short SMILES that might have failed
for i, s in enumerate(df['Substrate_SMILES']):
    if '[H]' in str(s) and len(str(s)) < 10:
        print(f'Substrate Index {i}: {s!r}')

# Morgan approach: load Morgan, it's definitely aligned to Substrate Index by row
morgan = np.load(f'{base}/features/morgan_fingerprint.npy')

# Load GROVER for each key, see if key N matches Substrate Index N or N+1
env = lmdb.open(f'{base}/features/grover_fingerprint.lmdb', subdir=False, readonly=True, lock=False, readahead=False)

# Trick: GROVER also stores a SMILES tag? Or we can use its shape/content as fingerprint
def get_grover(k):
    with env.begin() as txn:
        raw = txn.get(str(k).encode())
    if raw is None: return None
    return pickle.loads(raw)

g0 = get_grover(0)
print('\nGROVER key 0 type:', type(g0).__name__)
if isinstance(g0, dict):
    print('  keys:', list(g0.keys()))
    for k in ['smiles', 'smile', 'SMILES']:
        if k in g0: print(f'  {k}: {g0[k]}')

# We can't directly read SMILES from GROVER, so use GROVER total_embedding similarity
# Hypothesis: if no bug, GROVER key N embedding matches Substrate Index N
# Strategy: GROVER reads substrates sequentially. If one fails, all subsequent shift.
# Find the failure substrate by loading all GROVER and comparing count
print('\nGROVER n_keys:', sum(1 for _,_ in env.begin().cursor()))
print('Substrate rows:', len(df))
print('Diff:', len(df) - sum(1 for _,_ in env.begin().cursor()))

# Quick test: GROVER key 2123 (last). Is it Substrate Index 2123 or 2124?
# Load total_embedding and check matches via fingerprint.npz which should have SMILES metadata
try:
    npz = np.load(f'{base}/features/grover_fingerprint.npz', allow_pickle=True)
    print('\ngrover npz files:', npz.files)
    for k in npz.files:
        v = npz[k]
        print(f'  {k}: shape={v.shape if hasattr(v,"shape") else type(v).__name__}')
except Exception as e:
    print('npz load fail:', e)

# Check grover_substrates.csv (the one fed into GROVER)
try:
    gdf = pd.read_csv(f'{base}/features/grover_substrates.csv')
    print('\ngrover_substrates.csv rows:', len(gdf))
    # Find any SMILES with [H]
    for i, s in enumerate(gdf['Substrate_SMILES']):
        if '[H]' in str(s) and len(str(s)) < 10:
            print(f'  gCSV row {i}: {s!r}')
except Exception as e:
    print('gCSV load fail:', e)
