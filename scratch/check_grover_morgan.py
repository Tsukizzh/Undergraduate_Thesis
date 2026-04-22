import pandas as pd, lmdb, pickle, numpy as np

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
df = pd.read_csv(f'{base}/Substrates.csv')
print('Substrates.csv rows:', len(df))

# GROVER
env = lmdb.open(f'{base}/features/grover_fingerprint.lmdb', subdir=False, readonly=True, lock=False, readahead=False)
with env.begin() as txn:
    gkeys = sorted([k.decode() for k,_ in txn.cursor()])
print('grover keys:', len(gkeys), 'first 5:', gkeys[:5], 'last 5:', gkeys[-5:])
# Are keys Substrate Index or SMILES or counter?
try:
    int_keys = sorted(int(k) for k in gkeys)
    print('grover int keys range:', int_keys[0], '-', int_keys[-1])
    print('grover contiguous?', int_keys == list(range(len(int_keys))))
    print('missing int keys:', set(range(int_keys[-1]+1)) - set(int_keys))
except ValueError:
    print('grover keys are not ints. Sample:', gkeys[:3])

# Morgan
m = np.load(f'{base}/features/morgan_fingerprint.npy')
print('\nmorgan shape:', m.shape)
