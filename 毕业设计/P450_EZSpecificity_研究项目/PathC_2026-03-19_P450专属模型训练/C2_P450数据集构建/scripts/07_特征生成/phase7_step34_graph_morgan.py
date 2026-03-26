"""Phase 7 Step 3+4: Graph features + Morgan fingerprint on CPU."""
import sys, lmdb, pickle, numpy as np, pandas as pd
sys.path.insert(0, '/root/rivermind-data/EZSpecificity/src')
from Datasets.create_features import get_reaction_feature_single
from rdkit.Chem import AllChem, DataStructs

substrate_csv = '/root/rivermind-data/EZSpecificity/PathC/P450/Substrates.csv'
features_dir = '/root/rivermind-data/EZSpecificity/PathC/P450'

df = pd.read_csv(substrate_csv).dropna(subset=['Substrate_SMILES'])

# Graph features
print('Building graph features...', flush=True)
params = {}
for s in df['Substrate_SMILES']:
    if s not in params and len(s) < 275:
        params[s] = (len(params), (f'{s}>>{s}', s, s, True))
items = [v[1] for v in params.values()]
db = lmdb.open(f'{features_dir}/reaction_features.lmdb', map_size=600*(1024**3), create=True, subdir=False, readonly=False)
reactions = []; substrates = []
for full_reaction, substrate, tag, data in map(get_reaction_feature_single, items):
    if data is not None or tag is None:
        with db.begin(write=True, buffers=True) as txn:
            txn.put(str(len(reactions)).encode(), pickle.dumps(data))
        reactions.append(full_reaction); substrates.append(substrate)
    else:
        print(f'Skipped {full_reaction} {tag}', flush=True)
db.close()
print(f'Graph features: {len(reactions)}', flush=True)

# Morgan fingerprint
print('Building morgan fingerprint...', flush=True)
fps = []
for s in df['Substrate_SMILES']:
    mol = AllChem.MolFromSmiles(s)
    if mol is None:
        fps.append(np.zeros(1024, dtype=np.int8)); continue
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)
    arr = np.zeros(1024, dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps.append(arr)
np.save(f'{features_dir}/morgan_fingerprint.npy', np.array(fps))
print(f'Morgan: {len(fps)}. DONE.', flush=True)
