"""Try GetNumAtoms (includes dummy *) instead of GetNumHeavyAtoms."""
import pandas as pd, lmdb, pickle
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
orig = pd.read_csv(f'{base}/Substrates.csv')
gcsv = pd.read_csv(f'{base}/features/grover_substrates.csv')

# Quick sanity: what does GetNumAtoms vs GetNumHeavyAtoms give for orig[0]?
m = Chem.MolFromSmiles(orig.iloc[0]['Substrate_SMILES'])
print(f'orig[0] SMILES: {orig.iloc[0]["Substrate_SMILES"]}')
print(f'  GetNumAtoms: {m.GetNumAtoms()}')
print(f'  GetNumHeavyAtoms: {m.GetNumHeavyAtoms()}')
print(f'  explicit H: {sum(1 for a in m.GetAtoms() if a.GetSymbol()=="H")}')

def natoms(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
        return m.GetNumAtoms() if m is not None else None
    except Exception:
        return None

env = lmdb.open(f'{base}/features/grover_fingerprint.lmdb', subdir=False, readonly=True, lock=False, readahead=False)

print('\nExhaustive scan using GetNumAtoms():')
n_orig = n_shift = n_gcsv = total = 0
mismatches = []
with env.begin() as txn:
    for k_bytes, raw in txn.cursor():
        k = int(k_bytes.decode())
        d = pickle.loads(raw)
        lmdb_atoms = d['embedding'].shape[0]
        oN = natoms(orig.iloc[k]['Substrate_SMILES']) if k < len(orig) else None
        oN1 = natoms(orig.iloc[k+1]['Substrate_SMILES']) if k+1 < len(orig) else None
        gN = natoms(gcsv.iloc[k]['Substrate_SMILES']) if k < len(gcsv) else None
        total += 1
        if oN == lmdb_atoms: n_orig += 1
        if oN1 == lmdb_atoms: n_shift += 1
        if gN == lmdb_atoms: n_gcsv += 1
        if oN != lmdb_atoms and oN1 != lmdb_atoms:
            mismatches.append((k, lmdb_atoms, oN, oN1, gN, orig.iloc[k]['Substrate_SMILES'][:40]))

print(f'  total: {total}')
print(f'  lmdb == orig[N] GetNumAtoms:   {n_orig}/{total} ({100*n_orig/total:.2f}%)')
print(f'  lmdb == orig[N+1] GetNumAtoms: {n_shift}/{total} ({100*n_shift/total:.2f}%)')
print(f'  lmdb == gcsv[N] GetNumAtoms:   {n_gcsv}/{total} ({100*n_gcsv/total:.2f}%)')
print(f'  entries not matching orig[N] nor orig[N+1]: {len(mismatches)}')

# First mismatch point under GetNumAtoms (likely the same row 8)
print('\nFirst 20 rows: key -> (lmdb_atoms, orig[N], orig[N+1], gcsv[N])')
for k in range(20):
    with env.begin() as txn:
        raw = txn.get(str(k).encode())
    d = pickle.loads(raw)
    lmdb_atoms = d['embedding'].shape[0]
    oN = natoms(orig.iloc[k]['Substrate_SMILES'])
    oN1 = natoms(orig.iloc[k+1]['Substrate_SMILES'])
    gN = natoms(gcsv.iloc[k]['Substrate_SMILES'])
    tag = ''
    if lmdb_atoms == oN: tag += ' orig✓'
    if lmdb_atoms == oN1: tag += ' shift✓'
    print(f'  k={k:3d}: lmdb={lmdb_atoms:3d}  orig[{k}]={oN}  orig[{k+1}]={oN1}  gcsv[{k}]={gN}{tag}')
