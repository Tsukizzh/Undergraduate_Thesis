"""Empirical proof: does grover_fingerprint.lmdb key N encode Substrate Index N or N+1?

Strategy: GROVER stores per-atom embeddings in data["embedding"] with shape (N_atoms, D).
Compare LMDB key N's atom count to:
  (a) Substrates.csv[N] heavy atom count (what the index claims)
  (b) grover_substrates.csv[N] heavy atom count (what was actually fed to GROVER)
  (c) Substrates.csv[N+1] heavy atom count (shift hypothesis)

RDKit provides reliable heavy atom counting via Mol.GetNumHeavyAtoms().
"""
import pandas as pd, lmdb, pickle, numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
orig = pd.read_csv(f'{base}/Substrates.csv')
gcsv = pd.read_csv(f'{base}/features/grover_substrates.csv')
print(f'Substrates.csv rows: {len(orig)}')
print(f'grover_substrates.csv rows: {len(gcsv)}')

env = lmdb.open(f'{base}/features/grover_fingerprint.lmdb', subdir=False, readonly=True, lock=False, readahead=False)
def get_lmdb(key):
    with env.begin() as txn:
        raw = txn.get(str(key).encode())
    return pickle.loads(raw) if raw is not None else None

def n_heavy(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
        return m.GetNumHeavyAtoms() if m is not None else None
    except Exception:
        return None

# Also compute chemprop-style atom count: chemprop uses ALL atoms including H if explicit, but for fingerprint it uses heavy atoms.
# GROVER uses chemprop's MolGraph. Heavy atom count is the canonical check.

print('\n' + '='*90)
print('Empirical verification: LMDB key N atom count vs Substrate Index N/N+1')
print('='*90)
print(f'{"key":>6} | {"lmdb_atoms":>10} | {"orig[N]_atoms":>13} | {"orig[N+1]_atoms":>15} | {"gcsv[N]_atoms":>13} | {"match":>20}')
print('-'*110)

test_keys = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 500, 1000, 1500, 2000, 2120, 2122, 2123]
mismatches_orig_N = 0
matches_orig_Nplus1 = 0
matches_orig_N = 0

for k in test_keys:
    d = get_lmdb(k)
    if d is None:
        print(f'{k:>6} | MISSING')
        continue
    lmdb_atoms = d['embedding'].shape[0]

    orig_smi = orig.iloc[k]['Substrate_SMILES']
    orig_n = n_heavy(orig_smi)

    orig_next_smi = orig.iloc[k+1]['Substrate_SMILES'] if k+1 < len(orig) else None
    orig_next_n = n_heavy(orig_next_smi) if orig_next_smi is not None else None

    gcsv_smi = gcsv.iloc[k]['Substrate_SMILES'] if k < len(gcsv) else None
    gcsv_n = n_heavy(gcsv_smi) if gcsv_smi is not None else None

    # Determine match
    tags = []
    if orig_n == lmdb_atoms: tags.append('orig[N]✓'); matches_orig_N += 1
    if orig_next_n == lmdb_atoms: tags.append('orig[N+1]✓'); matches_orig_Nplus1 += 1
    if gcsv_n == lmdb_atoms: tags.append('gcsv[N]✓')
    tag = ','.join(tags) if tags else 'NONE'

    print(f'{k:>6} | {lmdb_atoms:>10} | {str(orig_n):>13} | {str(orig_next_n):>15} | {str(gcsv_n):>13} | {tag:>20}')

print('\n' + '='*90)
print('Summary across tested keys:')
print(f'  matches orig[N] (alignment theory):        {matches_orig_N}/{len(test_keys)}')
print(f'  matches orig[N+1] (shift theory):          {matches_orig_Nplus1}/{len(test_keys)}')

# Additional: exhaustive scan of ALL keys, compute match rate
print('\nExhaustive scan over all 2124 LMDB keys:')
n_match_orig = n_match_shift = n_match_gcsv = total = 0
mismatches_both = []
with env.begin() as txn:
    for k_bytes, raw in txn.cursor():
        k = int(k_bytes.decode())
        d = pickle.loads(raw)
        lmdb_atoms = d['embedding'].shape[0]
        oN = n_heavy(orig.iloc[k]['Substrate_SMILES']) if k < len(orig) else None
        oN1 = n_heavy(orig.iloc[k+1]['Substrate_SMILES']) if k+1 < len(orig) else None
        gN = n_heavy(gcsv.iloc[k]['Substrate_SMILES']) if k < len(gcsv) else None
        total += 1
        if oN == lmdb_atoms: n_match_orig += 1
        if oN1 == lmdb_atoms: n_match_shift += 1
        if gN == lmdb_atoms: n_match_gcsv += 1
        if oN != lmdb_atoms and oN1 != lmdb_atoms and gN != lmdb_atoms:
            mismatches_both.append((k, lmdb_atoms, oN, oN1, gN))

print(f'  total LMDB entries: {total}')
print(f'  lmdb_atoms == Substrates.csv[key] atoms:       {n_match_orig}/{total} ({100*n_match_orig/total:.2f}%)')
print(f'  lmdb_atoms == Substrates.csv[key+1] atoms:     {n_match_shift}/{total} ({100*n_match_shift/total:.2f}%)')
print(f'  lmdb_atoms == grover_substrates.csv[key] atoms:{n_match_gcsv}/{total} ({100*n_match_gcsv/total:.2f}%)')
print(f'  entries matching NONE (all three disagree):    {len(mismatches_both)}')
if mismatches_both[:5]:
    print(f'  first 5 unexplained mismatches (k, lmdb, orig[k], orig[k+1], gcsv[k]):')
    for m in mismatches_both[:5]:
        print(f'    {m}')
