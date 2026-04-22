import pandas as pd
base = '/root/autodl-tmp/EZSpecificity/PathC/P450/data'
orig = pd.read_csv(f'{base}/Substrates.csv')
gcsv = pd.read_csv(f'{base}/features/grover_substrates.csv')
print('orig Substrates.csv:', len(orig))
print('grover_substrates.csv:', len(gcsv))

# Check if grover CSV has *[H]
has_h = gcsv['Substrate_SMILES'].apply(lambda s: '[H]' in str(s) and len(str(s)) < 10)
print('*[H] in grover CSV:', has_h.sum())
if has_h.any():
    print(gcsv[has_h])

# Row-by-row compare first 20
print('\nFirst 20 rows comparison:')
for i in range(20):
    o = orig.iloc[i]['Substrate_SMILES']
    g = gcsv.iloc[i]['Substrate_SMILES'] if i < len(gcsv) else 'N/A'
    flag = '' if o == g else ' <-- DIFF'
    print(f'  {i}: orig={o[:40]!r:50s} grover={str(g)[:40]!r}{flag}')

# Find the first mismatch point
for i in range(min(len(orig), len(gcsv))):
    if orig.iloc[i]['Substrate_SMILES'] != gcsv.iloc[i]['Substrate_SMILES']:
        print(f'\nFIRST MISMATCH at row {i}')
        print(f'  orig[{i}]: {orig.iloc[i]["Substrate_SMILES"]!r}')
        print(f'  grover[{i}]: {gcsv.iloc[i]["Substrate_SMILES"]!r}')
        # check if orig[i+1] == grover[i]
        if i+1 < len(orig):
            if orig.iloc[i+1]['Substrate_SMILES'] == gcsv.iloc[i]['Substrate_SMILES']:
                print(f'  -> Shift by +1: orig[{i+1}] matches grover[{i}]. *[H] was removed at position {i}')
        break
else:
    print('\nAll matching rows identical. Only trailing difference.')
    print(f'  tail orig[-1]: {orig.iloc[-1]["Substrate_SMILES"]!r}')
    print(f'  tail grover[-1]: {gcsv.iloc[-1]["Substrate_SMILES"]!r}')
