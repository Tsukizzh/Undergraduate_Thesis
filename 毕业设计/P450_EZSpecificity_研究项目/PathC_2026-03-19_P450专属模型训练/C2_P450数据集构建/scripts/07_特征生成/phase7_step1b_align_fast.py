"""Phase 7 Step 1b FAST: Align ligands per-substrate (not per-dock). ~25x faster."""
import os
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS
RDLogger.DisableLog('rdApp.*')

RAW_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/raw_ligand')
LIG_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/ligand')
LIG_DIR.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv('/root/rivermind-data/EZSpecificity/PathC/P450/data.csv')
subs = pd.read_csv('/root/rivermind-data/EZSpecificity/PathC/P450/Substrates.csv')
dock_to_sub = {int(d): int(s) for d, s in zip(meta['Dock Index'], meta['Substrate Index'])}
smiles_arr = subs['Substrate_SMILES'].values

# Group by substrate
sub_to_docks = defaultdict(list)
for p in RAW_DIR.glob('*.sdf'):
    if p.name.startswith('tmp_'): continue
    idx = int(p.stem)
    if (LIG_DIR / f'{idx}.sdf').exists(): continue
    if idx in dock_to_sub:
        sub_to_docks[dock_to_sub[idx]].append(idx)

print(f'Unique substrates: {len(sub_to_docks)}, dock_indices: {sum(len(v) for v in sub_to_docks.values())}', flush=True)

def align_substrate(args):
    sub_idx, dock_indices, smile = args
    first = dock_indices[0]
    sdf_path = RAW_DIR / f'{first}.sdf'
    try:
        mol = next(iter(Chem.SDMolSupplier(str(sdf_path), sanitize=True)))
        mol = Chem.RemoveHs(mol)
    except:
        try:
            mol = next(iter(Chem.SDMolSupplier(str(sdf_path), sanitize=False)))
            mol = Chem.RemoveHs(mol, sanitize=False)
        except:
            return len(dock_indices)

    smile_mol = Chem.MolFromSmiles(smile)
    if smile_mol is None:
        return len(dock_indices)

    # MCS alignment ONCE per substrate
    try:
        mcs = rdFMCS.FindMCS([smile_mol, mol], timeout=120)
        patt = Chem.MolFromSmarts(mcs.smartsString)
        qmatch = mol.GetSubstructMatch(patt)
        tmatch = smile_mol.GetSubstructMatch(patt)
        mapping = [-1] * mol.GetNumAtoms()
        for q, t in zip(qmatch, tmatch):
            mapping[q] = t
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1 and mapping[atom.GetIdx()] == -1:
                raise ValueError('unmapped')
    except:
        try:
            mol2 = Chem.Mol(mol); refmol = Chem.Mol(smile_mol)
            for b in mol2.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
            for b in refmol.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
            for a in refmol.GetAtoms(): a.SetFormalCharge(0)
            for a in mol2.GetAtoms(): a.SetFormalCharge(0)
            matches = mol2.GetSubstructMatches(refmol, uniquify=False)
            if not matches: return len(dock_indices)
            matching = matches[0]
            mapping = [-1] * mol.GetNumAtoms()
            for ref_idx, mol_idx in enumerate(matching): mapping[mol_idx] = ref_idx
        except:
            return len(dock_indices)

    # Apply to all dock_indices
    bad = 0
    for didx in dock_indices:
        out = LIG_DIR / f'{didx}.sdf'
        if out.exists(): continue
        try:
            m = next(iter(Chem.SDMolSupplier(str(RAW_DIR / f'{didx}.sdf'), sanitize=False)))
            m = Chem.RemoveHs(m, sanitize=False)
            if m.GetNumAtoms() != len(mapping): bad += 1; continue
            for atom, idx in zip(m.GetAtoms(), mapping): atom.SetAtomMapNum(idx)
            w = Chem.SDWriter(str(out))
            w.write(m); w.close()
        except:
            bad += 1
    return bad

tasks = [(si, docks, smiles_arr[si]) for si, docks in sub_to_docks.items()]
print(f'Tasks: {len(tasks)}', flush=True)

total_bad = 0
with Pool(12) as pool:
    for i, bad in enumerate(pool.imap_unordered(align_substrate, tasks, chunksize=10), 1):
        total_bad += bad
        if i % 200 == 0 or i == len(tasks):
            print(f'[{i}/{len(tasks)}] bad={total_bad}', flush=True)

final = len(list(LIG_DIR.iterdir()))
print(f'DONE. Total ligands: {final}, bad: {total_bad}')
