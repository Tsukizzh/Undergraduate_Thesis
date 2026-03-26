"""Phase 7 Step 1b: Align raw ligand SDF to substrate SMILES numbering."""
import os, glob
from pathlib import Path
from multiprocessing import Pool
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdFMCS

RAW_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/raw_ligand')
LIG_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/ligand')
LIG_DIR.mkdir(parents=True, exist_ok=True)

data_csv = '/root/rivermind-data/EZSpecificity/PathC/P450/data.csv'
substrate_csv = '/root/rivermind-data/EZSpecificity/PathC/P450/Substrates.csv'

def AssignBondOrdersFromTemplate(refmol, mol):
    refmol2 = Chem.Mol(refmol); mol2 = Chem.Mol(mol)
    matching = mol2.GetSubstructMatch(refmol2)
    if not matching:
        for b in mol2.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
        for b in refmol2.GetBonds(): b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
        for a in refmol2.GetAtoms(): a.SetFormalCharge(0)
        for a in mol2.GetAtoms(): a.SetFormalCharge(0)
        matches = mol2.GetSubstructMatches(refmol2, uniquify=False)
        if not matches: raise ValueError("No matching found")
        for matching in matches:
            for b in refmol.GetBonds():
                a1, a2 = matching[b.GetBeginAtomIdx()], matching[b.GetEndAtomIdx()]
                b2 = mol2.GetBondBetweenAtoms(a1, a2)
                b2.SetBondType(b.GetBondType()); b2.SetIsAromatic(b.GetIsAromatic())
            for a in refmol.GetAtoms():
                a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
                a2.SetHybridization(a.GetHybridization()); a2.SetIsAromatic(a.GetIsAromatic())
                a2.SetNumExplicitHs(a.GetNumExplicitHs()); a2.SetFormalCharge(a.GetFormalCharge())
            try: Chem.SanitizeMol(mol2); break
            except: pass
    return mol2

def alignment_number_system(sdf, smile_mol):
    mcs = rdFMCS.FindMCS([smile_mol, sdf], timeout=120)
    patt = Chem.MolFromSmarts(mcs.smartsString)
    qmatch = sdf.GetSubstructMatch(patt); tmatch = smile_mol.GetSubstructMatch(patt)
    result = [-1] * sdf.GetNumAtoms()
    for q, t in zip(qmatch, tmatch): result[q] = t
    for atom in sdf.GetAtoms():
        assert atom.GetAtomicNum() == 1 or result[atom.GetIdx()] != -1
    return result

def single_match(item):
    idx, smile = item
    RDLogger.DisableLog('rdApp.*')
    sdf_path = RAW_DIR / f'{idx}.sdf'
    out_path = LIG_DIR / f'{idx}.sdf'
    if out_path.exists(): return 0
    try:
        mol = next(iter(Chem.SDMolSupplier(str(sdf_path), sanitize=True)))
        mol = Chem.RemoveHs(mol)
    except:
        try:
            mol = next(iter(Chem.SDMolSupplier(str(sdf_path), sanitize=False)))
            mol = Chem.RemoveHs(mol, sanitize=False)
        except: return 1
    smile_mol = Chem.MolFromSmiles(smile)
    if smile_mol is None: return 1
    try:
        aligned_idx = alignment_number_system(mol, smile_mol)
    except:
        try:
            mol = AssignBondOrdersFromTemplate(smile_mol, mol)
            aligned_idx = alignment_number_system(mol, smile_mol)
        except: return 1
    for atom, ai in zip(mol.GetAtoms(), aligned_idx): atom.SetAtomMapNum(ai)
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() != -1:
            if atom.GetAtomicNum() != smile_mol.GetAtomWithIdx(atom.GetAtomMapNum()).GetAtomicNum():
                return 1
    w = Chem.SDWriter(str(out_path))
    try: w.write(mol); w.close()
    except: w.close(); return 1
    return 0

if __name__ == '__main__':
    meta = pd.read_csv(data_csv)
    subs = pd.read_csv(substrate_csv)
    dock_to_sub = {int(d): int(s) for d, s in zip(meta['Dock Index'], meta['Substrate Index'])}
    smiles = subs['Substrate_SMILES'].values

    params = []
    for p in RAW_DIR.glob('*.sdf'):
        if p.name.startswith('tmp_'): continue
        idx = int(p.stem)
        if (LIG_DIR / f'{idx}.sdf').exists(): continue
        if idx in dock_to_sub:
            params.append((idx, smiles[dock_to_sub[idx]]))

    print(f'Aligning {len(params)} ligands', flush=True)
    bad = 0
    with Pool(10) as pool:
        for i, rc in enumerate(pool.imap_unordered(single_match, params, chunksize=50), 1):
            bad += rc
            if i % 1000 == 0 or i == len(params):
                print(f'[{i}/{len(params)}] bad={bad}', flush=True)
    print(f'DONE. processed={len(params)} bad={bad}')
