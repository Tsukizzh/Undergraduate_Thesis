"""Phase 7 Step 1: Extract pocket + raw ligand from complex PDBs. 10 CPU workers."""
import os, glob
from pathlib import Path
from multiprocessing import Pool
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

COMPLEX_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/complex')
POCKET_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/pocket')
RAW_DIR = Path('/root/rivermind-data/EZSpecificity/PathC/P450/structure/str_tmp_data/raw_ligand')
POCKET_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

def distance(x1,y1,z1,x2,y2,z2):
    return ((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)**0.5

def process_one(name):
    name = Path(name)
    idx = name.stem
    pocket_out = POCKET_DIR / f'{idx}.pdb'
    ligand_sdf = RAW_DIR / f'{idx}.sdf'
    if pocket_out.exists() and ligand_sdf.exists():
        return 0, 0  # skip
    lines = name.read_text(errors='replace').splitlines(True)
    sep = None
    for i, line in enumerate(lines):
        if 'COMPND' in line: sep = i; break
    if sep is None: return 1, 0
    protein_lines = lines[:sep]
    ligand_lines = lines[sep+1:]
    tmp_pdb = RAW_DIR / f'tmp_{idx}.pdb'
    with open(tmp_pdb, 'w') as f: f.writelines(ligand_lines)
    mol = Chem.MolFromPDBFile(str(tmp_pdb), flavor=1, sanitize=False)
    try: tmp_pdb.unlink()
    except: pass
    if mol is None: return 1, 0
    coords = []
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        coords.append((p.x, p.y, p.z))
    try:
        w = Chem.SDWriter(str(ligand_sdf))
        w.write(mol, confId=0); w.close()
    except: return 1, 0
    with open(pocket_out, 'w') as out:
        for line in protein_lines:
            rec = line[:6].strip()
            if rec == 'ATOM':
                try:
                    x=float(line[30:38]); y=float(line[38:46]); z=float(line[46:54])
                    if 'H' in line[12:16].strip(): continue
                    for lx,ly,lz in coords:
                        if distance(x,y,z,lx,ly,lz) < 10:
                            out.write(line); break
                except: continue
            elif rec in {'HETATM','ENDMDL','TER'}:
                out.write(line)
    return 0, 1

if __name__ == '__main__':
    files = sorted(glob.glob(str(COMPLEX_DIR / '*.pdb')))
    print(f'Processing {len(files)} complex PDBs with 10 workers', flush=True)
    with Pool(10) as pool:
        ok = bad = skip = 0
        for rc, is_ok in pool.imap_unordered(process_one, files, chunksize=50):
            if rc == 0 and is_ok: ok += 1
            elif rc == 0: skip += 1
            else: bad += 1
            if (ok+bad+skip) % 2000 == 0:
                print(f'  ok={ok} bad={bad} skip={skip}', flush=True)
    print(f'DONE. ok={ok} bad={bad} skip={skip}')
