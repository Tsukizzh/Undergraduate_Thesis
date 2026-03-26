"""
Phase 6.5: Uni-Dock post-processing -> EZSpecificity complex PDB
Converts docking outputs to: [receptor ATOM/HETATM] + COMPND + [RDKit-parseable ligand]
Uses original SDF templates to preserve bond orders.
Template normalization: strip dummy atoms, keep largest fragment, remove Hs.
Validates against input PDBQT (pre-docking) as authoritative atom count.
"""
from __future__ import annotations
import argparse, csv, os, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from rdkit import Chem
from rdkit.Geometry import Point3D

PTABLE = Chem.GetPeriodicTable()
PDBQT_TYPE_TO_ELEM = {
    "A":"C","C":"C","N":"N","NA":"N","OA":"O","O":"O","S":"S","SA":"S",
    "P":"P","F":"F","CL":"Cl","BR":"Br","I":"I","HD":"H","HS":"H",
    "FE":"Fe","ZN":"Zn","MG":"Mg","MN":"Mn","CA":"Ca","CU":"Cu","NI":"Ni","CO":"Co",
}
TWO_LETTER_ELEMENTS = {"CL","BR","FE","ZN","MG","MN","CA","CU","NI","CO"}

@dataclass(frozen=True)
class Atom:
    element: str; x: float; y: float; z: float

def infer_element(line, atom_name):
    atom_name = atom_name.strip().upper()
    if atom_name.startswith("H"):
        return "H"
    tokens = line.split()
    if tokens:
        m = PDBQT_TYPE_TO_ELEM.get(tokens[-1].strip().upper())
        if m: return m
    e = line[76:78].strip()
    if e:
        return e[0].upper() + e[1:].lower() if len(e) > 1 else e.upper()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if len(letters) >= 2 and letters[:2].upper() in TWO_LETTER_ELEMENTS:
        return letters[0].upper() + letters[1].lower()
    return letters[0].upper() if letters else "C"

def parse_pdbqt_heavy(path):
    atoms = []; in_first = False; saw = False
    with path.open("r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            rec = line[:6].strip()
            if rec == "MODEL":
                if not saw: saw = True; in_first = True; continue
                else: break
            if rec == "ENDMDL" and in_first: break
            if saw and not in_first: continue
            if rec not in {"ATOM","HETATM"}: continue
            atom_name = line[12:16].strip()
            elem = infer_element(line, atom_name)
            if elem.upper() == "H": continue
            atoms.append(Atom(elem, float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not atoms: raise ValueError(f"No heavy atoms in {path}")
    return atoms

def receptor_lines(path):
    lines = []
    with path.open("r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if line[:6].strip() in {"ATOM","HETATM","TER"}: lines.append(line)
    return lines

def strip_dummy_atoms(mol):
    dummy_idxs = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if not dummy_idxs: return mol
    rw = Chem.RWMol(mol)
    for idx in reversed(dummy_idxs):
        rw.RemoveAtom(idx)
    mol = rw.GetMol()
    try: Chem.SanitizeMol(mol)
    except: pass
    return mol

def largest_fragment(mol):
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) <= 1: return mol
    return max(frags, key=lambda m: sum(1 for a in m.GetAtoms() if a.GetAtomicNum() > 1))

def normalize_template(mol):
    mol = strip_dummy_atoms(Chem.Mol(mol))
    mol = largest_fragment(mol)
    mol = Chem.RemoveHs(mol, sanitize=True)
    mol = strip_dummy_atoms(mol)
    mol = largest_fragment(mol)
    if mol.GetNumAtoms() == 0:
        raise ValueError("Template normalization produced zero atoms")
    return mol

def build_ligand(template_mol, docked_atoms):
    n_tmpl = template_mol.GetNumAtoms()
    n_dock = len(docked_atoms)
    if n_tmpl != n_dock:
        raise ValueError(f"Atom count mismatch: template {n_tmpl} vs docked {n_dock}")
    conf = Chem.Conformer(n_tmpl)
    for i, a in enumerate(docked_atoms):
        conf.SetAtomPosition(i, Point3D(a.x, a.y, a.z))
    mol = Chem.Mol(template_mol)
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    ec = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol(); ec[sym] = ec.get(sym, 0) + 1
        info = Chem.AtomPDBResidueInfo()
        info.SetName(f" {sym}{ec[sym]:<2d}"[:4])
        info.SetResidueName("LIG"); info.SetResidueNumber(1)
        info.SetChainId("Z"); info.SetIsHeteroAtom(True)
        info.SetSerialNumber(atom.GetIdx()+1)
        atom.SetMonomerInfo(info)
    return mol

def ligand_pdb_lines(mol):
    block = Chem.MolToPDBBlock(mol)
    return [l for l in block.splitlines() if l[:6].strip() in {"HETATM","ATOM","CONECT"}]

def process_one(row, args_d):
    dk = int(row["dock_index"]); ei = int(row["enzyme_index"]); si = int(row["substrate_index"])
    dd = Path(args_d["docking_dir"]); sd = Path(args_d["sdf_dir"]); od = Path(args_d["output_dir"])
    out = od / f"{dk}.pdb"
    if out.exists() and not args_d["overwrite"] and out.stat().st_size > 0:
        return {"dock_index": str(dk), "status": "skipped"}

    rec_path = dd / "receptors_clean" / f"ENZ_G{ei+1:06d}.pdb"
    inp_pdbqt = dd / "ligands_pdbqt" / f"CMP_G{si+1:06d}.pdbqt"
    dock_pdbqt = dd / "results" / f"enz_{ei:06d}" / f"CMP_G{si+1:06d}_out.pdbqt"
    sdf_path = sd / f"CMP_G{si+1:06d}.sdf"

    for p, tag in [(rec_path,"receptor"),(inp_pdbqt,"input_ligand"),(dock_pdbqt,"docked")]:
        if not p.exists():
            return {"dock_index": str(dk), "status": "failed", "reason": f"missing_{tag}"}

    try:
        rec = receptor_lines(rec_path)
        inp_ha = parse_pdbqt_heavy(inp_pdbqt)
        dock_ha = parse_pdbqt_heavy(dock_pdbqt)
        if len(inp_ha) != len(dock_ha):
            raise ValueError(f"Dock/input mismatch: input {len(inp_ha)} vs docked {len(dock_ha)}")

        # Load and normalize template
        if sdf_path.exists():
            sup = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            tmpl = next((m for m in sup if m is not None), None)
            if tmpl is None: raise ValueError(f"Bad SDF: {sdf_path}")
            tmpl = normalize_template(tmpl)
        elif args_d["allow_smiles_fallback"]:
            smiles = row.get("Substrate_SMILES", "")
            if smiles and "*" in smiles:
                smiles = smiles.replace("*", "[H]")
            tmpl = Chem.MolFromSmiles(smiles) if smiles else None
            if tmpl is None: raise ValueError("Bad SMILES fallback")
            tmpl = normalize_template(tmpl)
        else:
            raise FileNotFoundError(f"No SDF: {sdf_path}")

        if tmpl.GetNumAtoms() != len(inp_ha):
            raise ValueError(
                f"Template/input mismatch: template {tmpl.GetNumAtoms()} vs input {len(inp_ha)}"
            )

        lig_mol = build_ligand(tmpl, dock_ha)
        lig_lines = ligand_pdb_lines(lig_mol)

        out.parent.mkdir(parents=True, exist_ok=True)
        with open(out, "w") as f:
            for l in rec: f.write(l + "\n")
            f.write(f"COMPND    {dock_pdbqt.name}\n")
            for l in lig_lines: f.write(l + "\n")
            f.write("END\n")
        return {"dock_index": str(dk), "status": "ok"}
    except Exception as e:
        return {"dock_index": str(dk), "status": "failed", "reason": f"{type(e).__name__}: {str(e)[:200]}"}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--docking-dir", type=Path, required=True)
    p.add_argument("--sdf-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--workers", type=int, default=14)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--allow-smiles-fallback", action="store_true")
    p.add_argument("--failure-csv", type=Path, default=None)
    args = p.parse_args()

    reg_dir = args.docking_dir / "registries"
    rows = []
    for csv_path in sorted(reg_dir.glob("pair_registry_batch_*.csv")):
        with open(csv_path, encoding="utf-8-sig") as f:
            rows.extend(csv.DictReader(f))
    print(f"Loaded {len(rows)} pairs. Workers={args.workers}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    wa = {"docking_dir": str(args.docking_dir), "sdf_dir": str(args.sdf_dir),
          "output_dir": str(args.output_dir), "overwrite": args.overwrite,
          "allow_smiles_fallback": args.allow_smiles_fallback}
    counts = {"ok": 0, "skipped": 0, "failed": 0}
    failures = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(process_one, r, wa) for r in rows]
        for i, fut in enumerate(as_completed(futs), 1):
            res = fut.result()
            counts[res["status"]] = counts.get(res["status"], 0) + 1
            if res["status"] == "failed": failures.append(res)
            if i % 500 == 0 or i == len(futs):
                print(f"[{i}/{len(futs)}] ok={counts['ok']} skip={counts['skipped']} fail={counts['failed']}")
    print(f"DONE. ok={counts['ok']} skipped={counts['skipped']} failed={counts['failed']}")
    if args.failure_csv and failures:
        with open(args.failure_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["dock_index","status","reason"])
            w.writeheader(); w.writerows(failures)
        print(f"Failures: {args.failure_csv}")

if __name__ == "__main__":
    main()
