"""
Phase 6 Stage B: Substrate SDF → PDBQT (Meeko)
- Converts 2,124 substrate SDF files to PDBQT format for AutoDock
- Multi-worker parallel processing
- Resume support: skips already completed files
- Updates substrate_registry.csv with PDBQT status

Codex-reviewed: 2026-03-25
"""

import pandas as pd
import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

# Use parent (not resolve) to preserve Chinese characters in paths
BASE = Path(__file__).parent.parent.parent
DATA = BASE / "data"
REGISTRIES = DATA / "registries"
SDF_DIR = DATA / "structures" / "ligands" / "sdf"
OUT_DIR = DATA / "structures" / "ligands_pdbqt"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MEEKO = "D:/anaconda3/envs/torch/Scripts/mk_prepare_ligand.exe"
PYTHON = "D:/anaconda3/envs/torch/python.exe"
MAX_WORKERS = 8


def validate_sdf(sdf_path):
    """Pre-validate SDF before Meeko."""
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if len(mols) == 0:
        return False, "no_valid_molecule"
    if len(mols) > 1:
        return False, f"multiple_molecules({len(mols)})"
    mol = mols[0]
    if mol.GetNumAtoms() == 0:
        return False, "zero_atoms"
    if mol.GetNumConformers() == 0:
        return False, "no_conformer"
    # Check for dummy atoms (wildcard *) — Meeko cannot handle query atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            return False, "contains_dummy_atom(*)"
    # Check for metals (would cause Meeko issues)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in {26, 30, 29, 25, 27, 28}:  # Fe, Zn, Cu, Mn, Co, Ni
            return False, f"contains_metal({atom.GetSymbol()})"
    return True, "ok"


def process_one_ligand(args):
    """Process a single SDF → PDBQT. Runs in subprocess worker."""
    substrate_index, sdf_path, out_path = args
    sdf_path = Path(sdf_path)
    out_path = Path(out_path)

    # Skip if already done
    ok_marker = out_path.with_suffix(".ok")
    if ok_marker.exists():
        return substrate_index, "skipped", "", ""

    if not sdf_path.exists():
        return substrate_index, "sdf_not_found", "", str(sdf_path)

    # All operations in temp dir to avoid Chinese path issues with RDKit/Meeko
    import tempfile, shutil, os
    try:
        with tempfile.TemporaryDirectory(prefix="meeko_lig_") as tmpdir:
            tmp_sdf = Path(tmpdir) / f"{sdf_path.stem}.sdf"
            tmp_pdbqt = Path(tmpdir) / f"{sdf_path.stem}.pdbqt"
            shutil.copy2(str(sdf_path), str(tmp_sdf))

            # Validate SDF from temp path (English path)
            valid, reason = validate_sdf(tmp_sdf)
            if not valid:
                return substrate_index, f"validation_failed({reason})", "", ""

            result = subprocess.run(
                [MEEKO, "-i", str(tmp_sdf), "-o", str(tmp_pdbqt), "--add_index_map"],
                capture_output=True, text=True, timeout=120
            )
            if result.returncode != 0:
                return substrate_index, "meeko_failed", result.stdout, result.stderr[:300]
            if not tmp_pdbqt.exists() or tmp_pdbqt.stat().st_size == 0:
                return substrate_index, "empty_output", result.stdout, result.stderr[:300]

            # Copy result to final location (copy2 handles Chinese paths better than move)
            shutil.copy2(str(tmp_pdbqt), str(out_path))

        # Write success marker
        ok_marker.write_text("ok")
        return substrate_index, "ok", "", ""
    except subprocess.TimeoutExpired:
        return substrate_index, "timeout", "", ""
    except Exception as e:
        return substrate_index, f"exception({type(e).__name__})", "", str(e)


def main():
    print("=" * 60)
    print("Phase 6 Stage B: Substrate SDF → PDBQT")
    print("=" * 60)

    # Load substrate registry
    registry = pd.read_csv(REGISTRIES / "substrate_registry.csv")
    ok_subs = registry[registry["ligand_status"] == "ok"].copy()
    print(f"Total substrates: {len(registry)}, OK for conversion: {len(ok_subs)}")

    # Build task list
    tasks = []
    for _, row in ok_subs.iterrows():
        sidx = int(row["substrate_index"])
        cid = row["global_compound_id"]
        sdf_path = SDF_DIR / f"{cid}.sdf"
        out_path = OUT_DIR / f"{cid}.pdbqt"
        tasks.append((sidx, str(sdf_path), str(out_path)))

    # Count already done
    already_done = sum(1 for _, _, out in tasks if Path(out).with_suffix(".ok").exists())
    remaining = len(tasks) - already_done
    print(f"Already done: {already_done}, Remaining: {remaining}")

    if remaining == 0:
        print("All ligands already converted!")
        return

    # Process in parallel
    print(f"Processing {remaining} ligands with {MAX_WORKERS} workers...")
    t0 = time.time()
    results = {"ok": 0, "skipped": 0, "failed": 0}
    errors = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_one_ligand, t): t[0] for t in tasks}
        for i, future in enumerate(as_completed(futures), 1):
            sidx, status, stdout, stderr = future.result()
            if status == "ok":
                results["ok"] += 1
            elif status == "skipped":
                results["skipped"] += 1
            else:
                results["failed"] += 1
                errors.append({"substrate_index": sidx, "status": status, "stderr": stderr[:200]})

            if i % 100 == 0 or i == len(tasks):
                elapsed = time.time() - t0
                print(f"  [{i}/{len(tasks)}] ok={results['ok']} skip={results['skipped']} fail={results['failed']} ({elapsed:.0f}s)")

    # Update registry
    pdbqt_status = {}
    pdbqt_paths = {}
    for _, sdf, out in tasks:
        cid = Path(sdf).stem
        out_path = Path(out)
        if out_path.with_suffix(".ok").exists():
            pdbqt_status[cid] = "ok"
            pdbqt_paths[cid] = str(out_path)
        else:
            pdbqt_status[cid] = "failed"
            pdbqt_paths[cid] = ""

    registry["pdbqt_status"] = registry["global_compound_id"].map(pdbqt_status).fillna("")
    registry["pdbqt_path"] = registry["global_compound_id"].map(pdbqt_paths).fillna("")
    registry.to_csv(REGISTRIES / "substrate_registry.csv", index=False)

    # Summary
    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"Newly converted: {results['ok']}")
    print(f"Skipped (already done): {results['skipped']}")
    print(f"Failed: {results['failed']}")
    print(f"Time: {elapsed:.1f}s ({elapsed/60:.1f}min)")
    if errors:
        print(f"\nFailed ligands:")
        for e in errors[:20]:
            print(f"  substrate_index={e['substrate_index']}: {e['status']} | {e['stderr']}")


if __name__ == "__main__":
    main()
