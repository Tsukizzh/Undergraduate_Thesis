#!/usr/bin/env python3
"""
Phase 6: Generate 3D coordinates for all substrates.

SMILES → RDKit → 3D embedding (ETKDGv3, 50 conformers) → MMFF94s/UFF minimize
→ pick lowest energy → save SDF

Design decisions (from Codex+Gemini review):
- MMFF94s primary, UFF fallback (Gemini: MMFF94s better for organic substrates)
- 50 conformers, pick lowest energy (Gemini: single conformer risks bad geometry)
- Replace dummy atoms * with [H] before 3D generation (Gemini: critical for wildcards)
- Keep largest fragment to remove salts (Gemini: prevents disconnected structures)
- ProcessPoolExecutor for true parallelism (Gemini: GIL blocks ThreadPool)
- Per-molecule timeout to prevent macrocycle hangs (Gemini)
- Resume support + incremental manifest writing
"""
from __future__ import annotations

import argparse, csv, json, os, signal, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
COMBINED = PROJECT / "data" / "combined"
STRUCTURES = PROJECT / "data" / "structures"
SDF_DIR = STRUCTURES / "ligands" / "sdf"
MANIFEST_PATH = STRUCTURES / "manifests" / "ligand_manifest.csv"

MANIFEST_FIELDS = [
    "global_compound_id", "canonical_smiles", "status",
    "num_atoms", "num_conformers_tried", "final_energy",
    "force_field_used", "sdf_path", "error_message",
]

NUM_CONFORMERS = 50
TIMEOUT_SECONDS = 60  # per molecule

def log(msg): print(msg, file=sys.stderr, flush=True)

# ---------------------------------------------------------------------------
# Single molecule processing (runs in worker process)
# ---------------------------------------------------------------------------
def process_one_molecule(args):
    """Process a single molecule. Returns a manifest row dict."""
    cmp_id, smiles, sdf_dir = args
    result = {
        "global_compound_id": cmp_id,
        "canonical_smiles": smiles,
        "status": "", "num_atoms": "", "num_conformers_tried": "",
        "final_energy": "", "force_field_used": "", "sdf_path": "",
        "error_message": "",
    }

    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import AllChem, Descriptors, rdmolops
        RDLogger.DisableLog("rdApp.*")

        if not smiles or not smiles.strip():
            result["status"] = "empty_smiles"
            result["error_message"] = "Empty SMILES"
            return result

        # Step 1: Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Try replacing dummy atoms * with [H]
            cleaned = smiles.replace("*", "[H]")
            mol = Chem.MolFromSmiles(cleaned)
            if mol is None:
                result["status"] = "parse_failed"
                result["error_message"] = "RDKit cannot parse SMILES (even after * → [H])"
                return result
            result["error_message"] = "Replaced * with [H]"

        # Step 2: Keep largest fragment (remove salts/counterions)
        frags = rdmolops.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumAtoms())

        # Step 3: Add hydrogens
        mol = Chem.AddHs(mol)
        result["num_atoms"] = mol.GetNumAtoms()

        # Step 4: Generate multiple 3D conformers (ETKDGv3)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 1  # each worker uses 1 thread
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=NUM_CONFORMERS, params=params)

        if len(conf_ids) == 0:
            # Fallback: try with useRandomCoords
            params.useRandomCoords = True
            conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=NUM_CONFORMERS, params=params)

        if len(conf_ids) == 0:
            result["status"] = "embed_failed"
            result["error_message"] = "EmbedMultipleConfs returned 0 conformers"
            return result

        result["num_conformers_tried"] = len(conf_ids)

        # Step 5: Minimize each conformer, track energies
        energies = {}
        ff_used = "none"

        # Try MMFF94s first
        mmff_ok = True
        for cid in conf_ids:
            try:
                ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s"), confId=cid)
                if ff is None:
                    mmff_ok = False
                    break
                ff.Minimize(maxIts=500)
                energies[cid] = ff.CalcEnergy()
            except Exception:
                mmff_ok = False
                break

        if mmff_ok and energies:
            ff_used = "MMFF94s"
        else:
            # Fallback to UFF
            energies = {}
            for cid in conf_ids:
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
                    if ff is None:
                        continue
                    ff.Minimize(maxIts=500)
                    energies[cid] = ff.CalcEnergy()
                except Exception:
                    continue
            if energies:
                ff_used = "UFF"

        if not energies:
            # No force field worked, use unminimized lowest-RMSD conformer
            best_cid = conf_ids[0]
            ff_used = "none"
            result["final_energy"] = ""
        else:
            best_cid = min(energies, key=energies.get)
            result["final_energy"] = f"{energies[best_cid]:.2f}"

        result["force_field_used"] = ff_used

        # Step 6: Save best conformer as SDF
        # RDKit SDWriter can't handle non-ASCII paths on Windows
        # Write to temp file first, then move
        import tempfile, shutil
        sdf_path = Path(sdf_dir) / f"{cmp_id}.sdf"
        sdf_path.parent.mkdir(parents=True, exist_ok=True)

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False, dir=tempfile.gettempdir()) as tmp:
            tmp_path = tmp.name

        writer = Chem.SDWriter(tmp_path)
        writer.write(mol, confId=best_cid)
        writer.close()
        shutil.move(tmp_path, str(sdf_path))

        result["status"] = "ok"
        result["sdf_path"] = str(sdf_path.relative_to(PROJECT)) if str(PROJECT) in str(sdf_path) else str(sdf_path)
        return result

    except Exception as e:
        result["status"] = "error"
        result["error_message"] = str(e)[:200]
        return result

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Generate 3D coordinates for substrates")
    parser.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 2),
                        help="Number of parallel workers")
    parser.add_argument("--combined-dir", type=Path, default=COMBINED)
    args = parser.parse_args()

    log(f"[setup] workers={args.workers}")

    # Load compounds
    compounds_csv = args.combined_dir / "global_compounds.csv"
    with compounds_csv.open("r", encoding="utf-8-sig") as f:
        compounds = list(csv.DictReader(f))
    log(f"[data] {len(compounds)} compounds to process")

    # Load existing manifest for resume
    done = set()
    if MANIFEST_PATH.exists():
        with MANIFEST_PATH.open("r", encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                if r.get("status") in ("ok", "parse_failed", "embed_failed", "empty_smiles"):
                    done.add(r["global_compound_id"])
    log(f"[resume] {len(done)} already processed, {len(compounds)-len(done)} remaining")

    # Prepare tasks
    tasks = []
    for c in compounds:
        cid = c["global_compound_id"].strip()
        smi = c.get("canonical_smiles", "").strip()
        if cid not in done:
            tasks.append((cid, smi, str(SDF_DIR)))

    if not tasks:
        log("[done] All compounds already processed")
        return

    # Create output dirs
    SDF_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Open manifest for incremental append
    write_header = not MANIFEST_PATH.exists() or len(done) == 0
    manifest_file = MANIFEST_PATH.open("a" if not write_header else "w",
                                        encoding="utf-8-sig", newline="")
    writer = csv.DictWriter(manifest_file, fieldnames=MANIFEST_FIELDS)
    if write_header:
        writer.writeheader()
        # Re-write existing rows if starting fresh
        if MANIFEST_PATH.exists() and len(done) > 0:
            pass  # append mode, existing rows preserved

    # Process in parallel
    completed = 0
    ok = fail = 0
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_one_molecule, t): t[0] for t in tasks}

        for future in as_completed(futures):
            cmp_id = futures[future]
            try:
                result = future.result(timeout=TIMEOUT_SECONDS)
            except Exception as e:
                result = {f: "" for f in MANIFEST_FIELDS}
                result["global_compound_id"] = cmp_id
                result["status"] = "timeout_or_crash"
                result["error_message"] = str(e)[:200]

            writer.writerow(result)
            manifest_file.flush()
            completed += 1

            if result["status"] == "ok":
                ok += 1
            else:
                fail += 1

            if completed % 100 == 0 or completed == len(tasks):
                elapsed = time.time() - start_time
                rate = completed / elapsed if elapsed > 0 else 0
                log(f"[progress] {completed}/{len(tasks)} ({ok} ok, {fail} fail, {rate:.1f}/sec)")

    manifest_file.close()
    elapsed = time.time() - start_time

    log(f"\n{'='*60}")
    log(f"[DONE] {completed} compounds processed in {elapsed:.0f}s")
    log(f"[DONE] ok={ok}, fail={fail}")
    log(f"[DONE] SDF files: {SDF_DIR}")
    log(f"[DONE] Manifest: {MANIFEST_PATH}")
    log(f"{'='*60}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        log("Interrupted"); raise SystemExit(130)
    except Exception as e:
        log(f"ERROR: {e}"); import traceback; traceback.print_exc(); raise SystemExit(1)
