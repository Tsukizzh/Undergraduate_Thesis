#!/usr/bin/env python3
"""
Process ColabFold outputs for Batch 2:
  1. Extract best model (rank_001) from ColabFold output
  2. Copy to alphafold/pdb/ with compatible naming
  3. Run HEM transplant (reuse transplant_heme_v2 logic)
  4. Run PDB cleanup → receptors_clean/
  5. Run MGLTools PDBQT conversion → receptors_pdbqt/
  6. Update enzyme_registry.csv with Fe coords

Usage:
  python process_colabfold_batch2.py
"""

import csv, os, re, shutil, sys, subprocess
from pathlib import Path
from collections import defaultdict

SCRIPT = Path(__file__).resolve()
PROJECT = SCRIPT.parents[2]
STRUCTURES = PROJECT / "data" / "structures"
REGISTRIES = PROJECT / "data" / "registries"

COLABFOLD_OUT = STRUCTURES / "colabfold_output"
AF_PDB_DIR = STRUCTURES / "alphafold" / "pdb"
HEME_DIR = STRUCTURES / "heme_transplant" / "pdb"
CLEAN_DIR = STRUCTURES / "receptors_clean"
PDBQT_DIR = STRUCTURES / "receptors_pdbqt"

def log(msg):
    print(msg, flush=True)


def step1_extract_best_models():
    """Extract rank_001 unrelaxed PDBs, copy to alphafold/pdb/ with AF-{UniProt}-F1.pdb naming."""
    log("=== Step 1: Extract best ColabFold models ===")
    AF_PDB_DIR.mkdir(parents=True, exist_ok=True)

    # Parse ColabFold filenames: ENZ_GXXXXXX_UniProtID_reason_unrelaxed_rank_001_...pdb
    copied = 0
    skipped = 0
    enz_to_uniprot = {}  # ENZ_GXXXXXX -> UniProtID

    for fname in sorted(os.listdir(COLABFOLD_OUT)):
        if not fname.endswith(".pdb"):
            continue
        if "unrelaxed" not in fname or "rank_001" not in fname:
            continue

        # Extract ENZ_ID and UniProt from filename
        # Format: ENZ_G000014_A0A067ZUX1_no_template_unrelaxed_rank_001_...pdb
        parts = fname.split("_")
        enz_id = parts[0] + "_" + parts[1]  # ENZ_GXXXXXX

        # Find UniProt ID (second meaningful part after ENZ_GXXXXXX)
        # It's between the ENZ_G part and the reason/unrelaxed part
        # Could be a real UniProt ID or SEQHASH or "nan"
        uniprot = parts[2] if len(parts) > 2 else ""

        # For SEQHASH or nan entries, use ENZ_ID as identifier
        if uniprot.startswith("SEQHASH") or uniprot == "nan":
            dest_name = f"CF-{enz_id}.pdb"
        else:
            dest_name = f"AF-{uniprot}-F1.pdb"

        src = COLABFOLD_OUT / fname
        dst = AF_PDB_DIR / dest_name

        if dst.exists():
            skipped += 1
            continue

        shutil.copy2(src, dst)
        enz_to_uniprot[enz_id] = (uniprot, dest_name)
        copied += 1

    log(f"  Copied {copied} models, skipped {skipped} (already exist)")
    log(f"  Total in {AF_PDB_DIR}: {len(list(AF_PDB_DIR.glob('*.pdb')))}")
    return enz_to_uniprot


def step2_run_heme_transplant():
    """Run transplant_heme_v2.py to add HEM to ColabFold structures."""
    log("\n=== Step 2: HEM transplant ===")
    transplant_script = PROJECT / "scripts" / "04_结构获取" / "transplant_heme_v2.py"
    if not transplant_script.exists():
        log(f"  ERROR: {transplant_script} not found!")
        return False

    result = subprocess.run(
        [sys.executable, str(transplant_script), "--workers", "6"],
        capture_output=True, text=True, timeout=3600
    )
    print(result.stdout)
    if result.stderr:
        print(result.stderr, file=sys.stderr)

    if result.returncode != 0:
        log(f"  WARNING: transplant_heme_v2.py exited with code {result.returncode}")
    return True


def step3_cleanup_and_pdbqt():
    """
    For each new heme-transplanted PDB:
      - Clean PDB → receptors_clean/
      - Convert to PDBQT via MGLTools (or just copy to PDBQT dir for later Cloud-2 processing)
    """
    log("\n=== Step 3: PDB cleanup ===")
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    PDBQT_DIR.mkdir(parents=True, exist_ok=True)

    # Read enzyme_registry to find Batch 2 enzymes
    registry_path = REGISTRIES / "enzyme_registry.csv"
    batch2_enzymes = {}
    with open(registry_path, encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            eidx = int(row["enzyme_index"])
            gid = row.get("global_enzyme_id", "")
            batch2_enzymes[gid] = eidx

    # Find new heme PDBs that don't have clean versions yet
    new_count = 0
    for pdb in sorted(HEME_DIR.glob("*_heme.pdb")):
        # Find corresponding ENZ_G ID
        stem = pdb.stem  # e.g. AF-A0A067ZUX1-F1_heme or CF-ENZ_G001404_heme
        clean_name = None

        # Try to match to enzyme_registry
        for gid, eidx in batch2_enzymes.items():
            enz_file = f"ENZ_G{eidx+1:06d}"
            clean_path = CLEAN_DIR / f"{enz_file}.pdb"
            pdbqt_path = PDBQT_DIR / f"{enz_file}.pdbqt"

            if clean_path.exists():
                continue

            # Check if this heme PDB matches this enzyme
            uid = gid.replace("ENZ_G", "").strip()
            if f"AF-{uid}" in stem or f"CF-{gid}" in stem:
                shutil.copy2(pdb, clean_path)
                new_count += 1
                break

    log(f"  {new_count} new clean PDBs created")


def step4_update_registry():
    """Update enzyme_registry.csv with Fe coordinates from new heme PDBs."""
    log("\n=== Step 4: Update enzyme_registry with Fe coords ===")
    registry_path = REGISTRIES / "enzyme_registry.csv"

    with open(registry_path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    updated = 0
    for row in rows:
        eidx = int(row["enzyme_index"])
        fe_x = row.get("heme_fe_x", "")
        if fe_x:  # already has Fe coords
            continue

        # Check if heme PDB exists for this enzyme
        enz_file = f"ENZ_G{eidx+1:06d}"
        clean_path = CLEAN_DIR / f"{enz_file}.pdb"
        if not clean_path.exists():
            continue

        # Extract Fe coords from PDB
        fe = extract_fe_from_pdb(clean_path)
        if fe:
            row["heme_fe_x"] = f"{fe[0]:.3f}"
            row["heme_fe_y"] = f"{fe[1]:.3f}"
            row["heme_fe_z"] = f"{fe[2]:.3f}"
            updated += 1

    if updated > 0:
        fields = list(rows[0].keys())
        with open(registry_path, "w", encoding="utf-8-sig", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)

    log(f"  Updated {updated} enzymes with Fe coordinates")


def extract_fe_from_pdb(pdb_path):
    """Extract Fe atom coordinates from a PDB file."""
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")):
                atom_name = line[12:16].strip().upper()
                res_name = line[17:20].strip().upper()
                if atom_name.startswith("FE") and res_name in ("HEM", "HEC", "HEA", "HEB", "HEO"):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    return (x, y, z)
    return None


def main():
    log("Processing ColabFold Batch 2 structures\n")

    # Step 1
    enz_map = step1_extract_best_models()

    # Step 2
    step2_run_heme_transplant()

    # Step 3
    step3_cleanup_and_pdbqt()

    # Step 4
    step4_update_registry()

    log("\n=== DONE ===")
    log("Next steps:")
    log("  1. Run batch_receptor_pdbqt.py for PDBQT conversion")
    log("  2. Upload to Cloud-2 for Batch 2 docking")


if __name__ == "__main__":
    main()
