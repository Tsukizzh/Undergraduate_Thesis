"""
Phase 6 Stage C: Batch receptor PDB → PDBQT using MGLTools prepare_receptor4.py
- Processes all 1,403 Batch 1 enzymes
- Uses MGLTools (Python 2.7) with temp English-path workaround
- Multi-worker parallel processing
- Resume support via .ok markers

Proven in Path B (282/292 = 96.6% success rate).
"""

import pandas as pd
import subprocess
import shutil
import time
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE = Path(__file__).parent.parent.parent
DATA = BASE / "data"
REGISTRIES = DATA / "registries"
STRUCTS = DATA / "structures"
MANIFESTS = STRUCTS / "manifests"
OUT_CLEAN = STRUCTS / "receptors_clean"
OUT_PDBQT = STRUCTS / "receptors_pdbqt"
OUT_CLEAN.mkdir(parents=True, exist_ok=True)
OUT_PDBQT.mkdir(parents=True, exist_ok=True)

MGLTOOLS_PYTHON = r"D:\autodock\MGLTools-1.5.7\python.exe"
PREPARE_RECEPTOR = r"D:\autodock\MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py"
TMP_DIR = Path(r"D:\autodock\tmp")
TMP_DIR.mkdir(parents=True, exist_ok=True)

MAX_WORKERS = 6  # MGLTools is CPU-bound, don't over-parallelize

# Import cleanup functions from pilot script
import sys
sys.path.insert(0, str(Path(__file__).parent))
from pilot_receptor_prep import cleanup_structure, normalize_pdb

# Source type → file path resolution
def resolve_source_path(row):
    """Resolve actual file path for each enzyme based on source type."""
    status = row["status"]
    eid = row["global_enzyme_id"]

    if status == "heme_transplanted":
        pdb_path = row.get("heme_transplant_pdb_path", "")
        if pd.notna(pdb_path) and pdb_path:
            return "v2_transplant", STRUCTS / "heme_transplant" / "pdb" / Path(str(pdb_path)).name
    elif status == "downloaded_with_heme":
        cif_path = row.get("alphafill_cif_path", "")
        if pd.notna(cif_path) and cif_path:
            return "alphafill_cif", STRUCTS / "alphafill" / "cif" / Path(str(cif_path)).name
    elif status == "pcpd_predicted":
        pdb_ref = row.get("existing_pdb_path", "")
        if pd.notna(pdb_ref) and pdb_ref:
            name = Path(str(pdb_ref)).name
            return "pcpd_pdb", BASE / "downloads" / "PCPD" / "PDB" / name
    elif status == "experimental_pdb":
        pdb_ref = row.get("existing_pdb_path", "")
        if pd.notna(pdb_ref) and pdb_ref:
            pdb_id = str(pdb_ref).replace("s1_rcsb/", "")
            for candidate in [
                BASE / "downloads" / "RCSB" / "PDB" / f"{pdb_id}.pdb",
                BASE.parent.parent / "PathA_2026-01-08_模型评估测试集构建" / "data" / "01_Step1_PDB文件" / f"{pdb_id}.pdb",
            ]:
                if candidate.exists():
                    return "experimental_pdb", candidate
            return "experimental_pdb", BASE / "downloads" / "RCSB" / "PDB" / f"{pdb_id}.pdb"

    return None, None


def process_one_receptor(args):
    """Process one enzyme: PDB cleanup → normalize → MGLTools → PDBQT."""
    enzyme_index, row_dict = args
    eid = row_dict["global_enzyme_id"]

    # Check skip
    ok_marker = OUT_PDBQT / f"{eid}.ok"
    if ok_marker.exists():
        return eid, "skipped", ""

    # Resolve source path
    source_type, source_path = resolve_source_path(pd.Series(row_dict))
    if source_type is None or source_path is None:
        return eid, "no_source_path", ""
    if not source_path.exists():
        return eid, "file_not_found", str(source_path)

    # Step 1: PDB cleanup
    try:
        result = cleanup_structure(source_path, source_type, row_dict, eid)
        if result["status"] != "ok":
            return eid, f"cleanup_failed({result['status']})", ""
    except Exception as e:
        return eid, f"cleanup_exception({type(e).__name__})", str(e)[:200]

    clean_pdb = Path(result["clean_pdb"])

    # Step 2: MGLTools in temp English-path directory
    try:
        import tempfile
        with tempfile.TemporaryDirectory(dir=str(TMP_DIR)) as tmpdir:
            tmp_pdb = Path(tmpdir) / f"{eid}.pdb"
            tmp_pdbqt = Path(tmpdir) / f"{eid}.pdbqt"

            # Fix altloc column (col 17 → space) before MGLTools
            with open(str(clean_pdb)) as f:
                lines = f.readlines()
            fixed = []
            for l in lines:
                if l.startswith(("ATOM", "HETATM")) and len(l) > 16:
                    l = l[:16] + " " + l[17:]
                fixed.append(l)
            with open(str(tmp_pdb), "w") as f:
                f.writelines(fixed)

            cmd = [MGLTOOLS_PYTHON, PREPARE_RECEPTOR,
                   "-r", str(tmp_pdb), "-o", str(tmp_pdbqt),
                   "-A", "hydrogens", "-U", "nphs_lps_waters"]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if not tmp_pdbqt.exists() or tmp_pdbqt.stat().st_size == 0:
                return eid, "mgltools_failed", proc.stderr[:300]

            # Validate Fe in PDBQT
            content = tmp_pdbqt.read_text()
            fe_lines = [l for l in content.split("\n") if "HEM" in l[17:20] and "FE" in l[12:16].upper()]
            if not fe_lines:
                return eid, "no_fe_in_pdbqt", ""

            # Copy to output
            shutil.copy2(str(tmp_pdbqt), str(OUT_PDBQT / f"{eid}.pdbqt"))

        # Write success marker
        ok_marker.write_text("ok")
        return eid, "ok", ""

    except subprocess.TimeoutExpired:
        return eid, "timeout", ""
    except Exception as e:
        return eid, f"mgltools_exception({type(e).__name__})", str(e)[:200]


def main():
    print("=" * 60)
    print("Phase 6 Stage C: Batch Receptor PDB → PDBQT")
    print("=" * 60)

    manifest = pd.read_csv(MANIFESTS / "receptor_manifest.csv")
    enz_reg = pd.read_csv(REGISTRIES / "enzyme_registry.csv")

    # Only Batch 1 enzymes
    batch1 = enz_reg[enz_reg["enzyme_batch"] == "batch_1"]
    print(f"Batch 1 enzymes: {len(batch1)}")

    # Merge with manifest for file paths
    merged = batch1.merge(manifest, on="global_enzyme_id", how="left", suffixes=("", "_manifest"))

    # Build task list
    tasks = []
    for _, row in merged.iterrows():
        tasks.append((int(row["enzyme_index"]), row.to_dict()))

    # Count already done
    already_done = sum(1 for eid in merged["global_enzyme_id"] if (OUT_PDBQT / f"{eid}.ok").exists())
    remaining = len(tasks) - already_done
    print(f"Already done: {already_done}, Remaining: {remaining}")

    if remaining == 0:
        print("All receptors already converted!")
        return

    # Process
    print(f"Processing {remaining} receptors with {MAX_WORKERS} workers...")
    t0 = time.time()
    results = {"ok": 0, "skipped": 0, "failed": 0}
    errors = []

    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(process_one_receptor, t): t[1]["global_enzyme_id"] for t in tasks}
        for i, future in enumerate(as_completed(futures), 1):
            eid, status, detail = future.result()
            if status == "ok":
                results["ok"] += 1
            elif status == "skipped":
                results["skipped"] += 1
            else:
                results["failed"] += 1
                errors.append({"enzyme": eid, "status": status, "detail": detail[:100]})

            if i % 50 == 0 or i == len(tasks):
                elapsed = time.time() - t0
                rate = i / elapsed if elapsed > 0 else 0
                eta = (len(tasks) - i) / rate if rate > 0 else 0
                print(f"  [{i}/{len(tasks)}] ok={results['ok']} skip={results['skipped']} fail={results['failed']} ({elapsed:.0f}s, ETA {eta:.0f}s)")

    elapsed = time.time() - t0
    print(f"\n{'='*60}")
    print(f"SUMMARY")
    print(f"{'='*60}")
    print(f"OK: {results['ok']}")
    print(f"Skipped: {results['skipped']}")
    print(f"Failed: {results['failed']}")
    print(f"Time: {elapsed:.1f}s ({elapsed/60:.1f}min)")

    if errors:
        print(f"\nFailed receptors ({len(errors)}):")
        for e in errors[:30]:
            print(f"  {e['enzyme']}: {e['status']} | {e['detail']}")


if __name__ == "__main__":
    main()
