"""
Phase 2a: Dock B6 inhibitor pairs with Vina.

Prepares the ABL-04 "bridge experiment" by docking B6 inhibitor-enzyme pairs
(Label=0 in B6) with AutoDock Vina, then generating structure features.

Pipeline:
  1. Read B6 inhibitor pairs → filter to dockable
  2. Assign new Dock Indices (starting from 4000)
  3. Parallel Vina docking → pocket PDB + ligand SDF
  4. Align ligands (step8_align_ligand.py)
  5. Generate structure_features.lmdb (step8_generate_structure_lmdb.py)

Usage:
    python step6_phase2a_dock_inhibitors.py \
        --b6_data      <data/00_shared/datasets/B6_v1/data.csv> \
        --enzymes_csv  <data/00_shared/datasets/B6_v1/Enzymes.csv> \
        --substrates_csv <data/00_shared/datasets/B6_v1/Substrates.csv> \
        --assets_dir   <data/03_Step3_对接预实验> \
        --pdb_dir      <PathA_.../data/01_Step1_PDB文件> \
        --output_dir   <data/06_Step6_消融实验/phase2_inhibitor_docking> \
        --step8_scripts <scripts/02_Step2_因子实验> \
        [--workers 12] [--resume]
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

DOCK_INDEX_START = 4000
DEFAULT_WORKERS = 12


# ═══════════════════════════════════════════════════════════════════════
#  Helpers (from step4_run_fullscale.py)
# ═══════════════════════════════════════════════════════════════════════

def _find_pdb(pdb_dir: Path, pdb_id: str) -> Optional[Path]:
    for name in (f"{pdb_id}.pdb", f"{pdb_id.upper()}.pdb"):
        p = pdb_dir / name
        if p.exists():
            return p
    return None


def _load_grid(path: Path):
    from lib.grid_locator import GridBox
    d = json.loads(path.read_text(encoding="utf-8"))
    return GridBox(
        float(d["center_x"]), float(d["center_y"]), float(d["center_z"]),
        float(d.get("size_x", 22.5)), float(d.get("size_y", 22.5)),
        float(d.get("size_z", 22.5)))


def run_single_pair(
    dock_index: int, pdb_id: str, enzyme_index: int, substrate_index: int,
    label: int, assets_dir: str, pdb_dir: str, output_dir: str,
    expected_heavy: Optional[int], step3_scripts_dir: str = "",
) -> Dict[str, Any]:
    """Dock one inhibitor pair. Runs in worker process."""
    if step3_scripts_dir and step3_scripts_dir not in sys.path:
        sys.path.insert(0, step3_scripts_dir)
    from lib.postprocess import postprocess_docked_pose
    from lib.vina_runner import run_vina_docking

    assets, out, pdb_root = Path(assets_dir), Path(output_dir), Path(pdb_dir)
    res: Dict[str, Any] = dict(
        dock_index=dock_index, status="failed", affinity=None,
        pocket_atoms=0, fe_distance=None, runtime_sec=0.0, error="")

    rec = assets / "receptors" / "pdbqt" / f"{pdb_id}.pdbqt"
    grd = assets / "receptors" / "grid_boxes" / f"{pdb_id}_grid.json"
    lig = assets / "ligands" / "pdbqt" / f"substrate_{substrate_index}.pdbqt"
    orig = _find_pdb(pdb_root, pdb_id)
    docked = out / "docking_output" / f"{dock_index}.pdbqt"

    for p, tag in [(rec, "receptor"), (grd, "grid"), (lig, "ligand")]:
        if not p.exists():
            res["status"], res["error"] = "skipped", f"Missing {tag}: {p.name}"
            return res
    if orig is None:
        res["status"], res["error"] = "skipped", f"PDB not found: {pdb_id}"
        return res

    try:
        grid = _load_grid(grd)
    except Exception as exc:
        res["error"] = f"Grid parse: {exc}"
        return res

    dock = run_vina_docking(rec, lig, docked, grid)
    res["runtime_sec"] = dock.runtime_sec
    res["affinity"] = dock.best_affinity
    if not dock.success:
        res["error"] = f"Docking: {dock.error}"
        return res

    post = postprocess_docked_pose(
        docked, orig, dock_index, out / "pocket", out / "raw_ligand",
        expected_heavy_atoms=expected_heavy)
    if not post.success:
        res["error"] = f"Postprocess: {post.error}"
        return res

    res.update(status="success", pocket_atoms=post.pocket_atoms,
               fe_distance=post.fe_ligand_distance)
    return res


# ═══════════════════════════════════════════════════════════════════════
#  Main
# ═══════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description="Phase 2a: Dock inhibitor pairs")
    p.add_argument("--b6_data", required=True)
    p.add_argument("--enzymes_csv", required=True)
    p.add_argument("--substrates_csv", required=True)
    p.add_argument("--assets_dir", required=True,
                   help="Step 3 assets (receptors/, ligands/, grid_boxes/)")
    p.add_argument("--pdb_dir", required=True,
                   help="Directory containing original PDB files")
    p.add_argument("--output_dir", required=True)
    p.add_argument("--step8_scripts", required=True,
                   help="Directory containing step8_align_ligand.py etc.")
    p.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    p.add_argument("--resume", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    b6_data = Path(args.b6_data).resolve()
    enzymes_csv = Path(args.enzymes_csv).resolve()
    substrates_csv = Path(args.substrates_csv).resolve()
    assets_dir = Path(args.assets_dir).resolve()
    pdb_dir = Path(args.pdb_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    step8_dir = Path(args.step8_scripts).resolve()

    # Add Step 3 scripts to path (for lib/ imports in workers)
    step3_scripts = assets_dir.parent.parent / "scripts" / "03_Step3_对接环境"
    if not step3_scripts.is_dir():
        step3_scripts = assets_dir.parent / "scripts" / "03_Step3_对接环境"
    if not step3_scripts.is_dir():
        print(f"[ERROR] Cannot find Step 3 scripts directory. Tried:")
        print(f"  {assets_dir.parent.parent / 'scripts' / '03_Step3_对接环境'}")
        print(f"  {assets_dir.parent / 'scripts' / '03_Step3_对接环境'}")
        return 1
    step3_scripts_str = str(step3_scripts)
    sys.path.insert(0, step3_scripts_str)

    # Create output directories
    for subdir in ["docking_output", "pocket", "raw_ligand", "ligand"]:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────
    print("Loading B6 data...")
    b6 = pd.read_csv(b6_data)
    enzymes = pd.read_csv(enzymes_csv)
    substrates = pd.read_csv(substrates_csv)

    enz_to_pdb = dict(zip(enzymes["Enzyme_Index"], enzymes["PDB_ID"]))

    # Receptor/ligand/grid readiness
    rec_summary = pd.read_csv(assets_dir / "receptor_prep_summary.csv")
    lig_summary = pd.read_csv(assets_dir / "ligand_prep_summary.csv")
    rec_ok = set(rec_summary[rec_summary["receptor_ok"] == True]["PDB_ID"])
    lig_ok = set(lig_summary[lig_summary["success"] == True]["Substrate_Index"])
    grid_dir = assets_dir / "receptors" / "grid_boxes"
    grids_ok = set(
        f.name.replace("_grid.json", "") for f in grid_dir.iterdir()
        if f.name.endswith("_grid.json"))

    # Expected heavy atoms per substrate (for QC)
    from rdkit import Chem
    from rdkit.RDLogger import DisableLog
    DisableLog("rdApp.*")
    smiles_map = dict(zip(substrates["Substrate_Index"], substrates["Substrate_SMILES"]))
    heavy_atoms = {}
    for si, smi in smiles_map.items():
        mol = Chem.MolFromSmiles(smi)
        heavy_atoms[si] = mol.GetNumHeavyAtoms() if mol else None

    # ── Build manifest ────────────────────────────────────────────
    manifest_path = output_dir / "manifest.csv"
    progress_path = output_dir / "progress.csv"
    mapping_path = output_dir / "dock_index_mapping.csv"

    if args.resume and manifest_path.exists():
        print("Resuming from existing manifest...")
        manifest = pd.read_csv(manifest_path).to_dict("records")
        if not manifest:
            print("[ERROR] Loaded manifest is empty.")
            return 1
    else:
        b6_neg = b6[b6["Label"] == 0]
        manifest = []
        skipped = 0
        for _, row in b6_neg.iterrows():
            ei = int(row["Enzyme Index"])
            si = int(row["Substrate Index"])
            b6_dock = int(row["Dock Index"])
            pdb_id = enz_to_pdb.get(ei)
            if not pdb_id:
                skipped += 1
                continue
            if pdb_id not in rec_ok or si not in lig_ok or pdb_id not in grids_ok:
                skipped += 1
                continue

            new_dock = DOCK_INDEX_START + len(manifest)
            manifest.append({
                "dock_index": new_dock,
                "b6_dock_index": b6_dock,
                "pdb_id": pdb_id,
                "enzyme_index": ei,
                "substrate_index": si,
                "label": 0,
            })

        print(f"Manifest: {len(manifest)} dockable inhibitor pairs ({skipped} skipped)")

        if not manifest:
            print("[ERROR] No dockable inhibitor pairs found.")
            return 1

        # Save manifest and mapping
        mf = pd.DataFrame(manifest)
        mf.to_csv(manifest_path, index=False)
        mf[["dock_index", "b6_dock_index", "enzyme_index",
            "substrate_index"]].to_csv(mapping_path, index=False)

    # ── Load progress (for resume) ────────────────────────────────
    progress: Dict[int, Dict] = {}
    if args.resume and progress_path.exists():
        prog_df = pd.read_csv(progress_path)
        for _, r in prog_df.iterrows():
            progress[int(r["dock_index"])] = r.to_dict()

    # Filter to pending pairs
    pending = [m for m in manifest if m["dock_index"] not in progress
               or progress[m["dock_index"]].get("status") != "success"]
    already_done = len(manifest) - len(pending)
    print(f"Pending: {len(pending)} pairs ({already_done} already done)")

    if not pending:
        print("All pairs already docked. Skipping to feature generation.")
    else:
        # ── Parallel docking ──────────────────────────────────────
        print(f"\nDocking {len(pending)} inhibitor pairs ({args.workers} workers)...")
        t0 = time.time()

        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {}
            for m in pending:
                fut = pool.submit(
                    run_single_pair,
                    dock_index=m["dock_index"],
                    pdb_id=m["pdb_id"],
                    enzyme_index=m["enzyme_index"],
                    substrate_index=m["substrate_index"],
                    label=0,
                    assets_dir=str(assets_dir),
                    pdb_dir=str(pdb_dir),
                    output_dir=str(output_dir),
                    expected_heavy=heavy_atoms.get(m["substrate_index"]),
                    step3_scripts_dir=step3_scripts_str,
                )
                futures[fut] = m["dock_index"]

            done_count = 0
            success_count = already_done
            for fut in as_completed(futures):
                dock_idx = futures[fut]
                try:
                    result = fut.result()
                except Exception as exc:
                    result = {"dock_index": dock_idx, "status": "error",
                              "error": str(exc), "affinity": None,
                              "pocket_atoms": 0, "fe_distance": None,
                              "runtime_sec": 0.0}
                progress[dock_idx] = result
                done_count += 1
                if result["status"] == "success":
                    success_count += 1

                if done_count % 20 == 0 or done_count == len(pending):
                    elapsed = time.time() - t0
                    print(f"  Progress: {done_count}/{len(pending)} "
                          f"({success_count} total success, {elapsed:.0f}s)")
                    # Checkpoint
                    _save_progress(progress_path, progress)

        wall = time.time() - t0
        _save_progress(progress_path, progress)
        print(f"\nDocking complete: {success_count}/{len(manifest)} success "
              f"({wall:.0f}s)")

    # Count successes (only within current manifest, not stale progress)
    manifest_ids = {m["dock_index"] for m in manifest}
    n_success = sum(1 for di, v in progress.items()
                    if di in manifest_ids and v.get("status") == "success")
    n_total = len(manifest)
    print(f"\nFinal: {n_success}/{n_total} successfully docked")

    if n_success == 0:
        print("[ERROR] No successful dockings. Cannot proceed.")
        return 1

    # ── Build data.csv for inhibitor negatives ────────────────────
    inhibitor_data_csv = output_dir / "inhibitor_data.csv"
    rows = []
    for m in manifest:
        di = m["dock_index"]
        if progress.get(di, {}).get("status") == "success":
            rows.append({
                "Dock Index": di,
                "Enzyme Index": m["enzyme_index"],
                "Substrate Index": m["substrate_index"],
                "Label": 0,
            })
    pd.DataFrame(rows).to_csv(inhibitor_data_csv, index=False)
    print(f"Inhibitor data.csv: {len(rows)} rows → {inhibitor_data_csv}")

    # ── Step 8.2: Align ligands ───────────────────────────────────
    print("\n" + "=" * 50)
    print("STEP 8.2: Align ligands")
    print("=" * 50)

    # Build mapping CSV for alignment
    mapping_for_align = output_dir / "mapping_for_align.csv"
    align_rows = []
    for m in manifest:
        if progress.get(m["dock_index"], {}).get("status") == "success":
            align_rows.append({
                "Dock_Index": m["dock_index"],
                "Substrate_Index": m["substrate_index"],
            })
    pd.DataFrame(align_rows).to_csv(mapping_for_align, index=False)

    align_script = step8_dir / "step8_align_ligand.py"
    if not align_script.exists():
        print(f"[ERROR] Cannot find {align_script}")
        return 1

    python_exe = sys.executable
    align_cmd = [
        python_exe, str(align_script),
        "--raw_ligand_dir", str(output_dir / "raw_ligand"),
        "--aligned_ligand_dir", str(output_dir / "ligand"),
        "--mapping_csv", str(mapping_for_align),
        "--substrates_csv", str(substrates_csv),
    ]
    print(f"Running: {' '.join(align_cmd[:3])}...")
    r = subprocess.run(align_cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(f"[ERROR] Alignment failed:\n{r.stderr[-2000:]}")
        return 1
    print(r.stdout[-1000:] if r.stdout else "(no output)")

    # ── Step 8.3: Generate structure LMDB ─────────────────────────
    print("\n" + "=" * 50)
    print("STEP 8.3: Generate structure features LMDB")
    print("=" * 50)

    lmdb_script = step8_dir / "step8_generate_structure_lmdb.py"
    if not lmdb_script.exists():
        print(f"[ERROR] Cannot find {lmdb_script}")
        return 1

    alignment_summary = output_dir / "alignment_summary.csv"
    if not alignment_summary.exists():
        print(f"[ERROR] alignment_summary.csv not found at {alignment_summary}")
        return 1

    lmdb_cmd = [
        python_exe, str(lmdb_script),
        "--pocket_dir", str(output_dir / "pocket"),
        "--ligand_dir", str(output_dir / "ligand"),
        "--alignment_summary", str(alignment_summary),
        "--output_dir", str(output_dir),
        "--dataset_csv", str(inhibitor_data_csv),
    ]
    print(f"Running: {' '.join(lmdb_cmd[:3])}...")
    r = subprocess.run(lmdb_cmd, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        print(f"[ERROR] LMDB generation failed:\n{r.stderr[-2000:]}")
        return 1
    print(r.stdout[-1000:] if r.stdout else "(no output)")

    # Verify outputs
    lmdb_path = output_dir / "structure_features.lmdb"
    hqid_path = output_dir / "high_quality_id.txt"
    if not lmdb_path.exists():
        print(f"[ERROR] Expected LMDB not found: {lmdb_path}")
        return 1

    if hqid_path.exists():
        hq_ids = [l.strip() for l in hqid_path.read_text().strip().split("\n") if l.strip()]
        print(f"\nHigh quality IDs: {len(hq_ids)}")
    else:
        print("[WARN] high_quality_id.txt not found")

    # ── Summary report ────────────────────────────────────────────
    report_path = output_dir / "phase2a_report.md"
    _write_report(report_path, manifest, progress, n_success,
                  len(hq_ids) if hqid_path.exists() else 0)

    print(f"\n{'=' * 50}")
    print("PHASE 2A COMPLETE")
    print(f"{'=' * 50}")
    print(f"  Docked: {n_success}/{n_total}")
    print(f"  LMDB: {lmdb_path}")
    print(f"  Report: {report_path}")
    return 0


def _save_progress(path: Path, progress: Dict[int, Dict]):
    if not progress:
        return
    rows = sorted(progress.values(), key=lambda r: r.get("dock_index", 0))
    keys = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _write_report(path: Path, manifest, progress, n_success, n_hq):
    lines = [
        "# Phase 2a: Inhibitor Docking Report",
        "",
        f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        "",
        f"## Summary",
        f"- Total inhibitor pairs: {len(manifest)}",
        f"- Successfully docked: {n_success}",
        f"- High quality IDs: {n_hq}",
        "",
        "## Status Breakdown",
    ]
    status_counts: Dict[str, int] = {}
    for v in progress.values():
        s = v.get("status", "unknown")
        status_counts[s] = status_counts.get(s, 0) + 1
    for s, c in sorted(status_counts.items()):
        lines.append(f"- {s}: {c}")

    # Affinity stats for successful dockings
    affinities = [v["affinity"] for v in progress.values()
                  if v.get("status") == "success" and v.get("affinity") is not None]
    if affinities:
        lines.extend([
            "",
            "## Docking Affinity (kcal/mol)",
            f"- Mean: {np.mean(affinities):.2f}",
            f"- Min: {np.min(affinities):.2f}",
            f"- Max: {np.max(affinities):.2f}",
            f"- Median: {np.median(affinities):.2f}",
        ])

    lines.append("")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"\n[ERROR] {exc}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
