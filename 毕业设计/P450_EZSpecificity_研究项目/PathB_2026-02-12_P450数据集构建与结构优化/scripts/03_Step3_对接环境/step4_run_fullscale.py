"""Step 4: Full-scale docking — all dockable positive pairs + random negatives.

Builds a manifest of enzyme-substrate pairs by:
1. Filtering B6 positives to only dockable receptor+ligand combinations
2. Generating random negatives from dockable enzyme pool (ratio ~9.36)
3. Running AutoDock Vina docking + postprocessing in parallel
4. Producing data.csv compatible with the model pipeline

Supports --resume to continue from checkpoint.

Usage:
    python step4_run_fullscale.py \\
        --data_csv      .../B6_v1/data.csv \\
        --enzymes_csv   .../B6_v1/Enzymes.csv \\
        --substrates_csv .../B6_v1/Substrates.csv \\
        --assets_dir    .../data/03_Step3_对接预实验/ \\
        --pdb_dir       .../PathA_.../data/01_Step1_PDB文件/ \\
        --output_dir    .../data/04_Step4_批量对接/ \\
        --results_dir   .../results/04_Step4_批量对接/ \\
        --num_workers 12
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import shutil
import statistics
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

LOG = logging.getLogger(__name__)

DOCK_INDEX_NEG_START = 1000
DEFAULT_NEG_RATIO = 9.36
MANIFEST_FIELDS = [
    "dock_index", "pdb_id", "enzyme_index", "substrate_index", "label", "source",
]
PROGRESS_FIELDS = [
    "dock_index", "status", "affinity", "pocket_atoms",
    "fe_distance", "runtime_sec", "error",
]


# ── helpers ──────────────────────────────────────────────────────────────────

def _nk(k: str) -> str:
    """Normalize CSV column key: strip, lowercase, spaces → underscores."""
    return k.strip().lower().replace(" ", "_")


def _csv_rows(path: Path) -> list[dict[str, str]]:
    """Read CSV returning rows with normalized keys."""
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return [{_nk(k): (v.strip() if v else "")
                 for k, v in row.items()} for row in csv.DictReader(f)]


def _safe_int(v: str) -> int:
    """Parse int from string, tolerating float-encoded ints like '12.0'."""
    return int(float(v)) if v else 0


def _opt_float(v: str) -> Optional[float]:
    return float(v) if v else None


# ── asset availability ───────────────────────────────────────────────────────

def load_dockable_pdb_ids(summary: Path) -> set[str]:
    """PDB IDs with receptor_ok == True in receptor_prep_summary.csv."""
    return {r["pdb_id"].upper() for r in _csv_rows(summary)
            if r.get("receptor_ok", "").lower() == "true"}


def load_dockable_substrates(summary: Path) -> set[int]:
    """Substrate indices with success == True in ligand_prep_summary.csv."""
    return {_safe_int(r["substrate_index"]) for r in _csv_rows(summary)
            if r.get("success", "").lower() == "true"}


# ── data loading ─────────────────────────────────────────────────────────────

def load_enzyme_map(path: Path) -> dict[int, str]:
    return {int(float(r["enzyme_index"])): r.get("pdb_id", "").upper()
            for r in _csv_rows(path)}


def load_all_positives(path: Path) -> list[tuple[int, int, int]]:
    """All (dock_index, enzyme_index, substrate_index) with Label==1."""
    return [(_safe_int(r["dock_index"]), _safe_int(r["enzyme_index"]),
             _safe_int(r["substrate_index"]))
            for r in _csv_rows(path) if _safe_int(r["label"]) == 1]


def load_substrate_heavy_atoms(path: Path) -> dict[int, int]:
    from rdkit import Chem
    out: dict[int, int] = {}
    for r in _csv_rows(path):
        mol = Chem.MolFromSmiles(r.get("substrate_smiles", ""))
        if mol:
            out[int(float(r["substrate_index"]))] = sum(
                1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    return out


# ── manifest ─────────────────────────────────────────────────────────────────

def build_manifest(
    data_csv: Path, enzymes_csv: Path,
    rec_summary: Path, lig_summary: Path,
    output_dir: Path, neg_ratio: float, seed: int,
) -> list[dict[str, Any]]:
    """Create manifest: filtered positives + asset-aware random negatives."""
    ok_pdbs = load_dockable_pdb_ids(rec_summary)
    ok_subs = load_dockable_substrates(lig_summary)
    enz_map = load_enzyme_map(enzymes_csv)
    positives = load_all_positives(data_csv)
    LOG.info("Dockable: %d receptors, %d substrates", len(ok_pdbs), len(ok_subs))

    # 1. Filter positives to dockable only
    tasks: list[dict[str, Any]] = []
    skip = 0
    for di, ei, si in positives:
        pdb = enz_map.get(ei, "")
        if pdb not in ok_pdbs or si not in ok_subs:
            skip += 1
            continue
        tasks.append(dict(dock_index=di, pdb_id=pdb, enzyme_index=ei,
                          substrate_index=si, label=1, source="positive"))
    n_pos = len(tasks)
    LOG.info("Positives: %d dockable, %d skipped", n_pos, skip)

    # 2. Prepare filtered CSVs for negative sampler (asset-aware)
    filt_pos = output_dir / "_filtered_positives.csv"
    with filt_pos.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Dock Index", "Enzyme Index", "Substrate Index", "Label"])
        for t in tasks:
            w.writerow([t["dock_index"], t["enzyme_index"], t["substrate_index"], 1])

    filt_enz = output_dir / "_filtered_enzymes.csv"
    with enzymes_csv.open("r", encoding="utf-8-sig", newline="") as fin:
        reader = csv.DictReader(fin)
        with filt_enz.open("w", encoding="utf-8", newline="") as fout:
            writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                pdb = {_nk(k): v.strip() for k, v in row.items()}.get("pdb_id", "").upper()
                if pdb in ok_pdbs:
                    writer.writerow(row)

    # 3. Generate negatives (only from dockable enzymes)
    from lib.negative_sampler import generate_negative_pairs

    neg_csv = output_dir / "_negatives.csv"
    nr = generate_negative_pairs(
        filt_pos, filt_enz, neg_csv,
        ratio=neg_ratio, random_seed=seed, dock_index_start=DOCK_INDEX_NEG_START,
    )
    if not nr.success:
        LOG.error("Negative sampling failed: %s", nr.error)
        sys.exit(1)
    LOG.info("Negatives generated: %d (ratio %.2f)", nr.negative_count, nr.achieved_ratio)

    # 4. Load negatives into manifest
    for r in _csv_rows(neg_csv):
        tasks.append(dict(
            dock_index=_safe_int(r["dock_index"]), pdb_id=r.get("pdb_id", "").upper(),
            enzyme_index=_safe_int(r["enzyme_index"]),
            substrate_index=_safe_int(r["substrate_index"]),
            label=0, source="negative"))

    n_neg = len(tasks) - n_pos
    LOG.info("Manifest: %d total (%d pos + %d neg, ratio %.2f)",
             len(tasks), n_pos, n_neg, n_neg / n_pos if n_pos else 0)

    # Verify dock_index uniqueness (positives 0-515, negatives start at 1000)
    all_di = [t["dock_index"] for t in tasks]
    if len(all_di) != len(set(all_di)):
        dupes = [di for di in all_di if all_di.count(di) > 1]
        LOG.error("dock_index collision detected: %s", sorted(set(dupes))[:10])
        sys.exit(1)

    # 5. Save manifest (immutable after creation)
    mpath = output_dir / "manifest.csv"
    with mpath.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        w.writeheader()
        w.writerows(tasks)
    return tasks


def load_manifest(output_dir: Path) -> list[dict[str, Any]]:
    return [dict(dock_index=_safe_int(r["dock_index"]), pdb_id=r["pdb_id"],
                 enzyme_index=_safe_int(r["enzyme_index"]),
                 substrate_index=_safe_int(r["substrate_index"]),
                 label=_safe_int(r["label"]), source=r["source"])
            for r in _csv_rows(output_dir / "manifest.csv")]


# ── progress tracking ────────────────────────────────────────────────────────

def load_progress(output_dir: Path) -> dict[int, dict[str, Any]]:
    path = output_dir / "progress.csv"
    if not path.exists():
        return {}
    out: dict[int, dict[str, Any]] = {}
    for r in _csv_rows(path):
        try:
            di = _safe_int(r["dock_index"])
            out[di] = dict(
                dock_index=di, status=r.get("status", "unknown"),
                affinity=_opt_float(r.get("affinity", "")),
                pocket_atoms=_safe_int(r.get("pocket_atoms", "0")),
                fe_distance=_opt_float(r.get("fe_distance", "")),
                runtime_sec=float(r.get("runtime_sec", 0) or 0),
                error=r.get("error", ""))
        except (ValueError, KeyError) as exc:
            LOG.warning("Skipping malformed progress row: %s", exc)
    return out


def save_progress(output_dir: Path, progress: dict[int, dict]) -> None:
    """Atomic write: write to temp file first, then replace target."""
    path = output_dir / "progress.csv"
    items = sorted(progress.values(), key=lambda x: x["dock_index"])
    fd, tmp = tempfile.mkstemp(dir=str(output_dir), suffix=".tmp", prefix="progress_")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=PROGRESS_FIELDS)
            w.writeheader()
            for it in items:
                aff = it.get("affinity")
                fed = it.get("fe_distance")
                w.writerow(dict(
                    dock_index=it["dock_index"], status=it["status"],
                    affinity=f"{aff:.3f}" if aff is not None else "",
                    pocket_atoms=it.get("pocket_atoms", 0),
                    fe_distance=f"{fed:.2f}" if fed is not None else "",
                    runtime_sec=f"{it.get('runtime_sec', 0):.1f}",
                    error=it.get("error", "")))
        shutil.move(tmp, str(path))
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def validate_outputs_on_resume(output_dir: Path, progress: dict[int, dict]) -> int:
    """Demote 'success' entries whose output files are missing."""
    demoted = 0
    for di, rec in progress.items():
        if rec["status"] != "success":
            continue
        pkt = output_dir / "pocket" / f"{di}.pdb"
        lig = output_dir / "raw_ligand" / f"{di}.sdf"
        if not pkt.exists() or not lig.exists():
            rec["status"] = "pending"
            rec["error"] = "output files missing on resume"
            demoted += 1
    return demoted


# ── single-pair docking ─────────────────────────────────────────────────────

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
    label: int, source: str,
    assets_dir: str, pdb_dir: str, output_dir: str,
    expected_heavy: Optional[int],
) -> dict[str, Any]:
    """Execute one docking + postprocess. Runs in worker process."""
    from lib.postprocess import postprocess_docked_pose
    from lib.vina_runner import run_vina_docking

    assets, out, pdb_root = Path(assets_dir), Path(output_dir), Path(pdb_dir)
    res: dict[str, Any] = dict(
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


# ── final outputs ────────────────────────────────────────────────────────────

def write_dataset(
    manifest: list[dict], progress: dict[int, dict],
    enzymes_csv: Path, substrates_csv: Path, results_dir: Path,
) -> int:
    """Write model-compatible data.csv with space-separated headers."""
    results_dir.mkdir(parents=True, exist_ok=True)
    rows = sorted(
        [{"Dock Index": t["dock_index"], "Enzyme Index": t["enzyme_index"],
          "Substrate Index": t["substrate_index"], "Label": t["label"]}
         for t in manifest
         if progress.get(t["dock_index"], {}).get("status") == "success"],
        key=lambda x: x["Dock Index"])

    with (results_dir / "data.csv").open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "Dock Index", "Enzyme Index", "Substrate Index", "Label"])
        w.writeheader()
        w.writerows(rows)

    for src in (enzymes_csv, substrates_csv):
        dest = results_dir / src.name
        if not dest.exists():
            shutil.copy2(str(src), str(dest))
    return len(rows)


def write_report(
    path: Path, manifest: list[dict], progress: dict[int, dict],
    wall_sec: float, workers: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    total = len(manifest)
    n_pos = sum(1 for t in manifest if t["label"] == 1)
    n_neg = total - n_pos
    cts = Counter(
        progress.get(t["dock_index"], {}).get("status", "missing")
        for t in manifest)
    ok = cts.get("success", 0)

    affs = [progress[t["dock_index"]]["affinity"]
            for t in manifest
            if progress.get(t["dock_index"], {}).get("status") == "success"
            and progress[t["dock_index"]].get("affinity") is not None]
    times = [progress[t["dock_index"]]["runtime_sec"]
             for t in manifest if t["dock_index"] in progress
             and progress[t["dock_index"]].get("runtime_sec")]
    pos_ok = sum(1 for t in manifest if t["label"] == 1
                 and progress.get(t["dock_index"], {}).get("status") == "success")
    neg_ok = ok - pos_ok

    ln = ["# Step 4 Full-Scale Docking Report", "",
          f"- **Total pairs**: {total} ({n_pos} positive, {n_neg} negative)"]
    if total:
        ln.append(f"- **Success**: {ok}/{total} ({ok / total * 100:.1f}%)")
    ln.append(f"- **Workers**: {workers}, **Wall time**: {wall_sec:.1f}s ({wall_sec / 3600:.1f}h)")
    if affs:
        ln.append(f"- **Affinity**: mean={statistics.mean(affs):.3f}, "
                   f"min={min(affs):.3f}, max={max(affs):.3f}")
    if times:
        ln.append(f"- **Per-pair time**: mean={statistics.mean(times):.1f}s, "
                   f"max={max(times):.1f}s")
        ln.append(f"- **Estimated serial time**: {sum(times) / 3600:.1f}h")

    ln += ["", "## Status Breakdown", ""]
    for s, n in sorted(cts.items()):
        ln.append(f"- {s}: {n}")

    ln += ["", "## Success by Label", ""]
    if n_pos:
        ln.append(f"- Positive: {pos_ok}/{n_pos} ({pos_ok / n_pos * 100:.1f}%)")
    if n_neg:
        ln.append(f"- Negative: {neg_ok}/{n_neg} ({neg_ok / n_neg * 100:.1f}%)")
    if pos_ok:
        ln.append(f"- Final ratio: 1:{neg_ok / pos_ok:.2f}")

    fails = [(t, progress.get(t["dock_index"], {}))
             for t in manifest
             if progress.get(t["dock_index"], {}).get("status")
             not in ("success", None, "missing")]
    if fails:
        ln += ["", f"## Failures ({len(fails)} total, showing first 50)", "",
               "| Dock_Index | PDB_ID | Label | Status | Error |",
               "|---:|---|---:|---|---|"]
        for t, r in sorted(fails, key=lambda x: x[0]["dock_index"])[:50]:
            ln.append(f"| {t['dock_index']} | {t['pdb_id']} | {t['label']} | "
                      f"{r.get('status', '?')} | {r.get('error', '')} |")
        if len(fails) > 50:
            ln.append(f"\n*...{len(fails) - 50} more (see progress.csv)*")

    path.write_text("\n".join(ln), encoding="utf-8")


# ── CLI ──────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data_csv", type=Path, required=True,
                   help="B6 data.csv")
    p.add_argument("--enzymes_csv", type=Path, required=True,
                   help="B6 Enzymes.csv")
    p.add_argument("--substrates_csv", type=Path, required=True,
                   help="B6 Substrates.csv")
    p.add_argument("--assets_dir", type=Path, required=True,
                   help="Step 3 assets (receptors/, ligands/)")
    p.add_argument("--pdb_dir", type=Path, required=True,
                   help="Raw PDB directory")
    p.add_argument("--output_dir", type=Path, required=True,
                   help="data/04_Step4_批量对接/")
    p.add_argument("--results_dir", type=Path, required=True,
                   help="results/04_Step4_批量对接/")
    p.add_argument("--receptor_summary", type=Path, default=None,
                   help="Default: {assets_dir}/receptor_prep_summary.csv")
    p.add_argument("--ligand_summary", type=Path, default=None,
                   help="Default: {assets_dir}/ligand_prep_summary.csv")
    p.add_argument("--num_workers", type=int, default=12)
    p.add_argument("--random_seed", type=int, default=2026)
    p.add_argument("--neg_ratio", type=float, default=DEFAULT_NEG_RATIO)
    p.add_argument("--resume", action="store_true",
                   help="Resume from existing manifest/progress")
    p.add_argument("--checkpoint_every", type=int, default=50,
                   help="Save progress every N completed pairs")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S")

    if args.checkpoint_every < 1:
        LOG.error("--checkpoint_every must be >= 1")
        return 2

    if args.receptor_summary is None:
        args.receptor_summary = args.assets_dir / "receptor_prep_summary.csv"
    if args.ligand_summary is None:
        args.ligand_summary = args.assets_dir / "ligand_prep_summary.csv"

    for tag, p in [("data_csv", args.data_csv), ("enzymes_csv", args.enzymes_csv),
                    ("substrates_csv", args.substrates_csv),
                    ("assets_dir", args.assets_dir), ("pdb_dir", args.pdb_dir),
                    ("receptor_summary", args.receptor_summary),
                    ("ligand_summary", args.ligand_summary)]:
        if not p.exists():
            LOG.error("Not found: %s (%s)", p, tag)
            return 2

    for d in (args.output_dir / "docking_output", args.output_dir / "pocket",
              args.output_dir / "raw_ligand", args.results_dir):
        d.mkdir(parents=True, exist_ok=True)

    # ── manifest ──
    manifest_exists = (args.output_dir / "manifest.csv").exists()
    if args.resume and manifest_exists:
        LOG.info("Resume: loading existing manifest")
        manifest = load_manifest(args.output_dir)
    else:
        if args.resume:
            LOG.warning("--resume but no manifest found; building fresh")
        manifest = build_manifest(
            args.data_csv, args.enzymes_csv,
            args.receptor_summary, args.ligand_summary,
            args.output_dir, args.neg_ratio, args.random_seed)

    # ── progress ──
    progress: dict[int, dict] = {}
    if args.resume:
        progress = load_progress(args.output_dir)
        dem = validate_outputs_on_resume(args.output_dir, progress)
        if dem:
            LOG.warning("Demoted %d entries with missing outputs", dem)
        done = sum(1 for v in progress.values() if v["status"] == "success")
        LOG.info("Resume: %d/%d already success", done, len(manifest))

    # ── docking ──
    pending = [t for t in manifest
               if progress.get(t["dock_index"], {}).get("status") != "success"]
    wall_sec = 0.0

    if not pending:
        LOG.info("All %d pairs already completed", len(manifest))
    else:
        LOG.info("Pending: %d pairs (%d workers)", len(pending), args.num_workers)
        heavy = load_substrate_heavy_atoms(args.substrates_csv)
        t0 = time.perf_counter()
        done_count = 0

        with ProcessPoolExecutor(max_workers=args.num_workers) as pool:
            futs = {
                pool.submit(
                    run_single_pair,
                    t["dock_index"], t["pdb_id"], t["enzyme_index"],
                    t["substrate_index"], t["label"], t["source"],
                    str(args.assets_dir), str(args.pdb_dir),
                    str(args.output_dir), heavy.get(t["substrate_index"]),
                ): t["dock_index"]
                for t in pending
            }

            for fut in as_completed(futs):
                di = futs[fut]
                try:
                    r = fut.result()
                except Exception as exc:
                    r = dict(dock_index=di, status="exception", affinity=None,
                             pocket_atoms=0, fe_distance=None,
                             runtime_sec=0, error=str(exc))

                progress[di] = {k: r.get(k) for k in PROGRESS_FIELDS}
                done_count += 1

                if done_count % 10 == 0 or done_count == len(pending):
                    elapsed = time.perf_counter() - t0
                    ok_n = sum(1 for v in progress.values()
                               if v["status"] == "success")
                    LOG.info("[%d/%d] %.0fs elapsed, %d success",
                             done_count, len(pending), elapsed, ok_n)

                if done_count % args.checkpoint_every == 0:
                    save_progress(args.output_dir, progress)

        wall_sec = time.perf_counter() - t0
        save_progress(args.output_dir, progress)
        LOG.info("Docking complete: %.1fs (%.1fh)", wall_sec, wall_sec / 3600)

    # ── outputs ──
    n_ds = write_dataset(manifest, progress, args.enzymes_csv,
                         args.substrates_csv, args.results_dir)
    report_path = args.results_dir / "fullscale_report.md"
    write_report(report_path, manifest, progress, wall_sec, args.num_workers)

    ok = sum(1 for v in progress.values() if v["status"] == "success")
    pct = ok / len(manifest) * 100 if manifest else 0
    print(f"\n{'=' * 50}")
    print(f"Step 4 Complete")
    print(f"  Manifest: {len(manifest)} pairs")
    print(f"  Success:  {ok} ({pct:.1f}%)")
    print(f"  Dataset:  {n_ds} rows -> {args.results_dir / 'data.csv'}")
    print(f"  Report:   {report_path}")
    print(f"{'=' * 50}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
