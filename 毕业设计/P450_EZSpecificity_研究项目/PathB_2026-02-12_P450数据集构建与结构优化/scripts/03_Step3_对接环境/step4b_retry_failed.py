"""Step 4b: Retry failed docking pairs with adjusted parameters.

Three-tier retry strategy:
  Tier 0: FE distance QC failures → re-postprocess with lower FE_DISTANCE_MIN (no re-dock)
  Tier 1: Timeout failures → re-dock with timeout=600s, exhaustiveness=8
  Tier 2: Remaining timeouts → re-dock with timeout=600s, exhaustiveness=4

Usage:
    python step4b_retry_failed.py \
        --assets_dir    .../data/03_Step3_对接预实验/ \
        --pdb_dir       .../PathA_.../data/01_Step1_PDB文件/ \
        --output_dir    .../data/04_Step4_批量对接/ \
        --results_dir   .../results/04_Step4_批量对接/ \
        --enzymes_csv   .../B6_v1/Enzymes.csv \
        --substrates_csv .../B6_v1/Substrates.csv \
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

PROGRESS_FIELDS = [
    "dock_index", "status", "affinity", "pocket_atoms",
    "fe_distance", "runtime_sec", "error",
]

# ── helpers ──────────────────────────────────────────────────────────────────

def _nk(k: str) -> str:
    return k.strip().lower().replace(" ", "_")


def _csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return [{_nk(k): (v.strip() if v else "")
                 for k, v in row.items()} for row in csv.DictReader(f)]


def _safe_int(v: str) -> int:
    return int(float(v)) if v else 0


def _opt_float(v: str) -> Optional[float]:
    return float(v) if v else None


# ── progress I/O ─────────────────────────────────────────────────────────────

def load_progress(output_dir: Path) -> dict[int, dict[str, Any]]:
    path = output_dir / "progress.csv"
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
        os.replace(tmp, str(path))
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def load_manifest(output_dir: Path) -> list[dict[str, Any]]:
    return [dict(dock_index=_safe_int(r["dock_index"]), pdb_id=r["pdb_id"],
                 enzyme_index=_safe_int(r["enzyme_index"]),
                 substrate_index=_safe_int(r["substrate_index"]),
                 label=_safe_int(r["label"]), source=r["source"])
            for r in _csv_rows(output_dir / "manifest.csv")]


# ── Tier 0: FE distance re-postprocess ───────────────────────────────────────

def repostprocess_fe_failure(
    dock_index: int, pdb_id: str, substrate_index: int,
    output_dir: str, pdb_dir: str,
    fe_dist_min: float,
    expected_heavy: Optional[int],
) -> dict[str, Any]:
    """Re-run postprocess only (docked PDBQT already exists) with relaxed FE threshold."""
    from lib.postprocess import postprocess_docked_pose

    out, pdb_root = Path(output_dir), Path(pdb_dir)
    docked = out / "docking_output" / f"{dock_index}.pdbqt"

    # Find original PDB
    orig = None
    for name in (f"{pdb_id}.pdb", f"{pdb_id.upper()}.pdb"):
        p = pdb_root / name
        if p.exists():
            orig = p
            break

    res: dict[str, Any] = dict(
        dock_index=dock_index, status="failed", affinity=None,
        pocket_atoms=0, fe_distance=None, runtime_sec=0.0, error="")

    if not docked.exists():
        res["error"] = f"Docked PDBQT missing: {docked.name}"
        return res
    if orig is None:
        res["error"] = f"Original PDB not found: {pdb_id}"
        return res

    post = postprocess_docked_pose(
        docked, orig, dock_index, out / "pocket", out / "raw_ligand",
        expected_heavy_atoms=expected_heavy,
        fe_dist_min=fe_dist_min)

    if not post.success:
        res["error"] = f"Postprocess: {post.error}"
        res["fe_distance"] = post.fe_ligand_distance
        return res

    res.update(status="success", pocket_atoms=post.pocket_atoms,
               fe_distance=post.fe_ligand_distance)
    return res


# ── Tier 1/2: re-dock with adjusted parameters ──────────────────────────────

def _load_grid(path: Path):
    from lib.grid_locator import GridBox
    d = json.loads(path.read_text(encoding="utf-8"))
    return GridBox(
        float(d["center_x"]), float(d["center_y"]), float(d["center_z"]),
        float(d.get("size_x", 22.5)), float(d.get("size_y", 22.5)),
        float(d.get("size_z", 22.5)))


def redock_pair(
    dock_index: int, pdb_id: str, enzyme_index: int, substrate_index: int,
    assets_dir: str, pdb_dir: str, output_dir: str,
    expected_heavy: Optional[int],
    timeout_sec: int, exhaustiveness: int,
    fe_dist_min: float,
    tier_tag: str,
) -> dict[str, Any]:
    """Re-dock a single pair with adjusted Vina parameters."""
    from lib.postprocess import postprocess_docked_pose
    from lib.vina_runner import run_vina_docking

    assets, out, pdb_root = Path(assets_dir), Path(output_dir), Path(pdb_dir)
    res: dict[str, Any] = dict(
        dock_index=dock_index, status="failed", affinity=None,
        pocket_atoms=0, fe_distance=None, runtime_sec=0.0, error="")

    rec = assets / "receptors" / "pdbqt" / f"{pdb_id}.pdbqt"
    grd = assets / "receptors" / "grid_boxes" / f"{pdb_id}_grid.json"
    lig = assets / "ligands" / "pdbqt" / f"substrate_{substrate_index}.pdbqt"
    docked = out / "docking_output" / f"{dock_index}.pdbqt"

    orig = None
    for name in (f"{pdb_id}.pdb", f"{pdb_id.upper()}.pdb"):
        p = pdb_root / name
        if p.exists():
            orig = p
            break

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

    dock = run_vina_docking(
        rec, lig, docked, grid,
        exhaustiveness=exhaustiveness,
        timeout_sec=timeout_sec)
    res["runtime_sec"] = dock.runtime_sec
    res["affinity"] = dock.best_affinity
    if not dock.success:
        res["error"] = f"Docking: {dock.error} [{tier_tag}]"
        return res

    post = postprocess_docked_pose(
        docked, orig, dock_index, out / "pocket", out / "raw_ligand",
        expected_heavy_atoms=expected_heavy,
        fe_dist_min=fe_dist_min)
    if not post.success:
        res["error"] = f"Postprocess: {post.error} [{tier_tag}]"
        res["fe_distance"] = post.fe_ligand_distance
        return res

    res.update(status="success", pocket_atoms=post.pocket_atoms,
               fe_distance=post.fe_ligand_distance)
    return res


# ── dataset regeneration ─────────────────────────────────────────────────────

def write_dataset(
    manifest: list[dict], progress: dict[int, dict],
    enzymes_csv: Path, substrates_csv: Path, results_dir: Path,
) -> int:
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

    ln = ["# Step 4 Full-Scale Docking Report (After Retry)", "",
          f"- **Total pairs**: {total} ({n_pos} positive, {n_neg} negative)"]
    if total:
        ln.append(f"- **Success**: {ok}/{total} ({ok / total * 100:.1f}%)")
    ln.append(f"- **Workers**: {workers}, **Retry wall time**: {wall_sec:.1f}s ({wall_sec / 3600:.1f}h)")
    if affs:
        ln.append(f"- **Affinity**: mean={statistics.mean(affs):.3f}, "
                   f"min={min(affs):.3f}, max={max(affs):.3f}")
    if times:
        ln.append(f"- **Per-pair time**: mean={statistics.mean(times):.1f}s, "
                   f"max={max(times):.1f}s")

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


# ── main ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--assets_dir", type=Path, required=True)
    p.add_argument("--pdb_dir", type=Path, required=True)
    p.add_argument("--output_dir", type=Path, required=True)
    p.add_argument("--results_dir", type=Path, required=True)
    p.add_argument("--enzymes_csv", type=Path, required=True)
    p.add_argument("--substrates_csv", type=Path, required=True)
    p.add_argument("--num_workers", type=int, default=12)
    p.add_argument("--fe_dist_min", type=float, default=1.5,
                   help="Relaxed FE distance minimum (default 1.5, was 2.0)")
    p.add_argument("--tier1_timeout", type=int, default=600,
                   help="Tier 1 timeout in seconds (default 600)")
    p.add_argument("--tier1_exhaustiveness", type=int, default=8)
    p.add_argument("--tier2_timeout", type=int, default=600,
                   help="Tier 2 timeout in seconds (default 600)")
    p.add_argument("--tier2_exhaustiveness", type=int, default=4)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S")

    for tag, p in [("assets_dir", args.assets_dir), ("pdb_dir", args.pdb_dir),
                    ("output_dir", args.output_dir)]:
        if not p.exists():
            LOG.error("Not found: %s (%s)", p, tag)
            return 2

    manifest = load_manifest(args.output_dir)
    progress = load_progress(args.output_dir)
    manifest_map = {t["dock_index"]: t for t in manifest}

    # Backup progress before retry
    backup = args.output_dir / "progress_before_retry.csv"
    if not backup.exists():
        shutil.copy2(str(args.output_dir / "progress.csv"), str(backup))
        LOG.info("Backed up progress.csv -> progress_before_retry.csv")

    # Load substrate heavy atom counts
    from rdkit import Chem
    heavy_map: dict[int, int] = {}
    for r in _csv_rows(args.substrates_csv):
        mol = Chem.MolFromSmiles(r.get("substrate_smiles", ""))
        if mol:
            heavy_map[_safe_int(r["substrate_index"])] = sum(
                1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)

    if not heavy_map:
        LOG.warning("No substrate heavy atom counts loaded — QC check disabled")

    # Classify failures
    fe_failures = []  # docked PDBQT exists, just needs re-postprocess
    timeout_failures = []  # needs full re-dock
    # Warn about never-attempted pairs (step4 interrupted mid-run)
    missing = [t["dock_index"] for t in manifest
               if t["dock_index"] not in progress]
    if missing:
        LOG.warning("%d pairs never attempted (not in progress.csv). "
                    "Run step4 with --resume first.", len(missing))

    for di, rec in progress.items():
        if rec["status"] not in ("failed", "exception"):
            continue
        err = rec.get("error", "")
        m = manifest_map.get(di)
        if m is None:
            continue
        if "FE-ligand distance" in err and "outside" in err:
            fe_failures.append(di)
        elif "timed out" in err:
            timeout_failures.append(di)
        else:
            timeout_failures.append(di)  # treat other failures as re-dock

    LOG.info("Failures: %d FE distance + %d timeout = %d total",
             len(fe_failures), len(timeout_failures),
             len(fe_failures) + len(timeout_failures))

    ok_before = sum(1 for v in progress.values() if v["status"] == "success")
    total_wall = 0.0

    # ── Tier 0: FE distance re-postprocess ───────────────────────────────
    if fe_failures:
        LOG.info("=== Tier 0: Re-postprocess %d FE distance failures (min=%.1fÅ) ===",
                 len(fe_failures), args.fe_dist_min)
        t0 = time.perf_counter()
        recovered = 0

        for di in fe_failures:
            m = manifest_map[di]
            res = repostprocess_fe_failure(
                dock_index=di,
                pdb_id=m["pdb_id"],
                substrate_index=m["substrate_index"],
                output_dir=str(args.output_dir),
                pdb_dir=str(args.pdb_dir),
                fe_dist_min=args.fe_dist_min,
                expected_heavy=heavy_map.get(m["substrate_index"]))

            # Preserve original affinity if we have it from the docking output
            old_aff = progress[di].get("affinity")
            if old_aff is not None and res["affinity"] is None:
                res["affinity"] = old_aff

            progress[di] = {k: res.get(k) for k in PROGRESS_FIELDS}
            if res["status"] == "success":
                recovered += 1

        elapsed = time.perf_counter() - t0
        total_wall += elapsed
        save_progress(args.output_dir, progress)
        LOG.info("Tier 0 done: %d/%d recovered in %.1fs", recovered, len(fe_failures), elapsed)

    # ── Tier 1: Re-dock timeouts with extended timeout ───────────────────
    if timeout_failures:
        LOG.info("=== Tier 1: Re-dock %d timeouts (timeout=%ds, exhaust=%d) ===",
                 len(timeout_failures), args.tier1_timeout, args.tier1_exhaustiveness)
        t0 = time.perf_counter()
        done_count = 0
        tier1_recovered = 0

        with ProcessPoolExecutor(max_workers=args.num_workers) as pool:
            futs = {}
            for di in timeout_failures:
                m = manifest_map[di]
                futs[pool.submit(
                    redock_pair,
                    di, m["pdb_id"], m["enzyme_index"], m["substrate_index"],
                    str(args.assets_dir), str(args.pdb_dir), str(args.output_dir),
                    heavy_map.get(m["substrate_index"]),
                    args.tier1_timeout, args.tier1_exhaustiveness,
                    args.fe_dist_min, "tier1",
                )] = di

            for fut in as_completed(futs):
                di = futs[fut]
                try:
                    r = fut.result()
                except Exception as exc:
                    r = dict(dock_index=di, status="failed", affinity=None,
                             pocket_atoms=0, fe_distance=None,
                             runtime_sec=0, error=f"Exception: {exc} [tier1]")

                progress[di] = {k: r.get(k) for k in PROGRESS_FIELDS}
                done_count += 1
                if r.get("status") == "success":
                    tier1_recovered += 1

                if done_count % 10 == 0 or done_count == len(timeout_failures):
                    elapsed = time.perf_counter() - t0
                    LOG.info("[Tier1 %d/%d] %.0fs, recovered=%d",
                             done_count, len(timeout_failures), elapsed, tier1_recovered)

                if done_count % 50 == 0:
                    save_progress(args.output_dir, progress)

        elapsed = time.perf_counter() - t0
        total_wall += elapsed
        save_progress(args.output_dir, progress)
        LOG.info("Tier 1 done: %d/%d recovered in %.1fs (%.1fh)",
                 tier1_recovered, len(timeout_failures), elapsed, elapsed / 3600)

        # ── Tier 2: Remaining timeouts with reduced exhaustiveness ───────
        still_failed = [di for di in timeout_failures
                        if progress[di]["status"] != "success"]

        if still_failed:
            LOG.info("=== Tier 2: Re-dock %d remaining (timeout=%ds, exhaust=%d) ===",
                     len(still_failed), args.tier2_timeout, args.tier2_exhaustiveness)
            t0 = time.perf_counter()
            done_count = 0
            tier2_recovered = 0

            with ProcessPoolExecutor(max_workers=args.num_workers) as pool:
                futs = {}
                for di in still_failed:
                    m = manifest_map[di]
                    futs[pool.submit(
                        redock_pair,
                        di, m["pdb_id"], m["enzyme_index"], m["substrate_index"],
                        str(args.assets_dir), str(args.pdb_dir), str(args.output_dir),
                        heavy_map.get(m["substrate_index"]),
                        args.tier2_timeout, args.tier2_exhaustiveness,
                        args.fe_dist_min, "tier2",
                    )] = di

                for fut in as_completed(futs):
                    di = futs[fut]
                    try:
                        r = fut.result()
                    except Exception as exc:
                        r = dict(dock_index=di, status="failed", affinity=None,
                                 pocket_atoms=0, fe_distance=None,
                                 runtime_sec=0, error=f"Exception: {exc} [tier2]")

                    progress[di] = {k: r.get(k) for k in PROGRESS_FIELDS}
                    done_count += 1
                    if r.get("status") == "success":
                        tier2_recovered += 1

                    if done_count % 10 == 0 or done_count == len(still_failed):
                        elapsed = time.perf_counter() - t0
                        LOG.info("[Tier2 %d/%d] %.0fs, recovered=%d",
                                 done_count, len(still_failed), elapsed, tier2_recovered)

                    if done_count % 50 == 0:
                        save_progress(args.output_dir, progress)

            elapsed = time.perf_counter() - t0
            total_wall += elapsed
            save_progress(args.output_dir, progress)
            LOG.info("Tier 2 done: %d/%d recovered in %.1fs (%.1fh)",
                     tier2_recovered, len(still_failed), elapsed, elapsed / 3600)

    # ── Regenerate outputs ───────────────────────────────────────────────
    ok_after = sum(1 for v in progress.values() if v["status"] == "success")
    n_ds = write_dataset(manifest, progress, args.enzymes_csv,
                         args.substrates_csv, args.results_dir)
    report_path = args.results_dir / "fullscale_report.md"
    write_report(report_path, manifest, progress, total_wall, args.num_workers)

    pct_before = ok_before / len(manifest) * 100
    pct_after = ok_after / len(manifest) * 100
    pos_ok = sum(1 for t in manifest if t["label"] == 1
                 and progress.get(t["dock_index"], {}).get("status") == "success")
    neg_ok = ok_after - pos_ok

    print(f"\n{'=' * 60}")
    print(f"Step 4b Retry Complete")
    print(f"  Before:  {ok_before}/{len(manifest)} ({pct_before:.1f}%)")
    print(f"  After:   {ok_after}/{len(manifest)} ({pct_after:.1f}%)")
    print(f"  Gained:  +{ok_after - ok_before} pairs")
    print(f"  Pos/Neg: {pos_ok} / {neg_ok} (ratio 1:{neg_ok / pos_ok:.2f})" if pos_ok else "")
    print(f"  Dataset: {n_ds} rows -> {args.results_dir / 'data.csv'}")
    print(f"  Report:  {report_path}")
    print(f"  Wall:    {total_wall:.0f}s ({total_wall / 3600:.1f}h)")
    print(f"{'=' * 60}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
