"""Step3D: Run 50-pair pilot docking (5 positives + 45 negatives, parallel)."""
from __future__ import annotations

import argparse
import csv
import json
import logging
import statistics
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

LOG = logging.getLogger(__name__)


def _norm_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_")


def load_enzyme_map(enzymes_csv: Path) -> dict[int, str]:
    mapping: dict[int, str] = {}
    with enzymes_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            mapping[int(float(nr["enzyme_index"]))] = nr.get("pdb_id", "").upper()
    return mapping


def load_first_n_positives(data_csv: Path, n: int) -> list[tuple[int, int, int]]:
    out: list[tuple[int, int, int]] = []
    with data_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            if int(float(nr["label"])) != 1:
                continue
            out.append((int(nr["dock_index"]), int(nr["enzyme_index"]), int(nr["substrate_index"])))
            if len(out) >= n:
                break
    return out


def load_substrate_heavy_atoms(substrates_csv: Path) -> dict[int, int]:
    from rdkit import Chem
    mapping: dict[int, int] = {}
    with substrates_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            idx = int(float(nr["substrate_index"]))
            mol = Chem.MolFromSmiles(nr.get("substrate_smiles", ""))
            if mol:
                mapping[idx] = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    return mapping


def find_pdb(pdb_dir: Path, pdb_id: str) -> Optional[Path]:
    for name in [f"{pdb_id}.pdb", f"{pdb_id.upper()}.pdb"]:
        p = pdb_dir / name
        if p.exists():
            return p
    return None


def load_grid(grid_json: Path):
    from lib.grid_locator import GridBox
    d = json.loads(grid_json.read_text(encoding="utf-8"))
    return GridBox(float(d["center_x"]), float(d["center_y"]), float(d["center_z"]),
                   float(d.get("size_x", 22.5)), float(d.get("size_y", 22.5)), float(d.get("size_z", 22.5)))


def run_single_pair(
    dock_index: int, pdb_id: str, enzyme_index: int, substrate_index: int,
    label: int, source: str,
    assets_dir: str, pdb_dir: str, output_dir: str,
    expected_heavy: Optional[int],
) -> dict[str, Any]:
    """Execute one docking + postprocess. Designed for ProcessPoolExecutor."""
    # Lazy imports for Windows spawn compatibility
    from lib.postprocess import postprocess_docked_pose
    from lib.vina_runner import run_vina_docking

    assets = Path(assets_dir)
    out = Path(output_dir)
    pdb_root = Path(pdb_dir)

    result: dict[str, Any] = {
        "dock_index": dock_index, "pdb_id": pdb_id,
        "enzyme_index": enzyme_index, "substrate_index": substrate_index,
        "label": label, "source": source,
        "status": "failed", "affinity": None, "pocket_atoms": 0,
        "fe_distance": None, "runtime_sec": 0.0, "error": "",
    }

    rec_pdbqt = assets / "receptors" / "pdbqt" / f"{pdb_id}.pdbqt"
    grid_json = assets / "receptors" / "grid_boxes" / f"{pdb_id}_grid.json"
    lig_pdbqt = assets / "ligands" / "pdbqt" / f"substrate_{substrate_index}.pdbqt"
    orig_pdb = find_pdb(pdb_root, pdb_id)
    docked_out = out / "docking_output" / f"{dock_index}.pdbqt"

    for p, name in [(rec_pdbqt, "receptor"), (grid_json, "grid"), (lig_pdbqt, "ligand")]:
        if not p.exists():
            result["error"] = f"Missing {name}: {p.name}"
            return result
    if orig_pdb is None:
        result["error"] = f"Original PDB not found for {pdb_id}"
        return result

    try:
        grid = load_grid(grid_json)
    except Exception as exc:
        result["error"] = f"Grid parse: {exc}"
        return result

    dock = run_vina_docking(rec_pdbqt, lig_pdbqt, docked_out, grid)
    result["runtime_sec"] = dock.runtime_sec
    result["affinity"] = dock.best_affinity
    if not dock.success:
        result["error"] = f"Docking: {dock.error}"
        return result

    post = postprocess_docked_pose(
        docked_out, orig_pdb, dock_index,
        out / "pocket", out / "raw_ligand",
        expected_heavy_atoms=expected_heavy,
    )
    if not post.success:
        result["error"] = f"Postprocess: {post.error}"
        return result

    result["status"] = "success"
    result["pocket_atoms"] = post.pocket_atoms
    result["fe_distance"] = post.fe_ligand_distance
    return result


def write_report(path: Path, results: list[dict], wall_sec: float, num_workers: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = Counter(r["status"] for r in results)
    affs = [r["affinity"] for r in results if r["status"] == "success" and r["affinity"] is not None]
    times = [r["runtime_sec"] for r in results if r["runtime_sec"]]

    lines = [
        "# Step3 Pilot 50 Report", "",
        f"- Total: {len(results)}, Success: {counts.get('success', 0)}",
        f"- Success rate: {counts.get('success', 0) / len(results) * 100:.1f}%",
        f"- Workers: {num_workers}, Wall time: {wall_sec:.1f}s",
    ]
    if affs:
        lines.append(f"- Affinity: mean={statistics.mean(affs):.3f}, "
                      f"min={min(affs):.3f}, max={max(affs):.3f}")
    if times:
        lines.append(f"- Docking time: mean={statistics.mean(times):.1f}s, max={max(times):.1f}s")

    lines.extend(["", "## Status Breakdown", ""])
    for status, n in sorted(counts.items()):
        lines.append(f"- {status}: {n}")

    lines.extend(["", "## Results", "",
        "| Dock_Index | PDB_ID | Label | Source | Affinity | Pocket | FE_Dist | Status | Error |",
        "|---:|---|---:|---|---:|---:|---:|---|---|"])
    for r in sorted(results, key=lambda x: x["dock_index"]):
        aff = "" if r["affinity"] is None else f"{r['affinity']:.3f}"
        fed = "" if r["fe_distance"] is None else f"{r['fe_distance']:.2f}"
        lines.append(
            f"| {r['dock_index']} | {r['pdb_id']} | {r['label']} | {r['source']} | "
            f"{aff} | {r['pocket_atoms']} | {fed} | {r['status']} | {r['error']} |"
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data_csv", type=Path, required=True, help="B6 data.csv")
    p.add_argument("--enzymes_csv", type=Path, required=True, help="B6 Enzymes.csv")
    p.add_argument("--substrates_csv", type=Path, required=True, help="B6 Substrates.csv")
    p.add_argument("--assets_dir", type=Path, required=True, help="Prepared assets from Step3B")
    p.add_argument("--pdb_dir", type=Path, required=True, help="PathA original PDB directory")
    p.add_argument("--output_dir", type=Path, required=True, help="Pilot output (data/.../pilot_50/)")
    p.add_argument("--results_dir", type=Path, required=True, help="Report output (results/.../)")
    p.add_argument("--num_workers", type=int, default=12, help="Parallel workers")
    p.add_argument("--random_seed", type=int, default=2026)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    for path in [args.data_csv, args.enzymes_csv, args.substrates_csv, args.assets_dir, args.pdb_dir]:
        if not path.exists():
            LOG.error("Not found: %s", path)
            return 2

    for d in [args.output_dir / "docking_output", args.output_dir / "pocket",
              args.output_dir / "raw_ligand", args.results_dir]:
        d.mkdir(parents=True, exist_ok=True)

    enzyme_map = load_enzyme_map(args.enzymes_csv)
    positives_raw = load_first_n_positives(args.data_csv, n=5)
    if len(positives_raw) < 5:
        LOG.error("Need >= 5 positives, found %d", len(positives_raw))
        return 1

    # Build positive tasks
    tasks: list[dict[str, Any]] = []
    for di, ei, si in positives_raw:
        pdb_id = enzyme_map.get(ei, "").upper()
        tasks.append({"dock_index": di, "pdb_id": pdb_id, "enzyme_index": ei,
                       "substrate_index": si, "label": 1, "source": "positive"})

    # Generate negatives via negative_sampler
    from lib.negative_sampler import generate_negative_pairs

    mini_csv = args.output_dir / "_pilot5_positives.csv"
    neg_csv = args.output_dir / "_pilot45_negatives.csv"
    with mini_csv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Dock Index", "Enzyme Index", "Substrate Index", "Label"])
        for t in tasks:
            w.writerow([t["dock_index"], t["enzyme_index"], t["substrate_index"], 1])

    neg_res = generate_negative_pairs(mini_csv, args.enzymes_csv, neg_csv,
                                       ratio=9.0, random_seed=args.random_seed, dock_index_start=1000)
    if not neg_res.success:
        LOG.error("Negative sampling failed: %s", neg_res.error)
        return 1

    with neg_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            tasks.append({
                "dock_index": int(nr["dock_index"]), "pdb_id": nr.get("pdb_id", "").upper(),
                "enzyme_index": int(nr["enzyme_index"]), "substrate_index": int(nr["substrate_index"]),
                "label": 0, "source": "negative",
            })

    if len(tasks) < 50:
        LOG.error("Need 50 tasks (5 pos + 45 neg), only got %d. "
                  "Check negative sampling and enzyme pool size.", len(tasks))
        return 1
    tasks = tasks[:50]
    LOG.info("Tasks: %d (%d positive, %d negative)", len(tasks),
             sum(1 for t in tasks if t["label"] == 1),
             sum(1 for t in tasks if t["label"] == 0))

    # Save pair manifest
    manifest = args.output_dir / "pilot_50_pairs.csv"
    with manifest.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Dock_Index", "PDB_ID", "Enzyme_Index",
                                           "Substrate_Index", "Label", "Source"])
        w.writeheader()
        for t in tasks:
            w.writerow({k.replace("dock_index", "Dock_Index").replace("pdb_id", "PDB_ID")
                        .replace("enzyme_index", "Enzyme_Index").replace("substrate_index", "Substrate_Index")
                        .replace("label", "Label").replace("source", "Source"): v
                        for k, v in t.items()})

    expected_heavy = load_substrate_heavy_atoms(args.substrates_csv)

    # Parallel docking
    t0 = time.perf_counter()
    results: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {
            executor.submit(
                run_single_pair,
                t["dock_index"], t["pdb_id"], t["enzyme_index"], t["substrate_index"],
                t["label"], t["source"],
                str(args.assets_dir), str(args.pdb_dir), str(args.output_dir),
                expected_heavy.get(t["substrate_index"]),
            ): t["dock_index"]
            for t in tasks
        }
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                res = fut.result()
            except Exception as exc:
                di = futures[fut]
                res = {"dock_index": di, "pdb_id": "", "enzyme_index": -1,
                       "substrate_index": -1, "label": -1, "source": "error",
                       "status": "exception", "affinity": None, "pocket_atoms": 0,
                       "fe_distance": None, "runtime_sec": 0, "error": str(exc)}
            results.append(res)
            if i % 10 == 0 or i == len(tasks):
                LOG.info("Completed %d/%d", i, len(tasks))

    wall_sec = time.perf_counter() - t0
    results.sort(key=lambda r: r["dock_index"])

    # Save results CSV
    results_csv = args.output_dir / "pilot_50_results.csv"
    with results_csv.open("w", encoding="utf-8", newline="") as f:
        fields = ["dock_index", "pdb_id", "enzyme_index", "substrate_index",
                  "label", "source", "status", "affinity", "pocket_atoms",
                  "fe_distance", "runtime_sec", "error"]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    report = args.results_dir / "pilot_50_report.md"
    write_report(report, results, wall_sec, args.num_workers)

    ok = sum(1 for r in results if r["status"] == "success")
    print(f"\n=== Pilot 50 Summary ===")
    print(f"Tasks: {len(results)}, Success: {ok} ({ok / len(results) * 100:.1f}%)")
    print(f"Wall time: {wall_sec:.1f}s")
    print(f"Report: {report}")

    return 0 if ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
