"""Step3C: Run 5-pair positive docking pilot (serial, for debugging)."""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from pathlib import Path
from typing import Any, Optional

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.grid_locator import GridBox
from lib.postprocess import postprocess_docked_pose
from lib.vina_runner import run_vina_docking

LOG = logging.getLogger(__name__)


def _norm_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_")


def load_positive_pairs(data_csv: Path, limit: int = 5) -> list[dict[str, int]]:
    pairs: list[dict[str, int]] = []
    with data_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            if int(float(nr["label"])) != 1:
                continue
            pairs.append({
                "dock_index": int(nr["dock_index"]),
                "enzyme_index": int(nr["enzyme_index"]),
                "substrate_index": int(nr["substrate_index"]),
            })
            if len(pairs) >= limit:
                break
    return pairs


def load_enzyme_map(enzymes_csv: Path) -> dict[int, str]:
    mapping: dict[int, str] = {}
    with enzymes_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            mapping[int(float(nr["enzyme_index"]))] = nr.get("pdb_id", "").upper()
    return mapping


def load_grid(grid_json: Path) -> GridBox:
    d = json.loads(grid_json.read_text(encoding="utf-8"))
    return GridBox(float(d["center_x"]), float(d["center_y"]), float(d["center_z"]),
                   float(d.get("size_x", 22.5)), float(d.get("size_y", 22.5)), float(d.get("size_z", 22.5)))


def find_original_pdb(pdb_dir: Path, pdb_id: str) -> Optional[Path]:
    for name in [f"{pdb_id}.pdb", f"{pdb_id.upper()}.pdb", f"{pdb_id.lower()}.pdb"]:
        p = pdb_dir / name
        if p.exists():
            return p
    return None


def save_report(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    ok = sum(1 for r in rows if r["status"] == "success")
    lines = [
        "# Step3 Pilot 5 Report", "",
        f"- Total: {len(rows)}, Success: {ok}, Rate: {ok / len(rows) * 100:.1f}%", "",
        "| Dock_Index | PDB_ID | Affinity | Pocket_Atoms | FE_Dist | Status | Error |",
        "|---:|---|---:|---:|---:|---|---|",
    ]
    for r in rows:
        aff = "" if r["affinity"] is None else f"{r['affinity']:.3f}"
        fed = "" if r["fe_distance"] is None else f"{r['fe_distance']:.2f}"
        lines.append(
            f"| {r['dock_index']} | {r['pdb_id']} | {aff} | "
            f"{r['pocket_atoms']} | {fed} | {r['status']} | {r['error']} |"
        )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data_csv", type=Path, required=True, help="B6 data.csv")
    p.add_argument("--enzymes_csv", type=Path, required=True, help="B6 Enzymes.csv")
    p.add_argument("--assets_dir", type=Path, required=True, help="Prepared assets from Step3B")
    p.add_argument("--pdb_dir", type=Path, required=True, help="PathA original PDB directory")
    p.add_argument("--output_dir", type=Path, required=True, help="Pilot output (data/.../pilot_5/)")
    p.add_argument("--results_dir", type=Path, required=True, help="Report output (results/.../)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    for path in [args.data_csv, args.enzymes_csv, args.assets_dir, args.pdb_dir]:
        if not path.exists():
            LOG.error("Not found: %s", path)
            return 2

    docking_dir = args.output_dir / "docking_output"
    pocket_dir = args.output_dir / "pocket"
    ligand_dir = args.output_dir / "raw_ligand"
    for d in [docking_dir, pocket_dir, ligand_dir]:
        d.mkdir(parents=True, exist_ok=True)

    pairs = load_positive_pairs(args.data_csv, limit=5)
    enzyme_map = load_enzyme_map(args.enzymes_csv)

    if len(pairs) < 5:
        LOG.error("Need >= 5 positives, found %d", len(pairs))
        return 1

    rows: list[dict[str, Any]] = []

    for i, pair in enumerate(pairs, 1):
        pdb_id = enzyme_map.get(pair["enzyme_index"], "").upper()
        row: dict[str, Any] = {
            "dock_index": pair["dock_index"], "pdb_id": pdb_id,
            "affinity": None, "pocket_atoms": 0, "fe_distance": None,
            "status": "failed", "error": "",
        }
        LOG.info("[%d/5] Dock_Index=%d PDB=%s Substrate=%d",
                 i, pair["dock_index"], pdb_id, pair["substrate_index"])

        rec_pdbqt = args.assets_dir / "receptors" / "pdbqt" / f"{pdb_id}.pdbqt"
        grid_json = args.assets_dir / "receptors" / "grid_boxes" / f"{pdb_id}_grid.json"
        lig_pdbqt = args.assets_dir / "ligands" / "pdbqt" / f"substrate_{pair['substrate_index']}.pdbqt"
        orig_pdb = find_original_pdb(args.pdb_dir, pdb_id)
        docked_out = docking_dir / f"{pair['dock_index']}.pdbqt"

        # Pre-flight checks
        missing = []
        if not pdb_id:
            missing.append("PDB_ID mapping")
        if not rec_pdbqt.exists():
            missing.append(f"receptor {rec_pdbqt.name}")
        if not grid_json.exists():
            missing.append(f"grid {grid_json.name}")
        if not lig_pdbqt.exists():
            missing.append(f"ligand {lig_pdbqt.name}")
        if orig_pdb is None:
            missing.append(f"original PDB for {pdb_id}")
        if missing:
            row["error"] = f"Missing: {', '.join(missing)}"
            rows.append(row)
            continue

        # Dock
        try:
            grid = load_grid(grid_json)
        except Exception as exc:
            row["error"] = f"Grid parse failed: {exc}"
            rows.append(row)
            continue
        dock_res = run_vina_docking(rec_pdbqt, lig_pdbqt, docked_out, grid)
        if not dock_res.success:
            row["error"] = f"Docking: {dock_res.error}"
            rows.append(row)
            continue

        row["affinity"] = dock_res.best_affinity
        LOG.info("  Affinity=%.3f kcal/mol (%.1fs)", dock_res.best_affinity, dock_res.runtime_sec)

        # Postprocess
        post = postprocess_docked_pose(docked_out, orig_pdb, pair["dock_index"],
                                       pocket_dir, ligand_dir)
        if not post.success:
            row["error"] = f"Postprocess: {post.error}"
            rows.append(row)
            continue

        row["pocket_atoms"] = post.pocket_atoms
        row["fe_distance"] = post.fe_ligand_distance
        row["status"] = "success"
        LOG.info("  Pocket=%d atoms, FE_dist=%.2fÅ", post.pocket_atoms, post.fe_ligand_distance)
        rows.append(row)

    # Console summary
    print(f"\n=== Pilot 5 Results ===")
    print(f"{'DockIdx':>8} {'PDB':>6} {'Affinity':>10} {'Pocket':>8} {'FE_Dist':>8} {'Status':>8}")
    for r in rows:
        aff = "" if r["affinity"] is None else f"{r['affinity']:.3f}"
        fed = "" if r["fe_distance"] is None else f"{r['fe_distance']:.2f}"
        print(f"{r['dock_index']:>8} {r['pdb_id']:>6} {aff:>10} "
              f"{r['pocket_atoms']:>8} {fed:>8} {r['status']:>8}")

    report = args.results_dir / "pilot_5_report.md"
    save_report(report, rows)
    print(f"\nReport: {report}")

    return 0 if all(r["status"] == "success" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
