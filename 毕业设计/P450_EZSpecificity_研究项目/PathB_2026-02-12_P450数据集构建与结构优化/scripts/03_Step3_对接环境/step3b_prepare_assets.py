"""Step3B: Batch-prepare receptor PDBQT and ligand PDBQT/SDF for docking."""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.grid_locator import locate_grid_box
from lib.ligand_prep import prepare_ligand_from_smiles
from lib.receptor_prep import prepare_receptor_pdbqt

LOG = logging.getLogger(__name__)


def _norm_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_")


def load_unique_pdb_ids(enzymes_csv: Path) -> dict[str, list[int]]:
    """Return {PDB_ID: [enzyme_indices]} with deduplication by PDB_ID."""
    pdb_map: dict[str, list[int]] = {}
    with enzymes_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            pdb_id = nr.get("pdb_id", "").upper()
            enz_idx = int(float(nr["enzyme_index"]))
            if pdb_id:
                pdb_map.setdefault(pdb_id, []).append(enz_idx)
    return pdb_map


def load_substrates(substrates_csv: Path) -> list[tuple[int, str]]:
    pairs: list[tuple[int, str]] = []
    with substrates_csv.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            idx = int(float(nr["substrate_index"]))
            smiles = nr.get("substrate_smiles", "")
            if smiles:
                pairs.append((idx, smiles))
    return sorted(pairs)


def index_pdb_files(pdb_dir: Path) -> dict[str, Path]:
    return {p.stem.upper(): p for p in pdb_dir.rglob("*.pdb")}


def write_summary(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--enzymes_csv", type=Path, required=True, help="B6 Enzymes.csv")
    p.add_argument("--substrates_csv", type=Path, required=True, help="B6 Substrates.csv")
    p.add_argument("--pdb_dir", type=Path, required=True, help="PathA PDB directory")
    p.add_argument("--output_dir", type=Path, required=True, help="Output base (data/03_Step3_对接预实验/)")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    for label, path in [("enzymes_csv", args.enzymes_csv),
                         ("substrates_csv", args.substrates_csv),
                         ("pdb_dir", args.pdb_dir)]:
        if not path.exists():
            LOG.error("%s not found: %s", label, path)
            return 2

    rec_pdbqt_dir = args.output_dir / "receptors" / "pdbqt"
    rec_grid_dir = args.output_dir / "receptors" / "grid_boxes"
    lig_pdbqt_dir = args.output_dir / "ligands" / "pdbqt"
    lig_sdf_dir = args.output_dir / "ligands" / "sdf_3d"
    for d in [rec_pdbqt_dir, rec_grid_dir, lig_pdbqt_dir, lig_sdf_dir]:
        d.mkdir(parents=True, exist_ok=True)

    pdb_index = index_pdb_files(args.pdb_dir)
    pdb_map = load_unique_pdb_ids(args.enzymes_csv)
    substrates = load_substrates(args.substrates_csv)

    LOG.info("Unique PDB_IDs: %d, Substrates: %d, PDB files indexed: %d",
             len(pdb_map), len(substrates), len(pdb_index))

    # === Prepare receptors (deduplicated by PDB_ID) ===
    rec_rows: list[dict[str, Any]] = []
    rec_ok = 0
    for i, (pdb_id, enz_indices) in enumerate(sorted(pdb_map.items()), 1):
        row: dict[str, Any] = {"PDB_ID": pdb_id, "Enzyme_Indices": str(enz_indices),
                                "grid_ok": False, "receptor_ok": False, "error": ""}

        if i % 20 == 0:
            LOG.info("Receptor %d/%d: %s", i, len(pdb_map), pdb_id)

        pdb_path = pdb_index.get(pdb_id)
        if pdb_path is None:
            row["error"] = "PDB file not found"
            rec_rows.append(row)
            continue

        grid = locate_grid_box(pdb_path, output_dir=rec_grid_dir, pdb_id=pdb_id)
        if grid is None:
            row["error"] = "HEM FE not found"
            rec_rows.append(row)
            continue
        row["grid_ok"] = True

        out_pdbqt = rec_pdbqt_dir / f"{pdb_id}.pdbqt"
        prep = prepare_receptor_pdbqt(pdb_path, out_pdbqt)
        row["receptor_ok"] = prep.success
        row["atom_count"] = prep.atom_count
        row["error"] = prep.error
        rec_rows.append(row)
        if prep.success:
            rec_ok += 1

    # === Prepare ligands ===
    lig_rows: list[dict[str, Any]] = []
    lig_ok = 0
    for i, (sub_idx, smiles) in enumerate(substrates, 1):
        if i % 20 == 0:
            LOG.info("Ligand %d/%d: substrate_%d", i, len(substrates), sub_idx)

        prep = prepare_ligand_from_smiles(
            smiles=smiles,
            ligand_id=f"substrate_{sub_idx}",
            pdbqt_dir=lig_pdbqt_dir,
            sdf_dir=lig_sdf_dir,
        )
        lig_rows.append({
            "Substrate_Index": sub_idx, "success": prep.success,
            "atom_count": prep.atom_count, "error": prep.error,
        })
        if prep.success:
            lig_ok += 1

    # === Summaries ===
    write_summary(args.output_dir / "receptor_prep_summary.csv", rec_rows,
                  ["PDB_ID", "Enzyme_Indices", "grid_ok", "receptor_ok", "atom_count", "error"])
    write_summary(args.output_dir / "ligand_prep_summary.csv", lig_rows,
                  ["Substrate_Index", "success", "atom_count", "error"])

    rec_fail = len(pdb_map) - rec_ok
    lig_fail = len(substrates) - lig_ok
    print(f"\n=== Step3B Summary ===")
    print(f"Receptors: {rec_ok}/{len(pdb_map)} OK ({rec_fail} failed)")
    print(f"Ligands:   {lig_ok}/{len(substrates)} OK ({lig_fail} failed)")
    print(f"Output:    {args.output_dir}")

    return 0 if (rec_fail + lig_fail) == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
