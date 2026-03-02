"""Locate the docking grid box from HEM FE atom coordinates in a PDB file."""
from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from Bio.PDB import PDBParser

logger = logging.getLogger(__name__)

HEME_RESNAMES = frozenset({"HEM", "HEC", "HEA", "HEO", "DHE", "HEB"})


@dataclass(frozen=True)
class GridBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float = 22.5
    size_y: float = 22.5
    size_z: float = 22.5


def find_heme_fe_center(pdb_path: Path) -> Optional[tuple[float, float, float]]:
    """Parse PDB and return (x, y, z) of the first HEM FE atom found."""
    try:
        structure = PDBParser(QUIET=True).get_structure("receptor", str(pdb_path))
    except Exception as exc:
        logger.warning("Failed to parse PDB %s: %s", pdb_path, exc)
        return None

    for model in structure:
        for residue in model.get_residues():
            if residue.get_resname().strip().upper() not in HEME_RESNAMES:
                continue
            for atom in residue.get_atoms():
                if atom.get_name().strip().upper() == "FE":
                    x, y, z = atom.get_coord()
                    return float(x), float(y), float(z)
    return None


def locate_grid_box(
    pdb_path: Path,
    output_dir: Optional[Path] = None,
    pdb_id: Optional[str] = None,
) -> Optional[GridBox]:
    """Locate docking grid box centered on HEM FE atom.

    Returns GridBox if FE found, None otherwise.
    Optionally saves grid parameters as JSON to output_dir.
    """
    fe_center = find_heme_fe_center(pdb_path)
    if fe_center is None:
        logger.warning("HEM FE atom not found in %s", pdb_path)
        return None

    grid = GridBox(*fe_center)
    resolved_id = (pdb_id or pdb_path.stem).upper()
    logger.info("Grid box for %s: center=(%.2f, %.2f, %.2f)", resolved_id, *fe_center)

    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        json_path = output_dir / f"{resolved_id}_grid.json"
        json_path.write_text(
            json.dumps({"pdb_id": resolved_id, **asdict(grid)}, indent=2),
            encoding="utf-8",
        )

    return grid
