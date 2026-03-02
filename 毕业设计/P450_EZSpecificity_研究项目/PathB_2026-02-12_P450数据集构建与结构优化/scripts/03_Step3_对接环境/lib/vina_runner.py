"""Execute a single AutoDock Vina docking run via subprocess."""
from __future__ import annotations

import logging
import re
import shutil
import subprocess
import time
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .grid_locator import GridBox

logger = logging.getLogger(__name__)

DEFAULT_VINA_EXE = Path(r"D:\autodock\vina.exe")
DEFAULT_TMP_DIR = Path(r"D:\autodock\tmp")

# Matches first docking mode line: "   1    -7.3    0.000    0.000"
AFFINITY_RE = re.compile(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s", re.MULTILINE)


@dataclass
class DockingResult:
    output_path: Optional[Path]
    best_affinity: Optional[float]
    success: bool
    error: str = ""
    runtime_sec: float = 0.0


def _needs_ascii_copy(path: Path) -> bool:
    return not str(path).isascii()


def _copy_to_tmp(path: Path, tmp_dir: Path, tag: str) -> Path:
    dest = tmp_dir / f"{tag}_{uuid.uuid4().hex[:8]}{path.suffix}"
    shutil.copy2(path, dest)
    return dest


def run_vina_docking(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    output_pdbqt: Path,
    grid: GridBox,
    vina_exe: Path = DEFAULT_VINA_EXE,
    tmp_dir: Path = DEFAULT_TMP_DIR,
    exhaustiveness: int = 8,
    num_modes: int = 9,
    energy_range: int = 3,
    cpu: int = 1,
    timeout_sec: int = 300,
) -> DockingResult:
    """Run Vina docking and return the best binding affinity."""
    result = DockingResult(output_path=output_pdbqt, best_affinity=None, success=False)
    start = time.perf_counter()

    for label, path in [("Receptor", receptor_pdbqt), ("Ligand", ligand_pdbqt),
                         ("Vina", vina_exe)]:
        if not path.exists():
            result.error = f"{label} not found: {path}"
            return result

    tmp_dir.mkdir(parents=True, exist_ok=True)
    temp_files: list[Path] = []
    receptor_arg, ligand_arg, output_arg = receptor_pdbqt, ligand_pdbqt, output_pdbqt
    copy_back = False

    try:
        # Ensure output directory exists before Vina runs (ASCII path case)
        output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

        if _needs_ascii_copy(receptor_pdbqt):
            receptor_arg = _copy_to_tmp(receptor_pdbqt, tmp_dir, "rec")
            temp_files.append(receptor_arg)
        if _needs_ascii_copy(ligand_pdbqt):
            ligand_arg = _copy_to_tmp(ligand_pdbqt, tmp_dir, "lig")
            temp_files.append(ligand_arg)
        if _needs_ascii_copy(output_pdbqt):
            output_arg = tmp_dir / f"docked_{uuid.uuid4().hex[:8]}.pdbqt"
            temp_files.append(output_arg)
            copy_back = True

        cmd = [
            str(vina_exe),
            "--receptor", str(receptor_arg),
            "--ligand", str(ligand_arg),
            "--center_x", f"{grid.center_x:.3f}",
            "--center_y", f"{grid.center_y:.3f}",
            "--center_z", f"{grid.center_z:.3f}",
            "--size_x", f"{grid.size_x:.3f}",
            "--size_y", f"{grid.size_y:.3f}",
            "--size_z", f"{grid.size_z:.3f}",
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(num_modes),
            "--energy_range", str(energy_range),
            "--cpu", str(cpu),
            "--out", str(output_arg),
        ]

        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
        result.runtime_sec = time.perf_counter() - start

        if proc.returncode != 0:
            result.error = (proc.stderr or proc.stdout or "Unknown Vina error").strip()
            return result

        if not output_arg.exists():
            result.error = f"Vina finished but output missing: {output_arg}"
            return result

        affinity = AFFINITY_RE.search(proc.stdout)
        if affinity is None:
            result.error = "Failed to parse binding affinity from Vina output"
            return result

        if copy_back:
            shutil.copy2(output_arg, output_pdbqt)

        result.best_affinity = float(affinity.group(1))
        result.success = True
        return result

    except subprocess.TimeoutExpired:
        result.runtime_sec = time.perf_counter() - start
        result.error = f"Vina timed out after {timeout_sec}s"
        return result
    except Exception as exc:
        result.runtime_sec = time.perf_counter() - start
        result.error = f"Unexpected docking error: {exc}"
        return result
    finally:
        for path in temp_files:
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass
