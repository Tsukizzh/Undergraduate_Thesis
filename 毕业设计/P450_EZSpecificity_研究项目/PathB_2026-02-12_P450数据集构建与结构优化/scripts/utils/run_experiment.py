"""
PathB experiment runner: orchestrates the structure feature → inference → evaluation
pipeline for a single factorial experiment.

Usage:
  # Single experiment
  python run_experiment.py --name EXP01_B6_10A_noHeme --pocket_radius 10.0

  # With heme
  python run_experiment.py --name EXP02_B6_10A_Heme --pocket_radius 10.0 --include_heme

  # Run all 4 from config
  python run_experiment.py --run_all --config ../configs/base.yaml

Pipeline per experiment:
  Step 8.1: extract_pocket_ligand.py  (pocket + ligand extraction)
  Step 8.2: align_ligand.py           (ligand atom alignment)
  Step 8.3: generate_structure_lmdb.py (structure feature LMDB)
  Step 9:   inference                  (model prediction)
  Step 10:  analysis                   (evaluation metrics)
"""

import argparse
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import yaml

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def resolve_paths() -> Dict[str, Path]:
    """Auto-resolve project paths from script location.

    Layout: scripts/utils/run_experiment.py
            ^^^^^^ ^^^^^ <- two .parent calls to reach PathB root
    """
    script_dir = Path(__file__).resolve().parent      # scripts/utils/
    pathb_root = script_dir.parent.parent             # PathB_2026-02-12_.../
    research_root = pathb_root.parent                 # P450_EZSpecificity_研究项目/
    project_root = research_root.parent.parent        # EZSpecificity_Project/

    # Find PathA directory (deterministic: sorted, pick first)
    patha_candidates = sorted(research_root.glob("PathA_*"))
    if not patha_candidates:
        raise FileNotFoundError("Cannot find PathA directory under " + str(research_root))
    patha_root = patha_candidates[0]

    return {
        "project_root": project_root,
        "research_root": research_root,
        "pathb_root": pathb_root,
        "patha_root": patha_root,
        "patha_data": patha_root / "data",
        "patha_scripts": patha_root / "scripts",
        "pathb_data": pathb_root / "data",
        "pathb_scripts": pathb_root / "scripts",
        "shared_features": pathb_root / "data" / "00_shared" / "features",
        "shared_datasets": pathb_root / "data" / "00_shared" / "datasets",
        "mapping_csv": patha_root / "data" / "04_Step4_格式修正后数据" / "dock_index_mapping.csv",
        "pdb_dir": patha_root / "data" / "01_Step1_PDB文件",
        "checkpoint_dir": project_root / "saved_model" / "model" / "run_0",
    }


def create_experiment_dir(pathb_root: Path, exp_name: str) -> Path:
    """Create isolated experiment directory under data/02_Step2_因子实验/."""
    exp_dir = pathb_root / "data" / "02_Step2_因子实验" / exp_name
    exp_dir.mkdir(parents=True, exist_ok=True)
    (exp_dir / "structure_features").mkdir(exist_ok=True)
    return exp_dir


def save_config_snapshot(exp_dir: Path, config: dict):
    """Save config snapshot for reproducibility."""
    with open(exp_dir / "config.yaml", "w", encoding="utf-8") as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)


def run_step(name: str, cmd: List[str], cwd: Optional[Path] = None) -> bool:
    """Run a pipeline step as subprocess."""
    logger.info("--- [%s] ---", name)
    logger.info("CMD: %s", " ".join(str(c) for c in cmd))
    t0 = time.time()
    try:
        result = subprocess.run(
            cmd, cwd=str(cwd) if cwd else None,
            capture_output=True, text=True, timeout=7200  # 2h timeout
        )
        elapsed = time.time() - t0
        if result.returncode != 0:
            logger.error("[%s] FAILED (%.1fs)\nSTDERR:\n%s", name, elapsed, result.stderr[-2000:])
            return False
        logger.info("[%s] OK (%.1fs)", name, elapsed)
        if result.stdout.strip():
            # Log last few lines
            lines = result.stdout.strip().split("\n")
            for line in lines[-5:]:
                logger.info("  %s", line)
        return True
    except subprocess.TimeoutExpired:
        logger.error("[%s] TIMEOUT (>2h)", name)
        return False
    except Exception as e:
        logger.error("[%s] ERROR: %s", name, e)
        return False


def run_single_experiment(
    exp_name: str,
    pocket_radius: float,
    include_heme: bool,
    paths: Dict[str, Path],
    dry_run: bool = False,
) -> bool:
    """Run complete pipeline for one experiment."""
    logger.info("=" * 60)
    logger.info("EXPERIMENT: %s", exp_name)
    logger.info("  pocket_radius=%.1f  include_heme=%s", pocket_radius, include_heme)
    logger.info("=" * 60)

    pathb_root = paths["pathb_root"]
    exp_dir = create_experiment_dir(pathb_root, exp_name)
    str_tmp = exp_dir / "structure_features"

    # Save config snapshot
    config = {
        "experiment_name": exp_name,
        "pocket_radius": pocket_radius,
        "include_heme": include_heme,
        "timestamp": datetime.now().isoformat(),
        "paths": {k: str(v) for k, v in paths.items()},
    }
    save_config_snapshot(exp_dir, config)

    if dry_run:
        logger.info("[DRY RUN] Would execute pipeline for %s", exp_name)
        logger.info("  Experiment dir: %s", exp_dir)
        return True

    python = sys.executable

    # --- Preflight checks ---
    required_files = [
        paths["mapping_csv"],
        paths["pdb_dir"],
        paths["shared_features"] / "enzyme_features.lmdb",
        paths["shared_features"] / "reaction_features.lmdb",
        paths["shared_features"] / "grover_fingerprint.lmdb",
        paths["shared_features"] / "morgan_fingerprint.npy",
        paths["checkpoint_dir"] / "models" / "best-checkpoint.ckpt",
    ]
    for f in required_files:
        if not f.exists():
            logger.error("PREFLIGHT FAIL: missing %s", f)
            return False
    logger.info("Preflight: all required files present")

    # --- Step 8.1: Extract pocket & ligand ---
    extract_script = pathb_root / "scripts" / "02_Step2_因子实验" / "extract_pocket_ligand.py"
    cmd_81 = [
        python, str(extract_script),
        "--mapping_csv", str(paths["mapping_csv"]),
        "--pdb_dir", str(paths["pdb_dir"]),
        "--output_dir", str(str_tmp),
        "--pocket_radius", str(pocket_radius),
        "--experiment_name", exp_name,
    ]
    if include_heme:
        cmd_81.append("--include_heme")

    if not run_step("Step 8.1 Extract", cmd_81):
        return False

    # --- Step 8.2: Align ligand ---
    # Use PathA's align script (no structural parameter changes needed)
    align_script = paths["patha_scripts"] / "08_Step8_结构特征生成" / "step8_align_ligand.py"
    if align_script.exists():
        # The align script uses hardcoded paths, so we need a PathB wrapper
        # For now, copy the raw_ligand files and run alignment
        logger.info("Step 8.2: Ligand alignment (using PathA script as reference)")
        # TODO: Create PathB version of align_ligand.py with configurable paths
        # For the factorial experiment, alignment is the same regardless of
        # pocket_radius/include_heme, so we can share alignment results
        logger.warning("Step 8.2: Skipping alignment (will be addressed in Step 2 execution)")
    else:
        logger.warning("Step 8.2: align_ligand.py not found, skipping")

    # --- Step 8.3: Generate structure LMDB ---
    # TODO: Create PathB version with configurable paths and HETATM support
    logger.warning("Step 8.3: generate_structure_lmdb not yet parameterized for PathB")

    # --- Step 9: Inference ---
    # TODO: Create PathB version with configurable feature paths
    logger.warning("Step 9: inference not yet parameterized for PathB")

    # --- Step 10: Analysis ---
    # TODO: Create PathB version with configurable input paths
    logger.warning("Step 10: analysis not yet parameterized for PathB")

    logger.info("Experiment %s: Step 8.1 completed. Steps 8.2-10 pending parameterization.", exp_name)
    # Partial success: only Step 8.1 is implemented so far
    return True  # Step 8.1 succeeded; caller should note TODO steps


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="PathB experiment runner")
    p.add_argument("--name", help="Single experiment name (e.g. EXP01_B6_10A_noHeme)")
    p.add_argument("--pocket_radius", type=float, default=10.0)
    p.add_argument("--include_heme", action="store_true")
    p.add_argument("--run_all", action="store_true", help="Run all experiments from config")
    p.add_argument("--config", help="Config YAML path (for --run_all)")
    p.add_argument("--dry_run", action="store_true", help="Print commands without executing")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    paths = resolve_paths()

    if args.run_all:
        config_path = Path(args.config) if args.config else paths["pathb_root"] / "configs" / "base.yaml"
        if not config_path.exists():
            logger.error("Config not found: %s", config_path)
            return 2
        with open(config_path, encoding="utf-8") as f:
            cfg = yaml.safe_load(f)

        matrix = cfg.get("factorial_matrix", [])
        if not matrix:
            logger.error("No factorial_matrix in config")
            return 2

        results = {}
        for exp in matrix:
            ok = run_single_experiment(
                exp_name=exp["name"],
                pocket_radius=exp["pocket_radius"],
                include_heme=exp["include_heme"],
                paths=paths,
                dry_run=args.dry_run,
            )
            results[exp["name"]] = ok

        logger.info("\n" + "=" * 60)
        logger.info("SUMMARY")
        for name, ok in results.items():
            logger.info("  %s: %s", name, "OK" if ok else "FAILED")
        return 0 if all(results.values()) else 1

    if not args.name:
        logger.error("Specify --name for single experiment or --run_all for batch")
        return 2

    ok = run_single_experiment(
        exp_name=args.name,
        pocket_radius=args.pocket_radius,
        include_heme=args.include_heme,
        paths=paths,
        dry_run=args.dry_run,
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
