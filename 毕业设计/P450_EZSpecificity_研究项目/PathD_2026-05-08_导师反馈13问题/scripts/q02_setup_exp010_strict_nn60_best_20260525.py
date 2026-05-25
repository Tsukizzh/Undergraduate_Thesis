#!/usr/bin/env python3
"""Set up Q2 EXP010 from one audited EXP009 strict NN60 candidate.

The script copies only the selected candidate split metadata, builds a fresh
pt-cache overlay from the existing PathD baseline cache, then creates an
independent GDTable training directory. Existing outputs are never replaced.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", default="/root/autodl-tmp/EZSpecificity/PathD/P450")
    p.add_argument("--candidate-dir", required=True, help="Selected EXP009 candidate split directory.")
    p.add_argument("--candidate-name", default=None, help="Optional name to store in manifests.")
    p.add_argument("--source-exp", default="EXP008_random_gdtable_b80_retry_20260525")
    p.add_argument("--data-name", default="exp010_strict_nn60_best")
    p.add_argument("--cache-name", default="strict_nn60_best_main")
    p.add_argument("--exp-name", default="EXP010_strict_nn60_best_gdtable")
    p.add_argument("--batch-size", type=int, default=80)
    p.add_argument("--max-epochs", type=int, default=200)
    p.add_argument("--num-workers", type=int, default=6)
    p.add_argument("--file-mode", choices=("hardlink", "copy", "symlink"), default="hardlink")
    return p.parse_args()


def ensure_abs_under(path: Path, root: Path) -> Path:
    path = path.resolve()
    root = root.resolve()
    if path != root and root not in path.parents:
        raise RuntimeError(f"Refusing to write outside PathD root: {path}")
    return path


def ensure_missing(path: Path) -> None:
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing path: {path}")


def ensure_exists(path: Path, kind: str) -> None:
    if kind == "dir" and not path.is_dir():
        raise FileNotFoundError(f"Missing directory: {path}")
    if kind == "file" and not path.is_file():
        raise FileNotFoundError(f"Missing file: {path}")


def copy_tree(src: Path, dst: Path) -> None:
    ensure_exists(src, "dir")
    ensure_missing(dst)
    shutil.copytree(src, dst)


def write_run_script(exp_dir: Path, cache_dir: Path, exp_name: str, batch_size: int, max_epochs: int, num_workers: int) -> None:
    script = exp_dir / "scripts/run_train_gdtable.sh"
    text = f"""#!/bin/bash
# Q2 {exp_name}: strict NN60 selected split GDTable training.

set -euo pipefail

ulimit -n 65536 2>/dev/null || true
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${{PATH}}

EXP={exp_dir}
CACHE={cache_dir}

BATCH_SIZE=${{BATCH_SIZE:-{batch_size}}}
MAX_EPOCHS=${{MAX_EPOCHS:-{max_epochs}}}
NUM_WORKERS=${{NUM_WORKERS:-{num_workers}}}
RUN_NAME=${{RUN_NAME:-Q2_{exp_name}_b{batch_size}_full}}
SHUTDOWN=${{SHUTDOWN:-false}}

ARGS=(
  --config ${{EXP}}/configs/config.yml
  --cache-dir ${{CACHE}}
  --edge-mode fixed
  --batch-size ${{BATCH_SIZE}}
  --max-epochs ${{MAX_EPOCHS}}
  --devices 1
  --num-workers ${{NUM_WORKERS}}
  --run-name ${{RUN_NAME}}
  --gdtable
  --gdtable-dense-dtype fp16
  --gdtable-graph-fp16
  --train-in-order false
)

case "${{SHUTDOWN}}" in
  true|TRUE|1|yes|YES)
    ARGS+=(--shutdown)
    ;;
esac

python ${{EXP}}/scripts/main_training_gdtable.py "${{ARGS[@]}}"
"""
    script.write_text(text, encoding="utf-8")
    os.chmod(script, 0o755)


def write_forward_script(exp_dir: Path) -> None:
    script = exp_dir / "scripts/run_train.sh"
    text = """#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/run_train_gdtable.sh" "$@"
"""
    script.write_text(text, encoding="utf-8")
    os.chmod(script, 0o755)


def run_build_cache(root: Path, candidate_copy: Path, output_root: Path, cache_name: str, file_mode: str) -> None:
    script = root / "scripts/q02_build_pt_cache_overlay_20260522.py"
    ensure_exists(script, "file")
    cmd = [
        sys.executable,
        str(script),
        "--root",
        str(root),
        "--threshold",
        "strict_nn60",
        "--split-dir",
        str(candidate_copy),
        "--output-root",
        str(output_root),
        "--name",
        cache_name,
        "--file-mode",
        file_mode,
        "--check-content",
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    candidate_dir = Path(args.candidate_dir).resolve()
    q2_data_root = root / "data/q02_sequence_similarity_split"
    q2_exp_root = root / "experiments/q02_sequence_similarity_split"
    data_dir = ensure_abs_under(q2_data_root / args.data_name, root)
    exp_dir = ensure_abs_under(q2_exp_root / args.exp_name, root)
    source_exp = q2_exp_root / args.source_exp
    cache_root = data_dir / "pt_cache"
    cache_dir = cache_root / args.cache_name
    candidate_name = args.candidate_name or candidate_dir.name
    candidate_copy = data_dir / "splits" / candidate_name

    for required in [
        candidate_dir,
        candidate_dir / "train_samples_strict_nn60.csv",
        candidate_dir / "val_samples_strict_nn60.csv",
        candidate_dir / "test_samples_strict_nn60.csv",
        candidate_dir / "split_summary.json",
        source_exp / "configs",
        source_exp / "scripts",
        source_exp / "src",
    ]:
        ensure_exists(required, "dir" if required.suffix == "" else "file")

    ensure_missing(data_dir)
    ensure_missing(exp_dir)
    summary = json.loads((candidate_dir / "split_summary.json").read_text(encoding="utf-8"))
    if not summary.get("all_passes"):
        raise RuntimeError(f"Candidate did not pass EXP009 audits: {candidate_dir}")

    data_dir.mkdir(parents=True)
    (data_dir / "splits").mkdir()
    (data_dir / "manifests").mkdir()
    (data_dir / "reports").mkdir()
    copy_tree(candidate_dir, candidate_copy)

    run_build_cache(
        root=root,
        candidate_copy=candidate_copy,
        output_root=cache_root,
        cache_name=args.cache_name,
        file_mode=args.file_mode,
    )

    exp_dir.mkdir(parents=True)
    copy_tree(source_exp / "configs", exp_dir / "configs")
    copy_tree(source_exp / "scripts", exp_dir / "scripts")
    copy_tree(source_exp / "src", exp_dir / "src")
    (exp_dir / "logs").mkdir()
    (exp_dir / "results/checkpoints").mkdir(parents=True)
    (exp_dir / "results/archive").mkdir(parents=True)
    write_run_script(exp_dir, cache_dir, args.exp_name, args.batch_size, args.max_epochs, args.num_workers)
    write_forward_script(exp_dir)

    manifest = {
        "experiment": args.exp_name,
        "created_by": "q02_setup_exp010_strict_nn60_best_20260525.py",
        "purpose": "Train one selected strict NN60 split from EXP009.",
        "selected_candidate_name": candidate_name,
        "selected_candidate_source": str(candidate_dir),
        "selected_candidate_copy": str(candidate_copy),
        "source_experiment": str(source_exp),
        "data_dir": str(data_dir),
        "pt_cache": str(cache_dir),
        "experiment_dir": str(exp_dir),
        "training_defaults": {
            "batch_size": args.batch_size,
            "max_epochs": args.max_epochs,
            "num_workers": args.num_workers,
            "devices": 1,
            "shutdown": False,
            "gdtable": True,
            "gdtable_dense_dtype": "fp16",
            "gdtable_graph_fp16": True,
        },
        "candidate_summary": summary,
    }
    (data_dir / "manifests/exp010_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    (data_dir / "README.md").write_text(
        "# EXP010 strict NN60 selected split\n\n"
        f"- selected candidate: `{candidate_name}`\n"
        f"- copied split dir: `{candidate_copy}`\n"
        f"- pt cache: `{cache_dir}`\n"
        f"- experiment dir: `{exp_dir}`\n"
        "\nThe pt cache is an overlay from the existing PathD baseline cache. Feature files are not recomputed.\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
