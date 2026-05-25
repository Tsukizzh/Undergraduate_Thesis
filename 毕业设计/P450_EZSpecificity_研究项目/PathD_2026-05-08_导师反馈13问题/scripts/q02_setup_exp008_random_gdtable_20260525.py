#!/usr/bin/env python3
"""Set up Q2 EXP008 random GDTable baseline on the PathD server.

The script creates a new experiment directory that reuses the existing PathD
random pt cache. It does not rebuild features or modify existing experiments.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path

import pandas as pd
import torch


SPLITS = ("train", "val", "test")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--root",
        default="/root/autodl-tmp/EZSpecificity/PathD/P450",
        help="PathD P450 root on the server.",
    )
    p.add_argument(
        "--source-exp",
        default="EXP007_strict_nn80_gdtable",
        help="Existing GDTable experiment to copy configs/src/scripts from.",
    )
    p.add_argument(
        "--exp-name",
        default="EXP008_random_gdtable",
        help="New experiment directory name.",
    )
    p.add_argument(
        "--data-name",
        default="exp008_random_baseline",
        help="Small documentation/audit directory name under data/q02_sequence_similarity_split.",
    )
    p.add_argument(
        "--cache",
        default=None,
        help="Existing random pt cache. Defaults to PathD base_from_PathC random cache.",
    )
    p.add_argument("--force", action="store_true", help="Allow replacing only this script's generated README/manifest files.")
    return p.parse_args()


def copy_tree(src: Path, dst: Path) -> None:
    if dst.exists():
        raise FileExistsError(f"Refusing to overwrite existing directory: {dst}")
    shutil.copytree(src, dst)


def load_index_count(cache: Path, split: str) -> int:
    index_path = cache / split / "index.pt"
    obj = torch.load(index_path, map_location="cpu")
    return len(obj["sample_ids"])


def split_label_stats(samples_csv: Path) -> dict[str, dict[str, int | float]]:
    df = pd.read_csv(samples_csv)
    out: dict[str, dict[str, int | float]] = {}
    for split in SPLITS:
        sub = df[df["split"] == split]
        n = int(len(sub))
        pos = int(sub["label"].sum())
        neg = int(n - pos)
        out[split] = {
            "n_samples": n,
            "n_positive": pos,
            "n_negative": neg,
            "positive_rate": round(pos / n, 6) if n else 0.0,
            "n_enzymes": int(sub["enzyme_index"].nunique()),
            "n_substrates": int(sub["substrate_index"].nunique()),
        }
    return out


def write_run_script(exp_dir: Path, cache: Path, exp_name: str) -> None:
    script = exp_dir / "scripts" / "run_train_gdtable.sh"
    text = f"""#!/bin/bash
# Q2 {exp_name}: PathD random split GDTable 1-GPU baseline.

set -euo pipefail

ulimit -n 65536 2>/dev/null || true
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${{PATH}}

EXP={exp_dir}
CACHE={cache}

BATCH_SIZE=${{BATCH_SIZE:-88}}
MAX_EPOCHS=${{MAX_EPOCHS:-200}}
NUM_WORKERS=${{NUM_WORKERS:-6}}
RUN_NAME=${{RUN_NAME:-Q2_{exp_name}_b88_full}}
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
    script = exp_dir / "scripts" / "run_train.sh"
    text = """#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "${SCRIPT_DIR}/run_train_gdtable.sh" "$@"
"""
    script.write_text(text, encoding="utf-8")
    os.chmod(script, 0o755)


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    q2_exp_root = root / "experiments/q02_sequence_similarity_split"
    q2_data_root = root / "data/q02_sequence_similarity_split"
    source_exp = q2_exp_root / args.source_exp
    exp_dir = q2_exp_root / args.exp_name
    data_dir = q2_data_root / args.data_name
    cache = (
        Path(args.cache).resolve()
        if args.cache
        else root / "data/base_from_PathC/cache_best_baseline/pt_cache_allfix_unified/random"
    )
    actual_used = root / "data/actual_used_baseline/tables/samples_actual_used.csv"

    required = [
        source_exp / "configs",
        source_exp / "scripts",
        source_exp / "src",
        cache / "manifest.pt",
        cache / "train/index.pt",
        cache / "val/index.pt",
        cache / "test/index.pt",
        actual_used,
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required inputs:\n" + "\n".join(missing))
    if exp_dir.exists():
        raise FileExistsError(f"Experiment directory already exists: {exp_dir}")

    exp_dir.mkdir(parents=True)
    copy_tree(source_exp / "configs", exp_dir / "configs")
    copy_tree(source_exp / "scripts", exp_dir / "scripts")
    copy_tree(source_exp / "src", exp_dir / "src")
    (exp_dir / "logs").mkdir()
    (exp_dir / "results/checkpoints").mkdir(parents=True)
    (exp_dir / "results/archive").mkdir(parents=True)

    write_run_script(exp_dir, cache, args.exp_name)
    write_forward_script(exp_dir)

    data_dir.mkdir(parents=True, exist_ok=True)
    (data_dir / "manifests").mkdir(exist_ok=True)
    (data_dir / "reports").mkdir(exist_ok=True)

    index_counts = {split: load_index_count(cache, split) for split in SPLITS}
    label_stats = split_label_stats(actual_used)
    mismatches = {
        split: {"index_count": index_counts[split], "table_count": label_stats[split]["n_samples"]}
        for split in SPLITS
        if index_counts[split] != label_stats[split]["n_samples"]
    }
    if mismatches:
        raise RuntimeError(f"Random cache index counts do not match actual_used table: {mismatches}")

    manifest = {
        "experiment": args.exp_name,
        "purpose": "PathD GDTable random split baseline for Q2 comparison.",
        "created_by": "q02_setup_exp008_random_gdtable_20260525.py",
        "source_experiment": str(source_exp),
        "experiment_dir": str(exp_dir),
        "random_cache": str(cache),
        "actual_used_samples_csv": str(actual_used),
        "split_stats": label_stats,
        "training_defaults": {
            "batch_size": 88,
            "max_epochs": 200,
            "num_workers": 6,
            "devices": 1,
            "shutdown": False,
            "gdtable": True,
            "gdtable_dense_dtype": "fp16",
            "gdtable_graph_fp16": True,
        },
    }

    (data_dir / "manifests/exp008_random_baseline_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    rows = [
        "# EXP008 random baseline audit",
        "",
        "EXP008 reuses the existing PathD random pt cache. It does not rebuild features or create a new pt cache.",
        "",
        f"- experiment dir: `{exp_dir}`",
        f"- random cache: `{cache}`",
        f"- source GDTable code/config: `{source_exp}`",
        "",
        "| split | samples | enzymes | substrates | positives | negatives | positive rate |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for split in SPLITS:
        st = label_stats[split]
        rows.append(
            f"| {split} | {st['n_samples']} | {st['n_enzymes']} | {st['n_substrates']} | "
            f"{st['n_positive']} | {st['n_negative']} | {st['positive_rate']:.6f} |"
        )
    rows += [
        "",
        "Validation:",
        "",
        "- random cache index counts match `samples_actual_used.csv` split counts.",
        "- original random cache is read-only input; EXP008 writes only its own experiment directory and this small audit directory.",
    ]
    (data_dir / "reports/exp008_random_baseline_audit.md").write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )

    readme = [
        "# EXP008 random baseline",
        "",
        "This directory documents the PathD GDTable random split baseline used by Q2 EXP008.",
        "",
        "The actual pt cache is reused from:",
        "",
        f"`{cache}`",
        "",
        "Generated files here are manifests and reports only.",
    ]
    (data_dir / "README.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
