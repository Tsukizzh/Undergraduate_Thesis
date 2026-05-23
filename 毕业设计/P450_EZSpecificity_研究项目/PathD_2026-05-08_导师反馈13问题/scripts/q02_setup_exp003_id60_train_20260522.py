#!/usr/bin/env python3
"""Create the PathD Q2 EXP003 id60 baseline training directory.

This script only writes inside the new Q2 experiment directory. It refuses to
run if the destination already exists, so earlier results are not overwritten.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path


ROOT = Path("/root/autodl-tmp/EZSpecificity/PathD/P450")
BASELINE = ROOT / "baselines/EXP001_allfix_unified_best"
DEST = ROOT / "experiments/q02_sequence_similarity_split/EXP003_id60_baseline_train"
CACHE = (
    ROOT
    / "data/q02_sequence_similarity_split/exp002_actual_used_cache_valid/pt_cache/id60_main"
)


RUN_TRAIN = f"""#!/bin/bash
# Q2 EXP003: train the original baseline architecture on the actual-used id60 split.
# This launcher writes only inside EXP003_id60_baseline_train and does not shut down the server.

set -euo pipefail

ulimit -n 65536 2>/dev/null || true
export OMP_NUM_THREADS=1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export TORCH_NCCL_ENABLE_EAGER_CONNECT=0
export PATH=/root/miniconda3/bin:${{PATH}}

EXP={DEST}
CACHE={CACHE}

python ${{EXP}}/scripts/main_training_pt.py \\
  --config ${{EXP}}/configs/config.yml \\
  --cache-dir ${{CACHE}} \\
  --edge-mode fixed \\
  --batch-size 88 \\
  --max-epochs 200 \\
  --devices 1 \\
  --num-workers 6 \\
  --run-name Q2_EXP003_id60_baseline_train
"""


MANIFEST = f"""# Q2 EXP003 id60 baseline train

Purpose: train the original PathC/PathD baseline architecture on the Q2 EXP002
actual-used id60 sequence-similarity split.

## Source

- Baseline code/config source: `{BASELINE}`
- Training cache: `{CACHE}`
- Split family: EXP002 actual-used, id60 cluster-held-out split

## Output boundary

All training logs, checkpoints, metrics and test results should stay under:

`{DEST}`

This directory is created as an isolated PathD experiment. It should not write
into `baselines/EXP001_allfix_unified_best`.

## Main command

`bash scripts/run_train.sh`

The launcher uses one GPU, batch size 88, `edge-mode=fixed`, and does not pass
`--shutdown`.
"""


def ignore_scripts(_dir: str, names: list[str]) -> set[str]:
    ignored: set[str] = set()
    for name in names:
        if name == "__pycache__" or name.endswith(".pyc"):
            ignored.add(name)
        if name.startswith("run_train") and name.endswith(".sh"):
            ignored.add(name)
    return ignored


def main() -> None:
    if not BASELINE.is_dir():
        raise SystemExit(f"Missing baseline directory: {BASELINE}")
    if not CACHE.is_dir():
        raise SystemExit(f"Missing id60 cache: {CACHE}")
    if DEST.exists():
        raise SystemExit(f"Destination already exists, refusing to overwrite: {DEST}")

    DEST.mkdir(parents=True)
    for dirname in ("scripts", "configs", "src"):
        shutil.copytree(
            BASELINE / dirname,
            DEST / dirname,
            ignore=ignore_scripts,
        )

    (DEST / "logs").mkdir()
    (DEST / "results").mkdir()
    (DEST / "manifest.md").write_text(MANIFEST, encoding="utf-8")

    run_path = DEST / "scripts/run_train.sh"
    run_path.write_text(RUN_TRAIN, encoding="utf-8")
    os.chmod(run_path, 0o755)

    print(f"CREATED {DEST}")
    print(f"RUN_SCRIPT {run_path}")
    print(f"CACHE {CACHE}")


if __name__ == "__main__":
    main()
