"""Phase 2: Build shared substrates flatbin from grover_fingerprint_fixed.lmdb.

Reuses `build_substrate_shards()` and `convert_substrate_shards_to_flatbin()`
from the existing build_pt_cache.py in EXP003_fixed.

Output: /root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_shared/substrates/
    grover_0000.pt               (shard files, raw output)
    index.pt                     (shard index, raw output)
    substrates_grover.bin        (fp16 flatbin, consumed by training)
    substrates_meta.pt           (GROVER mean + Morgan, consumed by training)
    substrates_index.pt          (substrate_id -> offset)

Codex-reviewed (session 019d88a7, round 5). Config format:
    grover_path and morgan_path must be LISTS (ds_id = list index).
"""
import os
import shutil
import sys
from pathlib import Path

# Import build_pt_cache from EXP003_fixed (latest version)
BUILD_PT_CACHE_DIR = "/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP003_fixed/scripts"
sys.path.insert(0, BUILD_PT_CACHE_DIR)

# Also ensure src/ is on path (build_pt_cache may import from Datasets/Models)
SRC_DIR = "/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP003_fixed/src"
sys.path.insert(0, SRC_DIR)

from build_pt_cache import build_substrate_shards, convert_substrate_shards_to_flatbin

# --- Paths ---
BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
GROVER_FIXED = BASE / "data/features/grover_fingerprint_fixed.lmdb"
MORGAN = BASE / "data/features/morgan_fingerprint.npy"
OUT_DIR = BASE / "data/pt_cache_allfix_shared"

assert GROVER_FIXED.exists(), f"missing {GROVER_FIXED}"
assert MORGAN.exists(), f"missing {MORGAN}"

# --- Clean target dir to avoid stale shard skip ---
substrates_dir = OUT_DIR / "substrates"
if substrates_dir.exists():
    print(f"Removing existing {substrates_dir}")
    shutil.rmtree(substrates_dir)

OUT_DIR.mkdir(parents=True, exist_ok=True)

# --- Build config (lists, not strings) ---
config = {
    "data": {
        "grover_path": [str(GROVER_FIXED)],
        "morgan_path": [str(MORGAN)],
    }
}

# --- Phase 2a: shard build ---
print("\n=== Phase 2a: build_substrate_shards ===")
print(f"  output_dir: {OUT_DIR}")
print(f"  grover: {GROVER_FIXED}")
print(f"  morgan: {MORGAN}")
print(f"  shard_size: 4096 (2124 substrates -> 1 shard)")

build_substrate_shards(
    output_dir=OUT_DIR,
    config=config,
    shard_size=4096,
    num_workers=0,   # single-process to avoid surprises
)

print("\nShard files after build:")
for p in sorted(substrates_dir.iterdir()):
    sz = p.stat().st_size
    print(f"  {p.name}: {sz/1024:.1f} KB")

# --- Phase 2b: flatbin conversion ---
print("\n=== Phase 2b: convert_substrate_shards_to_flatbin ===")
convert_substrate_shards_to_flatbin(
    output_dir=OUT_DIR,
    overwrite=True,
)

print("\nFinal files:")
for p in sorted(substrates_dir.iterdir()):
    sz = p.stat().st_size
    print(f"  {p.name}: {sz/1024:.1f} KB")

print("\n=== Phase 2 DONE ===")
