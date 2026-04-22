"""
Build the paperfilter cache overlay on the server.

Non-destructive: all paths except test/index.pt are symlinks to the source
pt_cache_allfix_unified/random. The new test/index.pt has rows removed where
enzyme_ids is in the blacklist.

Run on server:
  cd /root/autodl-tmp/EZSpecificity/PathC/P450
  /root/miniconda3/bin/python /root/autodl-tmp/.../build_paperfilter_cache.py
"""
from __future__ import annotations
import json, os
from pathlib import Path

import torch

# Paths (hard-coded for server, intentional)
ROOT = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
SRC_CACHE = ROOT / "data/pt_cache_allfix_unified/random"
DST_CACHE = ROOT / "data/pt_cache_allfix_unified_paperfilter/random"
BLACKLIST_JSON = ROOT / "experiments/EXP004_paper_baseline_unified_workdir/paper_blacklist.json"
REPORT_OUT = ROOT / "experiments/EXP004_paper_baseline_unified_workdir/filter_report.json"


def symlink(src: Path, dst: Path):
    """Create a symlink dst -> src (absolute). Skip if exists and correct."""
    assert src.exists(), f"source does not exist: {src}"
    if dst.is_symlink() or dst.exists():
        existing = os.readlink(dst) if dst.is_symlink() else None
        if existing == str(src):
            print(f"  skip (exists): {dst.name}")
            return
        raise RuntimeError(f"destination exists but not symlink to source: {dst}")
    dst.symlink_to(src, target_is_directory=src.is_dir())
    print(f"  symlink: {dst.name} -> {src}")


def main():
    assert SRC_CACHE.exists(), f"source cache missing: {SRC_CACHE}"
    assert BLACKLIST_JSON.exists(), f"blacklist json missing: {BLACKLIST_JSON}"

    # Load blacklist
    with open(BLACKLIST_JSON) as f:
        bl = json.load(f)
    blacklist = set(bl["blacklisted_enzyme_global_ids"])
    print(f"Blacklist: {len(blacklist)} enzyme_global_ids")

    # Create destination root
    DST_CACHE.mkdir(parents=True, exist_ok=True)
    print(f"\nBuilding destination: {DST_CACHE}")

    # Symlink shared dirs + unused splits (keep structure complete)
    print("\n[1/3] Symlinking shared dirs and unused splits:")
    for name in ["enzymes", "substrates", "manifest.pt", "train", "val"]:
        symlink(SRC_CACHE / name, DST_CACHE / name)

    # Build new test/ with filtered index.pt
    print("\n[2/3] Building filtered test/:")
    test_dst = DST_CACHE / "test"
    test_dst.mkdir(exist_ok=True)
    symlink(SRC_CACHE / "test" / "samples", test_dst / "samples")

    # Load source test index
    src_idx_path = SRC_CACHE / "test" / "index.pt"
    src_idx = torch.load(src_idx_path, map_location="cpu", weights_only=False)
    print(f"  Source test index keys: {list(src_idx.keys())}")
    src_n = len(src_idx["sample_ids"])
    print(f"  Source test samples: {src_n}")

    # Verify all arrays same length
    for k, v in src_idx.items():
        assert len(v) == src_n, f"array length mismatch: {k} = {len(v)} vs {src_n}"
        assert hasattr(v, "shape"), f"expected tensor for {k}, got {type(v)}"

    # Build mask
    enzyme_ids_src = src_idx["enzyme_ids"]
    blacklist_t = torch.tensor(sorted(blacklist), dtype=enzyme_ids_src.dtype)
    hit_mask = torch.isin(enzyme_ids_src, blacklist_t)
    keep_mask = ~hit_mask

    # Apply mask to every field (positional alignment preserved)
    new_idx = {k: v[keep_mask].clone() for k, v in src_idx.items()}

    # Stats
    n_src = int(src_n)
    n_dropped = int(hit_mask.sum())
    n_kept = int(keep_mask.sum())
    assert n_dropped + n_kept == n_src
    unique_src_enz = set(enzyme_ids_src.tolist())
    unique_kept_enz = set(new_idx["enzyme_ids"].tolist())
    unique_dropped_enz = unique_src_enz - unique_kept_enz

    # Label breakdown (sample_id -> need to load samples to get label; skip here,
    # let verify script compute it. We only have the index-level fields here.)

    print(f"\n  Source samples:  {n_src}")
    print(f"  Dropped (in blacklist enzymes): {n_dropped} ({100*n_dropped/n_src:.1f}%)")
    print(f"  Kept:            {n_kept} ({100*n_kept/n_src:.1f}%)")
    print(f"  Unique enzymes src:     {len(unique_src_enz)}")
    print(f"  Unique enzymes dropped: {len(unique_dropped_enz)}")
    print(f"  Unique enzymes kept:    {len(unique_kept_enz)}")

    # Save
    dst_idx_path = test_dst / "index.pt"
    torch.save(new_idx, dst_idx_path)
    print(f"\n  Wrote: {dst_idx_path}")

    # Save report
    print("\n[3/3] Saving report:")
    report = {
        "src_cache": str(SRC_CACHE),
        "dst_cache": str(DST_CACHE),
        "blacklist_enzyme_count": len(blacklist),
        "src_test_samples": n_src,
        "dropped_samples": n_dropped,
        "kept_samples": n_kept,
        "dropped_fraction": round(n_dropped / n_src, 4),
        "src_unique_enzymes_in_test": len(unique_src_enz),
        "dropped_unique_enzymes": len(unique_dropped_enz),
        "kept_unique_enzymes": len(unique_kept_enz),
        "index_keys_verified": list(new_idx.keys()),
        "new_index_array_lengths": {k: len(v) for k, v in new_idx.items()},
    }
    with open(REPORT_OUT, "w") as f:
        json.dump(report, f, indent=2)
    print(f"  {REPORT_OUT}")

    print("\n" + "=" * 70)
    print("DONE. Non-destructive overlay cache ready:")
    print(f"  {DST_CACHE}")
    print("=" * 70)


if __name__ == "__main__":
    main()
