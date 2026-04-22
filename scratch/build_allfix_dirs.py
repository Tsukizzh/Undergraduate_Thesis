"""Phase 4: Build 6 _allfix cache directories with symlinks.

For each name in:
    bare_allfix, bare_allfix_unified,
    heme_allfix, heme_allfix_unified,
    geom_allfix, geom_allfix_unified,

Create directory structure matching its corresponding base _fixed cache layout:
    - bare/heme have random/ sub-level (random/{train,val,test})
    - geom is flat ({train,val,test})

Each _allfix cache contains:
    - manifest.pt      (copy from _fixed)
    - enzymes/         symlink -> _fixed's enzymes dir
    - substrates/      symlink -> pt_cache_allfix_shared/substrates (shared)
    - {train,val,test}/
        - index.pt      (from Phase 3 staging)
        - samples/ OR graph_*.pt  symlink(s) -> base cache heavy data

Codex-reviewed round 7.
"""
import os
import shutil
import torch
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
STAGING = Path("/tmp/allfix_indices")
SPLITS = ['train', 'val', 'test']
SHARED_SUBSTRATES = BASE / "data/pt_cache_allfix_shared/substrates"

# Each config:
#   base_cache:    root of base cache's heavy data (graph_*.pt or samples/)
#   fixed_cache:   root of _fixed cache (source for manifest.pt and enzymes/)
#   has_random:    True if layout is random/{train,val,test}/, False if flat
#   mode:          'per_sample' (symlink samples/) or 'shard' (symlink graph_*.pt)
CACHES = {
    'bare': {
        'base':      BASE / "data/pt_cache/random",
        'fixed':     BASE / "data/pt_cache_fixed/random",
        'has_random': True,
        'mode':       'per_sample',
    },
    'heme': {
        'base':      BASE / "data/pt_cache_heme/random",
        'fixed':     BASE / "data/pt_cache_heme_fixed/random",
        'has_random': True,
        'mode':       'shard',
    },
    'geom': {
        'base':      BASE / "data/pt_cache_geom",
        'fixed':     BASE / "data/pt_cache_geom_fixed",
        'has_random': False,
        'mode':       'shard',
    },
}


def build_one_allfix(kind: str, variant: str) -> None:
    """kind in {'bare','heme','geom'}; variant in {'natural','unified'}."""
    cfg = CACHES[kind]
    variant_suffix = "_allfix" if variant == 'natural' else "_allfix_unified"

    if kind == 'bare':
        target_root = BASE / f"data/pt_cache{variant_suffix}"
    elif kind == 'heme':
        target_root = BASE / f"data/pt_cache_heme{variant_suffix}"
    else:  # geom
        target_root = BASE / f"data/pt_cache_geom{variant_suffix}"

    if cfg['has_random']:
        target_inner = target_root / "random"
    else:
        target_inner = target_root

    # Wipe and recreate target
    if target_root.exists():
        print(f"  removing existing {target_root}")
        shutil.rmtree(target_root)
    target_inner.mkdir(parents=True)

    # 1. manifest.pt (copy)
    src_manifest = cfg['fixed'] / "manifest.pt"
    assert src_manifest.exists(), f"missing {src_manifest}"
    shutil.copy2(src_manifest, target_inner / "manifest.pt")

    # 2. enzymes/ symlink
    src_enzymes = cfg['fixed'] / "enzymes"
    assert src_enzymes.exists(), f"missing {src_enzymes}"
    os.symlink(src_enzymes, target_inner / "enzymes")

    # 3. substrates/ symlink -> shared
    assert SHARED_SUBSTRATES.exists()
    os.symlink(SHARED_SUBSTRATES, target_inner / "substrates")

    # 4. For each split: index.pt + heavy data symlinks
    for split in SPLITS:
        split_dir = target_inner / split
        split_dir.mkdir()

        # 4a. index.pt from staging
        staging_name = f"{kind}_natural" if variant == 'natural' else f"{kind}_unified"
        src_index = STAGING / staging_name / split / "index.pt"
        assert src_index.exists(), f"missing staging {src_index}"
        shutil.copy2(src_index, split_dir / "index.pt")

        # 4b. heavy data symlinks
        base_split = cfg['base'] / split
        if cfg['mode'] == 'per_sample':
            # symlink entire samples/ dir
            src_samples = base_split / "samples"
            assert src_samples.exists() and src_samples.is_dir(), f"missing {src_samples}"
            os.symlink(src_samples, split_dir / "samples")
        else:
            # symlink all graph_*.pt shards
            shards = sorted(base_split.glob("graph_*.pt"))
            assert len(shards) > 0, f"no shards in {base_split}"
            for sh in shards:
                os.symlink(sh, split_dir / sh.name)


print("=== Phase 4: Building 6 _allfix cache directories ===\n")

variants_to_build = [('bare', 'natural'), ('heme', 'natural'), ('geom', 'natural'),
                     ('bare', 'unified'), ('heme', 'unified'), ('geom', 'unified')]

for kind, variant in variants_to_build:
    print(f"\n--- {kind}_{variant} ---")
    build_one_allfix(kind, variant)

    # Show what got built
    variant_suffix = "_allfix" if variant == 'natural' else "_allfix_unified"
    if kind == 'bare':
        root = BASE / f"data/pt_cache{variant_suffix}"
    elif kind == 'heme':
        root = BASE / f"data/pt_cache_heme{variant_suffix}"
    else:
        root = BASE / f"data/pt_cache_geom{variant_suffix}"

    # Walk and report
    inner = root / "random" if CACHES[kind]['has_random'] else root
    files = sorted(inner.iterdir())
    for f in files:
        if f.is_symlink():
            target = os.readlink(f)
            print(f"  {f.name} -> {target}")
        elif f.is_dir():
            subfiles = list(f.iterdir())
            n_symlinks = sum(1 for s in subfiles if s.is_symlink())
            n_files = sum(1 for s in subfiles if not s.is_symlink())
            print(f"  {f.name}/: {n_files} files + {n_symlinks} symlinks")
        else:
            print(f"  {f.name}: file ({f.stat().st_size} bytes)")

    # Verify all symlinks are valid (targets exist)
    broken = []
    for p in inner.rglob("*"):
        if p.is_symlink() and not p.exists():
            broken.append(p)
    assert not broken, f"broken symlinks: {broken}"
    print(f"  [OK] all symlinks valid")

print("\n=== Phase 4 DONE ===")
print("\nFinal directory sizes:")
import subprocess
for v in ['pt_cache_allfix', 'pt_cache_heme_allfix', 'pt_cache_geom_allfix',
          'pt_cache_allfix_unified', 'pt_cache_heme_allfix_unified', 'pt_cache_geom_allfix_unified']:
    d = BASE / "data" / v
    if d.exists():
        r = subprocess.run(['du', '-sh', str(d)], capture_output=True, text=True)
        print(f"  {r.stdout.strip()}")
