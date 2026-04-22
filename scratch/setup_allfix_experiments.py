"""Build 6 _allfix experiment directories by copying from _fixed and patching paths.

For each (source _fixed, target _allfix[_unified]):
1. cp -a source target  (preserves symlinks, ownership, etc)
2. rm -rf target/logs target/results  (stale run data)
3. Clean __pycache__ everywhere
4. Patch scripts/run_train.sh: EXP path, CACHE path, --run-name, comment header
5. Patch configs/config.yml: data.tag
6. Verify diff vs source
"""
import os
import shutil
import subprocess
import re
from pathlib import Path

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/experiments")

# (source_fixed, target_allfix, target_cache_rel, tag)
SETUPS = [
    # Natural group
    ('EXP001_fixed', 'EXP001_allfix',
     'pt_cache_allfix/random', 'EXP001-p450-allfix-natural'),
    ('EXP002a_fixed', 'EXP002a_allfix',
     'pt_cache_heme_allfix/random', 'EXP002a-p450-allfix-natural'),
    ('EXP003_fixed', 'EXP003_allfix',
     'pt_cache_geom_allfix', 'EXP003-p450-allfix-natural'),
    # Unified group
    ('EXP001_fixed', 'EXP001_allfix_unified',
     'pt_cache_allfix_unified/random', 'EXP001-p450-allfix-unified'),
    ('EXP002a_fixed', 'EXP002a_allfix_unified',
     'pt_cache_heme_allfix_unified/random', 'EXP002a-p450-allfix-unified'),
    ('EXP003_fixed', 'EXP003_allfix_unified',
     'pt_cache_geom_allfix_unified', 'EXP003-p450-allfix-unified'),
]


def sh(cmd: str):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd failed: {cmd}\nstderr: {r.stderr}")
    return r.stdout


def patch_file(path: Path, replacements: list):
    """Apply list of (old_substr, new_substr) substitutions. old_substr must exist."""
    text = path.read_text()
    for old, new in replacements:
        assert old in text, f"pattern not found in {path}: {old!r}"
        text = text.replace(old, new)
    path.write_text(text)


for src_name, tgt_name, cache_rel, tag in SETUPS:
    src = BASE / src_name
    tgt = BASE / tgt_name
    print(f"\n=== {src_name} -> {tgt_name} ===")
    assert src.exists(), f"source missing: {src}"

    # 1. Remove existing target (safety)
    if tgt.exists():
        print(f"  removing existing {tgt}")
        shutil.rmtree(tgt)

    # 2. Copy (cp -a to preserve symlinks + dir perms)
    sh(f"cp -a '{src}' '{tgt}'")

    # 3. Clean stale logs/results and __pycache__
    for sub in ['logs', 'results']:
        p = tgt / sub
        if p.exists():
            shutil.rmtree(p)
    # Recursive __pycache__ cleanup
    sh(f"find '{tgt}' -type d -name __pycache__ -exec rm -rf {{}} + 2>/dev/null || true")

    # 4. Patch run_train.sh
    run_sh = tgt / "scripts/run_train.sh"
    assert run_sh.exists()

    src_cache_rel_map = {
        'EXP001_fixed': 'pt_cache_fixed/random',
        'EXP002a_fixed': 'pt_cache_heme_fixed/random',
        'EXP003_fixed': 'pt_cache_geom_fixed',
    }
    src_cache_rel = src_cache_rel_map[src_name]

    patch_file(run_sh, [
        (f"experiments/{src_name}", f"experiments/{tgt_name}"),
        (f"data/{src_cache_rel}", f"data/{cache_rel}"),
        (f"--run-name {src_name}", f"--run-name {tgt_name}"),
    ])

    # 5. Patch config.yml data.tag
    cfg = tgt / "configs/config.yml"
    assert cfg.exists()
    text = cfg.read_text()
    # Replace the `tag: ...` line under data:
    new_text, n = re.subn(r'tag:\s*[^\n]+', f'tag: {tag}', text, count=1)
    assert n == 1, f"tag line not found in {cfg}"
    cfg.write_text(new_text)

    # 6. Verify: diff scripts/run_train.sh briefly
    print(f"  run_train.sh lines (post-patch):")
    for line in run_sh.read_text().splitlines():
        if any(k in line for k in ['EXP=', 'CACHE=', '--run-name']):
            print(f"    {line}")
    print(f"  config.yml tag: {tag}")
    print(f"  stale logs/results removed")

    # Verify cache target exists
    data_root = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data")
    cache_path = data_root / cache_rel
    assert cache_path.exists(), f"cache target missing: {cache_path}"
    print(f"  cache verified at {cache_path}")

print("\n=== All 6 allfix experiment dirs built ===")
for src_name, tgt_name, _, _ in SETUPS:
    tgt = BASE / tgt_name
    size = sh(f"du -sh '{tgt}' | cut -f1").strip()
    print(f"  {tgt_name}: {size}")
