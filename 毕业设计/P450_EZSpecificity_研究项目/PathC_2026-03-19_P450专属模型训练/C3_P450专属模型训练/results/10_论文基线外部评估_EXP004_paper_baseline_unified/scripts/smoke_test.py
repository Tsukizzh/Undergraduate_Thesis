"""
Smoke test before full grid. Runs paper ckpt on 1 batch of filtered test,
dumps logit statistics, catches silent failure modes (NaN / near-constant /
output head misselection) BEFORE burning GPU time on 4 full runs.

Usage (on server, after opening GPU):
  cd /root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP004_paper_baseline_unified
  python scripts/smoke_test.py

Exits 0 on clean sanity, 1 on any red flag.
"""
from __future__ import annotations
import os, sys, json
from pathlib import Path

import numpy as np
import torch
from torch_geometric.loader import DataLoader

EXP_DIR = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP004_paper_baseline_unified")
CACHE_DIR = Path("/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_unified_paperfilter/random")
CKPT = EXP_DIR / "results/checkpoints/paper_best-checkpoint.ckpt"
CFG = EXP_DIR / "configs/config.yml"

sys.path.insert(0, str(EXP_DIR / "src"))
sys.path.insert(0, str(EXP_DIR / "scripts"))


def main():
    # Load config
    sys.path.insert(0, str(EXP_DIR / "src"))
    from utils import load_config
    config = load_config(str(CFG))

    # Dataset + loader
    from pt_dataset import PtCacheDataset
    print(f"[1/5] Building PtCacheDataset (legacy_bug edge mode)...")
    ds = PtCacheDataset(
        cache_dir=str(CACHE_DIR),
        split="test",
        edge_mode="legacy_bug",
        dist_noise=False,
        cutoff=config.transform.cutoff,
        num_r_gaussian=config.transform.num_r_gaussian,
        max_enzyme_len=config.data.max_enzyme_length,
        max_substrate_len=config.data.max_substrate_length,
        preload=False,
    )
    print(f"  dataset size: {len(ds)}")
    loader = DataLoader(
        ds, batch_size=88, shuffle=False, num_workers=2,
        follow_batch=["ligand_index"],
    )

    # Model + ckpt
    print(f"[2/5] Loading SS model + paper ckpt...")
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    from Models.ss import SS
    model = SS(config)
    ckpt = torch.load(str(CKPT), map_location="cpu", weights_only=False)
    sd = ckpt.get("state_dict", ckpt)
    res = model.load_state_dict(sd, strict=True)  # strict=True: we already proved this works
    print(f"  strict load OK")
    model.to(device).eval()

    # Forward pass on ONE batch
    print(f"[3/5] Running 1 batch forward...")
    batch = next(iter(loader))
    batch = batch.to(device, non_blocking=True)
    with torch.inference_mode():
        output = model(batch)
        logits = output[0] if isinstance(output, tuple) else output
        logits = logits.squeeze(-1).detach().cpu().numpy()
    labels = batch.label.detach().cpu().numpy()
    try:
        tags = batch.tag
        if isinstance(tags, torch.Tensor):
            tags = tags.cpu().numpy().tolist()
    except Exception:
        tags = ["?"] * len(labels)

    # Sanity stats
    print(f"\n[4/5] Sanity stats on batch")
    n = len(logits)
    n_pos = int((labels == 1).sum())
    n_neg = int((labels == 0).sum())
    finite = bool(np.isfinite(logits).all())
    l_min, l_max = float(logits.min()), float(logits.max())
    l_mean, l_std = float(logits.mean()), float(logits.std())

    print(f"  batch size    : {n}")
    print(f"  pos / neg     : {n_pos} / {n_neg}")
    print(f"  all finite    : {finite}")
    print(f"  logit min/max : {l_min:.4f} / {l_max:.4f}")
    print(f"  logit mean/std: {l_mean:.4f} / {l_std:.4f}")

    # Red flags
    flags = []
    if not finite:
        flags.append("NON_FINITE_LOGITS")
    if l_std < 1e-4:
        flags.append(f"NEAR_CONSTANT_LOGITS (std={l_std:.2e})")
    if n_pos == 0 or n_neg == 0:
        flags.append(f"SINGLE_CLASS_BATCH (pos={n_pos}, neg={n_neg})")

    # First 10 (label, logit, tag) tuples
    print(f"\n[5/5] First 10 (label, logit, tag)")
    for i in range(min(10, n)):
        t = tags[i] if i < len(tags) else "?"
        print(f"  {i:3d}  label={int(labels[i])}  logit={logits[i]:+.4f}  tag={t}")

    print()
    print("=" * 60)
    if flags:
        print(f"SMOKE TEST FAILED: {flags}")
        print("Do NOT run the full grid until this is debugged.")
        sys.exit(1)
    else:
        print("SMOKE TEST PASSED - pipeline healthy, safe to run full grid.")
        print("=" * 60)
        sys.exit(0)


if __name__ == "__main__":
    main()
