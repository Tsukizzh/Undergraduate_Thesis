"""Quick scan: extract n_atoms per sample from .pt files, save to sizes.pt"""
import torch, time
from pathlib import Path

cache_dir = Path("/root/rivermind-data/EZSpecificity/data/10_Step10_pt训练/ezspec_pt_v1")

for split in ["train", "val"]:
    idx = torch.load(cache_dir / split / "index.pt", map_location="cpu", weights_only=False)
    sample_ids = idx["sample_ids"]
    n = len(sample_ids)
    sizes = torch.zeros(n, dtype=torch.int32)

    t0 = time.time()
    for i in range(n):
        sid = int(sample_ids[i])
        p = cache_dir / split / "samples" / f"{sid // 1000:03d}" / f"sample_{sid:06d}.pt"
        if p.exists():
            s = torch.load(p, map_location="cpu", weights_only=False)
            sizes[i] = s["ligand_pos"].shape[0] + s["protein_pos"].shape[0]
        if (i + 1) % 20000 == 0:
            print(f"  [{split}] {i+1}/{n} ({time.time()-t0:.0f}s)")

    out = cache_dir / split / "sizes.pt"
    torch.save(sizes, out)
    print(f"[{split}] Done: {n} samples, min={sizes.min().item()}, max={sizes.max().item()}, "
          f"mean={sizes.float().mean().item():.0f}, saved to {out} ({time.time()-t0:.0f}s)")
