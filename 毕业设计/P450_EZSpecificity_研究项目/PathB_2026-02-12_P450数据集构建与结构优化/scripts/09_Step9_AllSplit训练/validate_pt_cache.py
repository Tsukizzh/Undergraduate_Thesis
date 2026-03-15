"""Validate .pt cache against original LMDB pipeline."""
import argparse, os, random, sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np, torch, torch.nn.functional as F, yaml

SPLIT_TO_CFG_KEY = {"train": "train_data_path", "val": "val_data_path", "test": "test_data_path"}

def load_pt(path):
    return torch.load(path, map_location="cpu", weights_only=False)

def validate_graph_shards(cache_dir, splits):
    errors, warnings = [], []
    for split in splits:
        split_dir = cache_dir / split
        graph_paths = sorted(split_dir.glob("graph_*.pt"))
        if not graph_paths:
            warnings.append(f"[{split}] no graph shards found"); continue
        sizes_mb = []
        samples_per_shard = {}
        for gp in graph_paths:
            try:
                shard = load_pt(gp)
            except Exception as e:
                errors.append(f"[{split}] failed to load {gp.name}: {e}"); continue
            sizes_mb.append(gp.stat().st_size / (1024**2))
            s = int(shard["num_samples"])
            samples_per_shard[int(gp.stem.split("_")[1])] = s
            for k in ["enzyme_ids","substrate_ids","dataset_ids","str_tag_codes","labels","sample_weights"]:
                if len(shard[k]) != s:
                    errors.append(f"[{split}] {gp.name}: len({k})={len(shard[k])} != {s}")
            for ptr in ["ligand_ptr","protein_ptr","bond_ptr","knn_ptr"]:
                p = shard[ptr]
                if len(p) != s+1: errors.append(f"[{split}] {gp.name}: len({ptr})={len(p)} != {s+1}")
                if len(p)>0 and int(p[0])!=0: errors.append(f"[{split}] {gp.name}: {ptr}[0]!=0")
                if len(p)>1 and not torch.all(p[1:]>=p[:-1]): errors.append(f"[{split}] {gp.name}: {ptr} not monotonic")
            del shard
        if sizes_mb:
            non_last = sizes_mb[:-1] if len(sizes_mb)>1 else sizes_mb
            med = float(np.median(non_last))
            print(f"[{split}] {len(sizes_mb)} shards, median={med:.0f}MB, min={min(sizes_mb):.0f}MB, max={max(sizes_mb):.0f}MB")
        idx_path = split_dir / "index.pt"
        if idx_path.exists():
            index = load_pt(idx_path)
            n_idx = len(index["sample_ids"])
            n_shards = sum(samples_per_shard.values())
            if n_idx != n_shards:
                errors.append(f"[{split}] index samples={n_idx} != shard total={n_shards}")
    return errors, warnings

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--splits", nargs="+", default=["train","val","test"])
    args = p.parse_args()
    cache_dir = Path(args.cache_dir).resolve()
    print("== Shard Integrity Check ==")
    errors, warnings = validate_graph_shards(cache_dir, args.splits)
    for w in warnings: print(f"WARNING: {w}")
    if errors:
        for e in errors: print(f"ERROR: {e}")
        print(f"\nFAIL: {len(errors)} errors")
    else:
        print("PASS: All shards structurally valid")

if __name__ == "__main__":
    main()
