#!/usr/bin/env python3
"""Validate that a Q10 PT cache is readable by the EXP008 dataset code."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from torch_geometric.loader import DataLoader


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--exp008-scripts", required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--batch-size", type=int, default=2)
    args = parser.parse_args()

    sys.path.insert(0, str(Path(args.exp008_scripts).resolve()))
    from pt_dataset import PtCacheDataset

    dataset = PtCacheDataset(
        cache_dir=args.cache_dir,
        split=args.split,
        edge_mode="fixed",
        dist_noise=False,
        cutoff=10.0,
        num_r_gaussian=32,
        max_enzyme_len=1450,
        max_substrate_len=280,
        preload=False,
    )
    first = dataset[0]
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, follow_batch=["ligand_index"])
    batch = next(iter(loader))
    report = {
        "cache_dir": args.cache_dir,
        "split": args.split,
        "dataset_len": len(dataset),
        "first": {
            "embedding": list(first.embedding.shape),
            "grover": list(first.grover.shape),
            "protein_x": list(first.protein_x.shape),
            "ligand_x": list(first.ligand_x.shape),
            "complex_edge_index": list(first.complex_edge_index.shape),
            "complex_edge_attr": list(first.complex_edge_attr.shape),
            "label": float(first.label.item()),
            "tag": str(first.tag),
            "str_tag": str(first.str_tag),
        },
        "batch": {
            "num_graphs": int(batch.num_graphs),
            "embedding": list(batch.embedding.shape),
            "grover": list(batch.grover.shape),
            "protein_x": list(batch.protein_x.shape),
            "ligand_x": list(batch.ligand_x.shape),
            "complex_edge_index": list(batch.complex_edge_index.shape),
            "complex_edge_attr": list(batch.complex_edge_attr.shape),
        },
    }
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
