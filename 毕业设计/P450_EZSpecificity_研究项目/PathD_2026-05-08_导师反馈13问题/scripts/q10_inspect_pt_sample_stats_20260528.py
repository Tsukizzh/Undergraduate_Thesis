#!/usr/bin/env python3
"""Read-only PT cache sample-size stats for Q10 compatibility checks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch


def sample_path(cache_dir: Path, split: str, sample_id: int) -> Path:
    return cache_dir / split / "samples" / f"{sample_id // 1000:03d}" / f"sample_{sample_id:06d}.pt"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--limit", type=int, default=64)
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    index = torch.load(cache_dir / args.split / "index.pt", map_location="cpu", weights_only=False)
    sample_ids = [int(x) for x in index["sample_ids"][: args.limit].tolist()]
    rows = []
    for sid in sample_ids:
        p = sample_path(cache_dir, args.split, sid)
        if p.exists():
            s = torch.load(p, map_location="cpu", weights_only=False)
        else:
            shard_id = int(index["graph_shards"][sample_ids.index(sid)])
            row_id = int(index["graph_rows"][sample_ids.index(sid)])
            shard = torch.load(cache_dir / args.split / f"graph_{shard_id:04d}.pt", map_location="cpu", weights_only=False)
            protein_ptr = shard["protein_ptr"]
            ligand_ptr = shard["ligand_ptr"]
            bond_ptr = shard["bond_ptr"]
            knn_ptr = shard["knn_ptr"]
            rows.append({
                "sample_id": sid,
                "n_ligand": int(ligand_ptr[row_id + 1] - ligand_ptr[row_id]),
                "n_protein": int(protein_ptr[row_id + 1] - protein_ptr[row_id]),
                "n_bond": int(bond_ptr[row_id + 1] - bond_ptr[row_id]),
                "n_knn": int(knn_ptr[row_id + 1] - knn_ptr[row_id]),
                "storage": "shard",
            })
            continue
        rows.append({
            "sample_id": sid,
            "n_ligand": int(s["ligand_pos"].shape[0]),
            "n_protein": int(s["protein_pos"].shape[0]),
            "n_bond": int(s["bond_index"].shape[1]),
            "n_knn": int(s["knn_edge_index"].shape[1]),
            "storage": "sample",
        })

    def stats(name: str):
        vals = [r[name] for r in rows]
        return {
            "min": min(vals) if vals else None,
            "max": max(vals) if vals else None,
            "mean": round(sum(vals) / len(vals), 3) if vals else None,
        }

    print(json.dumps({
        "cache_dir": str(cache_dir),
        "split": args.split,
        "checked": len(rows),
        "stats": {
            "n_ligand": stats("n_ligand"),
            "n_protein": stats("n_protein"),
            "n_bond": stats("n_bond"),
            "n_knn": stats("n_knn"),
        },
        "preview": rows[:10],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
