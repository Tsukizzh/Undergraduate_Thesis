#!/usr/bin/env python3
"""Read-only inspection helper for Q10 PT cache compatibility."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch


def describe_value(value):
    out = {"type": type(value).__name__}
    if hasattr(value, "shape"):
        out["shape"] = list(value.shape)
    if hasattr(value, "dtype"):
        out["dtype"] = str(value.dtype)
    if isinstance(value, (list, tuple)):
        out["len"] = len(value)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--split", default="test")
    args = parser.parse_args()

    cache_dir = Path(args.cache_dir)
    index_path = cache_dir / args.split / "index.pt"
    manifest_path = cache_dir / "manifest.pt"
    substrate_index_path = cache_dir / "substrates" / "substrates_index.pt"
    enzyme_index_path = cache_dir / "enzymes" / "enzymes_index.pt"

    report = {
        "cache_dir": str(cache_dir),
        "split": args.split,
        "index_path": str(index_path),
        "index_exists": index_path.exists(),
    }

    if index_path.exists():
        index = torch.load(index_path, map_location="cpu", weights_only=False)
        report["index_keys"] = {k: describe_value(v) for k, v in index.items()}

    if manifest_path.exists():
        manifest = torch.load(manifest_path, map_location="cpu", weights_only=False)
        report["manifest"] = describe_value(manifest)
        if hasattr(manifest, "keys"):
            report["manifest_keys"] = list(manifest.keys())

    if substrate_index_path.exists():
        sub_index = torch.load(substrate_index_path, map_location="cpu", weights_only=False)
        report["substrates_index"] = describe_value(sub_index)
        if hasattr(sub_index, "keys"):
            report["substrates_index_keys"] = list(sub_index.keys())
        idx = sub_index.get("index") if isinstance(sub_index, dict) else None
        if isinstance(idx, dict):
            report["substrates_index_count"] = len(idx)
            report["substrates_index_preview"] = list(idx.items())[:5]

    if enzyme_index_path.exists():
        enz_index = torch.load(enzyme_index_path, map_location="cpu", weights_only=False)
        report["enzymes_index"] = describe_value(enz_index)
        if hasattr(enz_index, "keys"):
            report["enzymes_index_keys"] = list(enz_index.keys())
        idx = enz_index.get("index") if isinstance(enz_index, dict) else None
        if isinstance(idx, dict):
            report["enzymes_index_count"] = len(idx)
            report["enzymes_index_preview"] = list(idx.items())[:5]

    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
