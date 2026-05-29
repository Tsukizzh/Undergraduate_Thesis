#!/usr/bin/env python3
"""Build an EXP003 trainable P450 subset cache from the EXP001 PT cache.

The script creates a new cache directory. It hard-links existing per-sample
files and shared entity feature files from the base cache, then writes subset
index.pt files for train/val/test. It does not modify the base cache.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import tempfile
import time
from collections import Counter
from pathlib import Path

import torch


def sample_path(cache_dir: Path, split: str, sample_id: int) -> Path:
    return cache_dir / split / "samples" / f"{sample_id // 1000:03d}" / f"sample_{sample_id:06d}.pt"


def hardlink_or_copy(src: Path, dst: Path) -> str:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return "exists"
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        shutil.copy2(src, dst)
        return "copy"


def atomic_torch_save(obj, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=str(path.parent))
    os.close(fd)
    tmp_path = Path(tmp_name)
    try:
        torch.save(obj, tmp_path)
        os.replace(tmp_path, path)
    finally:
        if tmp_path.exists():
            tmp_path.unlink()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build_rows_by_split(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    by_split = {"train": [], "val": [], "test": []}
    for row in rows:
        split = str(row["split"]).strip()
        if split not in by_split:
            raise SystemExit(f"unexpected split: {split!r}")
        if str(row.get("cache_sample_id_present", "")).strip().lower() not in {"true", "1", "yes"}:
            raise SystemExit(f"manifest contains non-trainable sample: {row}")
        if str(row.get("cache_file_exists", "")).strip().lower() not in {"true", "1", "yes"}:
            raise SystemExit(f"manifest contains missing cache file: {row}")
        by_split[split].append(row)
    return by_split


def link_shared_dir(base_cache: Path, out_cache: Path, name: str) -> dict[str, int]:
    src_dir = base_cache / name
    dst_dir = out_cache / name
    if not src_dir.exists():
        raise SystemExit(f"missing base shared dir: {src_dir}")
    counts = Counter()
    for src in sorted(src_dir.rglob("*")):
        if src.is_dir():
            continue
        rel = src.relative_to(src_dir)
        dst = dst_dir / rel
        status = hardlink_or_copy(src, dst)
        counts[status] += 1
    return dict(counts)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-cache", required=True)
    parser.add_argument("--trainable-manifest", required=True)
    parser.add_argument("--out-cache", required=True)
    parser.add_argument("--tag", default=None)
    args = parser.parse_args()

    tag = args.tag or time.strftime("%Y%m%d_%H%M%S")
    base_cache = Path(args.base_cache)
    trainable_manifest = Path(args.trainable_manifest)
    out_cache = Path(args.out_cache)

    if not base_cache.exists():
        raise SystemExit(f"base cache does not exist: {base_cache}")
    if not trainable_manifest.exists():
        raise SystemExit(f"trainable manifest does not exist: {trainable_manifest}")
    if out_cache.exists() and any(out_cache.iterdir()):
        raise SystemExit(f"refusing to write into non-empty out cache: {out_cache}")
    out_cache.mkdir(parents=True, exist_ok=True)

    rows = read_manifest(trainable_manifest)
    by_split = build_rows_by_split(rows)

    base_manifest_path = base_cache / "manifest.pt"
    manifest = torch.load(base_manifest_path, map_location="cpu", weights_only=False)
    manifest = dict(manifest)
    manifest["exp003_full_p450_subset"] = True
    manifest["exp003_source_cache"] = str(base_cache)
    manifest["exp003_trainable_manifest"] = str(trainable_manifest)
    manifest["exp003_tag"] = tag
    manifest["exp003_created_at"] = time.strftime("%Y-%m-%d %H:%M:%S")
    atomic_torch_save(manifest, out_cache / "manifest.pt")

    shared_counts = {
        "enzymes": link_shared_dir(base_cache, out_cache, "enzymes"),
        "substrates": link_shared_dir(base_cache, out_cache, "substrates"),
    }

    split_summaries = {}
    all_status = Counter()
    for split, split_rows in by_split.items():
        base_index = torch.load(base_cache / split / "index.pt", map_location="cpu", weights_only=False)
        sample_ids = [int(row["sample_id"]) for row in split_rows]
        wanted = set(sample_ids)
        positions = [
            i for i, sid in enumerate(base_index["sample_ids"].tolist())
            if int(sid) in wanted
        ]
        found_ids = [int(base_index["sample_ids"][i]) for i in positions]
        if set(found_ids) != wanted:
            missing = sorted(wanted - set(found_ids))
            raise SystemExit(f"{split} index missing sample_ids: {missing[:20]} total={len(missing)}")
        if len(found_ids) != len(sample_ids):
            raise SystemExit(f"{split} duplicate or mismatch: manifest={len(sample_ids)} index={len(found_ids)}")

        subset_index = {}
        pos_tensor = torch.tensor(positions, dtype=torch.long)
        for key, value in base_index.items():
            if torch.is_tensor(value):
                subset_index[key] = value.index_select(0, pos_tensor)
            else:
                subset_index[key] = value
        atomic_torch_save(subset_index, out_cache / split / "index.pt")

        link_counts = Counter()
        missing_files = []
        for sid in found_ids:
            src = sample_path(base_cache, split, sid)
            dst = sample_path(out_cache, split, sid)
            if not src.exists():
                missing_files.append(str(src))
                continue
            link_counts[hardlink_or_copy(src, dst)] += 1
        if missing_files:
            raise SystemExit(f"{split} missing sample files: {missing_files[:20]} total={len(missing_files)}")
        all_status.update(link_counts)
        labels = Counter(str(row.get("label", "")) for row in split_rows)
        split_summaries[split] = {
            "manifest_rows": len(split_rows),
            "index_rows": int(len(subset_index["sample_ids"])),
            "first_sample_id": int(subset_index["sample_ids"][0]) if len(subset_index["sample_ids"]) else None,
            "last_sample_id": int(subset_index["sample_ids"][-1]) if len(subset_index["sample_ids"]) else None,
            "label_counts": dict(labels),
            "sample_file_link_counts": dict(link_counts),
        }

    summary = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "tag": tag,
        "base_cache": str(base_cache),
        "trainable_manifest": str(trainable_manifest),
        "out_cache": str(out_cache),
        "total_rows": len(rows),
        "split_summaries": split_summaries,
        "sample_file_link_counts_total": dict(all_status),
        "shared_file_link_counts": shared_counts,
    }
    summary_path = out_cache / f"exp003_p450_subset_cache_summary_{tag}.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
