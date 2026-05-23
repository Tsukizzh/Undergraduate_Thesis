#!/usr/bin/env python3
"""Build a Q2 pt-cache overlay from existing baseline per-sample pt files.

This script does not recompute ESM, GROVER, Morgan, docking graphs, or kNN
features. It only reorganizes existing sample_*.pt files according to a Q2
cluster split and writes fresh train/val/test index.pt files.

Default behavior uses hard links. A hard link gives the new cache its own file
path while reusing the same disk blocks as the baseline cache. Deleting either
path does not remove the other path until all links are gone.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import torch


SPLITS = ("train", "val", "test")
INDEX_KEYS = ("sample_ids", "enzyme_ids", "substrate_ids", "graph_shards", "graph_rows")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create Q2 train/val/test pt-cache overlay from existing baseline samples."
    )
    p.add_argument(
        "--root",
        default="/root/autodl-tmp/EZSpecificity/PathD/P450",
        help="PathD/P450 root on the server.",
    )
    p.add_argument("--threshold", default="id60", help="Q2 split threshold, e.g. id60.")
    p.add_argument("--name", default="id60_main", help="Output cache name under pt_cache/.")
    p.add_argument(
        "--source-cache",
        default=None,
        help="Existing baseline pt cache. Defaults to PathD base_from_PathC best baseline cache.",
    )
    p.add_argument(
        "--split-dir",
        default=None,
        help="Directory containing train_samples_<threshold>.csv etc.",
    )
    p.add_argument(
        "--output-root",
        default=None,
        help="Directory that contains output cache directories. Defaults to Q2 pt_cache/.",
    )
    p.add_argument(
        "--file-mode",
        choices=("hardlink", "copy", "symlink"),
        default="hardlink",
        help="How to place files in the new cache. hardlink is fast and space efficient.",
    )
    p.add_argument("--dry-run", action="store_true", help="Validate and print plan only.")
    p.add_argument("--force", action="store_true", help="Allow replacing an existing output dir.")
    p.add_argument(
        "--check-content",
        action="store_true",
        help="Load a small sample from each split after build and check IDs.",
    )
    return p.parse_args()


def default_paths(args: argparse.Namespace) -> dict[str, Path]:
    root = Path(args.root).resolve()
    source_cache = (
        Path(args.source_cache).resolve()
        if args.source_cache
        else root
        / "data/base_from_PathC/cache_best_baseline/pt_cache_allfix_unified/random"
    )
    q2_base = root / "data/q02_sequence_similarity_split/exp002_actual_used_cache_valid"
    split_dir = Path(args.split_dir).resolve() if args.split_dir else q2_base / "splits" / args.threshold
    output_root = Path(args.output_root).resolve() if args.output_root else q2_base / "pt_cache"
    output_dir = output_root / args.name
    return {
        "root": root,
        "source_cache": source_cache,
        "split_dir": split_dir,
        "output_root": output_root,
        "output_dir": output_dir,
    }


def fail(msg: str) -> None:
    raise RuntimeError(msg)


def ensure_exists(path: Path, kind: str) -> None:
    if kind == "dir" and not path.is_dir():
        fail(f"Missing directory: {path}")
    if kind == "file" and not path.is_file():
        fail(f"Missing file: {path}")


def place_file(src: Path, dst: Path, mode: str) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists() or dst.is_symlink():
        fail(f"Destination already exists: {dst}")
    if mode == "hardlink":
        os.link(src, dst)
    elif mode == "copy":
        shutil.copy2(src, dst)
    elif mode == "symlink":
        rel_src = os.path.relpath(src, start=dst.parent)
        dst.symlink_to(rel_src)
    else:
        fail(f"Unsupported file mode: {mode}")


def iter_files(root: Path) -> Iterable[Path]:
    for p in sorted(root.rglob("*")):
        if p.is_file():
            yield p


def read_split_csvs(split_dir: Path, threshold: str) -> dict[str, pd.DataFrame]:
    result: dict[str, pd.DataFrame] = {}
    required_cols = {
        "split",
        "sample_id_from_cache",
        "enzyme_id_from_cache",
        "substrate_id_from_cache",
        "label",
        "graph_shard",
        "graph_row",
    }
    for split in SPLITS:
        path = split_dir / f"{split}_samples_{threshold}.csv"
        ensure_exists(path, "file")
        df = pd.read_csv(path)
        missing = sorted(required_cols - set(df.columns))
        if missing:
            fail(f"{path} missing required columns: {missing}")
        result[split] = df
    return result


def source_sample_path(source_cache: Path, old_split: str, old_sample_id: int) -> Path:
    return (
        source_cache
        / old_split
        / "samples"
        / f"{old_sample_id // 1000:03d}"
        / f"sample_{old_sample_id:06d}.pt"
    )


def validate_inputs(paths: dict[str, Path], threshold: str) -> dict[str, pd.DataFrame]:
    source_cache = paths["source_cache"]
    split_dir = paths["split_dir"]
    ensure_exists(source_cache, "dir")
    ensure_exists(source_cache / "manifest.pt", "file")
    ensure_exists(source_cache / "enzymes", "dir")
    ensure_exists(source_cache / "substrates", "dir")
    for split in SPLITS:
        ensure_exists(source_cache / split / "index.pt", "file")
        ensure_exists(source_cache / split / "samples", "dir")
    ensure_exists(split_dir, "dir")

    dfs = read_split_csvs(split_dir, threshold)

    all_source_keys: list[tuple[str, int]] = []
    for new_split, df in dfs.items():
        source_key_dups = df.duplicated(["split", "sample_id_from_cache"]).sum()
        if source_key_dups:
            fail(f"{new_split} has duplicated original split+sample_id rows: {source_key_dups}")
        for row in df.itertuples(index=False):
            old_split = str(getattr(row, "split"))
            if old_split not in SPLITS:
                fail(f"Unexpected original split {old_split!r} in {new_split}")
            old_sid = int(getattr(row, "sample_id_from_cache"))
            all_source_keys.append((old_split, old_sid))
            src = source_sample_path(source_cache, old_split, old_sid)
            if not src.is_file():
                fail(f"Missing source sample file for {new_split}: {src}")

    if len(set(all_source_keys)) != len(all_source_keys):
        fail("Q2 split CSVs contain duplicated original samples across train/val/test.")

    return dfs


def copy_shared_files(source_cache: Path, tmp_dir: Path, file_mode: str) -> dict[str, int]:
    counts = {"shared_files": 0}
    place_file(source_cache / "manifest.pt", tmp_dir / "manifest.pt", file_mode)
    counts["shared_files"] += 1
    for dirname in ("enzymes", "substrates"):
        src_dir = source_cache / dirname
        dst_dir = tmp_dir / dirname
        for src in iter_files(src_dir):
            dst = dst_dir / src.relative_to(src_dir)
            place_file(src, dst, file_mode)
            counts["shared_files"] += 1
    return counts


def build_one_split(
    source_cache: Path,
    tmp_dir: Path,
    new_split: str,
    df: pd.DataFrame,
    file_mode: str,
) -> dict[str, object]:
    split_dir = tmp_dir / new_split
    samples_dir = split_dir / "samples"
    split_dir.mkdir(parents=True, exist_ok=True)

    n = len(df)
    sample_ids: list[int] = []
    enzyme_ids: list[int] = []
    substrate_ids: list[int] = []
    graph_shards: list[int] = []
    graph_rows: list[int] = []
    source_map_rows: list[dict[str, object]] = []

    pos = int((df["label"] == 1).sum())
    old_split_counts = {k: int(v) for k, v in df["split"].value_counts().sort_index().items()}

    for new_sid, row in enumerate(df.itertuples(index=False)):
        old_split = str(getattr(row, "split"))
        old_sid = int(getattr(row, "sample_id_from_cache"))
        src = source_sample_path(source_cache, old_split, old_sid)
        dst = samples_dir / f"{new_sid // 1000:03d}" / f"sample_{new_sid:06d}.pt"
        place_file(src, dst, file_mode)

        enzyme_id = int(getattr(row, "enzyme_id_from_cache"))
        substrate_id = int(getattr(row, "substrate_id_from_cache"))
        sample_ids.append(new_sid)
        enzyme_ids.append(enzyme_id)
        substrate_ids.append(substrate_id)
        graph_shards.append(int(getattr(row, "graph_shard")))
        graph_rows.append(int(getattr(row, "graph_row")))
        source_map_rows.append(
            {
                "new_split": new_split,
                "new_sample_id": new_sid,
                "source_split": old_split,
                "source_sample_id": old_sid,
                "enzyme_id": enzyme_id,
                "substrate_id": substrate_id,
                "label": int(getattr(row, "label")),
                "cluster_id": getattr(row, "cluster_id", ""),
            }
        )

    index = {
        "sample_ids": torch.tensor(sample_ids, dtype=torch.int32),
        "enzyme_ids": torch.tensor(enzyme_ids, dtype=torch.int32),
        "substrate_ids": torch.tensor(substrate_ids, dtype=torch.int32),
        "graph_shards": torch.tensor(graph_shards, dtype=torch.int16),
        "graph_rows": torch.tensor(graph_rows, dtype=torch.int32),
    }
    torch.save(index, split_dir / "index.pt")
    pd.DataFrame(source_map_rows).to_csv(split_dir / "source_map.csv", index=False)

    return {
        "samples": n,
        "positives": pos,
        "negatives": n - pos,
        "unique_enzymes": int(df["enzyme_id_from_cache"].nunique()),
        "unique_substrates": int(df["substrate_id_from_cache"].nunique()),
        "source_split_counts": old_split_counts,
    }


def verify_index_lengths(cache_dir: Path) -> dict[str, int]:
    counts: dict[str, int] = {}
    for split in SPLITS:
        idx_path = cache_dir / split / "index.pt"
        idx = torch.load(idx_path, map_location="cpu", weights_only=False)
        if sorted(idx.keys()) != sorted(INDEX_KEYS):
            fail(f"{idx_path} has unexpected keys: {sorted(idx.keys())}")
        lens = {k: len(v) for k, v in idx.items()}
        if len(set(lens.values())) != 1:
            fail(f"{idx_path} array length mismatch: {lens}")
        n = lens["sample_ids"]
        expected_ids = torch.arange(n, dtype=torch.int32)
        if not torch.equal(idx["sample_ids"], expected_ids):
            fail(f"{idx_path} sample_ids are not sequential 0..N-1")
        counts[split] = n
    return counts


def check_content(cache_dir: Path) -> dict[str, list[dict[str, object]]]:
    checks: dict[str, list[dict[str, object]]] = {}
    for split in SPLITS:
        idx = torch.load(cache_dir / split / "index.pt", map_location="cpu", weights_only=False)
        n = len(idx["sample_ids"])
        candidates = sorted({0, 1, 2, max(0, n // 2), max(0, n - 3), max(0, n - 2), max(0, n - 1)})
        split_checks = []
        for i in candidates:
            sid = int(idx["sample_ids"][i])
            sample_path = cache_dir / split / "samples" / f"{sid // 1000:03d}" / f"sample_{sid:06d}.pt"
            sample = torch.load(sample_path, map_location="cpu", weights_only=False)
            got = {
                "row": i,
                "sample_id": sid,
                "enzyme_id_index": int(idx["enzyme_ids"][i]),
                "enzyme_id_sample": int(sample["enzyme_id"]),
                "substrate_id_index": int(idx["substrate_ids"][i]),
                "substrate_id_sample": int(sample["substrate_id"]),
                "label_sample": int(sample["label"]),
            }
            if got["enzyme_id_index"] != got["enzyme_id_sample"]:
                fail(f"{split} sample {sid} enzyme mismatch: {got}")
            if got["substrate_id_index"] != got["substrate_id_sample"]:
                fail(f"{split} sample {sid} substrate mismatch: {got}")
            split_checks.append(got)
        checks[split] = split_checks
    return checks


def build_cache(args: argparse.Namespace) -> None:
    paths = default_paths(args)
    source_cache = paths["source_cache"]
    output_dir = paths["output_dir"]
    output_root = paths["output_root"]
    print(f"Source cache : {source_cache}")
    print(f"Split dir    : {paths['split_dir']}")
    print(f"Output dir   : {output_dir}")
    print(f"Threshold    : {args.threshold}")
    print(f"File mode    : {args.file_mode}")

    dfs = validate_inputs(paths, args.threshold)
    summary = {
        split: {
            "samples": int(len(df)),
            "positives": int((df["label"] == 1).sum()),
            "negatives": int((df["label"] == 0).sum()),
            "unique_enzymes": int(df["enzyme_id_from_cache"].nunique()),
            "unique_substrates": int(df["substrate_id_from_cache"].nunique()),
            "source_split_counts": {
                k: int(v) for k, v in df["split"].value_counts().sort_index().items()
            },
        }
        for split, df in dfs.items()
    }
    total = sum(v["samples"] for v in summary.values())
    print("\nPlan:")
    for split, stats in summary.items():
        print(
            f"  {split:5s}: {stats['samples']:6d} samples, "
            f"pos={stats['positives']:4d}, enzymes={stats['unique_enzymes']:4d}, "
            f"from={stats['source_split_counts']}"
        )
    print(f"  total: {total} samples")

    if args.dry_run:
        print("\nDRY-RUN complete. No files written.")
        return

    if output_dir.exists():
        if not args.force:
            fail(f"Output already exists. Use --force only after reviewing: {output_dir}")
        shutil.rmtree(output_dir)

    output_root.mkdir(parents=True, exist_ok=True)
    tmp_dir = output_root / f".{args.name}.tmp_{time.strftime('%Y%m%d_%H%M%S')}"
    if tmp_dir.exists():
        fail(f"Temporary directory already exists: {tmp_dir}")

    try:
        tmp_dir.mkdir(parents=True)
        shared_counts = copy_shared_files(source_cache, tmp_dir, args.file_mode)
        build_summary = {}
        for split in SPLITS:
            build_summary[split] = build_one_split(
                source_cache=source_cache,
                tmp_dir=tmp_dir,
                new_split=split,
                df=dfs[split],
                file_mode=args.file_mode,
            )
        index_counts = verify_index_lengths(tmp_dir)
        content_checks = check_content(tmp_dir) if args.check_content else {}

        report = {
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
            "threshold": args.threshold,
            "name": args.name,
            "file_mode": args.file_mode,
            "source_cache": str(source_cache),
            "split_dir": str(paths["split_dir"]),
            "output_dir": str(output_dir),
            "shared_counts": shared_counts,
            "splits": build_summary,
            "index_counts": index_counts,
            "content_checks": content_checks,
            "notes": [
                "No ESM/GROVER/Morgan/docking/kNN features were recomputed.",
                "New sample_ids are sequential within each new split to avoid collisions.",
                "source_map.csv preserves the original baseline split and sample_id for every row.",
            ],
        }
        reports_dir = tmp_dir / "reports"
        reports_dir.mkdir(exist_ok=True)
        with open(reports_dir / "build_report.json", "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2, ensure_ascii=False)

        os.replace(tmp_dir, output_dir)
    except Exception:
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
        raise

    print("\nDONE")
    print(f"Created: {output_dir}")
    print(f"Report : {output_dir / 'reports/build_report.json'}")


if __name__ == "__main__":
    build_cache(parse_args())
