"""
Convert enzyme_features.lmdb embedding arrays from float32 to float16.

Creates NEW *_fp16.lmdb files alongside originals. Does NOT overwrite.

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe "../毕业设计/P450_EZSpecificity_研究项目/PathB_2026-02-12_P450数据集构建与结构优化/scripts/09_Step9_AllSplit训练/convert_enzyme_fp16.py"

    # Dry-run (scan only, no conversion):
    D:/anaconda3/envs/torch/python.exe "...convert_enzyme_fp16.py" --dry-run
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
import time

import lmdb
import numpy as np

STEP9_DATA = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..",
    "data", "09_Step9_AllSplit训练",
))

FAMILIES = ["brenda", "Duf", "Esterase", "Gt_acceptor",
            "Nitrilase", "Phosphatase", "Thiolase"]

FP16_MAX = np.finfo(np.float16).max  # 65504.0


def scan_lmdb(src_path: str) -> dict:
    """Scan LMDB for embedding stats without writing anything."""
    env = lmdb.open(src_path, readonly=True, lock=False, readahead=True,
                    map_size=128 * 1024**3, subdir=False)
    stats = {"n_entries": 0, "global_max_abs": 0.0,
             "overflow_count": 0, "overflow_keys": [],
             "total_elements": 0, "nan_count": 0}
    with env.begin() as txn:
        cursor = txn.cursor()
        for key, value in cursor:
            data = pickle.loads(value)
            emb = data["embedding"]
            stats["n_entries"] += 1
            stats["total_elements"] += emb.size
            max_abs = float(np.max(np.abs(emb)))
            if max_abs > stats["global_max_abs"]:
                stats["global_max_abs"] = max_abs
            if max_abs > FP16_MAX:
                stats["overflow_count"] += 1
                stats["overflow_keys"].append(key.decode())
            if np.any(np.isnan(emb)):
                stats["nan_count"] += 1
    env.close()
    return stats


def convert_lmdb(src_path: str, dst_path: str) -> dict:
    """Convert embedding fields from float32 to float16, write to new LMDB."""
    src_env = lmdb.open(src_path, readonly=True, lock=False, readahead=True,
                        map_size=128 * 1024**3, subdir=False)
    # Estimate target size: ~60% of source (embedding halved + overhead)
    src_size = os.path.getsize(src_path)
    dst_map_size = max(int(src_size * 0.7), 1024**3)  # at least 1GB

    dst_env = lmdb.open(dst_path, map_size=dst_map_size, subdir=False,
                        readonly=False, meminit=False, writemap=True)
    result = {"n_converted": 0, "n_skipped_overflow": 0,
              "src_bytes": src_size, "dst_bytes": 0}
    t0 = time.time()
    batch_size = 500
    with src_env.begin() as src_txn:
        cursor = src_txn.cursor()
        dst_txn = dst_env.begin(write=True)
        for key, value in cursor:
            data = pickle.loads(value)
            emb = data["embedding"]
            # Safety: skip if overflow would occur
            if np.max(np.abs(emb)) > FP16_MAX:
                # Keep as float32 for this entry
                dst_txn.put(key, value)
                result["n_skipped_overflow"] += 1
            else:
                data["embedding"] = emb.astype(np.float16)
                dst_txn.put(key, pickle.dumps(data, protocol=pickle.HIGHEST_PROTOCOL))
                result["n_converted"] += 1

            if (result["n_converted"] + result["n_skipped_overflow"]) % batch_size == 0:
                dst_txn.commit()
                dst_txn = dst_env.begin(write=True)
                elapsed = time.time() - t0
                total = result["n_converted"] + result["n_skipped_overflow"]
                print(f"  [{total} entries] {elapsed:.0f}s", flush=True)

        dst_txn.commit()

    dst_env.sync()
    dst_env.close()
    src_env.close()
    result["dst_bytes"] = os.path.getsize(dst_path)
    result["time_s"] = time.time() - t0
    return result


def verify_lmdb(src_path: str, dst_path: str, n_samples: int = 20) -> bool:
    """Spot-check converted LMDB: compare a few entries for correctness."""
    src_env = lmdb.open(src_path, readonly=True, lock=False,
                        map_size=128 * 1024**3, subdir=False)
    dst_env = lmdb.open(dst_path, readonly=True, lock=False,
                        map_size=128 * 1024**3, subdir=False)
    ok = True
    with src_env.begin() as src_txn, dst_env.begin() as dst_txn:
        src_cursor = src_txn.cursor()
        checked = 0
        for key, src_value in src_cursor:
            if checked >= n_samples:
                break
            dst_value = dst_txn.get(key)
            if dst_value is None:
                print(f"  FAIL: key {key.decode()} missing in dst")
                ok = False
                continue
            src_data = pickle.loads(src_value)
            dst_data = pickle.loads(dst_value)
            # Check non-embedding fields are identical
            if src_data["active_site"] != dst_data["active_site"]:
                print(f"  FAIL: key {key.decode()} active_site mismatch")
                ok = False
            if not np.array_equal(src_data["sequence"], dst_data["sequence"]):
                print(f"  FAIL: key {key.decode()} sequence mismatch")
                ok = False
            # Check embedding: dst should be float16
            if dst_data["embedding"].dtype != np.float16:
                # Could be float32 if overflow was skipped
                if np.max(np.abs(src_data["embedding"])) <= FP16_MAX:
                    print(f"  FAIL: key {key.decode()} expected float16 but got {dst_data['embedding'].dtype}")
                    ok = False
            else:
                # Check numerical closeness
                src_fp16 = src_data["embedding"].astype(np.float16)
                if not np.array_equal(dst_data["embedding"], src_fp16):
                    print(f"  FAIL: key {key.decode()} embedding values differ")
                    ok = False
            # Check entry count in dst matches src
            if dst_data["embedding"].shape != src_data["embedding"].shape:
                print(f"  FAIL: key {key.decode()} shape mismatch {src_data['embedding'].shape} vs {dst_data['embedding'].shape}")
                ok = False
            checked += 1
        # Check total entry count
        src_count = src_txn.stat()["entries"]
        dst_count = dst_txn.stat()["entries"]
        if src_count != dst_count:
            print(f"  FAIL: entry count mismatch src={src_count} dst={dst_count}")
            ok = False

    src_env.close()
    dst_env.close()
    return ok


def main():
    parser = argparse.ArgumentParser(description="Convert enzyme LMDB float32→float16")
    parser.add_argument("--dry-run", action="store_true",
                        help="Scan only, report stats, do not convert")
    parser.add_argument("--families", nargs="+", default=FAMILIES,
                        help="Families to convert (default: all)")
    parser.add_argument("--data-dir", default=STEP9_DATA,
                        help="Base data directory")
    args = parser.parse_args()

    print("=" * 60)
    print("Enzyme LMDB float32 → float16 Converter")
    print("=" * 60)
    print(f"Data dir: {args.data_dir}")
    print(f"Families: {args.families}")
    print(f"Mode: {'DRY-RUN (scan only)' if args.dry_run else 'CONVERT'}")
    print()

    for family in args.families:
        src_path = os.path.join(args.data_dir, family, "enzyme_features.lmdb")
        dst_path = os.path.join(args.data_dir, family, "enzyme_features_fp16.lmdb")

        if not os.path.exists(src_path):
            print(f"[{family}] SKIP: {src_path} not found")
            continue

        src_size_gb = os.path.getsize(src_path) / 1024**3
        print(f"[{family}] Source: {src_size_gb:.2f} GB")

        # Phase 1: Scan
        print(f"[{family}] Scanning...", flush=True)
        stats = scan_lmdb(src_path)
        print(f"[{family}] Entries: {stats['n_entries']}")
        print(f"[{family}] Max |embedding|: {stats['global_max_abs']:.4f} (fp16 max: {FP16_MAX})")
        print(f"[{family}] Overflow entries: {stats['overflow_count']}")
        print(f"[{family}] NaN entries: {stats['nan_count']}")

        if stats["overflow_count"] > 0:
            print(f"[{family}] WARNING: {stats['overflow_count']} entries exceed fp16 range!")
            print(f"[{family}]   Keys: {stats['overflow_keys'][:10]}")

        if args.dry_run:
            safe = "YES" if stats["overflow_count"] == 0 and stats["nan_count"] == 0 else "NO"
            print(f"[{family}] Safe to convert: {safe}")
            print()
            continue

        # Phase 2: Convert
        if os.path.exists(dst_path):
            print(f"[{family}] SKIP: {dst_path} already exists (delete to reconvert)")
            print()
            continue

        print(f"[{family}] Converting...", flush=True)
        result = convert_lmdb(src_path, dst_path)
        dst_size_gb = result["dst_bytes"] / 1024**3
        ratio = result["dst_bytes"] / result["src_bytes"] * 100
        print(f"[{family}] Done in {result['time_s']:.0f}s")
        print(f"[{family}] Result: {dst_size_gb:.2f} GB ({ratio:.1f}% of original)")
        print(f"[{family}] Converted: {result['n_converted']}, Overflow kept float32: {result['n_skipped_overflow']}")

        # Phase 3: Verify
        print(f"[{family}] Verifying (spot-check 20 entries)...", flush=True)
        if verify_lmdb(src_path, dst_path):
            print(f"[{family}] Verification PASSED")
        else:
            print(f"[{family}] Verification FAILED — check output above")

        print()

    print("=" * 60)
    print("DONE. Next steps:")
    print("  1. Update config enzyme_lmdb_path: enzyme_features.lmdb → enzyme_features_fp16.lmdb")
    print("  2. Run a short training epoch to verify")
    print("  3. If OK, optionally delete original enzyme_features.lmdb to free space")
    print("=" * 60)


if __name__ == "__main__":
    main()
