"""
Copy all training data for all_split fold 0 to local project directory.
- brenda: only enzyme_features.lmdb (other files already local)
- 7 small families: all feature files + CSVs + structure

CSVs are also renamed (lowercase -> uppercase column names) for code compatibility.

Usage:
    D:/anaconda3/envs/torch/python.exe copy_all_data.py
"""

import os
import sys
import shutil
import time
import pandas as pd

# ============================================================
# Paths
# ============================================================
GDRIVE = "G:/.shortcut-targets-by-id/173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr/ESIBank"
LOCAL_BASE = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "data", "09_Step9_AllSplit训练"
))

FAMILIES = ["Duf", "Esterase", "Gt_acceptor", "Nitrilase",
            "Phosphatase", "Thiolase", "halogenase"]

FOLD = 0

# Column rename mapping (G: drive uses lowercase, code expects uppercase)
COL_RENAME = {
    "enzyme": "Enzyme Index",
    "reaction": "Substrate Index",
    "label": "Label",
    "structure_index": "Dock Index",
}


# ============================================================
# Helpers
# ============================================================
def copy_file(src, dst, label=""):
    """Copy a single file with progress reporting."""
    if not os.path.exists(src):
        print(f"  [SKIP] {label}: source not found: {src}")
        return False

    os.makedirs(os.path.dirname(dst), exist_ok=True)

    if os.path.exists(dst):
        src_size = os.path.getsize(src)
        dst_size = os.path.getsize(dst)
        if src_size == dst_size:
            print(f"  [SKIP] {label}: already exists ({dst_size/1024**2:.1f}MB)")
            return True
        else:
            print(f"  [REDO] {label}: size mismatch (src={src_size}, dst={dst_size})")

    src_size = os.path.getsize(src)
    size_str = f"{src_size/1024**3:.2f}GB" if src_size > 1024**3 else f"{src_size/1024**2:.1f}MB"
    print(f"  [COPY] {label}: {size_str} ...", end="", flush=True)

    t0 = time.time()
    dst_tmp = dst + ".tmp"
    try:
        with open(src, 'rb') as fsrc, open(dst_tmp, 'wb') as fdst:
            copied = 0
            last_report = t0
            while True:
                chunk = fsrc.read(16 * 1024 * 1024)  # 16MB chunks
                if not chunk:
                    break
                fdst.write(chunk)
                copied += len(chunk)
                now = time.time()
                if now - last_report >= 30:
                    pct = copied / src_size * 100
                    rate = copied / (now - t0) / 1024**2
                    print(f"\n    {pct:.1f}% ({rate:.1f}MB/s)", end="", flush=True)
                    last_report = now

        # Verify size
        tmp_size = os.path.getsize(dst_tmp)
        if tmp_size != src_size:
            print(f" FAILED (size mismatch: {tmp_size} vs {src_size})")
            os.remove(dst_tmp)
            return False

        os.replace(dst_tmp, dst)

        elapsed = time.time() - t0
        rate = src_size / elapsed / 1024**2 if elapsed > 0 else 0
        print(f" OK ({elapsed:.0f}s, {rate:.1f}MB/s)")
        return True

    except Exception as e:
        print(f" ERROR: {e}")
        if os.path.exists(dst_tmp):
            os.remove(dst_tmp)
        return False


def copy_and_rename_csv(src, dst_dir, label=""):
    """Copy CSV and create a renamed version with uppercase columns."""
    fname = os.path.basename(src)
    dst_orig = os.path.join(dst_dir, fname)
    dst_renamed = os.path.join(dst_dir, fname.replace(".csv", "_renamed.csv"))

    if os.path.exists(dst_renamed):
        print(f"  [SKIP] {label} CSV: already exists")
        return True

    if not os.path.exists(src):
        print(f"  [SKIP] {label} CSV: source not found")
        return False

    os.makedirs(dst_dir, exist_ok=True)
    shutil.copy2(src, dst_orig)

    df = pd.read_csv(src)
    df.rename(columns=COL_RENAME, inplace=True)
    df.to_csv(dst_renamed, index=False)
    print(f"  [CSV]  {label}: {len(df)} rows -> {fname} + {os.path.basename(dst_renamed)}")
    return True


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 70)
    print("Copy all_split fold 0 data to local project directory")
    print(f"Source: {GDRIVE}")
    print(f"Destination: {LOCAL_BASE}")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    errors = []

    # ----------------------------------------------------------
    # 1. BRENDA: only enzyme_features.lmdb is missing
    # ----------------------------------------------------------
    print(f"\n[1/8] brenda (enzyme_features.lmdb only)")
    src = f"{GDRIVE}/brenda/enzyme_features.lmdb"
    dst = os.path.join(LOCAL_BASE, "brenda", "enzyme_features.lmdb")
    if not copy_file(src, dst, "brenda/enzyme_features"):
        errors.append("brenda/enzyme_features.lmdb")

    # ----------------------------------------------------------
    # 2-8. Small families: all files
    # ----------------------------------------------------------
    for i, fam in enumerate(FAMILIES, start=2):
        print(f"\n[{i}/8] {fam}")
        src_base = f"{GDRIVE}/small_family/{fam}"
        dst_base = os.path.join(LOCAL_BASE, fam)

        # Feature files (at family root level)
        for feat_file in ["enzyme_features.lmdb", "grover_fingerprint.lmdb",
                          "reaction_features.lmdb", "morgan_fingerprint.npy"]:
            src = f"{src_base}/{feat_file}"
            dst = os.path.join(dst_base, feat_file)
            if not copy_file(src, dst, f"{fam}/{feat_file}"):
                errors.append(f"{fam}/{feat_file}")

        # Structure files (from structure/af2/ -> local structure/)
        for struct_file in ["structure_features.lmdb", "sequence_features.lmdb"]:
            src = f"{src_base}/structure/af2/{struct_file}"
            dst = os.path.join(dst_base, "structure", struct_file)
            if not copy_file(src, dst, f"{fam}/structure/{struct_file}"):
                errors.append(f"{fam}/structure/{struct_file}")

        # CSVs (fold 0 only: train, val, test)
        csv_dir = os.path.join(dst_base, "all_split")
        for split in ["training", "val", "testing"]:
            csv_name = f"{split}_datas_{FOLD}.csv"
            src = f"{src_base}/all_split/{csv_name}"
            copy_and_rename_csv(src, csv_dir, f"{fam}/{split}")

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    print("\n" + "=" * 70)
    print(f"Done: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    if errors:
        print(f"ERRORS ({len(errors)}):")
        for e in errors:
            print(f"  - {e}")
        sys.exit(1)
    else:
        print("All files copied successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
