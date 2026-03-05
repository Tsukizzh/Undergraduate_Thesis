"""
E6 Prep: Copy small family data files from G: drive to local tmp_lmdb.
Copies data.csv + feature LMDBs for each family.
"""
import os
import shutil

G_SMALL_FAMILY = r"G:\我的云端硬盘\ESIBank\small_family"
LOCAL_DIR = r"D:\EZSpecificity_Project\tmp_lmdb"

FAMILIES = ["Duf", "Esterase", "Gt_acceptor", "Halogenase", "Nitrilase", "Phosphatase", "Thiolase"]

NEEDED_FILES = [
    "data.csv",
    "enzyme_features.lmdb",
    "reaction_features.lmdb",
    "grover_fingerprint.lmdb",
    "morgan_fingerprint.npy",
    "structure_features.lmdb",
]

os.makedirs(LOCAL_DIR, exist_ok=True)

for fam in FAMILIES:
    src_dir = os.path.join(G_SMALL_FAMILY, fam)
    if not os.path.isdir(src_dir):
        print(f"SKIP {fam}: directory not found at {src_dir}")
        continue

    print(f"\n=== {fam} ===")
    files = os.listdir(src_dir)
    print(f"  Available: {files}")

    for fname in NEEDED_FILES:
        src = os.path.join(src_dir, fname)
        dst = os.path.join(LOCAL_DIR, f"{fam}_{fname}")

        if os.path.exists(dst):
            size_mb = os.path.getsize(dst) / 1024 / 1024
            print(f"  EXISTS {fname} -> {dst} ({size_mb:.1f} MB)")
            continue

        if not os.path.exists(src):
            print(f"  MISSING {fname}")
            continue

        size_mb = os.path.getsize(src) / 1024 / 1024
        if size_mb > 500:
            print(f"  TOO LARGE {fname} ({size_mb:.0f} MB) - skipping")
            continue

        print(f"  Copying {fname} ({size_mb:.1f} MB)...", end="", flush=True)
        shutil.copy2(src, dst)
        print(" done")

# Also read data.csv headers for each family
print("\n\n=== Data CSV Column Check ===")
import pandas as pd
for fam in FAMILIES:
    csv_path = os.path.join(LOCAL_DIR, f"{fam}_data.csv")
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path, nrows=3)
        print(f"{fam}: cols={list(df.columns)}, shape_full=", end="")
        df_full = pd.read_csv(csv_path)
        print(f"{df_full.shape}")
        print(f"  n_unique_enzyme={df_full['enzyme'].nunique()}, n_unique_reaction={df_full['reaction'].nunique()}")
    else:
        print(f"{fam}: no data.csv")
