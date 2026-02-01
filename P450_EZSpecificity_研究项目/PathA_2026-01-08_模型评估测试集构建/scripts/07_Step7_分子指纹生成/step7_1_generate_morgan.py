"""
Step 7.1: Morgan Fingerprint Generation
Generates 1024-bit Morgan fingerprints (radius=2) for all substrates.
Output: morgan_fingerprint.npy of shape (N, 1024), dtype=int8
Index alignment: row i corresponds to Substrate_Index=i from Substrates.csv
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PATHA_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

DEFAULT_INPUT = os.path.join(PATHA_ROOT, "data", "04_Step4_格式修正后数据", "Substrates.csv")
DEFAULT_OUTPUT = os.path.join(PATHA_ROOT, "data", "07_Step7_分子指纹", "morgan_fingerprint.npy")

MORGAN_RADIUS = 2
MORGAN_NBITS = 1024


def generate_morgan_fingerprints(input_csv, output_npy, radius, nbits):
    df = pd.read_csv(input_csv)

    if "Substrate_SMILES" not in df.columns:
        sys.exit("ERROR: Column 'Substrate_SMILES' not found.")
    if "Substrate_Index" not in df.columns:
        sys.exit("ERROR: Column 'Substrate_Index' not found.")

    n = len(df)
    expected_indices = set(range(n))
    actual_indices = set(df["Substrate_Index"].values)
    if actual_indices != expected_indices:
        sys.exit(f"ERROR: Substrate_Index not contiguous 0..{n-1}. "
                 f"Missing: {expected_indices - actual_indices}")

    df = df.sort_values("Substrate_Index").reset_index(drop=True)

    results = np.zeros((n, nbits), dtype=np.int8)
    failed = []

    for i, row in df.iterrows():
        idx = int(row["Substrate_Index"])
        smiles = row["Substrate_SMILES"]

        if pd.isna(smiles) or str(smiles).strip() == "":
            failed.append((idx, smiles, "Empty SMILES"))
            continue

        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            failed.append((idx, smiles, "RDKit parse failed"))
            continue

        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)
        arr = np.zeros((nbits,), dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        results[idx] = arr

    if failed:
        print(f"\nERROR: {len(failed)} SMILES failed fingerprint generation:")
        for idx, smi, reason in failed:
            print(f"  Substrate_Index={idx}, SMILES={smi}, Reason={reason}")
        sys.exit(1)

    os.makedirs(os.path.dirname(output_npy), exist_ok=True)
    np.save(output_npy, results)

    # ---- validation ----
    loaded = np.load(output_npy)
    assert loaded.shape == (n, nbits), f"Shape mismatch: {loaded.shape} != ({n}, {nbits})"
    assert loaded.dtype == np.int8, f"Dtype mismatch: {loaded.dtype}"
    nonzero_counts = np.count_nonzero(loaded, axis=1)
    all_zero_rows = np.sum(nonzero_counts == 0)

    # ---- summary ----
    print(f"\n{'='*60}")
    print(f"Morgan Fingerprint Generation Complete")
    print(f"{'='*60}")
    print(f"Input:         {input_csv}")
    print(f"Output:        {output_npy}")
    print(f"Substrates:    {n}")
    print(f"Fingerprint:   radius={radius}, nBits={nbits}")
    print(f"Shape:         {loaded.shape}")
    print(f"Dtype:         {loaded.dtype}")
    print(f"File size:     {os.path.getsize(output_npy) / 1024:.1f} KB")
    print(f"All-zero rows: {all_zero_rows}")
    print(f"Avg bits set:  {nonzero_counts.mean():.1f}")
    print(f"Min bits set:  {nonzero_counts.min()}")
    print(f"Max bits set:  {nonzero_counts.max()}")
    print(f"{'='*60}")

    # spot-check first 3
    print(f"\nSpot-check (first 3 substrates):")
    for i in range(min(3, n)):
        smi = df.loc[df["Substrate_Index"] == i, "Substrate_SMILES"].values[0]
        bits_on = np.count_nonzero(loaded[i])
        print(f"  [{i}] {smi[:50]:50s}  bits_on={bits_on}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Step 7.1: Morgan Fingerprint Generation")
    parser.add_argument("--input-csv", default=DEFAULT_INPUT, help="Path to Substrates.csv")
    parser.add_argument("--output-npy", default=DEFAULT_OUTPUT, help="Output .npy path")
    parser.add_argument("--radius", type=int, default=MORGAN_RADIUS)
    parser.add_argument("--nbits", type=int, default=MORGAN_NBITS)
    args = parser.parse_args()

    t0 = time.time()
    generate_morgan_fingerprints(args.input_csv, args.output_npy, args.radius, args.nbits)
    print(f"\nElapsed: {time.time() - t0:.1f}s")
