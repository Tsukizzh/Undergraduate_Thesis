#!/usr/bin/env python3
"""
Merge P450 enzyme information with substrate and reaction data from ESIBank.

This script:
1. Reads the 389 verified P450 enzymes
2. Maps them to enzyme indices in ESIBank
3. Extracts all enzyme-substrate pairs involving these P450s
4. Creates a comprehensive output with substrate SMILES, EC numbers, and reaction labels

Author: Claude Code
Date: 2026-01-02
"""

import pandas as pd
import os
from collections import defaultdict

# File paths
ESIBANK_PATH = r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\brenda"
P450_LIST_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-02_01-46_P450精确验证\P450酶列表_最终版389个.csv"
OUTPUT_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-02_01-46_P450精确验证"


def main():
    print("=" * 60)
    print("P450 Substrate Information Merger")
    print("=" * 60)

    # 1. Load P450 list
    print("\n[1] Loading verified P450 enzyme list...")
    p450_df = pd.read_csv(P450_LIST_PATH)
    p450_uniprots = set(p450_df["uniprot_id"].str.upper().tolist())
    print(f"    Loaded {len(p450_uniprots)} P450 enzymes")

    # 2. Load ESIBank enzymes.csv
    print("\n[2] Loading ESIBank enzymes.csv...")
    enzymes_path = os.path.join(ESIBANK_PATH, "enzymes.csv")
    enzymes_df = pd.read_csv(enzymes_path)
    print(f"    Loaded {len(enzymes_df)} enzymes")

    # 3. Create UniProt to enzyme index mapping
    print("\n[3] Mapping UniProt IDs to enzyme indices...")
    uniprot_to_idx = {}
    idx_to_uniprot = {}

    for idx, row in enzymes_df.iterrows():
        uniprots = str(row.get("uniprots", "")).upper()
        for uniprot in uniprots.split(","):
            uniprot = uniprot.strip()
            if uniprot:
                uniprot_to_idx[uniprot] = idx
                idx_to_uniprot[idx] = uniprot

    # Find P450 enzyme indices
    p450_indices = set()
    p450_uniprot_to_idx = {}
    for uniprot in p450_uniprots:
        if uniprot in uniprot_to_idx:
            idx = uniprot_to_idx[uniprot]
            p450_indices.add(idx)
            p450_uniprot_to_idx[uniprot] = idx

    print(f"    Found {len(p450_indices)} P450 enzymes in ESIBank")
    print(f"    Missing: {len(p450_uniprots) - len(p450_indices)} P450s not in enzymes.csv")

    # 4. Load reaction_smiles.csv (substrate SMILES)
    print("\n[4] Loading reaction_smiles.csv (substrate SMILES)...")
    reaction_smiles_path = os.path.join(ESIBANK_PATH, "reaction_smiles.csv")
    reaction_df = pd.read_csv(reaction_smiles_path)
    print(f"    Loaded {len(reaction_df)} substrates")

    # Create reaction index to SMILES mapping
    reaction_to_smiles = {}
    for idx, row in reaction_df.iterrows():
        reaction_to_smiles[idx] = row.get("substrates", "")

    # 5. Load data.csv (enzyme-reaction pairs)
    print("\n[5] Loading data.csv (enzyme-reaction pairs)...")
    data_path = os.path.join(ESIBANK_PATH, "data.csv")
    data_df = pd.read_csv(data_path)
    print(f"    Loaded {len(data_df)} enzyme-substrate pairs")

    # 6. Filter for P450 enzymes
    print("\n[6] Filtering for P450 enzyme-substrate pairs...")
    p450_data = data_df[data_df["enzyme"].isin(p450_indices)].copy()
    print(f"    Found {len(p450_data)} pairs involving P450 enzymes")

    # 7. Add UniProt ID, substrate SMILES, and sequence
    print("\n[7] Enriching data with UniProt IDs, SMILES, and sequences...")

    # Add UniProt ID
    p450_data["uniprot_id"] = p450_data["enzyme"].apply(
        lambda x: idx_to_uniprot.get(x, "")
    )

    # Add substrate SMILES
    p450_data["substrate_smiles"] = p450_data["reaction"].apply(
        lambda x: reaction_to_smiles.get(x, "")
    )

    # Add enzyme sequence
    p450_data["sequence"] = p450_data["enzyme"].apply(
        lambda x: enzymes_df.loc[x, "sequences"] if x < len(enzymes_df) else ""
    )

    # 8. Merge with P450 annotation data
    print("\n[8] Merging with P450 annotation data...")
    p450_data = p450_data.merge(
        p450_df[["uniprot_id", "protein_families", "has_heme"]],
        on="uniprot_id",
        how="left"
    )

    # 9. Create summary statistics
    print("\n[9] Generating statistics...")
    unique_p450s = p450_data["uniprot_id"].nunique()
    unique_substrates = p450_data["reaction"].nunique()
    positive_pairs = len(p450_data[p450_data["label"] == 1])
    negative_pairs = len(p450_data[p450_data["label"] == 0])

    print(f"    Unique P450 enzymes: {unique_p450s}")
    print(f"    Unique substrates: {unique_substrates}")
    print(f"    Positive pairs (label=1): {positive_pairs}")
    print(f"    Negative pairs (label=0): {negative_pairs}")

    # 10. Sort and select columns
    output_cols = [
        "uniprot_id",
        "protein_families",
        "ecnumber",
        "substrate_smiles",
        "label",
        "has_heme",
        "enzyme",
        "reaction",
        "difficulty",
        "sequence"
    ]

    # Ensure all columns exist
    for col in output_cols:
        if col not in p450_data.columns:
            p450_data[col] = ""

    p450_data = p450_data[output_cols].sort_values(
        ["uniprot_id", "label"],
        ascending=[True, False]
    )

    # 11. Save full detailed table
    print("\n[10] Saving outputs...")

    # Full table with all pairs
    full_output_path = os.path.join(OUTPUT_PATH, "P450酶底物反应详表_完整版.csv")
    p450_data.to_csv(full_output_path, index=False)
    print(f"    Saved: {full_output_path}")
    print(f"    Total rows: {len(p450_data)}")

    # Summary table (one row per P450, aggregating substrates)
    print("\n[11] Creating summary table...")
    summary_data = []

    for uniprot_id, group in p450_data.groupby("uniprot_id"):
        positive_substrates = group[group["label"] == 1]["substrate_smiles"].tolist()
        negative_substrates = group[group["label"] == 0]["substrate_smiles"].tolist()
        ec_numbers = group["ecnumber"].unique().tolist()

        summary_data.append({
            "uniprot_id": uniprot_id,
            "protein_families": group["protein_families"].iloc[0],
            "has_heme": group["has_heme"].iloc[0],
            "ec_numbers": "; ".join([str(ec) for ec in ec_numbers if pd.notna(ec)]),
            "positive_substrate_count": len(positive_substrates),
            "negative_substrate_count": len(negative_substrates),
            "total_pairs": len(group),
            "positive_substrates_smiles": " | ".join(positive_substrates[:10]),  # Limit to first 10
            "sample_positive_substrates": " | ".join(positive_substrates[:3]),  # First 3 for readability
            "sequence": group["sequence"].iloc[0]
        })

    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values("total_pairs", ascending=False)

    summary_output_path = os.path.join(OUTPUT_PATH, "P450酶底物反应详表_汇总版.csv")
    summary_df.to_csv(summary_output_path, index=False)
    print(f"    Saved: {summary_output_path}")
    print(f"    Total P450 enzymes with data: {len(summary_df)}")

    # 12. Statistics report
    print("\n" + "=" * 60)
    print("FINAL STATISTICS")
    print("=" * 60)
    print(f"Total P450 enzymes in list:         {len(p450_uniprots)}")
    print(f"P450 enzymes found in ESIBank:      {unique_p450s}")
    print(f"P450 enzymes not in ESIBank:        {len(p450_uniprots) - unique_p450s}")
    print(f"Total enzyme-substrate pairs:       {len(p450_data)}")
    print(f"  - Positive pairs (label=1):       {positive_pairs}")
    print(f"  - Negative pairs (label=0):       {negative_pairs}")
    print(f"Unique substrates:                  {unique_substrates}")
    print("=" * 60)

    # Print sample
    print("\n[Sample Data - First 5 rows]")
    print(p450_data[["uniprot_id", "ecnumber", "substrate_smiles", "label"]].head().to_string())

    return p450_data, summary_df


if __name__ == "__main__":
    main()
