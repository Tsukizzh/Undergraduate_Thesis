#!/usr/bin/env python3
"""
Merge species classification data with substrate reaction data.
Creates integrated files for each species category.

Author: Claude Code
Date: 2026-01-03
"""

import pandas as pd
import os
import json

# Paths
SPECIES_CLASSIFICATION_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-03_03-10_P450物种分类\P450物种分类_完整版389个.csv"
SUBSTRATE_DATA_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-02_23-00_底物数据整合\P450酶底物反应详表_汇总版.csv"
OUTPUT_DIR = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-03_03-10_P450物种分类"


def main():
    print("=" * 70)
    print("P450 Species Classification + Substrate Data Integration")
    print("=" * 70)

    # 1. Load species classification data
    print("\n[1] Loading species classification data...")
    species_df = pd.read_csv(SPECIES_CLASSIFICATION_PATH)
    print(f"    Loaded {len(species_df)} P450 enzymes with species info")

    # 2. Load substrate data
    print("\n[2] Loading substrate reaction data...")
    substrate_df = pd.read_csv(SUBSTRATE_DATA_PATH)
    print(f"    Loaded {len(substrate_df)} P450 enzymes with substrate data")

    # 3. Merge the two datasets
    print("\n[3] Merging datasets on uniprot_id...")
    merged_df = pd.merge(
        species_df,
        substrate_df,
        on="uniprot_id",
        how="left"
    )

    # Count how many have substrate data
    has_substrate = merged_df["positive_substrate_count"].notna()
    print(f"    Total P450 enzymes: {len(merged_df)}")
    print(f"    With substrate data: {has_substrate.sum()}")
    print(f"    Without substrate data: {(~has_substrate).sum()}")

    # 4. Save full merged data
    print("\n[4] Saving merged data...")
    full_output_path = os.path.join(OUTPUT_DIR, "P450_物种分类_含底物数据_完整版.csv")
    merged_df.to_csv(full_output_path, index=False, encoding="utf-8-sig")
    print(f"    Saved: {full_output_path}")

    # 5. Generate per-category files
    print("\n[5] Generating per-category files with substrate data...")

    # Category mapping for Chinese filenames
    category_names = {
        "Plants_BIA": "BIA生产植物",
        "Plants_Other": "其他植物",
        "Mammals": "哺乳动物",
        "Other_Animals": "其他动物",
        "Bacteria": "细菌",
        "Fungi": "真菌",
        "Unknown": "未知"
    }

    category_stats = []

    for category in merged_df["category"].unique():
        cat_df = merged_df[merged_df["category"] == category].copy()

        # Calculate stats
        has_substrate_in_cat = cat_df["positive_substrate_count"].notna()
        total_positive = cat_df["positive_substrate_count"].fillna(0).sum()
        total_negative = cat_df["negative_substrate_count"].fillna(0).sum()
        total_pairs = cat_df["total_pairs"].fillna(0).sum()

        cn_name = category_names.get(category, category)

        stats = {
            "category": category,
            "category_cn": cn_name,
            "enzyme_count": len(cat_df),
            "with_substrate_data": int(has_substrate_in_cat.sum()),
            "without_substrate_data": int((~has_substrate_in_cat).sum()),
            "total_positive_pairs": int(total_positive),
            "total_negative_pairs": int(total_negative),
            "total_pairs": int(total_pairs)
        }
        category_stats.append(stats)

        # Save category file
        cat_filename = f"P450_{cn_name}_{len(cat_df)}个_含底物数据.csv"
        cat_path = os.path.join(OUTPUT_DIR, cat_filename)
        cat_df.to_csv(cat_path, index=False, encoding="utf-8-sig")

        print(f"    {cn_name}: {len(cat_df)} enzymes, {int(has_substrate_in_cat.sum())} with substrate data")
        print(f"        -> Positive pairs: {int(total_positive)}, Negative pairs: {int(total_negative)}")

    # 6. Save statistics
    print("\n[6] Saving statistics...")

    stats_df = pd.DataFrame(category_stats)
    stats_df = stats_df.sort_values("enzyme_count", ascending=False)
    stats_path = os.path.join(OUTPUT_DIR, "物种分类_底物统计汇总.csv")
    stats_df.to_csv(stats_path, index=False, encoding="utf-8-sig")
    print(f"    Saved: {stats_path}")

    # JSON stats
    json_stats = {
        "total_enzymes": len(merged_df),
        "with_substrate_data": int(has_substrate.sum()),
        "without_substrate_data": int((~has_substrate).sum()),
        "total_positive_pairs": int(merged_df["positive_substrate_count"].fillna(0).sum()),
        "total_negative_pairs": int(merged_df["negative_substrate_count"].fillna(0).sum()),
        "categories": category_stats
    }
    json_path = os.path.join(OUTPUT_DIR, "物种分类_底物统计.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_stats, f, indent=2, ensure_ascii=False)
    print(f"    Saved: {json_path}")

    # 7. Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nTotal P450 enzymes: {len(merged_df)}")
    print(f"With substrate data: {has_substrate.sum()} ({100*has_substrate.sum()/len(merged_df):.1f}%)")
    print(f"\nTotal enzyme-substrate pairs: {int(merged_df['total_pairs'].fillna(0).sum())}")
    print(f"  - Positive (label=1): {int(merged_df['positive_substrate_count'].fillna(0).sum())}")
    print(f"  - Negative (label=0): {int(merged_df['negative_substrate_count'].fillna(0).sum())}")

    print("\nPer-category breakdown:")
    print("-" * 70)
    print(f"{'Category':<20} {'Enzymes':>10} {'With Data':>12} {'Pos Pairs':>12} {'Neg Pairs':>12}")
    print("-" * 70)
    for stat in sorted(category_stats, key=lambda x: x["enzyme_count"], reverse=True):
        print(f"{stat['category_cn']:<20} {stat['enzyme_count']:>10} {stat['with_substrate_data']:>12} {stat['total_positive_pairs']:>12} {stat['total_negative_pairs']:>12}")
    print("-" * 70)

    print("\n" + "=" * 70)
    print("DONE!")
    print("=" * 70)

    return merged_df, category_stats


if __name__ == "__main__":
    merged_df, stats = main()
