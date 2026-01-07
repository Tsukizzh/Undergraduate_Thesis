#!/usr/bin/env python3
"""
Classify P450 enzymes by species/organism using UniProt taxonomy data.

Classification scheme (based on Codex + Gemini discussion):
1. Bacteria - soluble enzymes, structural contrast
2. Archaea - rare P450s
3. Fungi - eukaryotic microorganisms
4. Mammals - drug metabolism (human, mouse, rat, etc.)
5. Other_Animals - insects, fish, etc.
6. Plants_BIA - Ranunculales (poppy, Coptis) - BIA-producing plants
7. Plants_Other - other plants (Arabidopsis, rice, etc.)

Author: Claude Code + Codex + Gemini collaboration
Date: 2026-01-03
"""

import pandas as pd
import requests
import time
import os
import json
from collections import defaultdict

# Paths
P450_LIST_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-02_01-46_P450精确验证\P450酶列表_最终版389个.csv"
OUTPUT_PATH = r"C:\Users\Administrator\Desktop\EZSpecificity_Project\提取P450过程日志\2026-01-03_03-10_P450物种分类"

# BIA-producing plant orders (Ranunculales and related)
BIA_ORDERS = ["Ranunculales"]
BIA_FAMILIES = [
    "Papaveraceae",      # Poppy family (opium poppy)
    "Ranunculaceae",     # Buttercup family (Coptis)
    "Berberidaceae",     # Barberry family
    "Menispermaceae",    # Moonseed family
    "Fumariaceae",       # Fumitory family (sometimes in Papaveraceae)
]


def fetch_taxonomy_batch(uniprot_ids, batch_size=100):
    """
    Fetch taxonomy information from UniProt for a batch of IDs.
    Returns a dict mapping uniprot_id -> taxonomy info.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    fields = "accession,organism_name,organism_id,lineage"

    results = {}

    for i in range(0, len(uniprot_ids), batch_size):
        batch = uniprot_ids[i:i+batch_size]
        query = " OR ".join([f"accession:{uid}" for uid in batch])

        params = {
            "query": query,
            "format": "tsv",
            "fields": fields
        }

        try:
            response = requests.get(base_url, params=params, timeout=60)
            response.raise_for_status()

            lines = response.text.strip().split("\n")
            if len(lines) > 1:  # Has header + data
                header = lines[0].split("\t")
                for line in lines[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 4:
                        entry = {
                            "accession": parts[0],
                            "organism_name": parts[1],
                            "organism_id": parts[2],
                            "lineage": parts[3] if len(parts) > 3 else ""
                        }
                        results[parts[0]] = entry

            print(f"  Fetched {len(batch)} IDs, got {len(results)} total results")
            time.sleep(0.5)  # Rate limiting

        except Exception as e:
            print(f"  Error fetching batch: {e}")
            time.sleep(2)

    return results


def classify_by_lineage(lineage, organism_name):
    """
    Classify organism based on taxonomic lineage.
    Returns (category, subcategory, is_bia_related)
    """
    lineage_lower = lineage.lower() if lineage else ""
    organism_lower = organism_name.lower() if organism_name else ""

    # Check for Bacteria
    if "bacteria" in lineage_lower:
        return "Bacteria", "Bacteria", False

    # Check for Archaea
    if "archaea" in lineage_lower:
        return "Archaea", "Archaea", False

    # Check for Fungi
    if "fungi" in lineage_lower:
        return "Fungi", "Fungi", False

    # Check for Plants (Viridiplantae)
    if "viridiplantae" in lineage_lower or "streptophyta" in lineage_lower:
        # Check if BIA-producing plant
        is_bia = False
        for order in BIA_ORDERS:
            if order.lower() in lineage_lower:
                is_bia = True
                break
        for family in BIA_FAMILIES:
            if family.lower() in lineage_lower:
                is_bia = True
                break

        if is_bia:
            return "Plants_BIA", "Ranunculales/BIA-producing", True
        else:
            return "Plants_Other", "Other plants", False

    # Check for Animals (Metazoa)
    if "metazoa" in lineage_lower:
        # Check for mammals
        if "mammalia" in lineage_lower:
            # Check for humans specifically
            if "homo sapiens" in organism_lower:
                return "Mammals", "Human", False
            else:
                return "Mammals", "Other mammals", False
        # Check for insects
        elif "insecta" in lineage_lower:
            return "Other_Animals", "Insects", False
        # Check for fish
        elif "actinopterygii" in lineage_lower or "fish" in organism_lower:
            return "Other_Animals", "Fish", False
        else:
            return "Other_Animals", "Other animals", False

    # Unknown
    return "Unknown", "Unknown", False


def main():
    print("=" * 60)
    print("P450 Species Classification")
    print("Based on Codex + Gemini collaborative design")
    print("=" * 60)

    # 1. Load P450 list
    print("\n[1] Loading P450 enzyme list...")
    df = pd.read_csv(P450_LIST_PATH)
    uniprot_ids = df["uniprot_id"].tolist()
    print(f"    Loaded {len(uniprot_ids)} P450 enzymes")

    # 2. Fetch taxonomy data from UniProt
    print("\n[2] Fetching taxonomy data from UniProt...")
    taxonomy_data = fetch_taxonomy_batch(uniprot_ids)
    print(f"    Retrieved taxonomy for {len(taxonomy_data)} enzymes")

    # 3. Classify each enzyme
    print("\n[3] Classifying enzymes by species...")
    classifications = []
    missing_taxonomy = []

    for uid in uniprot_ids:
        if uid in taxonomy_data:
            tax = taxonomy_data[uid]
            category, subcategory, is_bia = classify_by_lineage(
                tax["lineage"], tax["organism_name"]
            )
            classifications.append({
                "uniprot_id": uid,
                "organism_name": tax["organism_name"],
                "organism_id": tax["organism_id"],
                "lineage": tax["lineage"],
                "category": category,
                "subcategory": subcategory,
                "is_bia_related": is_bia
            })
        else:
            missing_taxonomy.append(uid)
            classifications.append({
                "uniprot_id": uid,
                "organism_name": "",
                "organism_id": "",
                "lineage": "",
                "category": "Unknown",
                "subcategory": "Missing data",
                "is_bia_related": False
            })

    # 4. Create classification DataFrame
    result_df = pd.DataFrame(classifications)

    # 5. Generate statistics
    print("\n[4] Generating statistics...")
    category_counts = result_df["category"].value_counts()
    subcategory_counts = result_df["subcategory"].value_counts()
    bia_count = result_df["is_bia_related"].sum()

    print("\n=== Category Distribution ===")
    for cat, count in category_counts.items():
        print(f"  {cat}: {count} ({100*count/len(result_df):.1f}%)")

    print(f"\n=== BIA-related P450s ===")
    print(f"  BIA-producing plants: {bia_count}")

    # 6. Save results
    print("\n[5] Saving results...")

    # Full classification table
    full_path = os.path.join(OUTPUT_PATH, "P450_species_classification.csv")
    result_df.to_csv(full_path, index=False)
    print(f"    Saved: {full_path}")

    # Summary by category
    summary_data = []
    for cat in category_counts.index:
        cat_df = result_df[result_df["category"] == cat]
        summary_data.append({
            "category": cat,
            "count": len(cat_df),
            "percentage": f"{100*len(cat_df)/len(result_df):.1f}%",
            "sample_organisms": ", ".join(cat_df["organism_name"].head(3).tolist())
        })
    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(OUTPUT_PATH, "P450_species_summary.csv")
    summary_df.to_csv(summary_path, index=False)
    print(f"    Saved: {summary_path}")

    # Statistics JSON
    stats = {
        "total_enzymes": len(result_df),
        "missing_taxonomy": len(missing_taxonomy),
        "category_counts": category_counts.to_dict(),
        "subcategory_counts": subcategory_counts.to_dict(),
        "bia_related_count": int(bia_count),
        "missing_ids": missing_taxonomy
    }
    stats_path = os.path.join(OUTPUT_PATH, "classification_statistics.json")
    with open(stats_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"    Saved: {stats_path}")

    # Per-category files
    print("\n[6] Creating per-category files...")
    for cat in category_counts.index:
        cat_df = result_df[result_df["category"] == cat]
        cat_path = os.path.join(OUTPUT_PATH, f"P450_{cat}.csv")
        cat_df.to_csv(cat_path, index=False)
        print(f"    {cat}: {len(cat_df)} enzymes -> {cat_path}")

    print("\n" + "=" * 60)
    print("CLASSIFICATION COMPLETE")
    print("=" * 60)

    return result_df, stats


if __name__ == "__main__":
    result_df, stats = main()
