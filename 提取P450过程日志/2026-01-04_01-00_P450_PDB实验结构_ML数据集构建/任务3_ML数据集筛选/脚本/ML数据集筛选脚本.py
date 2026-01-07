#!/usr/bin/env python3
"""
P450 ML Dataset Filter - Task 3

Filter high-quality P450 structures for machine learning training:
1. Resolution <= 3.0 A (high quality X-ray)
2. Has bound ligand (non-cofactor molecule, may be substrate/inhibitor/drug)
3. Sequence length within reasonable P450 range
4. Remove redundant entries (same UniProt, keep best resolution)

NOTE: "bound_ligands" = non-cofactor ligands, NOT confirmed substrates.
      PDB structures may contain inhibitors, drugs, or other ligands.
      True substrate identification requires functional validation.

Author: Claude Code
Date: 2026-01-04 (Updated: terminology fix - substrates → bound_ligands)
"""

import csv
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple

# Configuration
BASE_DIR = Path(r"C:\Users\Administrator\Desktop\EZSpecificity_Project")
INPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "01_全局搜索"
OUTPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "02_ML筛选"
LOG_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "日志"

# Quality thresholds
MAX_RESOLUTION = 3.0
MIN_SEQ_LENGTH = 300
MAX_SEQ_LENGTH = 600

# Common cofactors and solvents to exclude when identifying bound ligands
COFACTORS_SOLVENTS = {
    "HEM", "HEC", "HEA", "HEB", "HAS",  # Heme variants
    "FMN", "FAD", "NAD", "NAP", "NDP",  # Flavins/NAD
    "HOH", "DOD", "H2O",  # Water
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "NI", "CO",  # Ions
    "EDO", "GOL", "PEG", "MPD", "DMS", "ACT", "FMT", "ACE", "BME",  # Solvents/additives
    "TRS", "EPE", "MES", "HEP",  # Buffers
    "PGE", "PG4", "1PE", "P6G",  # PEG variants
    "OXY", "O2", "O",  # Oxygen
    "CYN", "CO",  # Cyanide, CO
    "FE2", "FE3",  # Iron ions
}


def log_message(log_file: Path, message: str):
    """Write message to log file with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"[{timestamp}] {message}\n")
    try:
        print(f"[{timestamp}] {message}")
    except UnicodeEncodeError:
        print(f"[{timestamp}] {message.encode('ascii', 'replace').decode()}")


def parse_ligands(ligands_str: str) -> Tuple[List[str], List[str]]:
    """
    Parse ligands string and separate into bound ligands vs cofactors.

    NOTE: "bound_ligands" here means non-cofactor molecules. They may be:
    - True substrates (enzyme can catalyze)
    - Inhibitors (bind but not catalyzed)
    - Drugs (for P450 drug metabolism studies)
    - Crystallization additives

    Returns: (bound_ligand_list, cofactor_list)
    """
    if not ligands_str:
        return [], []

    bound_ligands = []
    cofactors = []

    for lig in ligands_str.split("; "):
        if ":" in lig:
            comp_id = lig.split(":")[0].strip()
            name = lig.split(":", 1)[1].strip() if ":" in lig else ""
        else:
            comp_id = lig.strip()
            name = ""

        if comp_id.upper() in COFACTORS_SOLVENTS:
            cofactors.append(f"{comp_id}:{name}" if name else comp_id)
        else:
            bound_ligands.append(f"{comp_id}:{name}" if name else comp_id)

    return bound_ligands, cofactors


def load_data(input_csv: Path) -> List[Dict]:
    """Load P450 data from CSV."""
    data = []
    with open(input_csv, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Parse resolution
            try:
                row["resolution"] = float(row["resolution"]) if row["resolution"] else None
            except ValueError:
                row["resolution"] = None

            # Parse sequence length
            try:
                row["sequence_length"] = int(row["sequence_length"]) if row["sequence_length"] else 0
            except ValueError:
                row["sequence_length"] = 0

            # Parse UniProt IDs
            try:
                row["uniprot_ids"] = json.loads(row["uniprot_ids"]) if row["uniprot_ids"] else []
            except json.JSONDecodeError:
                row["uniprot_ids"] = []

            data.append(row)

    return data


def filter_quality(data: List[Dict], log_file: Path) -> List[Dict]:
    """Apply quality filters."""
    filtered = []

    stats = {
        "no_resolution": 0,
        "low_resolution": 0,
        "short_sequence": 0,
        "long_sequence": 0,
        "passed": 0,
    }

    for row in data:
        # Resolution filter
        if row["resolution"] is None:
            stats["no_resolution"] += 1
            continue
        if row["resolution"] > MAX_RESOLUTION:
            stats["low_resolution"] += 1
            continue

        # Sequence length filter
        if row["sequence_length"] < MIN_SEQ_LENGTH:
            stats["short_sequence"] += 1
            continue
        if row["sequence_length"] > MAX_SEQ_LENGTH:
            stats["long_sequence"] += 1
            continue

        stats["passed"] += 1
        filtered.append(row)

    log_message(log_file, f"Quality filtering:")
    log_message(log_file, f"  - No resolution data: {stats['no_resolution']}")
    log_message(log_file, f"  - Resolution > {MAX_RESOLUTION}A: {stats['low_resolution']}")
    log_message(log_file, f"  - Sequence < {MIN_SEQ_LENGTH}aa: {stats['short_sequence']}")
    log_message(log_file, f"  - Sequence > {MAX_SEQ_LENGTH}aa: {stats['long_sequence']}")
    log_message(log_file, f"  - Passed quality: {stats['passed']}")

    return filtered


def filter_with_ligands(data: List[Dict], log_file: Path) -> Tuple[List[Dict], List[Dict]]:
    """
    Separate entries with bound ligands from those with only cofactors.

    NOTE: "bound_ligands" = non-cofactor molecules (may include inhibitors, drugs, etc.)
          This does NOT confirm they are true enzymatic substrates.

    Returns: (with_ligands, without_ligands)
    """
    with_ligands = []
    without_ligands = []

    for row in data:
        bound_ligands, cofactors = parse_ligands(row.get("ligands", ""))
        row["bound_ligands"] = bound_ligands
        row["cofactors_only"] = cofactors
        row["has_ligand"] = len(bound_ligands) > 0

        if bound_ligands:
            with_ligands.append(row)
        else:
            without_ligands.append(row)

    log_message(log_file, f"Ligand filtering:")
    log_message(log_file, f"  - With bound ligands: {len(with_ligands)}")
    log_message(log_file, f"  - Cofactors only: {len(without_ligands)}")

    return with_ligands, without_ligands


def remove_redundancy(data: List[Dict], log_file: Path) -> List[Dict]:
    """
    Remove redundant entries: keep best resolution for each UniProt ID.
    Entries without UniProt ID are kept separately.
    """
    uniprot_best = {}  # uniprot_id -> best entry
    no_uniprot = []

    for row in data:
        uniprot_ids = row.get("uniprot_ids", [])

        if not uniprot_ids:
            no_uniprot.append(row)
            continue

        # Use first UniProt ID as key
        uniprot_id = uniprot_ids[0]

        if uniprot_id not in uniprot_best:
            uniprot_best[uniprot_id] = row
        else:
            # Keep better resolution
            current_res = uniprot_best[uniprot_id].get("resolution", 99)
            new_res = row.get("resolution", 99)
            if new_res and (current_res is None or new_res < current_res):
                uniprot_best[uniprot_id] = row

    unique = list(uniprot_best.values()) + no_uniprot

    log_message(log_file, f"Redundancy removal:")
    log_message(log_file, f"  - Unique UniProt IDs: {len(uniprot_best)}")
    log_message(log_file, f"  - Entries without UniProt: {len(no_uniprot)}")
    log_message(log_file, f"  - Total non-redundant: {len(unique)}")

    return unique


def save_csv(data: List[Dict], output_path: Path, extra_cols: List[str] = None):
    """Save data to CSV with optional extra columns."""
    if not data:
        return

    base_cols = [
        "pdb_id", "entity_id", "title", "resolution", "experimental_method",
        "organism_name", "ncbi_taxonomy_id", "uniprot_ids", "sequence",
        "sequence_length", "ligands", "species_category"
    ]

    if extra_cols:
        all_cols = base_cols + extra_cols
    else:
        all_cols = base_cols

    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=all_cols, extrasaction="ignore")
        writer.writeheader()
        for row in data:
            # Convert lists back to strings for CSV
            row_copy = row.copy()
            if isinstance(row_copy.get("uniprot_ids"), list):
                row_copy["uniprot_ids"] = json.dumps(row_copy["uniprot_ids"])
            if isinstance(row_copy.get("bound_ligands"), list):
                row_copy["bound_ligands"] = "; ".join(row_copy["bound_ligands"])
            writer.writerow(row_copy)


def main():
    """Main filtering pipeline."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOG_DIR / f"ml_filter_log_{timestamp}.txt"

    log_message(log_file, "=" * 60)
    log_message(log_file, "P450 ML Dataset Filtering - Task 3")
    log_message(log_file, f"Started at: {datetime.now()}")
    log_message(log_file, "=" * 60)

    # Load data
    input_csv = INPUT_DIR / "all_p450_structures.csv"
    log_message(log_file, f"\nLoading data from: {input_csv}")
    data = load_data(input_csv)
    log_message(log_file, f"Total entries loaded: {len(data)}")

    # Step 1: Quality filtering
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 1: Quality Filtering")
    log_message(log_file, "=" * 60)
    quality_data = filter_quality(data, log_file)

    # Step 2: Ligand filtering (non-cofactor molecules)
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 2: Ligand Filtering (excluding cofactors/solvents)")
    log_message(log_file, "=" * 60)
    with_ligands, apo_structures = filter_with_ligands(quality_data, log_file)

    # Step 3: Remove redundancy (for ligand-bound)
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 3: Redundancy Removal (ligand-bound)")
    log_message(log_file, "=" * 60)
    unique_with_ligands = remove_redundancy(with_ligands, log_file)

    # Step 4: Species distribution
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 4: Species Distribution (ML dataset)")
    log_message(log_file, "=" * 60)
    species_counts = defaultdict(int)
    for row in unique_with_ligands:
        species_counts[row.get("species_category", "unknown")] += 1

    for cat, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        log_message(log_file, f"  - {cat}: {count}")

    # Save results
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 5: Saving Results")
    log_message(log_file, "=" * 60)

    # 1. Full ML dataset (with ligands, non-redundant)
    ml_dataset_path = OUTPUT_DIR / "p450_ml_dataset.csv"
    save_csv(unique_with_ligands, ml_dataset_path, ["bound_ligands", "has_ligand"])
    log_message(log_file, f"ML dataset: {ml_dataset_path} ({len(unique_with_ligands)} entries)")

    # 2. All high-quality structures (including apo)
    all_quality_path = OUTPUT_DIR / "p450_high_quality_all.csv"
    save_csv(quality_data, all_quality_path)
    log_message(log_file, f"All high-quality: {all_quality_path} ({len(quality_data)} entries)")

    # 3. By species (for ML dataset)
    for species in species_counts.keys():
        species_data = [r for r in unique_with_ligands if r.get("species_category") == species]
        if species_data:
            species_path = OUTPUT_DIR / f"p450_ml_{species}.csv"
            save_csv(species_data, species_path, ["bound_ligands", "has_ligand"])
            log_message(log_file, f"  ML {species}: {species_path} ({len(species_data)} entries)")

    # Summary statistics
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "FINAL SUMMARY")
    log_message(log_file, "=" * 60)
    log_message(log_file, f"Input entries: {len(data)}")
    log_message(log_file, f"After quality filter: {len(quality_data)}")
    log_message(log_file, f"With bound ligands: {len(with_ligands)}")
    log_message(log_file, f"Non-redundant ML dataset: {len(unique_with_ligands)}")

    # Resolution statistics for ML dataset
    resolutions = [r["resolution"] for r in unique_with_ligands if r["resolution"]]
    if resolutions:
        avg_res = sum(resolutions) / len(resolutions)
        min_res = min(resolutions)
        max_res = max(resolutions)
        log_message(log_file, f"\nML dataset resolution stats:")
        log_message(log_file, f"  - Average: {avg_res:.2f} A")
        log_message(log_file, f"  - Best: {min_res:.2f} A")
        log_message(log_file, f"  - Worst: {max_res:.2f} A")

    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, f"Pipeline completed at: {datetime.now()}")
    log_message(log_file, "=" * 60)

    return unique_with_ligands


if __name__ == "__main__":
    main()
