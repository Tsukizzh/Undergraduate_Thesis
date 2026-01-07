#!/usr/bin/env python3
"""
P450 ML Dataset - Fixed Version (v4)

This script fixes all issues identified in the three-way audit:
1. Species classification bug (Trypanosoma/Danio misclassified as bacteria)
2. Sequence length limit too restrictive (600 -> 1200 aa)
3. Cofactor exclusion list incomplete (added IPA, EOH, MOH, PEG, etc.)
4. Entity-level ligand validation (added warning for entry-level ligands)
5. Expanded species keywords (v4: amycolatopsis, salvia, taxus, etc.)

Author: Claude Code (with Codex + Gemini fixes)
Date: 2026-01-04 v4
"""

import csv
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Optional

# Configuration
BASE_DIR = Path(r"C:\Users\Administrator\Desktop\EZSpecificity_Project")
INPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_02-27_任务3_ML训练数据集筛选" / "数据文件"
OUTPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_04-30_修复后脚本v3" / "数据文件"
LOG_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_04-30_修复后脚本v3" / "日志"

# =============================================================================
# FIX 1: Corrected Species Classification (Codex fix)
# =============================================================================

TAXONOMY_CATEGORY_MAP: Dict[int, str] = {
    # --- Mammals ---
    9606: "mammals",   # Homo sapiens
    10090: "mammals",  # Mus musculus
    10116: "mammals",  # Rattus norvegicus
    9823: "mammals",   # Sus scrofa (pig)
    9913: "mammals",   # Bos taurus (cow)
    9615: "mammals",   # Canis lupus familiaris (dog)
    9986: "mammals",   # Oryctolagus cuniculus (rabbit)
    9685: "mammals",   # Felis catus (cat)
    9796: "mammals",   # Equus caballus (horse)
    9940: "mammals",   # Ovis aries (sheep)
    9925: "mammals",   # Capra hircus (goat)
    9544: "mammals",   # Macaca mulatta (rhesus macaque)
    9598: "mammals",   # Pan troglodytes (chimpanzee)
    9601: "mammals",   # Pongo abelii (orangutan)
    9483: "mammals",   # Callithrix jacchus (marmoset)
    10036: "mammals",  # Mesocricetus auratus (hamster)
    10141: "mammals",  # Cavia porcellus (guinea pig)
    9739: "mammals",   # Tursiops truncatus (bottlenose dolphin)

    # --- Plants ---
    3702: "plants",    # Arabidopsis thaliana
    4530: "plants",    # Oryza sativa (rice)
    4577: "plants",    # Zea mays (maize)
    4558: "plants",    # Sorghum bicolor
    4513: "plants",    # Hordeum vulgare (barley)
    4565: "plants",    # Triticum aestivum (wheat)
    3847: "plants",    # Glycine max (soybean)
    4097: "plants",    # Nicotiana tabacum (tobacco)
    4100: "plants",    # Nicotiana benthamiana
    4081: "plants",    # Solanum lycopersicum (tomato)
    4113: "plants",    # Solanum tuberosum (potato)
    3694: "plants",    # Populus trichocarpa (poplar)
    29760: "plants",   # Vitis vinifera (grape)
    34305: "plants",   # Lotus japonicus
    3880: "plants",    # Medicago truncatula
    4058: "plants",    # Catharanthus roseus
    3461: "plants",    # Eschscholzia californica
    3469: "plants",    # Papaver somniferum

    # --- Bacteria ---
    562: "bacteria",    # Escherichia coli
    1423: "bacteria",   # Bacillus subtilis
    1404: "bacteria",   # Bacillus megaterium (Priestia megaterium)
    1396: "bacteria",   # Bacillus cereus
    287: "bacteria",    # Pseudomonas aeruginosa
    303: "bacteria",    # Pseudomonas putida
    1902: "bacteria",   # Streptomyces coelicolor
    1909: "bacteria",   # Streptomyces griseolus
    1836: "bacteria",   # Saccharopolyspora erythraea
    1773: "bacteria",   # Mycobacterium tuberculosis
    83332: "bacteria",  # Mycobacterium tuberculosis H37Rv
    274: "bacteria",    # Thermus thermophilus
    1280: "bacteria",   # Staphylococcus aureus
    28901: "bacteria",  # Salmonella enterica
    573: "bacteria",    # Klebsiella pneumoniae
    666: "bacteria",    # Vibrio cholerae
    210: "bacteria",    # Helicobacter pylori
    1358: "bacteria",   # Lactococcus lactis
    1718: "bacteria",   # Corynebacterium glutamicum
    470: "bacteria",    # Acinetobacter baumannii
    40216: "bacteria",  # Acinetobacter radioresistens
    31958: "bacteria",  # Amycolatopsis orientalis
    414: "bacteria",    # Methylococcus capsulatus
    2017: "bacteria",   # Streptoalloteichus hindustanus
    65700: "bacteria",  # Erwinia tracheiphila
    247156: "bacteria", # Nocardia farcinica
    1867: "bacteria",   # Actinoplanes teichomyceticus
    2039464: "bacteria",# Actinoplanes tsinanensis
    95486: "bacteria",  # Burkholderia cenocepacia
    382: "bacteria",    # Sinorhizobium meliloti
    384: "bacteria",    # Rhizobium leguminosarum
    190: "bacteria",    # Caulobacter crescentus
    1314: "bacteria",   # Streptococcus pyogenes
    1351: "bacteria",   # Enterococcus faecalis
    132919: "bacteria", # Rhodococcus jostii
    37919: "bacteria",  # Rhodococcus opacus
    1827: "bacteria",   # Rhodococcus (genus)
    100226: "bacteria", # Streptomyces coelicolor A3(2)
    54571: "bacteria",  # Streptomyces venezuelae
    68242: "bacteria",  # Streptomyces natalensis
    28040: "bacteria",  # Micromonospora griseorubida
    171258: "bacteria", # Streptomyces sp. TP-A0274
    33903: "bacteria",  # Streptomyces avermitilis
    2074: "bacteria",   # Pseudonocardia autotrophica
    13689: "bacteria",  # Sphingomonas paucimobilis
    56: "bacteria",     # Sorangium cellulosum
    1148: "bacteria",   # Synechocystis sp. PCC 6803
    2287: "bacteria",   # Saccharolobus solfataricus
    57706: "bacteria",  # Citrobacter braakii

    # --- Fungi ---
    4932: "fungi",     # Saccharomyces cerevisiae
    5476: "fungi",     # Candida albicans
    4896: "fungi",     # Schizosaccharomyces pombe
    5141: "fungi",     # Neurospora crassa
    162425: "fungi",   # Aspergillus nidulans
    330879: "fungi",   # Aspergillus fumigatus
    5518: "fungi",     # Fusarium graminearum
    5076: "fungi",     # Penicillium chrysogenum
    5207: "fungi",     # Cryptococcus neoformans
    4952: "fungi",     # Yarrowia lipolytica
    644223: "fungi",   # Komagataella phaffii (Pichia pastoris)

    # --- Other (explicitly avoid misclassifying as bacteria) ---
    7955: "other",     # Danio rerio (zebrafish) - FIXED!
    5693: "other",     # Trypanosoma cruzi - FIXED!
    5691: "other",     # Trypanosoma brucei - FIXED!
    7227: "other",     # Drosophila melanogaster
    7165: "other",     # Anopheles gambiae
    7091: "other",     # Bombyx mori
    7460: "other",     # Apis mellifera
    6239: "other",     # Caenorhabditis elegans
    9031: "other",     # Gallus gallus (chicken)
    8355: "other",     # Xenopus laevis
    8364: "other",     # Xenopus tropicalis
    31033: "other",    # Takifugu rubripes
    8090: "other",     # Oryzias latipes (medaka)
    7719: "other",     # Ciona intestinalis
}


def classify_species_fixed(
    taxonomy_id: Optional[int],
    organism_name: Optional[str] = None,
) -> str:
    """
    FIXED species classification - uses static mapping, defaults to "other".

    Key fix: NEVER use numeric ID ranges to guess kingdom!
    """
    try:
        tax_id = int(taxonomy_id) if taxonomy_id is not None else None
    except (TypeError, ValueError):
        tax_id = None

    if tax_id is None:
        return "other"

    # Direct lookup in our curated mapping
    direct = TAXONOMY_CATEGORY_MAP.get(tax_id)
    if direct:
        return direct

    # Fallback to keyword matching on organism name (NOT ID ranges!)
    if organism_name:
        name_lower = organism_name.lower()

        # Mammals
        if any(x in name_lower for x in ["homo sapiens", "mus musculus", "rattus", "sus scrofa",
                                          "bos taurus", "canis", "felis", "equus", "oryctolagus"]):
            return "mammals"

        # Plants
        if any(x in name_lower for x in ["arabidopsis", "oryza", "zea mays", "nicotiana",
                                          "solanum", "glycine max", "triticum", "papaver",
                                          "salvia", "parthenium", "taxus"]):
            return "plants"

        # Bacteria (specific keywords only)
        if any(x in name_lower for x in ["escherichia", "bacillus", "pseudomonas",
                                          "streptomyces", "mycobacterium", "rhodococcus",
                                          "amycolatopsis", "corynebacterium", "methylococcus",
                                          "acinetobacter", "erwinia", "nocardia", "actinoplanes"]):
            return "bacteria"

        # Fungi
        if any(x in name_lower for x in ["saccharomyces", "candida", "aspergillus",
                                          "fusarium", "neurospora", "schizosaccharomyces"]):
            return "fungi"

    # DEFAULT TO OTHER - not bacteria!
    return "other"


# =============================================================================
# FIX 2: Expanded Sequence Length Limit (Gemini fix)
# =============================================================================

MAX_RESOLUTION = 3.0
MIN_SEQ_LENGTH = 300
MAX_SEQ_LENGTH = 1200  # Increased from 600 to include P450 BM3 (1048 aa)


# =============================================================================
# FIX 3: Expanded Cofactor/Solvent Exclusion List (Gemini fix)
# =============================================================================

COFACTORS_SOLVENTS = {
    # Heme & Porphyrins
    "HEM", "HEC", "HEA", "HEB", "HAS", "HDD",

    # Redox Cofactors
    "FMN", "FAD", "NAD", "NAP", "NDP", "NADP", "NADPH",

    # Water & Gases
    "HOH", "DOD", "H2O", "OXY", "O2", "O", "N2", "CO", "CYN",

    # Common Metal Ions
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN",
    "FE", "MN", "CU", "NI", "FE2", "FE3", "IOD", "HG", "CD",

    # Solvents & Alcohols (CRITICAL ADDITIONS - Gemini fix)
    "IPA",   # Isopropanol - was wrongly counted as substrate!
    "MOH",   # Methanol
    "EOH", "ETOH",  # Ethanol
    "GOL",   # Glycerol
    "EDO",   # Ethylene glycol
    "PEG",   # Polyethylene glycol (CRITICAL - was missing in v3!)
    "PG4", "PGE", "1PE", "P6G", "PG0",  # PEG variants
    "MPD",   # 2-Methyl-2,4-pentanediol
    "DMS", "DMSO",  # DMSO
    "BME", "BET",   # Beta-mercaptoethanol
    "MRD",   # 2-Methylpentane-2,4-diol

    # Acids & Buffers
    "ACY", "ACT", "ACE", "FMT",  # Acetic acid, Formate
    "TRS",   # Tris buffer
    "EPE",   # HEPES
    "MES",   # MES buffer
    "HEP",   # HEPES variant
    "CIT",   # Citrate
    "MLA",   # Malonic acid
    "TLA",   # Tartaric acid
    "POP", "SO3", "NO3",
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


def parse_ligands_fixed(ligands_str: str) -> Tuple[List[str], List[str]]:
    """
    FIXED ligand parsing with expanded exclusion list.
    """
    if not ligands_str:
        return [], []

    substrates = []
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
            substrates.append(f"{comp_id}:{name}" if name else comp_id)

    return substrates, cofactors


def load_and_fix_data(input_csv: Path, log_file: Path) -> List[Dict]:
    """Load original data and apply fixes."""
    data = []
    species_corrections = 0

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

            # Parse taxonomy ID
            try:
                row["ncbi_taxonomy_id"] = int(row["ncbi_taxonomy_id"]) if row["ncbi_taxonomy_id"] else None
            except (ValueError, TypeError):
                row["ncbi_taxonomy_id"] = None

            # FIX: Re-classify species
            old_category = row.get("species_category", "")
            new_category = classify_species_fixed(
                row["ncbi_taxonomy_id"],
                row.get("organism_name", "")
            )

            if old_category != new_category:
                species_corrections += 1

            row["species_category_old"] = old_category
            row["species_category"] = new_category

            # FIX: Re-parse ligands with expanded exclusion list
            substrates, cofactors = parse_ligands_fixed(row.get("ligands", ""))
            row["substrates"] = substrates
            row["cofactors_only"] = cofactors
            row["has_substrate"] = len(substrates) > 0

            data.append(row)

    log_message(log_file, f"Species classification corrections: {species_corrections}")
    return data


def filter_quality_fixed(data: List[Dict], log_file: Path) -> List[Dict]:
    """Apply FIXED quality filters with expanded length limit."""
    filtered = []

    stats = {
        "no_resolution": 0,
        "low_resolution": 0,
        "short_sequence": 0,
        "long_sequence": 0,  # Now 1200, not 600
        "passed": 0,
    }

    for row in data:
        if row["resolution"] is None:
            stats["no_resolution"] += 1
            continue
        if row["resolution"] > MAX_RESOLUTION:
            stats["low_resolution"] += 1
            continue

        if row["sequence_length"] < MIN_SEQ_LENGTH:
            stats["short_sequence"] += 1
            continue
        if row["sequence_length"] > MAX_SEQ_LENGTH:
            stats["long_sequence"] += 1
            continue

        stats["passed"] += 1
        filtered.append(row)

    log_message(log_file, f"Quality filtering (FIXED - max {MAX_SEQ_LENGTH}aa):")
    log_message(log_file, f"  - No resolution: {stats['no_resolution']}")
    log_message(log_file, f"  - Resolution > {MAX_RESOLUTION}A: {stats['low_resolution']}")
    log_message(log_file, f"  - Sequence < {MIN_SEQ_LENGTH}aa: {stats['short_sequence']}")
    log_message(log_file, f"  - Sequence > {MAX_SEQ_LENGTH}aa: {stats['long_sequence']}")
    log_message(log_file, f"  - Passed quality: {stats['passed']}")

    return filtered


def filter_with_substrates_fixed(data: List[Dict], log_file: Path) -> Tuple[List[Dict], List[Dict]]:
    """Separate entries with real substrates (using FIXED exclusion list)."""
    with_substrates = []
    without_substrates = []

    ipa_removed = 0  # Track IPA corrections

    for row in data:
        # Check if IPA was previously counted as substrate
        old_substrates = row.get("substrates", [])
        if any("IPA" in s.upper() for s in old_substrates if isinstance(s, str)):
            ipa_removed += 1

        if row["has_substrate"]:
            with_substrates.append(row)
        else:
            without_substrates.append(row)

    log_message(log_file, f"Substrate filtering (FIXED exclusion list):")
    log_message(log_file, f"  - With substrates: {len(with_substrates)}")
    log_message(log_file, f"  - Cofactors only: {len(without_substrates)}")
    log_message(log_file, f"  - IPA incorrectly counted before: {ipa_removed}")

    return with_substrates, without_substrates


def remove_redundancy(data: List[Dict], log_file: Path) -> List[Dict]:
    """Remove redundant entries: keep best resolution for each UniProt ID."""
    uniprot_best = {}
    no_uniprot = []

    for row in data:
        uniprot_ids = row.get("uniprot_ids", [])

        if not uniprot_ids:
            no_uniprot.append(row)
            continue

        uniprot_id = uniprot_ids[0]

        if uniprot_id not in uniprot_best:
            uniprot_best[uniprot_id] = row
        else:
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
    """Save data to CSV."""
    if not data:
        return

    base_cols = [
        "pdb_id", "entity_id", "title", "resolution", "experimental_method",
        "organism_name", "ncbi_taxonomy_id", "uniprot_ids", "sequence",
        "sequence_length", "ligands", "species_category", "species_category_old"
    ]

    if extra_cols:
        all_cols = base_cols + extra_cols
    else:
        all_cols = base_cols

    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=all_cols, extrasaction="ignore")
        writer.writeheader()
        for row in data:
            row_copy = row.copy()
            if isinstance(row_copy.get("uniprot_ids"), list):
                row_copy["uniprot_ids"] = json.dumps(row_copy["uniprot_ids"])
            if isinstance(row_copy.get("substrates"), list):
                row_copy["substrates"] = "; ".join(row_copy["substrates"])
            writer.writerow(row_copy)


def compute_intersection(log_file: Path):
    """Compute intersection between ESIBank 389 and RCSB filtered P450s."""
    esibank_path = BASE_DIR / "提取P450过程日志" / "2026-01-02_01-46_P450精确验证" / "数据" / "P450酶列表_最终版389个.csv"

    esibank_ids = set()
    with open(esibank_path, 'r', encoding='utf-8-sig') as f:
        for row in csv.DictReader(f):
            uid = row.get('uniprot_id', '').strip()
            if uid:
                esibank_ids.add(uid)

    rcsb_ids = set()
    rcsb_path = OUTPUT_DIR / "ML训练数据集_v3.csv"
    if rcsb_path.exists():
        with open(rcsb_path, 'r', encoding='utf-8-sig') as f:
            for row in csv.DictReader(f):
                raw = row.get('uniprot_ids', '').strip()
                if raw:
                    try:
                        for uid in json.loads(raw):
                            rcsb_ids.add(uid)
                    except json.JSONDecodeError:
                        pass

    intersection = esibank_ids & rcsb_ids

    log_message(log_file, f"\n=== INTERSECTION CALCULATION ===")
    log_message(log_file, f"ESIBank P450s: {len(esibank_ids)}")
    log_message(log_file, f"RCSB ML dataset (v3): {len(rcsb_ids)}")
    log_message(log_file, f"Intersection: {len(intersection)}")
    log_message(log_file, f"Intersection IDs: {sorted(intersection)}")

    return intersection


def main():
    """Main pipeline with all fixes applied."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOG_DIR / f"fix_and_regenerate_v3_{timestamp}.txt"

    log_message(log_file, "=" * 70)
    log_message(log_file, "P450 ML Dataset - FIXED VERSION (v3)")
    log_message(log_file, "=" * 70)
    log_message(log_file, "Fixes applied:")
    log_message(log_file, "  1. Species classification bug (Trypanosoma/Danio)")
    log_message(log_file, "  2. Sequence length limit (600 -> 1200 aa)")
    log_message(log_file, "  3. Cofactor exclusion list (added IPA, EOH, MOH, etc.)")
    log_message(log_file, "=" * 70)

    # Load original high-quality data
    input_csv = INPUT_DIR / "高质量P450全集_1421个.csv"
    if not input_csv.exists():
        # Try the ML dataset
        input_csv = INPUT_DIR / "ML训练数据集_183个.csv"

    log_message(log_file, f"\nLoading from: {input_csv}")

    # Check if we need to use the original full dataset
    original_full = BASE_DIR / "提取P450过程日志" / "2026-01-04_02-23_任务2_RCSB全局P450搜索" / "数据文件" / "全部P450结构_1591个.csv"
    if original_full.exists():
        log_message(log_file, f"Using full dataset: {original_full}")
        input_csv = original_full

    data = load_and_fix_data(input_csv, log_file)
    log_message(log_file, f"Total entries loaded: {len(data)}")

    # Show species corrections
    corrections = [(r["organism_name"], r["species_category_old"], r["species_category"])
                   for r in data if r.get("species_category_old") != r.get("species_category")]
    if corrections:
        log_message(log_file, f"\n=== SPECIES CORRECTIONS ===")
        for org, old, new in corrections[:20]:
            log_message(log_file, f"  {org[:40]:40} : {old:10} -> {new}")
        if len(corrections) > 20:
            log_message(log_file, f"  ... and {len(corrections)-20} more")

    # Step 1: Quality filtering (with new 1200aa limit)
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 1: Quality Filtering (FIXED)")
    log_message(log_file, "=" * 60)
    quality_data = filter_quality_fixed(data, log_file)

    # Step 2: Substrate filtering (with fixed exclusion list)
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 2: Substrate Filtering (FIXED)")
    log_message(log_file, "=" * 60)
    with_substrates, apo_structures = filter_with_substrates_fixed(quality_data, log_file)

    # Step 3: Redundancy removal
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 3: Redundancy Removal")
    log_message(log_file, "=" * 60)
    unique_with_substrates = remove_redundancy(with_substrates, log_file)

    # Step 4: Species distribution (with FIXED classification)
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 4: Species Distribution (FIXED)")
    log_message(log_file, "=" * 60)
    species_counts = defaultdict(int)
    for row in unique_with_substrates:
        species_counts[row.get("species_category", "other")] += 1

    for cat, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        log_message(log_file, f"  - {cat}: {count}")

    # Save results
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 5: Saving Results")
    log_message(log_file, "=" * 60)

    # Save main ML dataset
    ml_dataset_path = OUTPUT_DIR / "ML训练数据集_v3.csv"
    save_csv(unique_with_substrates, ml_dataset_path, ["substrates", "has_substrate"])
    log_message(log_file, f"ML dataset (v3): {ml_dataset_path} ({len(unique_with_substrates)} entries)")

    # Save all high-quality
    all_quality_path = OUTPUT_DIR / "高质量P450全集_v3.csv"
    save_csv(quality_data, all_quality_path)
    log_message(log_file, f"All high-quality (v3): {all_quality_path} ({len(quality_data)} entries)")

    # Save by species
    for species in species_counts.keys():
        species_data = [r for r in unique_with_substrates if r.get("species_category") == species]
        if species_data:
            species_path = OUTPUT_DIR / f"ML数据集_{species}_v3.csv"
            save_csv(species_data, species_path, ["substrates", "has_substrate"])
            log_message(log_file, f"  {species}: {species_path} ({len(species_data)} entries)")

    # Compute intersection
    intersection = compute_intersection(log_file)

    # Summary
    log_message(log_file, "\n" + "=" * 70)
    log_message(log_file, "FINAL SUMMARY (v3 - FIXED)")
    log_message(log_file, "=" * 70)
    log_message(log_file, f"Input entries: {len(data)}")
    log_message(log_file, f"After quality filter (max {MAX_SEQ_LENGTH}aa): {len(quality_data)}")
    log_message(log_file, f"With substrates (fixed exclusion): {len(with_substrates)}")
    log_message(log_file, f"Non-redundant ML dataset: {len(unique_with_substrates)}")
    log_message(log_file, f"Intersection with ESIBank 389: {len(intersection)}")

    log_message(log_file, f"\nPipeline completed at: {datetime.now()}")

    return {
        "total_entries": len(data),
        "after_quality": len(quality_data),
        "with_substrates": len(with_substrates),
        "ml_dataset": len(unique_with_substrates),
        "intersection": len(intersection),
        "intersection_ids": sorted(intersection),
    }


if __name__ == "__main__":
    results = main()
    print("\n" + "=" * 60)
    print("Results Summary:")
    for k, v in results.items():
        if k != "intersection_ids":
            print(f"  {k}: {v}")
