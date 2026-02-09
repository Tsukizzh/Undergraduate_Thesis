#!/usr/bin/env python3
"""
P450 PDB Structure Query and Filtering Pipeline

This script performs three main tasks:
1. Query PDB structures for our 389 verified P450 enzymes
2. Search RCSB PDB for ALL P450s with structures, classify by species
3. Filter high-quality P450s for ML training

Author: Claude Code
Date: 2026-01-04

Technical recommendations from:
- Codex: Use SIFTS/PDBe API for UniProt→PDB mapping, PFAM/InterPro for broad search
- Gemini: Double-lock validation, InterPro domain filtering

Bug fixes applied based on Codex review:
- Fixed PDBe mapping parsing for different response formats
- Changed RCSB search to return polymer_entity instead of entry
- Added retry logic with exponential backoff
- Fixed search_rcsb_by_text to use search_term parameter
- Added proper error logging
"""

import requests
import pandas as pd
import json
import time
import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any, Set
from pathlib import Path

# Configuration
BASE_DIR = Path(r"C:\Users\Administrator\Desktop\EZSpecificity_Project")
OUTPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "数据文件"
LOG_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "日志"

# API Endpoints
PDBE_API = "https://www.ebi.ac.uk/pdbe/api"
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_DATA_API = "https://data.rcsb.org/rest/v1/core"
RCSB_GRAPHQL_API = "https://data.rcsb.org/graphql"

# P450 identifiers for searching
P450_PFAM = "PF00067"  # Cytochrome P450 PFAM domain
P450_INTERPRO = ["IPR001128", "IPR036396", "IPR002401"]  # InterPro domains

# Quality thresholds
MAX_RESOLUTION = 3.0  # Angstroms
MIN_COVERAGE = 0.5    # 50% minimum sequence coverage

# Session with retry logic
SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "EZSpecificity-P450Pipeline/2026-01-04",
    "Accept": "application/json",
})


def request_json(
    method: str,
    url: str,
    *,
    log_file: Optional[Path] = None,
    timeout: int = 30,
    max_retries: int = 5,
    backoff_base: float = 1.5,
    **kwargs: Any,
) -> Optional[Dict[str, Any]]:
    """
    Small retry wrapper for transient HTTP failures (429/5xx/timeouts).
    Returns parsed JSON dict on success, else None.
    """
    for attempt in range(max_retries):
        try:
            resp = SESSION.request(method, url, timeout=timeout, **kwargs)
        except requests.RequestException as e:
            sleep_s = backoff_base ** attempt
            if log_file:
                log_message(log_file, f"Request error ({method} {url}): {e}; retry in {sleep_s:.1f}s")
            time.sleep(sleep_s)
            continue

        if resp.status_code == 429:
            retry_after = resp.headers.get("Retry-After")
            try:
                sleep_s = float(retry_after) if retry_after else (backoff_base ** attempt)
            except ValueError:
                sleep_s = backoff_base ** attempt
            if log_file:
                log_message(log_file, f"429 rate-limited ({method} {url}); retry in {sleep_s:.1f}s")
            time.sleep(sleep_s)
            continue

        if 500 <= resp.status_code <= 599:
            sleep_s = backoff_base ** attempt
            if log_file:
                log_message(log_file, f"{resp.status_code} from {url}; retry in {sleep_s:.1f}s")
            time.sleep(sleep_s)
            continue

        if resp.status_code != 200:
            if log_file:
                log_message(log_file, f"HTTP {resp.status_code} ({method} {url}): {resp.text[:300]}")
            return None

        try:
            return resp.json()
        except ValueError:
            if log_file:
                log_message(log_file, f"Non-JSON response from {url}: {resp.text[:300]}")
            return None

    return None


def setup_logging():
    """Setup logging directory and file."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_file = LOG_DIR / f"pipeline_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    return log_file


def log_message(log_file: Path, message: str):
    """Write message to log file and print to console."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    full_message = f"[{timestamp}] {message}"
    print(full_message)
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(full_message + "\n")


# =============================================================================
# TASK 1: Query PDB structures for our 389 P450s
# =============================================================================

def load_our_p450s() -> pd.DataFrame:
    """Load our verified 389 P450 enzymes."""
    p450_file = BASE_DIR / "提取P450过程日志" / "2026-01-02_01-46_P450精确验证" / "数据" / "P450酶列表_最终版389个.csv"
    df = pd.read_csv(p450_file)
    return df


def query_pdbe_mapping(uniprot_id: str, log_file: Optional[Path] = None) -> Dict[str, List[Dict[str, Any]]]:
    """
    Query PDBe SIFTS mapping for a UniProt ID.
    Returns dict mapping PDB ID to list of chain/mapping info.

    Handles both response formats:
    - {uniprot_id: {"PDB": {pdb_id: {"mappings": [...]}}}}
    - {uniprot_id: {"PDB": {pdb_id: [...]}}}
    """
    uniprot_id = (uniprot_id or "").strip().upper()
    if not uniprot_id:
        return {}

    url = f"{PDBE_API}/mappings/uniprot/{uniprot_id}"

    try:
        data = request_json("GET", url, log_file=log_file, timeout=30)
        if not data:
            return {}

        # Handle case-insensitive key lookup
        entry = data.get(uniprot_id) or data.get(uniprot_id.lower())
        if not entry:
            return {}

        pdb_block = entry.get("PDB") or {}

        normalized: Dict[str, List[Dict[str, Any]]] = {}
        for pdb_id, pdb_info in pdb_block.items():
            # Handle format: {pdb_id: {"mappings": [...]}}
            if isinstance(pdb_info, dict) and isinstance(pdb_info.get("mappings"), list):
                normalized[pdb_id.upper()] = pdb_info["mappings"]
            # Handle format: {pdb_id: [...]}
            elif isinstance(pdb_info, list):
                normalized[pdb_id.upper()] = pdb_info
            # Handle format: {pdb_id: {chain_id: [...]}}
            elif isinstance(pdb_info, dict):
                all_mappings = []
                for chain_id, chain_mappings in pdb_info.items():
                    if isinstance(chain_mappings, list):
                        for m in chain_mappings:
                            if isinstance(m, dict):
                                m["chain_id"] = chain_id
                                all_mappings.append(m)
                if all_mappings:
                    normalized[pdb_id.upper()] = all_mappings

        return normalized
    except Exception as e:
        if log_file:
            log_message(log_file, f"PDBe mapping error for {uniprot_id}: {e}")
        return {}


def query_pdb_for_our_p450s(df: pd.DataFrame, log_file: Path) -> pd.DataFrame:
    """
    Query PDB structures for all our 389 P450 enzymes.
    Uses PDBe SIFTS mapping for authoritative UniProt→PDB mapping.
    """
    log_message(log_file, "=" * 60)
    log_message(log_file, "TASK 1: Querying PDB structures for our 389 P450 enzymes")
    log_message(log_file, "=" * 60)

    results = []
    total = len(df)

    for idx, row in df.iterrows():
        uniprot_id = row["uniprot_id"]

        if idx % 50 == 0:
            pct = (idx / total * 100.0) if total else 100.0
            log_message(log_file, f"Progress: {idx}/{total} ({pct:.1f}%)")

        pdb_mappings = query_pdbe_mapping(uniprot_id, log_file=log_file)

        if pdb_mappings:
            for pdb_id, chain_infos in pdb_mappings.items():
                for chain_info in (chain_infos or []):
                    if isinstance(chain_info, dict):
                        # Handle different start/end formats
                        start_info = chain_info.get("start", {})
                        end_info = chain_info.get("end", {})

                        if isinstance(start_info, dict):
                            start_res = start_info.get("author_residue_number") or start_info.get("residue_number")
                        else:
                            start_res = start_info

                        if isinstance(end_info, dict):
                            end_res = end_info.get("author_residue_number") or end_info.get("residue_number")
                        else:
                            end_res = end_info

                        results.append({
                            "uniprot_id": uniprot_id,
                            "pdb_id": pdb_id,
                            "chain_id": chain_info.get("chain_id", chain_info.get("struct_asym_id", "")),
                            "start": start_res,
                            "end": end_res,
                            "coverage": chain_info.get("coverage"),
                            "has_structure": True
                        })
        else:
            results.append({
                "uniprot_id": uniprot_id,
                "pdb_id": None,
                "chain_id": None,
                "start": None,
                "end": None,
                "coverage": None,
                "has_structure": False
            })

        # Rate limiting
        time.sleep(0.15)

    result_df = pd.DataFrame(results)

    # Statistics
    with_structure = result_df[result_df["has_structure"] == True]["uniprot_id"].nunique()
    without_structure = result_df[result_df["has_structure"] == False]["uniprot_id"].nunique()
    total_pdb_entries = result_df[result_df["pdb_id"].notna()]["pdb_id"].nunique()

    log_message(log_file, f"\nTask 1 Results:")
    log_message(log_file, f"  - P450s with PDB structure: {with_structure}")
    log_message(log_file, f"  - P450s without PDB structure: {without_structure}")
    log_message(log_file, f"  - Total unique PDB entries: {total_pdb_entries}")

    return result_df


# =============================================================================
# TASK 2: Search RCSB for ALL P450s with structures
# =============================================================================

def parse_polymer_entity_identifier(identifier: str) -> Optional[Tuple[str, str]]:
    """
    RCSB polymer_entity identifiers are typically like '1ABC_1'.
    Returns (pdb_id, entity_id) or None if unparseable.
    """
    if not identifier or "_" not in identifier:
        return None
    pdb_id, entity_id = identifier.split("_", 1)
    pdb_id = pdb_id.strip().upper()
    entity_id = entity_id.strip()
    if not pdb_id or not entity_id:
        return None
    return pdb_id, entity_id


def search_rcsb_by_interpro(interpro_id: str, log_file: Optional[Path] = None) -> List[str]:
    """
    Search RCSB PDB for structures containing a specific InterPro domain.
    Returns list of polymer_entity identifiers (e.g., 1ABC_1).
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match",
                "value": interpro_id
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "polymer_entity"
    }

    try:
        data = request_json("POST", RCSB_SEARCH_API, json=query, timeout=60, log_file=log_file)
        if data:
            return [hit["identifier"] for hit in data.get("result_set", [])]
        return []
    except Exception as e:
        if log_file:
            log_message(log_file, f"RCSB InterPro search error for {interpro_id}: {e}")
        return []


def search_rcsb_by_pfam(pfam_id: str, log_file: Optional[Path] = None) -> List[str]:
    """
    Search RCSB PDB for structures containing a specific PFAM domain.
    Returns list of polymer_entity identifiers.
    """
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match",
                "value": pfam_id
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "polymer_entity"
    }

    try:
        data = request_json("POST", RCSB_SEARCH_API, json=query, timeout=60, log_file=log_file)
        if data:
            return [hit["identifier"] for hit in data.get("result_set", [])]
        return []
    except Exception as e:
        if log_file:
            log_message(log_file, f"RCSB PFAM search error for {pfam_id}: {e}")
        return []


def search_rcsb_by_text(search_term: str, log_file: Optional[Path] = None) -> List[str]:
    """
    Fallback text search for P450 structures.
    Returns list of polymer_entity identifiers.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "or",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": search_term  # Fixed: now uses search_term parameter
                    }
                },
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": "CYP"
                    }
                }
            ]
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "polymer_entity"
    }

    try:
        data = request_json("POST", RCSB_SEARCH_API, json=query, timeout=60, log_file=log_file)
        if data:
            return [hit["identifier"] for hit in data.get("result_set", [])]
        return []
    except Exception as e:
        if log_file:
            log_message(log_file, f"RCSB text search error for {search_term}: {e}")
        return []


def get_entity_source_organism(pdb_id: str, log_file: Optional[Path] = None) -> Dict[str, Any]:
    """
    Get source organism information for a PDB entry using GraphQL.
    """
    query = """
    query ($pdb_id: String!) {
      entry(entry_id: $pdb_id) {
        polymer_entities {
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
          rcsb_entity_source_organism {
            ncbi_taxonomy_id
            ncbi_scientific_name
            ncbi_common_names
          }
          rcsb_polymer_entity_container_identifiers {
            entity_id
            uniprot_ids
          }
        }
        rcsb_entry_info {
          resolution_combined
          experimental_method
        }
      }
    }
    """

    try:
        pdb_id = (pdb_id or "").strip().upper()
        data = request_json(
            "POST",
            RCSB_GRAPHQL_API,
            json={"query": query, "variables": {"pdb_id": pdb_id}},
            timeout=30,
            log_file=log_file,
        )
        if data and data.get("errors") and log_file:
            log_message(log_file, f"GraphQL errors for {pdb_id}: {data['errors']}")
        return data or {}
    except Exception as e:
        if log_file:
            log_message(log_file, f"GraphQL error for {pdb_id}: {e}")
        return {}


def classify_species(scientific_name: str) -> str:
    """
    Classify organism into broad categories.
    """
    if not scientific_name:
        return "Unknown"

    name_lower = scientific_name.lower()

    # Human
    if "homo sapiens" in name_lower:
        return "Human"

    # Other mammals
    mammal_keywords = ["mus ", "rattus", "bos ", "sus ", "ovis ", "canis ", "felis ",
                       "oryctolagus", "macaca", "pan ", "equus"]
    if any(kw in name_lower for kw in mammal_keywords):
        return "Mammal (non-human)"

    # Plants
    plant_keywords = ["arabidopsis", "oryza", "zea ", "nicotiana", "solanum",
                      "glycine ", "triticum", "papaver", "eschscholzia", "coptis",
                      "catharanthus", "vitis", "lotus", "medicago"]
    if any(kw in name_lower for kw in plant_keywords):
        return "Plant"

    # Bacteria
    bacteria_keywords = ["escherichia", "bacillus", "pseudomonas", "streptomyces",
                         "mycobacterium", "rhodococcus", "thermus"]
    if any(kw in name_lower for kw in bacteria_keywords):
        return "Bacteria"

    # Fungi
    fungi_keywords = ["saccharomyces", "candida", "aspergillus", "fusarium",
                      "penicillium", "neurospora", "schizosaccharomyces"]
    if any(kw in name_lower for kw in fungi_keywords):
        return "Fungi"

    # Insects
    insect_keywords = ["drosophila", "bombyx", "apis ", "anopheles"]
    if any(kw in name_lower for kw in insect_keywords):
        return "Insect"

    return "Other"


def search_all_p450_structures(log_file: Path) -> pd.DataFrame:
    """
    Search RCSB for ALL P450 structures worldwide.
    Uses PFAM + InterPro + text search for comprehensive coverage.
    Now properly returns polymer_entity level results.
    """
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "TASK 2: Searching RCSB for ALL P450 structures")
    log_message(log_file, "=" * 60)

    all_polymer_entities: Set[str] = set()

    # Search by PFAM
    log_message(log_file, f"Searching by PFAM {P450_PFAM}...")
    pfam_results = search_rcsb_by_pfam(P450_PFAM, log_file)
    log_message(log_file, f"  Found {len(pfam_results)} polymer entities")
    all_polymer_entities.update(pfam_results)

    # Search by InterPro domains
    for interpro in P450_INTERPRO:
        log_message(log_file, f"Searching by InterPro {interpro}...")
        interpro_results = search_rcsb_by_interpro(interpro, log_file)
        log_message(log_file, f"  Found {len(interpro_results)} polymer entities")
        all_polymer_entities.update(interpro_results)

    # Text search as fallback
    log_message(log_file, "Searching by text (cytochrome P450, CYP)...")
    text_results = search_rcsb_by_text("cytochrome P450", log_file)
    log_message(log_file, f"  Found {len(text_results)} polymer entities")
    all_polymer_entities.update(text_results)

    # Parse polymer entity identifiers to group by PDB entry
    parsed = [parse_polymer_entity_identifier(x) for x in all_polymer_entities]
    parsed = [p for p in parsed if p is not None]
    entry_to_entities: Dict[str, Set[str]] = {}
    for pdb_id, entity_id in parsed:
        entry_to_entities.setdefault(pdb_id, set()).add(entity_id)

    log_message(log_file, f"\nTotal unique P450 polymer entities found: {len(all_polymer_entities)}")
    log_message(log_file, f"Total unique PDB entries found: {len(entry_to_entities)}")

    # Fetch metadata for each entry
    log_message(log_file, "\nFetching metadata for all entries...")
    results = []
    total = len(entry_to_entities)

    for idx, (pdb_id, wanted_entity_ids) in enumerate(entry_to_entities.items()):
        if idx % 100 == 0:
            pct = (idx / total * 100.0) if total else 100.0
            log_message(log_file, f"Progress: {idx}/{total} ({pct:.1f}%)")

        entity_data = get_entity_source_organism(pdb_id, log_file)

        if entity_data and "data" in entity_data and entity_data["data"] and entity_data["data"].get("entry"):
            entry = entity_data["data"]["entry"]

            # Get resolution and method
            entry_info = entry.get("rcsb_entry_info", {}) or {}
            resolution = entry_info.get("resolution_combined")
            if resolution and isinstance(resolution, list):
                resolution = resolution[0] if resolution else None
            method = entry_info.get("experimental_method")
            if method and isinstance(method, list):
                method = method[0] if method else None

            # Process each polymer entity (only the ones we searched for)
            for entity in entry.get("polymer_entities", []) or []:
                container = entity.get("rcsb_polymer_entity_container_identifiers", {}) or {}
                entity_id = str(container.get("entity_id") or "").strip()

                # Only process entities that matched our P450 search
                if entity_id not in wanted_entity_ids:
                    continue

                source_organisms = entity.get("rcsb_entity_source_organism", []) or []
                uniprot_ids = container.get("uniprot_ids", []) or []

                for org in source_organisms:
                    if org:
                        scientific_name = org.get("ncbi_scientific_name", "")
                        taxonomy_id = org.get("ncbi_taxonomy_id")
                        species_category = classify_species(scientific_name)

                        results.append({
                            "pdb_id": pdb_id,
                            "entity_id": entity_id,
                            "scientific_name": scientific_name,
                            "taxonomy_id": taxonomy_id,
                            "species_category": species_category,
                            "resolution": resolution,
                            "method": method,
                            "uniprot_ids": ";".join(uniprot_ids) if uniprot_ids else None
                        })
                        break  # Only take first organism per entity

        # Rate limiting - increased to avoid 429 errors
        time.sleep(0.2)

    result_df = pd.DataFrame(results)

    # Remove duplicates (same polymer entity can appear multiple times from different sources)
    if not result_df.empty and "entity_id" in result_df.columns:
        result_df = result_df.drop_duplicates(subset=["pdb_id", "entity_id"])

    # Statistics
    log_message(log_file, f"\nTask 2 Results:")
    log_message(log_file, f"  - Total P450 polymer entities: {len(result_df)}")
    unique_pdb = result_df["pdb_id"].nunique() if not result_df.empty else 0
    log_message(log_file, f"  - Unique PDB entries: {unique_pdb}")

    if not result_df.empty:
        species_counts = result_df["species_category"].value_counts()
        log_message(log_file, f"\n  Species distribution:")
        for species, count in species_counts.items():
            log_message(log_file, f"    - {species}: {count}")

    return result_df


# =============================================================================
# TASK 3: Filter high-quality P450s for ML
# =============================================================================

def load_substrate_data() -> pd.DataFrame:
    """Load substrate data for our P450s."""
    substrate_file = BASE_DIR / "提取P450过程日志" / "2026-01-02_23-00_底物数据整合" / "P450酶底物反应详表_完整版.csv"
    if substrate_file.exists():
        return pd.read_csv(substrate_file)
    return pd.DataFrame()


def filter_high_quality_p450s(
    our_pdb_df: pd.DataFrame,
    all_pdb_df: pd.DataFrame,
    log_file: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Filter high-quality P450s for ML training.

    Criteria (based on Codex recommendations):
    - Resolution ≤ 3.0 Å (must have resolution data)
    - X-ray or cryo-EM method (must have method data)
    - Has UniProt mapping
    - Preferably has substrate data
    """
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "TASK 3: Filtering high-quality P450s for ML")
    log_message(log_file, "=" * 60)

    # Load substrate data
    substrate_df = load_substrate_data()
    substrate_uniprots = set()
    if not substrate_df.empty and "uniprot_id" in substrate_df.columns:
        substrate_uniprots = set(substrate_df["uniprot_id"].dropna().unique())

    log_message(log_file, f"P450s with substrate data: {len(substrate_uniprots)}")

    if all_pdb_df.empty:
        log_message(log_file, "No PDB structures found in Task 2, skipping Task 3")
        return pd.DataFrame(), pd.DataFrame()

    # Filter all_pdb_df for quality
    quality_filtered = all_pdb_df.copy()

    # Resolution filter - require resolution data to be present
    initial_count = len(quality_filtered)
    quality_filtered = quality_filtered[
        quality_filtered["resolution"].notna() &
        (quality_filtered["resolution"] <= MAX_RESOLUTION)
    ]
    log_message(log_file, f"After resolution filter (≤{MAX_RESOLUTION}Å, must have data): {len(quality_filtered)} (removed {initial_count - len(quality_filtered)})")

    # Method filter (X-ray or cryo-EM) - require method data to be present
    valid_methods = ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
    initial_count = len(quality_filtered)
    quality_filtered = quality_filtered[
        quality_filtered["method"].notna() &
        quality_filtered["method"].isin(valid_methods)
    ]
    log_message(log_file, f"After method filter (X-ray/cryo-EM only): {len(quality_filtered)} (removed {initial_count - len(quality_filtered)})")

    # UniProt mapping filter
    initial_count = len(quality_filtered)
    quality_filtered = quality_filtered[quality_filtered["uniprot_ids"].notna()]
    log_message(log_file, f"After UniProt mapping filter: {len(quality_filtered)} (removed {initial_count - len(quality_filtered)})")

    # Check for substrate data
    def has_substrate(uniprot_ids_str):
        if pd.isna(uniprot_ids_str):
            return False
        uniprots = str(uniprot_ids_str).split(";")
        return any(u.strip() in substrate_uniprots for u in uniprots)

    quality_filtered["has_substrate_data"] = quality_filtered["uniprot_ids"].apply(has_substrate)

    # Split into two sets
    with_substrate = quality_filtered[quality_filtered["has_substrate_data"] == True].copy()
    without_substrate = quality_filtered[quality_filtered["has_substrate_data"] == False].copy()

    log_message(log_file, f"\nFinal Results:")
    log_message(log_file, f"  - High-quality with substrate data: {len(with_substrate)}")
    log_message(log_file, f"  - High-quality without substrate data (potential expansion): {len(without_substrate)}")

    return with_substrate, without_substrate


# =============================================================================
# Main Pipeline
# =============================================================================

def main():
    """Run the complete pipeline."""
    # Setup
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    log_file = setup_logging()

    log_message(log_file, "=" * 60)
    log_message(log_file, "P450 PDB Structure Query and Filtering Pipeline")
    log_message(log_file, f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_message(log_file, "=" * 60)

    # Task 1: Query our 389 P450s
    our_p450s_df = load_our_p450s()
    log_message(log_file, f"\nLoaded {len(our_p450s_df)} verified P450 enzymes")

    our_pdb_results = query_pdb_for_our_p450s(our_p450s_df, log_file)
    our_pdb_results.to_csv(OUTPUT_DIR / "task1_our_p450s_pdb_mapping.csv", index=False, encoding="utf-8-sig")

    # Save summary of P450s with and without structures
    with_structure = our_pdb_results[our_pdb_results["has_structure"] == True]
    without_structure = our_pdb_results[our_pdb_results["has_structure"] == False]

    # Get unique UniProt IDs for each category
    with_structure_unique = with_structure.drop_duplicates(subset=["uniprot_id"])
    without_structure_unique = without_structure.drop_duplicates(subset=["uniprot_id"])

    with_structure.to_csv(OUTPUT_DIR / "task1_p450s_with_pdb_structure.csv", index=False, encoding="utf-8-sig")
    without_structure.to_csv(OUTPUT_DIR / "task1_p450s_without_pdb_structure.csv", index=False, encoding="utf-8-sig")
    with_structure_unique.to_csv(OUTPUT_DIR / "task1_p450s_with_pdb_structure_unique.csv", index=False, encoding="utf-8-sig")
    without_structure_unique.to_csv(OUTPUT_DIR / "task1_p450s_without_pdb_structure_unique.csv", index=False, encoding="utf-8-sig")

    # Task 2: Search all P450 structures
    all_p450_structures = search_all_p450_structures(log_file)
    all_p450_structures.to_csv(OUTPUT_DIR / "task2_all_pdb_p450_structures.csv", index=False, encoding="utf-8-sig")

    # Save species classification summary
    if not all_p450_structures.empty:
        species_summary = all_p450_structures.groupby("species_category").agg({
            "pdb_id": "count",
            "resolution": "mean"
        }).rename(columns={"pdb_id": "count", "resolution": "avg_resolution"})
        species_summary.to_csv(OUTPUT_DIR / "task2_species_classification_summary.csv", encoding="utf-8-sig")

        # Save by species category
        for species in all_p450_structures["species_category"].unique():
            species_df = all_p450_structures[all_p450_structures["species_category"] == species]
            safe_name = species.replace(" ", "_").replace("(", "").replace(")", "")
            species_df.to_csv(OUTPUT_DIR / f"task2_p450_{safe_name}.csv", index=False, encoding="utf-8-sig")

    # Task 3: Filter high-quality P450s
    high_quality_with_substrate, high_quality_without_substrate = filter_high_quality_p450s(
        our_pdb_results, all_p450_structures, log_file
    )

    high_quality_with_substrate.to_csv(
        OUTPUT_DIR / "task3_high_quality_with_substrate.csv",
        index=False, encoding="utf-8-sig"
    )
    high_quality_without_substrate.to_csv(
        OUTPUT_DIR / "task3_high_quality_potential_expansion.csv",
        index=False, encoding="utf-8-sig"
    )

    # Final summary
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "PIPELINE COMPLETE")
    log_message(log_file, "=" * 60)
    log_message(log_file, f"\nOutput files saved to: {OUTPUT_DIR}")
    log_message(log_file, f"\nFiles generated:")
    log_message(log_file, "  Task 1:")
    log_message(log_file, "    - task1_our_p450s_pdb_mapping.csv (all mappings)")
    log_message(log_file, "    - task1_p450s_with_pdb_structure.csv")
    log_message(log_file, "    - task1_p450s_without_pdb_structure.csv")
    log_message(log_file, "    - task1_p450s_with_pdb_structure_unique.csv (unique UniProt IDs)")
    log_message(log_file, "    - task1_p450s_without_pdb_structure_unique.csv (unique UniProt IDs)")
    log_message(log_file, "  Task 2:")
    log_message(log_file, "    - task2_all_pdb_p450_structures.csv")
    log_message(log_file, "    - task2_species_classification_summary.csv")
    log_message(log_file, "    - task2_p450_[Species].csv (per species)")
    log_message(log_file, "  Task 3:")
    log_message(log_file, "    - task3_high_quality_with_substrate.csv")
    log_message(log_file, "    - task3_high_quality_potential_expansion.csv")

    log_message(log_file, f"\nCompleted at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Return summary stats
    return {
        "our_p450s_with_structure": with_structure_unique["uniprot_id"].nunique() if not with_structure_unique.empty else 0,
        "our_p450s_without_structure": without_structure_unique["uniprot_id"].nunique() if not without_structure_unique.empty else 0,
        "total_pdb_p450_structures": len(all_p450_structures),
        "high_quality_with_substrate": len(high_quality_with_substrate),
        "high_quality_potential_expansion": len(high_quality_without_substrate)
    }


if __name__ == "__main__":
    stats = main()
    print("\n" + "=" * 60)
    print("Summary Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
