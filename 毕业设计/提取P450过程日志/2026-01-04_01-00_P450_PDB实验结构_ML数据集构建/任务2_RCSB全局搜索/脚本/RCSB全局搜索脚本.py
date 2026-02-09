#!/usr/bin/env python3
"""
P450 Global Search Pipeline v2

Based on Codex + Gemini recommendations:
1. Use PFAM PF00067 + InterPro domains only (skip text search to avoid false positives)
2. Use RCSB GraphQL API for efficient batch metadata retrieval
3. Use SQLite for incremental saving (resume-safe)
4. Classify by species using NCBI taxonomy ID

Author: Claude Code (with Codex + Gemini collaboration)
Date: 2026-01-04
"""

import requests
import sqlite3
import json
import time
import csv
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Set, Any

# Configuration
BASE_DIR = Path(r"C:\Users\Administrator\Desktop\EZSpecificity_Project")
OUTPUT_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "01_全局搜索"
LOG_DIR = BASE_DIR / "提取P450过程日志" / "2026-01-04_PDB结构查询与筛选" / "日志"

# API Endpoints
RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_API = "https://data.rcsb.org/graphql"

# P450 identifiers (high precision - no text search!)
P450_DOMAINS = {
    "PFAM": ["PF00067"],
    "InterPro": ["IPR001128", "IPR036396", "IPR002401"]
}

# Species classification based on NCBI taxonomy
SPECIES_CATEGORIES = {
    "mammals": {9606, 10090, 10116, 9823, 9913, 9986, 9615},  # Human, Mouse, Rat, Pig, Cow, Rabbit, Dog
    "plants": set(),  # Will check lineage for Viridiplantae (33090)
    "bacteria": set(),  # Will check lineage for Bacteria (2)
    "fungi": set(),  # Will check lineage for Fungi (4751)
}

# Quality thresholds
MAX_RESOLUTION = 3.0  # Angstroms

# Session setup
SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "EZSpecificity-P450Pipeline-v2/2026-01-04",
    "Accept": "application/json",
    "Content-Type": "application/json",
})


def log_message(log_file: Path, message: str):
    """Write message to log file with timestamp."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(log_file, "a", encoding="utf-8") as f:
        f.write(f"[{timestamp}] {message}\n")
    # Use ASCII-safe print for Windows console
    try:
        print(f"[{timestamp}] {message}")
    except UnicodeEncodeError:
        print(f"[{timestamp}] {message.encode('ascii', 'replace').decode()}")


def init_database(db_path: Path) -> sqlite3.Connection:
    """Initialize SQLite database for incremental saving."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create tables
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS polymer_entities (
            pdb_id TEXT,
            entity_id TEXT,
            source TEXT,
            fetched INTEGER DEFAULT 0,
            PRIMARY KEY (pdb_id, entity_id)
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS p450_metadata (
            pdb_id TEXT,
            entity_id TEXT,
            title TEXT,
            resolution REAL,
            experimental_method TEXT,
            organism_name TEXT,
            ncbi_taxonomy_id INTEGER,
            uniprot_ids TEXT,
            sequence TEXT,
            sequence_length INTEGER,
            ligands TEXT,
            species_category TEXT,
            PRIMARY KEY (pdb_id, entity_id)
        )
    """)

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS progress (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)

    conn.commit()
    return conn


def search_by_annotation(annotation_type: str, annotation_id: str, log_file: Path) -> List[str]:
    """Search RCSB for polymer entities by annotation (PFAM or InterPro)."""

    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_annotation.annotation_id",
                "operator": "exact_match",
                "value": annotation_id
            }
        },
        "request_options": {
            "return_all_hits": True
        },
        "return_type": "polymer_entity"
    }

    try:
        resp = SESSION.post(RCSB_SEARCH_API, json=query, timeout=60)
        if resp.status_code == 200:
            data = resp.json()
            results = data.get("result_set", [])
            entity_ids = [r["identifier"] for r in results]
            log_message(log_file, f"  {annotation_type} {annotation_id}: Found {len(entity_ids)} polymer entities")
            return entity_ids
        else:
            log_message(log_file, f"  {annotation_type} {annotation_id}: HTTP {resp.status_code} - {resp.text[:200]}")
            return []
    except Exception as e:
        log_message(log_file, f"  {annotation_type} {annotation_id}: Error - {e}")
        return []


def fetch_metadata_batch(pdb_ids: List[str], log_file: Path) -> Dict[str, Any]:
    """Fetch metadata for a batch of PDB entries using GraphQL."""

    # GraphQL query to get essential metadata
    query = """
    query($ids: [String!]!) {
        entries(entry_ids: $ids) {
            rcsb_id
            struct {
                title
            }
            rcsb_entry_info {
                resolution_combined
                experimental_method
            }
            polymer_entities {
                rcsb_id
                entity_poly {
                    pdbx_seq_one_letter_code_can
                }
                rcsb_polymer_entity_container_identifiers {
                    uniprot_ids
                }
                rcsb_entity_source_organism {
                    ncbi_scientific_name
                    ncbi_taxonomy_id
                }
            }
            nonpolymer_entities {
                pdbx_entity_nonpoly {
                    name
                    comp_id
                }
            }
        }
    }
    """

    try:
        resp = SESSION.post(
            RCSB_GRAPHQL_API,
            json={"query": query, "variables": {"ids": pdb_ids}},
            timeout=120
        )
        if resp.status_code == 200:
            return resp.json()
        else:
            log_message(log_file, f"GraphQL error: HTTP {resp.status_code}")
            return {}
    except Exception as e:
        log_message(log_file, f"GraphQL error: {e}")
        return {}


def classify_species(taxonomy_id: Optional[int], organism_name: Optional[str]) -> str:
    """Classify species into categories based on NCBI taxonomy ID."""
    if not taxonomy_id:
        return "unknown"

    # Direct mapping for common organisms
    mammal_ids = {9606, 10090, 10116, 9823, 9913, 9986, 9615, 9685, 9796}  # Human, Mouse, Rat, Pig, Cow, Rabbit, Dog, Cat, Horse

    if taxonomy_id in mammal_ids:
        return "mammals"

    # Use organism name for broader classification
    if organism_name:
        name_lower = organism_name.lower()
        if any(x in name_lower for x in ["arabidopsis", "oryza", "zea", "nicotiana", "solanum", "glycine"]):
            return "plants"
        if any(x in name_lower for x in ["escherichia", "bacillus", "pseudomonas", "streptomyces", "mycobacterium"]):
            return "bacteria"
        if any(x in name_lower for x in ["saccharomyces", "candida", "aspergillus", "fusarium", "neurospora"]):
            return "fungi"
        if any(x in name_lower for x in ["mus ", "rattus", "homo ", "sus ", "bos ", "canis", "equus", "oryctolagus"]):
            return "mammals"

    # Default classification based on taxonomy ID ranges (approximate)
    if 1 <= taxonomy_id <= 100000:
        return "bacteria"  # Very rough approximation

    return "other"


def main():
    """Main pipeline function."""

    # Setup
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = LOG_DIR / f"global_search_log_{timestamp}.txt"
    db_path = OUTPUT_DIR / "p450_global_search.db"

    log_message(log_file, "=" * 60)
    log_message(log_file, "P450 Global Search Pipeline v2")
    log_message(log_file, f"Started at: {datetime.now()}")
    log_message(log_file, "Strategy: PFAM + InterPro only (no text search)")
    log_message(log_file, "=" * 60)

    # Initialize database
    conn = init_database(db_path)
    cursor = conn.cursor()

    # Check if search already completed
    cursor.execute("SELECT value FROM progress WHERE key = 'search_completed'")
    search_completed = cursor.fetchone()

    all_entities: Set[str] = set()

    if not search_completed:
        # Step 1: Search for P450 polymer entities
        log_message(log_file, "\n" + "=" * 60)
        log_message(log_file, "STEP 1: Searching for P450 polymer entities")
        log_message(log_file, "=" * 60)

        # Search by PFAM
        for pfam_id in P450_DOMAINS["PFAM"]:
            entities = search_by_annotation("PFAM", pfam_id, log_file)
            all_entities.update(entities)
            time.sleep(0.5)

        # Search by InterPro
        for interpro_id in P450_DOMAINS["InterPro"]:
            entities = search_by_annotation("InterPro", interpro_id, log_file)
            all_entities.update(entities)
            time.sleep(0.5)

        log_message(log_file, f"\nTotal unique polymer entities: {len(all_entities)}")

        # Parse entity IDs and save to database
        for entity_full_id in all_entities:
            parts = entity_full_id.split("_")
            pdb_id = parts[0]
            entity_id = parts[1] if len(parts) > 1 else "1"
            cursor.execute(
                "INSERT OR IGNORE INTO polymer_entities (pdb_id, entity_id, source) VALUES (?, ?, ?)",
                (pdb_id, entity_id, "PFAM+InterPro")
            )

        cursor.execute("INSERT OR REPLACE INTO progress (key, value) VALUES (?, ?)",
                      ("search_completed", "true"))
        conn.commit()

        log_message(log_file, f"Saved {len(all_entities)} entities to database")
    else:
        # Load existing entities
        cursor.execute("SELECT pdb_id, entity_id FROM polymer_entities")
        for row in cursor.fetchall():
            all_entities.add(f"{row[0]}_{row[1]}")
        log_message(log_file, f"Loaded {len(all_entities)} entities from database (search already completed)")

    # Get unique PDB IDs
    pdb_ids = list(set(e.split("_")[0] for e in all_entities))
    log_message(log_file, f"Total unique PDB entries: {len(pdb_ids)}")

    # Step 2: Fetch metadata
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 2: Fetching metadata via GraphQL")
    log_message(log_file, "=" * 60)

    # Check which PDBs still need metadata
    cursor.execute("SELECT DISTINCT pdb_id FROM p450_metadata")
    fetched_pdbs = set(row[0] for row in cursor.fetchall())
    pending_pdbs = [p for p in pdb_ids if p not in fetched_pdbs]

    log_message(log_file, f"Already fetched: {len(fetched_pdbs)}, Pending: {len(pending_pdbs)}")

    # Process in batches of 50
    batch_size = 50
    total_batches = (len(pending_pdbs) + batch_size - 1) // batch_size

    for i in range(0, len(pending_pdbs), batch_size):
        batch = pending_pdbs[i:i + batch_size]
        batch_num = i // batch_size + 1

        log_message(log_file, f"Processing batch {batch_num}/{total_batches} ({len(batch)} entries)")

        result = fetch_metadata_batch(batch, log_file)

        if result and "data" in result and result["data"]["entries"]:
            for entry in result["data"]["entries"]:
                if not entry:
                    continue

                pdb_id = entry.get("rcsb_id", "")
                title = entry.get("struct", {}).get("title", "") if entry.get("struct") else ""

                entry_info = entry.get("rcsb_entry_info", {}) or {}
                resolution = entry_info.get("resolution_combined")
                # Handle if resolution is a list
                if isinstance(resolution, list):
                    resolution = resolution[0] if resolution else None
                method = entry_info.get("experimental_method", "")
                # Handle if method is a list
                if isinstance(method, list):
                    method = method[0] if method else ""

                # Get ligands
                ligands = []
                for ne in (entry.get("nonpolymer_entities") or []):
                    pdbx = ne.get("pdbx_entity_nonpoly", {}) or {}
                    comp_id = pdbx.get("comp_id", "")
                    name = pdbx.get("name", "")
                    if comp_id and comp_id not in ["HOH", "SO4", "PO4", "CL", "NA", "MG", "CA", "ZN"]:
                        ligands.append(f"{comp_id}:{name}")

                ligands_str = "; ".join(ligands) if ligands else ""

                # Process each polymer entity
                for pe in (entry.get("polymer_entities") or []):
                    entity_full_id = pe.get("rcsb_id", "")
                    if not entity_full_id:
                        continue

                    entity_id = entity_full_id.split("_")[1] if "_" in entity_full_id else "1"

                    # Get sequence
                    seq = ""
                    if pe.get("entity_poly"):
                        seq = pe["entity_poly"].get("pdbx_seq_one_letter_code_can", "") or ""

                    # Get UniProt IDs
                    uniprot_ids = []
                    container_ids = pe.get("rcsb_polymer_entity_container_identifiers", {}) or {}
                    if container_ids.get("uniprot_ids"):
                        uniprot_ids = container_ids["uniprot_ids"]

                    # Get organism info
                    organism_name = ""
                    taxonomy_id = None
                    for org in (pe.get("rcsb_entity_source_organism") or []):
                        if org:
                            organism_name = org.get("ncbi_scientific_name", "") or ""
                            taxonomy_id = org.get("ncbi_taxonomy_id")
                            break

                    # Classify species
                    species_cat = classify_species(taxonomy_id, organism_name)

                    # Save to database
                    cursor.execute("""
                        INSERT OR REPLACE INTO p450_metadata
                        (pdb_id, entity_id, title, resolution, experimental_method,
                         organism_name, ncbi_taxonomy_id, uniprot_ids, sequence,
                         sequence_length, ligands, species_category)
                        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """, (
                        pdb_id, entity_id, title, resolution, method,
                        organism_name, taxonomy_id, json.dumps(uniprot_ids), seq,
                        len(seq) if seq else 0, ligands_str, species_cat
                    ))

        conn.commit()

        # Progress report
        cursor.execute("SELECT COUNT(DISTINCT pdb_id) FROM p450_metadata")
        done_count = cursor.fetchone()[0]
        progress_pct = done_count / len(pdb_ids) * 100
        log_message(log_file, f"  Progress: {done_count}/{len(pdb_ids)} ({progress_pct:.1f}%)")

        time.sleep(0.3)  # Rate limiting

    # Step 3: Generate summary statistics
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 3: Summary Statistics")
    log_message(log_file, "=" * 60)

    cursor.execute("SELECT COUNT(*) FROM p450_metadata")
    total_entries = cursor.fetchone()[0]

    cursor.execute("SELECT species_category, COUNT(*) FROM p450_metadata GROUP BY species_category")
    species_counts = dict(cursor.fetchall())

    cursor.execute("SELECT COUNT(*) FROM p450_metadata WHERE resolution <= ?", (MAX_RESOLUTION,))
    high_quality = cursor.fetchone()[0]

    cursor.execute("SELECT COUNT(*) FROM p450_metadata WHERE ligands != ''")
    with_ligands = cursor.fetchone()[0]

    log_message(log_file, f"Total P450 entities: {total_entries}")
    log_message(log_file, f"Species distribution:")
    for cat, count in sorted(species_counts.items(), key=lambda x: -x[1]):
        log_message(log_file, f"  - {cat}: {count}")
    log_message(log_file, f"High quality (resolution <= {MAX_RESOLUTION}A): {high_quality}")
    log_message(log_file, f"With ligands/substrates: {with_ligands}")

    # Step 4: Export to CSV
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "STEP 4: Exporting to CSV")
    log_message(log_file, "=" * 60)

    csv_path = OUTPUT_DIR / "all_p450_structures.csv"
    cursor.execute("SELECT * FROM p450_metadata ORDER BY species_category, pdb_id")
    columns = [desc[0] for desc in cursor.description]

    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(columns)
        writer.writerows(cursor.fetchall())

    log_message(log_file, f"Exported to: {csv_path}")

    # Export by species
    for category in species_counts.keys():
        if category:
            cat_path = OUTPUT_DIR / f"p450_{category}.csv"
            cursor.execute("SELECT * FROM p450_metadata WHERE species_category = ?", (category,))
            with open(cat_path, "w", newline="", encoding="utf-8-sig") as f:
                writer = csv.writer(f)
                writer.writerow(columns)
                writer.writerows(cursor.fetchall())
            log_message(log_file, f"  Exported {category}: {cat_path}")

    conn.close()

    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, f"Pipeline completed at: {datetime.now()}")
    log_message(log_file, "=" * 60)


if __name__ == "__main__":
    main()
