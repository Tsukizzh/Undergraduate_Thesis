#!/usr/bin/env python3
"""Build a resumable UniProt -> RCSB candidate structure manifest for Q1 EXP002.

This script does not download structure files. It queries official RCSB APIs to
find candidate PDB entries for each UniProt accession, then annotates whether
those entries contain heme/iron-like non-polymer components.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
from pathlib import Path
from typing import Iterable

import pandas as pd
import requests


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"

HEME_OR_FE_COMPONENTS = {
    "HEM",
    "HEC",
    "HEA",
    "HEB",
    "HEO",
    "HAS",
    "HDD",
    "DHE",
    "FE",
    "FES",
    "SF4",
    "F3S",
    "FEO",
}

UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def chunks(items: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def normalize_accession(value: str) -> tuple[str, str] | None:
    """Return (original, accession_for_rcsb) or None for clearly invalid values."""
    raw = str(value).strip()
    if not raw:
        return None
    # RCSB indexes UniProt accessions, not isoform suffixes such as P12345-2.
    query = raw.split("-")[0].upper()
    if not UNIPROT_ACCESSION_RE.match(query):
        return None
    return raw, query


def post_json(session: requests.Session, url: str, payload: dict, retries: int = 4) -> dict:
    last_error = None
    for attempt in range(retries):
        try:
            resp = session.post(url, json=payload, timeout=60)
            if resp.status_code == 204:
                return {}
            if resp.status_code == 400:
                raise RuntimeError(f"HTTP 400 Bad Request: {resp.text[:500]}")
            resp.raise_for_status()
            return resp.json()
        except Exception as exc:  # noqa: BLE001 - record retry context for manifest building.
            last_error = exc
            time.sleep(2 + attempt * 3)
    raise RuntimeError(f"POST failed after {retries} retries: {url}: {last_error}")


def search_polymer_entities(session: requests.Session, accessions: list[str]) -> list[str]:
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "in",
                "value": accessions,
            },
        },
        "return_type": "polymer_entity",
        "request_options": {"return_all_hits": True},
    }
    data = post_json(session, SEARCH_URL, query)
    return [r["identifier"] for r in data.get("result_set", [])]


def search_polymer_entities_resilient(
    session: requests.Session, accessions: list[str]
) -> tuple[list[str], list[str]]:
    """Search RCSB and isolate bad accession values instead of failing a full chunk."""
    if not accessions:
        return [], []
    try:
        return search_polymer_entities(session, accessions), []
    except Exception:
        if len(accessions) == 1:
            return [], accessions
        mid = len(accessions) // 2
        left_hits, left_bad = search_polymer_entities_resilient(session, accessions[:mid])
        right_hits, right_bad = search_polymer_entities_resilient(session, accessions[mid:])
        return left_hits + right_hits, left_bad + right_bad


def fetch_polymer_entities(session: requests.Session, entity_ids: list[str]) -> list[dict]:
    query = """
    query($ids:[String!]!) {
      polymer_entities(entity_ids:$ids) {
        rcsb_id
        rcsb_polymer_entity {
          pdbx_description
        }
        rcsb_polymer_entity_container_identifiers {
          entry_id
          entity_id
          asym_ids
          auth_asym_ids
          uniprot_ids
          reference_sequence_identifiers {
            database_accession
            database_name
            entity_sequence_coverage
            reference_sequence_coverage
          }
        }
      }
    }
    """
    rows: list[dict] = []
    for block in chunks(entity_ids, 100):
        data = post_json(session, GRAPHQL_URL, {"query": query, "variables": {"ids": block}})
        rows.extend(data.get("data", {}).get("polymer_entities", []) or [])
    return rows


def fetch_entry_ligands(session: requests.Session, entry_ids: list[str]) -> dict[str, dict]:
    query = """
    query($ids:[String!]!) {
      entries(entry_ids:$ids) {
        rcsb_id
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            nonpolymer_comp_id
          }
          nonpolymer_comp {
            chem_comp {
              id
              name
              formula
            }
          }
        }
      }
    }
    """
    out: dict[str, dict] = {}
    for block in chunks(entry_ids, 50):
        data = post_json(session, GRAPHQL_URL, {"query": query, "variables": {"ids": block}})
        for entry in data.get("data", {}).get("entries", []) or []:
            comp_ids: list[str] = []
            comp_names: list[str] = []
            has_heme_or_fe = False
            for nonpoly in entry.get("nonpolymer_entities") or []:
                ident = nonpoly.get("rcsb_nonpolymer_entity_container_identifiers") or {}
                comp_id = str(ident.get("nonpolymer_comp_id") or "").upper()
                comp = ((nonpoly.get("nonpolymer_comp") or {}).get("chem_comp") or {})
                formula = str(comp.get("formula") or "")
                name = str(comp.get("name") or "")
                if comp_id:
                    comp_ids.append(comp_id)
                    comp_names.append(name)
                if comp_id in HEME_OR_FE_COMPONENTS or " Fe" in f" {formula} " or formula.startswith("Fe"):
                    has_heme_or_fe = True
            out[entry["rcsb_id"]] = {
                "ligand_comp_ids": ";".join(sorted(set(comp_ids))),
                "ligand_names": ";".join(sorted(set(comp_names))),
                "has_heme_or_fe_component": has_heme_or_fe,
            }
    return out


def load_processed(output_csv: Path) -> set[str]:
    if not output_csv.exists() or output_csv.stat().st_size == 0:
        return set()
    processed: set[str] = set()
    with output_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get("query_uniprot")
            if acc:
                processed.add(acc)
    return processed


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure-map", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--progress-json", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=80)
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    df = pd.read_csv(args.structure_map, usecols=["uniprot"])
    raw_values = sorted({str(v).strip() for v in df["uniprot"].dropna() if str(v).strip()})
    if args.limit is not None:
        raw_values = raw_values[: args.limit]
    invalid_accessions = [v for v in raw_values if normalize_accession(v) is None]
    normalized_pairs = [normalize_accession(v) for v in raw_values]
    normalized_pairs = [p for p in normalized_pairs if p is not None]
    original_to_query = {raw: query for raw, query in normalized_pairs}
    query_to_originals: dict[str, list[str]] = {}
    for raw, query in normalized_pairs:
        query_to_originals.setdefault(query, []).append(raw)
    accessions = sorted(query_to_originals)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    processed = load_processed(args.output_csv)
    remaining = [
        acc for acc in accessions
        if not set(query_to_originals.get(acc, [acc])).issubset(processed)
    ]

    fieldnames = [
        "query_uniprot",
        "rcsb_query_uniprot",
        "status",
        "entry_id",
        "entity_id",
        "rcsb_polymer_entity_id",
        "asym_ids",
        "auth_asym_ids",
        "description",
        "entity_sequence_coverage",
        "reference_sequence_coverage",
        "ligand_comp_ids",
        "has_heme_or_fe_component",
    ]
    write_header = not args.output_csv.exists() or args.output_csv.stat().st_size == 0
    session = requests.Session()

    with args.output_csv.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        for acc in invalid_accessions:
            if acc not in processed:
                writer.writerow({
                    "query_uniprot": acc,
                    "rcsb_query_uniprot": "",
                    "status": "invalid_uniprot_format",
                    "entry_id": "",
                    "entity_id": "",
                    "rcsb_polymer_entity_id": "",
                    "asym_ids": "",
                    "auth_asym_ids": "",
                    "description": "",
                    "entity_sequence_coverage": "",
                    "reference_sequence_coverage": "",
                    "ligand_comp_ids": "",
                    "has_heme_or_fe_component": False,
                })
        f.flush()

        for block in chunks(remaining, args.chunk_size):
            block_set = set(block)
            try:
                entity_ids, bad_query_accessions = search_polymer_entities_resilient(session, block)
                for query_acc in bad_query_accessions:
                    for raw_acc in query_to_originals.get(query_acc, [query_acc]):
                        writer.writerow({
                            "query_uniprot": raw_acc,
                            "rcsb_query_uniprot": query_acc,
                            "status": "query_error",
                            "entry_id": "",
                            "entity_id": "",
                            "rcsb_polymer_entity_id": "",
                            "asym_ids": "",
                            "auth_asym_ids": "",
                            "description": "",
                            "entity_sequence_coverage": "",
                            "reference_sequence_coverage": "",
                            "ligand_comp_ids": "",
                            "has_heme_or_fe_component": False,
                        })
                entities = fetch_polymer_entities(session, entity_ids) if entity_ids else []
                entry_ids = sorted({
                    (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("entry_id")
                    for e in entities
                    if (e.get("rcsb_polymer_entity_container_identifiers") or {}).get("entry_id")
                })
                ligand_info = fetch_entry_ligands(session, entry_ids) if entry_ids else {}

                hit_accessions: set[str] = set()
                for entity in entities:
                    ident = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                    entry_id = ident.get("entry_id") or ""
                    refs = ident.get("reference_sequence_identifiers") or []
                    for ref in refs:
                        rcsb_acc = str(ref.get("database_accession") or "").upper()
                        if rcsb_acc not in block_set:
                            continue
                        hit_accessions.add(rcsb_acc)
                        lig = ligand_info.get(entry_id, {})
                        for raw_acc in query_to_originals.get(rcsb_acc, [rcsb_acc]):
                            writer.writerow({
                                "query_uniprot": raw_acc,
                                "rcsb_query_uniprot": rcsb_acc,
                                "status": "hit",
                                "entry_id": entry_id,
                                "entity_id": ident.get("entity_id") or "",
                                "rcsb_polymer_entity_id": entity.get("rcsb_id") or "",
                                "asym_ids": ";".join(ident.get("asym_ids") or []),
                                "auth_asym_ids": ";".join(ident.get("auth_asym_ids") or []),
                                "description": (entity.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "",
                                "entity_sequence_coverage": ref.get("entity_sequence_coverage"),
                                "reference_sequence_coverage": ref.get("reference_sequence_coverage"),
                                "ligand_comp_ids": lig.get("ligand_comp_ids", ""),
                                "has_heme_or_fe_component": lig.get("has_heme_or_fe_component", False),
                            })

                for acc in block:
                    if acc not in hit_accessions and acc not in set(bad_query_accessions):
                        for raw_acc in query_to_originals.get(acc, [acc]):
                            writer.writerow({
                                "query_uniprot": raw_acc,
                                "rcsb_query_uniprot": acc,
                                "status": "no_rcsb_hit",
                                "entry_id": "",
                                "entity_id": "",
                                "rcsb_polymer_entity_id": "",
                                "asym_ids": "",
                                "auth_asym_ids": "",
                                "description": "",
                                "entity_sequence_coverage": "",
                                "reference_sequence_coverage": "",
                                "ligand_comp_ids": "",
                                "has_heme_or_fe_component": False,
                            })
                f.flush()
            finally:
                processed = load_processed(args.output_csv)
                progress = {
                    "total_accessions": len(raw_values),
                    "valid_rcsb_query_accessions": len(accessions),
                    "invalid_accessions": len(invalid_accessions),
                    "processed_accessions": len(processed),
                    "remaining_accessions": max(len(raw_values) - len(processed), 0),
                    "output_csv": str(args.output_csv),
                    "updated_at": time.strftime("%F %T"),
                }
                args.progress_json.write_text(json.dumps(progress, indent=2), encoding="utf-8")
                print(json.dumps(progress), flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
