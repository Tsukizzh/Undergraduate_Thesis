#!/usr/bin/env python3
"""
Identify true Cytochrome P450 enzymes from ESIBank enzymes.csv

Correct methodology (2026-01-02):
  - Query UniProt protein_families field for "Cytochrome P450" (not just EC number!)
  - Check InterPro domains: IPR001128, IPR036396, IPR002401
  - Use Heme cofactor as supporting evidence (not hard filter)

Input: enzymes.csv with columns: sequences, uniprots
Output:
  - true_p450_enzymes.csv: confirmed P450 enzymes
  - non_p450_enzymes.csv: confirmed non-P450 enzymes
  - p450_identification_summary.json: statistics

Based on three-way verification (Claude + Codex + Gemini):
  - Codex: Provided robust batch UniProt query logic
  - Gemini: Refined criteria (broader string match, multiple InterPro IDs)
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import os
import re
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict
from typing import Dict, List, Optional, Sequence, Set, Tuple

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# P450-related InterPro domains
P450_INTERPRO_IDS = {"IPR001128", "IPR036396", "IPR002401"}

# UniProt accession regex
ACCESSION_RE = re.compile(
    r"\b(?:"
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")\b"
)


def eprint(*args):
    print(*args, file=sys.stderr)


def now_utc_iso() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class UniProtRecord:
    accession: str
    protein_families: str
    cc_cofactor: str
    xref_interpro: str
    status: str  # "ok" | "not_found" | "error"
    retrieved_at: str


def parse_accessions_from_cell(cell: str) -> List[str]:
    if not cell:
        return []
    return list(dict.fromkeys(ACCESSION_RE.findall(cell.upper())))


def read_uniprot_ids_from_enzymes_csv(path: str) -> List[str]:
    ids: List[str] = []
    with open(path, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ids.extend(parse_accessions_from_cell(row.get("uniprots", "") or ""))
    seen: Set[str] = set()
    out: List[str] = []
    for acc in ids:
        if acc not in seen:
            seen.add(acc)
            out.append(acc)
    return out


def ensure_cache_schema(conn: sqlite3.Connection) -> None:
    conn.execute("""
        CREATE TABLE IF NOT EXISTS uniprot_cache (
          accession TEXT PRIMARY KEY,
          protein_families TEXT NOT NULL DEFAULT '',
          cc_cofactor TEXT NOT NULL DEFAULT '',
          xref_interpro TEXT NOT NULL DEFAULT '',
          status TEXT NOT NULL,
          retrieved_at TEXT NOT NULL
        )
    """)


def load_cached_records(conn: sqlite3.Connection, accessions: Sequence[str]) -> Dict[str, UniProtRecord]:
    cached: Dict[str, UniProtRecord] = {}
    if not accessions:
        return cached
    chunk_size = 900
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i:i + chunk_size]
        placeholders = ",".join(["?"] * len(chunk))
        rows = conn.execute(
            f"SELECT accession, protein_families, cc_cofactor, xref_interpro, status, retrieved_at "
            f"FROM uniprot_cache WHERE accession IN ({placeholders})",
            chunk
        ).fetchall()
        for r in rows:
            cached[r[0]] = UniProtRecord(
                accession=r[0], protein_families=r[1], cc_cofactor=r[2],
                xref_interpro=r[3], status=r[4], retrieved_at=r[5]
            )
    return cached


def upsert_records(conn: sqlite3.Connection, records: List[UniProtRecord]) -> None:
    conn.executemany("""
        INSERT INTO uniprot_cache (accession, protein_families, cc_cofactor, xref_interpro, status, retrieved_at)
        VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(accession) DO UPDATE SET
          protein_families=excluded.protein_families,
          cc_cofactor=excluded.cc_cofactor,
          xref_interpro=excluded.xref_interpro,
          status=excluded.status,
          retrieved_at=excluded.retrieved_at
    """, [(r.accession, r.protein_families, r.cc_cofactor, r.xref_interpro, r.status, r.retrieved_at) for r in records])


def is_p450_by_family(protein_families: str) -> bool:
    """Check if protein_families contains 'Cytochrome P450' (case-insensitive)"""
    if not protein_families:
        return False
    return "cytochrome p450" in protein_families.lower()


def is_p450_by_interpro(xref_interpro: str) -> bool:
    """Check if any P450-related InterPro domain is present"""
    if not xref_interpro:
        return False
    return any(ipr in xref_interpro.upper() for ipr in P450_INTERPRO_IDS)


def has_heme(cc_cofactor: str) -> bool:
    return bool(cc_cofactor) and "heme" in cc_cofactor.lower()


def http_get_with_retries(url: str, headers: Dict[str, str],
                          timeout_s: int, max_retries: int) -> str:
    for attempt in range(max_retries + 1):
        try:
            req = urllib.request.Request(url, headers=headers, method="GET")
            with urllib.request.urlopen(req, timeout=timeout_s) as resp:
                charset = resp.headers.get_content_charset() or "utf-8"
                return resp.read().decode(charset, errors="replace")
        except urllib.error.HTTPError as e:
            if e.code == 429:
                retry_after = e.headers.get("Retry-After") if e.headers else None
                sleep_s = float(retry_after) if retry_after and retry_after.isdigit() else 2 ** attempt
                eprint(f"Rate limited (429), sleeping {sleep_s:.1f}s...")
                time.sleep(min(sleep_s, 120.0))
                continue
            if 500 <= e.code < 600:
                time.sleep(2 ** attempt)
                continue
            raise
        except Exception:
            time.sleep(2 ** attempt)
    raise RuntimeError(f"Failed after {max_retries} retries")


def uniprot_fetch_batch(accessions: Sequence[str], contact: str, timeout_s: int,
                        max_retries: int) -> Dict[str, UniProtRecord]:
    retrieved_at = now_utc_iso()
    if not accessions:
        return {}

    # Build GET request URL
    fields = "accession,protein_families,cc_cofactor,xref_interpro"
    query = "(" + " OR ".join([f"accession:{a}" for a in accessions]) + ")"

    params = {
        "query": query,
        "format": "tsv",
        "fields": fields,
        "size": str(len(accessions)),
    }
    url = UNIPROT_SEARCH_URL + "?" + urllib.parse.urlencode(params)
    headers = {
        "Accept": "text/plain",
        "User-Agent": f"P450Identifier/2.0 ({contact})",
    }

    tsv = http_get_with_retries(url, headers, timeout_s, max_retries)

    lines = [ln for ln in tsv.splitlines() if ln.strip()]
    if not lines:
        return {acc: UniProtRecord(acc, "", "", "", "not_found", retrieved_at) for acc in accessions}

    reader = csv.DictReader(lines, delimiter="\t")
    results: Dict[str, UniProtRecord] = {}

    for row in reader:
        acc = (row.get("Entry") or "").strip().upper()
        if not acc:
            continue
        results[acc] = UniProtRecord(
            accession=acc,
            protein_families=(row.get("Protein families") or "").strip(),
            cc_cofactor=(row.get("Cofactor") or "").strip(),
            xref_interpro=(row.get("InterPro") or "").strip(),
            status="ok",
            retrieved_at=retrieved_at
        )

    for acc in accessions:
        if acc not in results:
            results[acc] = UniProtRecord(acc, "", "", "", "not_found", retrieved_at)

    return results


def fetch_with_auto_split(accessions: Sequence[str], contact: str, timeout_s: int,
                          max_retries: int) -> Dict[str, UniProtRecord]:
    try:
        return uniprot_fetch_batch(accessions, contact, timeout_s, max_retries)
    except urllib.error.HTTPError as e:
        if e.code in (400, 414) and len(accessions) > 1:
            mid = len(accessions) // 2
            left = fetch_with_auto_split(accessions[:mid], contact, timeout_s, max_retries)
            right = fetch_with_auto_split(accessions[mid:], contact, timeout_s, max_retries)
            left.update(right)
            return left
        raise


def main() -> int:
    ap = argparse.ArgumentParser(description="Identify true P450 enzymes via UniProt annotations")
    ap.add_argument("--input", required=True, help="Path to enzymes.csv")
    ap.add_argument("--output-dir", required=True, help="Output directory")
    ap.add_argument("--cache", required=True, help="SQLite cache file path")
    ap.add_argument("--contact", default="p450-identifier@example.com", help="Contact email for UniProt")
    ap.add_argument("--batch-size", type=int, default=150, help="Batch size for UniProt queries")
    ap.add_argument("--sleep", type=float, default=0.3, help="Sleep between batches")
    ap.add_argument("--timeout", type=int, default=60, help="HTTP timeout")
    ap.add_argument("--max-retries", type=int, default=5, help="Max retries per request")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Read all UniProt IDs
    eprint(f"Reading UniProt IDs from {args.input}...")
    accessions = read_uniprot_ids_from_enzymes_csv(args.input)
    eprint(f"Found {len(accessions)} unique UniProt IDs")

    # Setup cache
    conn = sqlite3.connect(args.cache)
    ensure_cache_schema(conn)

    # Load cached records
    cached = load_cached_records(conn, accessions)
    to_fetch = [a for a in accessions if a not in cached]
    eprint(f"Cached: {len(cached)}, To fetch: {len(to_fetch)}")

    # Fetch missing records
    fetched_count = 0
    for i in range(0, len(to_fetch), args.batch_size):
        batch = to_fetch[i:i + args.batch_size]
        eprint(f"Fetching batch {i // args.batch_size + 1}/{(len(to_fetch) - 1) // args.batch_size + 1} ({len(batch)} accessions)...")

        batch_results = fetch_with_auto_split(batch, args.contact, args.timeout, args.max_retries)
        upsert_records(conn, list(batch_results.values()))
        conn.commit()
        fetched_count += len(batch_results)

        if i + args.batch_size < len(to_fetch):
            time.sleep(args.sleep)

    eprint(f"Fetched {fetched_count} records from UniProt")

    # Load all records
    all_records = load_cached_records(conn, accessions)
    conn.close()

    # Classify enzymes
    p450_list = []
    non_p450_list = []
    stats = {
        "total_enzymes": len(accessions),
        "p450_by_family": 0,
        "p450_by_interpro": 0,
        "p450_with_heme": 0,
        "not_found_in_uniprot": 0,
    }

    for acc in accessions:
        rec = all_records.get(acc)
        if not rec or rec.status != "ok":
            stats["not_found_in_uniprot"] += 1
            continue

        by_family = is_p450_by_family(rec.protein_families)
        by_interpro = is_p450_by_interpro(rec.xref_interpro)
        heme = has_heme(rec.cc_cofactor)

        is_p450 = by_family or by_interpro

        row = {
            "uniprot_id": acc,
            "protein_families": rec.protein_families,
            "cc_cofactor": rec.cc_cofactor,
            "xref_interpro": rec.xref_interpro,
            "is_p450_by_family": by_family,
            "is_p450_by_interpro": by_interpro,
            "has_heme": heme,
        }

        if is_p450:
            p450_list.append(row)
            if by_family:
                stats["p450_by_family"] += 1
            if by_interpro:
                stats["p450_by_interpro"] += 1
            if heme:
                stats["p450_with_heme"] += 1
        else:
            non_p450_list.append(row)

    # Write outputs
    p450_path = os.path.join(args.output_dir, "true_p450_enzymes.csv")
    non_p450_path = os.path.join(args.output_dir, "non_p450_enzymes.csv")
    summary_path = os.path.join(args.output_dir, "p450_identification_summary.json")

    fieldnames = ["uniprot_id", "protein_families", "cc_cofactor", "xref_interpro",
                  "is_p450_by_family", "is_p450_by_interpro", "has_heme"]

    with open(p450_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(sorted(p450_list, key=lambda x: x["uniprot_id"]))

    with open(non_p450_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(sorted(non_p450_list, key=lambda x: x["uniprot_id"]))

    stats["total_p450"] = len(p450_list)
    stats["total_non_p450"] = len(non_p450_list)
    stats["identification_date"] = now_utc_iso()

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)

    eprint(f"\n=== P450 Identification Results ===")
    eprint(f"Total enzymes analyzed: {len(accessions)}")
    eprint(f"True P450 enzymes: {len(p450_list)}")
    eprint(f"  - Identified by protein_families: {stats['p450_by_family']}")
    eprint(f"  - Identified by InterPro: {stats['p450_by_interpro']}")
    eprint(f"  - With Heme cofactor: {stats['p450_with_heme']}")
    eprint(f"Non-P450 enzymes: {len(non_p450_list)}")
    eprint(f"Not found in UniProt: {stats['not_found_in_uniprot']}")
    eprint(f"\nOutputs written to: {args.output_dir}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
