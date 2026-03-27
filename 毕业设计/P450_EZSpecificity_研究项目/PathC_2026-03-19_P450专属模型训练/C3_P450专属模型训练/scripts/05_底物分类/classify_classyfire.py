#!/usr/bin/env python3
"""
Batch ClassyFire classification via Wishart Lab API.

Workflow:
1. Load global_compounds.csv
2. Generate InChIKey via RDKit (skip wildcards)
3. Query http://classyfire.wishartlab.com/entities/{InChIKey}.json
4. Cache results and save to CSV

Usage:
    python classify_classyfire.py [--resume]
"""

import csv
import json
import time
import urllib.request
import argparse
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "combined"
OUTPUT_DIR = PROJECT_DIR / "data" / "classifications"
CACHE_FILE = OUTPUT_DIR / "classyfire_cache.json"
OUTPUT_FILE = OUTPUT_DIR / "classyfire_results.csv"

CLASSYFIRE_URL = "http://classyfire.wishartlab.com/entities"
REQUEST_DELAY = 0.3  # seconds between requests


def load_compounds(csv_path):
    compounds = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            compounds.append({
                "id": row["global_compound_id"],
                "name": row["canonical_name"].strip(),
                "smiles": row["canonical_smiles"].strip(),
            })
    return compounds


def smiles_to_inchikey(smiles):
    """Convert SMILES to InChIKey via RDKit. Returns None on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        inchi = MolToInchi(mol)
        if inchi is None:
            return None
        return InchiToInchiKey(inchi)
    except Exception:
        return None


def query_classyfire(inchikey, retries=3):
    """Query ClassyFire by InChIKey. Returns dict or None."""
    url = f"{CLASSYFIRE_URL}/{inchikey}.json"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "Mozilla/5.0",
                "Accept": "application/json",
            })
            resp = urllib.request.urlopen(req, timeout=20)
            data = json.loads(resp.read().decode("utf-8"))
            result = {}
            for level in ["kingdom", "superclass", "class", "subclass", "direct_parent"]:
                v = data.get(level, {})
                result[level] = v.get("name", "") if isinstance(v, dict) else ""
            result["molecular_framework"] = data.get("molecular_framework", "")
            return result
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"error": "not_found"}
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return {"error": str(e)}
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return {"error": str(e)}
    return None


def load_cache(path):
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(cache, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    compounds = load_compounds(DATA_DIR / "global_compounds.csv")
    print(f"Loaded {len(compounds)} compounds")

    # Generate InChIKeys
    print("Generating InChIKeys...")
    for c in compounds:
        if "*" in c["smiles"]:
            c["inchikey"] = None
        else:
            c["inchikey"] = smiles_to_inchikey(c["smiles"])

    has_key = sum(1 for c in compounds if c["inchikey"])
    no_key = sum(1 for c in compounds if not c["inchikey"])
    print(f"InChIKey generated: {has_key}, failed/wildcard: {no_key}")

    # Load cache
    cache = load_cache(CACHE_FILE) if args.resume else {}
    print(f"Cached: {len(cache)}")

    # Query ClassyFire
    to_query = [c for c in compounds if c["inchikey"] and c["inchikey"] not in cache]
    print(f"To query: {len(to_query)}")

    if to_query:
        print("Querying ClassyFire API...")
        for i, c in enumerate(to_query):
            result = query_classyfire(c["inchikey"])
            cache[c["inchikey"]] = result
            if (i + 1) % 50 == 0 or i == len(to_query) - 1:
                save_cache(cache, CACHE_FILE)
                status = "ok" if result and "error" not in result else "error"
                print(f"  Progress: {i+1}/{len(to_query)} ({status})")
            time.sleep(REQUEST_DELAY)

    # Write output
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "global_compound_id", "canonical_name", "canonical_smiles", "inchikey",
            "cf_kingdom", "cf_superclass", "cf_class", "cf_subclass",
            "cf_direct_parent", "cf_molecular_framework",
        ])
        for c in compounds:
            ik = c.get("inchikey", "")
            r = cache.get(ik, {}) if ik else {}
            writer.writerow([
                c["id"], c["name"], c["smiles"], ik or "",
                r.get("kingdom", ""), r.get("superclass", ""),
                r.get("class", ""), r.get("subclass", ""),
                r.get("direct_parent", ""), r.get("molecular_framework", ""),
            ])

    # Summary
    from collections import Counter
    superclass_counts = Counter()
    errors = 0
    not_found = 0
    for c in compounds:
        ik = c.get("inchikey")
        if not ik:
            continue
        r = cache.get(ik, {})
        if not r or "error" in r:
            if r and r.get("error") == "not_found":
                not_found += 1
            else:
                errors += 1
        else:
            superclass_counts[r.get("superclass", "(empty)")] += 1

    print(f"\n{'='*60}")
    print(f"ClassyFire complete! Output: {OUTPUT_FILE}")
    print(f"Classified: {sum(superclass_counts.values())}, Not found: {not_found}, Errors: {errors}")
    print(f"\n--- ClassyFire superclass distribution ---")
    for sc, count in superclass_counts.most_common(20):
        print(f"  {sc:<45} {count:>5}")


if __name__ == "__main__":
    main()
