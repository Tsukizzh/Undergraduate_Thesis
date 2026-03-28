#!/usr/bin/env python3
"""
Batch query ClassyFire via Wishart Lab InChIKey endpoint for all compounds.
Uses RDKit to convert SMILES → InChIKey, then queries ClassyFire.
Supports resume via cache file.

Usage:
    python query_classyfire_batch.py [--workers 8]
"""

import csv
import json
import sys
import time
import argparse
import urllib.request
import urllib.parse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey
from rdkit import RDLogger
RDLogger.logger().setLevel(RDLogger.ERROR)

sys.stdout.reconfigure(encoding="utf-8")

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "05_底物分类"
INPUT_CSV = DATA_DIR / "substrate_classifications_corrected.csv"
CACHE_FILE = DATA_DIR / "classyfire_inchikey_cache.json"
OUTPUT_CSV = DATA_DIR / "classyfire_full_results.csv"

CF_URL = "http://classyfire.wishartlab.com/entities/{}.json"


def load_cache():
    if CACHE_FILE.exists():
        return json.load(open(CACHE_FILE, "r", encoding="utf-8"))
    return {}


def save_cache(cache):
    with open(CACHE_FILE, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=1)


def smiles_to_inchikey(smiles):
    if not smiles or "*" in smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        inchi = MolToInchi(mol)
        return InchiToInchiKey(inchi) if inchi else None
    except Exception:
        return None


def query_classyfire(inchikey, retries=3):
    url = CF_URL.format(inchikey)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            resp = urllib.request.urlopen(req, timeout=20)
            data = json.loads(resp.read().decode("utf-8"))
            return {
                "kingdom": data.get("kingdom", {}).get("name", ""),
                "superclass": data.get("superclass", {}).get("name", ""),
                "class": data.get("class", {}).get("name", ""),
                "subclass": data.get("subclass", {}).get("name", "") if data.get("subclass") else "",
                "direct_parent": data.get("direct_parent", {}).get("name", ""),
                "description": data.get("description", ""),
                "status": "found",
            }
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return {"status": "not_found"}
            if e.code == 429:
                time.sleep(5 * (attempt + 1))
            elif attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return {"status": "error", "error": str(e)}
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return {"status": "error", "error": str(e)}
    return {"status": "error"}


def query_one(args):
    cid, smiles, inchikey, cache = args
    if inchikey in cache:
        return cid, inchikey, cache[inchikey]
    result = query_classyfire(inchikey)
    return cid, inchikey, result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    rows = list(csv.DictReader(open(INPUT_CSV, "r", encoding="utf-8-sig")))
    cache = load_cache()
    print(f"Loaded {len(rows)} compounds, cache has {len(cache)} entries")

    # Generate InChIKeys
    tasks = []
    skip_wildcard = 0
    skip_fail = 0
    for r in rows:
        smi = r["canonical_smiles"]
        ik = smiles_to_inchikey(smi)
        if ik is None:
            if "*" in smi:
                skip_wildcard += 1
            else:
                skip_fail += 1
            continue
        if ik in cache:
            continue
        tasks.append((r["global_compound_id"], smi, ik, cache))

    already = len(cache)
    print(f"Already cached: {already}, To query: {len(tasks)}, Wildcards: {skip_wildcard}, Failed: {skip_fail}")

    if not tasks:
        print("Nothing to query.")
    else:
        print(f"Querying ClassyFire with {args.workers} workers...")
        done = 0
        found = 0
        not_found = 0
        errors = 0

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(query_one, t): t for t in tasks}
            for future in as_completed(futures):
                cid, ik, result = future.result()
                cache[ik] = result
                done += 1
                if result.get("status") == "found":
                    found += 1
                elif result.get("status") == "not_found":
                    not_found += 1
                else:
                    errors += 1

                if done % 100 == 0 or done == len(tasks):
                    save_cache(cache)
                    print(f"  Progress: {done}/{len(tasks)} (found={found}, not_found={not_found}, error={errors})")

        save_cache(cache)
        print(f"\nDone. found={found}, not_found={not_found}, error={errors}")

    # Write full results CSV
    print(f"\nWriting results to {OUTPUT_CSV}")
    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "global_compound_id", "canonical_smiles", "inchikey",
            "cf_status", "cf_kingdom", "cf_superclass", "cf_class",
            "cf_subclass", "cf_direct_parent",
        ])
        found_total = 0
        for r in rows:
            smi = r["canonical_smiles"]
            ik = smiles_to_inchikey(smi)
            entry = cache.get(ik, {}) if ik else {}
            if entry.get("status") == "found":
                found_total += 1
            writer.writerow([
                r["global_compound_id"], smi, ik or "",
                entry.get("status", "no_inchikey"),
                entry.get("kingdom", ""),
                entry.get("superclass", ""),
                entry.get("class", ""),
                entry.get("subclass", ""),
                entry.get("direct_parent", ""),
            ])

    print(f"Total found: {found_total}/{len(rows)}")


if __name__ == "__main__":
    main()
