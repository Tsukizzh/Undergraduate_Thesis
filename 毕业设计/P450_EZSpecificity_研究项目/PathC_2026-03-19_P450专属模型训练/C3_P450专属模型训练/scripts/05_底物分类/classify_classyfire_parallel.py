#!/usr/bin/env python3
"""Parallel ClassyFire batch classification. Resumes from existing cache."""

import csv
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent
DATA_DIR = PROJECT_DIR / "data" / "combined"
OUTPUT_DIR = PROJECT_DIR / "data" / "classifications"
CACHE_FILE = OUTPUT_DIR / "classyfire_cache.json"
OUTPUT_FILE = OUTPUT_DIR / "classyfire_results.csv"

CLASSYFIRE_URL = "http://classyfire.wishartlab.com/entities"
WORKERS = 8

cache_lock = threading.Lock()
progress_lock = threading.Lock()
completed = 0


def smiles_to_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        inchi = MolToInchi(mol)
        return InchiToInchiKey(inchi) if inchi else None
    except Exception:
        return None


def query_classyfire(inchikey):
    url = f"{CLASSYFIRE_URL}/{inchikey}.json"
    for attempt in range(3):
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
            return inchikey, result
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return inchikey, {"error": "not_found"}
            time.sleep(1)
        except Exception:
            time.sleep(1)
    return inchikey, {"error": "failed_after_retries"}


def worker_task(inchikey):
    global completed
    result = query_classyfire(inchikey)
    with progress_lock:
        completed += 1
    return result


def main():
    sys.stdout.reconfigure(encoding="utf-8")

    # Load compounds
    compounds = []
    with open(DATA_DIR / "global_compounds.csv", "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            compounds.append({
                "id": row["global_compound_id"],
                "name": row["canonical_name"].strip(),
                "smiles": row["canonical_smiles"].strip(),
            })

    # Generate InChIKeys
    print("Generating InChIKeys...")
    for c in compounds:
        if "*" in c["smiles"]:
            c["inchikey"] = None
        else:
            c["inchikey"] = smiles_to_inchikey(c["smiles"])

    has_key = sum(1 for c in compounds if c["inchikey"])
    print(f"InChIKey generated: {has_key}")

    # Load existing cache
    cache = {}
    if CACHE_FILE.exists():
        with open(CACHE_FILE, "r", encoding="utf-8") as f:
            cache = json.load(f)
    print(f"Existing cache: {len(cache)}")

    # Find what still needs querying
    to_query = []
    for c in compounds:
        ik = c.get("inchikey")
        if ik and ik not in cache:
            to_query.append(ik)

    # Deduplicate
    to_query = list(set(to_query))
    print(f"To query: {len(to_query)} (deduplicated InChIKeys)")
    print(f"Using {WORKERS} parallel workers")

    if to_query:
        global completed
        completed = 0
        total = len(to_query)

        with ThreadPoolExecutor(max_workers=WORKERS) as executor:
            futures = {executor.submit(worker_task, ik): ik for ik in to_query}
            batch_results = {}
            for future in as_completed(futures):
                ik, result = future.result()
                batch_results[ik] = result
                # Save cache every 100
                if completed % 100 == 0 or completed == total:
                    with cache_lock:
                        cache.update(batch_results)
                        with open(CACHE_FILE, "w", encoding="utf-8") as f:
                            json.dump(cache, f, ensure_ascii=False, indent=1)
                        batch_results = {}
                    print(f"  Progress: {completed}/{total}")

            # Final save
            cache.update(batch_results)
            with open(CACHE_FILE, "w", encoding="utf-8") as f:
                json.dump(cache, f, ensure_ascii=False, indent=1)
            print(f"  Done: {completed}/{total}")

    # Write output CSV
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
    sc_counts = Counter()
    ok = nf = err = 0
    for c in compounds:
        ik = c.get("inchikey")
        if not ik:
            continue
        r = cache.get(ik, {})
        if not r or "error" in r:
            if r and r.get("error") == "not_found":
                nf += 1
            else:
                err += 1
        else:
            ok += 1
            sc_counts[r.get("superclass", "(empty)")] += 1

    print(f"\n{'='*60}")
    print(f"ClassyFire complete!")
    print(f"Classified: {ok}, Not found: {nf}, Errors: {err}")
    print(f"\nTop 15 ClassyFire superclasses:")
    for sc, cnt in sc_counts.most_common(15):
        print(f"  {sc:<45} {cnt:>5}")


if __name__ == "__main__":
    main()
