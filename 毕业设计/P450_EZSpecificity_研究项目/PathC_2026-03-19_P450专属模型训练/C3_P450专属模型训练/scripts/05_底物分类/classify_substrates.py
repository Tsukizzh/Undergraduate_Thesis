#!/usr/bin/env python3
"""
Classify P450 substrates using NPClassifier API.

Pipeline:
1. Load global_compounds.csv
2. Split into exact SMILES vs wildcard/template SMILES
3. For wildcards: use curated name→class dictionary
4. For exact: query NPClassifier API (with caching + resume)
5. Remap raw NPClassifier labels to controlled vocabulary (~20 classes)
6. Output: substrate_classifications.csv

Usage:
    python classify_substrates.py [--resume] [--dry-run]
"""

import csv
import json
import time
import urllib.request
import urllib.parse
import argparse
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent.parent  # C2_P450数据集构建/
DATA_DIR = PROJECT_DIR / "data" / "combined"
OUTPUT_DIR = PROJECT_DIR / "data" / "classifications"
CACHE_FILE = OUTPUT_DIR / "npclassifier_cache.json"
OUTPUT_FILE = OUTPUT_DIR / "substrate_classifications.csv"

NPCLASSIFIER_URL = "https://npclassifier.gnps2.org/classify"
REQUEST_DELAY = 0.15  # seconds between API calls (stay under rate limit)

# ── Wildcard name→class mapping (15 template compounds from P450Rdb) ──
WILDCARD_NAME_MAP = {
    "A Steroid": ("Terpenoids", "Steroids", "Steroid"),
    "A 14Alpha-Methyl Steroid": ("Terpenoids", "Steroids", "Steroid"),
    "A 14alpha-formyl steroid": ("Terpenoids", "Steroids", "Steroid"),
    "A 14alpha-hydroxymethyl steroid": ("Terpenoids", "Steroids", "Steroid"),
    "A 1,2-Saturated Fatty Acid": ("Fatty acids", "Fatty acyls", "Fatty_acid"),
    "2,3-Saturated Fatty Acid": ("Fatty acids", "Fatty acyls", "Fatty_acid"),
    "A Jasmonyl-L-Amino Acid": ("Terpenoids", "Jasmonoids", "Jasmonoid"),
    "A 12-Hydroxyjasmonyl-L-Alpha-Amino Acid": ("Terpenoids", "Jasmonoids", "Jasmonoid"),
    "An Organic Molecule": ("Other", "Other", "Unclassified"),
    "A Flavanone": ("Shikimates and Phenylpropanoids", "Flavonoids", "Flavonoid"),
    "An Omega-Methyl-Long-Chain Fatty Acid": ("Fatty acids", "Fatty acyls", "Fatty_acid"),
    "An (Omega-1)-Ethyl Fatty Acid": ("Fatty acids", "Fatty acyls", "Fatty_acid"),
}

# Wildcards with mojibake names (encoding-corrupted apostrophes in CSV)
# These won't match exact dict keys, so we handle them via keyword fallback.
# Keyword order matters: "isoflavon" must be checked BEFORE "flavon".

# ── NPClassifier Superclass → Controlled label remapping ──────────────
# NPClassifier has 70 superclasses. We remap to ~20 ML-friendly labels.
SUPERCLASS_REMAP = {
    # Terpenoids
    "Monoterpenoids": "Monoterpenoid",
    "Sesquiterpenoids": "Sesquiterpenoid",
    "Diterpenoids": "Diterpenoid",
    "Sesterterpenoids": "Sesterterpenoid",
    "Triterpenoids": "Triterpenoid",
    "Steroids": "Steroid",
    "Jasmonoids": "Jasmonoid",
    "Prenol lipids": "Terpenoid_other",
    # Polyketides / Phenylpropanoids
    "Flavonoids": "Flavonoid",
    "Isoflavonoids": "Isoflavonoid",
    "Chalcones and dihydrochalcones": "Chalcone",
    "Stilbenes": "Phenylpropanoid",
    "Lignans": "Phenylpropanoid",
    "Coumarins": "Coumarin",
    "Chromones": "Phenylpropanoid",
    "Cinnamic acids and derivatives": "Phenylpropanoid",
    "Phenylpropanoids": "Phenylpropanoid",
    # Alkaloids
    "Tyrosine alkaloids": "Alkaloid",
    "Tryptophan alkaloids": "Alkaloid",
    "Pseudoalkaloids": "Alkaloid",
    "Histidine alkaloids": "Alkaloid",
    "Nicotinic acid alkaloids": "Alkaloid",
    "Ornithine alkaloids": "Alkaloid",
    "Lysine alkaloids": "Alkaloid",
    "Purine alkaloids": "Alkaloid",
    "Steroidal alkaloids": "Alkaloid",
    "Terpenoid alkaloids": "Alkaloid",
    # Fatty acids
    "Fatty acyls": "Fatty_acid",
    "Fatty esters": "Fatty_acid",
    "Fatty amides": "Fatty_acid",
    "Fatty alcohols": "Fatty_acid",
    "Glycerolipids": "Lipid",
    "Glycerophospholipids": "Lipid",
    "Sphingolipids": "Lipid",
    "Saccharolipids": "Lipid",
    # Amino acids
    "Amino acids and peptides": "Amino_acid",
    "Small peptides": "Amino_acid",
    # Others
    "Carbohydrates": "Carbohydrate",
    "Nucleosides": "Nucleoside",
    "Polyketides": "Polyketide",
    "Macrolides": "Macrolide",
    "Polyphenols": "Polyphenol",
    "Porphyrins": "Porphyrin",
}

# Pathway-level fallback (when superclass is empty or unmapped)
PATHWAY_REMAP = {
    "Terpenoids": "Terpenoid_other",
    "Shikimates and Phenylpropanoids": "Phenylpropanoid",
    "Fatty acids": "Fatty_acid",
    "Alkaloids": "Alkaloid",
    "Polyketides": "Polyketide",
    "Amino acids and Peptides": "Amino_acid",
    "Carbohydrates": "Carbohydrate",
}


def load_compounds(csv_path):
    """Load global_compounds.csv, return list of dicts."""
    compounds = []
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            compounds.append({
                "id": row["global_compound_id"],
                "name": row["canonical_name"].strip(),
                "smiles": row["canonical_smiles"].strip(),
                "sources": row["sources"],
            })
    return compounds


def load_cache(cache_path):
    """Load cached NPClassifier results."""
    if cache_path.exists():
        with open(cache_path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def save_cache(cache, cache_path):
    """Save NPClassifier cache to disk."""
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as f:
        json.dump(cache, f, ensure_ascii=False, indent=1)


def query_npclassifier(smiles, retries=3):
    """Query NPClassifier API. Returns dict or None on failure."""
    encoded = urllib.parse.quote(smiles, safe="")
    url = f"{NPCLASSIFIER_URL}?smiles={encoded}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            resp = urllib.request.urlopen(req, timeout=20)
            data = json.loads(resp.read().decode("utf-8"))
            return {
                "pathway": data.get("pathway_results", []),
                "superclass": data.get("superclass_results", []),
                "class": data.get("class_results", []),
                "isglycoside": data.get("isglycoside", False),
            }
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return {"pathway": [], "superclass": [], "class": [], "error": str(e)}
    return None


def derive_controlled_label(npc_result):
    """Map NPClassifier output to controlled vocabulary label."""
    if npc_result is None or "error" in npc_result:
        return "Unclassified", "error"

    superclasses = npc_result.get("superclass", [])
    pathways = npc_result.get("pathway", [])

    # Try superclass-level mapping first
    for sc in superclasses:
        if sc in SUPERCLASS_REMAP:
            return SUPERCLASS_REMAP[sc], "superclass"

    # Fallback to pathway-level
    for pw in pathways:
        if pw in PATHWAY_REMAP:
            return PATHWAY_REMAP[pw], "pathway"

    # If NPClassifier returned something but we can't map it
    if superclasses:
        return superclasses[0], "superclass_unmapped"
    if pathways:
        return pathways[0], "pathway_unmapped"

    return "Unclassified", "empty"


def classify_wildcard(compound):
    """Classify wildcard compounds by name."""
    name = compound["name"]
    # Try exact match first
    if name in WILDCARD_NAME_MAP:
        pw, sc, label = WILDCARD_NAME_MAP[name]
        return {
            "pathway": [pw], "superclass": [sc], "class": [label],
            "isglycoside": False,
        }, label, "name_exact"

    # Try fuzzy keyword match (order matters: specific before general)
    name_lower = name.lower()
    keyword_rules = [
        ("isoflavon", "Isoflavonoid"),  # must be before "flavon"
        ("flavon", "Flavonoid"),
        ("flavanon", "Flavonoid"),
        ("steroid", "Steroid"),
        ("sterol", "Steroid"),
        ("fatty acid", "Fatty_acid"),
        ("fatty acyl", "Fatty_acid"),
        ("terpen", "Terpenoid_other"),
        ("alkaloid", "Alkaloid"),
        ("coumarin", "Coumarin"),
        ("jasmon", "Jasmonoid"),
    ]
    for kw, label in keyword_rules:
        if kw in name_lower:
            return {
                "pathway": [], "superclass": [], "class": [label],
                "isglycoside": False,
            }, label, "name_keyword"

    return {
        "pathway": [], "superclass": [], "class": [],
        "isglycoside": False,
    }, "Unclassified", "wildcard_unknown"


def main():
    parser = argparse.ArgumentParser(description="Classify P450 substrates")
    parser.add_argument("--resume", action="store_true", help="Resume from cache")
    parser.add_argument("--dry-run", action="store_true", help="Only show stats, don't query API")
    args = parser.parse_args()

    # Load data
    compounds_path = DATA_DIR / "global_compounds.csv"
    print(f"Loading compounds from {compounds_path}")
    compounds = load_compounds(compounds_path)
    print(f"Total compounds: {len(compounds)}")

    # Split wildcard vs exact
    wildcards = [c for c in compounds if "*" in c["smiles"]]
    exacts = [c for c in compounds if "*" not in c["smiles"]]
    print(f"Exact SMILES: {len(exacts)}, Wildcard: {len(wildcards)}")

    # Load cache
    cache = load_cache(CACHE_FILE) if args.resume else {}
    print(f"Cached results: {len(cache)}")

    if args.dry_run:
        to_query = [c for c in exacts if c["smiles"] not in cache]
        print(f"Would query NPClassifier for {len(to_query)} compounds")
        return

    # Classify wildcards by name
    results = {}
    for c in wildcards:
        npc, label, source = classify_wildcard(c)
        results[c["id"]] = {
            "npc_raw": npc,
            "controlled_label": label,
            "label_source": source,
        }
    print(f"Classified {len(wildcards)} wildcard compounds by name")

    # Query NPClassifier for exact SMILES
    to_query = [c for c in exacts if c["smiles"] not in cache]
    already_cached = [c for c in exacts if c["smiles"] in cache]
    print(f"Already cached: {len(already_cached)}, To query: {len(to_query)}")

    if to_query:
        print(f"Querying NPClassifier API for {len(to_query)} compounds...")
        for i, c in enumerate(to_query):
            npc = query_npclassifier(c["smiles"])
            cache[c["smiles"]] = npc
            if (i + 1) % 50 == 0 or i == len(to_query) - 1:
                save_cache(cache, CACHE_FILE)
                print(f"  Progress: {i+1}/{len(to_query)} "
                      f"({'error' if npc and 'error' in npc else 'ok'})")
            time.sleep(REQUEST_DELAY)

    # Derive controlled labels for exact compounds
    for c in exacts:
        npc = cache.get(c["smiles"], {})
        label, source = derive_controlled_label(npc)
        results[c["id"]] = {
            "npc_raw": npc,
            "controlled_label": label,
            "label_source": source,
        }

    # Write output CSV
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "global_compound_id", "canonical_name", "canonical_smiles",
            "npc_pathway", "npc_superclass", "npc_class", "is_glycoside",
            "controlled_label", "label_source",
        ])
        for c in compounds:
            r = results.get(c["id"], {})
            npc = r.get("npc_raw", {})
            writer.writerow([
                c["id"],
                c["name"],
                c["smiles"],
                "|".join(npc.get("pathway", [])),
                "|".join(npc.get("superclass", [])),
                "|".join(npc.get("class", [])),
                npc.get("isglycoside", False),
                r.get("controlled_label", "Unclassified"),
                r.get("label_source", "none"),
            ])

    # Print summary
    from collections import Counter
    label_counts = Counter(r["controlled_label"] for r in results.values())
    source_counts = Counter(r["label_source"] for r in results.values())

    print(f"\n{'='*60}")
    print(f"Classification complete! Output: {OUTPUT_FILE}")
    print(f"Total classified: {len(results)}")
    print(f"\n--- Label distribution ---")
    for label, count in label_counts.most_common():
        print(f"  {label:<25} {count:>5}  ({100*count/len(results):.1f}%)")
    print(f"\n--- Classification source ---")
    for src, count in source_counts.most_common():
        print(f"  {src:<25} {count:>5}")


if __name__ == "__main__":
    main()
