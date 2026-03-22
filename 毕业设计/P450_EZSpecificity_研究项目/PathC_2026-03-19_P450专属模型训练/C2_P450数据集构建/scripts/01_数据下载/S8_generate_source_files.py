"""
Generate S8 PlantP450DB standardized source files:
  enzymes.csv, compounds.csv, interactions.csv, unresolved.csv
Only include pairs where enzyme has sequence AND substrate has SMILES.
"""
import json, csv, hashlib, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
PLANT_DIR = os.path.join(BASE, "downloads", "PlantP450DB")
OUT_DIR = os.path.join(BASE, "data", "sources", "Source_PlantP450DB")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Load data ──
# 1. SMILES cache
smiles_map = {}  # name -> {smiles, inchikey, cid}
with open(os.path.join(PLANT_DIR, "smiles_cache.jsonl"), encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        if e.get("smiles"):
            smiles_map[e["name"]] = e

# 2. UniProt results
uniprot_map = {}  # (cyp_name, species) -> {uniprot_id, sequence, ...}
with open(os.path.join(PLANT_DIR, "uniprot_results.jsonl"), encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        if e.get("sequence"):
            uniprot_map[(e["cyp_name"], e["species"])] = e

# 3. Final merged results
results = []
with open(os.path.join(PLANT_DIR, "final_merged_results.jsonl"), encoding="utf-8") as f:
    for line in f:
        results.append(json.loads(line))

print(f"SMILES available: {len(smiles_map)}")
print(f"Enzymes with sequence: {len(uniprot_map)}")
print(f"Total entries: {len(results)}")

# ── Build enzymes and compounds ──
# Unique enzymes: keyed by (cyp_name, species) with sequence
enzyme_key_to_id = {}  # (cyp_name, species) -> enzyme_id
enzyme_rows = []

# Unique compounds: keyed by SMILES
smiles_to_compound_id = {}  # smiles -> compound_id
compound_rows = []

# Interactions
interaction_rows = []

# Unresolved
unresolved_rows = []

for entry in results:
    cyp = entry.get("cyp_name", "")
    species = entry.get("species", "")
    subs = entry.get("substrates", [])
    status = entry.get("status", "")
    doi = entry.get("evidence_doi", "")
    notes = entry.get("notes", "")
    has_pmid = 1 if doi else 0

    # Check QC status
    qc = "none"
    if "OPUS QC: VERIFIED" in notes:
        qc = "verified"
    elif "OPUS QC: PLAUSIBLE" in notes:
        qc = "plausible"

    # Check if enzyme has sequence
    enzyme_key = (cyp, species)
    uni = uniprot_map.get(enzyme_key)

    if not subs:
        unresolved_rows.append({
            "cyp_name": cyp, "species": species,
            "reason": f"no_substrate ({status})", "num_pairs_lost": 0
        })
        continue

    if not uni:
        # Count how many usable substrates this enzyme has
        usable_count = sum(1 for s in subs
                          if (s if isinstance(s, str) else s.get("name", "")) in smiles_map)
        if usable_count > 0:
            unresolved_rows.append({
                "cyp_name": cyp, "species": species,
                "reason": "no_sequence", "num_pairs_lost": usable_count
            })
        continue

    # Enzyme has sequence - register it
    if enzyme_key not in enzyme_key_to_id:
        eid = len(enzyme_rows)
        enzyme_key_to_id[enzyme_key] = eid

        uniprot_id = uni.get("uniprot_id", "")
        sequence = uni.get("sequence", "")
        if not uniprot_id and sequence:
            h = hashlib.sha256(sequence.encode()).hexdigest()[:16]
            uniprot_id = f"SEQHASH_{h}"

        enzyme_rows.append({
            "enzyme_id": eid,
            "uniprot_id": uniprot_id,
            "p450_symbol": cyp,
            "species": species,
            "species_class": "Plant",
            "ec_number": "",
            "sequence": sequence,
            "sequence_length": len(sequence)
        })

    eid = enzyme_key_to_id[enzyme_key]

    # Process each substrate
    for s in subs:
        sub_name = s if isinstance(s, str) else s.get("name", "")
        if sub_name not in smiles_map:
            continue

        smiles = smiles_map[sub_name]["smiles"]

        # Register compound (deduplicate by SMILES)
        if smiles not in smiles_to_compound_id:
            cid = len(compound_rows)
            smiles_to_compound_id[smiles] = cid
            compound_rows.append({
                "compound_id": cid,
                "smiles": smiles,
                "name": sub_name
            })

        cid = smiles_to_compound_id[smiles]

        # Check for duplicate interaction (same enzyme + same compound)
        pair_key = (eid, cid)
        # Simple dedup: skip if already exists
        existing = [r for r in interaction_rows if r["enzyme_id"] == eid and r["compound_id"] == cid]
        if existing:
            continue

        interaction_rows.append({
            "interaction_id": len(interaction_rows),
            "enzyme_id": eid,
            "compound_id": cid,
            "label": 1,
            "source": "S8_PlantP450DB",
            "quality_tier": "B",
            "num_reactions": 1,
            "has_pmid": has_pmid,
            "has_products": 0
        })

# ── Write files ──
# enzymes.csv
with open(os.path.join(OUT_DIR, "enzymes.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["enzyme_id", "uniprot_id", "p450_symbol", "species",
                                       "species_class", "ec_number", "sequence", "sequence_length"])
    w.writeheader()
    w.writerows(enzyme_rows)

# compounds.csv
with open(os.path.join(OUT_DIR, "compounds.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["compound_id", "smiles", "name"])
    w.writeheader()
    w.writerows(compound_rows)

# interactions.csv
with open(os.path.join(OUT_DIR, "interactions.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["interaction_id", "enzyme_id", "compound_id", "label",
                                       "source", "quality_tier", "num_reactions", "has_pmid", "has_products"])
    w.writeheader()
    w.writerows(interaction_rows)

# unresolved.csv
with open(os.path.join(OUT_DIR, "unresolved.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["cyp_name", "species", "reason", "num_pairs_lost"])
    w.writeheader()
    w.writerows(unresolved_rows)

print(f"\n=== S8 Source Files Generated ===")
print(f"enzymes.csv:      {len(enzyme_rows)} enzymes")
print(f"compounds.csv:    {len(compound_rows)} compounds")
print(f"interactions.csv: {len(interaction_rows)} pairs")
print(f"unresolved.csv:   {len(unresolved_rows)} dropped entries")
print(f"\nLost pairs (no sequence): {sum(r['num_pairs_lost'] for r in unresolved_rows if r['reason']=='no_sequence')}")
print(f"Output: {OUT_DIR}")
