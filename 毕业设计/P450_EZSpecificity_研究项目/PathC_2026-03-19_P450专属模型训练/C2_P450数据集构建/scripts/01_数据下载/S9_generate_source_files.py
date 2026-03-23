"""
Generate S9 PCPD standardized source files:
  enzymes.csv, compounds.csv, interactions.csv, unresolved.csv
"""
import json, csv, os, re

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PCPD")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "sources", "Source_PCPD")
os.makedirs(OUT_DIR, exist_ok=True)

# Load data
with open(os.path.join(BASE, "fasta_parsed.json"), encoding="utf-8") as f:
    fasta_data = json.load(f)

with open(os.path.join(BASE, "resource_new.json"), encoding="utf-8") as f:
    json_entries = {e["ID"]: e for e in json.load(f)}

smiles_map = {}
with open(os.path.join(BASE, "smiles_cache.jsonl"), encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        if e.get("smiles"):
            smiles_map[e["name"]] = e

extractions = {}
with open(os.path.join(BASE, "substrate_extraction_cleaned.jsonl"), encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        extractions[e["cyp_id"]] = e

# Build enzymes, compounds, interactions
enzyme_key_to_id = {}
enzyme_rows = []
smiles_to_cid = {}
compound_rows = []
interaction_rows = []
unresolved_rows = []

for cyp_id, ext in extractions.items():
    if ext["status"] != "specific":
        unresolved_rows.append({
            "cyp_id": cyp_id, "kingdom": ext.get("kingdom", ""),
            "reason": ext["status"], "num_pairs_lost": 0
        })
        continue

    # Check FASTA
    fasta = fasta_data.get(cyp_id)
    if not fasta or not fasta.get("sequence"):
        usable = sum(1 for s in ext.get("substrates_specific", []) if s in smiles_map)
        unresolved_rows.append({
            "cyp_id": cyp_id, "kingdom": ext.get("kingdom", ""),
            "reason": "no_sequence", "num_pairs_lost": usable
        })
        continue

    # Find matching substrates with SMILES
    subs = ext.get("substrates_specific", [])
    matched_subs = []
    for s in subs:
        if s in smiles_map:
            matched_subs.append((s, smiles_map[s]))
        else:
            # Try splitting
            for part in re.split(r'\s+and\s+|\|', s):
                part = part.strip()
                if part in smiles_map:
                    matched_subs.append((part, smiles_map[part]))

    if not matched_subs:
        unresolved_rows.append({
            "cyp_id": cyp_id, "kingdom": ext.get("kingdom", ""),
            "reason": "no_smiles", "num_pairs_lost": len(subs)
        })
        continue

    # Register enzyme
    species = fasta.get("species", "")
    enzyme_key = (cyp_id, species)
    if enzyme_key not in enzyme_key_to_id:
        eid = len(enzyme_rows)
        enzyme_key_to_id[enzyme_key] = eid
        enzyme_rows.append({
            "enzyme_id": eid,
            "uniprot_id": fasta.get("uniprot", ""),
            "p450_symbol": cyp_id,
            "species": species,
            "species_class": json_entries.get(cyp_id, {}).get("Kindom", ""),
            "ec_number": "",
            "sequence": fasta["sequence"],
            "sequence_length": fasta.get("seq_len", len(fasta["sequence"]))
        })

    eid = enzyme_key_to_id[enzyme_key]

    # Register compounds and interactions
    for sub_name, sub_data in matched_subs:
        smiles = sub_data["smiles"]
        if smiles not in smiles_to_cid:
            cid = len(compound_rows)
            smiles_to_cid[smiles] = cid
            compound_rows.append({
                "compound_id": cid,
                "smiles": smiles,
                "name": sub_name
            })
        cid = smiles_to_cid[smiles]

        # Deduplicate interactions
        existing = any(r["enzyme_id"] == eid and r["compound_id"] == cid for r in interaction_rows)
        if existing:
            continue

        interaction_rows.append({
            "interaction_id": len(interaction_rows),
            "enzyme_id": eid,
            "compound_id": cid,
            "label": 1,
            "source": "S9_PCPD",
            "quality_tier": "B",
            "num_reactions": 1,
            "has_pmid": 0,
            "has_products": 0
        })

# Write files
with open(os.path.join(OUT_DIR, "enzymes.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["enzyme_id", "uniprot_id", "p450_symbol", "species",
                                       "species_class", "ec_number", "sequence", "sequence_length"])
    w.writeheader()
    w.writerows(enzyme_rows)

with open(os.path.join(OUT_DIR, "compounds.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["compound_id", "smiles", "name"])
    w.writeheader()
    w.writerows(compound_rows)

with open(os.path.join(OUT_DIR, "interactions.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["interaction_id", "enzyme_id", "compound_id", "label",
                                       "source", "quality_tier", "num_reactions", "has_pmid", "has_products"])
    w.writeheader()
    w.writerows(interaction_rows)

with open(os.path.join(OUT_DIR, "unresolved.csv"), "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["cyp_id", "kingdom", "reason", "num_pairs_lost"])
    w.writeheader()
    w.writerows(unresolved_rows)

print(f"=== S9 PCPD Source Files Generated ===")
print(f"enzymes.csv:      {len(enzyme_rows)} enzymes")
print(f"compounds.csv:    {len(compound_rows)} compounds")
print(f"interactions.csv: {len(interaction_rows)} pairs")
print(f"unresolved.csv:   {len(unresolved_rows)} dropped entries")
print(f"\nOutput: {OUT_DIR}")
