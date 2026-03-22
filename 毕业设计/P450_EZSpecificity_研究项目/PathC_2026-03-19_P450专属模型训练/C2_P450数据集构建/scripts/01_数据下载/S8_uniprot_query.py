"""
Query UniProt for S8 PlantP450DB enzymes.
For each CYP name + species, search UniProt REST API for protein ID + sequence.
Results written incrementally to uniprot_results.jsonl.
"""
import json, re, time, urllib.request, urllib.parse, os, sys

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PlantP450DB")
RESULTS_FILE = os.path.join(BASE, "final_merged_results.jsonl")
SMILES_FILE = os.path.join(BASE, "smiles_cache.jsonl")
OUTPUT_FILE = os.path.join(BASE, "uniprot_results.jsonl")

# Load SMILES set
smiles_set = set()
with open(SMILES_FILE, encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        if e.get("smiles"):
            smiles_set.add(e["name"])

# Get unique CYP+species with usable substrates
cyp_species = {}
with open(RESULTS_FILE, encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        subs = e.get("substrates", [])
        has_usable = any(
            (s if isinstance(s, str) else s.get("name", "")) in smiles_set
            for s in subs
        )
        if has_usable:
            key = (e.get("cyp_name", ""), e.get("species", ""))
            if key not in cyp_species:
                cyp_species[key] = e

print(f"Total CYP+species to query: {len(cyp_species)}")

# Load existing results to support resume
existing = set()
if os.path.exists(OUTPUT_FILE):
    with open(OUTPUT_FILE, encoding="utf-8") as f:
        for line in f:
            e = json.loads(line)
            existing.add((e.get("cyp_name", ""), e.get("species", "")))
    print(f"Already queried: {len(existing)}, resuming...")

# Extract CYP number for query (e.g., "AtCYP51G1" -> "CYP51G1")
def extract_cyp_id(name):
    """Extract the CYP identifier from various naming formats."""
    # Standard CYP format: CYP + family + subfamily + number
    m = re.search(r'(CYP\d+[A-Z]+\d*)', name, re.IGNORECASE)
    if m:
        return m.group(1)
    return name

def query_uniprot(cyp_name, species):
    """Query UniProt REST API. Returns (uniprot_id, protein_name, sequence, gene_name) or None."""
    cyp_id = extract_cyp_id(cyp_name)

    # Try multiple query strategies
    queries = []

    # Strategy 1: gene name + organism (most specific)
    queries.append(f'(gene:{cyp_id}) AND (organism_name:"{species}")')

    # Strategy 2: protein name contains CYP ID + organism
    queries.append(f'(protein_name:{cyp_id}) AND (organism_name:"{species}")')

    # Strategy 3: just the CYP ID + organism (broader)
    queries.append(f'{cyp_id} AND (organism_name:"{species}")')

    # Strategy 4: original name + organism
    if cyp_id != cyp_name:
        queries.append(f'(gene:{cyp_name}) AND (organism_name:"{species}")')

    for q in queries:
        encoded = urllib.parse.quote(q)
        url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded}&format=json&size=5&fields=accession,protein_name,gene_names,sequence,organism_name"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PlantP450DB-UniProt/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode())
                results = data.get("results", [])
                if results:
                    # Pick best match (prefer reviewed)
                    best = results[0]
                    uniprot_id = best.get("primaryAccession", "")
                    protein_name = ""
                    pn = best.get("proteinDescription", {})
                    if pn.get("recommendedName"):
                        protein_name = pn["recommendedName"].get("fullName", {}).get("value", "")
                    elif pn.get("submissionNames"):
                        protein_name = pn["submissionNames"][0].get("fullName", {}).get("value", "")
                    gene_name = ""
                    genes = best.get("genes", [])
                    if genes:
                        gene_name = genes[0].get("geneName", {}).get("value", "") or genes[0].get("orderedLocusNames", [{}])[0].get("value", "") if genes[0].get("orderedLocusNames") else ""
                    sequence = best.get("sequence", {}).get("value", "")
                    seq_len = best.get("sequence", {}).get("length", 0)
                    return uniprot_id, protein_name, gene_name, sequence, seq_len
        except Exception:
            pass
        time.sleep(0.3)

    return None

# Query UniProt
out_f = open(OUTPUT_FILE, "a", encoding="utf-8")
total = len(cyp_species)
queried = 0
found = 0
skipped = 0

for (cyp_name, species) in sorted(cyp_species.keys()):
    if (cyp_name, species) in existing:
        skipped += 1
        continue

    result = query_uniprot(cyp_name, species)
    queried += 1

    if result:
        uniprot_id, protein_name, gene_name, sequence, seq_len = result
        entry = {
            "cyp_name": cyp_name,
            "species": species,
            "uniprot_id": uniprot_id,
            "protein_name": protein_name,
            "gene_name": gene_name,
            "sequence": sequence,
            "sequence_length": seq_len,
            "error": ""
        }
        found += 1
        status = "FOUND"
    else:
        entry = {
            "cyp_name": cyp_name,
            "species": species,
            "uniprot_id": "",
            "protein_name": "",
            "gene_name": "",
            "sequence": "",
            "sequence_length": 0,
            "error": "not_found"
        }
        status = "MISS"

    out_f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    out_f.flush()

    if queried % 20 == 0 or status == "MISS":
        print(f"  [{skipped + queried}/{total}] {status}: {cyp_name} ({species})" +
              (f" -> {entry['uniprot_id']}" if status == "FOUND" else ""))

    time.sleep(0.5)  # Rate limit

out_f.close()

print(f"\n=== DONE ===")
print(f"Total: {total}")
print(f"Skipped (already done): {skipped}")
print(f"Queried this run: {queried}")
print(f"Found: {found}")
print(f"Not found: {queried - found}")
print(f"Success rate: {found}/{queried} = {found/max(queried,1)*100:.1f}%")
