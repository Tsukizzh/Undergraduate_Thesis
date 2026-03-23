"""
Query PubChem for SMILES of all substrate names from PCPD extraction.
Writes results incrementally to smiles_cache.jsonl. Supports resume.
"""
import json, re, time, urllib.request, urllib.parse, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PCPD")
INPUT = os.path.join(BASE, "substrate_extraction_cleaned.jsonl")
OUTPUT = os.path.join(BASE, "smiles_cache.jsonl")

# Collect all unique substrate names
all_names = set()
with open(INPUT, encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        all_names.update(e.get("substrates_specific", []))

print(f"Unique substrate names to query: {len(all_names)}")

# Load existing cache
existing = {}
if os.path.exists(OUTPUT):
    with open(OUTPUT, encoding="utf-8") as f:
        for line in f:
            e = json.loads(line)
            existing[e["name"]] = e
    print(f"Already cached: {len(existing)}")

to_query = [n for n in sorted(all_names) if n not in existing]
print(f"Need to query: {len(to_query)}")

def clean_for_pubchem(name):
    """Clean substrate name for PubChem query. Returns list of query strings."""
    queries = []
    clean = re.sub(r'\s*\([^()]*\)\s*$', '', name).strip()
    clean = re.sub(r'\s*\([^()]*\)\s*$', '', clean).strip()
    if clean:
        queries.append(clean)
    parens = re.findall(r'\(([^()]+)\)', name)
    for p in parens:
        p = p.strip()
        if len(p) > 3 and not any(kw in p.lower() for kw in [
            'unspecified', 'broad', 'implied', 'trace', 'minor', 'poor',
            'tested', 'primary', 'pathway', 'derived', 'specificity',
            'e.g.', 'likely', 'class'
        ]):
            queries.append(p)
    fallback = re.sub(r'^\(\+\)-|^\(-\)-|^\(S\)-|^\(R\)-|^\(E\)-|^\(Z\)-', '', clean).strip()
    if fallback and fallback != clean:
        queries.append(fallback)
    return queries

def query_pubchem(name_queries):
    for q in name_queries:
        encoded = urllib.parse.quote(q)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/IsomericSMILES,InChIKey/JSON"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PCPD-SMILES/1.0"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
                props = data["PropertyTable"]["Properties"][0]
                smiles = props.get("IsomericSMILES", "") or props.get("SMILES", "")
                inchikey = props.get("InChIKey", "")
                cid = str(props.get("CID", ""))
                if smiles:
                    return smiles, inchikey, cid, ""
        except urllib.error.HTTPError as e:
            if e.code == 404:
                continue
            elif e.code == 503:
                time.sleep(3)
                continue
        except Exception:
            continue
        time.sleep(0.2)
    return "", "", "", f"not_found_after_{len(name_queries)}_queries"

out_f = open(OUTPUT, "a", encoding="utf-8")
queried = 0
found = 0

for name in to_query:
    queries = clean_for_pubchem(name)
    smiles, inchikey, cid, error = query_pubchem(queries)
    entry = {"name": name, "smiles": smiles, "inchikey": inchikey, "cid": cid, "error": error}
    out_f.write(json.dumps(entry, ensure_ascii=False) + "\n")
    out_f.flush()
    queried += 1
    if smiles:
        found += 1
    if queried % 50 == 0:
        print(f"  [{queried}/{len(to_query)}] found={found}")
    time.sleep(0.3)

out_f.close()

print(f"\n=== DONE ===")
print(f"Queried: {queried}")
print(f"Found SMILES: {found}")
print(f"Not found: {queried - found}")
print(f"Success rate: {found}/{queried} = {found/max(queried,1)*100:.1f}%")

# Final stats including existing cache
total_cached = len(existing) + queried
total_with_smiles = sum(1 for e in existing.values() if e.get("smiles")) + found
print(f"\nTotal cached: {total_cached}")
print(f"Total with SMILES: {total_with_smiles}")
