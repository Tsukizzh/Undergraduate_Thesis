"""
Retry UniProt query for 370 misses using GenBank accession numbers from all_entries.json.
"""
import json, re, time, urllib.request, urllib.parse, os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "downloads", "PlantP450DB")
RESULTS_FILE = os.path.join(BASE, "uniprot_results.jsonl")
ALL_ENTRIES = os.path.join(BASE, "all_entries.json")

# Load accession numbers from all_entries
acc_map = {}
with open(ALL_ENTRIES, encoding="utf-8") as f:
    data = json.load(f)
    for e in data:
        name = e.get("Name", "")
        species = e.get("Species", "")
        acc = e.get("Accession number", "").strip()
        if acc:
            acc_map[(name, species)] = acc

print(f"Accession numbers available: {len(acc_map)}")

# Load current results, find misses
entries = []
misses = []
with open(RESULTS_FILE, encoding="utf-8") as f:
    for line in f:
        e = json.loads(line)
        entries.append(e)
        if not e.get("sequence"):
            misses.append(e)

print(f"Total misses to retry: {len(misses)}")

# Count misses with accession
misses_with_acc = [(e, acc_map.get((e["cyp_name"], e["species"]), "")) for e in misses]
has_acc = sum(1 for _, a in misses_with_acc if a)
print(f"Misses with GenBank accession: {has_acc}")

def query_uniprot_by_xref(accession):
    """Query UniProt by EMBL/GenBank cross-reference."""
    q = urllib.parse.quote(f'(xref:embl-{accession})')
    url = f"https://rest.uniprot.org/uniprotkb/search?query={q}&format=json&size=3&fields=accession,protein_name,gene_names,sequence,organism_name"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "PlantP450DB-UniProt-R2/1.0"})
        with urllib.request.urlopen(req, timeout=20) as resp:
            data = json.loads(resp.read().decode())
            results = data.get("results", [])
            if results:
                best = results[0]
                uid = best.get("primaryAccession", "")
                seq = best.get("sequence", {}).get("value", "")
                seq_len = best.get("sequence", {}).get("length", 0)
                pn = best.get("proteinDescription", {})
                protein_name = ""
                if pn.get("recommendedName"):
                    protein_name = pn["recommendedName"].get("fullName", {}).get("value", "")
                elif pn.get("submissionNames"):
                    protein_name = pn["submissionNames"][0].get("fullName", {}).get("value", "")
                gene_name = ""
                genes = best.get("genes", [])
                if genes and genes[0].get("geneName"):
                    gene_name = genes[0]["geneName"].get("value", "")
                return uid, protein_name, gene_name, seq, seq_len
    except:
        pass
    return None

def query_uniprot_by_gene_organism(cyp_name, species):
    """Additional query strategies for plant P450s."""
    # Try with species genus only (e.g., "Arabidopsis" instead of "Arabidopsis thaliana")
    genus = species.split()[0] if species else ""

    cyp_id = re.search(r'(CYP\d+[A-Z]+\d*)', cyp_name, re.IGNORECASE)
    cyp_id = cyp_id.group(1) if cyp_id else cyp_name

    queries = []
    if genus:
        queries.append(f'(gene:{cyp_id}) AND (organism_name:"{genus}")')
        # Try without prefix (e.g., "CYP51G1" from "AtCYP51G1")
        bare = re.sub(r'^[A-Za-z]{2,4}(?=CYP)', '', cyp_name)
        if bare != cyp_name:
            queries.append(f'(gene:{bare}) AND (organism_name:"{genus}")')

    for q in queries:
        encoded = urllib.parse.quote(q)
        url = f"https://rest.uniprot.org/uniprotkb/search?query={encoded}&format=json&size=3&fields=accession,protein_name,gene_names,sequence,organism_name"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "PlantP450DB-UniProt-R2/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read().decode())
                results = data.get("results", [])
                if results:
                    best = results[0]
                    uid = best.get("primaryAccession", "")
                    seq = best.get("sequence", {}).get("value", "")
                    seq_len = best.get("sequence", {}).get("length", 0)
                    pn = best.get("proteinDescription", {})
                    protein_name = ""
                    if pn.get("recommendedName"):
                        protein_name = pn["recommendedName"].get("fullName", {}).get("value", "")
                    elif pn.get("submissionNames"):
                        protein_name = pn["submissionNames"][0].get("fullName", {}).get("value", "")
                    gene_name = ""
                    genes = best.get("genes", [])
                    if genes and genes[0].get("geneName"):
                        gene_name = genes[0]["geneName"].get("value", "")
                    return uid, protein_name, gene_name, seq, seq_len
        except:
            pass
        time.sleep(0.3)
    return None

# Retry
found_xref = 0
found_gene = 0
updated = {}

for e, acc in misses_with_acc:
    key = (e["cyp_name"], e["species"])

    # Strategy 1: GenBank accession cross-reference
    if acc:
        result = query_uniprot_by_xref(acc)
        time.sleep(0.4)
        if result:
            uid, pname, gname, seq, slen = result
            updated[key] = {
                "cyp_name": e["cyp_name"], "species": e["species"],
                "uniprot_id": uid, "protein_name": pname, "gene_name": gname,
                "sequence": seq, "sequence_length": slen, "error": ""
            }
            found_xref += 1
            if (found_xref + found_gene) % 20 == 0:
                print(f"  [{found_xref+found_gene}] XREF: {e['cyp_name']} -> {uid}")
            continue

    # Strategy 2: genus-level gene search
    result = query_uniprot_by_gene_organism(e["cyp_name"], e["species"])
    time.sleep(0.4)
    if result:
        uid, pname, gname, seq, slen = result
        updated[key] = {
            "cyp_name": e["cyp_name"], "species": e["species"],
            "uniprot_id": uid, "protein_name": pname, "gene_name": gname,
            "sequence": seq, "sequence_length": slen, "error": ""
        }
        found_gene += 1
        if (found_xref + found_gene) % 20 == 0:
            print(f"  [{found_xref+found_gene}] GENE: {e['cyp_name']} -> {uid}")

print(f"\nRetry results: xref={found_xref}, gene={found_gene}, total new={found_xref+found_gene}")

# Rewrite results file with updates
with open(RESULTS_FILE, "w", encoding="utf-8") as f:
    for e in entries:
        key = (e["cyp_name"], e["species"])
        if key in updated:
            f.write(json.dumps(updated[key], ensure_ascii=False) + "\n")
        else:
            f.write(json.dumps(e, ensure_ascii=False) + "\n")

# Final stats
total_with_seq = sum(1 for e in entries if e.get("sequence") or (e["cyp_name"], e["species"]) in updated)
print(f"\nFinal: {total_with_seq}/{len(entries)} have sequences ({total_with_seq/len(entries)*100:.1f}%)")
print(f"Still missing: {len(entries) - total_with_seq}")
