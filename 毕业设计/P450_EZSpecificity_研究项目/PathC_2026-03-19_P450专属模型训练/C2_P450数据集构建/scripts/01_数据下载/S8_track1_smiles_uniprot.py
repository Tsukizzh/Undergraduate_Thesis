"""
S8 Track 1: 将 found 条目的底物名称转为 SMILES，获取 UniProt ID + 序列
逐条写入，支持断点续传
"""
import json, os, time, urllib.request, urllib.parse, csv

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DL_DIR = os.path.join(BASE, 'downloads', 'PlantP450DB')
OUT_DIR = os.path.join(DL_DIR, 'track1_output')
os.makedirs(OUT_DIR, exist_ok=True)

SMILES_CACHE_FILE = os.path.join(OUT_DIR, 'smiles_cache.jsonl')
UNIPROT_CACHE_FILE = os.path.join(OUT_DIR, 'uniprot_cache.jsonl')
RESULTS_FILE = os.path.join(OUT_DIR, 'track1_results.jsonl')

# ============================================================
# Step 1: Load found entries
# ============================================================
print("=" * 60)
print("S8 Track 1: SMILES + UniProt")
print("=" * 60)

with open(os.path.join(DL_DIR, 'all_results_merged.jsonl'), 'r', encoding='utf-8') as f:
    all_results = [json.loads(l.strip()) for l in f if l.strip()]

found = [r for r in all_results if r['status'] == 'found']
print(f"Found entries: {len(found)}")

# Load combined input for species/accession info
with open(os.path.join(DL_DIR, 'combined_input.json'), 'r', encoding='utf-8') as f:
    combined = json.load(f)
entry_map = {e['entry_id']: e for e in combined}

# ============================================================
# Step 2: Collect unique substrate names and query PubChem
# ============================================================
all_substrates = set()
for r in found:
    for s in r.get('substrates', []):
        all_substrates.add(s)

print(f"\nUnique substrate names: {len(all_substrates)}")

# Load existing cache
smiles_cache = {}
if os.path.exists(SMILES_CACHE_FILE):
    with open(SMILES_CACHE_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rec = json.loads(line.strip())
                smiles_cache[rec['name']] = rec
    print(f"SMILES cache loaded: {len(smiles_cache)}")

remaining_subs = [s for s in sorted(all_substrates) if s not in smiles_cache]
print(f"Remaining to query: {len(remaining_subs)}")

cache_f = open(SMILES_CACHE_FILE, 'a', encoding='utf-8')
success = 0
fail = 0

for i, name in enumerate(remaining_subs):
    rec = {'name': name, 'smiles': '', 'inchikey': '', 'cid': '', 'error': ''}
    try:
        encoded = urllib.parse.quote(name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
        req = urllib.request.Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = json.loads(resp.read().decode('utf-8'))
        props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
        rec['smiles'] = props.get('IsomericSMILES', '') or props.get('CanonicalSMILES', '')
        rec['inchikey'] = props.get('InChIKey', '')
        rec['cid'] = str(props.get('CID', ''))
        if rec['smiles']:
            success += 1
        else:
            rec['error'] = 'no_smiles_in_response'
            fail += 1
    except Exception as e:
        rec['error'] = str(e)[:80]
        fail += 1

    smiles_cache[name] = rec
    cache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
    cache_f.flush()

    if (i + 1) % 50 == 0:
        print(f"  SMILES: {i+1}/{len(remaining_subs)}, success={success}, fail={fail}")
    time.sleep(0.25)

cache_f.close()
print(f"SMILES done. Success: {success}, Fail: {fail}")

# ============================================================
# Step 3: Get UniProt ID + sequence for each CYP
# ============================================================
print(f"\n[Step 3] UniProt lookup...")

# Collect unique (cyp_name, species) pairs from ALL entries (not just found)
cyp_species_pairs = set()
for e in combined:
    cyp_species_pairs.add((e['cyp_name'], e['species']))

print(f"Unique (CYP, species) pairs: {len(cyp_species_pairs)}")

# Load cache
uniprot_cache = {}
if os.path.exists(UNIPROT_CACHE_FILE):
    with open(UNIPROT_CACHE_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rec = json.loads(line.strip())
                uniprot_cache[(rec['cyp_name'], rec['species'])] = rec
    print(f"UniProt cache loaded: {len(uniprot_cache)}")

remaining_cyps = [(c, s) for c, s in sorted(cyp_species_pairs) if (c, s) not in uniprot_cache]
print(f"Remaining to query: {len(remaining_cyps)}")

ucache_f = open(UNIPROT_CACHE_FILE, 'a', encoding='utf-8')
u_success = 0
u_fail = 0

for i, (cyp, species) in enumerate(remaining_cyps):
    rec = {'cyp_name': cyp, 'species': species, 'uniprot_id': '', 'sequence': '', 'seq_length': 0, 'error': ''}

    # Try gene name search
    query = f'gene_exact:{cyp}'
    if species:
        query += f' AND organism_name:"{species}"'
    try:
        search_url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&format=json&size=3&fields=accession,gene_names,organism_name,sequence"
        req = urllib.request.Request(search_url, headers={'User-Agent': 'Python/EZSpecificity'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = json.loads(resp.read().decode('utf-8'))
        results = data.get('results', [])

        if len(results) == 1:
            rec['uniprot_id'] = results[0]['primaryAccession']
            rec['sequence'] = results[0].get('sequence', {}).get('value', '')
            rec['seq_length'] = len(rec['sequence'])
            u_success += 1
        elif len(results) > 1:
            rec['uniprot_id'] = results[0]['primaryAccession']
            rec['sequence'] = results[0].get('sequence', {}).get('value', '')
            rec['seq_length'] = len(rec['sequence'])
            rec['error'] = f'multiple_hits:{len(results)}'
            u_success += 1
        else:
            # Try protein name search
            query2 = f'protein_name:{cyp}'
            if species:
                query2 += f' AND organism_name:"{species}"'
            search_url2 = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query2)}&format=json&size=3&fields=accession,gene_names,organism_name,sequence"
            req2 = urllib.request.Request(search_url2, headers={'User-Agent': 'Python/EZSpecificity'})
            resp2 = urllib.request.urlopen(req2, timeout=10)
            data2 = json.loads(resp2.read().decode('utf-8'))
            results2 = data2.get('results', [])
            if results2:
                rec['uniprot_id'] = results2[0]['primaryAccession']
                rec['sequence'] = results2[0].get('sequence', {}).get('value', '')
                rec['seq_length'] = len(rec['sequence'])
                rec['error'] = 'matched_by_protein_name'
                u_success += 1
            else:
                rec['error'] = 'not_found'
                u_fail += 1
    except Exception as e:
        rec['error'] = str(e)[:80]
        u_fail += 1

    uniprot_cache[(cyp, species)] = rec
    ucache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
    ucache_f.flush()

    if (i + 1) % 50 == 0:
        print(f"  UniProt: {i+1}/{len(remaining_cyps)}, success={u_success}, fail={u_fail}")
    time.sleep(0.4)

ucache_f.close()
print(f"UniProt done. Success: {u_success}, Fail: {u_fail}")

# ============================================================
# Step 4: Summary
# ============================================================
smiles_ok = sum(1 for v in smiles_cache.values() if v.get('smiles'))
smiles_fail = sum(1 for v in smiles_cache.values() if not v.get('smiles'))
uniprot_ok = sum(1 for v in uniprot_cache.values() if v.get('uniprot_id'))
uniprot_fail = sum(1 for v in uniprot_cache.values() if not v.get('uniprot_id'))

print(f"\n{'=' * 60}")
print(f"Track 1 Summary")
print(f"{'=' * 60}")
print(f"SMILES: {smiles_ok}/{len(smiles_cache)} resolved ({smiles_ok/max(len(smiles_cache),1)*100:.1f}%)")
print(f"UniProt: {uniprot_ok}/{len(uniprot_cache)} resolved ({uniprot_ok/max(len(uniprot_cache),1)*100:.1f}%)")
