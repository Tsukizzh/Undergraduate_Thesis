"""
S8 Track 1 (fixed): PubChem SMILES lookup with correct field names + name cleaning
"""
import json, os, urllib.request, urllib.parse, time, re

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DL_DIR = os.path.join(BASE, 'downloads', 'PlantP450DB')
CACHE_FILE = os.path.join(DL_DIR, 'track1_output', 'smiles_cache.jsonl')

# Load found results
with open(os.path.join(DL_DIR, 'all_results_merged.jsonl'), 'r', encoding='utf-8') as f:
    all_results = [json.loads(l.strip()) for l in f if l.strip()]
found = [r for r in all_results if r['status'] == 'found']
all_substrates = set()
for r in found:
    for s in r.get('substrates', []):
        all_substrates.add(s)

# Load existing cache
cache = {}
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                rec = json.loads(line.strip())
                cache[rec['name']] = rec

remaining = sorted(s for s in all_substrates if s not in cache)
print(f"Total substrates: {len(all_substrates)}, Cached: {len(cache)}, Remaining: {len(remaining)}")

def clean_name(name):
    cleaned = re.sub(r'^\(\+\)-', '', name)
    cleaned = re.sub(r'^\(-\)-', '', cleaned)
    cleaned = re.sub(r'^\(\+/-\)-', '', cleaned)
    cleaned = re.sub(r'\s*\(tested\)', '', cleaned)
    cleaned = re.sub(r'\s*\(secondary\)', '', cleaned)
    cleaned = re.sub(r'\s*\(primary\)', '', cleaned)
    cleaned = re.sub(r'\s*\(preferred\)', '', cleaned)
    cleaned = re.sub(r'\s*\(best\)', '', cleaned)
    cleaned = re.sub(r'\s*\(ABA\)', '', cleaned)
    cleaned = re.sub(r'\s*\(major\)', '', cleaned)
    cleaned = re.sub(r'\s*\(minor\)', '', cleaned)
    return cleaned.strip()

def query_pubchem_name(name):
    encoded = urllib.parse.quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/property/CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    req = urllib.request.Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
    resp = urllib.request.urlopen(req, timeout=10)
    data = json.loads(resp.read().decode('utf-8'))
    props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
    smiles = props.get('SMILES', '') or props.get('IsomericSMILES', '') or props.get('CanonicalSMILES', '')
    inchikey = props.get('InChIKey', '')
    cid = str(props.get('CID', ''))
    return smiles, inchikey, cid

cache_f = open(CACHE_FILE, 'a', encoding='utf-8')
success = 0
fail = 0

for i, name in enumerate(remaining):
    rec = {'name': name, 'smiles': '', 'inchikey': '', 'cid': '', 'error': ''}

    # Strategy 1: Try original name
    try:
        smiles, inchikey, cid = query_pubchem_name(name)
        if smiles:
            rec = {'name': name, 'smiles': smiles, 'inchikey': inchikey, 'cid': cid, 'error': ''}
            success += 1
            cache[name] = rec
            cache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
            cache_f.flush()
            if (i+1) % 50 == 0:
                print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
            time.sleep(0.2)
            continue
    except:
        pass

    # Strategy 2: Clean name
    cleaned = clean_name(name)
    if cleaned and cleaned != name:
        try:
            smiles, inchikey, cid = query_pubchem_name(cleaned)
            if smiles:
                rec = {'name': name, 'smiles': smiles, 'inchikey': inchikey, 'cid': cid, 'error': ''}
                success += 1
                cache[name] = rec
                cache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
                cache_f.flush()
                if (i+1) % 50 == 0:
                    print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
                time.sleep(0.2)
                continue
        except:
            pass

    # Strategy 3: Remove all special chars
    simple = re.sub(r'[^a-zA-Z0-9\s-]', '', name).strip()
    if simple and simple != name and simple != cleaned:
        try:
            smiles, inchikey, cid = query_pubchem_name(simple)
            if smiles:
                rec = {'name': name, 'smiles': smiles, 'inchikey': inchikey, 'cid': cid, 'error': ''}
                success += 1
                cache[name] = rec
                cache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
                cache_f.flush()
                if (i+1) % 50 == 0:
                    print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
                time.sleep(0.2)
                continue
        except:
            pass

    rec['error'] = 'not_found_after_3_strategies'
    fail += 1
    cache[name] = rec
    cache_f.write(json.dumps(rec, ensure_ascii=False) + '\n')
    cache_f.flush()

    if (i+1) % 50 == 0:
        print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
    time.sleep(0.2)

cache_f.close()

ok = sum(1 for v in cache.values() if v.get('smiles'))
print(f"\nDone. New success: {success}, New fail: {fail}")
print(f"Total SMILES resolved: {ok}/{len(cache)} ({ok/len(cache)*100:.1f}%)")
print(f"Failed names sample: {[v['name'] for v in cache.values() if not v.get('smiles')][:10]}")
