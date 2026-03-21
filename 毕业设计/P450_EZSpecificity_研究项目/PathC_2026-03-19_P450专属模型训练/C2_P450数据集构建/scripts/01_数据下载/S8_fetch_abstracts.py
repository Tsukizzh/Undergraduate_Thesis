"""
S8 Plant P450 DB: 批量获取论文摘要
逐条写入，支持断点续传
"""
import json, os, time, urllib.request

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DL_DIR = os.path.join(BASE, 'downloads', 'PlantP450DB')
ENTRIES_FILE = os.path.join(DL_DIR, 'all_entries.json')
ABSTRACTS_FILE = os.path.join(DL_DIR, 'abstracts.jsonl')  # JSONL: one JSON per line

# Load entries
with open(ENTRIES_FILE, 'r', encoding='utf-8') as f:
    entries = json.load(f)

# Collect unique DOIs
all_dois = set()
for e in entries:
    dois = e.get('DOIs', '')
    if dois:
        for doi in dois.split('; '):
            doi = doi.strip()
            if doi:
                all_dois.add(doi)

doi_list = sorted(all_dois)
print(f'Unique DOIs: {len(doi_list)}')

# Load already fetched (resume support)
done_dois = set()
if os.path.exists(ABSTRACTS_FILE):
    with open(ABSTRACTS_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    rec = json.loads(line)
                    done_dois.add(rec['doi'])
                except:
                    pass
    print(f'Already fetched: {len(done_dois)}')

remaining = [d for d in doi_list if d not in done_dois]
print(f'Remaining: {len(remaining)}')

# Fetch abstracts
out_f = open(ABSTRACTS_FILE, 'a', encoding='utf-8')
success = 0
fail = 0

for i, doi in enumerate(remaining):
    doi_id = doi.replace('https://doi.org/', '').replace('http://doi.org/', '')
    record = {'doi': doi, 'pmid': '', 'title': '', 'abstract': '', 'error': ''}

    try:
        # Step 1: DOI -> PMID
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={doi_id}[doi]&retmode=json'
        req = urllib.request.Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
        resp = urllib.request.urlopen(req, timeout=10)
        data = json.loads(resp.read().decode('utf-8'))
        ids = data.get('esearchresult', {}).get('idlist', [])

        if not ids:
            record['error'] = 'no_pmid'
            fail += 1
        else:
            pmid = ids[0]
            record['pmid'] = pmid
            time.sleep(0.3)

            # Step 2: PMID -> abstract
            url2 = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id={pmid}&rettype=abstract&retmode=text'
            req2 = urllib.request.Request(url2, headers={'User-Agent': 'Python/EZSpecificity'})
            resp2 = urllib.request.urlopen(req2, timeout=10)
            text = resp2.read().decode('utf-8')

            # Parse title and abstract from the text
            lines = text.strip().split('\n')
            # Title is usually after the citation line
            title_lines = []
            abstract_lines = []
            in_abstract = False
            for line in lines:
                if line.strip() == '' and title_lines and not in_abstract:
                    in_abstract = True
                    continue
                if not in_abstract and len(title_lines) < 3:
                    if not line.startswith('Author') and not line.startswith('(') and not line.strip().startswith('DOI') and len(line.strip()) > 5:
                        if not any(line.strip().startswith(x) for x in ['PMID', 'PMCID', 'Copyright', 'Published']):
                            title_lines.append(line.strip())
                elif in_abstract:
                    abstract_lines.append(line.strip())

            record['title'] = ' '.join(title_lines)
            record['abstract'] = ' '.join(abstract_lines)
            record['full_text'] = text
            success += 1

    except Exception as ex:
        record['error'] = str(ex)[:100]
        fail += 1

    # Write immediately
    out_f.write(json.dumps(record, ensure_ascii=False) + '\n')
    out_f.flush()

    if (i + 1) % 20 == 0:
        print(f'  Progress: {i+1}/{len(remaining)}, success={success}, fail={fail}')

    time.sleep(0.35)

out_f.close()
print(f'\nDone. Success: {success}, Fail: {fail}, Total: {success+fail}')
print(f'Saved to: {ABSTRACTS_FILE}')
