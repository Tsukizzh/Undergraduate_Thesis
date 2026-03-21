"""
S8 Track 2: 下载 OA 全文论文（HTML/text 形式）
从 Unpaywall 找到的 76 篇 OA 论文中获取可读文本
逐条写入，支持断点续传
"""
import json, os, urllib.request, time, re

BASE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DL_DIR = os.path.join(BASE, 'downloads', 'PlantP450DB')
OUT_DIR = os.path.join(DL_DIR, 'track2_fulltext')
os.makedirs(OUT_DIR, exist_ok=True)

FULLTEXT_FILE = os.path.join(OUT_DIR, 'fulltext.jsonl')

# Load Unpaywall results
with open(os.path.join(DL_DIR, 'unpaywall_results.json'), 'r', encoding='utf-8') as f:
    unpaywall = json.load(f)

oa_papers = unpaywall['oa']
print(f"OA papers to fetch: {len(oa_papers)}")

# Load already fetched (resume support)
done_dois = set()
if os.path.exists(FULLTEXT_FILE):
    with open(FULLTEXT_FILE, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                try:
                    rec = json.loads(line.strip())
                    done_dois.add(rec['doi'])
                except:
                    pass
    print(f"Already fetched: {len(done_dois)}")

remaining = [p for p in oa_papers if p['doi'] not in done_dois]
print(f"Remaining: {len(remaining)}")

# Also try PMC for papers with PMIDs
doi_to_pmid = {}
with open(os.path.join(DL_DIR, 'abstracts.jsonl'), 'r', encoding='utf-8') as f:
    for line in f:
        if line.strip():
            rec = json.loads(line.strip())
            if rec.get('pmid'):
                doi_to_pmid[rec['doi']] = rec['pmid']

out_f = open(FULLTEXT_FILE, 'a', encoding='utf-8')
success = 0
fail = 0

for i, paper in enumerate(remaining):
    doi = paper['doi']
    pdf_url = paper.get('pdf_url', '')
    host = paper.get('host', '')
    record = {'doi': doi, 'text': '', 'source': '', 'error': ''}

    # Strategy 1: Try PMC full text (best quality - structured text)
    pmid = doi_to_pmid.get(doi, '')
    if pmid:
        try:
            # Get PMCID from PMID
            url = f"https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids={pmid}&format=json&tool=ezspec&email=test@test.com"
            req = urllib.request.Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
            resp = urllib.request.urlopen(req, timeout=10)
            data = json.loads(resp.read().decode('utf-8'))
            records = data.get('records', [])
            pmcid = records[0].get('pmcid', '') if records else ''

            if pmcid:
                # Fetch full text from PMC
                pmc_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmcid}&rettype=text&retmode=text"
                req2 = urllib.request.Request(pmc_url, headers={'User-Agent': 'Python/EZSpecificity'})
                resp2 = urllib.request.urlopen(req2, timeout=30)
                text = resp2.read().decode('utf-8', errors='replace')
                if len(text) > 500:
                    record['text'] = text
                    record['source'] = f'PMC:{pmcid}'
                    success += 1
                    out_f.write(json.dumps(record, ensure_ascii=False) + '\n')
                    out_f.flush()
                    if (i+1) % 10 == 0:
                        print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
                    time.sleep(0.5)
                    continue
        except:
            pass

    # Strategy 2: Try Unpaywall URL (HTML page -> extract text)
    if pdf_url and not pdf_url.endswith('.pdf'):
        try:
            req3 = urllib.request.Request(pdf_url, headers={'User-Agent': 'Mozilla/5.0'})
            resp3 = urllib.request.urlopen(req3, timeout=20)
            html = resp3.read().decode('utf-8', errors='replace')
            # Basic HTML to text
            text = re.sub(r'<script[^>]*>.*?</script>', '', html, flags=re.DOTALL)
            text = re.sub(r'<style[^>]*>.*?</style>', '', text, flags=re.DOTALL)
            text = re.sub(r'<[^>]+>', ' ', text)
            text = re.sub(r'\s+', ' ', text).strip()
            if len(text) > 1000:
                record['text'] = text[:50000]  # Cap at 50K chars
                record['source'] = f'unpaywall_html:{pdf_url[:80]}'
                success += 1
                out_f.write(json.dumps(record, ensure_ascii=False) + '\n')
                out_f.flush()
                if (i+1) % 10 == 0:
                    print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
                time.sleep(1)
                continue
        except:
            pass

    # Strategy 3: For PDF URLs, we can't extract text easily - skip
    record['error'] = f'pdf_only_or_failed:{pdf_url[:80] if pdf_url else "no_url"}'
    fail += 1
    out_f.write(json.dumps(record, ensure_ascii=False) + '\n')
    out_f.flush()

    if (i+1) % 10 == 0:
        print(f"  {i+1}/{len(remaining)}: success={success}, fail={fail}")
    time.sleep(0.5)

out_f.close()
print(f"\nDone. Success: {success}, Fail: {fail}")
