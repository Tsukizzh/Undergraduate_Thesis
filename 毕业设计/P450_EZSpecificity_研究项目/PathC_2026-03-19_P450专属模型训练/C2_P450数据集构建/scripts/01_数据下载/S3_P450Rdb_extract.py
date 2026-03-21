"""
S3_P450Rdb: 完整提取脚本
修复所有已知问题后提取到统一schema
"""
import csv, sys, os, json, hashlib
from urllib.request import urlopen, Request
from urllib.parse import quote
from urllib.error import HTTPError
import time

csv.field_size_limit(10**7)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))  # scripts/01_数据下载/
C2_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))    # C2_P450数据集构建/
RAW_DIR = os.path.join(C2_DIR, 'downloads', 'P450Rdb')
OUT_DIR = os.path.join(C2_DIR, 'data', 'sources', 'Source_P450Rdb')
os.makedirs(OUT_DIR, exist_ok=True)

print("=" * 60)
print("S3_P450Rdb 完整提取")
print("=" * 60)

# ============================================================
# Step 1: Load validated reactions
# ============================================================
print("\n[Step 1] 加载验证反应数据...")
with open(os.path.join(RAW_DIR, 'Reactions.csv'), 'r', encoding='utf-8', errors='replace') as f:
    reader = csv.reader(f)
    header = next(reader)
    all_rows = list(reader)

validated = [r for r in all_rows if len(r) >= 85 and r[84] == 'Experimentally validated']
print(f"  验证反应: {len(validated)} 行")

# ============================================================
# Step 2: Load FASTA sequences
# ============================================================
print("\n[Step 2] 加载FASTA序列...")
sequences = {}
current_id = None
current_seq = []
with open(os.path.join(RAW_DIR, 'Sequences.fasta'), 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            if current_id:
                sequences[current_id] = ''.join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        sequences[current_id] = ''.join(current_seq)
print(f"  FASTA序列: {len(sequences)} 条")

# ============================================================
# Step 3 (Issue 3): Split compound UniProt IDs
# ============================================================
print("\n[Step 3] 拆分复合UniProt ID...")
expanded_rows = []
split_count = 0
for r in validated:
    uid = r[4]
    if not uid or uid == '/':
        expanded_rows.append(r)
        continue

    parts = []
    if ';' in uid:
        parts = [p.strip() for p in uid.split(';') if p.strip() and p.strip() != '/']
    elif '/' in uid:
        parts = [p.strip() for p in uid.split('/') if p.strip()]

    if len(parts) > 1:
        split_count += 1
        for p in parts:
            new_row = list(r)
            new_row[4] = p
            expanded_rows.append(new_row)
    else:
        expanded_rows.append(r)

print(f"  拆分了 {split_count} 条复合ID, 新增 {len(expanded_rows) - len(validated)} 行")
validated = expanded_rows

# ============================================================
# Step 4 (Issue 4): Fetch missing sequences from UniProt API
# ============================================================
print("\n[Step 4] 获取缺失的FASTA序列...")
val_uniprots = set(r[4] for r in validated if r[4] and r[4] != '/')
missing_in_fasta = val_uniprots - set(sequences.keys())
print(f"  需要获取: {len(missing_in_fasta)} 个")

fetched = 0
failed_fetch = []
for uid in sorted(missing_in_fasta):
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
        req = Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
        with urlopen(req, timeout=10) as resp:
            fasta_text = resp.read().decode('utf-8')
            seq_lines = [l for l in fasta_text.strip().split('\n') if not l.startswith('>')]
            seq = ''.join(seq_lines)
            if seq:
                sequences[uid] = seq
                fetched += 1
                print(f"    OK {uid}: {len(seq)} aa")
    except Exception as e:
        failed_fetch.append((uid, str(e)[:60]))
        print(f"    FAIL {uid}: {str(e)[:60]}")
    time.sleep(0.3)

# Fallback: use blast column for still-missing
for r in validated:
    uid = r[4]
    if uid and uid != '/' and uid not in sequences:
        blast_seq = r[80] if len(r) > 80 else ''
        if blast_seq and len(blast_seq) > 50 and blast_seq[0].isalpha():
            sequences[uid] = blast_seq
            fetched += 1
            print(f"    BLAST {uid}: {len(blast_seq)} aa (from blast column)")

print(f"  成功获取: {fetched}, 失败: {len(failed_fetch)}")

# ============================================================
# Step 5 (Issue 2): Resolve UniProt for rows with CYP symbol
# ============================================================
print("\n[Step 5] 解析缺失UniProt的CYP符号...")
no_uniprot = [r for r in validated if not r[4] or r[4] == '/']
print(f"  待解析: {len(no_uniprot)} 行")

cyp_groups = {}
for r in no_uniprot:
    key = (r[1], r[6])
    if key not in cyp_groups:
        cyp_groups[key] = []
    cyp_groups[key].append(r)

print(f"  唯一 (CYP, species) 组合: {len(cyp_groups)} 个")

resolved_count = 0
unresolved_cyps = []

for (cyp, species), rows in sorted(cyp_groups.items()):
    if not cyp or cyp == '/':
        unresolved_cyps.append((cyp, species, len(rows), "empty CYP symbol"))
        continue

    # Try gene name search first
    query = f'gene_exact:{cyp}'
    if species:
        query += f' AND organism_name:"{species}"'
    try:
        search_url = f"https://rest.uniprot.org/uniprotkb/search?query={quote(query)}&format=json&size=5&fields=accession,gene_names,organism_name,protein_name,sequence"
        req = Request(search_url, headers={'User-Agent': 'Python/EZSpecificity'})
        with urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            results = data.get('results', [])

            if len(results) == 1:
                uid = results[0]['primaryAccession']
                seq = results[0].get('sequence', {}).get('value', '')
                for r in rows:
                    r[4] = uid
                if seq and uid not in sequences:
                    sequences[uid] = seq
                resolved_count += len(rows)
                print(f"    OK {cyp} ({species[:20] if species else '?'}): -> {uid} ({len(rows)} rows)")
            elif len(results) > 1:
                unresolved_cyps.append((cyp, species, len(rows), f"{len(results)} candidates"))
            else:
                # Try protein name search
                query2 = f'protein_name:{cyp}'
                if species:
                    query2 += f' AND organism_name:"{species}"'
                search_url2 = f"https://rest.uniprot.org/uniprotkb/search?query={quote(query2)}&format=json&size=5&fields=accession,gene_names,organism_name,protein_name,sequence"
                req2 = Request(search_url2, headers={'User-Agent': 'Python/EZSpecificity'})
                with urlopen(req2, timeout=10) as resp2:
                    data2 = json.loads(resp2.read().decode('utf-8'))
                    results2 = data2.get('results', [])
                    if len(results2) == 1:
                        uid = results2[0]['primaryAccession']
                        seq = results2[0].get('sequence', {}).get('value', '')
                        for r in rows:
                            r[4] = uid
                        if seq and uid not in sequences:
                            sequences[uid] = seq
                        resolved_count += len(rows)
                        print(f"    OK {cyp} ({species[:20] if species else '?'}): -> {uid} (protein name, {len(rows)} rows)")
                    else:
                        reason = "not found" if not results2 else f"{len(results2)} candidates"
                        unresolved_cyps.append((cyp, species, len(rows), reason))
    except Exception as e:
        unresolved_cyps.append((cyp, species, len(rows), str(e)[:50]))
    time.sleep(0.5)

print(f"  解析成功: {resolved_count} 行")
print(f"  未解析: {sum(n for _,_,n,_ in unresolved_cyps)} 行 ({len(unresolved_cyps)} 个CYP)")

# ============================================================
# Step 5b: Use sequence hash for rows still without UniProt but with blast sequence
# ============================================================
print("\n[Step 5b] 用序列hash恢复仍无UniProt的行...")
still_no_uniprot = [r for r in validated if not r[4] or r[4] == '/']
recovered_by_hash = 0
for r in still_no_uniprot:
    blast_seq = r[80] if len(r) > 80 else ''
    if blast_seq and len(blast_seq) > 50 and blast_seq[0].isalpha():
        seq_hash = 'SEQHASH_' + hashlib.sha256(blast_seq.upper().encode()).hexdigest()[:16]
        r[4] = seq_hash
        if seq_hash not in sequences:
            sequences[seq_hash] = blast_seq
        recovered_by_hash += 1

print(f"  通过序列hash恢复: {recovered_by_hash} 行")
print(f"  仍然无法恢复: {sum(1 for r in validated if not r[4] or r[4] == '/')} 行")

# ============================================================
# Step 6 (Issue 5): Resolve missing SMILES from PubChem
# ============================================================
print("\n[Step 6] 解析缺失的底物SMILES...")
no_smiles_rows = [r for r in validated if not r[14] or r[14] == '/']
print(f"  待解析: {len(no_smiles_rows)} 行")

resolved_smiles = 0
for r in no_smiles_rows:
    sub_name = r[9]
    if not sub_name or sub_name == '/':
        continue
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{quote(sub_name)}/property/CanonicalSMILES/JSON"
        req = Request(url, headers={'User-Agent': 'Python/EZSpecificity'})
        with urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            props = data.get('PropertyTable', {}).get('Properties', [])
            if props and props[0].get('CanonicalSMILES'):
                r[14] = props[0]['CanonicalSMILES']
                resolved_smiles += 1
                print(f"    OK {sub_name[:30]}: {r[14][:50]}")
    except:
        pass
    time.sleep(0.3)

print(f"  解析成功: {resolved_smiles} / {len(no_smiles_rows)}")

# ============================================================
# Step 7: Build final output
# ============================================================
print("\n[Step 7] 构建最终输出...")

pairs = {}
skipped = {'no_uniprot': 0, 'no_smiles': 0}

for r in validated:
    uniprot = r[4]
    smiles1 = r[14]

    if not uniprot or uniprot == '/':
        skipped['no_uniprot'] += 1
        continue
    if not smiles1 or smiles1 == '/':
        skipped['no_smiles'] += 1
        continue

    key = (uniprot, smiles1)
    if key not in pairs:
        pairs[key] = {
            'reaction_ids': [],
            'p450_symbol': r[1],
            'species': r[6],
            'species_class': r[7],
            'ec_number': r[5],
            'substrate_name': r[9],
            'pmids': set(),
            'products': [],
        }

    pairs[key]['reaction_ids'].append(r[0])
    if r[72] and r[72] != '/':
        pairs[key]['pmids'].add(r[72])
    for pi in [35, 41, 46, 51, 56, 61, 66]:
        if pi < len(r) and r[pi] and r[pi] != '/':
            if r[pi] not in pairs[key]['products']:
                pairs[key]['products'].append(r[pi])

print(f"  跳过: no_uniprot={skipped['no_uniprot']}, no_smiles={skipped['no_smiles']}")

# Write enzymes.csv
enzyme_map = {}
enzyme_rows = []
for uniprot, smiles in sorted(pairs.keys()):
    if uniprot not in enzyme_map:
        eid = len(enzyme_map)
        enzyme_map[uniprot] = eid
        info = pairs[(uniprot, smiles)]
        seq = sequences.get(uniprot, '')
        enzyme_rows.append({
            'enzyme_id': eid,
            'uniprot_id': uniprot,
            'p450_symbol': info['p450_symbol'],
            'species': info['species'],
            'species_class': info['species_class'],
            'ec_number': info['ec_number'],
            'sequence': seq,
            'sequence_length': len(seq),
        })

with open(os.path.join(OUT_DIR, 'enzymes.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['enzyme_id', 'uniprot_id', 'p450_symbol', 'species', 'species_class', 'ec_number', 'sequence', 'sequence_length'])
    w.writeheader()
    w.writerows(enzyme_rows)

# Write compounds.csv
compound_map = {}
compound_rows = []
for uniprot, smiles in sorted(pairs.keys()):
    if smiles not in compound_map:
        cid = len(compound_map)
        compound_map[smiles] = cid
        info = pairs[(uniprot, smiles)]
        compound_rows.append({
            'compound_id': cid,
            'smiles': smiles,
            'name': info['substrate_name'],
        })

with open(os.path.join(OUT_DIR, 'compounds.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['compound_id', 'smiles', 'name'])
    w.writeheader()
    w.writerows(compound_rows)

# Write interactions.csv
interaction_rows = []
for (uniprot, smiles), info in sorted(pairs.items()):
    interaction_rows.append({
        'interaction_id': len(interaction_rows),
        'enzyme_id': enzyme_map[uniprot],
        'compound_id': compound_map[smiles],
        'label': 1,
        'source': 'S3_P450Rdb',
        'quality_tier': 'B',
        'num_reactions': len(info['reaction_ids']),
        'has_pmid': 1 if info['pmids'] else 0,
        'has_products': 1 if info['products'] else 0,
    })

with open(os.path.join(OUT_DIR, 'interactions.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['interaction_id', 'enzyme_id', 'compound_id', 'label', 'source', 'quality_tier', 'num_reactions', 'has_pmid', 'has_products'])
    w.writeheader()
    w.writerows(interaction_rows)

# Write reactions.csv
reaction_rows = []
for (uniprot, smiles), info in sorted(pairs.items()):
    for rid in info['reaction_ids']:
        reaction_rows.append({
            'reaction_id': rid,
            'enzyme_id': enzyme_map[uniprot],
            'compound_id': compound_map[smiles],
            'substrate_smiles': smiles,
            'product_smiles': ';'.join(info['products']) if info['products'] else '',
            'pmids': ';'.join(info['pmids']) if info['pmids'] else '',
        })

with open(os.path.join(OUT_DIR, 'reactions.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['reaction_id', 'enzyme_id', 'compound_id', 'substrate_smiles', 'product_smiles', 'pmids'])
    w.writeheader()
    w.writerows(reaction_rows)

# Write unresolved log
with open(os.path.join(OUT_DIR, 'unresolved.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['cyp_symbol', 'species', 'row_count', 'reason'])
    for cyp, sp, n, reason in unresolved_cyps:
        w.writerow([cyp, sp, n, reason])

# ============================================================
# Final Summary
# ============================================================
enzymes_with_seq = sum(1 for e in enzyme_rows if e['sequence'])
print(f"\n{'=' * 60}")
print(f"最终提取结果")
print(f"{'=' * 60}")
print(f"enzymes.csv:      {len(enzyme_rows)} 个酶 ({enzymes_with_seq} 个有序列)")
print(f"compounds.csv:    {len(compound_rows)} 个化合物")
print(f"interactions.csv: {len(interaction_rows)} 条正样本对")
print(f"reactions.csv:    {len(reaction_rows)} 条反应记录")
print(f"unresolved.csv:   {len(unresolved_cyps)} 个未解析的CYP")
print(f"")
print(f"对比首次提取:")
print(f"  酶: 833 -> {len(enzyme_rows)} ({len(enzyme_rows)-833:+d})")
print(f"  化合物: 1,437 -> {len(compound_rows)} ({len(compound_rows)-1437:+d})")
print(f"  正样本: 2,710 -> {len(interaction_rows)} ({len(interaction_rows)-2710:+d})")
print(f"  反应: 3,248 -> {len(reaction_rows)} ({len(reaction_rows)-3248:+d})")
