"""
S2_ESIBank: 从本地已有的P450提取文件中提取正样本到统一schema
"""
import csv, os, hashlib

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
C2_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
# Source file path - use absolute path to avoid Chinese path encoding issues
PROJECT_ROOT = os.path.abspath(os.path.join(C2_DIR, '..', '..', '..', '..', '..'))
SRC_FILE = os.path.join(PROJECT_ROOT, '毕业设计', 'P450_EZSpecificity_研究项目',
                         'PathA_2026-01-08_模型评估测试集构建',
                         'source_data', '02_底物数据', 'P450酶底物反应详表_完整版.csv')
if not os.path.exists(SRC_FILE):
    # Fallback: try direct D: path
    SRC_FILE = r'D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\source_data\02_底物数据\P450酶底物反应详表_完整版.csv'
OUT_DIR = os.path.join(C2_DIR, 'data', 'sources', 'Source_ESIBank')
os.makedirs(OUT_DIR, exist_ok=True)

TRIVIAL_SMILES = {'O', 'OO', 'O=O', '[O][O]', '[H]O[H]', '[OH2]', '[H][H]', 'N=N', '[H]O', 'O=[O]'}

print("=" * 60)
print("S2_ESIBank P450 提取")
print("=" * 60)

# Step 1: Load source data
print("\n[Step 1] 加载源数据...")
with open(SRC_FILE, 'r', encoding='utf-8') as f:
    all_rows = list(csv.DictReader(f))
print(f"  总行数: {len(all_rows)}")

positives = [r for r in all_rows if r['label'] == '1']
print(f"  正样本: {len(positives)}")

# Step 2: Filter trivial SMILES
print("\n[Step 2] 过滤无意义底物...")
excluded = []
valid = []
for r in positives:
    smi = r['substrate_smiles'].strip()
    if smi in TRIVIAL_SMILES or len(smi) < 3:
        excluded.append(r)
        print(f"  排除: UniProt={r['uniprot_id']}, SMILES='{smi}' (辅因子/水)")
    else:
        valid.append(r)

print(f"  有效正样本: {len(valid)}, 排除: {len(excluded)}")

# Step 3: Build enzyme registry
print("\n[Step 3] 构建酶注册表...")
enzyme_map = {}
enzyme_rows = []
for r in valid:
    uid = r['uniprot_id']
    if uid not in enzyme_map:
        eid = len(enzyme_map)
        enzyme_map[uid] = eid
        enzyme_rows.append({
            'enzyme_id': eid,
            'uniprot_id': uid,
            'p450_symbol': '',
            'species': '',
            'species_class': '',
            'ec_number': r['ecnumber'],
            'sequence': r['sequence'],
            'sequence_length': len(r['sequence']),
        })

print(f"  唯一酶: {len(enzyme_rows)}")
has_seq = sum(1 for e in enzyme_rows if e['sequence'] and len(e['sequence']) > 10)
print(f"  有序列: {has_seq}")

# Step 4: Build compound registry
print("\n[Step 4] 构建化合物注册表...")
compound_map = {}
compound_rows = []
for r in valid:
    smi = r['substrate_smiles'].strip()
    if smi not in compound_map:
        cid = len(compound_map)
        compound_map[smi] = cid
        compound_rows.append({
            'compound_id': cid,
            'smiles': smi,
            'name': '',
        })

print(f"  唯一化合物: {len(compound_rows)}")

# Step 5: Build interactions
print("\n[Step 5] 构建配对表...")
seen_pairs = set()
interaction_rows = []
for r in valid:
    uid = r['uniprot_id']
    smi = r['substrate_smiles'].strip()
    pair_key = (uid, smi)
    if pair_key in seen_pairs:
        continue
    seen_pairs.add(pair_key)
    interaction_rows.append({
        'interaction_id': len(interaction_rows),
        'enzyme_id': enzyme_map[uid],
        'compound_id': compound_map[smi],
        'label': 1,
        'source': 'S2_ESIBank',
        'quality_tier': 'B',
        'num_reactions': 1,
        'has_pmid': 0,
        'has_products': 0,
        'source_enzyme_index': r['enzyme'],
        'source_reaction_index': r['reaction'],
    })

print(f"  唯一配对: {len(interaction_rows)}")

# Step 5b: Build reactions table from G drive reaction.csv
print("\n[Step 5b] 从ESIBank提取反应SMILES...")
REACTION_FILE = r'G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\brenda\reaction.csv'
reaction_rows_out = []
if os.path.exists(REACTION_FILE):
    with open(REACTION_FILE, 'r', encoding='utf-8') as f:
        all_rxns = list(csv.DictReader(f))
    print(f"  ESIBank reaction.csv: {len(all_rxns)} 条")

    matched_rxn = 0
    for r in valid:
        uid = r['uniprot_id']
        smi = r['substrate_smiles'].strip()
        rxn_idx = int(r['reaction'])
        pair_key = (uid, smi)
        if 0 <= rxn_idx < len(all_rxns):
            rxn_smiles = all_rxns[rxn_idx]['reactions']
            if rxn_smiles:
                # Extract product from reaction SMILES (substrate>>product)
                parts = rxn_smiles.split('>>')
                product_smiles = parts[1] if len(parts) > 1 else ''
                reaction_rows_out.append({
                    'reaction_id': f'ESI_{rxn_idx}',
                    'enzyme_id': enzyme_map[uid],
                    'compound_id': compound_map[smi],
                    'substrate_smiles': smi,
                    'product_smiles': product_smiles,
                    'reaction_smiles': rxn_smiles,
                    'pmids': '',
                })
                matched_rxn += 1
    print(f"  匹配到反应SMILES: {matched_rxn}/{len(valid)}")
else:
    print(f"  G盘文件不可用，跳过反应SMILES提取")

# Step 6: Write output
print("\n[Step 6] 写入文件...")

with open(os.path.join(OUT_DIR, 'enzymes.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['enzyme_id', 'uniprot_id', 'p450_symbol', 'species', 'species_class', 'ec_number', 'sequence', 'sequence_length'])
    w.writeheader()
    w.writerows(enzyme_rows)

with open(os.path.join(OUT_DIR, 'compounds.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['compound_id', 'smiles', 'name'])
    w.writeheader()
    w.writerows(compound_rows)

with open(os.path.join(OUT_DIR, 'interactions.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['interaction_id', 'enzyme_id', 'compound_id', 'label', 'source', 'quality_tier', 'num_reactions', 'has_pmid', 'has_products', 'source_enzyme_index', 'source_reaction_index'])
    w.writeheader()
    w.writerows(interaction_rows)

if reaction_rows_out:
    with open(os.path.join(OUT_DIR, 'reactions.csv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=['reaction_id', 'enzyme_id', 'compound_id', 'substrate_smiles', 'product_smiles', 'reaction_smiles', 'pmids'])
        w.writeheader()
        w.writerows(reaction_rows_out)
    print(f"  reactions.csv: {len(reaction_rows_out)} 条")

# Step 7: Quick overlap check with S3_P450Rdb
print("\n[Step 7] 与S3_P450Rdb重叠检查...")
s3_enzyme_path = os.path.join(C2_DIR, 'data', 'sources', 'Source_P450Rdb', 'enzymes.csv')
s3_compound_path = os.path.join(C2_DIR, 'data', 'sources', 'Source_P450Rdb', 'compounds.csv')
s3_interaction_path = os.path.join(C2_DIR, 'data', 'sources', 'Source_P450Rdb', 'interactions.csv')

if os.path.exists(s3_enzyme_path):
    with open(s3_enzyme_path) as f:
        s3_enzymes = set(r['uniprot_id'] for r in csv.DictReader(f) if not r['uniprot_id'].startswith('SEQHASH'))
    with open(s3_compound_path) as f:
        s3_compounds = set(r['smiles'] for r in csv.DictReader(f))

    s2_enzymes = set(enzyme_map.keys())
    s2_compounds = set(compound_map.keys())

    enzyme_overlap = s2_enzymes & s3_enzymes
    compound_overlap = s2_compounds & s3_compounds

    print(f"  S2酶: {len(s2_enzymes)}, S3酶(UniProt): {len(s3_enzymes)}")
    print(f"  酶重叠: {len(enzyme_overlap)} ({len(enzyme_overlap)/len(s2_enzymes)*100:.1f}% of S2)")
    print(f"  S2化合物: {len(s2_compounds)}, S3化合物: {len(s3_compounds)}")
    print(f"  化合物重叠: {len(compound_overlap)} ({len(compound_overlap)/len(s2_compounds)*100:.1f}% of S2)")

    # Pair overlap (need to load S3 interactions and map back to UniProt+SMILES)
    with open(s3_enzyme_path) as f:
        s3_eid_to_uid = {r['enzyme_id']: r['uniprot_id'] for r in csv.DictReader(f)}
    with open(s3_compound_path) as f:
        s3_cid_to_smi = {r['compound_id']: r['smiles'] for r in csv.DictReader(f)}
    with open(s3_interaction_path) as f:
        s3_pairs = set()
        for r in csv.DictReader(f):
            uid = s3_eid_to_uid.get(r['enzyme_id'], '')
            smi = s3_cid_to_smi.get(r['compound_id'], '')
            if uid and smi:
                s3_pairs.add((uid, smi))

    pair_overlap = seen_pairs & s3_pairs
    print(f"  配对重叠: {len(pair_overlap)} ({len(pair_overlap)/len(seen_pairs)*100:.1f}% of S2)")
else:
    print("  S3数据未找到，跳过重叠检查")

# Summary
print(f"\n{'=' * 60}")
print(f"S2_ESIBank 提取结果")
print(f"{'=' * 60}")
print(f"enzymes.csv:      {len(enzyme_rows)} 个酶 (全部有序列)")
print(f"compounds.csv:    {len(compound_rows)} 个化合物")
print(f"interactions.csv: {len(interaction_rows)} 条正样本对 (Tier B)")
print(f"reactions.csv:    {len(reaction_rows_out)} 条反应记录 (含反应SMILES)")
print(f"排除: {len(excluded)} 条 (辅因子/水)")
