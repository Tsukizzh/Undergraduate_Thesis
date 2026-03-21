"""
S1_RCSB: 从PathA已有数据中提取P450正样本到统一schema
使用B1_仅底物_271pos（PathA审核后的底物正样本）
"""
import csv, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
C2_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
PATHA_DIR = os.path.join(C2_DIR, '..', '..', 'PathA_2026-01-08_模型评估测试集构建')
if not os.path.exists(PATHA_DIR):
    PATHA_DIR = r'D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建'

DATA_DIR = os.path.join(PATHA_DIR, 'data', '04_Step4_格式修正后数据')
SRC_DIR = os.path.join(PATHA_DIR, 'source_data', '01_核心数据', '修复后最终版')
OUT_DIR = os.path.join(C2_DIR, 'data', 'sources', 'Source_RCSB')
os.makedirs(OUT_DIR, exist_ok=True)

TRIVIAL_SMILES = {'O', 'OO', 'O=O', '[O][O]', '[H]O[H]', '[OH2]', '[H][H]', 'N=N', '[H]O', 'O=[O]'}

print("=" * 60)
print("S1_RCSB P450 提取")
print("=" * 60)

# Step 1: Load PathA data
print("\n[Step 1] 加载PathA数据...")
with open(os.path.join(DATA_DIR, 'Enzymes.csv'), 'r', encoding='utf-8-sig') as f:
    enz_data = {r['Enzyme_Index']: r for r in csv.DictReader(f)}
print(f"  Enzymes.csv: {len(enz_data)} 个酶")

with open(os.path.join(DATA_DIR, 'Substrates.csv'), 'r', encoding='utf-8-sig') as f:
    sub_data = {r['Substrate_Index']: r['Substrate_SMILES'] for r in csv.DictReader(f)}
print(f"  Substrates.csv: {len(sub_data)} 个底物")

# Use B1 (substrate-only positives, curated in Path A)
b1_data_file = os.path.join(DATA_DIR, 'B1_仅底物_271pos', 'data.csv')
with open(b1_data_file, 'r', encoding='utf-8-sig') as f:
    b1_rows = list(csv.DictReader(f))
positives = [r for r in b1_rows if r['Label'] == '1']
print(f"  B1 正样本: {len(positives)} 条")

# Load 682 source for extra metadata (classification, pdb_id, ligand info)
src_682 = os.path.join(SRC_DIR, '独立测试集_682条.csv')
with open(src_682, 'r', encoding='utf-8-sig') as f:
    src_rows = list(csv.DictReader(f))
# Build lookup by (uniprot_id, ligand_ccd) or by pdb_id
pdb_info = {}
for r in src_rows:
    pdb_info[r['pdb_id']] = r

# Step 2: Filter trivial SMILES
print("\n[Step 2] 过滤无意义底物...")
excluded = []
valid = []
for r in positives:
    smi = sub_data.get(r['Substrate Index'], '')
    if smi in TRIVIAL_SMILES or len(smi) < 3:
        excluded.append(r)
        print(f"  排除: enzyme={r['Enzyme Index']}, SMILES='{smi}'")
    else:
        valid.append(r)
print(f"  有效: {len(valid)}, 排除: {len(excluded)}")

# Step 3: Build enzyme registry
print("\n[Step 3] 构建酶注册表...")
enzyme_map = {}
enzyme_rows = []
for r in valid:
    eidx = r['Enzyme Index']
    enz = enz_data.get(eidx, {})
    uid = enz.get('UniProt_ID', '')
    if uid and uid not in enzyme_map:
        eid = len(enzyme_map)
        enzyme_map[uid] = eid
        enzyme_rows.append({
            'enzyme_id': eid,
            'uniprot_id': uid,
            'p450_symbol': '',
            'species': enz.get('Organism', ''),
            'species_class': '',
            'ec_number': '',
            'sequence': enz.get('Protein sequence', ''),
            'sequence_length': enz.get('Sequence_Length', ''),
            'pdb_id': enz.get('PDB_ID', ''),
        })

print(f"  唯一酶: {len(enzyme_rows)}")
has_seq = sum(1 for e in enzyme_rows if e['sequence'] and len(e['sequence']) > 10)
print(f"  有序列: {has_seq}")

# Step 4: Build compound registry
print("\n[Step 4] 构建化合物注册表...")
compound_map = {}
compound_rows = []
for r in valid:
    smi = sub_data.get(r['Substrate Index'], '')
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
    eidx = r['Enzyme Index']
    sidx = r['Substrate Index']
    enz = enz_data.get(eidx, {})
    uid = enz.get('UniProt_ID', '')
    smi = sub_data.get(sidx, '')
    if not uid or not smi:
        continue
    pair_key = (uid, smi)
    if pair_key in seen_pairs:
        continue
    seen_pairs.add(pair_key)
    interaction_rows.append({
        'interaction_id': len(interaction_rows),
        'enzyme_id': enzyme_map[uid],
        'compound_id': compound_map[smi],
        'label': 1,
        'source': 'S1_RCSB',
        'quality_tier': 'A',
        'num_reactions': 1,
        'has_pmid': 0,
        'has_products': 0,
        'pdb_id': enz.get('PDB_ID', ''),
    })

print(f"  唯一配对: {len(interaction_rows)}")

# Step 6: Write output
print("\n[Step 6] 写入文件...")

with open(os.path.join(OUT_DIR, 'enzymes.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['enzyme_id', 'uniprot_id', 'p450_symbol', 'species', 'species_class', 'ec_number', 'sequence', 'sequence_length', 'pdb_id'])
    w.writeheader()
    w.writerows(enzyme_rows)

with open(os.path.join(OUT_DIR, 'compounds.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['compound_id', 'smiles', 'name'])
    w.writeheader()
    w.writerows(compound_rows)

with open(os.path.join(OUT_DIR, 'interactions.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['interaction_id', 'enzyme_id', 'compound_id', 'label', 'source', 'quality_tier', 'num_reactions', 'has_pmid', 'has_products', 'pdb_id'])
    w.writeheader()
    w.writerows(interaction_rows)

# Step 7: Overlap check
print("\n[Step 7] 与S2/S3重叠检查...")
s2_path = os.path.join(C2_DIR, 'data', 'sources', 'Source_ESIBank', 'enzymes.csv')
s3_path = os.path.join(C2_DIR, 'data', 'sources', 'Source_P450Rdb', 'enzymes.csv')

s1_enzymes = set(enzyme_map.keys())
s1_compounds = set(compound_map.keys())

for name, path in [('S2_ESIBank', s2_path), ('S3_P450Rdb', s3_path)]:
    if os.path.exists(path):
        with open(path) as f:
            other_enz = set(r['uniprot_id'] for r in csv.DictReader(f) if not r['uniprot_id'].startswith('SEQHASH'))
        comp_path = path.replace('enzymes.csv', 'compounds.csv')
        with open(comp_path) as f:
            other_comp = set(r['smiles'] for r in csv.DictReader(f))

        int_path = path.replace('enzymes.csv', 'interactions.csv')
        eid_map = {}
        with open(path) as f:
            for r in csv.DictReader(f):
                eid_map[r['enzyme_id']] = r['uniprot_id']
        cid_map = {}
        with open(comp_path) as f:
            for r in csv.DictReader(f):
                cid_map[r['compound_id']] = r['smiles']
        with open(int_path) as f:
            other_pairs = set()
            for r in csv.DictReader(f):
                u = eid_map.get(r['enzyme_id'], '')
                s = cid_map.get(r['compound_id'], '')
                if u and s:
                    other_pairs.add((u, s))

        enz_overlap = s1_enzymes & other_enz
        comp_overlap = s1_compounds & other_comp
        pair_overlap = seen_pairs & other_pairs
        print(f"  vs {name}: 酶重叠={len(enz_overlap)}/{len(s1_enzymes)} ({len(enz_overlap)/len(s1_enzymes)*100:.1f}%), "
              f"化合物重叠={len(comp_overlap)}/{len(s1_compounds)} ({len(comp_overlap)/len(s1_compounds)*100:.1f}%), "
              f"配对重叠={len(pair_overlap)}/{len(seen_pairs)} ({len(pair_overlap)/len(seen_pairs)*100:.1f}%)")

# Summary
print(f"\n{'=' * 60}")
print(f"S1_RCSB 提取结果")
print(f"{'=' * 60}")
print(f"enzymes.csv:      {len(enzyme_rows)} 个酶 ({has_seq} 个有序列)")
print(f"compounds.csv:    {len(compound_rows)} 个化合物")
print(f"interactions.csv: {len(interaction_rows)} 条正样本对 (Tier A)")
print(f"排除: {len(excluded)} 条 (辅因子/水)")
