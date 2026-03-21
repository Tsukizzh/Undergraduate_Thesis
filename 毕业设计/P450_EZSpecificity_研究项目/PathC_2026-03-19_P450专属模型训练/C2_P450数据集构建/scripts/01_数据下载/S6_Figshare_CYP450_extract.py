"""
S6_Figshare_CYP450: 提取6个人类CYP的底物/非底物数据
用途：辅助数据源（化合物池扩展 + 确认生物学负样本），不进入主benchmark
"""
import csv, os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
C2_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DL_DIR = os.path.join(C2_DIR, 'downloads', 'Figshare_CYP450')
OUT_DIR = os.path.join(C2_DIR, 'data', 'sources', 'Source_Figshare_CYP450')
os.makedirs(OUT_DIR, exist_ok=True)

CYP_TO_UNIPROT = {
    'CYP1A2': 'P05177',
    'CYP2C9': 'P11712',
    'CYP2C19': 'P33261',
    'CYP2D6': 'P10635',
    'CYP2E1': 'P05181',
    'CYP3A4': 'P08684',
}

print("=" * 60)
print("S6_Figshare_CYP450 提取（辅助数据源）")
print("=" * 60)

# Step 1: Load all CSV files
print("\n[Step 1] 加载数据...")
all_pairs = []  # (cyp, name, smiles, label, source_paper)
for cyp in CYP_TO_UNIPROT:
    for split in ['training', 'testing']:
        fname = os.path.join(DL_DIR, f'{cyp}_{split}set.csv')
        if not os.path.exists(fname):
            print(f"  WARNING: {fname} not found")
            continue
        with open(fname, 'r', encoding='utf-8-sig') as f:
            for r in csv.DictReader(f):
                all_pairs.append({
                    'cyp': cyp,
                    'name': r.get('Name', ''),
                    'smiles': r.get('SMILES', ''),
                    'label': int(r.get('Label', '0')),
                    'source_paper': r.get('Source', ''),
                })

print(f"  总记录: {len(all_pairs)}")
positives = [p for p in all_pairs if p['label'] == 1]
negatives = [p for p in all_pairs if p['label'] == 0]
print(f"  底物 (label=1): {len(positives)}")
print(f"  非底物 (label=0): {len(negatives)}")

# Step 2: Build enzyme registry (6 CYPs, no sequences yet)
print("\n[Step 2] 构建酶注册表...")
enzyme_rows = []
enzyme_map = {}
for i, (cyp, uid) in enumerate(CYP_TO_UNIPROT.items()):
    enzyme_map[cyp] = i
    enzyme_rows.append({
        'enzyme_id': i,
        'uniprot_id': uid,
        'p450_symbol': cyp,
        'species': 'Homo sapiens',
        'species_class': 'Animal',
        'ec_number': '',
        'sequence': '',
        'sequence_length': 0,
    })
print(f"  酶数量: {len(enzyme_rows)}")

# Step 3: Build compound registry
print("\n[Step 3] 构建化合物注册表...")
compound_map = {}
compound_rows = []
for p in all_pairs:
    smi = p['smiles']
    if smi and smi not in compound_map:
        cid = len(compound_map)
        compound_map[smi] = cid
        compound_rows.append({
            'compound_id': cid,
            'smiles': smi,
            'name': p['name'],
        })
print(f"  唯一化合物: {len(compound_rows)}")

# Step 4: Build interactions (positives only)
print("\n[Step 4] 构建正样本配对表...")
seen_pos = set()
pos_rows = []
for p in positives:
    smi = p['smiles']
    cyp = p['cyp']
    key = (cyp, smi)
    if key in seen_pos:
        continue
    seen_pos.add(key)
    pos_rows.append({
        'interaction_id': len(pos_rows),
        'enzyme_id': enzyme_map[cyp],
        'compound_id': compound_map[smi],
        'label': 1,
        'source': 'S6_Figshare',
        'quality_tier': 'B',
        'usage': 'auxiliary',
        'source_paper': p['source_paper'],
    })
print(f"  唯一正样本配对: {len(pos_rows)}")

# Step 5: Build biological negatives (separate file)
print("\n[Step 5] 构建生物学负样本表...")
seen_neg = set()
neg_rows = []
for p in negatives:
    smi = p['smiles']
    cyp = p['cyp']
    key = (cyp, smi)
    if key in seen_neg:
        continue
    seen_neg.add(key)
    neg_rows.append({
        'interaction_id': len(neg_rows),
        'enzyme_id': enzyme_map[cyp],
        'compound_id': compound_map[smi],
        'label': 0,
        'source': 'S6_Figshare',
        'quality_tier': 'B',
        'usage': 'biological_negative',
        'source_paper': p['source_paper'],
    })
print(f"  唯一负样本配对: {len(neg_rows)}")

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

with open(os.path.join(OUT_DIR, 'interactions_positive.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['interaction_id', 'enzyme_id', 'compound_id', 'label', 'source', 'quality_tier', 'usage', 'source_paper'])
    w.writeheader()
    w.writerows(pos_rows)

with open(os.path.join(OUT_DIR, 'biological_negatives.csv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['interaction_id', 'enzyme_id', 'compound_id', 'label', 'source', 'quality_tier', 'usage', 'source_paper'])
    w.writeheader()
    w.writerows(neg_rows)

# Step 7: Overlap with S2/S3
print("\n[Step 7] 与S2/S3化合物重叠...")
s6_sub_smiles = set(p['smiles'] for p in positives)
for name, src in [('S2_ESIBank', 'Source_ESIBank'), ('S3_P450Rdb', 'Source_P450Rdb')]:
    comp_path = os.path.join(C2_DIR, 'data', 'sources', src, 'compounds.csv')
    if os.path.exists(comp_path):
        with open(comp_path) as f:
            other_smi = set(r['smiles'] for r in csv.DictReader(f))
        overlap = s6_sub_smiles & other_smi
        print(f"  vs {name}: {len(overlap)} 个化合物重叠")

new_compounds = len(compound_map) - len(s6_sub_smiles & set())
print(f"  S6带来的新化合物: {len(compound_rows)} 个（含正负样本）")

# Summary
print(f"\n{'=' * 60}")
print(f"S6_Figshare_CYP450 提取结果")
print(f"{'=' * 60}")
print(f"enzymes.csv:              6 个人类CYP（无序列，需从UniProt获取）")
print(f"compounds.csv:            {len(compound_rows)} 个化合物")
print(f"interactions_positive.csv: {len(pos_rows)} 条正样本对 (Tier B, auxiliary)")
print(f"biological_negatives.csv:  {len(neg_rows)} 条确认负样本 (Tier B, biological_negative)")
print(f"")
print(f"注意：此数据源标记为 auxiliary，不进入主benchmark合并")
