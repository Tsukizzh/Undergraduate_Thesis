"""
Step 3.5: 数据集重构脚本

根据导师要求（2026-01-20组会）：
1. 方案B拆分为B1/B2/B3三种子方案
2. 去冗余逻辑从按UniProt分组改为按(UniProt + 配体CCD)分组

Author: Claude Code
Date: 2026-01-20
"""

import pandas as pd
import json
import re
from pathlib import Path

# 路径配置（中文文件夹名称）
BASE_DIR = Path(__file__).parent.parent
INPUT_DIR = BASE_DIR / "data" / "03_Step3.5_数据集重构" / "输入_源数据"
INTERMEDIATE_DIR = BASE_DIR / "data" / "03_Step3.5_数据集重构" / "中间结果"
OUTPUT_DIR = BASE_DIR / "data" / "03_Step3.5_数据集重构" / "输出_最终数据"
SCHEMES_DIR = BASE_DIR / "data" / "04_Step3.5_方案数据集"

# 排除列表：辅因子、溶剂、缓冲液、离子等
EXCLUDE_LIGANDS = {
    # 血红素相关
    'HEM', 'HEC', 'HEA', 'HEB', 'HAS', 'HDD', 'HDM', 'HEG', 'HNI',
    # 辅酶
    'FAD', 'FMN', 'NAD', 'NAP', 'NDP', 'ADP', 'ATP', 'GTP', 'GDP',
    # 水和离子
    'HOH', 'H2O', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN', 'FE', 'MN', 'CU', 'NI', 'CO',
    # 缓冲盐
    'SO4', 'PO4', 'ACT', 'ACE', 'FMT', 'TRS', 'MES', 'HEP', 'CIT',
    # 实验添加剂/溶剂
    'GOL', 'PEG', 'PGE', 'MPD', 'DMS', 'DMF', 'P6G', 'PG4', '1PE', '2PE',
    'EDO', 'EOH', 'MOH', 'IPA', 'IPH', 'ETA', 'BME',
    # 常见缓冲组分
    'IMD', 'EPE', 'DTT', 'MLI', 'IOD', 'SCN', 'CYS', 'SER', 'GLY',
    # 洗涤剂
    'BOG', 'LDA', 'LMT', 'OLC', 'C8E', 'C10', 'C12', 'DDM', 'OG', 'DM',
    # 其他常见非底物
    'NO', 'CN', 'O2', 'NO2', 'HIS', 'LI', 'NH4', 'OCT', 'UNX', 'UNL',
    # 气体分子（特别排除）
    'CMO',  # Carbon Monoxide
    'PEO',  # Peroxide相关
}

# 已知P450抑制剂关键词
INHIBITOR_KEYWORDS = [
    'inhibitor', 'antagonist', 'blocker', 'ketoconazole', 'itraconazole',
    'fluconazole', 'voriconazole', 'clotrimazole', 'miconazole', 'econazole',
    'ritonavir', 'cobicistat', 'erythromycin', 'clarithromycin', 'cimetidine',
    'azole', 'imidazole', 'triazole'
]

# 底物关键词
SUBSTRATE_KEYWORDS = [
    'substrate', 'camphor', 'fatty acid', 'steroid', 'hormone',
    'testosterone', 'progesterone', 'cortisol', 'cholesterol',
    'arachidonic', 'lauric', 'palmitic', 'linoleic'
]

# 产物关键词（使用正则边界匹配）
PRODUCT_PATTERNS = [
    r'\bproduct\b',
    r'\bmetabolite\b',
    r'\bhydroxyl\b',
    r'\bepoxide\b',
    r'\boxidized\b',
    r'hydroxy-',
    r'-ol\b',
]

# 产物排除词（避免误匹配）
PRODUCT_EXCLUDE_PATTERNS = [
    r'\bmonoxide\b',
    r'\bperoxide\b',
    r'\bdioxide\b',
    r'\bhydroxide\b',
    r'carbon monoxide',
    r'nitric oxide',
    r'hydrogen peroxide',
]


def parse_ligands(ligand_str):
    """解析配体字符串，返回配体列表 [(ccd, name), ...]"""
    if pd.isna(ligand_str) or not ligand_str.strip():
        return []

    ligands = []
    for part in ligand_str.split(';'):
        part = part.strip()
        if ':' in part:
            ccd, name = part.split(':', 1)
            ccd = ccd.strip().upper()
            name = name.strip()
            ligands.append((ccd, name))
    return ligands


def parse_uniprot_ids(uniprot_str):
    """解析UniProt ID字符串（JSON格式）"""
    if pd.isna(uniprot_str) or not uniprot_str.strip():
        return []
    try:
        ids = json.loads(uniprot_str.replace("'", '"'))
        return ids if isinstance(ids, list) else [ids]
    except:
        return []


def classify_ligand(ccd, name, title=""):
    """
    对配体进行分类
    返回: 'substrate', 'product', 'inhibitor', 'unknown', 'exclude'
    """
    ccd_upper = ccd.upper()
    name_lower = name.lower() if name else ""
    title_lower = title.lower() if title else ""
    combined = f"{name_lower} {title_lower}"

    # 1. 检查排除列表
    if ccd_upper in EXCLUDE_LIGANDS:
        return 'exclude'

    # 2. 检查名称中的排除词（避免误分类）
    for pattern in PRODUCT_EXCLUDE_PATTERNS:
        if re.search(pattern, combined, re.IGNORECASE):
            return 'exclude'

    # 3. 检查抑制剂关键词
    for kw in INHIBITOR_KEYWORDS:
        if kw in combined:
            return 'inhibitor'

    # 4. 检查底物关键词
    for kw in SUBSTRATE_KEYWORDS:
        if kw in combined:
            return 'substrate'

    # 5. 检查产物关键词（使用正则边界）
    for pattern in PRODUCT_PATTERNS:
        if re.search(pattern, combined, re.IGNORECASE):
            return 'product'

    # 6. 默认unknown
    return 'unknown'


def save_scheme(df, scheme_dir, scheme_name):
    """保存方案数据集（仅正样本）"""
    scheme_dir = Path(scheme_dir)
    scheme_dir.mkdir(parents=True, exist_ok=True)

    # 保存完整数据
    df.to_csv(scheme_dir / f"dataset_{scheme_name}.csv", index=False)

    # 生成Enzymes.csv格式（用于EZSpecificity）
    enzymes = df[['pdb_id', 'uniprot_id', 'sequence']].drop_duplicates(subset=['pdb_id', 'uniprot_id'])
    enzymes = enzymes.rename(columns={'sequence': 'Protein sequence'})
    enzymes.to_csv(scheme_dir / "Enzymes.csv", index=False)

    # 生成Substrates.csv格式（需要后续补充SMILES）
    substrates = df[['ligand_ccd', 'ligand_name']].drop_duplicates()
    substrates['Substrate_SMILES'] = ''  # 需要后续从API获取
    substrates.to_csv(scheme_dir / "Substrates_template.csv", index=False)

    # 统计信息
    print(f"    {scheme_name}: {len(df)} 条记录 -> {scheme_dir.name}/")


def main():
    print("=" * 60)
    print("Step 3.5: 数据集重构 (修复版)")
    print("=" * 60)

    # 1. 读取输入数据
    print("\n[1] 读取输入数据...")
    df_structures = pd.read_csv(INPUT_DIR / "高质量P450全集_1426个.csv")
    df_esibank = pd.read_csv(INPUT_DIR / "P450酶列表_最终版389个.csv")

    print(f"    高质量结构数: {len(df_structures)}")
    print(f"    ESIBank P450数: {len(df_esibank)}")

    # 获取ESIBank UniProt ID集合
    esibank_uniprots = set(df_esibank['uniprot_id'].tolist())
    print(f"    ESIBank UniProt ID数: {len(esibank_uniprots)}")

    # 2. PDB级别严格排除（修复多UniProt映射数据泄露）
    print("\n[2] PDB级别严格排除检查...")

    def pdb_has_esibank_uniprot(uniprot_ids_str):
        """检查PDB是否包含任何ESIBank UniProt ID"""
        uids = parse_uniprot_ids(uniprot_ids_str)
        return any(uid in esibank_uniprots for uid in uids)

    # 找出所有包含ESIBank UniProt的PDB
    pdbs_with_esibank = set(
        df_structures.loc[
            df_structures['uniprot_ids'].apply(pdb_has_esibank_uniprot),
            'pdb_id'
        ].astype(str).str.upper().tolist()
    )
    print(f"    包含ESIBank UniProt的PDB数: {len(pdbs_with_esibank)}")

    # 3. 展开配体，创建(UniProt, 配体CCD)记录
    print("\n[3] 展开配体，创建酶-配体对记录...")
    records = []

    for _, row in df_structures.iterrows():
        pdb_id = row['pdb_id']
        resolution = row['resolution']
        title = row['title']
        sequence = row['sequence']
        seq_length = row['sequence_length']
        organism = row['organism_name']
        species_cat = row['species_category']

        # 解析UniProt IDs
        uniprot_ids = parse_uniprot_ids(row['uniprot_ids'])
        if not uniprot_ids:
            continue

        # 解析配体
        ligands = parse_ligands(row['ligands'])
        if not ligands:
            continue

        # 为每个(UniProt, 配体)组合创建记录
        for uniprot_id in uniprot_ids:
            for ccd, ligand_name in ligands:
                # 分类配体
                classification = classify_ligand(ccd, ligand_name, title)

                records.append({
                    'pdb_id': pdb_id,
                    'uniprot_id': uniprot_id,
                    'ligand_ccd': ccd,
                    'ligand_name': ligand_name,
                    'classification': classification,
                    'resolution': resolution,
                    'title': title,
                    'sequence': sequence,
                    'sequence_length': seq_length,
                    'organism': organism,
                    'species_category': species_cat,
                })

    df_all = pd.DataFrame(records)
    print(f"    总记录数（展开后）: {len(df_all)}")

    # 保存展开后的完整列表
    df_all.to_csv(INTERMEDIATE_DIR / "master_ligand_list.csv", index=False)
    print(f"    已保存: intermediate/master_ligand_list.csv")

    # 4. 统计分类结果
    print("\n[4] 配体分类统计...")
    class_counts = df_all['classification'].value_counts()
    for cls, cnt in class_counts.items():
        print(f"    {cls}: {cnt}")

    # 5. 过滤掉exclude类别
    print("\n[5] 过滤排除类别...")
    df_valid = df_all[df_all['classification'] != 'exclude'].copy()
    print(f"    有效记录数: {len(df_valid)}")

    # 6. 新去冗余逻辑：按(UniProt, 配体CCD)分组
    print("\n[6] 执行新去冗余逻辑（按UniProt+配体分组）...")

    # 按分辨率排序（升序，小的更好）
    df_valid = df_valid.sort_values('resolution')

    # 按(uniprot_id, ligand_ccd)分组，保留第一个（分辨率最好的）
    df_dedup = df_valid.drop_duplicates(subset=['uniprot_id', 'ligand_ccd'], keep='first')
    print(f"    去冗余后记录数: {len(df_dedup)}")

    # 保存去冗余日志
    dedup_log = df_valid.groupby(['uniprot_id', 'ligand_ccd']).agg({
        'pdb_id': list,
        'resolution': list
    }).reset_index()
    dedup_log['selected_pdb'] = dedup_log['pdb_id'].apply(lambda x: x[0])
    dedup_log['candidate_count'] = dedup_log['pdb_id'].apply(len)
    dedup_log.to_csv(INTERMEDIATE_DIR / "deduplication_log.csv", index=False)
    print(f"    已保存: intermediate/deduplication_log.csv")

    # 7. 排除ESIBank交集（PDB级别严格排除）
    print("\n[7] 排除ESIBank交集（PDB级别严格排除）...")
    df_before = df_dedup.copy()

    # PDB级别排除
    pdb_upper = df_before['pdb_id'].astype(str).str.upper()
    df_independent = df_before[~pdb_upper.isin(pdbs_with_esibank)].copy()

    # 二次安全检查：排除任何直接匹配ESIBank UniProt的记录
    df_independent = df_independent[~df_independent['uniprot_id'].isin(esibank_uniprots)].copy()

    excluded_count = len(df_before) - len(df_independent)
    excluded_pdbs = df_before[pdb_upper.isin(pdbs_with_esibank)]['pdb_id'].nunique()
    excluded_uniprots = df_before[df_before['uniprot_id'].isin(esibank_uniprots)]['uniprot_id'].nunique()
    print(f"    排除记录数: {excluded_count}")
    print(f"    排除PDB ID数（PDB级别）: {excluded_pdbs}")
    print(f"    排除UniProt ID数: {excluded_uniprots}")
    print(f"    独立测试集记录数: {len(df_independent)}")

    # 保存独立测试集
    df_independent.to_csv(OUTPUT_DIR / "deduplicated_pairs.csv", index=False)
    print(f"    已保存: output/deduplicated_pairs.csv")

    # 8. 统计最终结果
    print("\n[8] 最终统计...")
    print(f"    总酶-配体对数: {len(df_independent)}")
    print(f"    唯一UniProt ID数: {df_independent['uniprot_id'].nunique()}")
    print(f"    唯一配体CCD数: {df_independent['ligand_ccd'].nunique()}")
    print(f"    唯一PDB ID数: {df_independent['pdb_id'].nunique()}")

    final_class_counts = df_independent['classification'].value_counts()
    print("\n    分类分布:")
    for cls, cnt in final_class_counts.items():
        print(f"      {cls}: {cnt}")

    # 9. 生成B1/B2/B3数据集（仅正样本）
    print("\n[9] 生成B1/B2/B3数据集...")

    # B1: 仅底物
    df_b1 = df_independent[df_independent['classification'] == 'substrate'].copy()
    save_scheme(df_b1, SCHEMES_DIR / "方案B1_仅底物", "B1")

    # B2: 底物+产物
    df_b2 = df_independent[df_independent['classification'].isin(['substrate', 'product'])].copy()
    save_scheme(df_b2, SCHEMES_DIR / "方案B2_底物加产物", "B2")

    # B3: 底物+产物+抑制剂
    df_b3 = df_independent[df_independent['classification'].isin(['substrate', 'product', 'inhibitor'])].copy()
    save_scheme(df_b3, SCHEMES_DIR / "方案B3_全部配体", "B3")

    # 10. 生成需要人工复核的列表
    print("\n[10] 生成人工复核列表...")
    df_unknown = df_independent[df_independent['classification'] == 'unknown']
    df_unknown_review = df_unknown[['pdb_id', 'uniprot_id', 'ligand_ccd', 'ligand_name', 'title', 'resolution']].copy()
    df_unknown_review['suggested_class'] = 'unknown'
    df_unknown_review['your_decision'] = ''  # 留空给用户填写
    df_unknown_review.to_csv(INTERMEDIATE_DIR / "unknown_review_list.csv", index=False)
    print(f"    需要复核的unknown配体数: {len(df_unknown_review)}")
    print(f"    已保存: intermediate/unknown_review_list.csv")

    # 11. 总结
    print("\n" + "=" * 60)
    print("Step 3.5 执行完成!")
    print("=" * 60)

    print(f"""
总结：
├── 输入：{len(df_structures)} 个高质量P450结构
├── 展开后：{len(df_all)} 条酶-配体记录
├── 过滤exclude后：{len(df_valid)} 条
├── 去冗余后：{len(df_dedup)} 条（按UniProt+配体）
├── 排除ESIBank交集后：{len(df_independent)} 条（独立测试集）
│   └── PDB级别严格排除：{excluded_pdbs} 个PDB
│
├── B1（仅底物）：{len(df_b1)} 条
├── B2（底物+产物）：{len(df_b2)} 条
├── B3（底物+产物+抑制剂）：{len(df_b3)} 条
│
└── 需要人工复核：{len(df_unknown_review)} 个unknown配体

下一步：
1. 检查 intermediate/unknown_review_list.csv，确认/修改分类
2. Step 4: 构建索引表（data.csv）
3. Step 5-8: 特征生成
""")


if __name__ == "__main__":
    main()
