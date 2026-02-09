"""
统计B3数据集（539条）的物种分布 - 修正版
"""
import pandas as pd
from pathlib import Path

# 路径配置
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # PathA_2026-01-08_模型评估测试集构建

# 1. 读取源数据（682条，含物种信息）
source_682 = pd.read_csv(BASE_DIR / "source_data/01_核心数据/修复后最终版/独立测试集_682条.csv")
print(f"源数据: {len(source_682)}条")

# 2. 读取三方验证结果
verified = pd.read_csv(BASE_DIR / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/03_合并结果/三方深度验证结果_merged.csv")
print(f"三方验证结果: {len(verified)}条")

# 3. 用pdb_id关联（正确的方式！）
print(f"\nverified pdb_id示例: {verified['pdb_id'].head().tolist()}")
print(f"source_682 pdb_id示例: {source_682['pdb_id'].head().tolist()}")

# 用pdb_id作为关联键
merged = verified.merge(source_682[['pdb_id', 'organism', 'species_category', 'uniprot_id']],
                        on='pdb_id', how='left', suffixes=('_verified', '_source'))

print(f"合并后: {len(merged)}条")
print(f"species_category非空: {merged['species_category'].notna().sum()}条")

# 4. 筛选B3数据集条件: SUBSTRATE + INHIBITOR + PRODUCT
b3_classes = ['SUBSTRATE', 'INHIBITOR', 'PRODUCT']
b3_data = merged[merged['final_class'].isin(b3_classes)]
print(f"\nB3数据集（SUBSTRATE+INHIBITOR+PRODUCT）: {len(b3_data)}条")

# 5. 统计物种分布
species_stats = b3_data['species_category'].value_counts()
print(f"\n{'='*50}")
print(f"B3数据集（539条）物种分布")
print(f"{'='*50}")
for species, count in species_stats.items():
    pct = count / len(b3_data) * 100
    print(f"  {species:12}: {count:4}条 ({pct:5.1f}%)")
print(f"  {'总计':12}: {len(b3_data):4}条")

# 6. 按物种和分类统计
print(f"\n{'='*50}")
print(f"按物种和配体类型交叉统计")
print(f"{'='*50}")
cross_tab = pd.crosstab(b3_data['species_category'], b3_data['final_class'], margins=True)
print(cross_tab.to_string())

# 7. 唯一酶（UniProt）的物种分布
unique_enzymes = b3_data.drop_duplicates(subset=['uniprot_id_source'])
enzyme_species = unique_enzymes['species_category'].value_counts()
print(f"\n{'='*50}")
print(f"唯一P450酶物种分布（{len(unique_enzymes)}个）")
print(f"{'='*50}")
for species, count in enzyme_species.items():
    pct = count / len(unique_enzymes) * 100
    print(f"  {species:12}: {count:4}个 ({pct:5.1f}%)")

# 8. 保存详细结果
output_dir = BASE_DIR / "data/11_物种统计"
output_dir.mkdir(parents=True, exist_ok=True)

# 保存B3数据集的物种信息
b3_with_species = b3_data[['record_id', 'pdb_id', 'uniprot_id_source', 'enzyme_name',
                           'ligand_name', 'final_class', 'organism', 'species_category']].copy()
b3_with_species.columns = ['record_id', 'pdb_id', 'uniprot_id', 'enzyme_name',
                           'ligand_name', 'final_class', 'organism', 'species_category']
b3_with_species.to_csv(output_dir / "B3_539条_物种分布详情.csv", index=False, encoding='utf-8-sig')

# 保存统计摘要
summary = pd.DataFrame({
    '物种类别': species_stats.index,
    '记录数': species_stats.values,
    '占比(%)': (species_stats.values / len(b3_data) * 100).round(1)
})
summary.to_csv(output_dir / "B3_539条_物种统计摘要.csv", index=False, encoding='utf-8-sig')

# 保存酶的物种统计
enzyme_summary = pd.DataFrame({
    '物种类别': enzyme_species.index,
    '酶数量': enzyme_species.values,
    '占比(%)': (enzyme_species.values / len(unique_enzymes) * 100).round(1)
})
enzyme_summary.to_csv(output_dir / "B3_唯一酶_物种统计.csv", index=False, encoding='utf-8-sig')

# 保存交叉表
cross_tab.to_csv(output_dir / "B3_物种配体类型交叉表.csv", encoding='utf-8-sig')

print(f"\n结果已保存到: {output_dir}")
