import pandas as pd
import sys
import io
import os

# 设置输出编码
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 80)
print("生成P450 only数据集（740 → 682条）")
print("=" * 80)

# 路径
base_path = r"c:\Users\Administrator\Desktop\EZSpecificity_Project"
exp2_path = os.path.join(base_path, "提取P450过程日志", "2026-01-21_实验二v2_数据重构", "数据文件")
source_path = os.path.join(base_path, "P450_EZSpecificity_研究项目", "PathA_2026-01-08_模型评估测试集构建", "source_data")

# 1. 读取排除列表（58条Photosystem I）
exclude_df = pd.read_csv(
    os.path.join(base_path, "P450_EZSpecificity_研究项目", "PathA_2026-01-08_模型评估测试集构建",
                 "data", "03_Step3_配体分类审核", "05_排除记录", "非P450排除记录_58条_Photosystem_I.csv"),
    encoding='utf-8-sig'
)

exclude_pairs = set(zip(exclude_df['uniprot_id'], exclude_df['pdb_id'], exclude_df['ligand_ccd']))
print(f"\n[1/5] 读取排除列表")
print(f"  排除的非P450记录: {len(exclude_pairs)}条")
print(f"  UniProt IDs: {exclude_df['uniprot_id'].unique().tolist()}")

# 2. 过滤实验二v2的740条数据
print(f"\n[2/5] 过滤实验二v2主数据集")
df_740_exp2 = pd.read_csv(os.path.join(exp2_path, "交付物", "独立测试集_740条.csv"), encoding='utf-8-sig')
df_682_exp2 = df_740_exp2[~df_740_exp2.apply(
    lambda row: (row['uniprot_id'], row['pdb_id'], row['ligand_ccd']) in exclude_pairs, axis=1
)]
output_file_exp2 = os.path.join(exp2_path, "交付物", "独立测试集_v2_P450only_682条.csv")
df_682_exp2.to_csv(output_file_exp2, index=False, encoding='utf-8-sig')
print(f"  原始740条 → 过滤后{len(df_682_exp2)}条")
print(f"  已保存: 数据文件/交付物/独立测试集_v2_P450only_682条.csv")

# 3. 按物种过滤
print(f"\n[3/5] 按物种分类过滤")
species_files = {
    'bacteria': '独立测试集_bacteria_426条.csv',
    'fungi': '独立测试集_fungi_30条.csv',
    'mammals': '独立测试集_mammals_101条.csv',
    'other': '独立测试集_other_174条.csv',
    'plants': '独立测试集_plants_9条.csv'
}

species_stats_before = {}
species_stats_after = {}

for species, filename in species_files.items():
    # 读取原始文件
    df_species = pd.read_csv(os.path.join(exp2_path, "交付物", "按物种分类", filename), encoding='utf-8-sig')
    species_stats_before[species] = len(df_species)

    # 过滤
    df_species_filtered = df_species[~df_species.apply(
        lambda row: (row['uniprot_id'], row['pdb_id'], row['ligand_ccd']) in exclude_pairs, axis=1
    )]
    species_stats_after[species] = len(df_species_filtered)

    # 保存
    output_filename = filename.replace('.csv', '_P450only.csv')
    output_file = os.path.join(exp2_path, "交付物", "按物种分类", output_filename)
    df_species_filtered.to_csv(output_file, index=False, encoding='utf-8-sig')

    excluded_count = species_stats_before[species] - species_stats_after[species]
    print(f"  {species:10s}: {species_stats_before[species]:3d} → {species_stats_after[species]:3d} (排除{excluded_count:2d}条)")

# 4. 过滤PathA的source_data
print(f"\n[4/5] 过滤PathA source_data")
df_740_patha = pd.read_csv(os.path.join(source_path, "01_核心数据", "独立测试集_740条.csv"), encoding='utf-8-sig')
df_682_patha = df_740_patha[~df_740_patha.apply(
    lambda row: (row['uniprot_id'], row['pdb_id'], row['ligand_ccd']) in exclude_pairs, axis=1
)]
output_file_patha = os.path.join(source_path, "01_核心数据", "独立测试集_v2_P450only_682条.csv")
df_682_patha.to_csv(output_file_patha, index=False, encoding='utf-8-sig')
print(f"  原始740条 → 过滤后{len(df_682_patha)}条")
print(f"  已保存: source_data/01_核心数据/独立测试集_v2_P450only_682条.csv")

# 按物种过滤（PathA）
patha_species_path = os.path.join(source_path, "03_物种分类", "按物种分类")
for species, filename in species_files.items():
    df_species = pd.read_csv(os.path.join(patha_species_path, filename), encoding='utf-8-sig')
    df_species_filtered = df_species[~df_species.apply(
        lambda row: (row['uniprot_id'], row['pdb_id'], row['ligand_ccd']) in exclude_pairs, axis=1
    )]
    output_filename = filename.replace('.csv', '_P450only.csv')
    output_file = os.path.join(patha_species_path, output_filename)
    df_species_filtered.to_csv(output_file, index=False, encoding='utf-8-sig')

print(f"  PathA物种分类文件已更新")

# 5. 生成统计报告
print(f"\n[5/5] 生成统计报告")
print(f"\n物种分类统计（740 → 682）:")
print(f"  {'物种':12s}  {'原始':>6s}  {'过滤后':>6s}  {'排除':>6s}")
print(f"  {'-'*12}  {'-'*6}  {'-'*6}  {'-'*6}")
total_excluded = 0
for species in ['bacteria', 'fungi', 'mammals', 'plants', 'other']:
    before = species_stats_before[species]
    after = species_stats_after[species]
    excluded = before - after
    total_excluded += excluded
    print(f"  {species:12s}  {before:6d}  {after:6d}  {excluded:6d}")
print(f"  {'-'*12}  {'-'*6}  {'-'*6}  {'-'*6}")
print(f"  {'Total':12s}  {sum(species_stats_before.values()):6d}  {sum(species_stats_after.values()):6d}  {total_excluded:6d}")

print(f"\n验证:")
print(f"  实验二v2主数据集: {len(df_682_exp2)} 条")
print(f"  PathA主数据集: {len(df_682_patha)} 条")
print(f"  物种分类总和: {sum(species_stats_after.values())} 条")
print(f"  一致性检查: {'[OK]' if len(df_682_exp2) == len(df_682_patha) == sum(species_stats_after.values()) == 682 else '[ERROR]'}")

print(f"\n[DONE] 所有682条P450 only数据集已生成")
print("=" * 80)
