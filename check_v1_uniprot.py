import pandas as pd

# 读取v1数据
df = pd.read_csv(r'C:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\_archived_v1_按UniProt去重_Step1至3\source_data\01_核心数据\ML训练数据集_186个.csv')

# 25个与ESIBank交集的UniProt
overlap_uniprots = ['P19099', 'Q9UNU6', 'P08684', 'P11509', 'P05108', 'P11712', 'P05181', 'Q02318', 'P20813', 'P14779', 'P15128', 'P33268', 'P08183', 'P51587', 'P20852', 'Q16678', 'P48422', 'P33261', 'P15538', 'Q9HBI6', 'P33260', 'P15391', 'P79401', 'P24462', 'P00186']

# 排除交集后的数据
df_filtered = df[~df['uniprot_id'].isin(overlap_uniprots)]

print(f'v1版本原始数据:')
print(f'  总记录数: {len(df)}')
print(f'  唯一PDB: {df["pdb_id"].nunique()}')
print(f'  唯一UniProt: {df["uniprot_id"].nunique()}')
print()
print(f'v1版本排除{len(overlap_uniprots)}个交集后:')
print(f'  剩余记录数: {len(df_filtered)}')
print(f'  唯一PDB: {df_filtered["pdb_id"].nunique()}')
print(f'  唯一UniProt: {df_filtered["uniprot_id"].nunique()}')
