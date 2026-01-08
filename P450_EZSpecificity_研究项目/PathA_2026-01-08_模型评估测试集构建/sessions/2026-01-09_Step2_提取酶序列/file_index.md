# 文件索引

## Step 2 产出文件

| 文件路径 | 说明 |
|---------|------|
| data/01_processed/Enzymes.csv | 158条酶序列（来自源数据） |
| reports/tables/step2_pdb_sequence_mapping.csv | **PDB ID对应关系审计文件** |
| reports/tables/step2_sequence_extraction_report.csv | 原始提取报告（已弃用） |
| scripts/step2_extract_sequence.py | 原始提取脚本（已弃用） |

## 更新的文件

| 文件路径 | 说明 |
|---------|------|
| source_data/01_核心数据/ML训练数据集_186个.csv | 更新为184条（移除6ZZX） |
| source_data/README.md | 添加sequence列用途说明 |
| data/00_raw/pdb_files/ | 删除6ZZX.pdb |
| 进度日志.md | 更新Step 2状态和技术修正 |
| 全局进度日志.md | 添加Step 0-2执行记录 |
