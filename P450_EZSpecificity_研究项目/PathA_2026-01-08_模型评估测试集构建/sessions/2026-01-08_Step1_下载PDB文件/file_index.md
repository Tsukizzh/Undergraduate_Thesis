# 文件索引

## Step 1 产生的文件

| 文件路径 | 说明 |
|---------|------|
| `data/00_raw/pdb_files/*.pdb` | 154个PDB格式结构文件 |
| `data/00_raw/pdb_files/*.cif` | 5个mmCIF格式结构文件 |
| `logs/excluded_pdbs.csv` | 被排除的26条记录列表 |
| `logs/errors/failed_pdbs.txt` | 下载失败记录（已通过.cif解决） |
| `scripts/step1_download_pdb.py` | 下载脚本 |

## 文件统计

| 类型 | 数量 |
|------|------|
| .pdb文件 | 154 |
| .cif文件 | 5 |
| **总计** | **159** |
