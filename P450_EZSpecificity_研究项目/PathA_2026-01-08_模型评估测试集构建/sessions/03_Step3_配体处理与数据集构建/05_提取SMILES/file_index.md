# 任务05 - 提取SMILES 文件索引

## 目录结构

```
05_提取SMILES/
├── data/           # 输出数据
├── scripts/        # 脚本文件
└── sessions/       # 会话记录
```

## 文件列表

### data/03_Step3_配体处理与数据集构建/05_提取SMILES/
| 文件名 | 说明 | 行数 |
|--------|------|------|
| ccd_smiles.csv | CCD代码对应的SMILES映射表 | 544 (含表头) |

### scripts/03_Step3_配体处理与数据集构建/05_提取SMILES/
| 文件名 | 说明 |
|--------|------|
| extract_smiles_from_ccd.py | SMILES提取主脚本 (RCSB API + RDKit验证) |
| check_columns.py | 辅助脚本：检查CSV列名 |

### sessions/03_Step3_配体处理与数据集构建/05_提取SMILES/
| 文件名 | 说明 |
|--------|------|
| session_log.md | 任务执行日志 |
| file_index.md | 本文件 |

## 输出数据格式

### ccd_smiles.csv
```csv
ccd_code,smiles,source,validation_status,error_message
TRP,N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O,API,VALID,
...
```

- **ccd_code**: 3字母PDB配体代码
- **smiles**: 标准化的SMILES字符串 (canonical, isomeric)
- **source**: 数据来源 (API/PDB_file/FAILED)
- **validation_status**: RDKit验证结果 (VALID/INVALID)
- **error_message**: 错误信息（如有）
