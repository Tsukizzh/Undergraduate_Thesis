# Step 7.2: GROVER指纹 - 文件索引

> **最后更新**: 2026-01-30

## 脚本文件

| 文件 | 路径 | 描述 |
|------|------|------|
| step7_2_generate_grover.py | scripts/07_Step7_分子指纹生成/ | GROVER指纹编排脚本 |

## 数据文件

### 输入
| 文件 | 路径 | 描述 |
|------|------|------|
| Substrates.csv | data/04_Step4_格式修正后数据/ | 436个底物及其SMILES |
| grover_large.pt | data/pretrain_model/ | 预训练的GROVER-Large模型（约1GB）|

### 输出
| 文件 | 路径 | 大小 | 描述 |
|------|------|------|------|
| grover_fingerprint.lmdb | data/07_Step7_分子指纹/ | 约10 GB | GROVER指纹（436条记录）|
| grover_substrates.csv | data/07_Step7_分子指纹/ | 20 KB | GROVER输入CSV |
| grover_substrates.npz | data/07_Step7_分子指纹/ | 3.7 KB | 原子/键描述符 |
| test_atom_vocab.pkl | data/07_Step7_分子指纹/grover_vocab/ | <1 KB | 原子词汇表（150种类型）|
| test_bond_vocab.pkl | data/07_Step7_分子指纹/grover_vocab/ | <1 KB | 键词汇表（212种类型）|

### 修改的源文件
| 文件 | 路径 | 修改内容 |
|------|------|----------|
| fingerprint.py | src/other_softwares/grover_software/task/ | map_size: 600GB → 10GB |
| data_representer.py | src/Datasets/ | map_size: 600GB → 10GB（3处）|

## 会话文档
| 文件 | 路径 | 描述 |
|------|------|------|
| session_log.md | sessions/07_Step7_分子指纹生成/02_GROVER指纹/ | 本次会话日志 |
| file_index.md | sessions/07_Step7_分子指纹生成/02_GROVER指纹/ | 本文件索引 |
