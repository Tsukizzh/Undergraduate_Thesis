# 任务06 - 数据集切分 文件索引

## 脚本文件
| 文件 | 路径 | 说明 |
|------|------|------|
| build_datasets.py | scripts/.../06_数据集切分/ | 数据集构建主脚本 |

## 输入文件
| 文件 | 路径 | 说明 |
|------|------|------|
| chunk_*.csv | data/.../02_配体分类并行验证/01_输入数据_chunks/ | 源数据(record_id → pdb_id, ligand_ccd) |
| 三方深度验证结果_merged.csv | data/.../04_配体三方深度验证/03_合并结果/ | 分类结果(record_id, final_class) |
| ccd_smiles.csv | data/.../05_提取SMILES/ | CCD→SMILES映射 |
| Enzymes.csv | data/02_Step2_酶序列/ | 酶索引表 |
| record_enzyme_mapping.csv | data/02_Step2_酶序列/ | (pdb_id,ligand_ccd)→Enzyme_Index |

## 输出文件
| 文件 | 路径 | 正样本 | 负样本 | 说明 |
|------|------|--------|--------|------|
| Substrates.csv | data/.../06_数据集切分/ | - | - | 436个唯一SMILES |
| data.csv | .../B1_仅底物_272pos/ | 272 | 0 | 仅SUBSTRATE |
| data.csv | .../B2_底物正产物负_272pos23neg/ | 272 | 23 | SUBSTRATE正+PRODUCT负 |
| data.csv | .../B3_完整数据集_272pos267neg/ | 272 | 267 | SUBSTRATE正+(PRODUCT+INHIBITOR)负 |
| data.csv | .../B4_仅产物_23neg/ | 0 | 23 | 仅PRODUCT(负样本) |
| data.csv | .../B5_仅抑制剂_244neg/ | 0 | 244 | 仅INHIBITOR(负样本) |
| data.csv | .../B6_底物正抑制剂负_272pos244neg/ | 272 | 244 | SUBSTRATE正+INHIBITOR负 |

## 方案说明 (2026-01-29修正)

**重要**: PRODUCT被标记为**负样本**，因为：
- 产物是催化反应的**结果**，不是底物
- EZSpecificity原始训练数据中，只有底物是正样本，产物被完全排除
- 参考：Codex分析 `src/Datasets/create_features.py:27` 强制 `substrate_only=True`

## Session文件
| 文件 | 路径 | 说明 |
|------|------|------|
| session_log.md | sessions/.../06_数据集切分/ | 执行日志 |
| file_index.md | sessions/.../06_数据集切分/ | 本文件 |

---
**最后更新**: 2026-01-29 | **版本**: v2.2
