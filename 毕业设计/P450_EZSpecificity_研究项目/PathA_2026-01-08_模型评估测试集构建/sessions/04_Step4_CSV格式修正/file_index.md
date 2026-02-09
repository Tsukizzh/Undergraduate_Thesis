# Step 4 文件索引

> **更新时间**: 2026-01-29

---

## 脚本文件

| 文件 | 用途 |
|------|------|
| [step4_csv_format_fix.py](../../scripts/04_Step4_CSV格式修正/step4_csv_format_fix.py) | CSV格式修正主脚本 |

---

## 输出数据

### 数据集文件

| 文件 | 行数 | 说明 |
|------|------|------|
| [B1_仅底物_272pos/data.csv](../../data/04_Step4_格式修正后数据/B1_仅底物_272pos/data.csv) | 272 | 仅正样本（底物） |
| [B2_底物正产物负_272pos23neg/data.csv](../../data/04_Step4_格式修正后数据/B2_底物正产物负_272pos23neg/data.csv) | 295 | 底物正+产物负 |
| [B3_完整数据集_272pos267neg/data.csv](../../data/04_Step4_格式修正后数据/B3_完整数据集_272pos267neg/data.csv) | 539 | 完整数据集 |
| [B4_仅产物_23neg/data.csv](../../data/04_Step4_格式修正后数据/B4_仅产物_23neg/data.csv) | 23 | 仅负样本（产物） |
| [B5_仅抑制剂_244neg/data.csv](../../data/04_Step4_格式修正后数据/B5_仅抑制剂_244neg/data.csv) | 244 | 仅负样本（抑制剂） |
| [B6_底物正抑制剂负_272pos244neg/data.csv](../../data/04_Step4_格式修正后数据/B6_底物正抑制剂负_272pos244neg/data.csv) | 516 | 底物正+抑制剂负 |

### 辅助文件

| 文件 | 行数 | 说明 |
|------|------|------|
| [Enzymes.csv](../../data/04_Step4_格式修正后数据/Enzymes.csv) | 292 | 酶序列表（从Step 2复制） |
| [Substrates.csv](../../data/04_Step4_格式修正后数据/Substrates.csv) | 436 | 底物SMILES表（从Step 3复制） |

### 报告文件

| 文件 | 说明 |
|------|------|
| [format_fix_report.txt](../../data/04_Step4_格式修正后数据/format_fix_report.txt) | 格式修正报告 |

---

## 数据格式说明

### data.csv 列定义

| 列名 | 类型 | 范围 | 说明 |
|------|------|------|------|
| Dock Index | int64 | 0-538 | 对接复合物索引 |
| Enzyme Index | int64 | 0-291 | 酶索引 |
| Substrate Index | int64 | 0-435 | 底物索引 |
| Label | int64 | 0/1 | 标签 |

### 与其他步骤的关联

```
Enzymes.csv (Enzyme_Index) ← data.csv (Enzyme Index) → enzyme_features.lmdb (Step 5)
Substrates.csv (Substrate_Index) ← data.csv (Substrate Index) → reaction_features.lmdb (Step 6)
PDB文件 ← data.csv (Dock Index) → structure_features.lmdb (Step 8)
```

---

## 核心参考文档

> **Q6: 完整代码库分析** 位于 `reports/数据集格式与负样本生成完整分析报告.md`
>
> 包含：
> - 6.1 代码库完整文件清单
> - 6.2 端到端数据流程详解（含ASCII流程图）
> - 6.3 CSV列名的铁证（代码追溯）
> - 6.4 LMDB为什么不可跳过
> - 6.5 配置文件详解
> - 6.6 example.ipynb步骤详解
> - 6.7 结构文件格式要求
> - 6.8 Index的含义详解
> - 6.9 模型输出解读（score → 概率转换）
> - 6.10 常见错误与解决方案
> - 6.11 核心结论汇总

---

**文档版本**: v1.0
