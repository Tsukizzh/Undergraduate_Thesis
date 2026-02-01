# Step 6 文件索引

> **更新时间**: 2026-01-30
> **状态**: ✅ 已完成

---

## 执行摘要

| 项目 | 数值 |
|------|------|
| 输入底物数量 | 436 |
| 成功生成 | 436 (100%) |
| 输出文件 | reaction_features.lmdb |
| 报告文件 | reaction_features_report.csv |
| 处理时间 | ~2秒 |
| GPU使用 | 无（纯CPU计算） |

---

## 脚本文件

| 文件 | 用途 |
|------|------|
| [step6_generate_reaction_graph_features.py](../../scripts/06_Step6_反应图特征生成/step6_generate_reaction_graph_features.py) | 反应图特征生成主脚本 |

---

## 输入数据

| 文件 | 来源 | 说明 |
|------|------|------|
| [Substrates.csv](../../data/04_Step4_格式修正后数据/Substrates.csv) | Step 4 | 436条底物SMILES |

---

## 输出数据

| 文件 | 记录数 | 说明 |
|------|--------|------|
| [reaction_features.lmdb](../../data/06_Step6_反应图特征/reaction_features.lmdb) | 436 | 分子图特征（节点+边+原子特征） |
| [reaction_features_report.csv](../../data/06_Step6_反应图特征/reaction_features_report.csv) | 436 | 处理状态报告 |

---

## LMDB格式说明

### Key-Value结构

| 字段 | 类型 | Shape | 说明 |
|------|------|-------|------|
| key | bytes | - | `str(Substrate_Index).encode()` |
| element | np.int64 | (n,) | 每个原子的元素序号 |
| edge_index | np.int64 | (2, E) | 有向边的起点和终点索引 |
| edge_type | np.int64 | (E,) | 化学键类型编号 |
| atom_feature | np.int64 | (n, 7) | 每个原子的7维特征 |
| reaction_attention_label | np.float32 | (n,) | 全0（substrate-only模式） |
| num_nodes | np.int64 | 标量 | 原子总数 |

### atom_feature 7维定义

| 维度 | 含义 | 取值示例 |
|------|------|----------|
| [0] | 原子序数 | 6=碳, 7=氮, 8=氧 |
| [1] | 是否芳香 | 0或1 |
| [2] | 连接度(degree) | 1~4 |
| [3] | 相邻氢原子数 | 0~3 |
| [4] | SP杂化 | 0或1 |
| [5] | SP2杂化 | 0或1 |
| [6] | SP3杂化 | 0或1 |

### edge_type 编号

| 编号 | 键类型 |
|------|--------|
| 1 | 单键 (Single) |
| 2 | 双键 (Double) |
| 3 | 三键 (Triple) |
| 12 | 芳香键 (Aromatic) |

---

## 与其他步骤的关联

```
Substrates.csv (Substrate_Index, Substrate_SMILES)
    ↓ Step 6 (RDKit解析 → 图特征)
reaction_features.lmdb (key = Substrate_Index)
    ↓ 模型推理
data.csv (Substrate Index) → 查找对应图特征
```

---

## 技术规格

| 项目 | 值 |
|------|-----|
| 解析工具 | RDKit 2025.03.6 |
| 存储格式 | LMDB (map_size=4GB) |
| 序列化 | pickle (HIGHEST_PROTOCOL) |
| 最大原子数限制 | 280 |
| torch_scatter | 未使用（自定义numpy实现） |

---

**文档版本**: v1.0
