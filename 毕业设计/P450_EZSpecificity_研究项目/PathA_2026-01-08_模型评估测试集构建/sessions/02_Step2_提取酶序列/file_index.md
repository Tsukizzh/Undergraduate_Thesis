# Step 2 v2 文件索引

> **创建时间**: 2026-01-21
> **步骤**: Step 2 - 提取酶序列

---

## 输入文件

| 文件名 | 位置 | 说明 |
|--------|------|------|
| `独立测试集_740条.csv` | `source_data/01_核心数据/` | 实验二v2产出的核心清单，包含740条记录 |

### 输入文件关键字段

| 字段名 | 说明 | 用途 |
|--------|------|------|
| `sequence` | RCSB deposited sequence | **序列来源**（不是从PDB文件提取！） |
| `pdb_id` | PDB ID | 用于关联记录 |
| `uniprot_id` | UniProt ID | 用于统计和验证 |
| `organism` | 物种信息 | 传递到输出文件 |

---

## 输出文件

| 文件名 | 位置 | 说明 | 行数 |
|--------|------|------|------|
| `Enzymes.csv` | `data/02_Step2_酶序列/` | 唯一酶序列（EZSpecificity格式） | 294 |
| `record_enzyme_mapping.csv` | `data/02_Step2_酶序列/` | 记录与酶索引的映射 | 740 |

### Enzymes.csv 字段说明

| 字段 | 说明 | 示例 |
|------|------|------|
| `Enzyme_Index` | 酶索引（0-293） | 0, 1, 2, ... |
| `Protein sequence` | 蛋白质序列（EZSpecificity要求的列名） | MALDP... |
| `UniProt_ID` | UniProt ID | P00183 |
| `PDB_ID` | 代表性PDB ID（该序列的第一个出现） | 1AKD |
| `Organism` | 物种 | Bacillus megaterium |
| `Sequence_Length` | 序列长度 | 472 |

### record_enzyme_mapping.csv 字段说明

| 字段 | 说明 |
|------|------|
| `Record_Index` | 原740条记录的索引（0-739） |
| `pdb_id` | PDB ID |
| `Enzyme_Index` | 对应的酶索引（指向Enzymes.csv） |

---

## 数据关系说明

```
数据去重逻辑：

740条记录
    ↓
按 sequence 字段去重（完全相同的序列才视为同一个酶）
    ↓
294条唯一序列

为什么294 > 155（唯一UniProt数）？
┌─────────────────────────────────────────────────────────────┐
│ 原因1：突变体                                                │
│   - 同一个UniProt ID可能对应多个PDB结构                       │
│   - 不同PDB可能是野生型、突变体、或不同表达条件的结构          │
│   - 例如：P00183可能有5个不同的PDB，其中2个是突变体            │
│   - 野生型和突变体的序列不同，因此算作不同的酶                 │
├─────────────────────────────────────────────────────────────┤
│ 原因2：表达系统差异                                           │
│   - 不同表达系统（大肠杆菌 vs 酵母）可能产生略有差异的序列     │
│   - 如N端标签、C端截短等                                      │
├─────────────────────────────────────────────────────────────┤
│ 原因3：deposited sequence变体                                │
│   - RCSB的deposited sequence是实验者提交的实际表达序列        │
│   - 同一个蛋白在不同实验室可能有不同的构建体                  │
└─────────────────────────────────────────────────────────────┘
```

---

## 重要技术说明

### 为什么不从PDB文件提取序列？

```
v1阶段的验证结果：

在158个PDB文件中，比较了两种序列来源：
- SEQRES（PDB文件中的序列记录）
- CSV的sequence字段（RCSB deposited sequence）

结果：
- 150/158 个序列长度不一致！
- 仅 8 个完全一致

原因分析：
┌─────────────────────────────────────────────────────────────┐
│ SEQRES（PDB文件）         vs    deposited sequence（CSV）    │
├─────────────────────────────────────────────────────────────┤
│ - 可能只包含晶体中可见的残基      - 完整的表达序列           │
│ - 无序区域可能被省略             - 包含标签、接头等           │
│ - 可能被手动编辑过               - RCSB API直接获取           │
└─────────────────────────────────────────────────────────────┘

结论：使用CSV的sequence字段更可靠、更一致
```

---

## 序列统计

| 项目 | 数值 |
|------|------|
| 最短序列 | 368 aa |
| 最长序列 | 786 aa |
| 平均长度 | 463.7 aa |

### P450酶的典型序列长度

- 细菌P450：~400-500 aa
- 真核P450：~480-560 aa（含膜锚定区）
- 本数据集范围：368-786 aa，符合P450特征

---

## 处理流程

```python
# 伪代码展示处理逻辑

# 1. 读取输入
df = read_csv("独立测试集_740条.csv")

# 2. 按序列去重
unique_sequences = df.drop_duplicates(subset=['sequence'])
# 结果：294条唯一序列

# 3. 创建Enzymes.csv
enzymes_df = DataFrame({
    'Enzyme_Index': range(294),
    'Protein sequence': unique_sequences['sequence'],
    'UniProt_ID': unique_sequences['uniprot_id'],
    'PDB_ID': unique_sequences['pdb_id'],  # 取第一个出现的PDB
    'Organism': unique_sequences['organism'],
    'Sequence_Length': unique_sequences['sequence'].apply(len)
})

# 4. 创建映射表
mapping_df = df.merge(enzymes_df[['sequence', 'Enzyme_Index']],
                       left_on='sequence', right_on='sequence')
```

---

## v1 vs v2 对比

| 项目 | v1 | v2 | 变化 |
|------|----|----|------|
| 输入记录数 | 158 | 740 | +582 (+368%) |
| 唯一酶序列数 | 158 | 294 | +136 (+86%) |
| 序列来源 | CSV sequence | CSV sequence | 不变 |
| 去重方式 | 按UniProt（上游已去重） | 按sequence（本步骤去重） | 逻辑改变 |

### 为什么v2的唯一酶数量变化不成比例？

```
v1: 158条记录 → 158个唯一酶 (1:1)
v2: 740条记录 → 294个唯一酶 (2.52:1)

原因：
v1是按UniProt去重的结果，每个UniProt只保留一条记录
v2是按(UniProt+配体)去重的结果，同一个UniProt可能有多条记录

平均每个酶对应 740/294 ≈ 2.52 条记录
这意味着平均每个酶与2-3种不同的配体有共晶结构
```

---

**文档版本**: v1.0 | **创建时间**: 2026-01-21
