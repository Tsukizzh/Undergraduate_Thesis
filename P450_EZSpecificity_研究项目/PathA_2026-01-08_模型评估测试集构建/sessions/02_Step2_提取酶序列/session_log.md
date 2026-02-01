# Step 2 v2 执行日志：提取酶序列

> **执行时间**: 2026-01-21 (初版) → 2026-01-29 (净化更新)
> **版本**: v2.1 (基于682条净化后数据，排除58条Photosystem I)

---

## 一、执行结果摘要

| 项目 | 初版(740条) | 净化后(682条) | 变化 |
|------|-------------|---------------|------|
| 输入记录数 | 740 | **682** | -58 (PS-I) |
| 唯一酶序列数 | 294 | **292** | -2 (W8SUA3, W8SY74) |
| 平均每个酶的记录数 | 2.52 | 2.34 | - |

### 序列统计

| 项目 | 数值 |
|------|------|
| 最短序列 | 368 aa |
| 最长序列 | 786 aa |
| 平均长度 | 463.7 aa |

---

## 二、详细执行过程

### 2.1 输入分析

**输入文件**: `source_data/01_核心数据/修复后最终版/独立测试集_682条.csv`

⚠️ **数据净化更新** (2026-01-29):
- 原始版本：`独立测试集_740条.csv` (含58条Photosystem I污染)
- 净化版本：`独立测试集_682条.csv` (已排除W8SUA3、W8SY74共58条PS-I)
- Step 2数据已同步更新至净化版本

**关键字段**:
| 字段名 | 用途 |
|--------|------|
| `sequence` | **序列来源**（RCSB deposited sequence） |
| `pdb_id` | 用于关联和参考 |
| `uniprot_id` | 用于统计和参考 |
| `organism` | 传递到输出 |

### 2.2 处理逻辑

```python
# 伪代码展示处理过程

# Step 1: 读取输入
df = pd.read_csv("独立测试集_740条.csv")
print(f"总记录数: {len(df)}")  # 740

# Step 2: 检查唯一序列数
unique_sequences = df['sequence'].unique()
print(f"唯一序列数: {len(unique_sequences)}")  # 294

# Step 3: 创建Enzymes.csv
enzymes = df.drop_duplicates(subset=['sequence'], keep='first')
enzymes = enzymes.reset_index(drop=True)
enzymes['Enzyme_Index'] = enzymes.index
enzymes['Sequence_Length'] = enzymes['sequence'].apply(len)

# 重命名为EZSpecificity格式
enzymes = enzymes.rename(columns={
    'sequence': 'Protein sequence',
    'uniprot_id': 'UniProt_ID',
    'pdb_id': 'PDB_ID',
    'organism': 'Organism'
})

# Step 4: 创建映射表
sequence_to_index = dict(zip(enzymes['Protein sequence'], enzymes['Enzyme_Index']))
df['Enzyme_Index'] = df['sequence'].map(sequence_to_index)
mapping = df[['pdb_id', 'Enzyme_Index']].reset_index()
mapping.columns = ['Record_Index', 'pdb_id', 'Enzyme_Index']
```

### 2.3 关键决策记录

#### 决策1：序列来源选择

| 选项 | 描述 | 采用 |
|------|------|------|
| 从PDB文件提取SEQRES | 解析PDB文件的SEQRES记录 | ❌ |
| 使用CSV的sequence字段 | RCSB API返回的deposited sequence | ✅ |

**选择理由**:
- v1阶段验证显示：150/158个PDB的SEQRES与CSV的sequence不一致
- SEQRES可能被手动编辑、可能缺少未解析残基
- CSV的sequence是RCSB官方API返回的完整表达序列

#### 决策2：去重策略

| 选项 | 描述 | 采用 |
|------|------|------|
| 按UniProt去重 | 同一UniProt只保留一条 | ❌ |
| 按sequence去重 | 完全相同的序列才合并 | ✅ |

**选择理由**:
- 同一UniProt可能有不同的序列变体（突变体、不同构建体）
- 不同序列对模型来说是不同的输入，应该分别保留
- 这样保留了294个唯一序列，而不是155个UniProt

---

## 三、输出文件详情

### 3.1 Enzymes.csv

**位置**: `data/02_Step2_酶序列/Enzymes.csv`
**行数**: **292行**（不含表头）← 净化后（原294行，删除2条PS-I序列）

| 字段 | 类型 | 说明 | 示例 |
|------|------|------|------|
| `Enzyme_Index` | int | 酶索引（0-291） | 0 |
| `Protein sequence` | str | 蛋白质序列 | MTIKEMPQPKT... |
| `UniProt_ID` | str | UniProt ID | P00183 |
| `PDB_ID` | str | 代表性PDB ID | 1AKD |
| `Organism` | str | 物种 | Bacillus megaterium |
| `Sequence_Length` | int | 序列长度 | 472 |

**EZSpecificity兼容性**:
- 列名`Protein sequence`是EZSpecificity要求的格式
- `Enzyme_Index`与后续data.csv中的索引对应

### 3.2 record_enzyme_mapping.csv

**位置**: `data/02_Step2_酶序列/record_enzyme_mapping.csv`
**行数**: **682行**（不含表头）← 净化后（原740行，删除58条PS-I记录）

| 字段 | 说明 |
|------|------|
| `Record_Index` | 净化后682条记录的索引（0-681） |
| `pdb_id` | PDB ID |
| `Enzyme_Index` | 对应的酶索引（指向Enzymes.csv） |

**用途**: 在后续Step中建立data.csv时，可通过此映射找到每条记录对应的酶索引

---

## 四、数据关系说明

### 4.1 为什么292 > 155？

```
682条记录（净化后）
    │
    ├── 155个唯一UniProt ID
    │
    └── 292个唯一sequence
         │
         └── 原因：同一UniProt有多个序列变体
```

**具体原因**:

1. **突变体研究**
   ```
   UniProt P00183 (P450-BM3) 有多个PDB：
   - 1AKD: 野生型
   - 2HPD: F87A突变体
   - 4KEY: A82F/F87V双突变体
   这些序列不同，算作不同的酶
   ```

2. **不同表达构建体**
   ```
   不同实验室可能使用：
   - 不同的N端标签
   - 不同的截短版本
   - 不同的表达系统优化密码子
   ```

3. **deposited sequence的差异**
   ```
   RCSB记录的是实验者提交的实际表达序列
   同一蛋白在不同实验可能有微小差异
   ```

### 4.2 数据流示意

```
独立测试集_682条.csv (净化版)
         │
         │ 读取sequence字段
         │
         ▼
    682条序列记录
         │
         │ 按sequence去重（完全相同才合并）
         │
         ▼
    292条唯一序列 ───→ Enzymes.csv
         │
         │ 建立索引映射
         │
         ▼
    682条记录的映射 ───→ record_enzyme_mapping.csv
```

**Photosystem I排除说明** (2026-01-24~29):
- 发现58条记录为Photosystem I蛋白（W8SUA3、W8SY74），非P450
- 这些记录在740→682净化过程中被排除
- 对应的2条酶序列同步从Enzymes.csv中删除
- 详见：`reports/03_Step3_配体处理与数据集构建/01_配体分类审核/Photosystem_I污染根因分析.md`

---

## 五、v1 vs v2 对比

| 项目 | v1 | v2 (净化后) | 变化 |
|------|----|----|------|
| 输入记录数 | 158 | 682 | +524 (+332%) |
| 唯一酶序列数 | 158 | 292 | +134 (+85%) |
| 序列来源 | CSV sequence | CSV sequence | 相同 |
| 去重依据 | 上游已按UniProt去重 | 本步骤按sequence去重 | 不同 |

**为什么变化率不同？**
- 记录数变化332%（从158到682），但
- 唯一酶数量只变化85%（从158到292）
- 这是因为v2中平均每个酶对应2.34条记录

---

## 六、验证检查

### 6.1 数据完整性

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| Enzymes.csv行数 | 292 | 292 | ✅ |
| 所有序列非空 | 是 | 是 | ✅ |
| Enzyme_Index连续 | 0-291 | 0-291 | ✅ |
| 映射表覆盖所有记录 | 682 | 682 | ✅ |
| 不含PS-I (W8SUA3) | 0条 | 0条 | ✅ |
| 不含PS-I (W8SY74) | 0条 | 0条 | ✅ |

### 6.2 格式兼容性

| 检查项 | EZSpecificity要求 | 实际 | 状态 |
|--------|-------------------|------|------|
| 列名 | "Protein sequence" | "Protein sequence" | ✅ |
| 序列格式 | 仅标准氨基酸 | 仅标准氨基酸 | ✅ |
| 索引从0开始 | 是 | 是 | ✅ |

---

## 七、下一步

Step 3: 提取配体SMILES（等待用户确认后执行）

需要确认的问题：
1. 选择使用哪个方案（B1/B2/B3）？
2. 如何处理448条unknown配体？
3. SMILES获取失败时的处理策略？

---

**执行者**: Claude Code + Codex + Gemini
**日期**: 2026-01-21 (初版) → 2026-01-29 (净化更新)
**文档版本**: v2.1 (基于682条净化后数据)
