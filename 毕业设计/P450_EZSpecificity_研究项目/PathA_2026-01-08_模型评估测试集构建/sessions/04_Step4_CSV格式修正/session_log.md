# Step 4 执行日志：CSV格式修正

> **执行时间**: 2026-01-29
> **版本**: v1.0
> **执行者**: Claude Code + Codex + Gemini

---

## 一、执行结果摘要

| 项目 | 数值 |
|------|------|
| 处理数据集数 | 6 (B1-B6) |
| 处理辅助文件 | 2 (Enzymes.csv, Substrates.csv) |
| 修正问题 | UTF-8 BOM移除 |
| 验证状态 | ALL PASSED |

---

## 二、问题分析

### 2.1 初始问题识别

通过三方协作（Claude + Codex + Gemini）确认了以下格式问题：

| 问题 | 描述 | 影响 |
|------|------|------|
| UTF-8 BOM | 所有CSV文件开头含`EF BB BF`字节 | 可能导致列名解析错误 |

### 2.2 关键发现

经过与Codex和Gemini的多轮讨论，确认：

1. **列名已正确**：当前格式 `Dock Index, Enzyme Index, Substrate Index, Label` 与代码期望一致
2. **列顺序不重要**：代码通过列名访问，非位置访问
3. **`positive_reactions`列不需要**：代码中未使用该字段
4. **数据类型已正确**：所有索引列已是int64类型

### 2.3 格式对比

| 格式类型 | 列名 | 说明 |
|----------|------|------|
| ESIBank原始 | `enzyme, reaction, label, structure_index` | 原始训练数据格式 |
| 代码期望 | `Enzyme Index, Substrate Index, Label, Dock Index` | 模型加载格式 |
| 我们的数据 | `Dock Index, Enzyme Index, Substrate Index, Label` | 已兼容 |

---

## 三、修正方案

### 3.1 最终方案

仅需移除UTF-8 BOM，其他格式已正确。

### 3.2 执行脚本

位置：`scripts/04_Step4_CSV格式修正/step4_csv_format_fix.py`

功能：
- 读取源CSV（pandas自动处理BOM）
- 验证必需列存在
- 确保索引列为int64类型
- 写出无BOM的UTF-8编码文件

---

## 四、执行结果

### 4.1 处理统计

| 数据集 | 行数 | 原BOM | 修正后BOM | 列验证 | 类型验证 |
|--------|------|-------|-----------|--------|----------|
| B1_仅底物_272pos | 272 | True | False | OK | OK |
| B2_底物正产物负_272pos23neg | 295 | True | False | OK | OK |
| B3_完整数据集_272pos267neg | 539 | True | False | OK | OK |
| B4_仅产物_23neg | 23 | True | False | OK | OK |
| B5_仅抑制剂_244neg | 244 | True | False | OK | OK |
| B6_底物正抑制剂负_272pos244neg | 516 | True | False | OK | OK |
| Substrates.csv | 436 | True | False | - | - |
| Enzymes.csv | 292 | True | False | - | - |

### 4.2 Codex审核结论

> "I inspected the actual outputs... All `B1..B6/*/data.csv` headers are exactly: `Dock Index,Enzyme Index,Substrate Index,Label` (no BOM; first bytes are not `EF BB BF`). All four columns parse as `int64` and contain no NaNs. Index ranges are valid vs the copied reference tables."

**审核状态**：✅ PASSED

---

## 五、输出文件

### 5.1 输出目录结构

```
data/04_Step4_格式修正后数据/
├── B1_仅底物_272pos/
│   └── data.csv
├── B2_底物正产物负_272pos23neg/
│   └── data.csv
├── B3_完整数据集_272pos267neg/
│   └── data.csv
├── B4_仅产物_23neg/
│   └── data.csv
├── B5_仅抑制剂_244neg/
│   └── data.csv
├── B6_底物正抑制剂负_272pos244neg/
│   └── data.csv
├── Enzymes.csv
├── Substrates.csv
└── format_fix_report.txt
```

### 5.2 data.csv格式

| 列名 | 类型 | 说明 |
|------|------|------|
| Dock Index | int64 | 对接复合物索引（对应PDB文件） |
| Enzyme Index | int64 | 酶索引（指向Enzymes.csv） |
| Substrate Index | int64 | 底物索引（指向Substrates.csv） |
| Label | int64 | 标签（1=正样本, 0=负样本） |

---

## 六、验证检查

### 6.1 格式验证

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| BOM移除 | False | False | ✅ |
| 列名正确 | 4列存在 | 4列存在 | ✅ |
| 类型正确 | int64 | int64 | ✅ |

### 6.2 数据完整性

| 检查项 | 预期 | 实际 | 状态 |
|--------|------|------|------|
| Enzyme Index范围 | [0, 291] | [0, 291] | ✅ |
| Substrate Index范围 | [0, 435] | [0, 435] | ✅ |
| 无NaN值 | 是 | 是 | ✅ |

---

## 七、下一步

Step 5: 生成ESM蛋白质嵌入

**输入**：`data/04_Step4_格式修正后数据/Enzymes.csv`
**输出**：`enzyme_features.lmdb`
**关键列**：`Protein sequence`, `Enzyme_Index`

---

## 八、关键决策记录

### 决策1：不添加positive_reactions列

| 选项 | 描述 | 采用 |
|------|------|------|
| 添加positive_reactions | 与toy_example格式一致 | ❌ |
| 不添加 | 代码中未使用该字段 | ✅ |

**理由**：Codex和Gemini均确认代码中不读取此列。

### 决策2：保持列顺序不变

| 选项 | 描述 | 采用 |
|------|------|------|
| 重排列顺序 | 与toy_example完全一致 | ❌ |
| 保持当前顺序 | 代码按名访问，顺序无关 | ✅ |

**理由**：代码使用`df['Enzyme Index']`等方式访问，不依赖列位置。

---

## 九、完整代码库分析参考

> **重要**: 关于CSV格式的详细分析、代码证据、LMDB必要性等，请参阅：
> `reports/数据集格式与负样本生成完整分析报告.md` 中的 **Q6: 完整代码库分析** 章节

### 9.1 关键代码证据摘要

#### data.csv列名要求

**来源**: [src/Datasets/data_representer.py:45-47](../../src/Datasets/data_representer.py#L45-L47)

```python
reation_idx = self.df.loc[idx, 'Substrate Index']  # 必须精确匹配
enzyme_idx = self.df.loc[idx, 'Enzyme Index']      # 必须精确匹配
```

**来源**: [src/Datasets/Structure/structure.py:61-65](../../src/Datasets/Structure/structure.py#L61-L65)

```python
for index, (reaction, protein, label, structure_index) in enumerate(zip(
    df['Substrate Index'].values,   # 不是 'reaction'
    df['Enzyme Index'].values,      # 不是 'enzyme'
    df["Label"].values,             # 不是 'label'（大写L）
    df['Dock Index'].values         # 不是 'structure_index'
)):
```

### 9.2 数据流图

```
data.csv (Enzyme Index, Substrate Index, Dock Index, Label)
    ↓
代码读取Index值
    ↓
用Index作为key查询LMDB
    ↓
enzyme_features.lmdb[str(Enzyme Index)] → ESM嵌入
reaction_features.lmdb[str(Substrate Index)] → 底物图特征
grover_fingerprint.lmdb[str(Substrate Index)] → GROVER嵌入
morgan_fingerprint.npy[Substrate Index] → Morgan指纹
structure_features.lmdb[str(Dock Index)] → 结构特征
    ↓
模型推理
    ↓
输出score列
```

### 9.3 Index含义

| Index列 | 含义 | 对应文件 |
|---------|------|----------|
| Enzyme Index = N | Enzymes.csv第N行（0-based） | enzyme_features.lmdb key="N" |
| Substrate Index = M | Substrates.csv第M行 | reaction_features.lmdb key="M" |
| Dock Index = K | 对接结构K | pocket/K.pdb + ligand/K.sdf |

---

**文档版本**: v1.0
**最后更新**: 2026-01-29
