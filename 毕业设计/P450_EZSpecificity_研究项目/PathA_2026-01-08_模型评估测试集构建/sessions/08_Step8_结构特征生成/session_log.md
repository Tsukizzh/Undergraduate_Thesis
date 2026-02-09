# Step 8: 结构特征生成 - 详细操作日志

> **执行时间**: 2026-01-30 ~ 2026-01-31
> **版本**: v3.1 (修复RDKit flavor参数后)
> **执行者**: Claude Code + Codex + Gemini (三方辩论审核)

---

## 一句话总结

**从 PDB 文件中提取蛋白质口袋（10Å范围）和配体的 3D 坐标，生成 StructureComplexData 存入 LMDB，供 EGNN 图神经网络使用。**

---

## ⚠️ 关键问题：为什么 539 条只剩 517 条？

### 数据损耗全景（v3.1 修复后）

```
B3 数据集: 539 条
    ↓
口袋/配体提取: 532 成功，7 失败 (损耗 1.3%)
    ↓
配体对齐: 517 成功，15 失败 (损耗 2.8%)
    ↓
最终: 517 条 (总损耗 4.1%)
```

### 失败原因分类（v3.1 修复后基于 CSV 实际数据）

| 类型 | 数量 | 真正原因 | 可修复？ |
|------|------|---------|---------|
| **Raw ligand file not found** | 8 | CIF 格式配体标识符不兼容 RDKit（7条CIF+1条提取失败） | 需 CCD 方案 |
| **AssignBondOrdersFromTemplate failed** | 9 | PDB 3D 构象与 SMILES 拓扑不匹配 | 部分可改善 |
| **RemoveHs failed (valence issues)** | 3 | 化学价异常（Cl=2, C=5, O=3） | 需手工处理 |
| **No substructure match found** | 2 | 提取的配体结构与 SMILES 不一致 | 部分可改善 |

### v3.0→v3.1 修复了什么？

**v3.1 发现并修复了 RDKit flavor 参数问题**：

```python
# v3.0 (bug): RDKit 默认把 altLoc (column 17) 和残基名 (columns 18-20) 拼接
# 例如配体 5LW 的 altLoc=B，RDKit 读成残基名 "B5LW" 而非 "5LW"
mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)

# v3.1 (修复): 使用 flavor=1 让 RDKit 从 columns 77-78 读取元素符号
mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False, flavor=1)
```

**修复效果**：
- 对齐成功率：511/539 → **517/539**（恢复 6 条）
- 恢复的配体：5EAH/5LW, 5EAH/5LX, 5EAH/5LY, 6CPP/CAE, 2P85/IND, 6WGW/CAM
- 最终数量从 511 增至 **517**（+6 条，每条数据都很珍贵！）

### 这是前面步骤的问题吗？

**不是。** 前面步骤（Step 1-7）均正常。22 条失败的原因是：
- 8 条：CIF 文件的 RDKit 解析兼容性问题（工程限制，需 CCD 连接性方案）
- 14 条：PDB 实验结构与 SMILES 的化学表示差异（互变异构体/质子化/键序等）

### 会影响后续推理吗？

**不会**。`high_quality_id.txt` 包含 517 个有效 ID，下游推理自动过滤：
```python
valid_ids = set(open('high_quality_id.txt').read().split())
df = df[df['Dock Index'].isin(valid_ids)]  # 539 → 517
```

---

## 为什么需要这一步？

EZSpecificity 使用 **SE(3)-等变图神经网络（EGNN）** 处理酶-底物复合物的 3D 结构。它需要：
- 蛋白质口袋的原子坐标（配体周围 10Å 范围内的蛋白质原子）
- 配体的原子坐标
- 原子类型、键连接等

---

## 输入数据

### 1. PDB 结构文件

| 属性 | 值 |
|------|-----|
| 位置 | `data/01_Step1_PDB文件/` |
| 数量 | 627 个（对应 539 条 B3 记录，有些 PDB 有多个配体） |
| 格式 | .pdb (619个) + .cif (8个) |

### 2. dock_index_mapping.csv

| 属性 | 值 |
|------|-----|
| 位置 | `data/04_Step4_格式修正后数据/dock_index_mapping.csv` |
| 内容 | Dock_Index → (pdb_id, ligand_ccd, Substrate_Index) 映射 |

### 3. Substrates.csv

| 属性 | 值 |
|------|-----|
| 位置 | `data/04_Step4_格式修正后数据/Substrates.csv` |
| 内容 | Substrate_Index → Substrate_SMILES（用于原子序对齐） |

---

## 输出数据

### 主输出：structure_features.lmdb

每条记录是一个 `StructureComplexData` 对象：

| 字段 | 形状 | 含义 |
|------|------|------|
| protein_pos | (N, 3) | 口袋原子的 xyz 坐标 |
| protein_element | (N,) | 口袋原子的元素（原子序数） |
| ligand_pos | (M, 3) | 配体原子的 xyz 坐标 |
| ligand_element | (M,) | 配体原子的元素 |
| **ligand_index** | (M,) | **关键字段**：PDB原子→SMILES原子的映射 |

### ligand_index 为什么重要？

下游模型代码（[src/Models/ss.py:90](src/Models/ss.py#L90)）：
```python
graph_output = graph_output[batch.ligand_index]  # 用 ligand_index 重排
```

如果 ligand_index 全是 0，模型会把所有原子都当成第 0 个原子，**预测必然错误**！

---

## 三步处理流程

```
步骤 8.1: 口袋/配体提取
  输入: 627 个 PDB 文件 + dock_index_mapping.csv
  工具: BioPython NeighborSearch
  输出: str_tmp_data/pocket/*.pdb (539个, 含失败的空文件)
        str_tmp_data/raw_ligand/*.sdf (532个)
  结果: 532/539 成功 (98.7%)

步骤 8.1b: 配体 SMILES 对齐 (关键！)
  输入: raw_ligand/*.sdf + Substrates.csv
  工具: RDKit AssignBondOrdersFromTemplate + GetSubstructMatches
  输出: str_tmp_data/ligand/*.sdf (511个，含正确的 AtomMapNum)
  结果: 511/539 成功 (94.8%)

步骤 8.2: LMDB 生成
  输入: pocket/*.pdb + ligand/*.sdf
  输出: structure_features.lmdb (511条)
        high_quality_id.txt (511个ID)
  结果: 511/511 成功 (100%)
```

---

## 失败原因详解

### 阶段 1：口袋/配体提取失败 (7 条)

全部为 CIF 文件，RDKit 无法解析 PDBIO 输出的 PDB block：
- Dock 26 (9KPU/A1EGF)
- Dock 192 (8VK6/A1AB6)
- Dock 277 (9KW3/A1L6Q)
- Dock 278 (9KW4/A1L6X)
- Dock 294 (9MS2/A1BNX)
- Dock 335 (8S4M/A1H47)
- Dock 534 (9CV8/A1AZ7)

### 阶段 2：配体对齐失败 (21 条)

| 原因 | 数量 | 说明 |
|------|------|------|
| AssignBondOrdersFromTemplate failed | 15 | PDB 键序与 SMILES 拓扑不匹配 |
| No substructure match found | 2 | 提取的配体结构与 SMILES 不一致 |
| RemoveHs failed (valence) | 3 | RDKit 检测到无效化合价 |
| AssignBondOrdersFromTemplate (valence) | 1 | 模板匹配时化合价异常 |

---

## 遇到的 BUG 及修复

### BUG 1: ligand_index 全零（最严重）

**现象**：第一次运行后，LMDB 里的 ligand_index 全是 `[0, 0, 0, ...]`

**根因**：PDB 原子顺序是晶体学实验顺序，和 SMILES 不一致；`MolFromPDBBlock` 不会设置 `AtomMapNum`

**修复**：新增 `step8_align_ligand.py`，使用 RDKit 模板匹配对齐原子序

### BUG 2: altLoc 过滤过严（v3.0 修复）

**现象**：13 条提取失败，其中 6 条是因为配体只有 altLoc B/C/D

**根因**：`accept_atom()` 只接受 `""` 或 `"A"`，丢弃了 B/C/D 配体；mmCIF 的 `label_alt_id='.'` 也被当作非空 altLoc 拒绝

**修复**：
- 添加 `_normalize_altloc()` 处理 PDB/mmCIF 差异（`.`/`?` 视为空白）
- 添加 `_preferred_residue_altloc()` 选择最常见 altLoc（投票排除空白原子，tie-break 优先 'A'）
- 三方审核（Claude+Codex+Gemini）确认逻辑正确

### BUG 3: Pickle 序列化失败

**修复**：注入假模块到 `sys.modules`

### BUG 4: UTF-8 BOM

**修复**：`encoding='utf-8-sig'`

### BUG 5: Valence 异常

**修复**：`try-except` 包裹 `RemoveHs`

---

## 验证结果

| 检查项 | 结果 |
|--------|------|
| LMDB 记录数 | 511 条 |
| ligand_index 不全零 | ✅ 验证通过 |
| ligand_index 范围正确 | ✅ [0, N-1] |
| protein_pos 非空 | ✅ |
| ligand_pos 非空 | ✅ |
| Pickle 反序列化成功 | ✅ |

---

## 文件位置

```
data/08_Step8_结构特征/
├── str_tmp_data/
│   ├── pocket/*.pdb          ← 539 个（含失败的空文件）
│   ├── raw_ligand/*.sdf      ← 532 个原始配体
│   └── ligand/*.sdf          ← 511 个对齐配体（ligand_index 正确）
├── structure_features.lmdb   ← 主输出（511 条）
└── high_quality_id.txt       ← 有效 ID 列表（511 个）

scripts/08_Step8_结构特征生成/
├── step8_extract_pocket_ligand.py   ← 步骤 8.1（v3.0 修复 altLoc bug）
├── step8_align_ligand.py            ← 步骤 8.1b
└── step8_generate_structure_lmdb.py ← 步骤 8.2
```

---

## 常见问题

### Q: 为什么用 B3 数据集（539条）而不是 B1（271条）？

A: B3 是最大的数据集（包含 substrate + product + inhibitor），结构特征是共享的。其他数据集（B1/B2/B4/B5/B6）是 B3 的子集，使用同一份 structure_features.lmdb，按 high_quality_id.txt 过滤即可。

### Q: 28 条失败的数据怎么办？

A: 这些数据在下游推理时会被自动跳过（因为不在 high_quality_id.txt 中）。不影响其他 511 条数据的推理结果。后续版本可通过 CCD 连接性方案和更鲁棒的模板匹配策略召回部分数据。

### Q: 下游推理需要做什么特殊处理？

A: 需要在反序列化 LMDB 前注入假模块：
```python
import types, sys
sys.modules['Datasets.Structure.utils'] = types.ModuleType('...')
# 然后才能 pickle.loads
```

### Q: 为什么 Step 6/7 是 436 条，Step 8 是 511 条？

A: 不同的索引系统：
- Step 6/7 按 **Substrate_Index** 索引（436 个唯一底物）
- Step 8 按 **Dock_Index** 索引（539 个酶-底物对，最终 511 个有效）

一个底物可能和多个酶形成复合物，所以 Dock 数量 > Substrate 数量。

---

## 下一步

Step 8 完成后，所有特征都已准备好。接下来：
- **Step 9**: 模型推理 - 在 B1-B6 数据集上运行 EZSpecificity 预测
- **Step 10**: 结果分析 - 计算 AUC/F1 等指标

---

**文档版本**: v3.0 (修复 altLoc bug，三方审核)
**最后更新**: 2026-01-31
