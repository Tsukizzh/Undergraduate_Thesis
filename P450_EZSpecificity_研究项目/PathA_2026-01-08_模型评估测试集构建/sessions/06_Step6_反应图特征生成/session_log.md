# Step 6: 反应图特征生成 - 详细操作日志

> **执行时间**: 2026-01-30
> **版本**: v1.0 (新手友好版)
> **执行者**: Claude Code + Codex (三轮讨论)

---

## 一句话总结

**把底物分子的"文字描述"(SMILES) 转换成"图结构数据"，让图神经网络能够理解分子的形状和连接方式。**

---

## 为什么需要这一步？

### 问题：模型看不懂SMILES文字

我们的底物数据长这样（人类化学家能看懂）：
```
COc1ccccc1O
```
这叫**SMILES**，是一种用文字描述分子结构的方法。但神经网络**看不懂文字**，它需要数字化的结构信息。

### 解决方案：把SMILES转成"分子图"

**分子图**(Molecular Graph) 是数学中"图论"的概念：
- **节点(Node)** = 原子（碳、氧、氮...）
- **边(Edge)** = 化学键（单键、双键、芳香键...）

这就像把一段话（SMILES）翻译成一张地图（Graph），地图上标明了每个原子在哪里、跟谁相连。

### 简单类比

```
SMILES文字:      COc1ccccc1O    (像一段暗语)
       ↓ RDKit解析
分子图结构:      C─O─C═C─C═C─C═C─O    (像一张地图)
                 节点+边
       ↓ 提取数字
图数据:          element=[6,8,6,6,6,6,6,6,8]    (数字化的地图)
                 edge_index=[[0,1],[1,0],[1,2],...]
```

---

## 输入数据长什么样？

### 文件：`Substrates.csv`

| 字段 | 含义 | 示例 |
|------|------|------|
| Substrate_Index | 底物编号 | 0, 1, 2, ... |
| Substrate_SMILES | 分子的SMILES表达 | COc1ccccc1O |

### 具体示例（前5条）

**第0号底物：** 邻甲氧基苯酚
```
Substrate_Index: 0
Substrate_SMILES: COc1ccccc1O
含义: C(碳)-O(氧)-苯环-O(氧)
原子数: 9个
```

**第1号底物：**
```
Substrate_Index: 1
Substrate_SMILES: COc1cccc(OC)c1O
原子数: 11个
```

**第2号底物：** 一个含有色氨酸结构的大分子
```
Substrate_Index: 2
Substrate_SMILES: CN[C@H](C(=O)N[C@H](CO)Cc1c[nH]c2ccccc12)C(C)C
原子数: 22个
```

**总计：436条底物记录**

---

## 输出数据长什么样？

### 文件：`reaction_features.lmdb`

LMDB是一种数据库格式，每条记录包含：

| 字段 | 含义 | 形状 |
|------|------|------|
| element | 每个原子的元素序号 | (原子数,) |
| edge_index | 哪些原子之间有化学键 | (2, 键数×2) |
| edge_type | 化学键的类型编号 | (键数×2,) |
| atom_feature | 原子的7维特征 | (原子数, 7) |
| reaction_attention_label | 反应注意力标签(全0) | (原子数,) |
| num_nodes | 原子总数 | 标量 |

### 具体示例（前3条）

**Key = "0"（第0号底物 `COc1ccccc1O`，邻甲氧基苯酚）：**
```
element: [6, 8, 6, 6, 6, 6, 6, 6, 8]
  → 9个原子：C(6), O(8), C(6), C(6), C(6), C(6), C(6), C(6), O(8)

edge_index:
  起点: [0, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7, 8]
  终点: [1, 0, 2, 1, 3, 7, 2, 4, 3, 5, 4, 6, 5, 7, 2, 6, 8, 7]
  → 9条化学键（每条键正反两个方向 = 18条有向边）

edge_type: [1, 1, 1, 1, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 1, 1]
  → 1=单键, 12=芳香键

atom_feature (每个原子的7维特征):
  第0个原子(C): [6, 0, 1, 0, 0, 0, 1]
  → [原子序数=6, 是否芳香=否, 连接度=1, 相邻氢数=0, SP=否, SP2=否, SP3=是]

reaction_attention_label: [0, 0, 0, 0, 0, 0, 0, 0, 0]
  → 全0（substrate-only模式）

num_nodes: 9
```

**Key = "1"（第1号底物 `COc1cccc(OC)c1O`）：**
```
element: [6, 8, 6, 6, 6, 6, 6, 8, 6, 6, 8]
  → 11个原子

edge_index: shape=(2, 22)  → 11条键

atom_feature: shape=(11, 7)
```

**Key = "2"（第2号底物，含色氨酸结构）：**
```
element: [6, 7, 6, 6, 8, 7, 6, 6, 8, 6, 6, 6, 7, 6, 6, 6, 6, 6, 6, 6, 6, 6]
  → 22个原子（含氮N=7、氧O=8）

edge_index: shape=(2, 46)  → 23条键

atom_feature: shape=(22, 7)
```

**总计：436条记录，与输入完全对应**

---

## 元素序号对照表

SMILES中的原子字母和`element`数组中的数字对应关系：

```
H  → 1    (氢)         N  → 7    (氮)         Cl → 17   (氯)
B  → 5    (硼)         O  → 8    (氧)         Br → 35   (溴)
C  → 6    (碳)         F  → 9    (氟)         I  → 53   (碘)
                       S  → 16   (硫)         P  → 15   (磷)
```

---

## atom_feature 7维特征详解

每个原子用7个数字描述：

| 维度 | 含义 | 示例值 | 说明 |
|------|------|--------|------|
| [0] | 原子序数 | 6 | 碳=6, 氧=8, 氮=7 |
| [1] | 是否芳香 | 0或1 | 苯环上的原子=1 |
| [2] | 连接度(degree) | 1~4 | 与几个其他原子相连 |
| [3] | 相邻氢原子数 | 0~3 | 连着几个氢 |
| [4] | SP杂化 | 0或1 | 三键碳 |
| [5] | SP2杂化 | 0或1 | 双键/芳香碳 |
| [6] | SP3杂化 | 0或1 | 单键碳 |

### 转换示例

以 `COc1ccccc1O` 的第0个原子（甲基碳C）为例：
```
atom_feature[0] = [6, 0, 1, 0, 0, 0, 1]

含义:
  [6]  原子序数 = 6（碳原子）
  [0]  不在芳香环上
  [1]  连接度 = 1（只连着一个氧原子）
  [0]  没有相邻氢原子（在显式表示中）
  [0]  不是SP杂化
  [0]  不是SP2杂化
  [1]  是SP3杂化（甲基碳）
```

---

## 图解：整体转换流程

```
输入 Substrates.csv                输出 reaction_features.lmdb
┌─────────────────────┐            ┌──────────────────────────────────┐
│ Substrate_Index: 0  │            │ Key: "0"                         │
│ SMILES:             │   RDKit    │ element: [6,8,6,6,6,6,6,6,8]    │
│ COc1ccccc1O         │ ────────→  │   每个原子 → 元素序号            │
│ (9个原子)           │            │ edge_index: (2, 18)              │
│                     │            │   哪些原子之间有键               │
│                     │            │ atom_feature: (9, 7)             │
│                     │            │   每个原子 → 7维特征             │
└─────────────────────┘            └──────────────────────────────────┘

         ↓ 重复436次 ↓

┌─────────────────────┐            ┌──────────────────────────────────┐
│ Substrate_Index: 435│            │ Key: "435"                       │
│ SMILES:             │   RDKit    │ element: [6,6,6,6,6,8,8,...]    │
│ CC[C@H](C)C(=O)... │ ────────→  │ edge_index: (2, 60)              │
│ (28个原子)          │            │ atom_feature: (28, 7)            │
└─────────────────────┘            └──────────────────────────────────┘
```

---

## 所需工具与环境

| 工具 | 版本 | 作用 |
|------|------|------|
| Python | 3.12.5 | 运行环境 |
| RDKit | 2025.03.6 | 解析SMILES，提取原子和键 |
| LMDB | 0.9.33 | 存储最终的图数据 |
| NumPy | 1.26.4 | 数组处理 |
| Pandas | 2.2.2 | 读取CSV |

**注意**: 本脚本使用了自定义的 `parse_smile_no_scatter()` 函数，无需安装 `torch_scatter`。

---

## 三方讨论过程

| 轮次 | 参与者 | 讨论内容 | 结论 |
|------|--------|----------|------|
| 1 | Codex | 实施方案设计、边缘情况处理 | 单进程、LMDB key=Substrate_Index |
| 2 | Codex | 文件夹命名、输入路径、map_size | 中文命名、默认Step4输出、4GB |
| 3 | Codex | 兼容性验证、集成测试方案 | 输出格式与brenda.py完全匹配 |
| 审核 | Codex | 代码审核 | 通过，num_hs计算等价 |

### 发现并解决的问题

| 问题 | 说明 | 修复 |
|------|------|------|
| torch_scatter缺失 | torch环境未安装torch_scatter | 自定义get_ligand_atom_features实现 |
| SMILES映射为空 | 简单分子的_smilesAtomOutputOrder为空 | 添加空值回退处理 |

---

## 实际执行过程

### 执行命令
```bash
D:\anaconda3\envs\torch\python.exe step6_generate_reaction_graph_features.py --validate
```

### 第一次执行（失败 - 47个分子报错）
```
Stats: total=436 ok=389 skipped=0 failed=47
原因: parse_smile failed: ValueError: invalid literal for int()
问题: 简单分子（如苯甲基 Cc1ccccc1）的 _smilesAtomOutputOrder 属性为空
```

### 修复后第二次执行（成功）
```
Stats: total=436 ok=436 skipped=0 failed=0
Validation: 436 / 436 keys present.
处理速度: ~200 it/s
```

### 硬件使用
- CPU: 单进程运行
- GPU: 未使用（此步骤不需要GPU）
- 处理时间: 约2秒

---

## 验证结果

| 检查项 | 结果 |
|--------|------|
| LMDB记录数量 | 436条 |
| Key完整性 | 0-435无缺失 |
| element形状 | 全部为 (n,) |
| edge_index形状 | 全部为 (2, E) |
| atom_feature形状 | 全部为 (n, 7) |
| reaction_attention_label | 全部为零向量 |

### 抽样验证

| Key | SMILES片段 | 原子数 | 键数 | atom_feature |
|-----|-----------|--------|------|-------------|
| 0 | COc1ccccc1O | 9 | 18 | (9, 7) |
| 1 | COc1cccc(OC)c1O | 11 | 22 | (11, 7) |
| 2 | CN[C@H](C(=O)... | 22 | 46 | (22, 7) |
| 50 | O=C1N[C@@H](Cc... | 21 | 48 | (21, 7) |
| 100 | CC(C)(C)[C@H](O... | 20 | 42 | (20, 7) |
| 200 | N[C@@H]1CCCc2cc... | 11 | 24 | (11, 7) |
| 435 | CC[C@H](C)C(=O)... | 28 | 60 | (28, 7) |

---

## 文件位置

```
PathA_2026-01-08_模型评估测试集构建/
├── data/
│   ├── 04_Step4_格式修正后数据/
│   │   └── Substrates.csv                    ← 输入 (436条底物)
│   └── 06_Step6_反应图特征/
│       ├── reaction_features.lmdb            ← 主要输出 (436条图特征)
│       └── reaction_features_report.csv      ← 处理报告
└── scripts/
    └── 06_Step6_反应图特征生成/
        └── step6_generate_reaction_graph_features.py  ← 执行脚本
```

---

## 常见问题

### Q: 为什么不直接用SMILES字符串？
A: 因为图神经网络需要知道原子之间的连接关系（拓扑结构），SMILES只是一串文字，模型无法从中提取结构信息。

### Q: "图"和"图片"有什么区别？
A: 完全不同。这里的"图"是数学中的Graph（节点+边），不是Image（像素矩阵）。

### Q: edge_type里的12是什么意思？
A: 12 = 芳香键(Aromatic Bond)。1 = 单键(Single Bond)，2 = 双键(Double Bond)，3 = 三键(Triple Bond)。

### Q: reaction_attention_label为什么全是0？
A: 因为我们使用的是"substrate-only"模式（只有底物，没有产物信息）。如果有完整的反应方程式，这个标签会标记出反应中发生变化的原子。

### Q: 为什么不需要GPU？
A: 这一步只用RDKit做分子解析，是纯CPU计算，不涉及深度学习模型。与Step 5的ESM嵌入不同。

---

## 下一步

Step 6完成后，436条底物的图特征已经准备好。接下来需要：
- **Step 7**: 生成Morgan指纹（`morgan_fingerprint.npy`，1024维二进制向量）
- **Step 8**: 生成GROVER分子嵌入（`grover_fingerprint.lmdb`）

---

**文档版本**: v1.0 (新手友好版 + 真实数据展示)
**最后更新**: 2026-01-30
