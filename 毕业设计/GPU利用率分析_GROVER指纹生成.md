# GPU 利用率分析：GROVER 指纹生成中的 CPU-GPU 瓶颈问题

> **观测日期**: 2026-01-30
> **硬件环境**: NVIDIA GeForce RTX 4070 Super / 32 GB RAM / Windows
> **观测现象**: GROVER 指纹生成过程中 GPU 利用率仅约 1%
> **处理规模**: 436 个化学底物，总耗时 54.5 秒

---

## 一、现象描述

在使用 GROVER-Large 预训练模型为 436 个底物生成分子指纹时，通过 `nvidia-smi` 监控发现 GPU 利用率始终徘徊在 0-1%。直觉上，既然代码中明确使用了 CUDA（`args.cuda=True`），且显卡是性能强劲的 RTX 4070 Super，GPU 利用率应该较高才对。

**但事实是**：这并非异常，而是该数据管线架构下的**正常行为**。本文将从代码级别深入剖析其原因。

---

## 二、核心类比：理解 CPU-GPU 协作模式

为了帮助理解，我们用一个贯穿全文的类比：

> **CPU 是备菜大厨，GPU 是光速烤箱。**
>
> - **CPU（备菜大厨）** 负责读懂化学食谱（SMILES）、洗菜切菜（解析分子结构）、按复杂规则摆盘（提取特征、构建图数据结构）。这些是精细的手工活，需要逐步判断、逐项填写。
> - **GPU（光速烤箱）** 负责烤披萨（矩阵乘法/Transformer推理）。它一次能放进32个披萨，按下按钮后 0.1 秒就全部烤熟。
> - **DataLoader（传送带）** 负责把大厨做好的生披萨送到烤箱。
>
> 现在的状况：大厨花了53秒切菜备料，烤箱只在最后1.5秒"轰"了几下就完成了全部烘烤。**不是烤箱不努力，是菜备得太慢。**

---

## 三、数据管线全景图

GROVER 指纹生成的完整数据管线如下：

```
┌─────────────────────────────────────────────────────────────────┐
│                        CPU 领域（占 97% 时间）                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Phase 1: 数据加载                                               │
│  ┌─────────────────────────────────────────────────────┐        │
│  │ get_data() → 读取 CSV → 创建 MoleculeDatapoint 列表  │        │
│  └─────────────────────────────────────────────────────┘        │
│                          ↓                                      │
│  Phase 2: DataLoader 分发任务                                    │
│  ┌─────────────────────────────────────────────────────┐        │
│  │ DataLoader(batch_size=32, num_workers=4)             │        │
│  │   → 将 436 个分子分成 14 个 batch                      │        │
│  │   → 每个 batch 交给 MolCollator 处理                   │        │
│  └─────────────────────────────────────────────────────┘        │
│                          ↓                                      │
│  Phase 3: MolCollator（瓶颈所在！）                               │
│  ┌─────────────────────────────────────────────────────┐        │
│  │ 对 batch 中的每个 SMILES:                              │        │
│  │   ├─ Chem.MolFromSmiles()     → 解析分子              │        │
│  │   ├─ 4x GetSubstructMatches() → 化学模式匹配           │        │
│  │   ├─ GetRingInfo()            → 环结构检测             │        │
│  │   ├─ atom_features()          → 151维原子特征          │        │
│  │   ├─ bond_features()          → 14维键特征             │        │
│  │   └─ 邻接关系索引构建                                   │        │
│  │ 最后: BatchMolGraph()                                  │        │
│  │   → 拼接所有分子图 → 转为 PyTorch 张量                  │        │
│  └─────────────────────────────────────────────────────┘        │
│                                                                 │
├─────────────────────────────────────────────────────────────────┤
│                        GPU 领域（占 3% 时间）                    │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  Phase 4: 模型前向推理                                           │
│  ┌─────────────────────────────────────────────────────┐        │
│  │ GroverFpGeneration.forward()                         │        │
│  │   ├─ GROVEREmbedding → GTransEncoder                 │        │
│  │   │   └─ 多头注意力 + 前馈网络（纯矩阵运算）            │        │
│  │   ├─ Readout（均值池化）                               │        │
│  │   └─ 特征拼接 → 最终指纹                               │        │
│  └─────────────────────────────────────────────────────┘        │
│                          ↓                                      │
│  Phase 5: 写入 LMDB（回到 CPU）                                  │
│  ┌─────────────────────────────────────────────────────┐        │
│  │ pickle.dumps(data) → txn.put()                       │        │
│  └─────────────────────────────────────────────────────┘        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 四、逐层深度剖析

### 4.1 Phase 1: 数据加载 — `get_data()`

**代码位置**: [fingerprint.py:82-87](src/other_softwares/grover_software/task/fingerprint.py#L82-L87)

```python
test_data = get_data(path=args.data_path, args=args,
                     use_compound_names=False,
                     max_data_size=float("inf"),
                     skip_invalid_smiles=False)
test_data = MoleculeDataset(test_data)
```

**做了什么**：
1. 读取 CSV 文件（`grover_substrates.csv`，只有一列 `smiles`）
2. 为每一行创建一个 `MoleculeDatapoint` 对象（[moldataset.py:18-50](src/other_softwares/grover_software/grover/data/moldataset.py#L18-L50)），存储 SMILES 字符串和预计算的特征数组
3. 用 `MoleculeDataset` 包装成 PyTorch Dataset

**类比**：这一步相当于**把菜单上的 436 道菜名抄到小纸条上**。很快，不是瓶颈。

**耗时估计**：< 0.5 秒

---

### 4.2 Phase 2: DataLoader 配置

**代码位置**: [fingerprint.py:35-42](src/other_softwares/grover_software/task/fingerprint.py#L35-L42)

```python
mol_collator = MolCollator(args=args, shared_dict={})

num_workers = 4
mol_loader = DataLoader(data,
                        batch_size=32,
                        shuffle=False,
                        num_workers=num_workers,
                        collate_fn=mol_collator)
```

**关键参数解释**：

| 参数 | 值 | 含义 |
|------|-----|------|
| `batch_size` | 32 | 每次处理 32 个分子（一炉 32 个披萨） |
| `num_workers` | 4 | 4 个子进程并行预处理数据（4 个厨师） |
| `collate_fn` | MolCollator | 自定义的数据整理函数（大厨的工作流程） |
| `shared_dict` | `{}` | 空字典，用于缓存已解析的分子图 |

**批次划分**：436 个分子 / 32 = **14 个批次**（最后一批不足 32 个）

**类比**：这一步相当于**安排 4 个厨师，每人轮流处理一批 32 个订单**。

---

### 4.3 Phase 3: MolCollator — CPU 瓶颈的核心

**代码位置**: [molgraph.py:369-387](src/other_softwares/grover_software/grover/data/molgraph.py#L369-L387)

```python
class MolCollator(object):
    def __init__(self, shared_dict, args):
        self.args = args
        self.shared_dict = shared_dict

    def __call__(self, batch):
        smiles_batch = [d.smiles for d in batch]          # 提取 SMILES
        features_batch = [d.features for d in batch]      # 提取预计算特征
        target_batch = [d.targets for d in batch]         # 提取目标值
        batch_mol_graph = mol2graph(smiles_batch,         # ← 瓶颈！
                                    self.shared_dict,
                                    self.args)
        batch = batch_mol_graph.get_components()          # 转为张量元组
        mask = torch.Tensor([[...] for tb in target_batch])
        targets = torch.Tensor([[...] for tb in target_batch])
        return smiles_batch, batch, features_batch, mask, targets
```

**每当 DataLoader 需要一个 batch 时**，就会调用这个 `__call__` 方法。关键瓶颈在第 383 行的 `mol2graph()` 调用。

**类比**：每次烤箱需要一炉披萨时，大厨就要从头开始：读菜单 → 洗菜 → 切菜 → 摆盘 → 装盘。烤箱在这整个过程中都在等。

---

#### 4.3.1 `mol2graph()` — 分子图构建

**代码位置**: [molgraph.py:347-366](src/other_softwares/grover_software/grover/data/molgraph.py#L347-L366)

```python
def mol2graph(smiles_batch, shared_dict, args):
    mol_graphs = []
    for smiles in smiles_batch:          # 逐个处理 batch 中的 SMILES
        if smiles in shared_dict:        # 检查缓存
            mol_graph = shared_dict[smiles]
        else:
            mol_graph = MolGraph(smiles, args)    # ← 最耗时的操作
            if not args.no_cache:
                shared_dict[smiles] = mol_graph   # 缓存结果
        mol_graphs.append(mol_graph)

    return BatchMolGraph(mol_graphs, args)        # 合并为批次图
```

这个函数对 batch 中的每个 SMILES 字符串：
1. 先查缓存（`shared_dict`）
2. 如果没命中，调用 `MolGraph()` 从头构建分子图
3. 最后用 `BatchMolGraph()` 将所有图合并

---

#### 4.3.2 `MolGraph.__init__()` — 单分子图构建（最耗时）

**代码位置**: [molgraph.py:80-230](src/other_softwares/grover_software/grover/data/molgraph.py#L80-L230)

这是整个管线中**最耗 CPU 的部分**。每个分子要经历以下 6 个步骤：

##### 步骤 A: SMILES 解析（第 113 行）

```python
mol = Chem.MolFromSmiles(smiles)
```

**做了什么**：将 SMILES 文本字符串（如 `"COc1ccccc1O"`）转换为 RDKit 的分子对象（`Mol`）。

**内部过程**：
1. 词法分析：逐字符读取 SMILES 符号（`C`=碳, `O`=氧, `c`=芳香碳, `1`=环闭合编号...）
2. 构建初始原子-键图
3. 化学完整性验证（价态检查、化合价规则）
4. 芳香性感知（Kekulization，判断哪些键是芳香键）
5. 隐式氢原子补全
6. 立体化学指定（手性、顺反异构）

**类比**：这就像大厨**读一份用化学速记法写的食谱**。比如看到 `"COc1ccccc1O"` 要在脑子里还原出"一个甲氧基连在苯环的邻位有个羟基"的完整三维结构。这不是简单的文本读取，而是涉及大量化学规则的推理。

**耗时**：每个分子约 0.5-5 毫秒（经验估计，取决于分子复杂度）

##### 步骤 B: 化学模式匹配（第 115-127 行）

```python
# 定义 4 个 SMARTS 模式
self.hydrogen_donor = Chem.MolFromSmarts(
    "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]")
self.hydrogen_acceptor = Chem.MolFromSmarts(
    "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),"
    "$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);"
    "!$([o,s]:c:n)])]")
self.acidic = Chem.MolFromSmarts("[$([C,S](=[O,S,P])-[O;H1,-1])]")
self.basic = Chem.MolFromSmarts(
    "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),...]]")

# 对分子执行 4 次子结构匹配
self.hydrogen_donor_match = sum(mol.GetSubstructMatches(self.hydrogen_donor), ())
self.hydrogen_acceptor_match = sum(mol.GetSubstructMatches(self.hydrogen_acceptor), ())
self.acidic_match = sum(mol.GetSubstructMatches(self.acidic), ())
self.basic_match = sum(mol.GetSubstructMatches(self.basic), ())
```

**做了什么**：用 4 个化学查询模式（SMARTS）去"扫描"分子，找出哪些原子属于氢键供体、氢键受体、酸性基团、碱性基团。

**SMARTS 是什么**：SMARTS 是 SMILES 的扩展语法，用于描述"化学搜索模式"。比如 `[$([N;!H0;v3])]` 的意思是"找到一个不带氢原子、价态为3的氮原子"。

**额外开销**：注意代码中 `Chem.MolFromSmarts()` 是在 `MolGraph.__init__()` 内部调用的（第 115-122 行），这意味着**每个分子都会重新编译 4 个 SMARTS 模式**，而不是复用预编译的模式对象。这是一个设计上的低效点。

**内部过程**：
1. 将 SMARTS 编译为查询图（每个分子重复做 4 次！）
2. 使用 VF2 类算法在分子图上做子图同构匹配（Subgraph Isomorphism）
3. 返回所有匹配的原子索引

**类比**：这就像大厨拿着 4 份安检清单，逐个扫描每颗菜：
- 清单1："你是氢键供体吗？"
- 清单2："你是氢键受体吗？"
- 清单3："你是酸性基团吗？"
- 清单4："你是碱性基团吗？"

更糟的是，大厨每次处理新菜时都要**重新抄写**这 4 份清单，而不是用同一份。

**耗时**：每个分子约 1-10 毫秒（经验估计，4 次匹配的总和）

##### 步骤 C: 环结构检测（第 128 行）

```python
self.ring_info = mol.GetRingInfo()
```

**做了什么**：计算分子中的环结构信息——哪些原子在环中、环的大小是多少。

**内部过程**：使用最小环基（SSSR）算法检测所有环。后续在第 199-204 行判断每个原子是否在 3/4/5/6/7/8 元环中。

**耗时**：每个分子约 0.2-2 毫秒（经验估计）

##### 步骤 D: 原子特征提取（第 135-205 行）

```python
for _, atom in enumerate(mol.GetAtoms()):     # 遍历每个原子
    self.f_atoms.append(self.atom_features(atom))
```

对每个原子调用 `atom_features()` 方法（第 171-205 行），提取 **151 维**特征向量：

| 特征类别 | 特征名称 | 编码方式 | 维度 |
|---------|----------|---------|------|
| 原子类型 | `atomic_num` (0-99) | one-hot | 101 |
| 连接度 | `degree` (0-5) | one-hot | 7 |
| 形式电荷 | `formal_charge` (-2 to +2) | one-hot | 6 |
| 手性标签 | `chiral_tag` (0-3) | one-hot | 5 |
| 氢原子数 | `num_Hs` (0-4) | one-hot | 6 |
| 杂化类型 | `hybridization` (SP/SP2/SP3/SP3D/SP3D2) | one-hot | 6 |
| 芳香性 | `IsAromatic` | 布尔 | 1 |
| 原子质量 | `mass × 0.01` | 浮点数 | 1 |
| 隐式价态 | `implicit_valence` (0-6) | one-hot | 8 |
| 氢键受体 | SMARTS 匹配结果 | 布尔 | 1 |
| 氢键供体 | SMARTS 匹配结果 | 布尔 | 1 |
| 酸性基团 | SMARTS 匹配结果 | 布尔 | 1 |
| 碱性基团 | SMARTS 匹配结果 | 布尔 | 1 |
| 环大小 | 3/4/5/6/7/8 元环 | 6个布尔 | 6 |
| **合计** | | | **151** |

**one-hot 编码是什么**：

> 不用数字直接表示类别（"碳=6"），而是用一个长向量，在对应位置标"1"，其他全"0"。
>
> 比如原子类型有 100 种可能，如果当前原子是碳（原子序数 6），那么 one-hot 向量就是：
> `[0, 0, 0, 0, 0, 1, 0, 0, ..., 0]`（第 6 位为 1，共 101 位，含一个"未知类别"位）
>
> **类比**：就像填调查问卷。不是写"碳"两个字，而是在一张有 100 个选项的表格上勾选。每个原子都要填一张 151 个问题的问卷。

**为什么耗时**：
- Python 层面的循环（`for atom in mol.GetAtoms()`）有较大的解释器开销
- 每个 `onek_encoding_unk()` 调用都创建一个新列表并逐元素赋值
- 大量临时 Python 对象的创建和销毁导致内存碎片

**耗时**：每个分子约 0.1-1 毫秒（经验估计，取决于原子数量）

##### 步骤 E: 键特征提取（第 143-169 行）

```python
for a1 in range(self.n_atoms):
    for a2 in range(a1 + 1, self.n_atoms):
        bond = mol.GetBondBetweenAtoms(a1, a2)
        if bond is None:
            continue
        f_bond = self.bond_features(bond)         # 14 维键特征
        self.f_bonds.append(self.f_atoms[a1] + f_bond)  # 正向: 151 + 14 = 165
        self.f_bonds.append(self.f_atoms[a2] + f_bond)  # 反向: 151 + 14 = 165
```

**做了什么**：
1. 遍历所有原子对 `(a1, a2)`，检查它们之间是否有化学键
2. 对每条键，提取 14 维特征（键类型、共轭、环内/外、立体化学）
3. **创建双向边**：每条键生成 2 个方向的特征向量（a1→a2 和 a2→a1）
4. 每个方向的特征 = 源原子特征(151) + 键特征(14) = **165 维**

键特征向量的 14 维构成（第 207-230 行）：

| 特征 | 维度 | 说明 |
|------|------|------|
| 空键标记 | 1 | 该键是否存在 |
| 键类型 | 4 | 单键/双键/三键/芳香键 |
| 共轭性 | 1 | 是否共轭 |
| 环内键 | 1 | 是否在环中 |
| 立体化学 | 7 | 6种可能 + 未知 (one-hot) |

**类比**：大厨不仅要了解每颗菜（原子），还要搞清楚**每根牙签**（化学键）是怎么插的：是单根的还是双根的？是直的还是扭的？而且每根牙签要从两头各看一遍。

**注意嵌套循环**：`for a1 ... for a2 ...` 是 O(n²) 复杂度。虽然实际键数远小于原子对数（大部分原子对之间没有键），但遍历检查本身就是 O(n²) 的。

**耗时**：每个分子约 0.2-2 毫秒（经验估计）

##### 步骤 F: 批次图合并 — `BatchMolGraph`（第 233-296 行）

```python
class BatchMolGraph:
    def __init__(self, mol_graphs, args):
        f_atoms = [[0] * self.atom_fdim]     # 初始零填充行
        f_bonds = [[0] * self.bond_fdim]     # 初始零填充行

        for mol_graph in mol_graphs:         # 合并 32 个分子的图
            f_atoms.extend(mol_graph.f_atoms)
            f_bonds.extend(mol_graph.f_bonds)
            # 调整索引偏移量...

        # 转换为 PyTorch 张量（CPU 张量）
        self.f_atoms = torch.FloatTensor(f_atoms)    # 第 288 行
        self.f_bonds = torch.FloatTensor(f_bonds)    # 第 289 行
        self.a2b = torch.LongTensor([...])           # 第 290 行 (含 padding)
        self.b2a = torch.LongTensor(b2a)             # 第 291 行
        self.b2revb = torch.LongTensor(b2revb)       # 第 292 行
        self.a2a = self.b2a[self.a2b]                # 第 294 行
```

**做了什么**：
1. 将 32 个独立的分子图合并为一个大的批次图
2. 调整索引，使每个分子的原子/键编号在批次中唯一
3. 将所有 Python 列表转换为 PyTorch CPU 张量

**输出**（由 `get_components()` 返回的 8 元组）：

| 张量 | 形状 | 说明 |
|------|------|------|
| `f_atoms` | (总原子数+1, 151) | 所有原子特征（含零填充行） |
| `f_bonds` | (总键数+1, 165) | 所有键特征（含零填充行） |
| `a2b` | (总原子数+1, 最大入度) | 原子→入边索引映射 |
| `b2a` | (总键数+1,) | 边→源原子索引映射 |
| `b2revb` | (总键数+1,) | 边→反向边索引映射 |
| `a_scope` | (32, 2) | 每个分子的(起始原子, 原子数) |
| `b_scope` | (32, 2) | 每个分子的(起始键, 键数) |
| `a2a` | (总原子数+1, 最大入度) | 原子→邻居原子索引映射 |

**为什么张量转换耗时**：`torch.FloatTensor(python_list)` 需要遍历 Python 列表中的每个元素，将其复制到连续的 C 内存中。对于深嵌套的列表（列表中的列表），这涉及大量的 Python 对象解包。

**类比**：32 个菜做好后，大厨还要**把它们整整齐齐地摆到一个标准化托盘上**（张量化），确保每个菜的位置有据可查（索引映射）。这个"摆盘"过程本身也需要时间。

**耗时**：每个 batch 约 0.5-2 毫秒（经验估计）

---

### 4.4 Phase 4: GPU 前向推理

**代码位置**: [models.py:306-363](src/other_softwares/grover_software/grover/model/models.py#L306-L363)

```python
class GroverFpGeneration(nn.Module):
    def forward(self, batch, features_batch):
        _, _, _, _, _, a_scope, b_scope, _ = batch

        # ① GPU 推理：GROVEREmbedding（含 GTransEncoder）
        output = self.grover(batch)

        # ② GPU 池化：对原子级嵌入做均值聚合
        mol_atom_from_bond_output = self.readout(output["atom_from_bond"], a_scope)
        mol_atom_from_atom_output = self.readout(output["atom_from_atom"], a_scope)

        # ③ CPU→GPU：将预计算特征移到 GPU
        if features_batch[0] is not None:
            features_batch = torch.from_numpy(np.stack(features_batch)).float()
            if self.iscuda:
                features_batch = features_batch.cuda()

        # ④ GPU 拼接：组合所有嵌入
        fp = torch.cat([mol_atom_from_atom_output,
                        mol_atom_from_bond_output,
                        mol_bond_from_atom_output,
                        mol_bond_from_bodd_output], 1)
        if features_batch is not None:
            fp = torch.cat([fp, features_batch], 1)

        return output["atom_from_bond"], output["atom_from_atom"], fp
```

**GPU 做了什么**：
1. **GTransEncoder**：多层 Transformer（多头注意力 + 前馈网络），这是 GROVER 的核心。对批次图进行消息传递，生成原子/键嵌入。全部是矩阵乘法，GPU 最擅长。
2. **Readout**：对每个分子的原子嵌入取平均，得到分子级表示。
3. **特征拼接**：将 4 种嵌入 + 预计算特征拼接成最终指纹。

**为什么 GPU 这么快**：

GPU 做的核心操作是**矩阵乘法**（General Matrix Multiplication, GEMM）。RTX 4070 Super 有：
- **7168 个 CUDA 核心**，每个核心可以并行执行浮点运算
- **FP32 算力约 35 TFLOPS**（每秒 35 万亿次浮点运算）

对于一个 batch 的 32 个小分子（平均 ~15 个原子），需要处理的张量尺寸大约是：
- `f_atoms`: ~(480, 151) → 经过 Transformer → ~(480, hidden_size)
- 这种规模的矩阵乘法对 RTX 4070 Super 来说**微不足道**

**类比**：GPU 就像一台有 7168 个火眼的工业烤箱。你给它 32 个小披萨，它 0.1 秒就全烤好了。它的能力远超当前需求。

**每个 batch 的 GPU 推理耗时**：约 50-150 毫秒

**14 个 batch 的 GPU 总耗时**：约 **1-2 秒**

---

### 4.5 Phase 5: 写入 LMDB

**代码位置**: [fingerprint.py:51-65](src/other_softwares/grover_software/task/fingerprint.py#L51-L65)

```python
for item in mol_loader:
    _, batch, features_batch, _, _ = item
    with torch.no_grad():
        _, _, _, _, _, a_scope, _, _ = batch
        batch_preds = model(batch, features_batch)    # GPU 推理
        for id, i in enumerate(a_scope):
            with env.begin(write=True, buffers=True) as txn:
                data = {
                    'embedding': np.concatenate((...), axis=1),      # (n_atoms, 2400)
                    'total_embedding': batch_preds[2][id].cpu().numpy(),  # (4885,)
                }
                txn.put(key=str(id+st).encode(), value=pickle.dumps(data))
    st += 32
```

**做了什么**：将 GPU 输出的嵌入从显存拷回 CPU，序列化后写入 LMDB 数据库。这一步涉及：
- **GPU→CPU 同步点**：`.cpu().numpy()` 调用会强制等待 GPU 计算完成，阻止流水线重叠
- **pickle 序列化**：每个分子的嵌入数据独立序列化
- **LMDB 事务写入**：每分子一个写事务

这一阶段可能不是最大的时间瓶颈，但它确实会增加 CPU 开销并强制 GPU-CPU 同步，使得 GPU 计算更加"碎片化"。具体占比需要 profiler 验证。

---

## 五、时间线对比图

### 5.1 本次实际运行的时间分配

```
┌────────────────────────────────────────────────────────────────────┐
│  时间轴 (0s ────────────────────────────────── 54.5s)              │
├────────────────────────────────────────────────────────────────────┤
│                                                                    │
│  CPU:  ████████████████████████████████████████████████████████░░  │
│        ├── Step 1-3 (CSV/Features/Vocab) ──┤├── Step 4 循环 ──────┤│
│        │<────── ~20s (纯CPU) ──────────────>│<─── ~34s ──────────>││
│        │                                    │                      │
│  GPU:  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░█░█░░ │
│                                                             ↑ ↑    │
│                                       间歇性的短暂推理突发 ─────┘    │
│                                       (每次 ~100ms, 共 14 次)      │
│                                                                    │
│  ████ = 忙碌  ░░░░ = 空闲                                          │
└────────────────────────────────────────────────────────────────────┘
```

### 5.2 单个 batch 的微观时间线（示意图）

> **注意**：以下时间数值为经验估计，仅供理解相对关系，非 profiler 实测结果。

```
时间 →  |  CPU (Worker进程)          |  GPU (RTX 4070 Super)
────────┼────────────────────────────┼──────────────────────
0ms     | [解析 SMILES #1]           | [空闲 💤]
2ms     | [解析 SMILES #2]           | [空闲 💤]
4ms     | [解析 SMILES #3]           | [空闲 💤]
...     | [解析 SMILES #4-#31]       | [空闲 💤]
60ms    | [解析 SMILES #32]          | [空闲 💤]
70ms    | [SMARTS 模式匹配 ×4]       | [空闲 💤]
150ms   | [提取原子/键特征]          | [空闲 💤]
200ms   | [构建 BatchMolGraph]       | [空闲 💤]
230ms   | [Python列表→PyTorch张量]   | [空闲 💤]
250ms   | [通过 IPC 传输到主进程]    | [空闲 💤]
────────┼────────────────────────────┼──────────────────────
280ms   | [等待 GPU 完成...]         | [接收数据]
285ms   | [等待...]                  | [GTransEncoder 推理 🔥]
350ms   | [等待...]                  | [Readout + 拼接]
380ms   | [接收结果, 写入LMDB]       | [完成 ✅, 等下一批 💤]
────────┼────────────────────────────┼──────────────────────

CPU 工作时间: ~280ms (占 74%)
GPU 工作时间: ~100ms (占 26%)
但实际 GPU 利用率≠26%——因为 nvidia-smi 取的是滚动平均值，
而 GPU 推理是突发式的（burst），平均下来就只有 ~1-3%
```

### 5.3 为什么 `nvidia-smi` 显示 ~1%

`nvidia-smi` 的 GPU Utilization 指标是一个**采样式的滚动平均值**，采样周期通常为 1 秒。它反映的是"过去 1 秒内，GPU 有多少比例的时间在执行至少一个 kernel"。

在我们的场景中：
- GPU 在 1 秒内可能只被激活了 1-2 次（处理 1-2 个 batch）
- 每次激活持续 ~100ms
- 所以 1 秒内的活跃占比 ≈ 100-200ms / 1000ms = **10-20%**

但是考虑到很多时候连续好几秒 GPU 完全空闲（CPU 在处理上一步的 save_features 或 build_vocab），**跨整个 54.5 秒的平均 GPU 活跃时间 ≈ 1.5s / 54.5s ≈ 2.7%**。而 `nvidia-smi` 的采样可能正好落在空闲期，因此观测到 ~1% 是合理的。

---

## 六、缓存失效问题

### 6.1 设计意图

代码在 [fingerprint.py:35](src/other_softwares/grover_software/task/fingerprint.py#L35) 中创建了一个空字典 `shared_dict={}`，传给 `MolCollator`。设计意图是：如果同一个 SMILES 在数据集中出现多次，第二次遇到时可以直接从缓存取出已解析的 `MolGraph`，避免重复计算。

### 6.2 为什么失效

**问题在于 `num_workers=4`**。

当 DataLoader 使用多个 worker 时：
- 每个 worker 运行在**独立的子进程**中
- 在 Windows 上，Python 的多进程使用 **spawn** 方式（不是 Linux 的 fork）
- spawn 方式下，子进程不继承父进程的内存，而是通过 **pickle 序列化**重新创建所有 Python 对象
- 这意味着：每个 worker 都有自己**独立的** `shared_dict` 副本

结果：
```
主进程:  shared_dict = {}  (交给DataLoader)
         ↓ pickle 序列化 ↓
Worker0: shared_dict = {}  → 解析分子0-31,  缓存到自己的dict
Worker1: shared_dict = {}  → 解析分子32-63, 缓存到自己的dict
Worker2: shared_dict = {}  → 解析分子64-95, 缓存到自己的dict (看不到Worker0/1的结果)
Worker3: shared_dict = {}  → 解析分子96-127,缓存到自己的dict (看不到Worker0/1/2的结果)
```

**类比**：4 个厨师在 4 个隔音的独立厨房里工作。厨师 A 算出了"苯环有 6 个碳"写在小纸条上，但厨师 B 看不到这张纸条，必须自己重新算一遍。

### 6.3 实际影响

对本次 436 个底物的场景，影响不大：
- 436 个底物全部是**不同的分子**（没有重复 SMILES）
- 即使缓存完全生效，也不会减少计算量
- 缓存失效的问题只在**训练场景**（相同分子在不同 epoch 中反复出现）时才会造成显著性能损失

---

## 七、"我的 CPU 太差了吗？"

### 7.1 简短回答

**不是。** 这个瓶颈与 CPU 的绝对性能关系不大。

### 7.2 详细分析

CPU 预处理慢的根本原因是**架构层面的**，而非硬件层面的：

| 因素 | 说明 |
|------|------|
| **Python 解释器开销** | 原子/键特征的提取都是纯 Python 循环。CPython 解释器执行循环比 C/C++ 慢 10-100 倍。即使 RDKit 的底层是 C++，每次 Python→C++ 的函数调用仍有跨语言开销。 |
| **单分子内部串行** | 单个分子的图构建是串行的逻辑操作链：先解析 SMILES、再检测环、再逐原子提取特征、再逐键提取特征。这些步骤有严格的依赖关系。虽然可以跨分子并行（`num_workers=4`），但受 Windows spawn 开销、IPC 序列化、缓存不共享等因素限制，并行收益递减。 |
| **GPU 太快了** | RTX 4070 Super 的 FP32 算力是 35 TFLOPS。处理 32 个小分子的 Transformer 推理，只是九牛一毛。即使 CPU 速度翻倍，GPU 仍然会在绝大多数时间处于空闲状态。 |
| **小数据集效应** | 436 个分子太少了。如果有 10 万个分子，CPU 预处理可以用流水线方式和 GPU 推理重叠（当 GPU 在推理第 N 批时，CPU 在准备第 N+1 批）。但 14 个 batch 太少，流水线根本来不及填满。 |

### 7.3 升级 CPU 会有帮助吗？

| 升级方案 | 预期效果 | 说明 |
|---------|---------|------|
| 更高主频 CPU | 小幅改善（10-30%） | Python 循环和 RDKit 调用的延迟会降低，但瓶颈不会消除 |
| 更多核心 CPU | 收益递减 | 当前 `num_workers=4`，增加 worker 需要更多核心，但 Windows spawn 模式下有进程启动开销、IPC 序列化开销和缓存不共享问题，实际收益需实测 |
| 顶级服务器 CPU | 改善但仍受限 | 假设 CPU 处理速度翻 10 倍（5.45s），GPU 推理（1.5s）的占比提升到 1.5/(5.45+1.5)≈22%，利用率仍然偏低。说明架构瓶颈比硬件瓶颈更根本 |

**核心结论**：问题不在于"CPU 太慢"，而在于"GPU 太快 + 数据预处理是串行逻辑操作"。这是一个**架构瓶颈**，不是硬件瓶颈。

---

## 八、与训练场景的对比

理解推理（inference）和训练（training）场景下的差异，有助于全面理解这个问题。

### 8.1 推理场景（本次情况）

```
数据量: 436 个分子, 14 个 batch
经过次数: 1 次
CPU 总预处理: ~33 秒
GPU 总推理: ~1.5 秒
总耗时: 54.5 秒
```

每个分子只需要处理**一次**。CPU 瓶颈存在但对绝对耗时的影响可接受。

### 8.2 训练场景（假设）

```
数据量: 假设 10,000 个分子, 313 个 batch
epoch 数: 100
如果不缓存/不预计算:
  CPU 总预处理: ~750 秒/epoch × 100 = ~75,000 秒 (≈ 20.8 小时)
  GPU 总训练: ~30 秒/epoch × 100 = ~3,000 秒 (≈ 50 分钟)
  总耗时: ~21.7 小时 (GPU 空闲 96% 的时间!)

如果预计算图数据:
  CPU 预计算(一次): ~750 秒
  GPU 总训练: ~3,000 秒
  总耗时: ~1 小时 (GPU 利用率大幅提升)
```

在训练场景下，CPU 瓶颈的代价被 **epoch 次数放大**，影响从"可接受"变成"灾难性的"。这就是为什么严肃的训练管线都会**预计算图数据**。

---

## 九、优化策略（如果将来需要）

以下策略按**影响大小**排序。对于当前 436 个分子的一次性推理，**不需要任何优化**。仅在将来处理大规模数据或进行模型训练时参考。

### 9.1 预计算图数据（最推荐，影响最大）

**原理**：将 RDKit 解析和特征提取从 DataLoader 热路径中移出，预先计算并保存到磁盘。

```
离线阶段（只做一次）:
  SMILES → RDKit 解析 → 特征提取 → 保存为 .pt 或 .npz 文件

在线阶段（每次推理/每个epoch）:
  直接从磁盘加载预计算的张量 → 送入 GPU
```

**效果**：消除 CPU 热路径中的 RDKit 操作，显著提升吞吐。GPU 利用率能否达到 50-80% 取决于 batch_size、模型结构、是否还有 `.cpu()` 同步点等因素，需实测验证。

**类比**：不在客人来的时候才切菜，而是在开店前就把所有食材切好做成半成品放在冰箱里。客人来了直接拿出来加热。

### 9.2 设置 `num_workers=0`（简单，对小数据集有效）

**原理**：不使用多进程，所有操作在主进程中完成。

**好处**：
- 消除 Windows spawn 的进程启动开销
- `shared_dict` 缓存真正生效（同一进程内共享）
- 消除进程间 pickle 序列化/反序列化开销

**适用场景**：数据集较小（< 1000 个分子）

### 9.3 增大 `batch_size`（简单，提升 GPU 效率）

**原理**：更大的 batch 意味着 GPU 每次处理的矩阵更大，能更好地利用并行计算资源。

**当前**：`batch_size=32`（GPU 严重"吃不饱"）
**建议**：如果显存允许，尝试 `batch_size=128` 或 `256`

**效果**：不减少 CPU 预处理时间，但提高了 GPU 每次被唤醒时的工作效率。

### 9.4 使用 `persistent_workers=True`

**原理**：保持 worker 进程在整个推理过程中存活，避免每个 epoch 重新启动进程。

```python
mol_loader = DataLoader(data, batch_size=32, num_workers=4,
                        collate_fn=mol_collator,
                        persistent_workers=True)  # 新增
```

**效果**：减少 Windows 上的进程启动/销毁开销。对短任务效果显著。

### 9.5 使用 `pin_memory=True` + 非阻塞传输

**原理**：使用锁页内存加速 CPU→GPU 的数据传输。

**效果**：对本场景影响极小（数据传输不是瓶颈），但在大规模训练中有一定帮助。

---

## 十、总结

| 问题 | 回答 |
|------|------|
| GPU 利用率为什么只有 ~1%？ | 54.5 秒中 GPU 只工作了 ~1.5 秒。其余时间都在等 CPU 做数据预处理。 |
| CPU 在做什么？ | 解析 SMILES → 化学模式匹配 → 提取原子/键特征 → 构建图数据结构 → 转为张量。 |
| 是 CPU 太差了吗？ | 不是。这是架构层面的瓶颈：图构建是串行逻辑操作，天然比 GPU 的矩阵乘法慢。 |
| 代码有 bug 吗？ | 没有。代码逻辑正确，GPU 确实被使用了，只是利用率低。 |
| 需要优化吗？ | 对 436 个分子的一次性推理：**不需要**，54.5 秒完全可以接受。 |
| 什么时候需要优化？ | 处理 10 万+ 分子，或需要反复训练时，应采用预计算策略。 |

### 核心洞见

> **GPU 利用率低 ≠ 程序有问题。**
>
> 它反映的是数据管线的特征：CPU 做的是**逻辑密集型**工作（解析化学结构、提取特征），GPU 做的是**计算密集型**工作（矩阵乘法）。两者的速度差异导致了 GPU 的"等待饥饿"。
>
> 这在计算化学领域是一个非常普遍的现象，尤其是在使用图神经网络处理分子数据时。

---

## 附录 A: 代码文件引用索引

| 文件 | 路径 | 关键行号 | 角色 |
|------|------|---------|------|
| fingerprint.py | src/other_softwares/grover_software/task/ | 35-42, 51-65 | 入口脚本，创建 DataLoader 和推理循环 |
| molgraph.py | src/other_softwares/grover_software/grover/data/ | 80-230, 233-296, 347-387 | MolGraph/BatchMolGraph/MolCollator |
| moldataset.py | src/other_softwares/grover_software/grover/data/ | 18-50 | MoleculeDatapoint/MoleculeDataset |
| models.py | src/other_softwares/grover_software/grover/model/ | 306-363 | GroverFpGeneration 前向推理 |

## 附录 B: 关键数值汇总

| 数值 | 来源 | 含义 |
|------|------|------|
| 436 | 数据集 | 底物总数 |
| 32 | fingerprint.py:39 | batch_size |
| 14 | 436/32 | 总 batch 数 |
| 4 | fingerprint.py:37 | DataLoader worker 数 |
| 151 | ATOM_FDIM(133)+18 | 每原子特征维度 |
| 14 | BOND_FDIM | 每键特征维度 |
| 165 | 151+14 | 每有向边特征维度（f_bonds） |
| 2400 | 模型输出 | 原子级嵌入维度 |
| 4885 | 模型输出 | 分子级嵌入维度 |
| ~1.5s | 实测估算 | GPU 总活跃时间 |
| ~53s | 实测估算 | CPU 总预处理时间 |
| 54.5s | 实测 | 端到端总耗时 |

---

**文档版本**: v1.0
**撰写日期**: 2026-01-30
**多模型协作**: Claude Code (架构分析 + 整合) + Codex (技术深度分析) + Gemini (通俗化转译)
