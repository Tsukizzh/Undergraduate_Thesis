# Step 7.2: GROVER指纹生成 - 详细会话日志

> **执行日期**: 2026-01-30
> **版本**: v1.0
> **执行者**: Claude Code + Codex + Gemini (多模型协作)

---

## 概述

**使用预训练的GROVER-Large模型为全部436个底物生成GROVER神经指纹，产出原子级（2400维）和分子级（4885维）嵌入，存储在LMDB中。**

---

## 背景知识

### 什么是GROVER？

GROVER（Graph Representation frOm self-superVised messagE passing tRansformer，自监督消息传递Transformer的图表示）是一个预训练的分子表示模型。与使用手工规则的Morgan指纹不同，GROVER使用深度神经网络，该网络在1000万个分子上预训练，以学习有意义的分子表示。

**类比**：Morgan指纹就像用清单描述一个人（身高、体重、眼睛颜色）。GROVER则像是让一位AI艺术家画一幅肖像——它能捕捉到清单无法描述的微妙特征。

### 两种GROVER嵌入类型

| 类型 | 代码中的名称 | 形状 | 描述 |
|------|-------------|------|------|
| 原子级 | `embedding` | (n_atoms, 2400) | 每个原子得到一个2400维向量 |
| 分子级 | `total_embedding` | (4885,) | 整个分子汇总为4885维 |

### 为什么需要两个级别？

EZSpecificity模型同时使用两者：
- `atom_features: ["grover", 2400]` → 原子级用于图神经网络
- `features: ["grover_mean", 4885]` → 分子级作为全局上下文

---

## GROVER四步流水线

```
步骤1：创建输入CSV
   Substrates.csv → grover_substrates.csv（单SMILES列）

步骤2：save_features.py
   SMILES → 原子/键描述符（.npz）
   [化学特征如官能团标签]

步骤3：build_vocab.py（低内存模式）
   SMILES → 分子词汇表（.pkl）
   [原子词汇：150种类型，键词汇：212种类型]

步骤4：main.py fingerprint
   CSV + NPZ + grover_large.pt → LMDB
   [神经网络在GPU上前向传播]
```

---

## 输入

### 文件：`Substrates.csv`

| 字段 | 描述 | 示例 |
|------|------|------|
| Substrate_Index | 底物ID（从0开始）| 0, 1, 2, ..., 435 |
| Substrate_SMILES | 分子SMILES表示 | COc1ccccc1O |

**总计：436个底物**

### 预训练模型：`grover_large.pt`

| 属性 | 值 |
|------|-----|
| 位置 | `data/pretrain_model/grover_large.pt` |
| 大小 | 约1 GB |
| 架构 | GROVER-Large（消息传递Transformer）|
| 预训练数据 | 1000万个无标签分子 |

---

## 输出

### 文件：`grover_fingerprint.lmdb`

| 属性 | 值 |
|------|-----|
| 记录数 | 436 |
| 键格式 | "0", "1", ..., "435"（字符串编码）|
| 文件大小 | 约10 GB（LMDB map_size预分配）|

### 记录结构

每条记录（通过 `pickle.loads()` 访问）包含：

| 字段 | 形状 | 数据类型 | 描述 |
|------|------|----------|------|
| `embedding` | (n_atoms, 2400) | float32 | 原子级表示 |
| `total_embedding` | (4885,) | float32 | 分子级表示 |

### 抽样检查结果

| 键 | 底物 | embedding | total_embedding |
|-----|------|-----------|-----------------|
| "0" | COc1ccccc1O（9个原子）| (9, 2400) | (4885,) |
| "1" | COc1cccc(OC)c1O（11个原子）| (11, 2400) | (4885,) |
| "2" | CN[C@H](...)（22个原子）| (22, 2400) | (4885,) |

### 中间文件

| 文件 | 大小 | 描述 |
|------|------|------|
| grover_substrates.csv | 20 KB | GROVER输入（仅SMILES）|
| grover_substrates.npz | 3.7 KB | 原子/键特征描述符 |
| grover_vocab/ | 约1 KB | 原子词汇（150）+ 键词汇（212）|

---

## 模型如何加载GROVER指纹

来自 `src/Datasets/data_representer.py`：

```python
# 加载LMDB（第100-109行）
self.grover_dbs = [lmdb.open(path, ..., readonly=True, lock=False) for path in ...]

# 访问分子级（第185-186行）
grover_data = pickle.loads(txn.get(str(substrate_index).encode()))
substrate_features["grover_mean"] = grover_data['total_embedding']   # (4885,)

# 访问原子级（第186-187行）
substrate_features["grover"] = grover_data['embedding']              # (n_atoms, 2400)
```

---

## 执行

### 命令
```bash
D:\anaconda3\envs\torch\python.exe step7_2_generate_grover.py
```

### 结果
```
GROVER Fingerprint Generation Complete
============================================================
Input:       .../04_Step4_格式修正后数据/Substrates.csv
Output LMDB: .../07_Step7_分子指纹/grover_fingerprint.lmdb
Substrates:  436
Checkpoint:  .../data/pretrain_model/grover_large.pt
CUDA:        Auto (GPU used for inference)
Status:      PASS
============================================================

Total elapsed: 54.5s
```

### 硬件
- GPU：NVIDIA GeForce RTX 4070 Super（CUDA推理）
- 处理时间：总计54.5秒
  - 步骤1-3（CPU）：约20秒
  - 步骤4（GPU推理）：约35秒

---

## 遇到的问题与修复

### 问题1：缺少 `descriptastorus` 模块

```
ModuleNotFoundError: No module named 'descriptastorus'
```

**修复**：`pip install descriptastorus`

### 问题2：内存爆炸（98% RAM，31.2/31.8 GB）

**根本原因**：原始 `build_vocab.py` 使用 `Pool(100)` 创建100个进程。在Windows上（基于spawn的多进程），每个进程复制父进程的内存，导致OOM。

**修复**：创建自定义 `build_vocab_low_memory()` 函数，直接调用 `MolVocab` 并设置 `num_workers=1`。

### 问题3：LMDB中文路径错误

```
lmdb.Error: C:\...\P450_EZSpecificity_研究项目\...\grover_fingerprint.lmdb: 乱码
```

**根本原因**：Windows上的LMDB无法处理文件路径中的非ASCII字符。

**修复**：使用临时的纯ASCII路径（`C:\grover_temp\`）进行GROVER处理，然后将结果复制到最终目标。

### 问题4：LMDB 600GB map_size错误

```
lmdb.Error: C:\grover_temp\grover_fingerprint.lmdb: 乱码
```

**根本原因**：LMDB `map_size=600*(1024**3)`（600GB）对于Windows虚拟地址空间分配来说太大。

**修复**：在 `task/fingerprint.py` 中将 `map_size` 减少到 `10*(1024**3)`（10GB）。同时更新 `data_representer.py` 以保持一致。

### 代码修改汇总

| 文件 | 修改 | 原因 |
|------|------|------|
| `step7_2_generate_grover.py` | 添加 `build_vocab_low_memory()` | Windows上Pool(100)导致OOM |
| `step7_2_generate_grover.py` | 使用 `C:\grover_temp` 临时目录 | LMDB中文路径不兼容 |
| `src/.../fingerprint.py` 第45行 | map_size: 600GB → 10GB | Windows虚拟地址空间限制 |
| `src/Datasets/data_representer.py` 第102,117,129行 | map_size: 600GB → 10GB | 保持Windows兼容性一致 |

---

## 多模型审查

| 轮次 | 审查者 | 发现 | 结果 |
|------|--------|------|------|
| 规划 | Codex + Gemini | Morgan + GROVER实施策略 | 批准 |
| Codex审查1 | Codex | GROVER输出维度与模型配置匹配 | 批准 |
| Codex审查2 | Codex | LMDB索引键"0".."435"正确 | 批准 |
| Codex审查3 | Codex | 发现data_representer.py 600GB map_size问题 | 已修复 |

---

## 验证清单

| 检查项 | 结果 |
|--------|------|
| LMDB记录数 = 436 | 通过 |
| 所有键"0".."435"存在 | 通过 |
| embedding形状 = (n_atoms, 2400) | 通过 |
| total_embedding形状 = (4885,) | 通过 |
| 使用GPU进行推理 | 通过 |
| 无缺失/多余的键 | 通过 |
| 与模型配置一致 | 通过 |

---

## 工具与环境

| 工具 | 版本 | 用途 |
|------|------|------|
| Python | 3.12.5 (torch环境) | 运行时 |
| PyTorch | 2.1.0+cu121 | GROVER模型推理 |
| GROVER | Large模型 (grover_large.pt) | 预训练分子编码器 |
| RDKit | 2024.03+ | SMILES处理 |
| LMDB | 1.4.1 | 特征存储 |
| descriptastorus | latest | GROVER特征生成 |
| CUDA | 12.1 | GPU加速 |

---

## 文件位置

```
PathA_2026-01-08_模型评估测试集构建/
├── data/
│   ├── 04_Step4_格式修正后数据/
│   │   └── Substrates.csv                    ← 输入（436个底物）
│   └── 07_Step7_分子指纹/
│       ├── grover_fingerprint.lmdb           ← 主输出（436条记录）
│       ├── grover_substrates.csv             ← GROVER输入CSV
│       ├── grover_substrates.npz             ← 原子/键特征
│       ├── grover_vocab/                     ← 词汇文件
│       │   ├── test_atom_vocab.pkl
│       │   └── test_bond_vocab.pkl
│       └── morgan_fingerprint.npy            ← （来自Step 7.1）
└── scripts/
    └── 07_Step7_分子指纹生成/
        └── step7_2_generate_grover.py        ← 脚本
```

---

## 后续步骤

Step 7（Morgan + GROVER指纹）现已完成。剩余的特征生成步骤是：
- **Step 8**：结构特征生成（从PDB文件中提取口袋/配体）

所有特征生成完成后，进行模型推理。

---

**文档版本**: v1.0
**最后更新**: 2026-01-30
