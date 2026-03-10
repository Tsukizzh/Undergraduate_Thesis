# EZSpecificity 训练管线与模型内部详解

**创建日期**: 2026-03-09
**目的**: 面向训练新手的完整技术文档。从"零训练经验"出发，讲清楚模型是什么、怎么训练、怎么微调。
**前置知识**: 了解 Python 基础语法、知道什么是神经网络（大概就行）

---

## 一、一句话总结

EZSpecificity 是一个**二分类模型**，判断"给定的酶能否催化给定的底物"：

```
输入：一个酶的氨基酸序列 + 一个底物的化学结构 + 两者的3D对接构象
输出：一个分数（越高越可能催化）
```

模型用了三个"信息通道"来做判断，最后投票得出结论：

```
通道1: 序列通道 — "酶的氨基酸序列说了什么？"    （ESM嵌入）
通道2: 分子通道 — "底物的化学结构长什么样？"    （GROVER + Morgan指纹）
通道3: 结构通道 — "酶和底物在3D空间怎么靠近？"  （SE(3)-GNN）
```

---

## 二、PyTorch Lightning 训练框架

### 本章与模型架构的关系

本文档后续章节会详细讲解模型的核心组件：ESM 蛋白质嵌入、GNN（图神经网络，处理 3D 结构）、双向交叉注意力等（第四章）。这些组件构成了模型的**计算逻辑**——给定一对酶和底物的特征，经过这些组件的计算，最终输出一个特异性预测分数。

但仅有计算逻辑是不够的。训练一个模型还需要回答：数据怎么分批送进来？loss 怎么算？梯度怎么更新？学习率怎么调？什么时候保存最佳模型？**PyTorch Lightning 就是管理这些训练流程的框架**。

两者的关系：

```
PyTorch Lightning 框架（本章内容）
├─ Trainer         控制训练循环（epoch、梯度更新、验证、保存）
├─ Singledataset   把原始数据打包成 batch 送给模型
└─ SS 类           包含以下两部分：
    ├─ 模型计算逻辑（第四章内容）
    │   ├─ ESM 投影层
    │   ├─ GNN（图神经网络，处理 3D 结构）
    │   ├─ 双向交叉注意力
    │   └─ 预测头 MLP
    └─ 训练/验证逻辑（本章内容）
        ├─ training_step()     调用 forward + 算 loss
        ├─ validation_step()   调用 forward + 记录指标
        └─ configure_optimizers()  定义优化器和学习率策略
```

换句话说：**第四章讲的是 SS 类的 `forward()` 内部做了什么计算，本章讲的是 `forward()` 是怎么被调用、训练流程是怎么组织的**。

### 2.1 三个核心组件

PyTorch Lightning 是 PyTorch 的封装框架。传统 PyTorch 需要手动写训练循环（`for epoch`, `for batch`, `loss.backward()`, `optimizer.step()`...），Lightning 把这些重复代码封装了，你只需要定义"做什么"。

它有三个核心组件：

| 组件 | 本项目对应 | 代码文件 | 职责 | 涉及的方法 |
|------|-----------|---------|------|-----------|
| `LightningModule` | SS 类 | [ss.py](src/Models/ss.py) | 定义模型结构 + 训练/验证逻辑 | `forward`, `training_step`, `validation_step`, `get_loss`, `configure_optimizers` |
| `LightningDataModule` | Singledataset 类 | [brenda.py](src/Datasets/brenda.py) | 加载 CSV + LMDB → 构建 DataLoader | `train_dataloader`, `val_dataloader`, `test_dataloader` |
| `Trainer` | **需要新写** | 目前仅有推理脚本 [main_testing.py](src/main_testing.py) | 控制训练循环、调度验证、保存检查点 | `trainer.fit()` |

### 2.2 训练的整体过程（简化版）

训练的目标：反复调整模型参数，让模型的预测越来越准。

一轮训练（一个 epoch）做两件事：先**学习**，再**考试**。

```
一个 Epoch = 一轮"学习 + 考试"

┌─────────────── 学习（训练阶段）──────────────────────────────────┐
│                                                                   │
│  训练数据被分成很多批（每批 16 个酶-底物对）。                       │
│  对每一批，执行以下 4 步：                                         │
│                                                                   │
│  1. 取一批数据                                                     │
│     每个样本包含：一个酶的蛋白质序列、一个底物的分子式、             │
│     它们的对接 3D 结构，以及正确答案（1=能催化，0=不能催化）         │
│                                                                   │
│  2. 模型做预测                                                     │
│     数据经过三条通道分别处理：                                      │
│     - 序列通道：处理酶的蛋白质序列 → 酶的特征                      │
│     - 结构通道：处理酶-底物的 3D 对接结构 → 结构特征                │
│     - 分子通道：处理底物的分子指纹 → 分子特征                      │
│     三条通道的结果通过交叉注意力融合后，输出一个预测分数              │
│     （分数越高 = 模型越认为这个酶能催化这个底物）                    │
│                                                                   │
│  3. 计算误差（loss）                                               │
│     把预测分数和正确答案对比，计算差距有多大                         │
│     使用 BCEWithLogitsLoss（二分类交叉熵损失函数）：                │
│     预测越接近正确答案 → loss 越小；预测越离谱 → loss 越大           │
│                                                                   │
│  4. 调整模型参数                                                    │
│     根据 loss，计算"每个参数应该往哪个方向调、调多少"（即梯度）     │
│     然后沿着梯度方向微调所有参数，使 loss 减小                      │
│     注意：每 2 批数据的梯度累加后才调整一次，                       │
│     避免单批数据的噪声导致调整方向不稳定                             │
│                                                                   │
│  所有批次跑完 = 训练数据全部用过一遍                                 │
│                                                                   │
└───────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────── 考试（验证阶段）──────────────────────────────────┐
│                                                                   │
│  用一批模型在训练中从未见过的数据来测试：                             │
│                                                                   │
│  1. 模型对验证数据做预测                                           │
│     过程和学习阶段的第 2 步完全一样，但不计算梯度、不调整参数         │
│                                                                   │
│  2. 汇总所有验证 batch 的预测结果，计算整体指标                     │
│     - AUC-ROC：模型区分"能催化"和"不能催化"的能力                  │
│       （0.5 = 相当于随机猜，1.0 = 完美区分）                       │
│     - AUPR：模型预测"能催化"时有多准确                             │
│     还会按酶家族、底物类型等分组，计算各组的细分指标                  │
│                                                                   │
│  3. 根据考试成绩做两个决策                                         │
│     - AUPR 连续多个 epoch 没提升？                                  │
│       → 降低学习率（让参数调整的步幅更小，更精细地搜索最优解）       │
│     - AUPR 是历史最佳？                                            │
│       → 保存当前模型参数到 .ckpt 文件（之后可以随时加载这个最佳状态）│
│                                                                   │
└───────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                     开始下一个 Epoch
                重复直到训练结束（达到上限或成绩不再提升）
```

**为什么要分学习和考试？** 如果只看训练数据上的表现，模型可能只是"背答案"（过拟合）——在训练数据上很准，但遇到新数据就不行了。验证阶段用模型没见过的数据来测试，确保模型学到的是通用规律。

### 2.3 三个组件如何协作（详细版）

下面是同一个过程的完整代码级流程图。每一步标注了由哪个组件执行、对应哪个代码文件。如果上面的简化版已经看懂，可以先跳过这里，等需要查代码时再回来看。

```
一个完整 Epoch 的执行流程（从上到下阅读）
==========================================
标注格式：[文件 / 类名]

Step 1  [需新写 / Trainer] 启动第 N 个 epoch
        │
        ▼
Step 2  [brenda.py / Singledataset] 提供训练数据
        Trainer 调用 Singledataset.train_dataloader()，获得一个迭代器
        每次迭代产出 1 个 batch = 16 个酶-底物对
        │
        ▼
Step 3  ┌─── 对每个训练 batch 循环 ──────────────────────────────────┐
        │                                                             │
        │  3a  [需新写 / Trainer] 把 batch 传给 SS.training_step()   │
        │       │                                                     │
        │       ▼                                                     │
        │  3b  [ss.py:163 / SS] 执行 SS.training_step(batch)：       │
        │       ├─ SS.forward(batch)           产出 logits [B, 1]    │
        │       ├─ SS.get_loss(batch, logits)  产出 loss（标量）      │
        │       └─ return loss                                       │
        │       │                                                     │
        │       ▼                                                     │
        │  3c  [需新写 / Trainer] 收到 loss，执行反向传播：           │
        │       ├─ loss.backward()     计算所有参数的梯度              │
        │       ├─ gradient_clip(8.0)  裁剪过大的梯度                  │
        │       └─ optimizer.step()    用梯度更新参数                  │
        │                                                             │
        │  注：accumulate_grad_batches=2，每 2 个 batch 的梯度        │
        │      累加后才执行一次 optimizer.step()                       │
        │      等效 batch_size = 16 × 2 = 32                         │
        │                                                             │
        └──────────── 循环直到所有训练 batch 跑完 ────────────────────┘
        │
        ▼
Step 4  [brenda.py / Singledataset] 提供验证数据
        Trainer 调用 Singledataset.val_dataloader()
        │
        ▼
Step 5  ┌─── 对每个验证 batch 循环 ──────────────────────────────────┐
        │                                                             │
        │  5a  [需新写 / Trainer] 把 batch 传给 SS.validation_step() │
        │       │                                                     │
        │       ▼                                                     │
        │  5b  [ss.py:176 / SS] 执行 SS.validation_step(batch)：     │
        │       ├─ SS.forward(batch)      产出 logits                 │
        │       ├─ SS.get_loss(...)       记录 loss（仅监控）         │
        │       └─ return logits, labels, tags                       │
        │                                                             │
        │  注：验证阶段不计算梯度、不更新参数                          │
        │                                                             │
        └──────────── 循环直到所有验证 batch 跑完 ────────────────────┘
        │
        ▼
Step 6  [ss.py:234 / SS] 执行 SS.validation_epoch_end()，汇总验证结果
        ├─ 拼接所有 batch 的 logits 和 labels
        ├─ 计算整体 AUC-ROC 和 AUPR
        └─ 按 tag（酶家族/底物类型）分组计算细分指标
        │
        ▼
Step 7  [需新写 / Trainer] 根据验证指标做两个决策：
        │
        ├─ ReduceLROnPlateau（配置在 ss.py:252）：
        │   AUPR/val 连续 patience 个 epoch 没提升？
        │   是 → 学习率 × factor
        │
        └─ ModelCheckpoint：
            AUPR/val 是历史最佳？
            是 → 保存 .ckpt 文件
        │
        ▼
Step 8  回到 Step 1，开始第 N+1 个 epoch
        直到达到 max_epochs 或 early stopping
```

### 2.4 训练 vs 推理：代码差异

```python
# 推理（现有代码 main_testing.py:35-48）
model = SS.load_from_checkpoint("best-checkpoint.ckpt", config=config)
trainer.test(model, datamodule=dm)     # ← 只调用 test_step() → evaluate()

# 训练（需要新写的代码）
model = SS(config)                      # 或 load_from_checkpoint 做微调
trainer.fit(model, datamodule=dm)       # ← 调用 training_step() + validation_step()
trainer.test(ckpt_path="best", datamodule=dm)
```

---

## 三、数据管线：从 CSV 到模型输入

### 3.1 总览流程图

```
                        ┌─────────────────────────────────────────┐
                        │       离线预计算（只做一次，很慢）         │
                        │                                         │
  Enzymes.csv ─────────→│  ESM-2 模型 ──→ enzyme_features.lmdb    │
  (蛋白质序列)           │  (每条序列 → 1280维向量/残基)            │
                        │                                         │
  Substrates.csv ──┬───→│  RDKit ──────→ reaction_features.lmdb   │
  (SMILES)         │    │  (原子、键、图结构)                       │
                   │    │                                         │
                   ├───→│  GROVER ─────→ grover_fingerprint.lmdb  │
                   │    │  (原子级 + 分子级嵌入)                    │
                   │    │                                         │
                   └───→│  RDKit ──────→ morgan_fingerprint.npy   │
                        │  (1024-bit 分子指纹)                     │
                        │                                         │
  Structure/*.pdb ─────→│  口袋提取 ───→ structure_features.lmdb   │
  (对接复合物)           │  (蛋白质口袋 + 配体3D坐标)               │
                        └───────────────────┬─────────────────────┘
                                            │
                                            ▼
                        ┌─────────────────────────────────────────┐
                        │      在线加载（每个 batch 实时做）         │
                        │                                         │
  data.csv ────────────→│  Singledataset (brenda.py)              │
  (酶索引,底物索引,      │  ├─ 读 CSV → 找到酶/底物/结构索引        │
   对接索引,标签)        │  ├─ 查 LMDB → 取出预计算的特征            │
                        │  ├─ 填充到固定长度                       │
                        │  │   (酶→1450残基, 底物→280原子)          │
                        │  ├─ 构建分子图 (EdgeConnection)           │
                        │  │   (k=48 近邻, 距离RBF编码)            │
                        │  └─ 打包成 batch (16个样本/批)            │
                        └───────────────────┬─────────────────────┘
                                            │
                                            ▼
                                     模型 forward()
```

### 3.2 四种特征详解

| 特征 | 含义 | 维度 | 来源 | 存储 |
|------|------|------|------|------|
| **ESM 嵌入** | 蛋白质序列的上下文表示（每个氨基酸一个向量） | [序列长度, 1280] | ESM-2 预训练语言模型 | LMDB |
| **GROVER 原子嵌入** | 底物每个原子的化学环境表示 | [原子数, 2400] | GROVER 预训练分子模型 | LMDB |
| **GROVER 分子嵌入** | 底物整体的化学表示（所有原子嵌入的聚合） | [4885] | GROVER | LMDB |
| **Morgan 指纹** | 底物的拓扑子结构编码（固定长度的 bit 向量） | [1024] | RDKit 直接计算 | NumPy |
| **3D 结构特征** | 蛋白质口袋+底物的原子坐标和化学键信息 | 图结构（节点+边） | PDB 文件 | LMDB |

**各特征的具体作用**：
- **ESM 嵌入**：ESM-2 是一个在 2.5 亿条蛋白质序列上预训练的语言模型。它把每个氨基酸编码为 1280 维向量，向量中包含了该氨基酸在整条序列上下文中的功能信息（类似 BERT 对自然语言做的事）
- **Morgan 指纹**：一种经典的分子表示方法。对分子中的每个原子，统计半径 r 内的子结构，哈希到 1024-bit 向量。它是全局级别的，不区分单个原子
- **GROVER 嵌入**：GROVER 是在大规模分子数据上预训练的图 Transformer。它能输出原子级嵌入（每个原子 2400 维）和分子级嵌入（4885 维），比 Morgan 指纹包含更丰富的化学信息
- **3D 结构特征**：直接来自蛋白质-底物对接复合物的 PDB 文件。提取口袋原子（配体 10Å 内的蛋白质原子）和配体原子的 3D 坐标、元素类型、化学键等，构建成图

### 3.3 训练 vs 推理的数据差异

| 方面 | 训练 | 推理 |
|------|------|------|
| 数据打乱 | 是（shuffle=True） | 否 |
| 距离噪声 | 加（Laplace噪声，scale≈0.032Å） | 不加 |
| 纯配体回退 | 当前配置关闭（`full_data: False`），只用高质量复合物 | 同左 |
| 输出 | loss（用于反向传播） | 预测分数（用于评估） |

---

## 四、模型内部结构（核心章节）

### 4.1 全景流程图

#### 简化版：模型在做什么

模型接收三类原始数据，经过三条通道分别处理，最终合并成一个预测分数。

```
          蛋白质序列              对接复合物 PDB             底物 SMILES
              │                      │                        │
              ▼                      ▼                        ▼
         ① 序列通道             ② 结构通道               ③ 分子通道
         ESM 嵌入投影           图神经网络处理            GROVER + Morgan
              │                      │                        │
              ▼                      ▼                        ▼
          酶的特征              结构特征（原子级）         分子特征（原子级+全局级）
         x_pro                  x_str                    grover_proj
              │                      │                        │
              │                      └───────┬────────────────┘
              │                              ▼
              │                    ④ 合并底物原子特征
              │                    x_str + grover_proj → x_reaction
              │                              │
              └──────────┬───────────────────┘
                         ▼
                  ⑤ 双向交叉注意力
                  酶特征 ↔ 底物特征（互相交换信息）
                         │
                         ▼
                  ⑥ 池化 + 拼接
                  把所有特征压缩成 7 个固定长度向量，拼在一起
                         │
                         ▼
                  ⑦ 预测头 MLP
                  输出特异性分数（logits）
```

#### 详细版：SS.forward() 的完整数据流

下面是同一个过程的代码级详解。每一步标注了张量形状（B=batch_size=16）和对应的代码位置。

```
Step 1  序列通道（ss.py:96-97）
        输入: G.embedding [B, 1450, 1280]    ESM-2 预计算的蛋白质嵌入
        处理: protein_mlp (nn.Linear: 1280→128)
        输出: x_pro [B, 1450, 128]           酶的残基级特征
              （1450 = 最大酶长度，不足的位置用 0 填充）

Step 2  结构通道（structure.py + egnn.py）
        输入: 蛋白质原子坐标、底物原子坐标、边连接关系
        处理: Graph 模块（原子类型嵌入） → EGNN（3 层消息传递）
        输出:
          - x_str [B, 280, 128]     底物原子级结构特征
          - str_mean [B, 128]       整个复合物的全局结构特征（所有节点取平均）
              （280 = 最大底物原子数）

Step 3  分子通道 — 原子级（ss.py:107-110）
        输入: G.grover [B, 280, 2400]        GROVER 预计算的原子级嵌入
        处理: atom_feature_mlp (nn.Linear: 2400→128 + LayerNorm)
        输出: grover_proj [B, 280, 128]      底物原子级分子特征

Step 4  合并底物原子特征（ss.py:112）
        输入: x_str [B, 280, 128] + grover_proj [B, 280, 128]
        处理: atom_feature_aggregator（拼接成 [B, 280, 256] → MLP → 128）
        输出: x_reaction [B, 280, 128]       底物的综合原子级特征

        此时有两个关键变量：
          x_pro      [B, 1450, 128]  代表酶（来自 Step 1）
          x_reaction [B, 280, 128]   代表底物（来自 Step 2+3+4）

Step 5  双向交叉注意力（ss.py:116-119）
        方向1 enzyme_attention（酶 → 底物）：
          Q=x_pro, K=V=x_reaction
          输出: weighted_sum_pro [B, 1450, 128]
        方向2 reaction_attention（底物 → 酶）：
          Q=x_reaction, K=V=x_pro
          输出: weighted_sum_reaction [B, 280, 128]

Step 6  掩码均值池化（ss.py:127-137）
        把变长序列压缩成固定长度向量（忽略填充位置）：
          x_pro                 [B, 1450, 128] → 均值 → [B, 128]
          x_reaction            [B, 280, 128]  → 均值 → [B, 128]
          weighted_sum_pro      [B, 1450, 128] → 均值 → [B, 128]
          weighted_sum_reaction [B, 280, 128]  → 均值 → [B, 128]

Step 7  分子通道 — 全局级（ss.py:148-151）
        G.grover_mean [B, 4885] → feature_mlp → [B, 128]
        G.morgan      [B, 1024] → feature_mlp → [B, 128]

Step 8  拼接 + 预测（ss.py:153）
        拼接以下 7 个 [B, 128] 向量，得到 [B, 896]：
          ① x_pro              （酶均值）
          ② x_reaction          （底物均值）
          ③ weighted_sum_pro    （交叉注意力：酶视角）
          ④ weighted_sum_react  （交叉注意力：底物视角）
          ⑤ str_mean            （全局结构特征，来自 Step 2）
          ⑥ grover_mean_proj    （全局 GROVER 特征）
          ⑦ morgan_proj         （Morgan 指纹特征）

        specificity_header (MLP: 896→128→128→1)
        输出: logits [B, 1]（原始分数，未经 sigmoid）
```

### 4.2 序列通道详解：ESM 投影

```python
# ss.py:54
self.protein_mlp = nn.Linear(1280, 128)

# ss.py:96-97
x_pro = self.protein_mlp(G.embedding)    # [B*1450, 1280] → [B*1450, 128]
x_pro = x_pro.view(-1, 1450, 128)       # → [B, 1450, 128]
```

**做了什么**：ESM-2 给每个氨基酸一个 1280 维向量（太大了），用一个线性层压缩到 128 维。

**做了什么**：ESM-2 输出的每残基表示是 1280 维，但模型内部统一用 128 维（`hidden_dim`）。这个线性层做的就是降维：把 1280 维压缩到 128 维，丢弃冗余信息，保留对分类任务有用的部分。

**为什么是 1450**：配置中 `max_enzyme_length=1450`，不足的序列用 0 填充，超过 1000 氨基酸的蛋白质在特征生成阶段就被跳过了。

### 4.3 结构通道详解：SE(3)-GNN

#### 4.3.1 图是怎么构建的？

在数据进入模型之前，`EdgeConnection` 变换（[transforms.py:84-162](src/Datasets/Structure/transforms.py#L84-L162)）已经把蛋白质口袋+底物构建成了一张**图**：

```
图的构成：
├── 节点 = 原子（蛋白质口袋原子 + 底物原子）
│   ├── 蛋白质原子特征: [28维] = 元素(6) + 氨基酸类型(21) + 是否骨架(1)
│   └── 配体原子特征:   [~30维] = 元素(11) + 原子性质（维度取决于RDKit枚举值）
│
├── 边 = 原子之间的连接
│   ├── 共价键（底物内部的化学键）
│   └── k-NN 空间近邻（3D距离最近的 k=48 个原子对）
│       └── 包括跨分子边（蛋白质原子 ↔ 底物原子）——模型学习酶-底物交互的基础
│
└── 边特征 [39维]（当前配置）
    ├── 键类型 one-hot [6维]（单键/双键/三键/芳香键/共价其他/空间近邻）
    ├── 距离 RBF 编码 [32维]（num_r_gaussian=32 个高斯核）
    └── 是否跨分子 [1维]（蛋白质-底物之间的边标记为 1）
```

这里最关键的是**跨分子边**——蛋白质原子和底物原子之间的空间近邻边。正是这些边让 GNN 能学到蛋白质口袋和底物之间的 3D 相互作用。

#### 4.3.2 消息传递（Message Passing）

**核心思想**：每个节点（原子）初始只有自己的特征。通过多轮消息传递，每个节点逐步获取邻居的信息，从而"感知"自己所在的局部化学环境。

具体过程：
1. **计算消息**：对每条边，把两端节点特征和边特征拼接，通过 MLP 生成一条消息
2. **聚合消息**：每个节点把所有指向自己的消息求和
3. **更新特征**：用聚合后的消息和节点自身特征一起通过 MLP，更新节点特征

经过 3 层消息传递后，每个节点的特征向量中已经编码了 3 跳范围内所有邻居的信息。对于蛋白质-底物界面上的原子，这意味着它们的特征中包含了与对方分子原子的相互作用信息。

```
EGNN 单层消息传递（egnn.py:24-56）：

Step 1: 计算边消息
  对每条边 (src → dst)：
  mij = edge_mlp([edge_feat, h_dst, h_src])    # 拼接边特征+两端节点特征 → MLP
                                                 # [39+128+128=295] → [128]（当前配置）

Step 2: 注意力门控（可选）
  eij = mij * sigmoid(linear(mij))              # 给消息加权重（哪些邻居更重要）

Step 3: 聚合到节点
  mi = scatter_sum(eij, dst)                    # 每个节点收集所有指向它的消息并求和

Step 4: 更新节点
  h = h + node_mlp([mi, h])                     # 残差连接：旧特征 + 变换后的新信息
```

整个 EGNN 堆叠 **3 层**（[egnn.py:72-78](src/Models/Structure/egnn.py#L72-L78)），每层之间还有一个额外的残差连接（[egnn.py:84](src/Models/Structure/egnn.py#L84)）。

#### 4.3.3 "SE(3)-等变"到底是什么意思？

**SE(3)** = Special Euclidean group in 3D = 三维空间中的旋转 + 平移。

**等变（Equivariant）**的含义：如果你旋转输入分子，输出也跟着旋转（而不是完全变了）。

**诚实说明**（Codex 审计确认）：这个实现**实际上是不变的（invariant）而非等变的**。因为坐标在进入模型前已经被转换成了**距离 RBF 特征**（[transforms.py:109](src/Datasets/Structure/transforms.py#L109)），模型内部不处理坐标。距离在旋转下不变，所以输出也不变。这对我们的二分类任务来说**完全没问题**——旋转一个蛋白质-底物复合物，"能不能催化"的答案不应该改变。

#### 4.3.4 图输出怎么回到模型？

```python
# ss.py:87-91
def _get_graph_output(self, G, B):
    ligand_embedding = torch.zeros(B, 281, 128)          # 预分配空间（281是因为蛋白质节点用280做哨兵）
    output, (h, batch, index) = self.graph_net(G)         # GNN前向传播
    ligand_embedding[batch, index, :] = h                 # 把GNN输出的底物节点放回对应位置
    return output, ligand_embedding[:, :280, :]           # 截掉哨兵，返回 [B,280,128]
```

两个输出：
- `str_mean [B, 128]`：整张图所有节点的平均 → 代表"复合物整体印象"
- `x_str [B, 280, 128]`：底物每个原子的 GNN 编码 → 后面和 GROVER 融合

### 4.4 分子通道详解

#### 原子级：GROVER + GNN 融合

```python
# ss.py:107-112
grover = G.grover.view(-1, 280, 2400)                    # GROVER原子嵌入
grover_proj = atom_feature_mlp["grover"](grover)          # [B,280,2400] → [B,280,128]

x_reaction = atom_feature_aggregator([x_str, grover_proj]) # 拼接 → MLP
# [B,280,128] + [B,280,128] = [B,280,256] → MLP → [B,280,128]
```

**做了什么**：把 GNN 学到的"3D结构信息"和 GROVER 的"化学信息"融合在一起。每个底物原子现在同时知道"自己在3D空间中的环境"和"自己的化学性质"。

#### 全局级：GROVER mean + Morgan

```python
# ss.py:148-151
grover_mean_proj = feature_mlp["grover_mean"](G.grover_mean)  # [B,4885] → [B,128]
morgan_proj = feature_mlp["morgan"](G.morgan)                  # [B,1024] → [B,128]
```

这两个是**底物级别**的全局特征，不区分单个原子，只描述整个分子的化学性质。

### 4.5 交叉注意力详解（最核心的创新）

#### 4.5.1 什么是注意力？

注意力机制（Attention）的核心是**加权聚合**：不是对所有位置一视同仁地平均，而是根据相关性给不同位置分配不同权重。

它有三个输入：
- **Query（查询）**：当前位置想要获取什么信息（"我在找什么"）
- **Key（键）**：其他位置的索引信息（"我有什么"）
- **Value（值）**：其他位置的实际内容（"我的具体内容是什么"）

计算过程：Query 和每个 Key 做点积 → 得到相关性分数 → softmax 归一化为权重 → 用权重对 Value 做加权求和。结果是：与 Query 最相关的位置贡献最大，不相关的位置贡献接近 0。

#### 4.5.2 为什么是"双向"交叉注意力？

```
方向1: enzyme_attention（酶 → 底物）
  对酶的每个残基，计算它与底物所有原子的相关性权重，加权求和底物特征
  Query = x_pro [B, 1450, 128]              酶残基的特征向量
  Key = Value = x_reaction [B, 280, 128]    底物原子的特征向量
  输出: weighted_sum_pro [B, 1450, 128]     每个酶残基获得底物信息的加权聚合

方向2: reaction_attention（底物 → 酶）
  对底物的每个原子，计算它与酶所有残基的相关性权重，加权求和酶特征
  Query = x_reaction [B, 280, 128]          底物原子的特征向量
  Key = Value = x_pro [B, 1450, 128]        酶残基的特征向量
  输出: weighted_sum_reaction [B, 280, 128] 每个底物原子获得酶信息的加权聚合
```

**为什么需要双向？** 单向只有酶对底物的注意力权重，双向则同时建模底物对酶的注意力权重。催化反应依赖双向的特异性匹配——需要知道酶活性位点的哪些残基与底物相关，也需要知道底物的哪些原子与酶活性位点相关。

#### 4.5.3 多头注意力（8 heads）

多头的含义：把 128 维的向量拆成 8 组（每组 16 维），每组独立计算注意力。每个头可以学到不同的相关性模式——例如一个头可能关注电荷互补，另一个关注空间距离，另一个关注疏水性匹配。8 个头的结果拼接后再投影回 128 维，综合了多种角度的信息。

#### 4.5.4 填充掩码（Padding Mask）

因为不同样本的酶长度不同（有的500残基，有的1000残基），统一填充到1450。填充位置需要用**掩码**告诉注意力机制"这些位置是假的，别看"。

```python
# ss.py:127-128
reaction_mask = (~G.reaction_padding_mask).unsqueeze(-1).float()  # True=真实, False=填充
enzyme_mask = (~G.enzyme_padding_mask).unsqueeze(-1).float()
```

### 4.6 预测头

```python
# ss.py:140-153
embeddings = [
    x_pro,                    # [B, 128]  酶的均值表示
    x_reaction,               # [B, 128]  底物的均值表示
    weighted_sum_pro,         # [B, 128]  酶经过注意力后的表示
    weighted_sum_reaction,    # [B, 128]  底物经过注意力后的表示
    str_mean,                 # [B, 128]  3D结构的整体表示
    grover_mean_proj,         # [B, 128]  GROVER分子级表示
    morgan_proj               # [B, 128]  Morgan指纹表示
]
# 拼接: 7 × 128 = 896维

logits = specificity_header(embeddings)  # MLP: 896 → 128 → 128 → 1
```

这 7 个向量覆盖了模型所有信息源：原始酶表示、原始底物表示、交叉注意力增强后的酶/底物表示、3D 结构全局表示、以及两种分子指纹。拼接后通过 3 层 MLP 逐步压缩（896→128→128→1），最终输出一个标量。

**这个分数是 logit**（原始分数，可正可负），不是概率。要变成概率需要 `sigmoid(logit)`。训练时用 `BCEWithLogitsLoss`，它内部自动做 sigmoid。

---

## 五、训练循环详解

### 5.1 一个 epoch 发生了什么？

```
一个 epoch = 遍历整个训练集一次

┌─────────────────────────────────────────────────────────────┐
│ for batch in train_dataloader:    # 每次取16个样本           │
│   │                                                         │
│   ├─ training_step(batch):        # ss.py:163-166          │
│   │   logits, tag = self(batch)   # 前向传播 → 得到预测分数  │
│   │   loss = get_loss(batch, logits, "train")  # 计算损失   │
│   │   return loss                                           │
│   │                                                         │
│   ├─ loss.backward()              # Lightning 自动做        │
│   │   (反向传播：计算每个参数的梯度)                           │
│   │                                                         │
│   ├─ optimizer_step():            # ss.py:268-277          │
│   │   if epoch < 8:  lr = 线性插值(min_lr → lr)  # 预热    │
│   │   optimizer.step()            # 用梯度更新参数           │
│   │                                                         │
│   └─ (每2个batch才真正更新一次，因为 accumulate_grad=2)      │
│                                                             │
│ 每个epoch结束后:                                             │
│   for batch in val_dataloader:                              │
│     validation_step(batch)        # 前向传播 + 计算指标      │
│   validation_epoch_end()          # 汇总所有batch的AUC/AUPR │
│   ReduceLROnPlateau.step(aupr)    # 如果AUPR没提升，降低lr   │
│   ModelCheckpoint.check()         # 如果是最好的，保存模型    │
└─────────────────────────────────────────────────────────────┘
```

### 5.2 损失函数

```python
# ss.py:156-161
def get_loss(self, x, logits, stage):
    logits = logits.squeeze(-1)                              # [B,1] → [B]
    loss = (BCEWithLogitsLoss(reduction='none')(logits, label)  # 每个样本单独算loss
            * x.sample_weight                                   # 乘以权重
           ).mean()                                             # 取平均
    return loss
```

**BCEWithLogitsLoss**（Binary Cross-Entropy with Logits）：
- 二分类的标准损失函数
- 内部先做 `sigmoid(logits)` 变成概率，再和真实标签（0或1）比较
- 预测对了 → loss 小；预测错了 → loss 大

**sample_weight**：在当前配置中两个权重都是 1.0（[YAML:117-119](saved_model/model/run_0/complete-full-random-all-0-complex.yml#L117-L119)），实际上没有起作用。原始论文没有在 loss 层面处理类别不平衡（1:9.36），而是在数据层面处理。

### 5.3 优化器与学习率

**优化的基本概念**：

- **梯度（Gradient）**：损失函数对每个参数的偏导数。方向表示"增大这个参数 loss 会增还是减"，大小表示"影响有多强"
- **优化器（AdamW）**：根据梯度更新参数的算法。AdamW 维护每个参数的一阶矩（梯度均值）和二阶矩（梯度方差），让更新方向更稳定、自适应地调整每个参数的步长
- **学习率（lr=3e-4）**：控制每次参数更新的幅度。太大会导致 loss 震荡不收敛，太小会导致训练极慢
- **梯度裁剪（gradient_clip_val=8）**：当梯度范数超过 8 时，等比缩放所有梯度使其范数恰好为 8。防止某些 batch 产生异常大的梯度导致参数更新失控

#### 学习率变化过程

```
lr
│   3e-4 ───────────────────────╮
│                               │  sched_patience=10
│                               │  连续10个epoch AUPR没提升
│   1.5e-4 ─────────────╮      ▼  则 lr × factor（减半）
│                        │
│   7.5e-5 ────╮         ▼
│               ▼
│                  ... 直到 min_lr=5e-6
│  /
│ / warmup 阶段
│/  前 8 个 epoch 从 5e-6 线性增到 3e-4
│
└───────────────────────────────────────→ epoch
0   8
```

**为什么要 warmup（预热）？**

刚开始训练时，模型参数是随机的，梯度方向很不稳定。如果一上来就用大学习率，容易"暴走"。所以先用很小的学习率慢慢走（epoch 0 时 lr=5e-6），等模型"稳定下来"再提速。

**为什么要 ReduceLROnPlateau（高原降速）？**

训练到后期，模型已经接近最优了。如果还用大步子走，会在最低点附近来回震荡。自动检测"如果连续10个epoch都没改善"，就把学习率减半，让步子更小更精确。

### 5.4 梯度累积

```yaml
accumulate_grad_batches: 2
```

**原理**：不是每个 batch 都更新参数，而是累积 2 个 batch 的梯度再更新。效果等价于把 batch_size 翻倍（16×2=32），但不需要更多显存。

**原论文的有效batch size**：16（每GPU）× 4（GPU数）× 2（累积）= **128**。
我们在单卡上如果保持 batch_size=16，有效 batch 只有 16×1×2=**32**。这个差异需要注意。

---

## 六、配置文件解读

### 6.1 完整配置（带注释）

来源：[complete-full-random-all-0-complex.yml](saved_model/model/run_0/complete-full-random-all-0-complex.yml)

```yaml
# === 数据配置 ===
batch_size: 16              # 每批16个酶-底物对
max_substrate_length: 280   # 底物最多280个原子（超过跳过）
max_enzyme_length: 1450     # 酶最长1450个残基（超过跳过）

features:                   # 全局特征（底物级别）
  - 'morgan'               # Morgan指纹
  - 1024                   # 1024维
  - 'grover_mean'          # GROVER分子级嵌入
  - 4885                   # 4885维

atom_features:              # 原子级特征
  - 'grover'               # GROVER原子嵌入
  - 2400                   # 每原子2400维（经截断）

sample_weight: [1.0, 1.0]  # 复合物/纯配体的loss权重（当前相同）

# === 图变换配置 ===
transform:
  dist_noise: True          # 训练时给距离加噪声（数据增强）
  cutoff: 10                # RBF高斯核的最大距离值（非硬截断，代码中距离过滤已注释掉）
  num_r_gaussian: 32        # 距离RBF编码用32个高斯核
  k: 48                     # k近邻图，每个原子连接最近的48个

# === 模型配置 ===
model:
  hidden_dim: 128           # 所有隐藏层统一128维
  use_gnn: True             # 启用结构通道（GNN）
  feature_norm: False       # 特征投影后不做LayerNorm

  header:                   # 预测头
    act_fn: 'relu'          # 激活函数
    norm: True              # 注意：代码实际不识别True，等于无归一化
    num_layers: 3           # 3层MLP (896→128→128→1)

  graph:                    # GNN配置
    name: egnn              # 用EGNN（不是GAT）
    hidden_dim: 128         # GNN隐藏维度
    num_layers: 3           # 3层消息传递
    attention: True         # 边注意力门控
    residual: True          # 残差连接（注意：此标志在代码中未被检查，残差是硬编码的）
    norm: layer             # LayerNorm
    mode: complex           # 只用蛋白质-底物复合图

  cross_attention:          # 交叉注意力
    n_head: 8               # 8个注意力头
    dropout: 0.9            # 极高的dropout！见下方解释

# === 训练配置 ===
training:
  optimizer: 'adamW'            # AdamW优化器
  lr: 3.0e-04                   # 基础学习率
  min_lr: 5.0e-06               # 最小学习率（warmup起点 + plateau下限）
  warmup_epochs: 8              # 前8轮线性预热
  weight_decay: 0               # L2正则化（注意：代码未传入AdamW，此值被忽略）
  gradient_clip_val: 8          # 梯度裁剪阈值
  sched_factor: 0.5             # 高原时lr减半
  sched_patience: 10            # 10个epoch没提升才降lr
  val_frequency: 1              # 每个epoch验证一次
  accumulate_grad_batches: 2    # 梯度累积2步
```

### 6.2 异常值：dropout=0.9

这个值**非常高**（通常 0.1-0.5）。0.9 意味着交叉注意力中 **90% 的权重被随机丢弃**。

可能的原因：
- 训练数据 32 万条，模型仅 185 万参数 → 不太需要正则化，但作者可能实验发现高 dropout 对注意力有正则效果
- 也可能是超参搜索的结果，或者一个有意为之的"极端正则化"策略

**对微调的启示**：P450 数据量远小于原始数据（几千 vs 32万），但 0.9 的 dropout 可能过于激进。Codex 建议微调时降到 **0.1-0.3**。

### 6.3 Codex 发现的代码-配置不一致

| 问题 | 位置 | 影响 |
|------|------|------|
| `weight_decay: 0` 在 YAML 中定义，但 `AdamW` 创建时没传入 | [ss.py:249](src/Models/ss.py#L249) | AdamW 使用默认 weight_decay=0.01 |
| `header.norm: True` 不被 MLP 识别（只认 `'layer'`/`'batch'`） | [common.py:164](src/Models/common.py#L164) | 预测头实际无归一化 |
| `need_weights=True` 返回注意力矩阵但被丢弃 | [ss.py:116](src/Models/ss.py#L116) | 浪费约 52MB 显存 |

---

## 七、模型参数量与显存估算

### 7.1 参数量分解（Codex 估算）

| 模块 | 参数量 | 占比 |
|------|--------|------|
| `protein_mlp`（ESM投影） | ~164K | 8.9% |
| `graph_net`（GNN编码器） | ~320K | 17.3% |
| `atom_feature_mlp_dict`（GROVER投影） | ~307K | 16.6% |
| `atom_feature_aggregator`（GNN+GROVER融合） | ~33K | 1.8% |
| `enzyme_attention` + `reaction_attention` | ~132K | 7.1% |
| `feature_mlp_dict`（全局特征投影） | ~757K | 40.9% |
| `specificity_header`（预测头） | ~131K | 7.1% |
| **总计** | **~1.85M** | 100% |

**对比参考**：GPT-2 有 1.5 亿参数，我们的模型只有 185 万 —— 非常小。FP32 权重只占 **7.4 MB**。

### 7.2 单 batch 显存估算（batch_size=16）

| 张量 | 形状 | 大小 |
|------|------|------|
| ESM 嵌入输入 | [16, 1450, 1280] | ~119 MB |
| GROVER 原子输入 | [16, 280, 2400] | ~43 MB |
| 投影后酶张量 | [16, 1450, 128] | ~12 MB |
| 注意力权重（被浪费的） | 2 × [16, 1450, 280] | ~52 MB |
| 图张量（变长） | 通常 < 50 MB | ~50 MB |
| 梯度 + 优化器状态 | 约等于参数量 ×3 | ~22 MB |
| **估计总计** | | **~300 MB** |

**注意**：上表仅统计主要张量的静态大小，是**下限估计**。实际训练时还有大量 autograd 激活缓存、注意力工作空间等临时张量，峰值显存会显著高于此值。Codex 建议**从 batch_size=8 开始**，确认不 OOM 后再增大。用 `precision='16-mixed'`（半精度训练）可以进一步节省约 30-50% 显存。

---

## 八、微调可行性分析

### 8.1 什么是微调？

**从头训练 vs 微调的区别**：

| | 从头训练 | 微调 |
|---|---------|------|
| 参数初始值 | 随机初始化 | 加载预训练权重 |
| 数据需求 | 大量（32万条） | 少量（几千条） |
| 计算资源 | 多GPU、2-4周 | 单GPU、几小时到几天 |
| 学习率 | 较大（3e-4） | 较小（1e-5 ~ 3e-5） |
| 风险 | 数据不够会欠拟合 | 数据不够会过拟合 |

微调的核心思想是**迁移学习**（Transfer Learning）：预训练模型已经在大数据上学到了通用的生化知识（如"什么样的化学结构倾向于被酶催化"），微调只需要在此基础上学习 P450 特异的知识。

具体步骤：
1. 加载已训练好的模型权重（预训练检查点）
2. **冻结**大部分层的参数（设置 `requires_grad=False`，反向传播时不计算梯度、不更新）
3. 只训练少数任务特异层的参数
4. 用 P450 专属数据训练

### 8.2 哪些层冻结？哪些层训练？

```
模型各层与冻结策略：

┌─────────────────────────────────────────────────────────────┐
│  ❄️ 冻结区（通用知识，不需要改）                              │
│                                                             │
│  protein_mlp          "理解蛋白质序列的能力"     ~164K       │
│  graph_net (EGNN)     "理解3D结构的能力"         ~320K       │
│  atom_feature_mlp     "理解分子化学的能力"        ~307K      │
│  feature_mlp_dict     "理解全局分子特征的能力"     ~757K     │
│                                                             │
│  这些层学到的是"通用的生化知识"，对所有酶家族都适用           │
│  冻结它们 = 保留通用知识，防止小数据集把它们"遗忘"            │
├─────────────────────────────────────────────────────────────┤
│  🔥 训练区（任务特异，需要针对P450调整）                      │
│                                                             │
│  enzyme_attention     "酶问底物什么问题"          ~66K       │
│  reaction_attention   "底物问酶什么问题"          ~66K       │
│  specificity_header   "怎么做最终判断"            ~131K      │
│  atom_feature_aggr.   "怎么融合结构和化学信息"     ~33K      │
│                                                             │
│  这些层负责"酶-底物配对的具体判断逻辑"                       │
│  P450的催化机制和其他酶不同，需要重新学习这个判断逻辑         │
└─────────────────────────────────────────────────────────────┘

冻结后只需训练 ~296K 参数（原来的 16%）→ 显存需求大幅降低
```

### 8.3 微调需要什么数据？

| 需求 | 来源 | 状态 |
|------|------|------|
| P450 酶序列 + ESM 嵌入 | ESIBank 367个P450酶 | 需确认是否可从 Google Drive 获取 |
| 底物 SMILES + GROVER/Morgan | ESIBank ~12,329个酶-底物对 | 需确认 |
| 对接复合物 PDB | ESIBank 用 AlphaFold 结构 | 需确认格式兼容性 |
| 正/负标签 + 负样本设计 | 需重新设计（匹配难度的负样本） | **关键任务** |
| 训练/验证/测试划分 | 需新设计（避免酶泄漏） | **关键任务** |

### 8.4 微调脚本原型（Codex 提供的参考）

```python
# 这是伪代码/参考，不能直接运行
import pytorch_lightning as pl
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint
from Datasets.brenda import Singledataset
from Models.ss import SS
from utils import load_config

def freeze_backbone(model):
    """冻结通用知识层，只训练任务特异层"""
    for module in [model.graph_net, model.protein_mlp,
                   model.atom_feature_mlp_dict, model.feature_mlp_dict]:
        for param in module.parameters():
            param.requires_grad = False

config = load_config("config_p450.yml")
model = SS.load_from_checkpoint("best-checkpoint.ckpt", config=config)
freeze_backbone(model)

dm = Singledataset(config)

callbacks = [
    ModelCheckpoint(monitor="aupr/val", mode="max", save_top_k=1),
    EarlyStopping(monitor="aupr/val", mode="max", patience=10),
]

trainer = pl.Trainer(
    accelerator="gpu", devices=1,
    max_epochs=30,
    precision="16-mixed",        # 半精度训练，省显存
    gradient_clip_val=8,
    accumulate_grad_batches=2,
    callbacks=callbacks,
)

trainer.fit(model, datamodule=dm)
trainer.test(ckpt_path="best", datamodule=dm)
```

### 8.5 微调的关键风险

| 风险 | 说明 | 应对 |
|------|------|------|
| **过拟合** | P450 数据少（几千 vs 原来32万），容易记住训练集 | 冻结大部分层 + 早停 + 降低 dropout 到 0.1-0.3 |
| **酶泄漏** | 如果训练/测试用了相同的酶，AUC 会虚高 | 严格按酶划分训练/测试集 |
| **负样本设计** | 我们已验证负样本难度是 AUC 差距的主因 | 使用 EC 层级匹配的负样本（≈难度3） |
| **小验证集不稳定** | 少量样本的 AUC 波动大 | 使用多折交叉验证或 bootstrap |
| **有效 batch size 差异** | 原文128，我们32 | 可能需要调整 lr（线性缩放规则） |

### 8.6 论文先例：HaloS 微调

论文中对卤化酶（Halogenase）做过微调：
- 数据量：~3,300 个酶-底物对
- 结果：Top-1 准确率 91.7%
- 说明：小数据微调在这个模型上是可行的

---

## 九、关键文件索引

| 文件 | 角色 | 重要方法/类 |
|------|------|------------|
| [src/Models/ss.py](src/Models/ss.py) | **主模型** | `SS`, `forward()`, `training_step()`, `configure_optimizers()` |
| [src/Models/Structure/egnn.py](src/Models/Structure/egnn.py) | GNN 层 | `EGNN`, `EnBaseLayer` |
| [src/Models/Structure/structure.py](src/Models/Structure/structure.py) | 图构建 | `Graph` |
| [src/Models/common.py](src/Models/common.py) | 工具模块 | `MLP`, `GaussianSmearing`（注意：`Aggregator` 在 ss.py:17） |
| [src/Datasets/brenda.py](src/Datasets/brenda.py) | 数据模块 | `Singledataset` |
| [src/Datasets/data_representer.py](src/Datasets/data_representer.py) | 数据表示 | `StructureSequence`, `Reaction` |
| [src/Datasets/Structure/transforms.py](src/Datasets/Structure/transforms.py) | 图变换 | `FeaturizeProteinAtom`, `FeaturizeLigandAtom`, `EdgeConnection` |
| [src/Datasets/create_features.py](src/Datasets/create_features.py) | 特征生成 | `generate_esm_embedding` |
| [src/main_testing.py](src/main_testing.py) | 推理入口 | `generate_prediction` |
| [saved_model/.../complete-full-random-all-0-complex.yml](saved_model/model/run_0/complete-full-random-all-0-complex.yml) | 配置文件 | 所有超参数 |

---

## 十、术语表

| 术语 | 英文 | 通俗解释 |
|------|------|---------|
| 前向传播 | Forward Pass | 数据从输入到输出走一遍模型 |
| 反向传播 | Backpropagation | 从 loss 反推每个参数"该往哪调" |
| 梯度 | Gradient | 参数的调整方向和幅度 |
| 学习率 | Learning Rate | 每次调整参数走多远 |
| epoch | Epoch | 把整个训练集看一遍 |
| batch | Batch | 一次看的样本数（这里16个） |
| logit | Logit | 模型输出的原始分数（未经sigmoid） |
| sigmoid | Sigmoid | 把任意实数压缩到 0-1 之间的函数 |
| BCE Loss | Binary Cross-Entropy | 二分类损失函数 |
| dropout | Dropout | 训练时随机"关掉"一部分连接，防止过拟合 |
| 残差连接 | Residual Connection | 跳跃连接，让信息可以"绕过"某一层 |
| 消息传递 | Message Passing | GNN中节点之间交换信息的过程 |
| 注意力 | Attention | 让模型学会"看重点"的机制 |
| 掩码 | Mask | 标记"哪些位置是真实数据，哪些是填充" |
| 微调 | Fine-tuning | 在已训练好的模型上用新数据继续训练 |
| 冻结 | Freeze | 锁定某些层的参数不更新 |
| 检查点 | Checkpoint | 训练中间保存的模型权重文件 |
| 混合精度 | Mixed Precision | 用 FP16 加速计算 + 省显存 |

---

## 附录A：Codex 审核勘误

本文档经 Codex 技术审核（2026-03-09），以下是已修正的关键问题：

| 原始错误 | 修正 | 严重性 |
|---------|------|:---:|
| 边特征写成31维 | 实际39维（6+32+1，因为 num_r_gaussian=32） | 高 |
| EGNN消息输入287维 | 实际295维（39+128+128） | 高 |
| cutoff=10 描述为"空间截断距离" | 实际只是RBF高斯核的stop值，距离过滤代码已被注释 | 高 |
| 训练时"可能用纯配体回退" | 当前配置 `full_data: False`，训练也只用复合物 | 高 |
| 显存估算~300MB描述为"总计" | 改为"下限估计"，实际峰值显著更高 | 中 |
| Aggregator 标注在 common.py | 实际定义在 ss.py:17 | 低 |

**未修正但需注意的细节**：
- PyG 批处理后，`G.embedding` 的实际 runtime 形状是 `[B*1450, 1280]`（扁平化），在 `ss.py:97` 用 `.view()` 还原为 `[B, 1450, 1280]`。流程图中直接写 `[B, 1450, 1280]` 是概念性简化。
- `follow_batch=['ligand_index']` 是数据加载的关键设置（[brenda.py:40](src/Datasets/brenda.py#L40)），没有它 GNN 的 batch 索引会缺失导致崩溃。
- ESM 特征生成跳过 >1000 氨基酸的序列（[create_features.py:296](src/Datasets/create_features.py#L296)），但 runtime 填充到 1450 维——大部分 1450 维空间是 padding。
- `optimizer: 'rmsport'` 是代码中的拼写错误（[ss.py:247](src/Models/ss.py#L247)），如果在 YAML 中写正确的 `'rmsprop'` 会匹配不到。
- AUC/AUPR 是从 raw logits 计算的（不经过 sigmoid），这对排序指标是正确的。
