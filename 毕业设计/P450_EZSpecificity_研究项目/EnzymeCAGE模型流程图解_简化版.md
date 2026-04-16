# EnzymeCAGE 模型完整流程图解（新手版）

> 本文档面向已经读过《EZSpecificity 模型流程图解_简化版.md》的读者。
> 不涉及具体函数名和代码，只讲每一步在做什么、为什么这样做、和 EZSpecificity 的区别。
> 源代码位置：`d:/EZSpecificity_Project/毕业设计/EnzymeCAGE/`
> 论文：Liu et al., *Nature Catalysis*, 2026. "A geometric foundation model for enzyme retrieval with evolutionary insights"

---

## 一、模型要解决什么问题

```
已知：一个酶（蛋白质）  +  一个完整的化学反应（不是单个底物！）
      反应的写法：A + B >> C + D
                 （底物）    （产物）

问题：这个酶能不能催化这个反应？

模型输出：一个分数
  分数 → sigmoid → 概率（0~1）
  概率 > 0.5 → 模型认为"能催化"
  概率 < 0.5 → 模型认为"不能催化"
```

**和 EZSpecificity 的核心区别**：

| 对比项 | EZSpecificity | EnzymeCAGE |
|---|---|---|
| 预测目标 | 酶 + **一个底物** → 能否催化 | 酶 + **完整反应（底物→产物）** → 能否催化 |
| 输入 | enzyme, substrate | enzyme, substrate, product |
| 任务 | 特异性预测（单向底物筛查）| 酶检索（给反应找最匹配的酶）|

EnzymeCAGE 明确建模"反应是从什么变成什么"，而不只是"酶碰上这个分子"。对很多催化任务这很关键——同一个底物可能走不同反应通路，产物才决定催化方向。

---

## 二、模型需要什么输入数据

EnzymeCAGE 用了**六种**输入特征（EZSpecificity 只用四种）：

```
┌──────────────────────────────────────────────────────────────────┐
│                    一个样本的六种描述                               │
│                                                                  │
│  ① 酶序列（ESM 嵌入）                                             │
│     全序列逐残基编码 (≈1152 或 1280 维/残基)                       │
│     和 EZSpecificity 一样                                         │
│                                                                  │
│  ② 酶口袋（GVP 残基图）★★★ 和 EZSpecificity 最大区别             │
│     pocket = 预测结构中距活性位点较近的残基（默认 8Å）             │
│     不是原子级！以"残基"为节点                                     │
│     每个残基记录 φ/ψ/ω 二面角 + 侧链方向向量                      │
│     边 = 残基之间的 kNN 空间连接                                  │
│                                                                  │
│  ③ 底物 3D 构象（SchNet 原子图）                                  │
│     RDKit 生成的 3D 立体构象（不是对接！是底物本身的最低能量构象） │
│     每个原子有 3D 坐标 + 原子序数                                 │
│                                                                  │
│  ④ 产物 3D 构象（SchNet 原子图）★★★ EZSpecificity 没有           │
│     同③，但针对反应的产物分子                                     │
│                                                                  │
│  ⑤ 反应中心标注（reacting center）★★★ EZSpecificity 没有         │
│     通过原子映射 AAM 识别反应中真正变化的原子                     │
│     每个原子一个 0/1 标志：1=这个原子是反应中心                    │
│                                                                  │
│  ⑥ 反应指纹（DRFP）                                               │
│     differential reaction fingerprint，2048 位 0/1 向量           │
│     描述"从底物到产物发生了什么化学变化"（相对 Morgan 更高级）     │
│                                                                  │
└──────────────────────────────────────────────────────────────────┘
```

**为什么这么多输入？** EnzymeCAGE 的野心是做酶检索的"基础模型"（foundation model），所以把所有可能有用的信息都喂进去：
- 酶的序列 → 身份
- 酶的口袋结构 → 空间几何
- 底物 + 产物的 3D 结构 → 反应前后的分子形状
- 反应中心 → 哪些原子在反应中关键
- 反应指纹 → 反应的化学类型

---

## 三、数据预处理：把原始数据变成数字

EnzymeCAGE 的预处理比 EZSpecificity 复杂得多，因为多了产物、反应中心、3D 构象等。

```
原始数据                        工具                        数字化特征
────────────────────────────────────────────────────────────────────

① 酶序列                       ESM-C 600M 或              每个残基 → 1152 或 1280 维向量
"MATLK..."                     ESM-2 650M                  （和 EZSpecificity 方式一致）

                                   │
                                   ▼

② pocket 残基列表               AlphaFill（预测结构）+      每个 pocket 残基:
（从预测结构提取）              提取 pocket (8Å)           - node_s (标量) 6 维
                                                          - node_v (向量) 3×3 方向向量
                                                          - edge_s (标量) 32 维（距离RBF）
                                                          - edge_v (向量) 1×3（相对方向）

                                   │
                                   ▼

③④ 底物/产物 SMILES           RDKit 生成 3D 构象          每个分子:
                               (ETKDGv2 算法，MMFF优化)    - 原子类型 (integer, 9 类)
                                                          - 原子 3D 坐标 (N×3)
                                                          （底物和产物分开处理）

                                   │
                                   ▼

⑤ 反应 SMILES                  AAM 原子映射                每个原子一个 0/1 标志：
"A>>B"                         （RXNMapper）               - 是否是反应中心
                                                          （底物和产物各自标注）

                                   │
                                   ▼

⑥ 反应 SMILES                  DRFP 库                    2048 位 0/1 向量
"A>>B"                         (differential              描述"分子图的差异"
                                reaction fingerprint)     (产物减去底物的子结构变化)
```

**重要区别：pocket 怎么提取？**

- **EZSpecificity**: 从对接复合物 PDB 里，选底物原子 10Å 范围内的蛋白质原子。**原子级**，结果是一张"蛋白原子 + 底物原子"的混合图。
- **EnzymeCAGE**: 从预测结构（AlphaFill 或 AlphaFold）里，用预定义方法（如 radius 或 固定 N 个）圈出 pocket **残基**。**残基级**。然后为每个残基计算二面角 + 侧链方向作为 GVP 输入。

**什么是 GVP？**

GVP = Geometric Vector Perceptron，**几何向量感知机**。它是一种专为处理 3D 几何设计的神经网络：

- 普通线性层只处理**标量**（数字）：距离=5.0Å、角度=60°
- GVP 同时处理**标量**和**向量**：距离=5.0Å（标量）、方向=(0.3, 0.5, -0.8)（向量）

**GVP 的关键性质：SE(3) 等变性**。你旋转整个分子 90°，GVP 算出的向量也跟着旋转 90°，但标量不变。这让模型"理解几何"而不只是"记住坐标"。

```
举例：一个残基的 node 特征

标量部分 node_s (6 维)：
  [sin(φ), cos(φ), sin(ψ), cos(ψ), sin(ω), cos(ω)]
   —— 记录二面角

向量部分 node_v (3 个 3D 向量)：
  [ CA→N   方向单位向量  ]  ← 骨架 N 端方向
  [ CA→C   方向单位向量  ]  ← 骨架 C 端方向
  [ CA→Cβ  方向单位向量  ]  ← 侧链方向

这 9 个数字加起来描述了"残基在 3D 空间中的朝向"
```

---

## 四、模型内部：从特征到预测

### 4.0 总览：完整架构图

```
=============================================================================
  阶段一：四条通道并行处理
=============================================================================

  ① 酶序列通道              ② 酶口袋通道             ③ 底物分子通道      ④ 产物分子通道
  （残基级全序列）          （残基级 3D 几何）        （原子级 3D）      （原子级 3D）

  ESM 全序列嵌入            pocket 残基集合           底物原子+坐标      产物原子+坐标
  [N_seq, 1152]             [N_pocket, ...]                │                │
       │                         │                         ▼                ▼
       │                         ▼                    SchNet (6 层)    SchNet (6 层)
       │                    GVP encoder (3 层)         ┌─ 原子级嵌入  ┌─ 原子级嵌入
       │                         │                     │   [Ns, 128]  │   [Np, 128]
       │                         ▼                     │              │
       │                    输出 [N_pocket, 128]       │              │
       │                         │                     │              │
       │                         ▼                     │              │
       │                    ★ 关键融合                  │              │
       │                    concat(GVP, ESM_pocket)    │              │
       │                      → [N_pocket, 1280]       │              │
       │                    （即：每个 pocket 残基      │              │
       │                     同时有几何表示+序列表示）   │              │
       │                         │                     │              │
       │                         ▼                     │              │
       │                    pocket 自注意力（8 头）      │              │
       │                    ★ 距离偏置（近权大）          │              │
       │                    → enz_node [N_pocket, 512] │              │
       │                         │                     │              │
       │                         │                     ▼              ▼
       │                         │            ★★★ 4.3 反应内部注意力（可选）
       │                         │            substrate ↔ product 交叉注意力
       │                         │            用 reacting_center 加权
       │                         │                     │              │
       │                         └──────────┬──────────┘              │
       │                                    │                         │
       │                                    ▼                         │
=======│===================================================================
  阶段二：几何增强交叉注意力（4 路）
=======│===================================================================
       │                         │
       │                ⑤ 计算 interaction_weight
       │                   pocket_weight × substrate_weight
       │                   （pocket 几何中心距离 × reacting center 权重）
       │                         │
       │                         ▼
       │                ⑥ 4 路交叉注意力（每路都用 interaction_weight 做 attention bias）
       │                   ┌─ 酶→底物 (enz_subs_output)
       │                   ├─ 底物→酶 (subs_enz_output)
       │                   ├─ 酶→产物 (enz_prod_output)
       │                   └─ 产物→酶 (prod_enz_output)
       │                         │
       │                         ▼
       │                   4 个 [128] 向量 → 拼接成 [512]
       │                   (interaction_fused)
       │                         │
=======│===================================================================
  阶段三：拼接 + 预测
=======│===================================================================
       │                         │
       ▼                         │
  ESM 全序列均值 [1152]          │
       │                         │
       ▼                         ▼
  ⑦ 末端大拼接
  interaction_fused [512]  +  ESM 全局 [1152]  +  DRFP [2048]  = [3840 ~ 3712]
                         │
                         ▼
  ⑧ 预测头（3 层大 MLP）
  3840 → 2048 → 1024 → 1
  （每层 BatchNorm + Dropout + Linear + LeakyReLU）
                         │
                         ▼
  输出：logits（预测分数）
  → sigmoid → 概率
```

**对比 EZSpecificity 总览（回顾）**：
- EZSpecificity 末端拼接 7 个 128 维 = 896 维
- EnzymeCAGE 末端拼接 512 + 1152 + 2048 ≈ 3712 维
- EnzymeCAGE 的预测头大得多（3840→2048→1024→1），因为要处理高得多的输入维度

---

### 4.1 酶序列通道：和 EZSpecificity 基本一样

```
输入：ESM 全序列嵌入（每个残基 1152 或 1280 维）

这里 EnzymeCAGE 比 EZSpecificity 还简单——ESM 嵌入不经过任何线性层降维，
直接保留原始 1152/1280 维，在两个地方使用：

  （a）pocket 残基那部分的 ESM 嵌入 → 和 GVP 输出拼接（见 4.2）
       这是"让每个 pocket 残基同时有几何信息和序列信息"

  （b）整条序列的 ESM 嵌入 → 末端直接拼接到预测头（见 4.7）
       这是"给模型一个全局的酶身份信号"
```

EnzymeCAGE 没有像 EZSpecificity 那样把 ESM 1280 维压到 128 维。作者选择**保留完整维度**，把降维交给末端预测头去做。好处是信息不丢失，代价是预测头更大。

---

### 4.2 酶口袋通道：残基级 GVP（核心创新之一）

这是 EnzymeCAGE 和 EZSpecificity 最本质的区别：**EZSpecificity 把口袋当原子集合，EnzymeCAGE 把口袋当残基集合**。

#### 4.2.1 建图：pocket 残基图

```
把 pocket 残基（~30~50 个）建成一张图：

  节点 = pocket 残基本身（不是原子！）
  边 = 残基之间的 kNN 空间近邻（比如每个残基连接最近的 30 个其他残基）

每个节点的初始特征是 (标量, 向量) 的 tuple：
  标量 node_s [6 维]：φ/ψ/ω 三个二面角的 sin/cos 编码
  向量 node_v [3 个 3D 向量]：骨架和侧链方向

每条边的初始特征也是 (标量, 向量) tuple：
  标量 edge_s [32 维]：两残基 Cα 距离的 RBF（径向基函数）展开
  向量 edge_v [1 个 3D 向量]：两残基 Cα 之间的单位方向向量
```

**为什么用残基级而不是原子级？**

作者认为：
- 酶的活性位点本质上由**残基**定义（哪个 His 是催化三元组？哪个 Ser 是亲核基团？）
- 残基级粒度更大，计算更快（~30 节点 vs ~200 原子）
- GVP 能处理"残基的方向"，比"单个原子的位置"更符合生物学直觉

**代价**：失去了原子级的精细接触信息（例如某个底物原子离某个具体骨架氧原子多近）。

#### 4.2.2 GVP encoder：SE(3) 等变消息传递

GVP 在图上做消息传递，和 EGNN 类似但更强：

```
GVP 消息传递（一轮）：

  1. 收集邻居信息（和 EGNN 类似）
     对每个节点的每个邻居，拼接：
       (节点的 s, 节点的 V) + (邻居的 s, 邻居的 V) + (边的 s, 边的 V)

  2. GVP 处理（关键！）
     通过 GVP 层计算消息
     GVP 内部：
       - 对标量 s 做普通线性层 + ReLU
       - 对向量 V 做"向量线性层"（矩阵左乘向量，保持 SE(3) 等变性）
       - 标量 s 可以"看到"向量 V 的范数（长度）作为额外特征
         即：向量部分的信息可以"告知"标量部分
       - 向量 V 被标量 s 缩放（gating）
         即：标量部分可以"控制"向量部分的强度

  3. 汇总更新
     所有邻居的消息求和 → 新的 (s, V) → 残差连接到原节点特征
```

做 3 轮消息传递后，每个 pocket 残基的特征向量融合了周围残基的几何信息。

**GVP 相比 EGNN 的优势**：
- EGNN 只处理标量（距离）→ 只能学"残基 A 离残基 B 有多近"
- GVP 同时处理标量和向量 → 能学"残基 A 的侧链朝哪个方向、离残基 B 的口袋入口有多近且朝向正确"

#### 4.2.3 GVP 输出

```
GVP 最后一层输出：每个 pocket 残基一个 128 维向量（纯标量，向量部分被转换成标量）

然后 EnzymeCAGE 做一件关键的事：
  ★ 把每个 pocket 残基的 GVP 向量 [128] 和它对应的 ESM 嵌入 [1152] 拼接
  ★ 得到 [1280] 维 = GVP（几何）+ ESM（序列）

这样每个 pocket 残基同时有两份信息：
  - GVP 部分：我在口袋里朝哪个方向、和邻居几何关系如何
  - ESM 部分：我是哪个氨基酸、在全序列的什么上下文

此后的 pocket 自注意力和交叉注意力都用这个 1280 维的联合表示。
```

EZSpecificity 的 pocket 只有"结构通道的原子特征"这一份信息；EnzymeCAGE 的 pocket 残基**天然就是结构+序列融合的**。

---

### 4.3 底物/产物分子通道：SchNet

EnzymeCAGE 用 **SchNet**（一个经典的 3D 分子神经网络）处理底物和产物的 3D 构象。

#### 4.3.1 什么是 SchNet？

SchNet 和 EGNN/GVP 一样是消息传递网络，但有两个特点：
- **连续滤波卷积（Continuous-filter Convolution）**：边的权重是距离的连续函数（用 RBF 展开 + MLP），不是离散的化学键类型
- **轴向对称**：只看距离不看方向（像 EGNN 一样 SE(3) 不变，但不等变）

```
SchNet 处理流程：

  输入：每个原子的原子序数（1=H, 6=C, 7=N, 8=O...）+ 3D 坐标

  1. 原子嵌入：原子序数 → 128 维初始特征（查 embedding 表）
  2. radius graph：每个原子连接 10Å 以内的所有邻居（包括非键连的）
  3. 6 层 InteractionBlock（消息传递）
     每层内部：
       - 边权重 W = MLP(RBF(距离)) × 平滑 cutoff
       - 消息 = 邻居特征 × W
       - 汇总邻居消息 + 当前节点特征 → 新特征
  4. 输出：每个原子 128 维 + 整个分子 128 维（scatter mean 得到）
```

**为什么用 SchNet 而不是 GROVER / Morgan？**

EZSpecificity 用的是**预训练的分子编码器**（GROVER）和**规则指纹**（Morgan）——这些是 2D 分子图特征，不看 3D 坐标。

EnzymeCAGE 认为反应是在 3D 空间里发生的，底物和产物的 3D 构象（尤其是反应过渡态附近的形状）是重要信息。所以用 **SchNet（3D 分子编码器）** 而不是 2D 编码器。

#### 4.3.2 底物-产物反应内部注意力（可选）

如果配置 `rxn_inner_interaction=True`，模型会在进入酶-分子交叉注意力之前，**先让底物和产物互相看一眼**：

```
底物 SchNet 输出 [Ns, 128]      产物 SchNet 输出 [Np, 128]
        │                             │
        ▼                             ▼
  ┌───────────────────────────────────────────┐
  │  底物-产物交叉注意力（CrossAttention）       │
  │                                           │
  │  用 reacting_center 做 attention bias:    │
  │    反应中心原子 → 权重高                    │
  │    非反应中心原子 → 权重低                  │
  │                                           │
  │  方向1：底物 query，产物 key/value         │
  │   → 增强后的底物表示                        │
  │                                           │
  │  方向2（可选）：产物 query，底物 key/value  │
  │   → 增强后的产物表示                        │
  └───────────────────────────────────────────┘
```

**为什么要先做这一步？**

因为反应前后分子的**变化位置**比分子本身更重要。通过反应中心加权的交叉注意力，模型能学到"底物的哪些原子对应产物的哪些原子、这对原子在反应中发生了什么变化"。

这是 EZSpecificity 完全没有的建模——EZSpecificity 只看底物，不知道反应变化。

---

### 4.4 pocket 自注意力：口袋内部信息传递

pocket 残基在进入酶-底物交叉注意力之前，**先做一次内部自注意力**：

```
pocket 残基联合特征 [N_pocket, 1280]  （GVP + ESM 拼接）
        │
        ▼
  ┌───────────────────────────────────────────┐
  │  Multi-Head Attention（8 头，512 维）      │
  │                                           │
  │  Q = K = V = pocket 残基特征              │
  │    → 残基两两之间计算注意力               │
  │                                           │
  │  ★ 距离偏置（pocket_inner_interaction=True）│
  │    attn_bias = 1 - distance / 30          │
  │    距离近 → bias 大 → 注意力强             │
  │    距离远 → bias 小 → 注意力弱             │
  └───────────────────────────────────────────┘
        │
        ▼
  增强后的 pocket 表示 [N_pocket, 512]
```

**为什么需要距离偏置？** 原始的 Transformer 自注意力是"全连接"的（每个位置都能关注所有其他位置）。对蛋白质口袋来说，空间上近的残基往往功能相关（协同催化、形成氢键网络），用距离作为先验能引导模型优先关注这些残基对。

---

### 4.5 几何增强交叉注意力（核心创新之二）

这是 EnzymeCAGE 最独特的地方。和 EZSpecificity 的双向交叉注意力不同，EnzymeCAGE 做**4 路交叉注意力**，而且**每一路都注入 3D 几何偏置**。

#### 4.5.1 4 路交叉注意力

```
4 组 Query-Key 配对：

  路径 1：enzyme → substrate
    Q = pocket 残基  [N_pocket, 512]
    K,V = 底物原子    [Ns, 128]
    → 酶的每个残基获得"与它相关的底物原子信息"

  路径 2：substrate → enzyme
    Q = 底物原子     [Ns, 128]
    K,V = pocket 残基 [N_pocket, 512]
    → 底物的每个原子获得"与它相关的酶残基信息"

  路径 3：enzyme → product   （use_prods_info=True 时启用）
    Q = pocket 残基
    K,V = 产物原子
    → 酶的每个残基获得"与它相关的产物原子信息"

  路径 4：product → enzyme   （同上）
    Q = 产物原子
    K,V = pocket 残基
    → 产物的每个原子获得"与它相关的酶残基信息"

每条路径的输出都做 mean pooling，得到 1 个 128 维向量。
4 路共 4 个 128 维向量 → 拼接成 512 维（interaction_fused）
```

EZSpecificity 只做 enzyme→substrate 和 substrate→enzyme 两路，因为没有产物。

#### 4.5.2 几何偏置（interaction_weight）：3D 先验注入注意力

这是 EnzymeCAGE 的"秘方"。普通交叉注意力只靠特征相似度决定关注谁，EnzymeCAGE 在注意力分数上**加一个 3D 几何偏置**：

```
interaction_weight[i, j] = pocket_weight[i] × substrate_weight[j]

其中：
  pocket_weight[i]：
    基于残基 i 到 pocket 几何中心的距离
    距离近 → 权重大（在口袋中心的残基更重要）
    距离远 → 权重小
    范围大约 [0, 0.2]

  substrate_weight[j]：
    基于原子 j 的反应中心标志
    是反应中心 → 权重大
    不是反应中心 → 权重小
    范围大约 [0.1, 0.6]

然后在交叉注意力中：
  attention_score[i, j] = Q_i · K_j / sqrt(d) + interaction_weight[i, j]
                          └─ 普通内积相似度 ─┘   └─ 3D 几何先验 ─┘
```

**直觉**：
- pocket 几何中心的残基 + 底物反应中心的原子 → **配对最重要**
- 外围残基 + 非反应原子 → 配对几乎不重要

这等于告诉模型："你不必自己瞎猜谁和谁相关，我告诉你空间上哪对最可能参与反应。"

#### 4.5.3 与 EZSpecificity 双向交叉注意力对比

| 维度 | EZSpecificity | EnzymeCAGE |
|---|---|---|
| 路径数 | 2（enzyme↔substrate）| 4（加上 enzyme↔product）|
| 几何先验 | 无 | **pocket 距离 + reacting center 权重** |
| Query/Key 维度 | 128 | enzyme 侧 512, substrate/product 侧 128 |
| 注意力头数 | 8 | 8 |
| 输出 | 2 个 [128] 向量 | 4 个 [128] 向量拼成 [512] |

---

### 4.6 末端大拼接

EnzymeCAGE 的末端拼接比 EZSpecificity 简单得多（没有那么多独立向量），但每个向量维度都很大：

```
interaction_fused [512]  ← 4 路交叉注意力输出
        +
ESM 全序列均值 [1152]    ← 全局酶身份信号
        +
DRFP 反应指纹 [2048]     ← 反应类型先验
        =
[3712 ~ 3840] 拼接向量（use_prods_info=False 时较小）
```

对比 EZSpecificity：
- EZSpecificity: 7 个 128 = 896 维
- EnzymeCAGE: 3 个大向量 ≈ 3840 维

**设计哲学不同**：
- EZSpecificity 倾向"很多个同维度的小向量"——统一到 128 维再拼接
- EnzymeCAGE 倾向"保留原始高维"——让预测头自己去压缩

---

### 4.7 预测头：更深的 MLP

```
  输入：[3840] 拼接向量
        │
        ▼
  ┌───────────────────────────────────────────────────┐
  │  层 1: BatchNorm(3840) + Dropout + Linear(3840→2048) + LeakyReLU │
  │  层 2: BatchNorm(2048) + Dropout + Linear(2048→1024) + LeakyReLU │
  │  层 3: BatchNorm(1024) + Dropout + Linear(1024→1)                │
  └───────────────────────────────────────────────────┘
        │
        ▼
  输出：logits [1]（原始分数）
  → sigmoid → 概率
```

**区别 EZSpecificity 的预测头**：
- EZSpecificity: 3 层，896→128→128→1，**无 BatchNorm**，使用 ReLU
- EnzymeCAGE: 3 层，3840→2048→1024→1，**有 BatchNorm**，使用 LeakyReLU

BatchNorm 对 EnzymeCAGE 很重要，因为输入维度差异大（ESM 1152 + DRFP 2048 + 交互 512），BatchNorm 能强制每个维度均值=0、方差=1，让训练更稳定。

---

## 五、与 EZSpecificity 的完整对比速查表

| 层面 | EZSpecificity | EnzymeCAGE |
|---|---|---|
| **任务** | 酶-底物特异性二分类 | 酶-反应匹配（检索）二分类 |
| **反应建模** | 只看底物 | **完整反应（底物 + 产物）** |
| **酶输入** | ESM-2 1280 全序列 | ESM-C 1152 或 ESM-2 1280 全序列 |
| **pocket 粒度** | 原子级（对接 10Å）| **残基级（预测结构 8Å）** |
| **结构编码器** | EGNN（SE(3) 不变）| **GVP（SE(3) 等变）+ SchNet** |
| **分子编码器** | GROVER (2D 预训练) + Morgan (2D 规则) | **SchNet (3D 几何)** |
| **反应中心标注** | 无 | **有（AAM 原子映射）** |
| **反应指纹** | Morgan (1024-bit) + GROVER mean (4885) | **DRFP (2048-bit)** |
| **产物信息** | 无 | **可选（use_prods_info）** |
| **pocket 内部注意力** | 无 | **有（距离偏置）** |
| **反应内部注意力** | 无 | **可选（reacting_center 加权）** |
| **酶-分子交叉注意力** | 2 路 | **4 路（含产物）** |
| **几何偏置注入注意力** | 无 | **有（pocket × reacting_center）** |
| **末端拼接维度** | 896（7×128）| 3840（混合大向量）|
| **预测头** | 3 层 MLP (ReLU) | 3 层 MLP (BatchNorm + LeakyReLU) |
| **预训练** | 无 | **有（foundation model）** |
| **发表年份** | Nature 2025 | Nature Catalysis 2026 |

**一句话总结**：EnzymeCAGE 是"EZSpecificity 的几何增强+反应级升级版"：
- 把 EZSpecificity 的**原子级结构**升级成**残基级几何 (GVP)**
- 把 EZSpecificity 的**单底物**升级成**完整反应 (substrate + product)**
- 新增了**反应中心标注**和**3D 几何偏置**作为强先验
- 用了**更大的预测头**来处理更丰富的输入

---

## 六、训练与推理

### 6.1 训练

EnzymeCAGE 是**基础模型**（foundation model），训练流程分两步：

**Step 1（预训练）**：在 RHEA 数据库（几十万酶-反应对）上大规模预训练
- 每个正样本：(UniProt ID, 反应 SMILES) 对
- 负样本：对每个反应，随机采样错误的酶作为负样本（1:N 比例）
- 优化目标：标准二分类交叉熵

**Step 2（下游评估）**：在小的测试集上评估零样本性能或做少量微调
- Orphan-335：对孤儿反应做酶检索（给反应排序候选酶）
- Enzyme-405：对新酶做功能预测
- 外部测试集（P450 / Terpene / Phosphatase）：家族级评估

EZSpecificity 没有预训练概念，它直接在 ESIBank 上从头训练。

### 6.2 推理

```
对每一对 (enzyme, reaction)：

  1. 查表拿到 enzyme 的 ESM 特征、GVP 特征
  2. 查表拿到 reaction 的 DRFP、反应中心、3D 构象
  3. 走第四章的完整前向计算
  4. 输出 logits → sigmoid → 概率
```

---

## 七、EnzymeCAGE 对我们 P450 研究的启示

回顾本项目 EXP005（双图架构 Dualgraph 2+）的实验：我们试图**在 EZSpecificity 上加一条残基级 GVP 通路**，结果 Test AUC **下降 0.0067**。

### 7.1 为什么 EXP005 没达到 EnzymeCAGE 的效果？

| 维度 | EnzymeCAGE | 我们 EXP005 |
|---|---|---|
| 残基级 GVP | ✅ 作为**主要**结构编码（甚至替代原子级）| 仅作**补充**通路（原子级 EGNN 仍是主干）|
| 反应中心 | ✅ 显式标注并作为 attention bias | ❌ 没有 |
| 产物信息 | ✅ 完整反应 | ❌ 只有底物 |
| 几何偏置注入注意力 | ✅ pocket × reacting_center | ❌ 普通交叉注意力 |
| 预训练 | ✅ 几十万对 | ❌ P450 从头训练 |
| 数据量 | RHEA 数十万对 | 47,510 对 |

我们只移植了 GVP 这一部分，没有移植**反应级建模**、**反应中心加权**、**几何偏置注入**等真正让 EnzymeCAGE 起效的核心机制。相当于把 BMW 的发动机装到自行车上——发动机本身没问题，但缺少与之配套的车架、变速箱、轮胎。

### 7.2 如果要进一步靠拢 EnzymeCAGE，需要做什么？

以下三步按工作量排序：

1. **数据侧增加"产物信息"**（工作量中等）
   - 需要为 P450 数据集每个正样本找到对应的反应产物
   - 挑战：P450 的许多反应记录没有明确产物（特别是 plant P450 的区域选择性反应，产物可能是混合物）
   - 收益预估：可能不大，因为 P450 数据同质性高，反应类型有限

2. **引入几何偏置注入交叉注意力**（工作量小）
   - 用现有 pocket_residue_idx 数据 + 某个"催化位点先验"（比如 Fe 血红素附近的残基）
   - 把距离偏置加到现有交叉注意力的 logits 上
   - 收益预估：可能小幅提升，值得一试

3. **完整迁移 EnzymeCAGE 架构**（工作量大）
   - 换掉 EZSpecificity 主干为 EnzymeCAGE
   - 需要重做整个数据预处理（AAM、反应中心提取、SchNet 构象生成）
   - 不再是"改进 EZSpecificity"，而是"用 EnzymeCAGE 在 P450 上复现+微调"
   - 收益预估：可能显著，但毕设时间不够

### 7.3 更现实的下一步：**零样本用 EnzymeCAGE 评估我们的 P450 测试集**

EnzymeCAGE 提供了预训练 checkpoint 和专门的 P450 外部测试脚本（`shells/evaluate_p450.sh`）。我们可以：

1. 下载 EnzymeCAGE 预训练 ckpt
2. 准备好我们 P450 数据集的反应 SMILES + 预测结构
3. 运行 EnzymeCAGE 零样本推理
4. 对比三方：**我们的 EXP001 (0.9320) vs EnzymeCAGE 零样本 vs 论文 EZSpecificity (0.5586)**

如果 EnzymeCAGE 零样本就能接近我们的 0.9320，说明"残基级几何 + 反应级建模 + 预训练"这套路线**真有价值**，毕设可以讨论这个对比作为未来方向；如果差距很大（EnzymeCAGE << 我们），则进一步验证了"P450 家族内部的特异性比 EnzymeCAGE 的通用催化建模更需要家族专属训练"。

无论哪种结果，都是毕设答辩的强材料。

---

## 附录：术语快速对照

| 术语 | EZSpecificity 中的含义 | EnzymeCAGE 中的含义 |
|---|---|---|
| pocket | 底物原子 10Å 内的蛋白质原子集合 | 预测结构里 8Å 的 pocket 残基集合 |
| 节点 | 原子 | 残基（pocket 图）或 原子（SchNet 图）|
| 边 | 化学键 + kNN 空间近邻 | kNN 残基邻居 / radius 原子邻居 |
| 消息传递 | EGNN（3 层）| GVP（3 层）+ SchNet（6 层）|
| 结构编码器 | SE(3) 不变（只用距离）| SE(3) 等变（距离 + 向量）|
| 分子编码器 | GROVER（2D 图）+ Morgan（指纹）| SchNet（3D 图）|
| 反应 | 只有底物 | 底物 + 产物（双分子图）|
| 反应中心 | 不建模 | AAM 原子映射标注 |
| 反应指纹 | Morgan 1024 位 + GROVER 4885 维 | DRFP 2048 位 |
| 交叉注意力 | 2 路（双向）| 4 路（含产物）+ 3D 几何偏置 |
| 预测头 | 896 → 128 → 128 → 1 | 3840 → 2048 → 1024 → 1 |

---

## 八、代码入口索引

读者若想深入 EnzymeCAGE 代码，从以下入口开始（路径均相对 `d:/EZSpecificity_Project/毕业设计/EnzymeCAGE/`）：

| 入口 | 文件 | 内容 |
|---|---|---|
| 主模型 | `enzymecage/model.py` | EnzymeCAGE 主类 + GVP_embedding + 反应中心权重计算 |
| 基类 | `enzymecage/base.py` | BaseModel（评估指标计算）|
| GVP 模块 | `enzymecage/gvp/__init__.py` | GVP / GVPConv / GVPConvLayer 实现 |
| SchNet | `enzymecage/schnet.py` | SchNet / InteractionBlock / GaussianSmearing |
| 交叉注意力 | `enzymecage/interaction.py` | CrossAttention / EnzymeCompoundCrossAttention（4 路）|
| 多头注意力 | `enzymecage/attention.py` | MultiHeadAttention + attn_bias 支持 |
| 数据集 | `enzymecage/dataset/geometric.py` | PyG Dataset（HeteroData 格式）|
| 特征提取 | `feature/main.py` | ESM / GVP / DRFP / 反应中心 / 构象 全流程 |
| 训练入口 | `train.py` | 完整训练 loop |
| 推理入口 | `infer.py` | 推理 + 结果保存 |
| 评估入口 | `evaluate.py` | AUC / 检索指标计算 |
| 论文原型 | `README.md` | 官方使用指南 |

本文档最后更新：2026-04-16
