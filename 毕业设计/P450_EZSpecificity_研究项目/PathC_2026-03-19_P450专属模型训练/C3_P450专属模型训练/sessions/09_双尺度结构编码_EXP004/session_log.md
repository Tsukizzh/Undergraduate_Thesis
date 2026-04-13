# Session 09：EXP004 双尺度结构编码 — 规划与决策

**实验编号**：EXP004
**对应研究路径**：Path C3 Step 14
**起始日期**：2026-04-13
**基准实验**：EXP003（残基几何特征 φ/ψ/χ1，Test AUC=0.7914）
**目标**：在 EXP003 基础上新增一条残基级 GVP-GNN 通路，形成双尺度结构编码（原子级 EGNN + 残基级 GVP）

---

## 一、实验背景与动机

### 1.1 为什么要做 EXP004

从 Codex + Gemini 三方讨论确认的事实出发：

1. **现有"EGNN"其实不是真 EGNN**：`src/Models/Structure/egnn.py` 里的 `EnBaseLayer.forward` 签名是 `(h, edge_index, edge_attr)`，**没有坐标参数**，也不做坐标更新。它是一个基于预计算 Gaussian RBF 距离的标量消息传递 GNN，是 **SE(3)-不变**（invariant）而非 SE(3)-等变（equivariant）。
2. **现有结构通道从头到尾只用标量**：
   - 节点特征：原子类型 one-hot、氨基酸类型 one-hot、骨架/侧链标志、角度 sin/cos
   - 边特征：`bond_feat [6] + dist_feat [24 Gaussian RBF] + cross_bond [1]`
   - 在 `EdgeConnection.__call__` 结束时 `data.protein_pos = None` 被显式丢弃
   - **任何 3D 方向向量都没有进入网络**
3. **GVP 带来的是 100% 净新增信息**：
   - `node_v [N_res, 3, 3]`（主链前向、主链反向、虚拟侧链方向）
   - `edge_v [E, 1, 3]`（CA-CA 单位方向向量）
   - 这些向量通道是现有架构完全没有的

### 1.2 与 EXP003 的关系

| 维度 | EXP003 | EXP004 |
|---|---|---|
| 结构通路数 | 1（仅原子级 EGNN）| 2（原子级 EGNN + 残基级 GVP）|
| 几何信息 | 标量距离 + 标量角度 | 标量 + **显式方向向量** |
| 粒度 | 原子级 | 原子级 + 残基级 |
| feature_dim | 37（EXP003 修改版）| 37 + GVP 残基嵌入（128 维独立通路）|

---

## 二、决策记录

### ✅ 决策 1：融合方式 — 选项 ②+（残基回填 + g_res bypass 双出口设计）

**最终方案**：GVP 输出同时走两条路，对称 EZSpecificity 结构通道的设计哲学。

#### 2.1.1 具体架构

```
GVP 残基通路（新增）
        │
        ▼
  GVP-GNN 处理（残基级 k-NN 图，k=30）
        │
        ▼
  h_res [N_pocket, 128]  ← 每个口袋残基一个 128 维嵌入
        │
        ├──── 出口 1：详细版 ────┐
        │                       │
        │                       ▼
        │              Linear 投影
        │                       │
        │                       ▼
        │              x_pro_enhanced[pocket_idx] +=
        │                   alpha * LayerNorm(Linear(h_res))
        │              （alpha 可学习，小值初始化）
        │                       │
        │                       ▼
        │              进入现有双向交叉注意力
        │
        └──── 出口 2：总结版 ────┐
                                │
                                ▼
                        scatter_mean(h_res)
                                │
                                ▼
                        g_res [128]
                                │
                                ▼
                        不进注意力，直接 bypass
                                │
                                ▼
                    末端 concat 第 8 个向量
                                │
                                ▼
                    预测头输入：896 → 1024
```

#### 2.1.2 为什么选这个

1. **对称 EZSpecificity 原作者的设计**：
   - 原架构中，结构通道有两个出口：
     - `x_str [280, 128]` → 进交叉注意力
     - `str_mean [128]` → 直接 bypass 到末端 concat
   - 这个模式被 Codex 和 Gemini 共同标记为**原则性设计**（Principled，A 类）
   - 我们为新增的 GVP 通路采用同样的双出口模式，和原作者风格一致
2. **防信息丢失**：
   - 万一 h_res 回填到 x_pro 后被 ESM 的 1450 维序列稀释
   - 经过交叉注意力 + 平均池化后 GVP 信号可能严重衰减
   - g_res bypass 保证预测头总能直接拿到 GVP 的全局总结
3. **自带消融实验**：
   - 可以关闭出口 1 只保留出口 2 → 看 g_res 单独贡献
   - 可以关闭出口 2 只保留出口 1 → 看 h_res 单独贡献
   - 双开 → 看联合效果
4. **成本极小**：相比原版选项 ②，只多了 ~15 行代码和 ~16K 参数（预测头第一层）

#### 2.1.3 Codex 标记的陷阱（实施时必须规避）

1. **残基索引对齐是最大风险点**：
   - `pocket_residue_idx` 必须和 ESM 生成的 `x_pro` 里的残基位置严格对齐
   - 当前代码路径里这个字段**不存在**，需要在 `src/Datasets/Structure/protein_ligand.py` 或 pt_cache 构建阶段新增
   - PDB 残基编号、口袋提取顺序、ESM 序列位置——三者必须完全一致
2. **用 `index_add_` 而非 `x_pro[idx] +=`**：
   - 口袋残基可能因为被切成两段或其他原因出现重复
   - `index_add_` 有正确的 scatter 语义
3. **必须加可学习门控**：
   - `x_pro_enhanced[idx] += alpha * LayerNorm(Linear(h_res))`
   - `alpha` 初始化为很小的值（例如 0.01）
   - 防止 GVP 信号一开始就压过 ESM，给模型自适应空间
4. **padding 位置不能被污染**：
   - `x_pro` 被 padding 到 1450，真实序列长度可能只有 ~400
   - `pocket_idx` 必须在真实序列长度内，不能落到 padding 区域
5. **预期的增益主要来自 `weighted_sum_reaction`**：
   - 因为 `x_pro` 的平均池化会把注入信号稀释在 1450 位置上
   - 真正的增益通过底物方向的交叉注意力体现

---

### ✅ 决策 2：残基图节点范围 — 选项 b（仅口袋残基）

**最终方案**：GVP 残基图只包含配体 10Å 范围内的残基（~30-60 个，具体数量实施时实测）。

#### 2.2.1 为什么选这个

1. **和现有 EGNN 原子通路同范围**：
   - 现有 EGNN 看配体 10Å 范围内的原子
   - GVP 也看同一个区域，只是粒度不同（残基 vs 原子）
   - 实验变量极干净：差别只有"新增残基级通路"这一项
2. **数据已现成**：
   - 使用现有 `P450/data/structure/complex/*.pdb` 文件
   - 零额外数据工程，ddl 前时间最紧
3. **计算成本很小**：
   - GVP 在 30-60 残基图上 k-NN (k=30) 只产生 ~1200 条边
   - 远小于原子级 EGNN 的规模

#### 2.2.2 实施时必须验证

- 现有 `complex/*.pdb` 到底是"口袋+配体"格式还是"全蛋白+配体"格式
- 如果是后者，在 GVP 特征提取时仍然要做 10Å 残基过滤

#### 2.2.3 拒绝选项 a（全蛋白）的理由

- 需要重新准备 1577 个酶的完整结构文件，1-2 天数据工程
- 实验变量不干净（同时改粒度+改范围）

#### 2.2.4 拒绝选项 c（口袋+扩展环）的理由

- 既不是和 EGNN 干净对比，也不是全蛋白视角
- 要单独写"扩展残基"逻辑，工程成本不匹配收益

---

### ✅ 决策 3：GVP 代码来源 — 选项 A（原版移植 + 逐项代码核对）

**最终方案**：从 EnzymeCAGE 仓库直接移植 GVP 实现，**但必须逐项核对**而非盲拷。

#### 2.3.1 来源文件

从 `毕业设计/EnzymeCAGE/enzymecage/gvp/` 移植以下文件：
- `__init__.py` — GVP、GVPConv、GVPConvLayer、GVP_embedding 核心实现
- `data.py` — ProteinGraphDataset（node_s / node_v / edge_s / edge_v 构建）
- 可能需要的：`gvp_dataset.py`、`prepare_feature.py`（取决于移植范围）

#### 2.3.2 核对清单（实施时 Codex 要逐项验证）

1. **依赖兼容性**：
   - 是否依赖 `torchdrug`、`atom3d`、`torch_cluster` 等 EZSpecificity 现有环境没有的包？
   - 和我们 PyTorch 2.1.0 + CUDA 12.1 的版本兼容性
   - 和 `torch_geometric` 现有版本的兼容性
2. **API 一致性**：
   - `GVP_embedding.forward` 的输入输出 shape 是否确实是文档里说的 `(6, 3) → (128,)`
   - 是否有隐藏的必需参数（如 batch、edge_index 格式）
3. **数据格式对齐**：
   - EnzymeCAGE 的 `ProteinGraphDataset` 输入是什么格式？字典、Data 对象、还是 PDB 路径？
   - 我们需要在哪一层做格式转换，才能让现有 pt_cache 流程喂得进去
4. **残基编号体系**：
   - EnzymeCAGE 用 PDB 残基编号还是连续序号？
   - 这决定了**决策 1 里的 pocket_residue_idx 构建方式**
5. **是否有训练好的权重需要下载**：
   - 如果 GVP 模块有预训练权重，路径在哪，下载方式是什么
6. **许可证**：
   - EnzymeCAGE 是什么 License？能不能在我们项目里复用

#### 2.3.3 为什么选 A 而非 B（最小重写）

1. **时间最值钱**：ddl 就在眼前，省 2 天等于多一次失败重跑机会
2. **正确性有论文背书**：EnzymeCAGE 发表在 Nature Chem Eng 2025，代码经过同行评议
3. **Codex 审查成本低**：只要看"移植是否忠实 + 输入输出对齐"即可，不需要验证数学正确性

---

### ✅ 决策 4：数据存储与管线接入 — 选项 γ（新 LMDB + 新 pt_cache 文件夹）

**最终方案**：顺着现有两级数据管线走。

#### 2.4.1 数据管线设计

```
原有 LMDB 原料（全部不动）          🟥 新增 LMDB 原料
├── enzyme_features.lmdb           └── structure_features_gvp.lmdb
├── reaction_features.lmdb                 （1577 个酶的 GVP 特征，
├── grover_fingerprint.lmdb                  按 enzyme_id 为 key）
└── structure_features.lmdb
           │                              │
           └──────────┬───────────────────┘
                      ▼
           🟥 修改版 build_pt_cache.py
           （读所有 LMDB，按样本组装）
                      │
                      ▼
           🟥 新建：pt_cache_geom_gvp/
                      ├── sample_00001.pt  (原数据 + GVP)
                      ├── sample_00002.pt  (原数据 + GVP)
                      └── ...（47,807 个）

           原有：pt_cache_geom/（完全不动，EXP003 仍可跑）
```

#### 2.4.2 新 LMDB 的内容结构

```python
# structure_features_gvp.lmdb 每个条目
{
    'enzyme_id': str,
    'node_s': np.ndarray [N_pocket, 6],          # 标量节点特征
    'node_v': np.ndarray [N_pocket, 3, 3],        # 向量节点特征
    'edge_index': np.ndarray [2, E],              # k-NN 连接
    'edge_s': np.ndarray [E, 32],                 # 边标量特征
    'edge_v': np.ndarray [E, 1, 3],               # 边向量特征
    'pocket_residue_idx': np.ndarray [N_pocket],  # 口袋残基在酶全长序列的位置
}
```

#### 2.4.3 修改版 build_pt_cache.py（伪代码）

```python
for sample in all_47807_samples:
    # 原有逻辑（保持不变）
    enzyme_data = enzyme_lmdb[sample.enzyme_id]
    substrate_data = reaction_lmdb[sample.substrate_id]
    grover_data = grover_lmdb[sample.substrate_id]
    structure_data = structure_lmdb[sample.dock_id]

    # 🟥 新增：读 GVP LMDB
    gvp_data = gvp_lmdb[sample.enzyme_id]

    combined = {
        **enzyme_data,
        **substrate_data,
        **grover_data,
        **structure_data,
        **gvp_data,  # 🟥 新增字段：gvp_node_s, gvp_node_v, gvp_edge_*, pocket_residue_idx
    }

    # 🟥 输出到新文件夹
    torch.save(combined, f"pt_cache_geom_gvp/sample_{sample.id}.pt")
```

#### 2.4.4 训练时的接入

```bash
# EXP003 训练（不变）
python main_training_pt.py --cache-dir pt_cache_geom/ ...

# EXP004 训练（只改 cache 路径）
python main_training_pt.py --cache-dir pt_cache_geom_gvp/ ...
```

Dataset 代码几乎不用改，只在 `__getitem__` 里多读几个 GVP 字段即可。

#### 2.4.5 为什么选 γ 而非 α/β

- **α（改原 .pt）**：破坏 EXP003 可复现性 → 拒绝
- **β（独立 LMDB 旁路，Dataset 两处读）**：
  - 省磁盘（30MB vs 16GB），但 **AutoDL 磁盘充足**
  - 代码复杂度高（两处读取 + 拼接逻辑 + 新 Dataset 类）
  - bug 风险更大
  - 拒绝
- **γ（新建文件夹）**：
  - **代码改动最小**：Dataset 和 EXP003 几乎一样
  - **顺着现有两级管线走**：不是新发明
  - **失败回退零成本**：删文件夹
  - **保留 EXP003 可复现性**
  - **选中** ✅

#### 2.4.6 磁盘与时间估算

| 项目 | 估算 |
|---|---|
| 新 LMDB 大小 | ~30MB（1577 个酶 × ~20KB） |
| 新 pt_cache 大小 | ~16GB（和 pt_cache_geom 一样） |
| 构建新 LMDB 耗时 | 1577 次 GVP 特征提取，~10-20 分钟 |
| 重跑 build_pt_cache.py 耗时 | 30-60 分钟（47,807 次读写） |
| 磁盘总增加 | ~16GB（AutoDL 充足，无压力） |

---

## 三、待决策事项

### ✅ 决策 5：训练超参数 — 选项 i（完全沿用 EXP003）

**最终方案**：EXP004 所有超参数完全照搬 EXP003，零改动。

#### 2.5.1 沿用的参数清单

```yaml
training:
  lr: 4.0e-04
  batch_size: 88
  warmup_epochs: 12
  weight_decay: 1.0e-05
  min_lr: 5.0e-06
  max_epochs: 200
  accumulate_grad_batches: 1
  gradient_clip_val: 8
  optimizer: adamW
  sched_factor: 0.5
  sched_patience: 8

model:
  hidden_dim: 128
  cross_attention:
    n_head: 8
    dropout: 0.9
  graph:
    num_layers: 3
    attention: true
```

#### 2.5.2 为什么选这个

1. **单变量原则最干净**：和 EXP003 的差别只有"加了 GVP 通路"一项，论文里能清晰归因
2. **新增参数量小**：GVP 通路只增加约 +10% 参数，对最优 lr 的影响非常小
3. **时间成本最低**：一次跑完，不用搜超参数空间
4. **可迭代**：如果跑出来效果不好，再考虑调参是第二轮的事

#### 2.5.3 拒绝其他选项

- **选项 ii（降 lr 到 2e-4 或 3e-4）**：不再是严格单变量，论文归因会打折扣
- **选项 iii（搜 3-4 个 lr）**：训练时间 × 4，ddl 前来不及

### ✅ 决策 6：实验服务器 — 选项 P（AutoDL 4×RTX5090）

**最终方案**：在 AutoDL 开一台 4×RTX5090 实例跑 EXP004，和 EXP003 硬件完全一致。

#### 2.6.1 为什么选这个

- **硬件一致性**：EXP003 就是在 AutoDL 4×5090 上跑的（73 ep，~78 min）
- **训练曲线可直接对比**：同硬件同环境，排除"机器差异"变量
- **速度最快**：5090 > 4090，迭代快
- **选项 Q（Cloud-2 4×4090）**：也能跑，但比 5090 慢
- **选项 R（wanglab 导师机器，10.45.246.249）**：**不可行**，只有 16GB RAM，加载 pt_cache_geom_gvp ~16GB 会 OOM

#### 2.6.2 训练配置估算

| 参数 | 值 |
|---|---|
| GPU | 4 × RTX5090 |
| batch_size | 88 / GPU（和 EXP003 一致）|
| effective batch | 88 × 4 = 352 |
| num_workers | 6 / rank |
| preload | 是（--preload 把 pt_cache 全部读入 RAM）|
| 每 epoch 耗时 | ~85-100 秒（比 EXP003 的 64 秒略慢，因为多一条 GVP 通路）|
| 总 epoch | ~100-150（参考 EXP003 的 73 早停）|
| 总训练时长 | ~2-3 小时 |
| auto shutdown | 是（训练+测试完成后关机）|

---

## 🎯 所有决策完成

所有 6 个决策已全部敲定（2026-04-13）。下一步进入阶段 1：Codex 代码核对与依赖验证。

## 四、参考资料

### 文件
- **决策 1 详细图解**：`../../EXP004_决策1_融合方式图解.md`
- **EnzymeCAGE 架构分析**：`毕业设计/EnzymeCAGE_架构分析.md`
- **EnzymeCAGE 源码**：`毕业设计/EnzymeCAGE/`
- **现有架构分析**：`毕业设计/P450_EZSpecificity_研究项目/EZSpecificity模型流程图解_简化版.md`

### 讨论会话
- **Codex Session ID**：`019d8212-57f8-7cc3-91fa-a60c37cefa44`
- **Gemini Session ID**：`5ada727b-efc4-4dab-9509-fea51f51d626`

### 关键代码位置
- **当前主模型**：`src/Models/ss.py:SS`
- **当前结构通路**：`src/Models/Structure/structure.py:Graph`
- **当前"EGNN"**：`src/Models/Structure/egnn.py:EnBaseLayer`（实为距离不变 GNN）
- **当前特征构建**：`src/Datasets/Structure/transforms.py:FeaturizeProteinAtom, EdgeConnection`

---

## 五、下一步

- 等待决策 4/5/6 敲定
- 阶段 1（Codex 代码核对 + 依赖验证）启动前，必须完成所有决策
- 阶段 1 输出：可移植的 GVP 模块原型 + 依赖清单 + 集成 diff patch（但不执行）

---

## 六、⚠️ 阶段 1 期间的重大发现：LMDB 对齐 Bug（2026-04-13）

### 6.1 发现过程

进入 EXP004 阶段 1（Codex 代码核对）后，在审查 `pt_cache_geom` 构建流程时发现 `phase7_step2_esm.py` 中存在严重的 enzyme-sample 对齐 bug：

```python
# PathC/P450/archive/phase7_scripts/phase7_step2_esm.py
uniprot_dict[str(idx)] = (len(uniprot_dict), 1)  # ← 压缩键，非原 CSV 行号
```

真实 CSV 行号 `idx` 被替换为 `len(uniprot_dict)`（顺序压缩计数），导致 `enzymes.lmdb` 的 key 是"第 N 个通过过滤的酶"，而非 `Enzymes.csv` 里的实际行号。样本在 `build_pt_cache.py` 里用 CSV 行号查 LMDB → **绝大多数样本拿到的是错配的酶特征**。

### 6.2 影响面

- **错误覆盖范围**：EXP001 / EXP002a / EXP002b / EXP003 全部受影响
- **受影响的样本比例**：需要 skip 的酶（非标准氨基酸或 > 1000aa）之前的样本无影响，之后的样本全部错位
- **实际错配率**：样本级 ~95.8%（几乎全错）
- **为什么没被发现**：训练仍能收敛，AUC 也不难看（0.77-0.79），模型部分靠了底物/结构通道以及"任意两个 P450 酶都有相关性"的家族同质性

### 6.3 修复方案

**原则**：非破坏性，全部新增文件不覆盖。

| 阶段 | 操作 | 输出 |
|---|---|---|
| A1-A7 | `fix_enzyme_lmdb.py`：读原 LMDB，重放 Phase 7 过滤逻辑，用原 CSV 行号作为新 key | `structure/enzymes_fixed.lmdb` |
| B | `fix_flatbin_build.py`：按 `build_pt_cache.py` 654-663 行 schema 重建 flatbin | `pt_cache/shared_fixed/enzymes/enzymes.bin + enzymes_index.pt` |
| C | `fix_geom_cache.py`：创建 `pt_cache_geom_fixed/`，样本和 graph shards 软链接到原 `pt_cache_geom/`，新建过滤后的 index.pt | `pt_cache_geom_fixed/` |

**孤儿样本处理**：48 个被过滤掉的酶（45 个未过滤成功的 + 3 个 PDB 不匹配）对应的样本在新 index 里被剔除：
- train: 23916 → 22312（-1604）
- val: 11968 → 11111（-857）
- test: 11923 → 11107（-816）
- 总损失：6.9%

**验证**：
- 原 LMDB → 新 LMDB：1577 个条目全部字节相等（按新 key 检索）
- 新 flatbin 10 个采样与新 LMDB 字节相等
- 端到端随机抽样测试：ESM 嵌入和 CSV 序列匹配率从 4.2% → 97.6%

### 6.4 EXP003_fixed：修复后的 EXP003 重跑

**配置**：完全沿用 EXP003 所有超参数，只改 cache 路径为 `pt_cache_geom_fixed/`。

| 维度 | EXP003 (原) | EXP003_fixed | Δ |
|---|---|---|---|
| Best ckpt | ep58 Val AUC=0.7850 | **ep57 Val AUC=0.8891** | +0.1041 |
| **Test AUC-ROC** | 0.7914 | **0.8943** | **+0.1029** |
| Test AUPR | 0.3814 | 0.5358 | +0.1544 |
| Samples | 11108 | 11108 (pos=993, neg=10115) | — |
| 硬件 | 4×RTX5090 | 4×RTX5090 | — |

**含义**：
1. EXP001 → EXP002a → EXP002b → EXP003 的整条提升曲线全部建立在错配的酶-样本对上，真实的 P450 专属数据集威力被对齐 bug 压住了
2. 在正确对齐的 baseline 上，一次训练就从 0.7914 跳到 0.8943（+10.3 个点），远超之前全部优化叠加的幅度
3. EXP004 双尺度结构编码仍然要做，但**所有 baseline 都必须在 fixed cache 上重建**，否则对比不可信

### 6.5 对 EXP004 计划的影响

| 决策 | 原计划 | 修订 |
|---|---|---|
| 决策 4 cache 路径 | `pt_cache_geom_gvp/` 基于 `pt_cache_geom/` | **改为基于 `pt_cache_geom_fixed/` 构建** |
| 参照 baseline | EXP003 Test=0.7914 | **EXP003_fixed Test=0.8943** |
| 先决条件 | Codex 阶段 1 完成 | Codex 阶段 1 + 确认 fixed cache 已就绪 |
| EXP001/002a/002b 重跑 | — | **待定**：要不要在 fixed cache 上重跑 baseline 曲线用于论文 |

### 6.6 清理

服务器端：`/root/*.py` 的 11 个临时修复脚本（`fix_*.py`、`inspect_*.py`、`dry_run.py`）已删除。
本地：同名 9 个 `fix_*.py` 已从 `d:/EZSpecificity_Project/` 根目录删除。
EXP003 / EXP003_fixed 完整实验目录（configs/logs/results/scripts/src）已拉回 `results/08_残基几何特征_EXP003/` 和 `results/09_对齐bug修复_EXP003_fixed/`。
