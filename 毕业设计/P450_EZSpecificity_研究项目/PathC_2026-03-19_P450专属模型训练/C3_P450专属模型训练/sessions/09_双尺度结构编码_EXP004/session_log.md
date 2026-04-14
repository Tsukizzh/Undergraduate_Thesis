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

---

## 七、EXP001_fixed / EXP002a_fixed 基线重跑（2026-04-14）

### 7.1 决策

用户要求：只重跑 EXP001 和 EXP002a（放弃 EXP002b 的 lr tuning 变体），**全部使用 EXP003_fixed 的训练配方**。这样三个实验构成极干净的单变量 ablation：

| 实验 | feature_dim | 变化 | 训练配方 |
|---|---|---|---|
| EXP001_fixed | **28** | bare baseline (no Fe/HEM/residue_geom) | 同 EXP003_fixed |
| EXP002a_fixed | **31** | +Fe atom + HEM residue + is_hetero | 同 EXP003_fixed |
| EXP003_fixed（已完成）| **37** | +φ/ψ/χ1 残基几何 | — |

**EXP003_fixed 训练配方**（作为统一超参）：
- lr: 4.0e-04, warmup_epochs: 12, weight_decay: 1.0e-05, accumulate_grad_batches: 1
- sched_patience: 8, gradient_clip_val: 8, optimizer: adamW
- batch_size: 88, devices: 4, num_workers: 6, max_epochs: 200

### 7.2 pt_cache 对齐修复扩展到 base + heme

**背景**：之前只修了 pt_cache_geom（for EXP003_fixed），现在需要把 pt_cache 和 pt_cache_heme 也做 overlay。

**脚本**：`fix_cache_overlay.py`（单文件支持两种布局）
- per-sample 布局（pt_cache）：symlink samples/ 整个目录 + 过滤 index.pt
- shard-only 布局（pt_cache_heme，pt_cache_geom）：symlink graph_*.pt + 过滤 index.pt
- 每个 split 的 `enzymes/` 指向共享的 `pt_cache/shared_fixed/enzymes/`
- Codex 两轮审查通过

**产出**（全部 symlink overlay，只占 <2 MB 总空间）：
- `pt_cache_fixed/` (1.1M)
  - random: 47,411 → 44,284 (-6.60%, 3127 orphan samples)
  - all: 16,492 → 15,567 (-5.61%, 925 orphan samples)
- `pt_cache_heme_fixed/` (808K)
  - random: 47,807 → 44,680 (-6.54%, 3127 orphan samples)

### 7.3 端到端数据正确性验证

对 pt_cache_fixed 和 pt_cache_heme_fixed 做了 5 项验证：

| 项 | 结果 |
|---|---|
| 每个 index.pt 里 enzyme_ids ⊆ fixed flatbin keys | ✅ 12 个 index 全部 0 orphan |
| Sample 文件 enzyme_id 和 index.pt 完全一致 | ✅ 5 个随机抽样全部吻合 |
| Buggy vs Fixed flatbin 同 key 字节不同 | ✅ 10/10 全部字节不同（证明 fix 真的换了数据）|
| pt_cache 和 pt_cache_heme 过滤相同的 orphan 酶 | ✅ 32 个 orphan id 集完全一致 |
| 样本丢失率 | random ~6.6%, all ~5.6% |

**端到端 chain 验证**（深度 probe）：
- LMDB remap：`buggy[compressed] == fixed[csv_row]` 7 个抽样全部字节相等
- Flatbin ↔ fixed LMDB：5 个 key 字节完全相等（fp16 精度）
- Sample → flatbin：3 个随机 sample 的 enzyme embedding 与 LMDB 匹配，shape 和序列长度都对

### 7.4 PL 1.x → PL 2.x 迁移

**背景**：EXP001/002a 代码是 4090 + PL 1.x 时代写的，5090 服务器是 torch 2.8 + PL 2.6，API 破坏性变更。

**Codex 两轮讨论输出的补丁清单**（按代码家族）：

**家族 1：EXP001**（`src/Models/ss.py` + `cpi.py` + `scripts/main_training_pt.py`）
**家族 2：EXP002a**（同上，但 main_training_pt.py 已有 `_apply_runtime_patches` 辅助函数）

**具体补丁**：
1. `ss.py`/`cpi.py` 加 `self.validation_step_outputs = []` 和 `self.test_step_outputs = []`
2. `test_step`/`validation_step` append 到各自 list
3. `validation_epoch_end(outputs)` → `on_validation_epoch_end()` + 从 list 读
4. `test_epoch_end(outputs)` → `on_test_epoch_end()`
5. 新增 `on_validation_epoch_start` / `on_test_epoch_start` 清空 list
6. `main_training_pt.py` 里 `on_before_optimizer_step(..., opt_idx=0)` 去掉 opt_idx
7. `precision=16` → `"16-mixed"`
8. `strategy=None` → `"auto"`（单卡时）
9. 整个手写 DDP gather 块 → `dist.all_gather_object`
10. `num_sanity_val_steps=2` → `0`
11. `shutdown -h now` → `/usr/bin/shutdown`

**smoke test 漏掉的 bug**：`SS.optimizer_step` 签名是 PL 1.x 7 参数版本，PL 2.x 只传 4 个。第一次 EXP002b_fixed 启动跑了 1 个 batch 就崩。补丁：
```python
# OLD (PL 1.x)
def optimizer_step(self, epoch, batch_idx, optimizer, optimizer_idx,
               optimizer_closure, on_tpu, using_native_amp, using_lbfgs):

# NEW (PL 2.x)
def optimizer_step(self, epoch, batch_idx, optimizer, optimizer_closure):
```

**教训**：单纯的 `SS(cfg)` 实例化 smoke test 测不到 training loop 里被 PL 调用的 hooks。真正的 smoke test 应当 run 一个 mini fit（1 GPU, bs=16, max_epochs=1）让 PL 实际调用所有 hooks。

### 7.5 踩坑记录（为什么启动比想象的慢）

1. **EXP002b_fixed 第一次启动崩**：optimizer_step 签名 bug → 已修
2. **盲目"优化" --num-workers 16 + --preload**：
   - 用户指令：严格按 EXP003_fixed 参数跑
   - 我自作主张把 workers 从 6 改到 16（想用满 208 核），并开 --preload（想用满 754GB RAM）
   - 结果：第二次启动 `ERROR: Unexpected bus error encountered in worker. This might be caused by insufficient shared memory`
   - 4 GPU × 16 workers = 64 worker 进程同时访问 /dev/shm，爆了
   - 用户明确纠正："我求你就按照EXP003的参数配置来跑"
   - 修复：严格复制 EXP003_fixed 的 run_train.sh（workers=6，无 preload），只改 EXP/CACHE/run-name
3. **paramiko 后台启动问题**：`nohup bash script &` 通过 paramiko 起不来，exec_command 会等 channel 关闭。改用 `(bash script > log 2>&1 </dev/null &)` 子 shell fork 才成功
4. **shell 变量展开陷阱**：Python heredoc 里的 `${EXP}` 被外层 bash 展开成空串，导致 `--cache-dir ` 后面什么都没有。改用本地 Write 文件 + ssh_upload 上传
5. **本地 DNS/socket buffer 抽风**：连续 paramiko 调用耗尽 Windows ephemeral port，`getaddrinfo failed`。用户重启电脑后恢复

### 7.6 EXP002a_fixed 已启动

2026-04-14 01:51 启动，4×RTX5090：
- Batch 88 × devices 4 × accum 1 = effective **352**/step（和 EXP003_fixed 完全一致）
- Train 22,384 samples / Val 11,148 samples
- Trainable parameters: **1,847,044**（和 smoke test 验证值完全一致，证明是真 31 维模型）
- GPU 利用率稳定 93-99% × 4 卡，显存 27GB/32GB 每卡
- 预计 ~1-1.5 小时完成（参考 EXP003_fixed 约 78 分钟）
- 跑完自动 `/usr/bin/shutdown`

**EXP001_fixed 等 EXP002a_fixed 跑完后由用户重开 4 卡再启动**。

### 7.7 目录结构

服务器端：
```
PathC/P450/experiments/
├── EXP001_fixed/          ← 28 dim baseline, 训练配方同 EXP003_fixed
│   ├── configs/config.yml (lr=4e-4, warmup=12, wd=1e-5, accum=1)
│   ├── scripts/run_train.sh (bs=88, devices=4, workers=6)
│   ├── src/ (full tree, ss.py/cpi.py/main_training_pt.py 全部 PL 2.x 补丁过)
│   └── (logs/ 和 results/ 训练时自动生成)
├── EXP002a_fixed/         ← 31 dim Fe/HEM, 训练配方同 EXP003_fixed
│   └── (同上，cache 指向 pt_cache_heme_fixed)
└── EXP003_fixed/          ← 37 dim residue geom, Test AUC=0.8943（已完成）

PathC/P450/data/
├── pt_cache/              ← 原版（bug, 不动）
├── pt_cache_fixed/        ← 1.1M overlay（本次新增）
├── pt_cache_heme/         ← 原版（bug, 不动）
├── pt_cache_heme_fixed/   ← 808K overlay（本次新增）
├── pt_cache_geom/         ← 原版（bug, 不动）
└── pt_cache_geom_fixed/   ← 之前修过的 overlay
```

### 7.8 下一步

1. 等 EXP002a_fixed 训练完成（自动关机）
2. 用户重开 4 卡 → 启动 EXP001_fixed
3. 两个都完成后，下载结果到本地 `results/10_EXP001_fixed/` 和 `results/11_EXP002a_fixed/`
4. 更新所有文档（session log / MEMORY / README / 项目进度日志）
5. git commit + push
6. 回到 EXP004 Step 14 规划（双尺度结构编码）

---

## 八、GROVER 对齐 Bug 发现（2026-04-14 深夜）

**⚠️ 重大中断**：EXP002a_fixed 训练过程中，用户追问"EXP003 的角度是否是对的酶"，深挖全链路审计意外发现 **GROVER LMDB 存在与 ESM 同类的对齐 bug**。

### 8.1 触发

用户关注点原本只是角度特征。我第一轮回答"只有 ESM 错位，结构通道（含角度）正确"——但用户要求实证验证。Codex 第一轮读代码时已基本定性，我随后用原子数做硬实证：**100% 的 LMDB key 原子数匹配 `grover_substrates.csv[N]`，精确断点落在 k=8**（`*[H]` 被删位置），Codex 二次复核判为 "effectively airtight"。

### 8.2 Bug 本质

`phase7_step5_grover.sh` 预处理时 `.dropna(subset=['Substrate_SMILES'])` 本应保留全部 2125 行，但实际 `grover_substrates.csv` 只有 2124 行——当时因 GROVER 在 `*[H]`（Substrate Index 8）上崩溃（MEMORY 里的 "Fixed `*[H]` index-out-of-bounds"），处理方式是**直接删行未补位**。GROVER 的 LMDB 用顺序计数器作 key，导致：

```
key 0..7       → Substrate Index 0..7      ✓
key 8..2123    → Substrate Index 9..2124   ✗（错位 1 格）
(key 2124 缺失)
```

**99.6% 的底物加载了错误的 GROVER 嵌入。**

### 8.3 受影响实验

| 实验 | ESM | GROVER |
|---|---|---|
| EXP001/002a/002b/003 | ❌ Bug 1 | ❌ Bug 2 |
| EXP003_fixed (Test=0.8943) | ✅ | ❌ **未修** |
| EXP002a_fixed (训练中) | ✅ | ❌ **未修** |

**EXP003_fixed 的绝对数值仍带 GROVER 污染**。但 EXP003→EXP003_fixed 的 +0.1029 增量确实只是 ESM 修复贡献。

### 8.4 完整记录

详见独立文档：[GROVER对齐bug发现_2026-04-14.md](GROVER对齐bug发现_2026-04-14.md)

该文档包含：
- 完整数据链路（三份 CSV 词汇表 + 5 条特征管线 + pt_cache 装配）
- 两次 Bug 的共同根源（压缩计数器 key 模式）
- GROVER Bug 的硬证据（CSV 对齐 + 原子数 100% 匹配实证）
- 修复方案 A/B（推荐 B：纯 rekey）
- pt_cache 连锁重建计划
- 教训与规则（禁用压缩计数器 key、Morgan 填零模式推广）

### 8.5 待决策

1. 是否立即停 EXP002a_fixed
2. 方案 A 重跑 vs 方案 B rekey
3. 修复后重跑顺序：EXP001_fixed → EXP002a_fixed → EXP003_fixed

---

## 九、AllFix 系列修复与执行（2026-04-15）

### 9.1 决策

- 方案 B（纯文件 rekey）执行
- 中断 EXP002a_fixed，废弃 fixed 系列，改走 allfix 命名
- 同时构建 natural（各自 orphan 过滤）和 unified（三套 sample_id 交集）两套 cache，natural 跑最大数据量，unified 跑严格 feature_dim ablation

### 9.2 五阶段修复（每步多轮 codex + 字节级实证）

1. **Phase 1 rekey LMDB** (`scratch/fix_grover_lmdb.py`)：规则 `new_int = old_int if old_int < 8 else old_int + 1`。Codex 关键点：`txn.put(overwrite=False)` 返回 False 不抛异常，必须 `assert ok`。字节级全扫 0/2116 + 0/8 mismatch。
2. **Phase 2 rebuild flatbin** (`scratch/build_allfix_substrates.py`)：复用 `build_substrate_shards + convert_substrate_shards_to_flatbin`。Config 陷阱（codex 纠正）：`grover_path` / `morgan_path` 必须是 list 不是 str。全扫 2124×3 字节 0 mismatch。
3. **Phase 3 indices** (`scratch/build_allfix_indices.py`)：生成 6 套 index.pt。Natural 保留各自 orphan 过滤（bare 22178 / heme 22384 / geom 22312 train），0 samples dropped（无样本用 Substrate Index 8）。Unified 取三套 sample_id 交集：train 22083 / val 11008 / test 11000。事先用 `check_sample_id_consistency.py` 验证三套 cache 的 sample_id 含义完全一致（0 disagreements）。
4. **Phase 4 symlink dirs** (`scratch/build_allfix_dirs.py`)：建 6 个目录。Symlink 策略依 layout：bare（per-sample）symlink `samples/`；heme/geom（shard-only）symlink 所有 `graph_*.pt`。所有目录含 manifest 副本 + `enzymes/` 共享 symlink + `substrates/` 指向新 allfix flatbin。
5. **Phase 5 端到端验证**：`verify_phase5_v3.py` 通过真实 `PtCacheDataset` 加载，确认 `protein_x.shape[-1]` 匹配预期 feature_dim（28/31/37）；`verify_final_sub9.py` 对 sub_id in [100, 1000] 验证 `sample == fixed_lmdb[sub_id] == True` 且 `sample == old_lmdb[sub_id] == False`。6 caches × 2 sub_ids 全通过。

### 9.3 数据传输与实验骨架

- 260MB 新文件通过 rsync 从 西北→北京（反向 pull，因 西北→北京 SSH 被阻）
- `scratch/setup_allfix_experiments.py` 从 `_fixed` 实验 `cp -a` 派生 6 个 `_allfix{,_unified}` 目录，regex patch `run_train.sh` 路径与 `config.yml` data.tag

### 9.4 结果（4×RTX4090 DDP, bs=88 eff=352, lr=4e-4, warmup=12, wd=1e-5, dropout=0.9）

| 实验 | feature_dim | Best ep | **Test AUC** | Test AUPR | 备注 |
|---|---|---|---|---|---|
| EXP001_allfix_unified | 28 (bare) | ep43 | **0.9320** | **0.6749** | ✅ |
| EXP002a_allfix_unified | 31 (+Fe/HEM) | ep59 | **0.9270** | 0.6300 | ✅ **-0.005** |
| EXP003_allfix_unified | 37 (+φ/ψ/χ1) | 🔄 ep74+ | — | — | best Val=0.9183@ep62 |

### 9.5 震撼发现

1. **bare baseline 从 ~0.77 跳到 0.9320**，GROVER bug 单独贡献 +0.04（相对 EXP003_fixed=0.8943）。两个 LMDB bug 合计影响 +0.16。
2. **Fe/HEM 在干净数据上反而掉点**。此前 EXP002a > EXP001 的优势完全是 GROVER 错位嵌入对 Fe 特征的偶然补偿；bug 修复后真实效果反转。
3. **残基几何 37 维当前 Val AUC 低于 28/31 维**，EXP003 的 "+0.0025 增量" 也是 bug 污染产物，Step 13/14 的双尺度结构编码方向严重存疑。
4. **EXP001-003 + EXP003_fixed 的整条 ablation 链全部作废**，feature_dim 单变量消融结论在干净数据上不成立。

### 9.6 后续

- 等 EXP003_allfix_unified 完成
- 3 套 natural 变体（最大数据量）训练
- 视最终结论决定 Step 14 是否继续，还是转向其他创新方向（dropout 扫、增广、Stage 2 多标签）
