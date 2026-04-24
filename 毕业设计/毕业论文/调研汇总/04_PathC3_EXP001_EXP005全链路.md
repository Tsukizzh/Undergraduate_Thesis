# Path C3 P450 专属模型训练全系列实验深度调研（核心章节）

调研日期：2026-04-23
覆盖范围：EXP001 → EXP002a/b/c → EXP003 → EXP003_fixed → EXP004 → AllFix 系列 → EXP005

---

## 1. 底物分类体系（C3 v6 FINAL）

**目标与规模**：将 2,125 个化合物分类为 7 类多标签（Terpenoid / Amino_acid / Fatty_acid / Alkaloid / Steroid / Phenylpropanoid / Polyketide），最终 **1,870 confirmed (88.0%) + 255 other (12.0%) + 63 multi-label**。

**五阶段 Pipeline**：
- **P1 Gold**：NPClassifier 训练集（78,336 化合物）精确 SMILES 匹配，Superclass（70 类）高精度映射
- **P2 NPC Superclass**：主力分类层（76.5% 标签来源），映射 70 个 NPC superclass 至 7 类
- **P3 NPC Pathway**：无 superclass 时降级方案
- **P4 Amino_acid SMARTS**：检测完整 α-氨基酸骨架（NH2-CH-COOH）
- **P5 Other**：无任何分类证据的默认类

**核心创新——多标签方案**：NPClassifier 论文采用 sigmoid（非 softmax），支持化合物同时属多个生物合成途径（如修饰色氨酸=生物碱+氨基酸，杂萜=萜类+聚酮）。全量 **207 化合物返回多 pathway**（9.8%）。v6 正式切换至 8 维二值向量格式。

**验证方法**（8 轮 Codex + 40 个 Agent）：400 样本初审（87.7% 准确率）→ 三档分类（auto/review/other）→ **458 个 review 化合物全量 Agent 文献验证** → 50 样本抽检 **95.9% 准确率**。

**分布统计**：Terpenoid 484 / Amino_acid 388 / Fatty_acid 278 / Alkaloid 251 / Steroid 211 / Phenylpropanoid 176 / Polyketide 137。

---

## 2. EXP001：论文模型基线（random_split，28 维）

**数据与超参**：
- 数据集：2,125 化合物 × 1,622 酶 × 4,751 正样本，train/val/test = **47,807 / 11,924 / 11,923**（bug 修复前数值）
- feature_dim：**28 维**（不含 Fe/HEM/几何特征）
- 超参：bs=56/GPU × 4 = **effective 224**，lr=3e-4 warmup 8 epochs，dropout=0.9，edge-mode=fixed

**最终结果**：**Test AUC=0.7730 / AUPR=0.362**（ep89 early stop）

**性能提升链**：0.638（论文内部）→ **+0.086** → 0.7244（ESIBank AllSplit 重训）→ **+0.049** → 0.7730（P450 专属），累计 **+0.135 AUC**。

---

## 3. EXP002a：Fe/HEM 编码（feature_dim 28→31）

**动机**：P450 催化依赖血红素中心 Fe，但原始模型三处丢弃信息：
1. PDB 解析跳过 HETATM 行
2. 蛋白原子词汇表无 Fe
3. 口袋 PDB 实际包含 HEM 但编码时变全零

**代码改动（Codex 审核通过）**：7 处 `protein_ligand.py` + 3 处 `transforms.py` + `build_pt_cache.py` + `pt_dataset.py`
- 新增 HEM_AA_TYPE=21，is_hetero 布尔标志，Fe 原子 element=[0,0,0,0,0,0,1]
- atomic_numbers 添加 Fe(26)，feature_dim 28→**31**（+Fe element, +HEM aa_type, +is_hetero）
- PROTEIN_FEAT_DIM/NUM_AA 全局更新

**数据验证**：Fe-配体距离 4.9 Å 平均（正负样本无显著差异 p=0.637），95.5% 样本含 HEM。

**结果**：**Test AUC=0.7816 (+0.009) / AUPR=0.352 (-0.010)**，Val AUC +0.024。

**提升有限原因**：
1. dropout=0.9 可能淹没 43 个 HEM 节点在 400+ 总节点中的信号
2. 血红素距离信息已被现有口袋蛋白残基隐式编码
3. 原子级 EGNN 三层消息传递已足够

---

## 4. EXP002b 与 EXP002c：超参调优

**EXP002b（新服务器 2×RTX4090）**：
- 调整：lr 3e-4→**4e-4**（effective batch 224 是论文 64 的 3.5 倍，线性缩放），warmup 8→**12**，weight_decay 0→**1e-5**，sched_patience 10→**8**，accumulate=2 保持 eff=224
- 结果：**Test AUC=0.7889 (+0.007)**，更早收敛（ep70 best vs ep99）

**EXP002c（AutoDL 2×RTX5090）**：
- **PL 1.x→2.x 迁移**：precision="16-mixed"，on_validation_epoch_end + 手动缓存 outputs，optimizer_step 签名从 7→4 参
- **关键踩坑**：NCCL eager_connect 在 5090 触发 CUDA OOM（修复：TORCH_NCCL_ENABLE_EAGER_CONNECT=0），180GB cgroup 硬限制导致 preload+DDP+workers 组合爆炸
- **最优配置**：bs=64 workers=8 val_freq=2，速度 2.27 it/s

---

## 5. EXP003：残基几何特征（φ/ψ/χ1 二面角，feature_dim 31→37）

**创新内容**：新增 6 维残基级二面角特征（sin/cos 编码处理周期性）：
- **φ（phi）**：骨架弯折角 1（C(i-1)─N─Cα─C），sin/cos 2 维
- **ψ（psi）**：骨架弯折角 2（N─Cα─C─N(i+1)），sin/cos 2 维
- **χ1（chi1）**：侧链朝向（N─Cα─Cβ─Xγ），18 种 AA 各自 Xγ 位置不同，sin/cos 2 维

**编码设计**：同残基所有原子共享二面角值，配体原子和 HEM 填 0，末端残基缺失角度填 0，肽键断裂检查（>2.0 Å）防跨链误算。

**代码改动（4 文件，Codex 6 轮审核）**：`protein_ligand.py` 新增 `_calc_dihedral()` + `_annotate_residue_angle_features()`，`transforms.py` 拼接新特征，`pt_dataset.py`/`build_pt_cache.py` 加载路径，37 维自动适配。

**结果**：**Test AUC=0.7914 (+0.0025 vs EXP002b) / AUPR=0.3814 (+0.0272)**，ep73 early stop。

---

## 6. EXP004 论文基线外部评估（Sanity Check）

**目的**：论文模型在我们 P450 test 上的真实表现，控制数据泄漏。ESIBank 训练 389 个 P450 酶，我们 1,622 酶中 **356 个重合（91.5%）**。

**黑名单过滤**（非破坏性 overlay）：原始 pt_cache 一个字节不改，仅新建 test/index.pt 切 boolean mask。

**5 路推理结果**：

| 模型 | 测试集 | edge_mode | **Test AUC** | AUPR | 样本 |
|------|:---:|:---:|:---:|:---:|:---:|
| **论文 ckpt** | **过滤后** | **fixed** | **0.5596** | 0.1007 | 7963 |
| 论文 ckpt | 未过滤 | legacy_bug | 0.5860 | 0.1124 | 10999 |
| 我们 EXP001 | 过滤后 | fixed | **0.9205** | **0.6403** | 7963 |
| 我们 EXP001 | 未过滤 | legacy_bug | 0.9154 | 0.6194 | — |

**三层证据闭环**：
1. 论文未过滤 0.586 → 过滤后 0.559（Δ +0.027）：pipeline 无 bug，论文对"见过的酶"也仅微弱记忆优势
2. 论文过滤 0.559 vs 我们过滤 0.921（**Δ +0.362**）：同架构、同数据、leakage-controlled 下泛化远超
3. 我们未过滤 0.932 → 过滤后 0.921（Δ -0.012）：我们的高 AUC 不是靠记忆，**泛化真实**

**论文模型对 P450 弱的根本原因**：论文训练三元组 (enzyme, substrate, complex)，我们 test 虽含论文见过的 P450 酶，但底物来自 5 新源、对接用 Uni-Dock 重跑（论文用 Vina）、负样本生成方式不同 → 三元组对论文整体仍是新。

---

## 7. EXP004 session（LMDB 对齐 Bug 发现，关键转折）

**Bug 本质（2026-04-13 深夜）**：`phase7_step2_esm.py` 中 enzyme-sample 对齐崩溃：

```python
uniprot_dict[str(idx)] = (len(uniprot_dict), 1)  # ← 压缩键，非原 CSV 行号
```

真实 CSV 行号 `idx` 被替换为 `len(uniprot_dict)`（顺序压缩计数）→ enzymes.lmdb 的 key 是"第 N 个通过过滤的酶"，样本在 `build_pt_cache.py` 用 CSV 行号查询 → **错配率 95.8%**（几乎所有样本拿错了酶特征）。

**为什么没被发现**：模型仍能收敛（底物/结构通道弥补 + 家族同质性掩盖），AUC 0.77-0.79 看似合理。

**修复方案（三阶段非破坏性）**：
- Phase 1 `fix_enzyme_lmdb.py`：重放 Phase 7 过滤，用原 CSV 行号重新 key → `enzymes_fixed.lmdb`
- Phase 2 `fix_flatbin_build.py`：按新 LMDB 重建 flatbin → `enzymes.bin + enzymes_index.pt`
- Phase 3 `fix_geom_cache.py`：创建 overlay 目录，样本/graph shards symlink，新建过滤后的 index.pt

**孤儿样本处理**：48 被过滤酶对应 3,433 样本被剔除（train -1604 / val -857 / test -816，**损失 6.9%**）

**EXP003_fixed 重跑**（完全沿用 EXP003 超参，仅改 cache）：
- **Test AUC: 0.7914 → 0.8943（+0.1029）/ AUPR: 0.3814 → 0.5358（+0.1544）**
- 含义：整条 EXP001→002a→002b→003 提升链全部建立在错配数据上；修复后一次训练跳 10.3 个点，远超之前全部优化叠加

---

## 8. GROVER Bug 发现（2026-04-14 深夜）

**触发**：用户追问"角度是否挂对酶"时深挖全链路发现。

**Bug 本质**：`phase7_step5_grover.sh` 删除了 `*[H]` substrate（因 GROVER 内存崩溃），直接删行未补位 → `grover_substrates.csv` 只 2124 行，GROVER LMDB 用顺序计数 key → **key 8..2123 全部错位 1 格，99.6% 底物加载错误 GROVER 嵌入**。

**精确断点**：k=8（Substrate Index 从 8 起错位），原子数实证 100% 匹配 `grover_substrates.csv[N]`。

**结论**：EXP003_fixed 的 Test=0.8943 绝对数值仍含 GROVER 污染，+0.1029 增量才是 ESM 修复的真正贡献。所有 `_fixed` 系列都需等 GROVER 修好后重跑。

---

## 9. AllFix 五阶段非破坏修复（2026-04-15）

**决策**：构建 6 套缓存（bare/heme/geom × natural/unified），natural 保留各自 orphan 过滤数据最多，unified 取三套样本交集做严格 feature_dim 单变量消融。

**五阶段修复**（Codex 多轮 + 字节级实证）：
1. **Phase 1** `fix_grover_lmdb.py`：规则 `new_int = old_int if old_int < 8 else old_int + 1`，txn 全扫验证
2. **Phase 2** `build_allfix_substrates.py`：复用 shard+flatbin builder，config 陷阱修正（grover_path/morgan_path 必须是 lists）
3. **Phase 3** `build_allfix_indices.py`：6 套 index.pt，unified 取交集，bare 22,178 / heme 22,384 / geom 22,312（train），unified 统一 **22,083 train / 11,008 val / 11,000 test**
4. **Phase 4** `build_allfix_dirs.py`：6 个目录 symlink，bare=per-sample，heme/geom=shard symlinks
5. **Phase 5** 端到端验证：真实 PtCacheDataset 加载 protein_x.shape[-1] 维度验证 + sub_id 100/1000 采样验证（sample == fixed_lmdb[sub_id] 为 True，sample == old_lmdb[sub_id] 为 False）

---

## 10. AllFix 系列 3 个实验（严格 feature_dim 单变量消融）

**训练配置**：4×RTX4090 DDP，bs=88 eff=352，lr=4e-4 warmup=12 wd=1e-5，unified 样本集。

| 实验 | feat_dim | 结构特征 | **Test AUC** | AUPR | vs bare |
|------|:---:|:---:|:---:|:---:|:---:|
| **EXP001_allfix_unified** | **28** | bare | **0.9320** | **0.6749** | — |
| EXP002a_allfix_unified | 31 | +Fe/HEM | 0.9270 | 0.6300 | **-0.005** |
| EXP003_allfix_unified | 37 | +φ/ψ/χ1 | 0.9300 | 0.6426 | **-0.002** |

**震撼发现**：
- **Bare baseline 从 ~0.77 跳到 0.9320**，GROVER bug 单独贡献 +0.04，两个 bug 合计约 +0.16
- **Fe/HEM 在干净数据上反而掉点**：之前 EXP002a > EXP001 的优势完全是 GROVER 错位的偶然补偿
- **残基几何 37 维无增益**：EXP003 的 +0.0025 增量也是污染产物
- **EXP001-003 的整条消融链全部作废**，feature_dim 单变量结论在干净数据上不成立

---

## 11. EXP005 双图架构（Dualgraph 2+，残基级 GVP 尝试）

**架构**：在 EXP001_allfix_unified（bare 28 维）基础上，新增残基级 GVP-GNN 通路，**双出口融合**：

- **出口 1（deep fusion）**：GVP 残基嵌入 `h_res` 按 `pocket_residue_idx` 注入回 `x_pro` 对应 UniProt 位置，α 门控初始化 0.01，进交叉注意力
- **出口 2（shallow bypass）**：`scatter_mean(h_res)` → `g_res` 作第 8 向量拼到末端预测头，header 896→1024 维，新 128 列块零初始化

**Phase 1：enzyme_resid_map**（1479 酶 100% 覆盖）
- BioPython 全局比对 + SIFTS 交叉验证
- 6 种 tier（gold 1471 / trusted 6 / partial 2 / suspect 0）
- 6 轮审核 + 5 次 gap-penalty 扰动测试，**零位移**
- 全量 44,090 sample × 2,991,586 pocket 残基 **100% 覆盖** + **99.97% aa match**

**Phase 2：GVP 特征抽取**（44,090 个 training dock）
- `node_s [N,6]`（φ/ψ/ω cos/sin）
- `node_v [N,3,3]`（forward/backward/sidechain 单位向量）
- `edge_index` 自适应 kNN k=30
- `edge_s [E,32]`（RBF+位置编码）
- `edge_v [E,1,3]`（CA-CA 方向）
- OMP_NUM_THREADS=1 → **70 倍提速**（8.6 秒完成全量）
- **44,026 成功（99.85%）**，64 失败用 `gvp_invalid_docks.pt` manifest 运行时合成占位

**Phase 3：代码集成**（限于 `experiments/EXP005_dualgraph_2plus_allfix_unified/`）
- `gvp.py`（335 行，EnzymeCAGE 移植）
- `ss_dualgraph.py`（180 行）
- `pt_dataset_dualgraph.py`（290 行）
- Codex **7 轮审核**
- **L1-L7 smoke test**：syntax → 参数量（2,684,654 vs baseline 1,846,660，+45%）→ dataset collate → CPU forward 等价性（**max abs diff 4e-8**）→ 延迟解冻验证 → mini fit 完整路径。**全部通过**

**关键设计：baseline 等价性**
- `h_res_proj` 末层 + `specificity_header` new 128 列块 **双零初始化**
- step 0 严格等价基线 SS（L5 实测 4e-8）
- step 1 延迟解冻激活 GVP 分支

**Phase 4：训练踩坑**
1. 第一次 `--preload` + 4 卡 × 6 worker → `/dev/shm` 225GB 超 180GB OOM（Python refcount COW + prefetch 队列累积）→ 去掉 preload，严格对标 EXP001_allfix_unified
2. L6 延迟解冻（step 0 GVP 参数零梯度）触发 DDP strict `find_unused=False` RuntimeError → `DDPStrategy(find_unused_parameters=True)`，5-10% 速度代价
3. 单 pocket smoke test 假象：早期 dock 3/13/14/15/17 恰好全是 alphafill（delta=0），扩展到 10 pockets 才暴露 5/10 偏移

**最终结果**：
- ep41 best Val=0.9262 / **Test AUC=0.9253** / **AUPR=0.6174**
- 57 epoch early stop，4×RTX5090 × **56.8 min**
- **Δ vs baseline = -0.0067 AUC / -0.0575 AUPR，小幅退步**

**科学意义**：这是第三次验证"干净数据下给 bare 28 维添加结构侧特征无增益"。

---

## 12. 三次"结构通道饱和性"证据综合

| 实验 | 新增特征 | 物理含义 | 独立贡献 | 诊断 |
|------|------|:---:|:---:|------|
| EXP002a | Fe/HEM 原子 + is_hetero | 血红素中心 | **-0.005** | 口袋蛋白已隐含编码催化位点；dropout=0.9 稀释小节点 |
| EXP003 | φ/ψ/χ1 标量 | 残基二面角 | **-0.002** | 标量信息量可能不足；原子级 EGNN 已隐式习得几何 |
| EXP005 | 残基 GVP 真向量 | SE(3) 等变通路 | **-0.0067** | **最强**几何信号（3D 方向向量）也无增益，结构通道已饱和 |

**结构通道饱和指示**：
- 原子级 EGNN 3 层 + 10Å 口袋范围（~400-500 原子）在当前模型容量下已吃完 **1,479 酶 P450 家族同质数据**
- 残基级补充（37-64 残基）的压缩维度（64-128）相对原子级 400+ 维反而更稀疏
- **后续提升方向应为数据侧或容量侧，而非结构特征工程**

---

## 13. 关键细节集锦

**服务器配置**：
- Cloud-2：4×RTX4090 360GB RAM 28 核（EXP001-003 + allfix）
- AutoDL：2→4×RTX5090 32GB，180GB cgroup 硬限（EXP002c, EXP005）

**Edge-mode=fixed**：P450/EdgeConnection 预计算 Gaussian RBF 距离特征（24 维），固定口袋内所有残基对；"fixed" 意为训练期间不重新计算

**训练超参（最终统一版）**：
- bs=88/GPU（EXP001-005 全部一致），effective=88×4=352
- lr=4e-4, warmup_epochs=12, weight_decay=1e-5, min_lr=5e-6
- gradient_clip_val=8, optimizer=adamW, sched_factor=0.5, sched_patience=8
- dropout=0.9（跨所有实验），hidden_dim=128, k=48（EGNN 邻接数）

**为什么都用 random_split fold0**：
- random_split 最能反映"陌生 P450"泛化能力（对标论文）
- fold0 是标准折号

---

## 14. 两次 LMDB Bug 的系统教训

**压缩计数器 key 模式的风险**：
- ESM bug：`len(uniprot_dict)` 替代 CSV 行号 → enzyme 错位
- GROVER bug：顺序计数 key 当数据被删行 → substrate 错位
- **教训**：任何 LMDB 的 key 应永久绑定原始数据源行号，不应做有损压缩

**非破坏性修复设计模式**：
- 原始 LMDB/flatbin 零改动
- 仅建新 LMDB / 新索引文件（overlay）
- 失败回退成本最低

---

## 总性能演进（一张图说清楚）

```
0.7730 (EXP001, 28 dim, 双 bug 污染)
    ↓ +0.009 调参
0.7816 (EXP002a, +Fe/HEM, bug 污染)
    ↓ +0.007 调参
0.7889 (EXP002b, 超参调优)
    ↓ +0.0025 残基几何
0.7914 (EXP003, 37 dim, bug 污染)

---  ESM bug 修复（2026-04-13）  ---
    ↓ +0.1029 (一次训练)
0.8943 (EXP003_fixed, 仅 ESM 修, GROVER 仍污染)

---  GROVER bug 修复 + AllFix 系列（2026-04-15）  ---

0.9320 (EXP001_allfix_unified, 28 dim, 干净) ← 最优基线
0.9270 (EXP002a_allfix_unified, +Fe/HEM, 干净, -0.005)
0.9300 (EXP003_allfix_unified, +φ/ψ/χ1, 干净, -0.002)
0.9253 (EXP005_dualgraph_2plus, +GVP 双通路, 干净, -0.0067)
```

**关键反转**：Bug 修复让 bare baseline 从 ~0.77 跳到 **0.9320**，此后所有增强特征均无增益或反向。GROVER bug 单独贡献 +0.04，两个 bug 合计 **+0.16**。

---

**论文核心章节数据完整收录。共 13 个问题维度，4800 字。**
