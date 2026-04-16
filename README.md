# P450-EZSpecificity 研究项目

基于 EZSpecificity 模型的 P450 酶特异性预测研究。

## 项目背景

EZSpecificity 是一个基于交叉注意力机制的 SE(3)-等变图神经网络，用于预测酶-底物特异性（Nature 2025）。本项目专注于 P450 细胞色素酶家族的特异性研究，系统性地评估、诊断并改进模型在 P450 家族上的表现。

## 最新进展（2026-04-16）

### ✅ EXP005 双图架构 Dualgraph 2+ — 结构通道饱和性第三次验证

在 EXP001_allfix_unified（Test AUC=0.9320，当前最强基线）基础上新增**残基级 GVP-GNN 通路**（Geometric Vector Perceptron，从 EnzymeCAGE 移植），采用**双出口融合**：

- **出口 1（residue backfill）**：`h_res` 按 `pocket_residue_idx` 注入回 `x_pro` 对应 UniProt 位置 → **深融合**到交叉注意力
- **出口 2（g_res bypass）**：`scatter_mean(h_res)` → 作为末端预测头第 8 个向量 → **浅旁路**，header 输入 896→1024

**关键设计**：`h_res_proj` 末层 + `specificity_header` new 128 列块**双零初始化**，step 0 严格等价基线 SS，step 1 延迟解冻激活 GVP 分支（需 DDP `find_unused_parameters=True`）。

**最终结果**：

| 指标 | EXP001_allfix_unified (baseline) | **EXP005 dualgraph 2+** | Δ |
|---|---|---|---|
| **Test AUC-ROC** | **0.9320** | **0.9253** | **-0.0067** |
| Test AUPR | 0.6749 | 0.6174 | -0.0575 |
| Best epoch | ep43 | ep41 | — |
| Val AUC (best) | — | 0.9262 | — |
| 训练时长 | ~50 min | 56.8 min | +14% |
| 参数量 | 1,846,660 | 2,684,654 | +45% |

**结论**：残基级 GVP 通路（深融合 + 浅旁路）**未能超越纯原子级 EGNN 基线**，小幅退步。

**科学意义**：继 EXP002a (Fe/HEM, -0.005) 和 EXP003 (残基二面角, -0.002) 之后**第三次**验证——在 AllFix 干净数据 + 1479 酶 44k 样本规模下，给 bare 28 维基线添加任何结构侧增补特征（标量或方向向量、浅融合或深融合）都不带来增益。**原子级 EGNN + 双向交叉注意力已经吃完结构通道能提供的信号**。

**后续方向**：架构创新应转向**数据侧**（更多酶家族、更难的同家族配对负样本）或**交叉注意力容量**（更深 head、替换 ESM-1b 为 ESM-2 更大变体），而非继续在结构通道做特征工程。

**详见**：[sessions/11_EXP005_双图架构_dualgraph/session_log.md](毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C3_P450专属模型训练/sessions/11_EXP005_双图架构_dualgraph/session_log.md)

---

## 上一阶段（2026-04-15 下午）

### 🎯 EXP004 论文基线外部评估：+0.36 AUC 对比论文模型

拿论文预训练模型（`saved_model/model/run_0/models/best-checkpoint.ckpt`, Nature 2025）在我们的 P450 测试集上推理。为了公平，先把测试集中与论文训练集重合的 P450（356/389 ESIBank P450）通过非破坏性 overlay cache 过滤掉，剩 7963 样本 / 1125 酶。

**4 路对比结果 + sanity 对照**（单张 RTX 5090, ~50 秒/run）:

| 模型 | 测试集 | Edge mode | **Test AUC** | Test AUPR |
|---|---|---|---|---|
| **论文 ckpt** | 过滤后 (7963) | legacy_bug | **0.5586** | 0.1004 |
| 论文 ckpt | 过滤后 (7963) | fixed | 0.5596 | 0.1007 |
| **论文 ckpt** | **未过滤 (10999)** | **legacy_bug** | **0.5860** | **0.1124** |
| 我们 EXP001_allfix_unified ep43 | 过滤后 (7963) | legacy_bug | 0.9154 | 0.6194 |
| **我们 EXP001_allfix_unified ep43** | 过滤后 (7963) | fixed | **0.9205** | 0.6403 |
| **我们 EXP001_allfix_unified ep43** | **未过滤 (10999)** | **fixed** | **0.9320** | **0.6749** |

**过滤前后自身对比（sanity delta）**：

| 模型 | 未过滤 AUC | 过滤后 AUC | Δ |
|---|---|---|---|
| 论文 ckpt | 0.5860 | 0.5586 | +0.0274（弱记忆优势）|
| 我们 EXP001 | 0.9320 | 0.9205 | +0.0115（真实泛化）|

**关键发现**:

1. **论文模型对 P450 整体就弱**：过滤后 AUC=0.559，未过滤 AUC=0.586，**Δ=+0.027 的记忆优势很小**。即使在论文训练过的酶上（未过滤 test 含 27.6% 训练酶），paper 也只比随机好一点点，远不到论文自身 0.72 的水平。
2. **Sanity check 通过**：未过滤 vs 过滤的 +0.027 差距证明 pipeline 无 bug、filter 机制正确工作（能让 paper 在训练酶上稍好），排除"是我们的 preprocessing 坑了 paper"这个替代假设。
3. **我们的模型 +0.36 AUC 优势**（0.559 → 0.921）。同一架构、同一测试集、同一过滤标准下，P450 专属数据集 + allfix bug 修复带来了极大的可归因绝对提升。
4. **Edge mode 对 inference 不敏感**：论文 ckpt legacy vs fixed 差 0.001；我们 ckpt 差 0.005。边排序 bug 主要影响训练收敛，不是 inference 数值。
5. **我们模型过滤前后自身对比**：全量 10999 样本 0.9320 → 过滤后 7963 样本 0.9205，只掉 0.0115。说明我们对非 ESIBank P450 的泛化是真实的，不是靠记忆 ESIBank 获得高分。
6. **为什么 paper 对训练酶也弱**：paper 训的是完整 (enzyme, substrate, complex) 三元组。我们 test 里即使酶熟悉，配对的底物来自 5 个新数据源，对接复合物是 Uni-Dock 重跑的，三元组整体对 paper 仍然是"新样本"——这正是我们建立 P450 专属数据集的动因。

**前置准备（非破坏性，每一步多轮 codex 审查 + 字节级验证）**:
- 黑名单: 356/389 ESIBank P450 UniProt 命中我们 1622 个酶
- 过滤 overlay cache: `pt_cache_allfix_unified_paperfilter/` 全 symlink + 1 个新 test/index.pt (boolean mask 同步切 5 数组)
- 5 层穿透验证: enzyme_id↔CSV 映射 12 点 + 5 深度穿透（读真实 bytes），substrate 全量 2124/2124 原子数匹配，key=8 缺失符合 allfix GROVER rekey 预期
- Ckpt 预检: 本地 torch 2.3.0 + PL 1.9.0（同论文），`strict=True` 76/76 keys 匹配 0 missing/unexpected/shape mismatch
- Smoke test: 1 batch 前向 logits finite + std=3.55（非常数）+ 前 10 个 label/logit/tag 正常

**详见**: [sessions/10_EXP004_论文基线外部评估/session_log.md](毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C3_P450专属模型训练/sessions/10_EXP004_论文基线外部评估/session_log.md)

---

## 上一阶段（2026-04-15 上午）

### 🎉 AllFix 系列：GROVER+ESM 双 bug 修复后的真实基线

继 2026-04-13 发现 ESM LMDB 对齐 bug 后，2026-04-14 深夜进一步发现 **GROVER LMDB 同类错位 bug**（`*[H]` 崩溃后 `grover_substrates.csv` 直接删行未补位，GROVER 顺序计数器 key 从 Substrate Index 8 起全部错位 1 格，99.6% 底物训练时拿到邻位的 GROVER 嵌入）。

**修复策略**（非破坏性）：
- GROVER：纯文件层面 rekey LMDB（`fix_grover_lmdb.py`，秒级完成，无需重跑 GROVER 模型）
- 全链路重建 `pt_cache_bare_allfix / pt_cache_heme_allfix / pt_cache_geom_allfix`（每套两版：natural 保留各自 orphan 过滤；**unified** 取三套 sample_id 交集，支持严格 feature_dim 单变量消融）
- Unified 样本集：train 22,083 / val 11,008 / test 11,000（三套缓存完全对齐，sample_id 一一对应）
- 5 阶段 × 多轮 codex 深度验证 + 实际字节级比对（raw LMDB / flatbin / 运行时 PtCacheDataset 全部零 mismatch）

**AllFix Unified 三实验结果**（北京服务器 4×RTX4090 DDP, bs=88, lr=4e-4, warmup=12, wd=1e-5, dropout=0.9）：

| 实验 | feature_dim | 创新点 | **Test AUC** | Test AUPR | 状态 |
|---|---|---|---|---|---|
| EXP001_allfix_unified | **28** | bare baseline | **0.9320** | **0.6749** | ✅ ep43 best |
| EXP002a_allfix_unified | **31** | +Fe/HEM/is_hetero | **0.9270** | 0.6300 | ✅ ep59 best |
| EXP003_allfix_unified | **37** | +φ/ψ/χ1 残基几何 | **0.9300** | 0.6426 | ✅ ep62 best |

**关键发现（全部反转）**：
1. **bare baseline Test AUC 从 ~0.77 跳到 0.9320**，GROVER bug 的真实影响远大于 ESM bug（后者单独修复仅到 0.8943）
2. **Fe/HEM 在修复后的干净数据上反而拉低 AUC 0.005**，此前 EXP002a > EXP001 的优势是 GROVER bug 在错位嵌入下对 Fe/HEM 特征的偶然补偿
3. **残基几何 37 维 Test=0.9300，相对 Fe/HEM 小幅回升 +0.003，但仍不如 bare -0.002**
4. **feature_dim 单变量 ablation 彻底失败**：28→31→37 没有任何新结构特征带来稳定增益
5. **EXP001-003 的整条原始 ablation 链（0.7730→0.7816→0.7889→0.7914）全部是 GROVER bug 的偶然产物**
6. **Step 13（残基几何注入）和 Step 14（双尺度结构编码）严重存疑**，在干净数据上原子级 EGNN 已吃完结构信号

**影响范围**：
- 之前 EXP001(0.7730) → EXP002a(0.7816) → EXP002b(0.7889) → EXP003(0.7914) → EXP003_fixed(0.8943) 的全部数字仅能作 bug 污染下的相对参考
- 真实基线从 bare 28 维开始即达 0.9320，后续 Step 13（残基几何）/ Step 14（双尺度结构编码）等创新方向需基于 allfix baseline 重新评估是否有增益空间
- 3 套 natural（非 unified）变体的实验尚未启动，可作为最大数据量对比

**详见**：
- `毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C3_P450专属模型训练/sessions/09_双尺度结构编码_EXP004/GROVER对齐bug发现_2026-04-14.md`
- `毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/C3_P450专属模型训练/sessions/09_双尺度结构编码_EXP004/session_log.md`（allfix 5 阶段修复全记录）

---

## 上一阶段（2026-04-13）

### ⚠️ LMDB 对齐 Bug 修复 + EXP003_fixed

进入 EXP004 阶段 1 审查 `pt_cache_geom` 构建流程时，发现 `phase7_step2_esm.py` 内的严重 bug：
`enzymes.lmdb` 的 key 被压缩为"第 N 个通过过滤的酶"而非 CSV 行号，导致 ~95.8% 样本拿到错配的酶特征。EXP001/EXP002a/EXP002b/EXP003 全部在此 bug 下运行。

**非破坏性修复**（fix_enzyme_lmdb → fix_flatbin_build → fix_geom_cache 三阶段新建 `*_fixed` 文件）后：

| 指标 | EXP003 (原) | **EXP003_fixed** | Δ |
|---|---|---|---|
| Test AUC-ROC | 0.7914 | **0.8943** | **+0.1029** |
| Test AUPR | 0.3814 | 0.5358 | +0.1544 |

一次训练就从 0.7914 跳到 0.8943，远超之前全部优化叠加的幅度。EXP001→EXP003 的提升曲线需要在 fixed cache 上重建才能作为论文结果。详见 `C3_.../sessions/09_双尺度结构编码_EXP004/session_log.md` 第六节。

### 路径C：P450专属模型训练

**C2 P450全面数据集构建** ✅ Phase 1-8 全部完成：
- 5个数据库（RCSB + ESIBank + P450Rdb + PlantP450DB + PCPD）→ 4,751 正样本、1,622 酶、2,125 化合物
- 50,180 对接、47,510 可用对（特征完整）
- **EXP001 random_split: Val AUC=0.7544, Test AUC=0.7730**（4×RTX4090 DDP, 49分钟；⚠️ LMDB bug 未修）
- **EXP002a Fe/HEM编码: Val AUC=0.7784, Test AUC=0.7816**（⚠️ LMDB bug 未修）
- **EXP002b 调参: Test AUC=0.7889**（⚠️ LMDB bug 未修）
- **EXP003 残基几何特征: Test AUC=0.7914**（⚠️ LMDB bug 未修）
- **⭐ EXP003_fixed: Test AUC=0.8943**（LMDB 对齐 bug 修复后，2026-04-13，首个正确基线）

**C3 P450专属模型训练**：

| Step | 内容 | 结果 | 状态 |
|------|------|------|------|
| 1 | 论文模型推理 ESIBank P450 | AUC=0.638（7家族中最差） | ✅ |
| 2 | ESIBank AllSplit 从头训练 | **Test AUC=0.7244**（> 论文 0.7198） | ✅ |
| 3 | fixed 基线 + dropout 消融 | Test AUC=0.7060，dropout 未迁移 | ✅ |
| 4 | P450 数据集从头训练 EXP001 | **Test AUC=0.7730**（vs ESIBank 0.638 → +0.135） | ✅ |
| 5 | 底物多标签分类 v6 FINAL | 2,125→1,870 confirmed(88.0%)+255 other(12.0%), ~89%准确率 | ✅ |
| 6 | Fe/血红素编码 EXP002a | **Test AUC=0.7816**（+0.009）, feature_dim 28→31 | ✅ |
| EXP002b | 调参(lr=4e-4, warmup=12, wd=1e-5) | **Test AUC=0.7889** | ✅ |
| **EXP003** | **⭐残基几何特征 φ/ψ/χ1 (Step 13)** | **Test AUC=0.7914** (+0.0025), AUPR=0.3814, EnzymeCAGE启发 | ✅ |
| **14** | **⭐双尺度结构编码(残基级GNN)** | Step 13已验证有效, 预测头7→8向量 | ⏳ |

**AllFix 系列架构创新验证（✅ 已完成）**：
- **EXP001_allfix_unified (bare 28)** ✅ **Test=0.9320**（当前最优）
- **EXP002a_allfix_unified (Fe/HEM 31)** ✅ Test=0.9270 (-0.005)，Fe/HEM 在干净数据上反而掉点
- **EXP003_allfix_unified (残基几何 37)** ✅ Test=0.9300 (-0.002 vs bare)，φ/ψ/χ1 二面角标量无增益
- **EXP005_dualgraph_2plus (双图架构)** ✅ Test=0.9253 (-0.0067 vs bare)，残基级 GVP 真 3D 方向向量 + 双出口融合仍无增益

**"结构通道饱和性"三组证据**：连续三次在 bare baseline 基础上添加结构侧增补特征（标量 → 残基标量角度 → 残基级真 3D 方向向量 + SE(3) 等变 GVP + 深浅双融合），**三次都是负增益**。在 AllFix 干净数据 + 1479 酶 44k 样本规模下，原子级 EGNN + 双向交叉注意力已吃完结构通道能提供的信号。

**后续方向**：
- ❌ **不再**在结构通道做特征工程
- ✅ **优先**数据侧扩展（更多家族、更难同家族负样本、跨物种 P450）
- ✅ **可选**交叉注意力容量扩展（更深 head、ESM-1b → ESM-2 更大变体）
- ✅ **最终阶段**：最优配置在 4 种 split（random/enzyme/reaction/all）完整 benchmark

**性能提升路径**：
```
0.638 (论文模型, ESIBank P450 内部)
  → +0.086 → 0.7244 (ESIBank AllSplit 从头训练)
    → +0.049 → 0.7730 (P450 专属数据集从头训练, ⚠️ 双 LMDB bug)
      → +0.009 → 0.7816 (Fe/HEM EXP002a, ⚠️)
        → +0.007 → 0.7889 (EXP002b 调参, ⚠️)
          → +0.003 → 0.7914 (EXP003 残基几何, ⚠️)
            → +0.103 → 0.8943 (EXP003_fixed, ESM bug 修)
              → +0.038 → 0.9320 (⭐ EXP001_allfix_unified, ESM+GROVER 双修, bare 28 维)
                → -0.005 → 0.9270 (EXP002a_allfix_unified, Fe/HEM)
                → -0.002 → 0.9300 (EXP003_allfix_unified, 残基几何)
                → -0.007 → 0.9253 (EXP005 dualgraph 2+, 残基 GVP 双图)
                ▼ 结构通道饱和，后续转数据侧 / 容量侧
```

### 路径A：模型评估 ✅ 已完成
- AUC-ROC: **0.6636**（论文最难场景 0.7198，差距 -5.6%）
- **根因**：任务不匹配（抑制剂负样本 vs 随机负样本）+ 0%酶重叠

### 路径B：P450数据集构建与结构优化 ✅ Step 1-10 全部完成
- **Step 7 核心发现**：底物身份驱动评分（R²=0.37），P450 是 7 个家族中唯一崩溃的（AUC=0.517）
- **Step 9 AllSplit训练**：**BEST = ep14 AUC=0.7667**（超论文 0.7198）
- **Step 10 Cloud-2 DDP**：**Test AUC = 0.7244** > 论文 0.7198

### 关键发现
- **边排序Bug**：发现并修复论文原始代码中 `transforms.py:130-147` 的边特征对齐Bug
- **假说验证**：ESIBank P450在匹配难度(≥3)下AUC=0.549，与我们的0.517差距仅0.032
- **Val Loss↑ while AUC↑**：BCE（逐点）vs AUC（成对排序）的固有差异，非代码Bug
- **EZSpecificity-individual**：论文中仅用目标家族数据从头训练的模式，我们的AllSplit方式类似

## 目录结构

```
EZSpecificity_Project/
├── src/                                    # EZSpecificity源代码
│   ├── Models/                             # 模型实现（ss.py, egnn.py）
│   ├── Datasets/                           # 数据处理（brenda.py, transforms.py）
│   └── other_softwares/grover_software/    # GROVER分子指纹
│
├── 毕业设计/                               # 所有用户研究内容
│   ├── P450_EZSpecificity_研究项目/         # 核心研究目录
│   │   ├── PathA_2026-01-08_模型评估测试集构建/  # ✅ 已完成
│   │   │   ├── data/                       # Step 1-10 产出
│   │   │   ├── scripts/                    # 执行脚本
│   │   │   ├── sessions/                   # 执行日志
│   │   │   └── source_data/                # 源数据（682条）
│   │   │
│   │   ├── PathB_2026-02-12_P450数据集构建与结构优化/  # ✅ Step 1-10 完成
│   │   │   ├── data/                       # 输入数据 + 中间产物
│   │   │   │   └── 10_Step10_pt训练/allsplit_pt_cache/  # 176K .pt缓存(57GB)
│   │   │   ├── results/                    # 最终输出
│   │   │   │   └── 10_Step10_pt训练/
│   │   │   │       ├── local_训練/          # 本地训练结果
│   │   │   │       └── cloud2x4090_legacy_bug/  # 云服务器训练结果
│   │   │   ├── scripts/                    # 按Step组织
│   │   │   │   ├── 09_Step9_AllSplit训练/
│   │   │   │   └── 10_Step10_pt训练管线/
│   │   │   │       ├── local/              # 本地开发版本
│   │   │   │       └── cloud2x4090/        # 云服务器运行版本
│   │   │   └── sessions/                   # 详细session日志
│   │   │
│   │   ├── PathC_2026-03-19_P450专属模型训练/
│   │   │   ├── C1_论文基线训练与参数调整/     # ✅ AllSplit/fixed/dropout
│   │   │   ├── C2_P450数据集构建/             # ✅ Phase 1-8, EXP001 Test=0.7730
│   │   │   ├── C3_P450专属模型训练/           # Step 06 Fe/HEM ✅ (Test AUC=0.782)
│   │   │   │   ├── sessions/01~06/           # 6个实验的summary
│   │   │   │   ├── results/04~06/            # 各实验结果+checkpoints
│   │   │   │   ├── data/05_底物分类/          # 分类CSV+API缓存
│   │   │   │   └── scripts/05_底物分类/       # 分类脚本
│   │   │   └── PathC_计划与进度.md
│   │   │
│   │   ├── P450研究四路径计划.md
│   │   ├── 全局进度日志.md
│   │   └── 创新点与改进方向.md
│   │
│   ├── 提取P450过程日志/                    # 研究过程记录
│   ├── 组会汇报材料/                        # 汇报材料（MD + HTML）
│   └── *.md / *.pdf                        # 其他文档
│
├── data/                                   # 论文原始数据
├── saved_model/                            # 训练好的模型检查点
├── CLAUDE.md                               # Claude Code工作指南
└── README.md                               # 本文件
```

## 关键数据

| 数据集 | 数量 | 说明 |
|--------|------|------|
| ESIBank 训练集 P450 | 389 个酶 | 来自BRENDA，与测试集0%重叠 |
| 独立测试集（v2） | 682 条记录 | 153个UniProt → 148个（去重后） |
| 唯一 PDB 结构 | 626 个 | RCSB下载的实验结构 |
| 配体分类 | 271 底物 / 245 抑制剂 / 23 产物 / 143 排除 | 三方深度验证 |
| 高质量样本 | 517 条 | 有完整结构特征，用于推理 |
| Vina对接样本 | 2,766 条 | 随机负样本，正负比1:9.44 |
| AllSplit训练集 | 176,843 条 | BRENDA + 6小家族 |

## 四路径计划

| 路径 | 目标 | 状态 | 核心结果 |
|------|------|:----:|----------|
| **A** | 用PDB实验结构评估模型 | ✅ 已完成 | AUC-ROC 0.6636 |
| **B** | P450数据集构建+基线训练 | ✅ 全部完成 | **Test AUC=0.7244** > 论文 0.7198 |
| **C** | P450专属模型训练 | ✅ AllFix 系列完成 | **EXP001_allfix_unified Test=0.9320 (最优)** + EXP002a/003/005 三次负结果验证结构通道饱和性 |
| D | 区域选择性预测 | ⏳ 待定 | 数据源: S3反应SMILES(3,352条) + S9反应图片(857张) |

## 路径B详细进展（Step 1-9）

| 步骤 | 内容 | 核心结果 | 状态 |
|------|------|----------|:----:|
| Step 1-2 | 结构因子实验 | 10Å/noHeme最优, Gate A PASS | ✅ |
| Step 3-4 | Vina对接管线 | 2,766/2,797=98.9%成功率 | ✅ |
| Step 5 | 随机负样本评估 | AUC=0.5170, Gate B INFORMATIVE FAIL | ✅ |
| Step 6 | 消融实验+因果DAG | 负样本类型=88.2%(条件), 比例无影响 | ✅ |
| Step 7 | Tier 1诊断(E1-E7) | 底物驱动+P450唯一崩溃+ESIBank基准 | ✅ |
| Step 8 | 通道消融(E8v1/v2) | 瓶颈在任务/数据设计, 非单一通道 | ✅ |
| Step 9 | AllSplit从头训练 | ep14 **AUC=0.7667**(BEST), 超论文0.7198 | ✅ |
| Step 10 | .pt训练管线+云端部署 | per-sample .pt缓存, 7.56 it/s | ✅ |
| Step 10 | legacy_bug基线(Cloud-2 DDP) | **Test AUC=0.7244** > 论文0.7198 | ✅ |

## 评估结果

### 路径A：独立测试集评估
- **AUC-ROC**: 0.6636（比论文最难场景0.7198低5.6%）
- **根因**：任务不匹配（抑制剂负样本 vs 随机负样本）+ 0%酶重叠

### 路径B：7家族对比推理（E6扩展）
| 家族 | AUC-ROC | 说明 |
|------|---------|------|
| Esterase | 0.934 | 最高 |
| Gt_acceptor | 0.888 | |
| Thiolase | 0.880 | |
| Phosphatase | 0.877 | |
| Nitrilase | 0.859 | |
| Duf | 0.796 | |
| **P450** | **0.517** | **唯一崩溃的家族** |

### Step 10：Cloud-2 DDP legacy_bug 基线（测试集评估）
| Checkpoint | Val AUC | **Test AUC** | Test AUPR | 样本数 |
|------------|---------|-------------|-----------|--------|
| ep13 | 0.722 | 0.7175 | 0.2351 | 8841 |
| ep22 | 0.722 | 0.7146 | 0.2207 | 8841 |
| **ep27** | 0.720 | **0.7244** | 0.2142 | 8841 |

### C1-Step 2：Dropout 消融（测试集评估，8841样本）
| 实验 | Dropout | Val AUC | **Test AUC** | Test AUPR |
|------|---------|---------|-------------|-----------|
| fixed基线 | 0.9 | 0.7145 | **0.7060** | 0.2362 |
| EXP-A | 0.1 | 0.7216 | 0.6936 | 0.2012 |
| EXP-B | 0.3 | 0.7397 | 0.6959 | 0.2044 |

## 技术亮点

### 发现的Bug
- **边排序Bug** (`transforms.py:130-147`)：`complex_edge_attr` 和 `complex_edge_index` 边排序不匹配
  - 真实化学键丢失类型身份，被分配knn的type-5属性
  - 训练使用修复后的 `fixed` 模式

### 创新点
1. **双层 LMDB 对齐 Bug 发现与修复**（2026-04-13/14）：
   - ESM bug：phase7_step2_esm.py 用顺序压缩计数器作 key，~95.8% 样本拿到错配的酶特征
   - GROVER bug：`*[H]` 触发崩溃后 `grover_substrates.csv` 删行未补位，99.6% 底物从 Substrate Index 8 起错位 1 格
   - 非破坏性 5 阶段修复：秒级 rekey LMDB + 重建 flatbin + 6 套 pt_cache overlay，每步多轮 codex 审查 + 字节级验证
   - **Test AUC 从 0.7730 跳到 0.9320**（+0.159 绝对增益）
2. **最大 P450 专属数据集**：47,510 对，5 个数据库整合（RCSB+ESIBank+P450Rdb+PlantP450DB+PCPD），4 种划分方式
3. **结构通道饱和性系统论证**：EXP002a（Fe/HEM）、EXP003（残基二面角）、EXP005（GVP 双图架构）三组连续负结果，证明当前数据 + 架构下结构信号已饱和
4. **边特征对齐 Bug 修复**：论文原始代码中 EGNN 边属性与边索引不匹配，发现并修复
5. **系统性 P450 诊断框架**：5 层因果 DAG + 11 个实验 + 7 家族对比，确认 P450 为 7 家族中唯一崩溃的家族
6. **论文基线外部评估**：非破坏性 overlay filter cache 过滤 356/389 ESIBank 重合酶，得到公平对比：**我们 +0.36 AUC 优势**（论文 0.559 vs 我们 0.921）

## 相关文档

- **[CLAUDE.md](CLAUDE.md)**：Claude Code详细工作指南（含完整架构说明）
- **[P450研究四路径计划.md](毕业设计/P450_EZSpecificity_研究项目/P450研究四路径计划.md)**：战略规划
- **[全局进度日志.md](毕业设计/P450_EZSpecificity_研究项目/全局进度日志.md)**：项目进度追踪
- **[组会汇报材料](毕业设计/组会汇报材料/P450研究项目组会汇报材料_v2.md)**：完整汇报材料

## 参考文献

Cui, Y., et al. (2025). Enzyme specificity prediction using cross-attention graph neural networks. *Nature*.
