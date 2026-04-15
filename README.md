# P450-EZSpecificity 研究项目

基于 EZSpecificity 模型的 P450 酶特异性预测研究。

## 项目背景

EZSpecificity 是一个基于交叉注意力机制的 SE(3)-等变图神经网络，用于预测酶-底物特异性（Nature 2025）。本项目专注于 P450 细胞色素酶家族的特异性研究，系统性地评估、诊断并改进模型在 P450 家族上的表现。

## 最新进展（2026-04-15）

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

**架构创新方向（EnzymeCAGE启发）**：
- **Step 13**：从口袋PDB提取残基级几何特征（骨架二面角φ/ψ + 侧链扭转角χ1），sin/cos编码后拼接到EGNN节点特征。EGNN代码不动，改动模式与Fe/HEM编码一致
- **Step 14**：新增残基级GNN通道，与现有原子级EGNN互补（原子级看精细接触，残基级看口袋几何形状）
- **最终阶段**：最优配置确定后在4种split(random/enzyme/reaction/all)上完整benchmark

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
| **C** | P450专属模型训练 | 🔄 C3 Step 06✅ → Step 13⏳ | EXP002a Test=0.7816, EXP002b调参中, **Step 13/14架构创新(EnzymeCAGE启发)** |
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
1. **边特征对齐Bug修复**：论文原始代码中EGNN边属性与边索引不匹配
2. **Fe/血红素编码** ✅：将Fe(26)加入蛋白原子词汇表 + HEM残基类型 + is_hetero标志位，feature_dim 28→31，Test AUC +0.009
3. **残基级几何特征注入**（Step 13）：φ/ψ/χ1二面角→EGNN节点，捕获口袋腔体形状信息
4. **双尺度结构编码**（Step 14）：原子级EGNN + 残基级GNN双通道，EnzymeCAGE启发
5. **系统性P450诊断框架**：5层因果DAG + 11个实验 + 7家族对比
6. **最大P450专属数据集**：47,510对，5个数据库整合，4种划分方式

## 相关文档

- **[CLAUDE.md](CLAUDE.md)**：Claude Code详细工作指南（含完整架构说明）
- **[P450研究四路径计划.md](毕业设计/P450_EZSpecificity_研究项目/P450研究四路径计划.md)**：战略规划
- **[全局进度日志.md](毕业设计/P450_EZSpecificity_研究项目/全局进度日志.md)**：项目进度追踪
- **[组会汇报材料](毕业设计/组会汇报材料/P450研究项目组会汇报材料_v2.md)**：完整汇报材料

## 参考文献

Cui, Y., et al. (2025). Enzyme specificity prediction using cross-attention graph neural networks. *Nature*.
