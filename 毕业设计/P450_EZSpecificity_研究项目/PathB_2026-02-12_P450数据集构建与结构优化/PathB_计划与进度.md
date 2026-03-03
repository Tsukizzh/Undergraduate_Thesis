# Path B：P450 数据集构建与结构优化

## 一、背景与目标

### 1.1 前置条件

Path A 已完成（AUC-ROC 0.6636），根因分析识别出三个核心问题：

| 根因 | 描述 | 影响 |
|:---:|------|------|
| 1 | P450 泛化不足：测试集 148 个 P450 与训练集 389 个完全不重叠（0%） | AUC 从 0.89 降至 ~0.72 |
| 2 | 负样本语义差异：训练用随机配对，测试用抑制剂（更难区分） | AUC 从 ~0.72 降至 0.66 |
| 3 | 结构信息缺失：Heme 被排除、10Å 口袋可能不够精确 | 加剧泛化失败 |

### 1.2 Path B 的总目标

**从多个维度构建高质量 P450 数据集，为 Path C（模型微调）提供数据基础。**

具体包含三个方案：
- **方案 3**：结构特征改进（Heme、口袋半径）→ 解决根因 3
- **方案 1**：随机配对负样本 + 分子对接 → 验证根因 2
- **方案 2**：外部 P450 数据库数据收集 → 解决根因 1

### 1.3 导师反馈（2026-02 组会）

1. 去除产物记录，仅保留底物 + 抑制剂（B3 → B6 方案）
2. 尝试纳入 Heme 辅因子
3. 同时测试 6Å 和 10Å 口袋半径
4. 三个改进方案都要推进

---

## 二、总体架构

### 2.1 三方案线性流程

```
方案3（结构特征改进）──→ 方案1（随机负样本）──→ 方案2（外部数据库）──→ 交付 Path C
     Step 1-2              Step 3-5              Step 6-8
```

方案之间有明确的依赖关系：

```
方案3 的结果（最佳结构参数）
  → 供方案1 的对接使用（用最佳参数生成结构特征）
  → 供方案2 的大规模处理使用

方案1 的基础设施（对接流水线）
  → 供方案2 直接复用（批量对接新数据）

方案2 的产出（完整 P450 数据集）
  → 交付 Path C（模型微调）
```

### 2.2 决策门机制

| 决策门 | 位置 | 决策内容 | 后果 |
|:-----:|:----:|---------|------|
| **Gate A** | Step 2 完成后 | 选定最佳结构参数（Heme on/off、6Å/10Å） | 后续所有步骤使用该参数 |
| **Gate B** | Step 5 完成后 | 确认随机负样本策略的有效性 | 决定方案2是否沿用相同负样本策略 |
| **Gate C** | Step 8 完成后 | 最终数据质量验证 | 确认数据集可交付 Path C |

### 2.3 共享资源设计

**核心原则**：昂贵特征只算一次，新数据只算增量。

```
data/00_shared/
├── datasets/          # 数据集版本（互不覆盖）
│   ├── B6_v1/         # 去产物后的 516 条（Step 1 产出）
│   └── B6_expanded_v2/# 未来扩增后（Step 6-7 产出后更新）
│
└── features/          # 昂贵特征（从 PathA 复用，可追加）
    ├── enzyme_features.lmdb      # ESM 嵌入（按 enzyme_id 索引）
    ├── reaction_features.lmdb    # 反应图特征（按 substrate_id 索引）
    ├── morgan_fingerprint.npy    # Morgan 指纹
    ├── grover_fingerprint.lmdb   # GROVER 指纹
    └── README.md                 # 来源、版本、覆盖范围
```

- **数据集切换**：修改 `configs/base.yaml` 中的 `dataset_dir` 字段即可
- **新数据追加**：只为新增的酶/底物计算特征，append 到现有 LMDB

### 2.4 实验管理

含多组实验的 Step（如 Step 2、Step 5）内部按实验编号隔离：

```
data/02_Step2_因子实验/
├── EXP01_B6_10A_noHeme/    # 每组实验独立目录，永不覆盖
│   ├── config.yaml          # 参数快照（自动生成）
│   ├── structure_features/
│   ├── predictions.csv
│   └── evaluation/
├── EXP02_B6_10A_Heme/
├── EXP03_B6_6A_noHeme/
└── EXP04_B6_6A_Heme/
```

跑新实验：`python scripts/utils/run_experiment.py --name EXP05 --include_heme --pocket_radius 15`
→ 自动创建新目录，旧实验完整保留。

---

## 三、Step 详细计划

### Step 1：数据准备与代码改造

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 3（结构特征改进）|
| **前置依赖** | Path A 完成 |
| **目标** | 准备 B6 数据集 + 改造结构特征提取代码 |

**子任务**：

1. **B6 数据集迁移**
   - Path A Step 4 已生成 B6 数据集（516 条 = 272 底物正 + 244 抑制剂负）
   - 来源：`PathA_.../data/04_Step4_格式修正后数据/B6_底物正抑制剂负_271pos245neg/`
   - 复制 data.csv + Enzymes.csv + Substrates.csv 到 `data/00_shared/datasets/B6_v1/`
   - 验证记录数和标签分布正确

2. **共享特征迁移**
   - 从 PathA 复制 Step 5-7 的特征文件到 `data/00_shared/features/`
   - 验证 LMDB 完整性（key 数量、数据类型）
   - ~~编写 `features/README.md` 记录来源~~（已在 session_log 中记录）

3. **代码改造**
   - 修改 `extract_pocket_ligand.py`：
     - 新增 `--include_heme` 参数（布尔值）
     - 新增 `--pocket_radius` 参数（浮点数，默认 10.0）
     - Heme 纳入逻辑：从 HETATM 记录中提取 HEM 残基的非氢原子
   - 编写 `run_experiment.py`：读取 config → 创建实验目录 → 调用 Step 8-10 子流程

**产出**：
- `data/00_shared/datasets/B6_v1/`（B6 数据集）
- `data/00_shared/features/`（共享特征库）
- `scripts/utils/run_experiment.py`（实验运行器）
- `scripts/01_Step1_数据准备/`（数据过滤脚本）

---

### Step 2：结构特征 2×2 因子实验

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 3（结构特征改进）|
| **前置依赖** | Step 1 完成 |
| **目标** | 确定最佳结构参数组合 |

**实验矩阵**（2×2 因子设计）：

| 实验 | Heme | 口袋半径 | 目的 |
|:---:|:---:|:---:|------|
| EXP01 | 无 | 10Å | 基线（对照组，与 PathA 条件一致） |
| EXP02 | **有** | 10Å | 单独测试 Heme 效果 |
| EXP03 | 无 | **6Å** | 单独测试口袋缩小效果 |
| EXP04 | **有** | **6Å** | 测试组合效果 |

**每组实验的执行流程**：
```
Step 8.1：口袋/配体提取（使用对应 Heme/半径参数）
Step 8.2：配体原子对齐
Step 8.3：生成 structure_features.lmdb
Step 9：  模型推理（使用 00_shared/ 的共享特征 + 本实验的结构特征）
Step 10： 评估分析（AUC-ROC、AUC-PR、得分分布）
```

**⚠️ Step 2 前置修复（2026-02-20 Codex 四轮验证确认）**：

Step 1 的 `extract_pocket_ligand.py` 已支持 `--include_heme`，但下游 `step8_generate_structure_lmdb.py` 的 `PDBProtein` 解析器**只读 ATOM 记录，会静默丢弃 Heme 的 HETATM 行**。若不修复，Heme 因子实际无效，2×2 退化为仅测半径。

**修复 TODO**（Step 2 执行前必须完成）— **✅ 全部完成（2026-02-21）**：
1. ✅ 创建 PathB 版 `step8_generate_structure_lmdb.py`，`PDBProtein._parse()` 支持 `ATOM` + `HETATM`
2. ✅ atom dict 增加 `record_type` 字段（防极端情况下的 ID 冲突）
3. ✅ 增加 HETATM 统计日志（n_atom, n_hetatm 等）
4. ⏳ 门禁校验待实验运行后验证（EXP01 HETATM=0，EXP02/04 > 0）
5. ⏳ 同一 Dock_Index 对比待实验运行后验证

**额外完成**（Codex 三轮审核修复 7 个 bug）：详见 `sessions/02_Step2_.../session_log.md`

**科学定位**：修复后本实验为「推理敏感性/结构扰动实验」，非「Heme 生化机理验证」。原因：
- 模型训练时未见 HETATM（PathA 539 个 pocket 中 HETATM=0，已全量验证）
- Fe(26) 不在 `FeaturizeProteinAtom.atomic_numbers` 中（全零向量），Heme 的 C/N/O 可正常编码
- 结论措辞应限定为"模型对 Heme 几何节点注入的敏感性"

**分析重点**：

| 对比 | 分析目的 |
|------|---------|
| EXP02 vs EXP01 | Heme 的独立贡献 |
| EXP03 vs EXP01 | 口袋半径的独立贡献 |
| EXP04 vs EXP01 | 二者组合效果 |
| 4 组交叉对比 | 是否存在交互效应 |

**注意事项**：
- 模型在 10Å 无 Heme 条件下训练，改变输入分布可能导致 OOD（分布外）效应
- 如果 EXP02/EXP03/EXP04 都不如 EXP01，说明模型无法利用这些新信息（需要微调才能受益）
- 如果某组明显提升，后续所有步骤采用该参数

**⛳ Gate A 决策**：
- 选定后续使用的 `include_heme` 和 `pocket_radius` 值
- 将决策记录在 `sessions/02_Step2_结构特征因子实验/decision.md`
- 更新 `configs/base.yaml` 中的默认结构参数

---

### Step 3：对接环境搭建

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 1（随机负样本）|
| **前置依赖** | Gate A 通过 |
| **目标** | 搭建 AutoDock Vina 自动化对接流水线 |

**已有资源**：
- AutoDock Vina 1.2.7（`D:\autodock\vina.exe`）
- MGLTools 1.5.7（含 `prepare_receptor4.py`）
- NVIDIA 4070 Super GPU（可选 Vina-GPU 加速）

**子任务**：

1. **流水线开发**
   - 蛋白准备自动化：PDB → 去水/加氢 → PDBQT
   - 配体准备自动化：SMILES → 3D 构象 → PDBQT
   - Grid box 自动定位：利用 HEM 残基 Fe 原子坐标作为中心
   - Vina 调用封装：config 自动生成 + 批量调用

2. **预实验验证**（5 对）
   - 选 5 对已知的酶-底物对（从 PathA 正样本中选）
   - 端到端跑通：准备 → 对接 → 评分 → 提取位姿
   - 验证：对接位姿是否合理、结合能分布

3. **扩大验证**（50 对）
   - 抽检对接结果质量
   - 估算全量时间（2,537 对）
   - 决策：CPU 并行 vs Vina-GPU

**产出**：
- `scripts/03_Step3_对接环境/` 下完整的对接流水线脚本
- `data/03_Step3_对接预实验/` 预实验结果
- 对接 SOP 文档

---

### Step 4：随机负样本生成与批量对接

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 1（随机负样本）|
| **前置依赖** | Step 3 完成 |
| **目标** | 生成与训练集分布一致的随机配对负样本 |

**子任务**：

1. **负样本生成策略**
   - 保持 271 条正样本不变
   - 随机配对负样本：~2,537 对（正负比 1:9.36，与 ESIBank 训练集一致）
   - 采样规则：每个底物随机配对非催化酶（单向采样，避免假负样本）

2. **批量分子对接**
   - 使用 Step 3 的流水线对 ~2,537 对进行对接
   - 使用 Gate A 选定的结构参数（Heme、口袋半径）
   - 并行策略：多核 CPU 或 Vina-GPU

3. **质量控制**
   - 检查对接成功率
   - 分析结合能分布（排除异常高结合能的"假负样本"）

**产出**：
- `data/04_Step4_批量对接/` 下所有对接复合物 + 评分
- 更新 `data/00_shared/datasets/` 下的负样本数据集

---

### Step 5：重构评估分析

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 1（随机负样本）|
| **前置依赖** | Step 4 完成 |
| **目标** | 验证随机负样本策略对 AUC-ROC 的影响 |

**子任务**：

1. **特征生成**
   - 对新增的随机负样本生成所有特征（ESM、GROVER、Morgan、Structure）
   - 追加到 `data/00_shared/features/`

2. **模型推理**
   - 使用 `run_experiment.py` 运行实验
   - 对比配置：与 Step 2 的最佳基线对比

3. **根因验证分析**
   - 对比 AUC-ROC：抑制剂负样本 vs 随机负样本
   - 量化根因 2 的贡献（预期：AUC 从 ~0.66 升至 ~0.72）
   - 分析得分分布变化

**⛳ Gate B 决策**：

| 结果 | 解读 | 后续行动 |
|------|------|---------|
| AUC 升至 ~0.72 | 负样本语义差异确认为主因 | 方案 2 沿用随机负样本策略 |
| AUC 几乎不变 | 根因主要在泛化不足 | 方案 2 需更大规模训练数据 |

---

### Step 6：外部数据库收集与清洗

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 2（外部数据库）|
| **前置依赖** | Gate B 通过 |
| **目标** | 从专业 P450 数据库收集训练数据 |

**数据源**：

| 数据库 | 预期规模 | 数据类型 | 优先级 |
|--------|---------|---------|:------:|
| Plant P450 Database (erda.dk) | 874 条 | 酶-底物对（文献验证） | 最高 |
| P450Rdb (cellknowledge.com.cn) | 342 条 / 1,692 reactions | 反应数据 | 次优 |
| PCPD (p450.biodesign.ac.cn) | 181 条 | 序列 + 结构 + 配体 | 备选 |

**子任务**：

1. **数据下载与解析**
   - 编写各数据库的解析脚本
   - 提取酶序列（FASTA）+ 底物 SMILES

2. **数据清洗与去重**
   - 统一 CYP 命名格式
   - 标准化 SMILES（RDKit Canonicalize）
   - 与 ESIBank 训练集去重（确保独立性）
   - 与 PathA 测试集去重

3. **负样本生成**
   - 采用 Gate B 确认的负样本策略
   - 正负比 1:9.36（与训练集一致）

**产出**：
- `data/06_Step6_外部数据/` 下清洗后的合并数据集
- 数据来源追踪表

---

### Step 7：结构获取与质量控制

| 项目 | 内容 |
|------|------|
| **所属方案** | 方案 2（外部数据库）|
| **前置依赖** | Step 6 完成 |
| **目标** | 为新数据获取 3D 结构并完成对接 |

**子任务**：

1. **酶结构获取**
   - 优先查找 RCSB PDB 实验结构
   - 无实验结构时使用 AlphaFold DB 预测结构
   - 质量过滤：实验结构要求分辨率 ≤3.0Å，AlphaFold 要求 pLDDT ≥70

2. **批量分子对接**
   - **直接复用 Step 3-4 的对接流水线**
   - 使用 Gate A 选定的结构参数
   - 对所有新增酶-底物对进行对接

3. **质量控制**
   - 对接成功率统计
   - 结合能分布检查
   - 与 PathA 数据的结构质量对比

**产出**：
- `data/07_Step7_结构库/` 下完整的 PDB/PDBQT 文件
- 质量控制报告

---

### Step 8：全量特征生成与数据集交付

| 项目 | 内容 |
|------|------|
| **所属方案** | 整合（方案 1+2+3）|
| **前置依赖** | Step 7 完成 |
| **目标** | 生成最终数据集，交付 Path C |

**子任务**：

1. **增量特征生成**
   - 对新增酶序列生成 ESM 嵌入（追加到 LMDB）
   - 对新增底物生成 GROVER + Morgan 指纹
   - 使用 Gate A 参数生成结构特征

2. **数据集构建**
   - 合并 PathA 数据 + Step 4 随机负样本 + Step 6-7 外部数据
   - 训练集 / 验证集 / 测试集划分
   - 存入 `data/00_shared/datasets/` 新版本

3. **回归验证**
   - 在新数据集上重跑基线模型
   - 确认数据质量（无异常得分、特征完整性）

**⛳ Gate C 决策**：
- 数据集质量是否达标？
- 是否可以交付 Path C 进行微调？

**产出**：
- `data/08_Step8_最终数据集/` 下完整的训练数据集
- 数据集统计报告（物种分布、正负比、特征覆盖率）
- `data/00_shared/datasets/` 更新为最终版本

---

## 四、目录结构总览

```
PathB_2026-02-12_P450数据集构建与结构优化/
│
├── sessions/                            # 工作日志（写报告用）
│   ├── 01_Step1_数据准备与代码改造/
│   ├── 02_Step2_结构特征因子实验/
│   │   ├── overview.md
│   │   ├── EXP01_log.md ... EXP04_log.md
│   │   └── decision.md                 # Gate A 决策记录
│   ├── 03_Step3_对接环境搭建/
│   ├── 04_Step4_随机负样本对接/
│   ├── 05_Step5_重构评估分析/
│   │   └── decision.md                 # Gate B 决策记录
│   ├── 06_Step6_外部数据库收集/
│   ├── 07_Step7_结构获取质量控制/
│   └── 08_Step8_全量特征与数据集交付/
│       └── decision.md                 # Gate C 决策记录
│
├── data/
│   ├── 00_shared/                       # 跨步骤共享资源
│   │   ├── datasets/                    # 数据集版本
│   │   └── features/                    # 昂贵特征（ESM/GROVER/Morgan）
│   ├── 01_Step1_数据准备/
│   ├── 02_Step2_因子实验/               # 含 EXP01-EXP04 子目录
│   ├── 03_Step3_对接预实验/
│   ├── 04_Step4_批量对接/
│   ├── 05_Step5_重构评估/               # 含 EXP 子目录
│   ├── 06_Step6_外部数据/
│   ├── 07_Step7_结构库/
│   └── 08_Step8_最终数据集/
│
├── scripts/
│   ├── 01_Step1_数据准备/
│   ├── 02_Step2_因子实验/
│   ├── 03_Step3_对接环境/
│   ├── 04_Step4_批量对接/
│   ├── 05_Step5_评估/
│   ├── 06_Step6_数据库解析/
│   ├── 07_Step7_结构处理/
│   ├── 08_Step8_特征生成/
│   └── utils/
│       └── run_experiment.py            # 通用实验运行器
│
├── configs/
│   └── base.yaml                        # 公共默认参数
│
├── logs/                                # 运行日志（stdout/stderr）
├── results/                             # 汇总结果（跨实验对比表、图表）
├── reports/                             # 正式报告
└── PathB_计划与进度.md                   # 本文件
```

---

## 五、进度记录

| 日期 | Step | 内容 | 状态 |
|------|:----:|------|:----:|
| 2026-02-12 | - | PathB 规划完成，目录结构创建 | ✅ |
| 2026-02-12 | Step 1 | **数据准备与代码改造** | **✅ 已完成** |
| | | - B6 数据集迁移到 00_shared/datasets/B6_v1/ (516条) | ✅ |
| | | - 共享特征迁移到 00_shared/features/ (~19GB，仅 Windows) | ✅ |
| | | - extract_pocket_ligand.py 改造 (Heme+半径参数化) | ✅ |
| | | - run_experiment.py 实验运行器 + base.yaml 配置 | ✅ |
| 2026-02-20 | Step 2 | **代码准备：HETATM 问题验证** | ✅ |
| | | - Codex 四轮深度验证确认 HETATM 阻断风险 | ✅ |
| | | - 确认修复方案：PathB fork 版 PDBProtein 支持 ATOM+HETATM | ✅ |
| | | - 确认实验定位：推理敏感性实验（非生化机理验证） | ✅ |
| 2026-02-21 | Step 2 | **代码准备：编写 Step 8.2/8.3 脚本** | ✅ |
| | | - step8_align_ligand.py（配体原子对齐，CLI 可配置） | ✅ |
| | | - step8_generate_structure_lmdb.py（LMDB 生成，含 HETATM 修复） | ✅ |
| | | - run_experiment.py 集成 8.2/8.3 | ✅ |
| 2026-02-21 | Step 2 | **代码审核：Codex 三轮审核** | ✅ |
| | | - 修复 7 个 bug（2 HIGH + 3 MEDIUM + 2 LOW） | ✅ |
| | | - 详见 sessions/02_Step2_.../session_log.md | ✅ |
| 2026-02-22 | Step 2 | **环境准备：MacBook 依赖安装** | ✅ |
| | | - 安装 RDKit 2025.09.5 + lmdb 1.7.5 + PyG 2.7.0 | ✅ |
| | | - dry-run 验证通过（4 个实验命令正确） | ✅ |
| 2026-02-22 | Step 2 | **阻断：缺失数据文件** | **🔴 阻断中** |
| | | - Mac 上缺少：PDB 原始文件 + 4 个特征 LMDB | 🔴 |
| | | - 需从 Windows 迁移或在 Windows 上运行 | 🔴 |
| | Step 2 | **待做：运行 Steps 8.1-8.3** | ⏳ 等数据迁移 |
| | Step 2 | **待做：编写 Step 9/10 脚本** | ⏳ |
| 2026-03-02 | Step 3 | **对接环境搭建** | **✅ 已完成** |
| | | - 6 lib 模块 + 5 CLI 脚本，全部经 Codex 审核 | ✅ |
| | | - 受体 282/292 成功，配体 250/252 成功 | ✅ |
| | | - Pilot 5: 4/5 成功; Pilot 50: 48/50 = 96% 成功 | ✅ |
| | | - 下游 LMDB 兼容性: 48/48 = 100% PASS | ✅ |
| 2026-03-03 | Step 4 | **随机负样本生成与全量对接** | **✅ 已完成** |
| | | - Manifest: 270 正 + 2,527 负 = 2,797 对 | ✅ |
| | | - 首次运行: 2,530/2,797 = 90.5% 成功 (4.9h, 12 workers) | ✅ |
| | | - 三层重试: Tier0(FE阈值) + Tier1(延时) + Tier2(降exhaustiveness) | ✅ |
| | | - 最终: 2,766/2,797 = **98.9%** 成功 (265正 + 2,501负, 比1:9.44) | ✅ |
| | | - 产出: results/04_Step4_批量对接/data.csv (模型兼容) | ✅ |
| 2026-03-03 | Step 5 | **重构评估 + Gate B 决策** | **✅ 已完成** |
| | | - 特征复用: ESM/Reaction/GROVER/Morgan 全覆盖 (无需重新生成) | ✅ |
| | | - 配体对齐: 2,766/2,766 = 100% (8s) | ✅ |
| | | - 结构 LMDB: 2,766/2,766 = 100% (23s, HETATM=0) | ✅ |
| | | - 模型推理: 2,766/2,766 = 100% | ✅ |
| | | - **Gate B: INFORMATIVE FAIL** (AUC-ROC=0.5170 [0.4804, 0.5521]) | ✅ |
| | | - 根因: Dockability≠Catalysis + OOD enzymes (0% overlap) | ✅ |
| | | - Codex 2 轮审核 + bootstrap CI + 混杂变量分析 | ✅ |
| | | - 产出: analysis/ (报告+ROC+分布图+per-enzyme), predictions.csv | ✅ |
| | Step 6 | 外部数据库收集与清洗 | ⏳ 待开始 |
| | Step 7 | 结构获取与质量控制 | ⏳ 待开始 |
| | Step 8 | 全量特征生成与数据集交付 | ⏳ 待开始 |

### 缺失文件清单（待迁移）

以下文件被 `.gitignore` 排除，仅在 Windows 机器上存在：

| 文件 | 大小 | 原始位置 | 目标位置 | Step 需要 |
|------|------|---------|---------|----------|
| 627 个 PDB 对接文件 | ~数 GB | `PathA/data/01_Step1_PDB文件/` | 同 | 8.1 |
| enzyme_features.lmdb | 4.29 GB | `PathA/data/05_Step5_.../` | `PathB/data/00_shared/features/` | 9 |
| reaction_features.lmdb | 4.29 GB | `PathA/data/06_Step6_.../` | `PathB/data/00_shared/features/` | 9 |
| grover_fingerprint.lmdb | 10.74 GB | `PathA/data/07_Step7_.../` | `PathB/data/00_shared/features/` | 9 |
| morgan_fingerprint.npy | 447 KB | `PathA/data/07_Step7_.../` | `PathB/data/00_shared/features/` | 9 |
