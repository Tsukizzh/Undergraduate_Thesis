# Step 2: 结构特征 2×2 因子实验 — Session Log

## 状态：代码准备完毕，等待数据文件迁移后运行

## 时间线

| 日期 | 阶段 | 内容 |
|------|------|------|
| 2026-02-20 | 代码准备 | HETATM 问题 Codex 四轮深度验证 |
| 2026-02-21 | 代码准备 | 编写 step8_align_ligand.py、step8_generate_structure_lmdb.py |
| 2026-02-21 | 代码审核 | Codex 三轮代码审核，修复 7 个 bug |
| 2026-02-22 | 代码审核 | run_experiment.py 集成 8.2/8.3，Codex 验证 |
| 2026-02-22 | 环境准备 | MacBook 安装依赖、dry-run 验证通过 |
| 2026-02-22 | **阻断** | 发现 Mac 上缺少大文件（PDB + 特征 LMDB），无法运行 |

---

## 一、Step 2 的目标

运行 2×2 因子实验，测试两个结构变量对模型预测的影响：

| 实验 | 口袋半径 | Heme | 目的 |
|:---:|:---:|:---:|------|
| EXP01 | 10Å | 无 | 基线（与 PathA 条件一致） |
| EXP02 | 10Å | 有 | 单独测试 Heme 效果 |
| EXP03 | 6Å | 无 | 单独测试口袋缩小效果 |
| EXP04 | 6Å | 有 | 测试组合效果 |

每个实验的执行流程：
```
Step 8.1: extract_pocket_ligand.py   → 从 PDB 提取口袋 + 配体
Step 8.2: step8_align_ligand.py      → 对齐配体原子编号（4 个实验共享）
Step 8.3: step8_generate_structure_lmdb.py → 生成 structure LMDB
Step 9:   模型推理                    → 用模型预测（TODO：脚本未写）
Step 10:  结果分析                    → 对比 AUC-ROC 等（TODO：脚本未写）
```

---

## 二、代码准备阶段（已完成）

### 2.1 HETATM 阻断问题发现与验证（2026-02-20）

**问题**：Step 1 编写的 `extract_pocket_ligand.py` 支持 `--include_heme` 参数，能将 Heme 原子写入 pocket PDB。但下游 `PDBProtein._parse()` 只读 ATOM 行，**会静默丢弃 Heme 的 HETATM 行**。如果不修复，Heme 因子无效，2×2 实验退化为只测半径。

**Codex 四轮深度验证确认**：
1. BioPython PDBIO 根据残基的 `hetfield` 决定写 ATOM 还是 HETATM
2. Heme 残基的 `hetfield != " "`，因此写出为 HETATM 行
3. `PDBProtein._parse()` 只匹配 `line[0:6].strip() == 'ATOM'`，完全跳过 HETATM
4. 结果：Heme 原子被静默丢弃

**额外技术发现**：
- PDB 格式中 ATOM 与 HETATM 列布局完全一致（wwPDB 标准）
- PathA 的 539 个 pocket PDB 中 HETATM 行数为 0（全量扫描已验证）
- Heme 原子量级：单个 HEM 残基约 43 个非氢重原子，占 pocket 平均原子数的 ~10.5%
- Fe(26) 不在 `FeaturizeProteinAtom.atomic_numbers` 中（会得到全零向量）
- Heme 中的 C/N/O 可正确编码；HEM 残基会被标记为 UNK（index=20）
- residue key 冲突风险极低（实际数据中 ATOM 与 HETATM 的 chain+resSeq+iCode 无冲突）

**实验科学定位**：
- 修复后本实验为「推理敏感性/结构扰动实验」，**不是**「Heme 生化机理验证」
- 模型训练时未见 HETATM，且 Fe 特征缺失
- 若要机理层面结论，需在 Path C 重训时纳入 Heme 和 Fe 特征通道

### 2.2 新脚本编写（2026-02-21）

#### step8_align_ligand.py
- **位置**: `scripts/02_Step2_因子实验/step8_align_ligand.py`
- **作用**: 将 PDB 提取的配体与 SMILES 原子编号对齐
- **与 PathA 差异**: CLI 可配置路径（PathA 为硬编码）；核心对齐逻辑完全一致
- **关键设计**: 对齐结果与 pocket 参数无关，4 个实验共享同一份对齐数据

#### step8_generate_structure_lmdb.py
- **位置**: `scripts/02_Step2_因子实验/step8_generate_structure_lmdb.py`
- **作用**: 将 pocket PDB + aligned ligand SDF 打包成模型能读的 LMDB 格式
- **核心改动**: `PDBProtein._parse()` 支持 ATOM + HETATM（修复 Heme 被丢弃的问题）
  - `line[0:6].strip() == 'ATOM'` → `line[0:6].strip() in ('ATOM', 'HETATM')`
  - `is_backbone` 仅对 ATOM 记录为 True（HETATM 不算骨架）
  - 增加 record_type 字段到 atom dict（防止 residue key 冲突）
  - 增加 HETATM 统计日志（n_atom, n_hetatm 等）

#### run_experiment.py 更新
- 集成 Step 8.2 和 8.3 调用
- 共享对齐：Step 8.2 的 alignment 是 4 个实验共享的，只需跑一次
- 智能排序：先跑 10Å+Heme（最宽松条件），确保共享对齐数据最完整
- 缓存验证：alignment_summary.csv 存在且行数 > 0 才视为有效

### 2.3 Codex 三轮代码审核（2026-02-21）

**审核范围**: step8_generate_structure_lmdb.py、step8_align_ligand.py、run_experiment.py

**已修复的 7 个问题**:

| # | 严重度 | 文件 | 问题 | 修复 |
|---|--------|------|------|------|
| 1 | MEDIUM | step8_generate_structure_lmdb.py | `element_symb` 回退 `line[13:14]` 对 Fe 等双字母元素只取 1 字符 | → `line[12:14].strip().capitalize()` |
| 2 | MEDIUM | step8_generate_structure_lmdb.py | `r["success"] == "True"` 大小写敏感 | → `.strip().lower() == "true"` |
| 3 | LOW | step8_generate_structure_lmdb.py | AtomMapNum 全零检查误拒单原子配体 | → 加 `mol.GetNumAtoms() > 1` 守卫 |
| 4 | HIGH | run_experiment.py | 共享对齐实验顺序依赖（小半径先跑→对齐不完整） | → 排序 radius DESC + heme DESC |
| 5 | HIGH | run_experiment.py | `alignment_summary.exists()` 不检测空/损坏文件 | → 行数验证 + 重跑逻辑 |
| 6 | MEDIUM | step8_align_ligand.py | 大量 raw_ligand 缺失时仍静默成功 | → 成功率 < 50% 返回错误码 1 |
| 7 | LOW | step8_generate_structure_lmdb.py | 空 target_docks 仍生成空 LMDB | → fail fast 返回 2 |

**确认无问题的项目**:
- `is_backbone` 只对 ATOM 记录为 True
- PDB 列布局 ATOM/HETATM 一致
- pickle 兼容性保持（StructureComplexData.__module__ 设为 'Datasets.Structure.utils'）
- `project_root` 路径解析正确

**接受的残余风险**:
- alignment_summary 验证用 `row_count > 0`（完美验证需循环依赖，当前已足够实用）

---

## 三、环境准备（已完成）

### 3.1 MacBook 环境信息
- **机器**: macOS Darwin 25.3.0 (ARM)
- **conda 环境**: `torch`（`/opt/anaconda3/envs/torch`）
- **Python**: 3.12.4
- **PyTorch**: 2.5.1（MPS，无 CUDA/GPU）
- **无 GPU**: Steps 8.1-8.3 不需要 GPU（纯 CPU 数据处理）

### 3.2 依赖安装（2026-02-22）

| 包 | 版本 | 安装方式 |
|---|---|---|
| RDKit | 2025.09.5 | conda-forge |
| lmdb | 1.7.5 | pip |
| PyTorch Geometric | 2.7.0 | pip |
| BioPython | 1.85 | 已有 |
| PyYAML | 6.0.1 | 已有 |

### 3.3 dry-run 验证（2026-02-22）

```
python run_experiment.py --run_all --dry_run
```

**结果**: 4 个实验全部通过。排序正确：
1. EXP02_B6_10A_Heme（10Å + Heme）← 最先跑
2. EXP01_B6_10A_noHeme（10Å + 无 Heme）
3. EXP04_B6_6A_Heme（6Å + Heme）
4. EXP03_B6_6A_noHeme（6Å + 无 Heme）← 最后跑

---

## 四、当前阻断：缺失数据文件

dry-run 通过后，preflight 检查发现以下文件在这台 Mac 上不存在（被 `.gitignore` 排除，仅在原始 Windows 机器上）：

### 缺失文件清单

| 文件 | 大小 | 用途 | 哪一步需要 | 来源 |
|------|------|------|-----------|------|
| `PathA/data/01_Step1_PDB文件/` | ~数 GB | 原始 PDB 对接文件（627个） | Step 8.1 | PathA Step 1 下载 |
| `data/00_shared/features/enzyme_features.lmdb` | 4.29 GB | 酶 ESM 嵌入 | Step 9 | PathA Step 5 |
| `data/00_shared/features/reaction_features.lmdb` | 4.29 GB | 底物反应图 | Step 9 | PathA Step 6 |
| `data/00_shared/features/grover_fingerprint.lmdb` | 10.74 GB | GROVER 分子指纹 | Step 9 | PathA Step 7 |
| `data/00_shared/features/morgan_fingerprint.npy` | 447 KB | Morgan 指纹 | Step 9 | PathA Step 7 |

### 影响分析

| 步骤 | 能否在 Mac 上运行 | 原因 |
|------|:---:|------|
| Step 8.1（提取口袋/配体） | 不能 | 需要原始 PDB 文件 |
| Step 8.2（对齐配体） | 不能 | 依赖 8.1 产出的 raw_ligand |
| Step 8.3（生成 LMDB） | 不能 | 依赖 8.1/8.2 产出 |
| Step 9（推理） | 不能 | 需要特征 LMDB + GPU |
| Step 10（分析） | 不能 | 依赖 Step 9 产出 |

### Mac 上已有的 PathA 中间产物（可用于局部测试）

| 文件 | 数量 | 说明 |
|------|------|------|
| `PathA/data/08_Step8_.../str_tmp_data/pocket/*.pdb` | 539 个 | PathA 10Å 无 Heme 条件下的口袋 |
| `PathA/data/08_Step8_.../str_tmp_data/ligand/*.sdf` | 517 个 | PathA 已对齐的配体 |
| `PathA/data/08_Step8_.../str_tmp_data/raw_ligand/*.sdf` | 532 个 | PathA 原始配体 |

这些可用于测试 Step 8.3 的 LMDB 生成逻辑，但**不能替代完整的因子实验**（因为条件相同，无法测试 Heme/6Å 的差异）。

---

## 五、待完成事项

### 高优先级
1. **迁移缺失文件**：将 PDB 原始文件和特征 LMDB 从 Windows 机器复制到 Mac（或改在 Windows 上运行）
2. **运行 Steps 8.1-8.3**：4 个实验的结构特征生成
3. **编写 Step 9 脚本**：PathB 版推理脚本（可配置特征路径）
4. **编写 Step 10 脚本**：PathB 版分析脚本（4 组对比）

### 文件迁移后的执行命令

```bash
# 在有完整数据的机器上
conda activate torch  # 或 ezspecificity
cd PathB_2026-02-12_P450数据集构建与结构优化/

# 先跑单个实验验证
python scripts/utils/run_experiment.py --name EXP02_B6_10A_Heme --pocket_radius 10.0 --include_heme

# 验证通过后跑全部 4 个
python scripts/utils/run_experiment.py --run_all
```

---

## 六、产出文件清单（Step 2 代码准备阶段）

```
scripts/02_Step2_因子实验/
├── extract_pocket_ligand.py      (Step 1 产出，Step 8.1 用)
├── step8_align_ligand.py         (本阶段新增，Step 8.2 用)
└── step8_generate_structure_lmdb.py (本阶段新增，Step 8.3 用，含 HETATM 修复)

scripts/utils/
└── run_experiment.py             (本阶段更新，集成 8.2+8.3)

data/02_Step2_因子实验/
├── EXP01_B6_10A_noHeme/config.yaml  (dry-run 自动生成)
├── EXP02_B6_10A_Heme/config.yaml
├── EXP03_B6_6A_noHeme/config.yaml
└── EXP04_B6_6A_Heme/config.yaml
```
