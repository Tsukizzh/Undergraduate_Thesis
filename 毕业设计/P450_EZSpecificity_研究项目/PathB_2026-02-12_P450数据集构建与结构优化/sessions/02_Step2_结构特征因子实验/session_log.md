# Step 2: 结构特征 2×2 因子实验 — Session Log

## 状态：✅ 已完成（Gate A PASS: 10Å/noHeme, AUC-ROC=0.7115）

## 时间线

| 日期 | 阶段 | 内容 |
|------|------|------|
| 2026-02-20 | 代码准备 | HETATM 问题 Codex 四轮深度验证 |
| 2026-02-21 | 代码准备 | 编写 step8_align_ligand.py、step8_generate_structure_lmdb.py |
| 2026-02-21 | 代码审核 | Codex 三轮代码审核，修复 7 个 bug |
| 2026-02-22 | 代码审核 | run_experiment.py 集成 8.2/8.3，Codex 验证 |
| 2026-02-22 | 环境准备 | MacBook 安装依赖、dry-run 验证通过 |
| 2026-02-22 | 阻断 | 发现 Mac 上缺少大文件（PDB + 特征 LMDB），无法运行 |
| 2026-02-25 | Windows 执行 | 创建 pathb-step2 分支，4 个实验 8.1-8.3 全部成功 |
| 2026-02-25 | 推理 | 编写 step9_inference.py，4 个实验推理完成 |
| 2026-02-25 | 分析 | 编写 step10_comparative_analysis.py，生成 ROC/heatmap/Gate A |
| 2026-02-25 | 审核 | Codex 审核通过，确认 Heme OOD 是真实信号 |
| 2026-02-26 | 文档 | 更新所有进度文件，提交 docs commit |

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

## 四、Windows 执行阶段（2026-02-25）

由于 Mac 上缺少 22GB 的数据文件（PDB + 特征 LMDB），决定直接在 Windows 机器上执行。

### 4.0 前置操作：git pull + 分支创建

**拉取 Mac 上的代码**：
- 用户从 Mac push 了 2 个提交到 GitHub（`ce003e1` 和 `b3f6899`）
- 在 Windows 上 pull 回来（SSH 环境下需用 Device code 认证方式：`github.com/login/device`）
- 添加 `.vscode/` 到 `.gitignore`（SSH 连接 VSCode 自动生成的配置文件）

**创建实验分支**：
```bash
git checkout -b pathb-step2
```
此分支从 main 分出，用于隔离 Step 2 实验工作。

### 4.1 Windows 环境验证

使用 `D:/anaconda3/envs/torch/python.exe`（因 conda 环境无法通过 `conda run` 在 SSH 环境下工作）。

| 组件 | 版本 | 状态 |
|------|------|------|
| Python | 3.10 (torch env) | ✅ |
| PyTorch | 2.3.0+cu121 | ✅ |
| CUDA | 12.1 | ✅ |
| GPU | NVIDIA (≥8GB) | ✅ |
| RDKit | 2025.03.6 | ✅ |
| lmdb | 0.9.33 | ✅ |
| PyTorch Geometric | 2.7.0 | ✅ |
| BioPython | OK | ✅ |

**注意**：`conda run -n torch python -c "多行脚本"` 在 Windows SSH 下会报 `AssertionError: Support for scripts where arguments contain newlines not implemented`。解决方案：直接用 Python 解释器的绝对路径。

### 4.2 必需文件检查

| 文件 | 状态 | 说明 |
|------|------|------|
| PDB 文件 (627 expected) | ⚠️ 619 found | 8 个缺失，非阻断 |
| enzyme_features.lmdb | ✅ | `data/00_shared/features/` |
| reaction_features.lmdb | ✅ | 同上 |
| grover_fingerprint.lmdb | ✅ | 同上 |
| morgan_fingerprint.npy | ✅ | 同上 |
| B6 datasets (data.csv etc.) | ✅ | `data/00_shared/datasets/B6_v1/` |
| Model checkpoint | ✅ | `saved_model/model/run_0/models/best-checkpoint.ckpt` (22.3MB) |
| Model config YAML | ✅ | `saved_model/model/run_0/complete-full-random-all-0-complex.yml` |

**模型 checkpoint 小插曲**：初次检查时 `saved_model/model/run_0/` 目录为空（可能被误删或未恢复）。用户手动从备份恢复了 `best-checkpoint.ckpt` 及 4 个备份版本。

### 4.3 运行 4 个因子实验 Steps 8.1-8.3

**执行命令**：
```bash
python scripts/utils/run_experiment.py --run_all
```

**执行顺序**（智能排序：radius DESC + heme DESC，确保最宽松条件先跑以生成完整的共享对齐数据）：
1. EXP02_B6_10A_Heme（10Å + Heme）← 第一个跑，产生 shared_alignment
2. EXP01_B6_10A_noHeme（10Å + 无 Heme）← 复用 shared_alignment
3. EXP04_B6_6A_Heme（6Å + Heme）← 复用
4. EXP03_B6_6A_noHeme（6Å + 无 Heme）← 复用

**各步骤耗时**：

| 实验 | Step 8.1 (提取) | Step 8.2 (对齐) | Step 8.3 (LMDB) | 总计 |
|------|----------------|----------------|-----------------|------|
| EXP02 (首个) | 80s | 2.2s | 11.1s | 93.3s |
| EXP01 | 77.5s | cached | 10.9s | 88.4s |
| EXP04 | 75.9s | cached | 10.2s | 86.1s |
| EXP03 | 75.2s | cached | 9.7s | 84.9s |
| **合计** | | | | **~6 分钟** |

**共享对齐缓存机制验证**：
- Step 8.2（对齐）只在 EXP02 时实际运行（2.2s），生成 `shared_alignment/alignment_summary.csv`（539 条记录）
- 后续 3 个实验检测到 alignment_summary.csv 存在且行数 > 0，直接跳过
- 这验证了 run_experiment.py 的缓存逻辑正常工作

**每个实验的产出**：
```
EXP*/
├── config.yaml                    # 实验参数快照（含完整路径信息）
├── high_quality_id.txt            # 有效 Dock Index 列表（均为 495 条）
├── structure_build_summary.csv    # LMDB 构建详情（每条记录的成功/失败）
├── structure_features.lmdb        # 模型可读的结构特征（二进制）
├── structure_features.lmdb-lock   # LMDB 锁文件
└── structure_features/
    ├── pocket/                    # 提取的口袋 PDB 文件
    ├── raw_ligand/                # 提取的原始配体 SDF 文件
    ├── extraction_summary.csv     # 提取统计（540 条）
    └── extraction_failures.csv    # 提取失败记录（7-8 条）
```

**关键数据量**：
- 输入：B6 数据集 516 条记录 → mapping 后 540 个 dock index → 619 个 PDB 文件可用
- 经过 Step 8.1 提取 + Step 8.2 对齐 + Step 8.3 质量过滤
- 输出：**每个实验均为 495 条有效记录**（high_quality_id.txt）
- 过滤掉的 21 条：invalid structures（PDB 缺失或提取失败）+ 对齐失败

---

## 五、推理阶段 — Step 9（2026-02-25）

### 5.1 脚本编写

**文件**: `scripts/02_Step2_因子实验/step9_inference.py`

**编写过程**：
1. 读取 PathA 的 `step9_inference.py`（322 行）作为参考
2. 向 Codex 请求 PathB 版代码原型（要求 CLI 可配置，read-only sandbox）
3. Codex 返回 281 行代码原型
4. 我对 Codex 原型进行审核和重写，最终版本 ~160 行（去除冗余，精简结构）

**与 PathA 版的关键区别**：

| 方面 | PathA | PathB |
|------|-------|-------|
| 路径配置 | 全部硬编码在文件顶部 | CLI 参数（`--experiment_dir`, `--shared_features` 等） |
| 结构 LMDB | 固定路径 | 从 `--experiment_dir` 下自动发现 |
| 共享特征 | 指向 PathA data 子目录 | 指向 `--shared_features` 目录 |
| 数据集 CSV | PathA Step 4 的 data.csv | `--shared_datasets` 下的 B6 data.csv |
| 输出位置 | PathA data/09_Step9_模型推理/ | 直接写入 experiment_dir/predictions.csv |
| 复用性 | 单次运行 | 可循环调用（每个实验一次） |

**核心推理逻辑（与 PathA 完全一致）**：
```python
config = load_config(str(config_yaml))
config.num_cpus = 0          # Windows: num_workers=0
config.data.full_data = False # 不需要 str_features.lmdb
config.data.representer = "structure_sequence"

dm = Singledataset(config)
test_loader = dm.test_dataloader()
model = SS.load_from_checkpoint(checkpoint, config=config, map_location="cpu")
model.eval(); model.to(device)

for batch in test_loader:
    batch = batch.to(device)
    logits, _ = model(batch)  # [B, 1]
```

**`full_data=False` 的关键性**：这个标志让模型跳过加载 `str_features.lmdb`（序列级结构特征），该文件由 `seq_process()` 从 SMILES 生成，与实验参数无关但生成非常耗时。PathA 在 Step 9 中也使用了此标志。

### 5.2 执行

**执行命令**（4 个实验逐个运行，避免 GPU OOM）：
```bash
for exp in EXP01 EXP02 EXP03 EXP04; do
    python step9_inference.py \
        --experiment_dir  data/02_Step2_因子实验/$exp \
        --shared_features data/00_shared/features \
        --shared_datasets data/00_shared/datasets \
        --checkpoint_dir  saved_model/model/run_0
done
```

**运行日志摘录**：
```
[EXP01_B6_10A_noHeme] Device: cuda:0
[EXP01_B6_10A_noHeme] Running inference...
Invalid reaction data: 0
#Invalid structures: 21
516 → 495 → 495
  Processed 160/495 samples...
  Processed 320/495 samples...
  Processed 480/495 samples...
[EXP01_B6_10A_noHeme] Saved: .../EXP01_B6_10A_noHeme/predictions.csv
[EXP01_B6_10A_noHeme] Samples: 495, prob range=(0.0000, 0.9280)
```

**4 个实验均成功**，每个处理 495 个样本。

**数据流解释**：
```
B6 data.csv: 516 条记录
    ↓ Singledataset 加载
test_prediction_df: 516 条（原始）
    ↓ high_quality_id.txt 过滤（495 条有效 Dock Index）
    ↓ invalid structure 排除（21 条无效结构）
最终推理: 495 条
```

**输出格式**：每个实验的 `predictions.csv` 包含 7 列：
```csv
Dock Index, Enzyme Index, Substrate Index, Label, score, logit, prob
```
- `logit`：模型原始输出值（可正可负）
- `prob`：经过 sigmoid 转换后的概率（0~1）
- `score`：与 logit 相同（保持与 PathA 输出格式兼容）

### 5.3 推理结果速览

| 实验 | prob 范围 | 说明 |
|------|-----------|------|
| EXP01 (10Å/noHeme) | (0.0000, 0.9280) | 区分度好 |
| EXP02 (10Å/Heme) | (0.0000, 1.0000) | 极端概率出现 |
| EXP03 (6Å/noHeme) | (0.0000, 0.9855) | 区分度好 |
| EXP04 (6Å/Heme) | (0.0000, 1.0000) | 极端概率出现 |

Heme 实验（EXP02/04）出现了 prob=1.0（logit 极大值），说明 Heme 原子扰动了模型的输出分布。

---

## 六、分析阶段 — Step 10（2026-02-25）

### 6.1 脚本编写

**文件**: `scripts/02_Step2_因子实验/step10_comparative_analysis.py`

**编写过程**：
1. 读取 PathA 的 `step10_analysis.py`（693 行，单实验深度分析）
2. 向 Codex 请求 PathB 版代码原型（侧重 4 实验对比，而非单实验深度）
3. Codex 返回 422 行代码原型
4. 重写为 ~200 行精简版本

**功能**：
- 自动发现所有 `EXP*` 目录并读取 `predictions.csv`
- 从实验名称解析因子（例如 `EXP02_B6_10A_Heme` → heme=True, radius=10.0）
- 计算每个实验的 AUC-ROC、AUC-PR、accuracy
- 绘制 4 条 ROC 曲线叠加图
- 绘制 2×2 因子热力图
- 计算主效应（Heme effect、Radius effect）
- 生成 Gate A 决策报告

**Gate A 决策逻辑**：
```
baseline = PathA AUC-ROC = 0.6636
threshold = 0.005 (0.5%)
if best_experiment.auc_roc - baseline > threshold:
    Gate A = PASS → 采用最佳配置
else:
    Gate A = HOLD → 保持基线
```

### 6.2 执行

```bash
python step10_comparative_analysis.py \
    --experiments_dir data/02_Step2_因子实验 \
    --shared_datasets data/00_shared/datasets \
    --output_dir      data/02_Step2_因子实验/analysis
```

### 6.3 实验结果

#### 核心指标

| 实验 | Heme | Radius | AUC-ROC | AUC-PR | Accuracy | N |
|------|------|--------|---------|--------|----------|---|
| **EXP01** | No | 10Å | **0.7115** | 0.7356 | 0.5091 | 495 |
| EXP03 | No | 6Å | 0.6678 | 0.6990 | 0.5131 | 495 |
| EXP02 | Yes | 10Å | 0.4894 | 0.5492 | 0.4848 | 495 |
| EXP04 | Yes | 6Å | 0.4257 | 0.4973 | 0.4525 | 495 |

#### 排名
```
EXP01 (0.7115) > EXP03 (0.6678) > EXP02 (0.4894) > EXP04 (0.4257)
```

#### 主效应分析

| 效应 | 计算方式 | 值 | 含义 |
|------|---------|-----|------|
| **Heme 效应** | mean(Heme实验) - mean(noHeme实验) | **-0.2322** | Heme 严重损害性能 |
| **Radius 效应** | mean(10Å实验) - mean(6Å实验) | **+0.0537** | 更大口袋有帮助 |

#### 关于 accuracy ~50% 的说明

accuracy 在 ~50% 是正常的，不代表模型差。原因：
- B6 数据集接近平衡（272 pos + 244 neg）
- accuracy 使用固定阈值 0.5，对 logit 分布偏移敏感
- **AUC-ROC 不受阈值影响**，是更可靠的排序能力指标
- PathA 结果也是类似的 accuracy（~50%）配合较好的 AUC-ROC（0.6636）

#### Gate A 决策

```
Best experiment:     EXP01_B6_10A_noHeme
AUC-ROC:            0.7115
PathA baseline:     0.6636
Delta:              +0.0479 (+4.79%)
Threshold:          0.005

Gate A Decision:    PASS
Recommendation:     采用 10Å/noHeme 配置用于后续 PathB 步骤
```

---

## 七、Codex 审核（2026-02-25）

Step 9/10 执行完成后，将结果提交 Codex 审核。

### 审核问题

1. Heme 的 -23% AUC 下降是否可疑？是 bug 还是真实信号？
2. EXP01(0.7115) vs PathA(0.6636) 的 +4.79% 差异如何解释？
3. accuracy ~50% 是否正常？
4. Gate A PASS 决策是否正确？

### Codex 审核结论

**1. Heme 惩罚是真实的 OOD 效应（非 bug）**

Codex 检查了 3 个关键代码路径：
- PathA `extract_pocket_ligand.py` → 明确排除 HETATM
- PathB `extract_pocket_ligand.py` → 正确注入 Heme HETATM 原子
- PathB `step8_generate_structure_lmdb.py` → 正确解析 HETATM 行

Heme 造成性能下降的根本原因：
- **Fe 原子不在模型元素词汇表中**（`src/Datasets/Structure/transforms.py` 的 `atomic_numbers` 列表不含 Fe=26）→ Fe 得到全零特征向量
- **HEM 残基映射为 UNK**（residue index=20）→ 与模型训练时从未见过的模式
- 模型从未在训练数据中见过 HETATM 记录 → 完全 out-of-distribution
- 这些 OOD 特征注入到 GNN 的消息传递中，**污染了**整个口袋的表示

**结论**：这不是 bug，而是模型的固有 OOD 限制。如果要让 Heme 信息真正有用，需要在 Path C 重训时扩展元素词汇表并纳入 HETATM 训练数据。

**2. EXP01 vs PathA 的差异不可直接比较**

Codex 发现两者使用了**不同的数据集**：
- PathA：Step 4 的 data.csv → 517 个可用样本
- PathB：B6 的 data.csv → 495 个可用样本

这两个数据集有以下差异：
- 记录数不同（PathA 540 行 vs B6 516 行）
- B6 排除了产物（PRODUCT 类配体），PathA 包含
- 样本构成不同 → AUC-ROC 不可直接做差

**如要严格对比**，需在完全相同的 495 条 ID 上重新跑 PathA 流程。但这不影响 Step 2 的内部比较（4 个实验都用 B6 数据集）。

**3. accuracy ~50% 正常**（同 PathA 模式）

**4. Gate A PASS 正确**（在 4 个实验的内部比较中，EXP01 明确最优）

### Codex 建议的后续步骤

1. 在相同 495 个 ID 上重算 PathA 基线，才能宣称 "+4.79% 改进"
2. 可尝试 Heme 消融：(a) Heme 去掉 Fe 原子，(b) 只保留口袋半径内的 Heme 原子
3. 将 Heme 结果记录为"模型 OOD 限制"，而非"生化结论"

---

## 八、产出文件完整清单

### 脚本文件
```
scripts/02_Step2_因子实验/
├── extract_pocket_ligand.py         (Step 1 产出，Step 8.1 用)
├── step8_align_ligand.py            (Mac 阶段编写，Step 8.2 用)
├── step8_generate_structure_lmdb.py (Mac 阶段编写，Step 8.3 用，含 HETATM 修复)
├── step9_inference.py               (Windows 阶段编写，CLI 可配置推理)
└── step10_comparative_analysis.py   (Windows 阶段编写，4 实验对比 + Gate A)

scripts/utils/
└── run_experiment.py                (Mac 阶段更新，集成 8.1→8.2→8.3)
```

### 数据文件
```
data/02_Step2_因子实验/
├── EXP01_B6_10A_noHeme/
│   ├── config.yaml                  # 实验参数快照
│   ├── high_quality_id.txt          # 495 条有效 Dock Index
│   ├── predictions.csv              # 推理结果（7 列 × 495 行）
│   ├── structure_build_summary.csv  # LMDB 构建详情
│   ├── structure_features.lmdb      # 结构特征（gitignored）
│   └── structure_features/
│       ├── pocket/                  # 口袋 PDB（gitignored）
│       ├── raw_ligand/              # 原始配体 SDF（gitignored）
│       ├── extraction_summary.csv   # 提取统计
│       └── extraction_failures.csv  # 提取失败记录
├── EXP02_B6_10A_Heme/              (同结构)
├── EXP03_B6_6A_noHeme/             (同结构)
├── EXP04_B6_6A_Heme/               (同结构)
├── shared_alignment/                # 共享对齐数据（gitignored）
│   ├── ligand/                      # 对齐后配体 SDF
│   └── alignment_summary.csv        # 539 条对齐记录
└── analysis/
    ├── comparative_metrics.csv      # 4 行指标汇总
    ├── gate_a_recommendation.txt    # Gate A 决策报告
    └── figures/
        ├── comparative_roc.png      # 4 条 ROC 曲线叠加图
        └── factorial_heatmap.png    # 2×2 AUC-ROC 热力图
```

### Git 提交
```
pathb-step2 分支:
├── ba7c22a  PathB Step 2: execute 2x2 factorial experiments + inference + analysis
└── 48bdbcb  docs: update all progress trackers for PathB Step 2 completion
```

---

## 九、经验总结

1. **Windows SSH 下 conda 的坑**：`conda run` 不支持多行参数，需直接用 Python 解释器绝对路径
2. **Git 凭证在 SSH 下的坑**：`wincredman` 凭证存储在 SSH 会话下不可用，需用 Device code 认证
3. **`full_data=False` 是推理的关键**：跳过 `str_features.lmdb`（序列级特征），大幅简化依赖
4. **共享对齐缓存设计正确**：4 个实验只跑了一次对齐（2.2s），节省 ~6s
5. **Heme OOD 是预期内的**：Phase 2.1 的 Codex 分析已预警 Fe 不在词汇表中，实验证实了这一预判
6. **不同数据集的 AUC-ROC 不可直接比较**：这是 Codex 审核发现的重要提醒
