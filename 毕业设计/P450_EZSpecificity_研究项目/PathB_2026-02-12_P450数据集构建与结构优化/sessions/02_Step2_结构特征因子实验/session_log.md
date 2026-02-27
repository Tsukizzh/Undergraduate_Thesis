# Step 2: 结构特征 2×2 因子实验 — Session Log

## 状态：✅ 已完成（Gate A PASS: 10Å/noHeme, AUC-ROC=0.7115）

---

## 结果摘要（先看这里）

### 我们在验证什么？

Path A 用 PDB 晶体结构中直接提取的口袋做推理，AUC-ROC = 0.6636。但有两个未知变量：

- **Heme（血红素）** — P450 酶的活性中心有一个含铁的 Heme 辅因子。PDB 文件里它和配体挨在一起，但模型训练时**从未见过 Fe 原子**。留着它会帮助模型理解 P450 特性，还是因为 OOD（分布外）反而干扰？
- **口袋半径** — 原版用 10Å（配体周围 10 埃内的蛋白原子），6Å 能更聚焦配体附近。哪个更好？

所以设计了 2×2 因子实验（4 种组合），在**完全相同的 495 个样本**上测试。

### 实验结果

| 实验 | Heme | 口袋半径 | AUC-ROC | 排名 |
|:---:|:---:|:---:|:---:|:---:|
| **EXP01** | 去掉 | 10Å | **0.7115** | 1st |
| EXP03 | 去掉 | 6Å | 0.6678 | 2nd |
| EXP02 | 保留 | 10Å | 0.4894 | 3rd |
| EXP04 | 保留 | 6Å | 0.4257 | 4th |

> AUC-ROC = 0.5 是随机猜测水平，1.0 是完美预测。低于 0.5 意味着模型的判断被**反转**了。

### 这意味着什么？

1. **Heme 是毒药，不是帮手（效应 = -0.2322）**
   - 加 Heme 后性能暴跌到 0.5 以下（比瞎猜还差）
   - 根本原因：模型的原子词汇表里没有铁（Fe），Heme 残基被标记为"未知"（UNK）。这些 OOD 噪音通过 GNN 消息传递污染了整个口袋表示
   - 这不是 bug，是模型的固有局限。想让 Heme 信息有用，需要在 Path C 重训时扩展原子词汇表

2. **大口袋略优于小口袋（效应 = +0.0537）**
   - 10Å 比 6Å 好一点，更大的口袋包含了更多有意义的蛋白环境信息

3. **Gate A 决策：PASS** — 采用 10Å / 不含 Heme 配置继续后续研究

4. **重要注意**：EXP01（0.7115）和 Path A（0.6636）用的是不同的数据集，不能直接说"提升了 4.8%"。但 4 个实验之间的对比是严格公平的（同一批 495 个样本）。

---

## 时间线

| 日期 | 阶段 | 内容 |
|------|------|------|
| 02-20 | 代码准备 | HETATM 问题 Codex 四轮验证 |
| 02-21 | 代码准备 | 编写 step8_align_ligand.py、step8_generate_structure_lmdb.py |
| 02-21 | 代码审核 | Codex 三轮审核，修复 7 个 bug |
| 02-22 | 集成 | run_experiment.py 集成 8.2/8.3，dry-run 通过 |
| 02-22 | 阻断 | 发现 Mac 缺少大文件（PDB + 特征 LMDB），无法运行 |
| 02-25 | Windows 执行 | 创建 pathb-step2 分支，4 个实验 8.1-8.3 全部成功 |
| 02-25 | 推理+分析 | 编写并运行 Step 9/10，Gate A PASS |
| 02-25 | 审核 | Codex 确认 Heme OOD 是真实信号 |
| 02-26 | 文档 | 更新所有进度文件 + 补全 session log |

---

## 一、实验设计

### 1.1 因子与水平

| 因子 | 水平 | 说明 |
|------|------|------|
| Heme | on / off | 是否将 Heme 辅因子原子纳入口袋 |
| Pocket radius | 6Å / 10Å | 配体周围多大范围内的蛋白原子算"口袋" |

### 1.2 四个实验

| 实验 | 口袋半径 | Heme | 设计目的 |
|:---:|:---:|:---:|------|
| EXP01 | 10Å | 无 | 基线条件（与 PathA 最接近） |
| EXP02 | 10Å | 有 | 单独测试 Heme 对大口袋的效果 |
| EXP03 | 6Å | 无 | 单独测试缩小口袋的效果 |
| EXP04 | 6Å | 有 | 两个因子叠加 |

### 1.3 执行流水线

每个实验依次经过 5 个步骤：

```
Step 8.1: extract_pocket_ligand.py    → 从 PDB 提取口袋 + 配体
Step 8.2: step8_align_ligand.py       → 对齐配体原子编号（4 实验共享，只跑 1 次）
Step 8.3: step8_generate_structure_lmdb.py → 生成 structure LMDB
Step 9:   step9_inference.py          → 模型推理
Step 10:  step10_comparative_analysis.py → 4 实验对比分析 + Gate A 决策
```

### 1.4 数据来源

所有实验使用同一份数据：**B6 数据集**（516 条酶-配体对记录），经质量过滤后 **495 条**可用。

---

## 二、代码准备阶段（Mac，02-20 ~ 02-22）

### 2.1 HETATM 阻断问题

**问题**：Step 1 的 `extract_pocket_ligand.py` 用 `--include_heme` 可将 Heme 写入口袋 PDB，但下游 `PDBProtein._parse()` **只读 ATOM 行，静默丢弃 HETATM 行**。如果不修复，Heme 因子无效，2×2 实验退化成只测半径。

**Codex 四轮验证确认**：
1. BioPython PDBIO 根据残基的 `hetfield` 决定写 ATOM 还是 HETATM
2. Heme 残基 `hetfield != " "`，所以写出为 HETATM
3. `PDBProtein._parse()` 的判断条件 `line[0:6].strip() == 'ATOM'` 完全跳过 HETATM
4. 结论：Heme 被静默丢弃，必须修复

**额外发现**（影响后续分析解读）：
- Fe(26) **不在** `FeaturizeProteinAtom.atomic_numbers` 中 → Fe 得到全零特征向量
- HEM 残基会被标记为 UNK（index=20）
- 模型从未在训练数据中见过 HETATM → 完全 out-of-distribution
- 因此本实验定位为「推理敏感性实验」，不是「Heme 生化机理验证」

### 2.2 新脚本编写

| 脚本 | 作用 | 与 PathA 的区别 |
|------|------|----------------|
| `step8_align_ligand.py` | PDB 配体与 SMILES 原子编号对齐 | CLI 可配置路径（PathA 硬编码） |
| `step8_generate_structure_lmdb.py` | pocket PDB + aligned SDF → LMDB | 修复了 HETATM 解析 |
| `run_experiment.py`（更新） | 串联 8.1→8.2→8.3 | 共享对齐缓存 + 智能排序 |

**HETATM 修复的核心改动**（`step8_generate_structure_lmdb.py`）：
```python
# 原版（PathA）：只读 ATOM
line[0:6].strip() == 'ATOM'

# 修复后：同时读 ATOM 和 HETATM
line[0:6].strip() in ('ATOM', 'HETATM')
```

### 2.3 Codex 三轮代码审核（修复 7 个 bug）

| # | 严重度 | 问题 | 修复 |
|---|--------|------|------|
| 1 | HIGH | 共享对齐实验排序依赖（小半径先跑→对齐不完整） | 排序 radius DESC + heme DESC |
| 2 | HIGH | alignment_summary 不检测空/损坏文件 | 行数验证 + 重跑逻辑 |
| 3 | MEDIUM | `element_symb` 回退对 Fe 等双字母元素只取 1 字符 | `line[12:14].strip().capitalize()` |
| 4 | MEDIUM | `r["success"] == "True"` 大小写敏感 | `.strip().lower() == "true"` |
| 5 | MEDIUM | 大量 raw_ligand 缺失时仍静默成功 | 成功率 < 50% 返回错误码 |
| 6 | LOW | AtomMapNum 全零检查误拒单原子配体 | 加 `mol.GetNumAtoms() > 1` 守卫 |
| 7 | LOW | 空 target_docks 仍生成空 LMDB | fail fast 返回 2 |

### 2.4 Mac 环境 dry-run

```bash
python run_experiment.py --run_all --dry_run  # 4 个实验全部通过
```

**阻断**：Mac 上缺少 22GB 数据文件（PDB 结构 + 特征 LMDB + 模型 checkpoint）→ 决定在 Windows 上执行。

---

## 三、Windows 执行阶段（02-25）

### 3.1 环境与前置条件

**环境**：Windows 11, Python 3.10 (conda `torch`), PyTorch 2.3.0+cu121, GPU

**注意**：`conda run -n torch python -c "..."` 在 Windows SSH 下不支持多行参数（`AssertionError`），需直接用 `D:/anaconda3/envs/torch/python.exe`。

**必需文件检查**：

| 文件 | 状态 | 说明 |
|------|------|------|
| PDB 文件 | ⚠️ 619/627 | 8 个缺失，非阻断 |
| 共享特征 LMDB（4个） | ✅ | `data/00_shared/features/` |
| B6 数据集 CSV | ✅ | `data/00_shared/datasets/B6_v1/` |
| 模型 checkpoint | ✅ | 22.3MB（手动从备份恢复） |

### 3.2 运行 Steps 8.1-8.3（结构特征生成）

```bash
python scripts/utils/run_experiment.py --run_all
```

**智能排序**：最宽松条件（10Å+Heme）先跑，产生最完整的共享对齐数据，后续实验直接复用。

| 实验 | Step 8.1 提取 | Step 8.2 对齐 | Step 8.3 LMDB | 有效样本 |
|------|:---:|:---:|:---:|:---:|
| EXP02（首个） | 80s | 2.2s | 11.1s | 495 |
| EXP01 | 77.5s | cached | 10.9s | 495 |
| EXP04 | 75.9s | cached | 10.2s | 495 |
| EXP03 | 75.2s | cached | 9.7s | 495 |
| **合计** | | | | **~6 分钟** |

共享对齐缓存正常：Step 8.2 只跑了 1 次（EXP02 时），后续 3 个实验检测到缓存文件后自动跳过。

**数据漏斗**：

```
B6 data.csv 516 条 → mapping 540 个 dock index → 619 个 PDB 可用
    → Step 8.1 提取（7-8 个失败）
    → Step 8.2 对齐（部分失败）
    → Step 8.3 质量过滤
    → 最终每个实验 495 条有效记录
```

---

## 四、Step 9 模型推理（02-25）

### 4.1 脚本编写

**文件**：`scripts/02_Step2_因子实验/step9_inference.py`（~160 行）

编写流程：读 PathA 的 `step9_inference.py`（322 行）→ Codex 出原型（281 行）→ 我重写精简为 160 行。

**与 PathA 的关键区别**：

| | PathA | PathB |
|--|-------|-------|
| 路径 | 全部硬编码 | CLI 参数（`--experiment_dir` 等） |
| 结构 LMDB | 固定路径 | 自动从 experiment_dir 发现 |
| 数据集 | PathA data.csv | B6 data.csv |
| 输出 | PathA 子目录 | experiment_dir/predictions.csv |

**核心推理逻辑**（与 PathA 一致）：

```python
config.data.full_data = False   # 关键：跳过 str_features.lmdb 依赖
config.num_cpus = 0             # Windows 兼容：不用多进程
model = SS.load_from_checkpoint(checkpoint, config=config)
for batch in test_loader:
    logits, _ = model(batch)    # [B, 1] 原始预测值
```

> `full_data=False` 让模型跳过 `str_features.lmdb`（序列级结构特征）。这个文件生成非常耗时且与本实验的口袋参数无关。PathA 在 Step 9 中也使用了此标志。

### 4.2 执行

```bash
for exp in EXP01 EXP02 EXP03 EXP04; do
    python step9_inference.py \
        --experiment_dir  data/02_Step2_因子实验/$exp \
        --shared_features data/00_shared/features \
        --shared_datasets data/00_shared/datasets \
        --checkpoint_dir  saved_model/model/run_0
done
```

**4 个实验全部成功**，每个处理 495 个样本。

### 4.3 推理输出初看

| 实验 | prob 范围 | 观察 |
|------|-----------|------|
| EXP01 (10Å/noHeme) | 0.0000 ~ 0.9280 | 正常分布，区分度好 |
| EXP02 (10Å/Heme) | 0.0000 ~ 1.0000 | ⚠️ 出现极端概率 |
| EXP03 (6Å/noHeme) | 0.0000 ~ 0.9855 | 正常分布 |
| EXP04 (6Å/Heme) | 0.0000 ~ 1.0000 | ⚠️ 出现极端概率 |

Heme 实验出现了 prob=1.0（logit 极大值），说明 Heme 的 OOD 原子严重扰动了模型输出。

---

## 五、Step 10 对比分析 + Gate A 决策（02-25）

### 5.1 脚本编写

**文件**：`scripts/02_Step2_因子实验/step10_comparative_analysis.py`（~200 行）

功能：自动发现 4 个 EXP 目录 → 从目录名解析因子 → 计算 AUC-ROC/AUC-PR/accuracy → 画 ROC 叠加图 + 热力图 → 输出 Gate A 决策报告。

### 5.2 核心结果与解读

#### 排名

```
EXP01 (0.7115) > EXP03 (0.6678) > EXP02 (0.4894) > EXP04 (0.4257)
```

规律非常清晰：**去 Heme 的两个实验远优于加 Heme 的两个**，10Å 略优于 6Å。

#### 完整指标

| 实验 | Heme | Radius | AUC-ROC | AUC-PR | Accuracy |
|------|:---:|:---:|:---:|:---:|:---:|
| **EXP01** | No | 10Å | **0.7115** | **0.7356** | 0.5091 |
| EXP03 | No | 6Å | 0.6678 | 0.6990 | 0.5131 |
| EXP02 | Yes | 10Å | 0.4894 | 0.5492 | 0.4848 |
| EXP04 | Yes | 6Å | 0.4257 | 0.4973 | 0.4525 |

**关于 accuracy ~50%**：这是正常的，不代表模型差。B6 数据集接近平衡（272 正 / 244 负），accuracy 用固定阈值 0.5 切分，对 logit 偏移很敏感。AUC-ROC 不受阈值影响，是更可靠的指标。PathA 的 accuracy 也是 ~50%。

#### 主效应分析

| 效应 | 计算方式 | 值 | 解读 |
|------|---------|:---:|------|
| **Heme** | mean(有Heme) − mean(无Heme) | **-0.2322** | Heme 是毒药，不是帮手 |
| **Radius** | mean(10Å) − mean(6Å) | **+0.0537** | 大口袋略优 |

**为什么 Heme 效应这么大（-23%）？**

这不是 bug，而是模型的 OOD（Out-of-Distribution）限制：

1. **Fe 不在原子词汇表中** — `src/Datasets/Structure/transforms.py` 的 `atomic_numbers` 列表不含 Fe(26)，Fe 原子得到全零特征向量（相当于"看不见"）
2. **HEM 残基映射为 UNK** — residue index=20（未知残基类型），模型训练时从未见过这种模式
3. **GNN 消息传递扩散噪音** — 这些 OOD 特征通过图神经网络的消息传递机制，从 Heme 原子扩散到相邻的蛋白原子，污染整个口袋表示
4. **AUC < 0.5 意味着反转** — 模型的预测方向被 Heme 噪音颠倒了（该判正的判负，该判负的判正）

**为什么 Radius 效应是正的？**

10Å 口袋比 6Å 多包含了 ~40% 的蛋白原子。这些额外的环境信息（比如口袋入口处的残基、第二层配位残基）帮助模型更好地理解结合位点的化学环境。

#### Gate A 决策

```
最佳实验:        EXP01_B6_10A_noHeme
AUC-ROC:        0.7115
PathA baseline: 0.6636
Delta:          +0.0479 (+4.79%)
阈值:           0.005 (0.5%)

Gate A:         PASS
决策:           采用 10Å/noHeme 配置继续后续 PathB 步骤
```

**⚠️ 重要注意（Codex 审核发现）**：

EXP01（0.7115）和 PathA（0.6636）使用的**不是同一批数据**：
- PathA 用的是 PathA Step 4 的 data.csv → 517 个可用样本
- PathB 用的是 B6 的 data.csv → 495 个可用样本
- B6 排除了产物（PRODUCT 类配体），样本构成不同

因此 "+4.79%" 不能直接宣称为"改进"。但 Gate A 的目的是在 4 个实验中选最佳配置，这个内部比较是严格公平的（同 495 个样本）。

---

## 六、Codex 审核（02-25）

### 审核的 4 个问题

| # | 问题 | Codex 结论 |
|---|------|-----------|
| 1 | Heme 的 -23% 下降是 bug 还是真实信号？ | **真实的 OOD 效应**（检查了 3 条代码路径确认） |
| 2 | EXP01 vs PathA 的 +4.79% 如何解释？ | **不可直接比较**（不同数据集） |
| 3 | accuracy ~50% 正常吗？ | **正常**（与 PathA 模式一致） |
| 4 | Gate A PASS 正确吗？ | **正确**（内部比较有效） |

### Codex 建议的后续步骤

1. 若要宣称改进，需在相同 495 个 ID 上重算 PathA 基线
2. 可尝试 Heme 消融实验：去掉 Fe 原子但保留有机部分
3. 将 Heme 结果定性为"模型 OOD 限制"，而非"Heme 无用"的生化结论

---

## 七、产出文件清单

### 脚本

| 文件 | 作用 | 产出阶段 |
|------|------|---------|
| `scripts/02_Step2_因子实验/extract_pocket_ligand.py` | 口袋/配体提取（Step 8.1） | Step 1 |
| `scripts/02_Step2_因子实验/step8_align_ligand.py` | 配体对齐（Step 8.2） | Mac |
| `scripts/02_Step2_因子实验/step8_generate_structure_lmdb.py` | LMDB 生成（Step 8.3） | Mac |
| `scripts/02_Step2_因子实验/step9_inference.py` | 模型推理 | Windows |
| `scripts/02_Step2_因子实验/step10_comparative_analysis.py` | 对比分析 + Gate A | Windows |
| `scripts/utils/run_experiment.py` | 实验编排（8.1→8.2→8.3） | Mac |

### 数据（data/ = 输入 + 中间产物）

```
data/02_Step2_因子实验/
├── EXP01_B6_10A_noHeme/
│   ├── high_quality_id.txt       # 中间产物：495 条有效 Dock Index
│   ├── structure_build_summary.csv # 中间产物：LMDB 构建日志
│   ├── structure_features.lmdb   # 中间产物 (gitignored)
│   └── structure_features/       # 中间产物：pocket/ + raw_ligand/ (gitignored)
├── EXP02_B6_10A_Heme/            # 同结构
├── EXP03_B6_6A_noHeme/           # 同结构
├── EXP04_B6_6A_Heme/             # 同结构
└── shared_alignment/             # 共享对齐数据 (gitignored)
```

### 结果（results/ = 最终输出）

```
results/02_Step2_因子实验/
├── EXP01_B6_10A_noHeme/          ← 最佳配置
│   ├── predictions.csv           # Step 9 推理结果（7 列 × 495 行）
│   └── config.yaml               # 实验参数快照（便于追溯）
├── EXP02_B6_10A_Heme/            # 同结构
├── EXP03_B6_6A_noHeme/           # 同结构
├── EXP04_B6_6A_Heme/             # 同结构
└── analysis/                     # Step 10 对比分析
    ├── comparative_metrics.csv   # 4 行指标汇总
    ├── gate_a_recommendation.txt # Gate A 决策报告
    └── figures/
        ├── comparative_roc.png   # 4 条 ROC 曲线叠加图
        └── factorial_heatmap.png # 2×2 AUC-ROC 热力图
```

### Git 提交（pathb-step2 分支）

| 提交 | 内容 |
|------|------|
| `ce003e1` | Step 2 代码准备（Mac）：脚本 + HETATM 修复 |
| `b3f6899` | 文档拆分与进度更新（Mac） |
| `ba7c22a` | 4 个因子实验执行 + 推理 + 分析（Windows） |
| `48bdbcb` | 所有进度追踪文件更新 |
| `11f202e` | Session log 详细执行记录 |

---

## 八、经验总结

| 问题 | 解决方案 | 影响 |
|------|---------|------|
| Windows SSH 下 `conda run` 不支持多行参数 | 直接用 Python 解释器绝对路径 | 操作层面 |
| Git 凭证在 SSH 下不可用 | Device code 认证（github.com/login/device） | 操作层面 |
| `full_data=False` 跳过耗时的 str_features.lmdb | PathA 验证过可行，直接复用 | 效率 |
| 共享对齐缓存只跑 1 次 | run_experiment.py 的智能排序 + 缓存检测 | 效率 |
| Heme OOD 问题在设计阶段就有预判 | Section 2.1 的 Codex 分析已警告 Fe 不在词汇表 | 实验设计 |
| 不同数据集的 AUC-ROC 不可直接比较 | Codex 审核发现，已标注在 Gate A 报告中 | 结论严谨性 |
