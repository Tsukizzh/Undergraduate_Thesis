# Step 5: 重构评估 (Feature Generation + Inference + Gate B)

**日期**: 2026-03-03
**状态**: COMPLETED
**Gate B 结论**: INFORMATIVE FAIL (AUC-ROC = 0.5170)

## 目标

对 Step 4 生成的 2,766 个 Vina 对接复合物（265 正 + 2,501 随机负）运行 EZSpecificity 完整推理管线，
评估随机负样本策略是否改善 AUC-ROC（相比 Path A 的抑制剂负样本）。

## 执行记录

### Phase 0: 资产检查与 Codex 审核

**检查项**:
- Step 4 输出: data.csv (2,766 rows), Enzymes.csv (292), Substrates.csv (436)
- pocket/ (2,766 PDB), raw_ligand/ (2,766 SDF)
- 共享特征: ESM (292 entries), Reaction (436), GROVER (436), Morgan (436×1024)
- 确认酶/底物索引与 B6/PathA 完全一致，可复用全部 4 种特征

**Codex 发现的风险**:
- 列名分裂: alignment 需要 `Dock_Index`（下划线），LMDB/inference 需要 `Dock Index`（空格）
- step10_comparative_analysis.py 不适用于 Gate B（硬编码 EXP* 命名）→ 需要独立的 Gate B 脚本
- LMDB Unicode 路径问题
- BOM 编码风险

### Phase 1: 目录创建与 Mapping CSV

- 创建 `data/05_Step5_重构评估/`, `results/05_Step5_重构评估/`, `sessions/05_Step5_重构评估/`
- 生成 `dock_index_mapping.csv`: 从 data.csv 提取，列名重命名为下划线格式，UTF-8 无 BOM

### Phase 2: 配体对齐

- 脚本: `scripts/02_Step2_因子实验/step8_align_ligand.py`（复用）
- 结果: **2,766/2,766 = 100%** 成功，耗时 8 秒
- 输出: `data/05_Step5_重构评估/aligned_ligand/` (2,766 SDF)

### Phase 3: 结构 LMDB

- 脚本: `scripts/02_Step2_因子实验/step8_generate_structure_lmdb.py`（复用）
- LMDB Unicode 路径修复: `subdir=False`（与 src 代码一致）
- 结果: **2,766/2,766 = 100%** 成功，耗时 23 秒
- HETATM=0（正确：noHeme 配置，Gate A 决策）
- 输出: `data/05_Step5_重构评估/structure_features.lmdb` (~10GB)

### Phase 4: 模型推理

- 脚本: `scripts/02_Step2_因子实验/step9_inference.py`（复用）
- GPU 利用率 ~0%（正常：小模型 + CPU-bound 数据加载）
- 结果: **2,766/2,766 = 100%** 推理完成
- 输出: `results/05_Step5_重构评估/predictions.csv`
- prob 范围: [0.0000, 0.9910]

### Phase 5: AUC-ROC = 0.5170 — 诊断

**关键发现: AUC-ROC = 0.5170（近似随机）**

诊断分析:
- 正样本均分: -2.94, 负样本均分: -3.20 → **几乎无分离** (separation=0.26)
- EXP01 对比: 正样本 -3.00（一致），负样本 -5.51（明显分离 = 2.51）
- Per-enzyme AUC: mean=0.4769, 仅 42.1% > 0.5
- Codex 确认：非管线错误，是真实结果

**根因假设**:
1. Dockability ≠ Catalysis: Vina 给所有配对生成合理 pose → 模型无法区分
2. OOD 酶 (0% overlap): 模型从未见过这 152 个 P450 酶
3. 混杂变量: 同时改变了负样本定义 + 结构来源

### Phase 6: Gate B 分析脚本 + 报告

- 创建 `scripts/05_Step5_评估/step5_gate_b_analysis.py`
- Codex 审核 (2 轮):
  - Round 1: 6 个问题 (硬编码值、因果过强、缺 CI、低支持度、表格不完整、边界)
  - Round 2: 确认全部修复，评为 "chapter-draft ready"
- 修复: 添加 bootstrap CI (2000 iter)、最低支持度过滤 (n>=3)、软化因果语言、2×2 混杂表、方法学注释

**最终 Gate B 指标**:

| Dataset | AUC-ROC | 95% CI | AUC-PR | Prevalence |
|---------|---------|--------|--------|------------|
| Step 5 (Random, Vina) | 0.5170 | [0.4804, 0.5521] | 0.1116 | 0.096 |
| EXP01 (Inhibitor, Crystal) | 0.7115 | [0.6657, 0.7569] | 0.7356 | 0.539 |

## 产出文件

| 文件 | 说明 |
|------|------|
| `data/05_Step5_重构评估/dock_index_mapping.csv` | Dock Index 映射 |
| `data/05_Step5_重构评估/aligned_ligand/` | 2,766 对齐配体 SDF |
| `data/05_Step5_重构评估/structure_features.lmdb` | 结构特征 LMDB |
| `data/05_Step5_重构评估/high_quality_id.txt` | 2,766 高质量 ID |
| `results/05_Step5_重构评估/predictions.csv` | 模型预测结果 |
| `results/05_Step5_重构评估/analysis/gate_b_report.md` | Gate B 决策报告 |
| `results/05_Step5_重构评估/analysis/gate_b_metrics.csv` | 指标汇总 |
| `results/05_Step5_重构评估/analysis/gate_b_roc.png` | ROC 曲线对比图 |
| `results/05_Step5_重构评估/analysis/gate_b_score_dist.png` | 分数分布对比图 |
| `results/05_Step5_重构评估/analysis/gate_b_per_enzyme_auc.png` | 每酶 AUC 分布图 |
| `results/05_Step5_重构评估/analysis/per_enzyme_auc.csv` | 每酶 AUC 数据 |
| `scripts/05_Step5_评估/step5_gate_b_analysis.py` | Gate B 分析脚本 |

## 关键决策

| 决策 | 选择 | 原因 |
|------|------|------|
| 特征复用 | 复用 PathA/B6 的 ESM/Reaction/GROVER/Morgan | 酶/底物索引完全一致 |
| LMDB 路径 | `subdir=False` | 解决 Windows Unicode 路径问题 |
| Gate B 结论 | INFORMATIVE FAIL | AUC-ROC CI 包含 0.5 |
| 因果表述 | "假设"而非"证明" | 存在负样本类型/结构来源混杂 |

## 下一步

1. Path C: 使用累积数据进行 P450 专项微调
2. 可选控制实验: 晶体正样本 + Vina 随机负样本（分离结构来源效应）
3. 更新全局日志并提交推送
