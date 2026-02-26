# Step 2 状态

> 上次更新：2026-02-26
> 状态：**已完成** ✅

---

## 结果摘要

| 实验 | Heme | Radius | AUC-ROC | AUC-PR | 样本数 |
|------|------|--------|---------|--------|--------|
| **EXP01** | No | 10Å | **0.7115** | 0.7356 | 495 |
| EXP03 | No | 6Å | 0.6678 | 0.6990 | 495 |
| EXP02 | Yes | 10Å | 0.4894 | 0.5492 | 495 |
| EXP04 | Yes | 6Å | 0.4257 | 0.4973 | 495 |

**Gate A 决策**: PASS — 采用 10Å/noHeme 配置
**Heme 效应**: -0.2322（严重负面，OOD：Fe 不在原子词汇表）
**Radius 效应**: +0.0537（10Å 优于 6Å）

**Codex 审核备注**: EXP01(0.7115) vs PathA(0.6636) 使用不同数据集，不可直接对比

---

## 全部完成项

| 项目 | 状态 | 说明 |
|------|------|------|
| Step 1 数据迁移 | ✅ | B6 数据集 + 共享特征 |
| extract_pocket_ligand.py | ✅ | --pocket_radius + --include_heme |
| step8_align_ligand.py | ✅ | CLI 可配置 |
| step8_generate_structure_lmdb.py | ✅ | HETATM 修复 |
| run_experiment.py | ✅ | 集成 8.1→8.2→8.3 |
| Codex 三轮代码审核 | ✅ | 修复 7 个 bug |
| 4 个因子实验 (8.1-8.3) | ✅ | EXP01-04 全部成功 (Windows) |
| step9_inference.py | ✅ | CLI 可配置推理脚本 |
| Step 9 推理 (4个实验) | ✅ | 每个 495 样本 |
| step10_comparative_analysis.py | ✅ | ROC/heatmap/Gate A |
| Step 10 对比分析 | ✅ | Gate A: PASS |
| Codex 结果审核 | ✅ | Heme OOD 确认 |
| Git 提交 | ✅ | ba7c22a on pathb-step2 |

---

## 下一步: Step 3

AutoDock Vina 对接管线，生成随机负样本
