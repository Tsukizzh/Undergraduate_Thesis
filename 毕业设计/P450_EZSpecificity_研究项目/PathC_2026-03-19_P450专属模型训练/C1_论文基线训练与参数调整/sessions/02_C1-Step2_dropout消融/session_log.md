# C1-Step 2: Dropout 消融实验

## 实验概要

| 项目 | 值 |
|------|-----|
| 日期 | 2026-03-20 ~ 2026-03-21 |
| 服务器 | Cloud-2 (2×RTX4090, DDP) |
| 基线 | C1-Step 1 fixed (Val AUC=0.7145, ep16) |
| 改动 | cross_attention.dropout: 0.9 → 0.1 / 0.3 |
| 其他参数 | 与基线完全一致 (bs=56×2, lr=3e-4, warmup=8, edge=fixed) |

## 实验结果

| 实验 | Dropout | Best Val AUC | Best Epoch | 总Epochs | 早停? | vs 基线 |
|------|---------|-------------|------------|----------|-------|---------|
| fixed基线 | 0.9 | 0.7145 | ep16 | 32 (早停ep31) | 是 | - |
| **EXP-A** | **0.1** | **0.7216** | **ep17** | 32 (早停ep32) | 是 | **+0.0071** |
| **EXP-B** | **0.3** | **0.7397** | **ep49** | 50 (跑满) | 否 | **+0.0252** |

### 关键发现

1. **dropout=0.3 显著优于 0.1 和基线**：AUC 0.7397 vs 0.7216 vs 0.7145
2. **dropout=0.3 在ep49仍在上升**，未触发早停，说明还有提升空间 → 需要延长epochs
3. dropout=0.1 早停于ep32（patience=15, best=ep17），略好于基线但不显著
4. 验证了原始 dropout=0.9 确实是主要瓶颈之一

### AUC 趋势对比

**dropout=0.1**:
- ep0(0.597) → ep7(0.680) → **ep10(0.720)** → ep17(**0.722 BEST**) → ep23(0.717) → ep32(早停)
- 特点：快速收敛，但plateau明显，ep10后几乎不涨

**dropout=0.3**:
- ep0(0.597) → ep10(0.686) → ep20(0.704) → ep31(0.723) → ep34(0.727) → ep48(0.728) → **ep49(0.740 BEST)**
- 特点：收敛较慢但持续上升，ep30后仍在稳步提升
- **最后3个epoch: 0.698 → 0.728 → 0.740**，上升趋势明显

### Per-family AUC (Best Epoch)

| 家族 | dropout=0.1 (ep17) | dropout=0.3 (ep49) | fixed基线 (ep16) |
|------|-------------------|-------------------|-----------------|
| brenda | 0.736 | 0.778 | 0.736 |
| Esterase | 0.604 | 0.640 | 0.575 |
| Gt_acceptor | 0.646 | 0.292 | 0.448 |
| Nitrilase | 0.000 | 1.000 | 1.000 |
| Phosphatase | 0.509 | 0.531 | 0.487 |

### 训练曲线特征

- **Val Loss ↑ while AUC ↑**：两个实验都出现此现象（与Step 1一致，正常）
- **Gradient Norm**：随epoch增长（0.47→1.70），dropout=0.3略高于0.1
- **Train Loss**：dropout=0.3最终更低（0.144 vs 0.160），dropout降低允许更多信号通过

## 配置详情

```yaml
# 改动项
cross_attention.dropout: 0.1 / 0.3  # 原始: 0.9

# 不变项
edge_mode: fixed
batch_size: 56 (per GPU) × 2 GPU = 112 effective
lr: 3e-4
warmup_epochs: 8
max_epochs: 50
patience: 15
gradient_clip_val: 8
seed: 3407
```

## 文件产出

```
results/02_dropout_ablation/
├── dropout01/
│   ├── metrics_dropout01.csv
│   ├── fig1_training_dynamics.png
│   ├── fig2_family_performance.png
│   └── checkpoints/
│       └── pt-dropout01-ep17-auc0.7216.ckpt (22MB)
└── dropout03/
    ├── metrics_dropout03.csv
    ├── fig1_training_dynamics.png
    ├── fig2_family_performance.png
    └── checkpoints/
        └── pt-dropout03-ep49-auc0.7397.ckpt (22MB)
```

## Test 集评估 (2026-03-21)

所有 3 个 checkpoint 在 8841 样本测试集上评估。

### 总体 Test 结果

| 实验 | Dropout | Val AUC (best ep) | **Test AUC** | **Test AUPR** | vs 基线 Test |
|------|---------|-------------------|-------------|--------------|-------------|
| fixed基线 | 0.9 | 0.7145 (ep16) | **0.7060** | 0.2362 | - |
| EXP-A | 0.1 | 0.7216 (ep17) | 0.6936 | 0.2012 | **-0.012** |
| EXP-B | 0.3 | 0.7397 (ep49) | 0.6959 | 0.2044 | **-0.010** |

### 关键发现

1. **Val 改善未迁移到 Test**：dropout=0.3 Val +0.025 但 Test -0.010，dropout=0.1 Val +0.007 但 Test -0.012
2. **三个模型 Test AUC 均低于 legacy_bug baseline (0.7244)**
3. **Val-Test gap 随 dropout 降低而增大**：基线 gap=0.009, dropout=0.1 gap=0.028, dropout=0.3 gap=0.044
4. **结论**：降低 dropout 改善了验证集拟合但加剧了过拟合到验证集分布

### Per-family Test AUC

| 家族 | fixed (d=0.9) | dropout=0.1 | dropout=0.3 |
|------|--------------|-------------|-------------|
| brenda | 0.7279 | 0.7119 | 0.7199 |
| Duf | 0.3947 | 0.3891 | **0.6767** |
| Esterase | 0.7124 | **0.7312** | 0.6571 |
| Gt_acceptor | **0.7431** | **0.7896** | 0.7025 |
| Nitrilase | **0.8333** | 0.7778 | 0.6111 |
| Phosphatase | 0.4618 | **0.5047** | 0.4423 |
| Thiolase | 0.6364 | **0.7273** | 0.6364 |

- dropout=0.1 在多数小家族上优于基线（Esterase, Gt_acceptor, Phosphatase, Thiolase）
- dropout=0.3 在 Duf 家族上大幅提升（0.39→0.68），但在其他家族上普遍下降

## 数据泄露分析 (2026-03-21)

### ESIBank all_split 泄露验证

| 数据集 | 酶重叠 | 底物重叠 |
|--------|--------|---------|
| **ESIBank brenda single-family all_split** | **0%** | **0%** |
| ESIBank random_split | 99.9% | 98.1% |

- 所有 7 个家族单独检查：均为 0% 酶重叠、0% 底物重叠
- **跨家族 ID 重用**（合并后）：酶 3.2%、底物 1.6%（解释了 .pt cache 中看到的 8-10%）
- 结论：all_split 模式下无数据泄露，之前的 6.56% 底物泄露来自跨数据集 Morgan FP 重叠（非 ID 重叠）

## 代码变更 (2026-03-21)

- `run_test_eval.py` 功能合并到 `main_training_pt.py` 的 `--test-only` 模式
- `start_ablation.sh` 已删除（全部用 Python 脚本）
- 清理了服务器临时文件（/tmp/、/root/ablation.log）和本地临时文件

## 下一步

1. 继续其他消融项（lr、weight_decay、LR scheduler 等）
2. dropout 结论：**降低 dropout 对 Val 有效但未迁移到 Test，不作为 C1-Step 3 组合改动**
3. 需要重新审视优化方向：关注能改善 Test 泛化的改动（如 weight_decay、数据增强）
