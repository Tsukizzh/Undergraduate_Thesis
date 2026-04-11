# EXP002b: 学习率调参实验

> 基于EXP002a (Test AUC=0.7816)，调整训练超参数

## 实验配置

**与EXP002a的区别（仅调参，代码和数据完全相同）**：

| 参数 | EXP002a | EXP002b | 调整理由 |
|------|---------|---------|---------|
| lr | 3e-4 | **4e-4** | effective batch=224 是论文64的3.5倍，线性缩放 |
| warmup_epochs | 8 | **12** | 更高lr需要更长warmup |
| weight_decay | 0 | **1e-5** | 轻度正则化（原论文不用） |
| sched_patience | 10 | **8** | 更积极降lr |
| accumulate_grad_batches | 1 | **2** | 新服务器2卡，保持effective batch=224 |
| devices | 4 | **2** | 新服务器 |
| num-workers | 6 | **5** | 14核÷2GPU |

**不变的**：dropout=0.9, hidden_dim=128, k=48, max-epochs=200, gradient_clip=8, 数据=pt_cache_heme/random (feature_dim=31)

## 服务器

新Cloud-2: 2×RTX4090 24GB, 180GB RAM, 14核

## 结果

| 指标 | EXP002a | EXP002b | 差异 |
|------|---------|---------|------|
| Best epoch | 99 | **70** | 更早收敛 |
| Val AUC | 0.7784 | **0.7809** | +0.003 |
| **Test AUC** | 0.7816 | **0.7889** | **+0.007** |
| Test AUPR | 0.3524 | **0.3575** | +0.005 |
| Total epochs | 114 | 86 | 更快停止 |
| 速度 | ~1.4 min/ep (4卡) | ~2.4 min/ep (2卡) | |
| 总耗时 | ~2.5小时 | ~3.5小时 | |

## 性能提升路径

```
0.638 (论文模型, ESIBank P450 内部)
  → +0.086 → 0.7244 (ESIBank AllSplit 从头训练)
    → +0.049 → 0.7730 (P450 专属数据集 EXP001)
      → +0.009 → 0.7816 (Fe/HEM编码 EXP002a)
        → +0.007 → 0.7889 (调参 EXP002b)
```

## 分析

- lr从3e-4提到4e-4 + weight_decay=1e-5 + warmup=12的组合有效
- 收敛更快（ep70 best vs ep99），说明更高的lr加速了有效学习
- Test AUC和AUPR都有提升，说明不是单纯过拟合
- 但无法区分lr/wd/warmup/patience各自的贡献（一次改了4个参数）

## 文件

- 结果: `results/EXP002b_lr_tuning/`
- 服务器目录: `P450/experiments/EXP002b_lr_tuning/`
