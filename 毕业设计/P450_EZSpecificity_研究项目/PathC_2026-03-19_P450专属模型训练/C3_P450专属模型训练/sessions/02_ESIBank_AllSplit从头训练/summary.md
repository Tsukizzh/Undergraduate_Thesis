# Step 2: ESIBank AllSplit 从头训练

> **对应实验**: Path B Step 9-10
> **日期**: 2026-03-12 ~ 2026-03-16
> **详细 log**: `PathB_.../sessions/09_Step9_AllSplit训练/` + `10_Step10_pt训练管线/`

## 做了什么

在 ESIBank 原始数据上，使用 all_split（Unknown Enzyme+Substrate，0% 酶/底物重叠）从头训练 EZSpecificity 模型，复现论文基线。

## 为什么做

1. 验证我们的训练管线（.pt 缓存 + DDP）能否复现论文结果
2. 建立 ESIBank 数据的基线，为后续 P450 专属训练提供对比参考

## 关键结果

| 配置 | 环境 | Val AUC | Test AUC | 说明 |
|------|------|---------|----------|------|
| 本地 fixed | 单卡 4070S | **0.7667** (ep14) | 未测 | 27 epochs |
| Cloud-2 legacy_bug | 2×4090 DDP | 0.722 (ep22) | **0.7244** (ep27) | 32 epochs |
| **论文报告** | 4 GPU, ~256 ep | — | **0.7198** | |

## 核心发现

1. **我们 0.7244 > 论文 0.7198**（仅 32 epochs vs 论文 ~256 epochs）
2. .pt 缓存管线成功解决了 LMDB thrashing（7.56 it/s stable vs 2.09 it/s decay）
3. DDP 训练管线验证通过（发现并修复了 test_epoch_end all_gather 缺失的 bug）
4. 边修复（fixed vs legacy_bug）对 AUC 影响不显著

## 技术贡献

- **.pt 缓存 v3**: per-sample .pt 文件，~10GB 工作集，无 thrashing
- **三步构建流程**: graph shards → per-sample .pt → flatbin (enzymes.bin + substrates)
- **DDP bug 修复**: validation_epoch_end 有 all_gather 但 test_epoch_end 没有
