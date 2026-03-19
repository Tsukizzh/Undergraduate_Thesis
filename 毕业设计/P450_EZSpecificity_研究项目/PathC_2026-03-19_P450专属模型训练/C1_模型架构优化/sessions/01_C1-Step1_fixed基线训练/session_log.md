# C1-Step 1: Fixed 基线训练 — Session Log

> **日期**: 2026-03-19
> **目标**: 在 Cloud-2 上用修复边排序的 `--edge-mode fixed` 训练 AllSplit 基线，与 legacy_bug 基线（test AUC=0.7244）对比
> **结论**: Val AUC **0.7145** (ep16 best)，Test AUC **0.7060**。均低于 legacy_bug（Val 0.722, Test 0.7244）。边排序修复未带来提升，Δ=-0.018 (test)。

---

## 1. 实验配置

| 参数 | 值 |
|------|------|
| 服务器 | Cloud-2 (2×RTX 4090, 180GB RAM) |
| 训练模式 | DDP (2 GPU) |
| Edge mode | **fixed**（边排序修复版） |
| Batch size | 56/GPU × 2 = **112 effective** |
| Accumulate grad | 1 |
| 学习率 | 3e-4 (warmup 8 epochs) |
| LR scheduler | ReduceLROnPlateau(monitor=aupr, patience=10, factor=0.5) |
| Weight decay | 0 |
| Dropout | 0.9（交叉注意力层） |
| 精度 | FP16 混合精度 |
| 最大轮数 | 50 |
| 早停 | monitor=auc/val, patience=15 |
| Seed | 3407 |
| 自动关机 | ✅ `--shutdown` |

**脚本路径（服务器）**: `PathC/scripts/main_training_pt.py`
**配置文件**: `PathC/configs/base.yml`
**缓存目录**: `data/allsplit_pt_cache/` (57GB, 176K .pt 文件)

---

## 2. 服务器重组

本次训练前对 Cloud-2 服务器进行了目录重组：

**旧结构**（扁平）:
```
EZSpecificity/
├── scripts/10_Step10_pt训练管线/
├── results/10_Step10_pt训练/
├── logs/10_Step10_pt训练/
└── data/10_Step10_pt训练/allsplit_pt_cache/
```

**新结构**（Path 分离）:
```
EZSpecificity/
├── src/                        ← 论文源码
├── data/allsplit_pt_cache/     ← 共享缓存（提升至项目根）
├── PathB/                      ← 归档 legacy_bug 结果
│   ├── scripts/, configs/, results/, logs/
└── PathC/                      ← 当前工作
    ├── scripts/                ← 从 PathB 复制，路径已修正
    ├── configs/base.yml
    ├── results/                ← 训练输出
    └── logs/                   ← 训练日志
```

**路径修正**（PathC `main_training_pt.py`）:
- `PATHC_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, ".."))`
- `PROJECT_ROOT = os.path.normpath(os.path.join(PATHC_DIR, ".."))`
- `SRC_DIR = os.path.join(PROJECT_ROOT, "src")`
- `RESULTS_DIR = os.path.join(PATHC_DIR, "results")`
- `LOG_DIR = os.path.join(PATHC_DIR, "logs")`

`pt_dataset.py` 无需修改（已有多路径自适应逻辑）。

---

## 3. 训练过程

### 3.1 启动过程

训练经历了3次启动：
1. **第1次**: 无 `--shutdown`，手动停止
2. **第2次**: 加 `--shutdown --resume last`，但无 checkpoint 可恢复，脚本正常从头训练，完成后触发关机，服务器自动重启
3. **第3次**: 加 `--shutdown`，无 `--resume`，从头训练，成功完成

### 3.2 训练时间线

- **启动**: 2026-03-18 22:36 (第3次)
- **ep0 完成**: 22:46 (10:03)
- **ep16 (best)**: 2026-03-19 01:23
- **ep27 (LR decay)**: 03:10 (ReduceLROnPlateau 触发, lr 3e-4 → 1.5e-4)
- **ep31 (早停)**: 03:49
- **总时长**: ~5.2 小时 (32 epochs)
- **自动关机**: 03:50 确认

### 3.3 速度

| 指标 | 值 |
|------|------|
| 平均速度 | 2.7 it/s |
| 每 Epoch | ~9:45 |
| 总训练时间 | ~5.2h |
| GPU利用率 | 70-100% |
| 显存占用 | ~23.7 GB/卡 |

---

## 4. 训练结果

### 4.1 AUC 趋势

| Epoch | Val AUC | Val Loss | Train Loss | LR | 备注 |
|-------|---------|----------|------------|------|------|
| 0 | 0.5971 | 0.3195 | 0.3621 | 5e-6 | warmup |
| 3 | 0.6190 | 0.3420 | 0.3110 | 1.16e-4 | warmup |
| 4 | 0.6390 | 0.3378 | 0.2967 | 1.52e-4 | |
| 6 | 0.6466 | 0.3418 | 0.2721 | 2.26e-4 | |
| 7 | 0.6690 | 0.3557 | 0.2621 | 2.63e-4 | |
| 8 | 0.6892 | 0.3299 | 0.2615 | 3e-4 | warmup结束 |
| 11 | 0.6954 | 0.3485 | 0.2401 | 3e-4 | |
| 14 | 0.6996 | 0.3892 | 0.2246 | 3e-4 | |
| **16** | **0.7145** | **0.3794** | **0.2195** | **3e-4** | **BEST** |
| 19 | 0.7057 | 0.3750 | 0.2094 | 3e-4 | |
| 23 | 0.7138 | 0.3987 | 0.1974 | 3e-4 | |
| 27 | 0.7018 | 0.4164 | 0.1877 | 3e-4 | LR decay触发 |
| 28 | 0.7010 | 0.4673 | 0.1668 | 1.5e-4 | |
| 31 | 0.7038 | 0.5042 | 0.1552 | 1.5e-4 | 早停 |

### 4.2 过拟合信号

- **Val Loss ↑ while AUC 平台**: ep16后 val_loss 从 0.379 → 0.504（+33%），但 AUC 维持在 0.70-0.71
- **Train Loss 持续下降**: 0.362 → 0.155，典型过拟合
- 与 legacy_bug 训练的行为模式一致（BCE pointwise vs AUC pairwise ranking，已知现象）

### 4.3 Checkpoints（Top-3）

| 文件 | Epoch | Val AUC |
|------|-------|---------|
| `pt-fixed-ep16-auc0.7145.ckpt` | 16 | **0.7145** |
| `pt-fixed-ep23-auc0.7138.ckpt` | 23 | 0.7138 |
| `pt-fixed-ep24-auc0.7122.ckpt` | 24 | 0.7122 |

---

## 5. 与 Legacy Bug 基线对比

| 指标 | Legacy Bug (PathB) | Fixed (PathC C1-Step1) | 差值 |
|------|-------------------|----------------------|------|
| Val AUC (best) | 0.7224 (ep22) | 0.7145 (ep16) | **-0.0079** |
| Test AUC | 0.7244 (ep27) | **0.7060** (ep16) | **-0.0184** |
| Best epoch | 22 | 16 | -6 |
| 早停 epoch | 32 | 31 | -1 |
| 总训练时间 | ~5.3h | ~5.2h | ≈ |
| Effective BS | 112 | 112 | = |

### 5.1 关键发现

**Fixed 修复并未带来 AUC 提升**，Val 低 0.008，Test 低 0.018。可能的解释：

1. **随机波动**: 0.018 差距可能在单次实验的噪声范围内，需多 seed 验证
2. **"有益噪声"假说**: legacy_bug 的错误边特征分配可能意外引入了某种正则化效果
3. **Batch size 问题**: 本地 Step 9 fixed (bs=32) 达到了 0.7667，同一修复在大 batch (112) 下只有 0.7145 → 暗示 **batch size 才是主要瓶颈**

### 5.2 Test AUC 评估

使用 best checkpoint (ep16) 在 AllSplit fold 0 test set 上评估：

| 指标 | 值 |
|------|------|
| Test AUC-ROC | **0.7060** |
| Test AUPR | 0.2362 |
| Test 样本数 | 8,841 |
| 正/负样本 | 966 / 7,875 |
| 评估方式 | 单GPU推理，bs=56 |

### 5.2 Edge Mode 验证

确认训练确实使用了 fixed 模式：
- 日志输出: `Edge mode : fixed`
- 代码逻辑: `edge_mode == "fixed"` → `complex_edge_attr = torch.cat([bond_attr, knn_attr], dim=0)` ← 与 edge_index 顺序一致

---

## 6. 下一步

1. ~~**评估 Test AUC**~~: ✅ 已完成。Test AUC = **0.7060**，低于 legacy_bug 的 0.7244（Δ=-0.018）
2. **C1-Step 2 消融实验**: 重点调整大 batch 下的超参数
   - lr 线性缩放 (3e-4 → 4e-4 或更高)
   - warmup steps 增加 (8 epochs → 更多)
   - weight_decay 加入 (0 → 1e-5)
   - accumulate_grad_batches 增加（降低 effective BS）
3. **多 Seed 验证**: fixed vs legacy_bug 的差距是否显著

---

## 7. 文件索引

| 文件 | 路径 | 说明 |
|------|------|------|
| 训练日志 | `logs/01_fixed_baseline/train_fixed.log` | 9.1MB, 完整训练输出 |
| Metrics CSV | `results/01_fixed_baseline/metrics_fixed.csv` | 32行 (ep0-31) |
| Best Checkpoint | `results/01_fixed_baseline/checkpoints/pt-fixed-ep16-auc0.7145.ckpt` | 22MB |
| 训练曲线 | `results/01_fixed_baseline/fig1_training_dynamics.png` | Loss/AUC 曲线 |
| 家族性能 | `results/01_fixed_baseline/fig2_family_performance.png` | 7家族 AUC 分布 |
| 服务器脚本 | 服务器 `PathC/scripts/main_training_pt.py` | 路径已修正 |
| 服务器配置 | 服务器 `PathC/configs/base.yml` | 与 PathB 相同超参 |
