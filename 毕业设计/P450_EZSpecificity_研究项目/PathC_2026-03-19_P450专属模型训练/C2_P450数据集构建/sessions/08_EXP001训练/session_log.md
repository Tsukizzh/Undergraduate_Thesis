# Session Log: C2 Phase 8 — EXP001 P450 Random Split 基线训练

> **日期**: 2026-03-26
> **状态**: ✅ 已完成

---

## 一、实验目标

用我们从头构建的 P450 数据集（C2 Phase 1-7 产出），在 random_split 上训练 EZSpecificity 模型，验证数据集质量。

## 二、实验配置

| 项目 | 值 |
|------|-----|
| 服务器 | Cloud-2 (4×RTX4090, 360GB RAM, 28核/64vCPU) |
| Split | random_split fold0 |
| 训练集 | 23,710 样本 |
| 验证集 | 11,878 样本 |
| 测试集 | 11,823 样本（1,101正 / 10,722负） |
| Batch size | 56/GPU, effective=224 (4 GPU DDP) |
| num-workers | 6/GPU (24 total) |
| 预加载 | --preload (~149GB RAM) |
| edge-mode | fixed |
| max-epochs | 200 (EarlyStopping patience=15) |
| 监控指标 | auc/val |
| save_top_k | 5 |

## 三、训练结果

| 指标 | 值 |
|------|-----|
| 总 epochs | 89 (early stopped at ep88, best=ep73) |
| **Val AUC** | **0.7544** |
| **Test AUC** | **0.7730** |
| Test AUPR | 0.3620 |
| Test 样本数 | 11,823 (1,101正 / 10,722负) |
| 训练速度 | ~1.8 it/s, ~2 min/epoch |
| GPU利用率 | 92-100% |
| 总耗时 | **49 分钟** |

## 四、与基线对比

| 模型 | 数据 | Split | Test AUC | 说明 |
|------|------|-------|----------|------|
| 论文 checkpoint | ESIBank 全部 | random | — | 论文未公布 random test AUC |
| 论文 checkpoint | ESIBank P450 子集 | internal | **0.638** | Path B Step 7 E7 内部基准 |
| 我们从头训练 | ESIBank all_split | all_split | **0.7244** | Path B Step 10 Cloud-2 DDP |
| **我们从头训练** | **P450 数据集** | **random** | **0.7730** | **本实验 (EXP001)** |

**关键结论**: P450 专属数据集 Test AUC=0.7730，比 ESIBank 内部 P450 评估的 0.638 提升了 **+0.135**，确认数据集构建成功。

## 五、训练模板改动

main_training_pt.py 经 Codex 4 轮 review 重写：
- 输出路径扁平化（不再嵌套 results/run_tag/）
- max-epochs 50→200, save_top_k 3→5
- 训练后自动 test (trainer.test())
- DDP test_epoch_end 补丁（之前只有 validation 有 all_gather）
- 训练+测试完成后自动关机
- SRC_DIR fallback: experiment-local src/ → global src/

## 六、实验管理

每个实验为自包含目录：
```
P450/experiments/EXP001_random_baseline/
├── src/Models, Datasets, utils.py    # 完整代码副本
├── scripts/main_training_pt.py, pt_dataset.py, launch.sh
├── configs/config.yml
├── logs/train.log + TensorBoard
└── results/checkpoints/ + metrics.csv + test_eval.json
```

## 七、test_eval.json

```json
{
  "test_auc_roc": 0.7730,
  "test_aupr": 0.3620,
  "n_samples": 11823,
  "n_positive": 1101,
  "n_negative": 10722
}
```
