# Run01: best-checkpoint.ckpt

## 基本信息

| 项目 | 值 |
|------|-----|
| **运行编号** | Run01 |
| **检查点文件** | `best-checkpoint.ckpt` |
| **检查点路径** | `saved_model/model/run_0/models/best-checkpoint.ckpt` |
| **运行日期** | 2026-01-30 |
| **运行脚本** | `scripts/09_Step9_模型推理/step9_inference.py` |

## 模型说明

`best-checkpoint.ckpt` 是EZSpecificity训练过程中验证集性能最佳的最终模型。这是论文发布时使用的官方模型检查点。

## 输出文件

- `predictions.csv`: 517条预测结果
  - Dock Index, Enzyme Index, Substrate Index, Label
  - score (logit), logit, prob

## 核心指标

| 指标 | 值 |
|------|-----|
| AUC-ROC | 0.6636 |
| AUC-PR | 0.6360 |
| Accuracy | 0.5164 |
| Recall | 0.0996 |
| F1 Score | 0.1722 |

详细分析见 `data/10_Step10_结果分析/Run01_best-checkpoint/`
