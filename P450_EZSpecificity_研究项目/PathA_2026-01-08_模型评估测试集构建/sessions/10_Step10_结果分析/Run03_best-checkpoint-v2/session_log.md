# Run03 分析会话记录

## 基本信息
- **输入**: Run03_best-checkpoint-v2的推理结果
- **运行日期**: 2026-02-01
- **执行脚本**: step10_batch_full_analysis.py

## 核心指标

| 指标 | 值 |
|------|-----|
| AUC-ROC | **0.6674** (最佳) |
| AUC-PR | 0.6336 |
| Accuracy | 0.5145 |
| Recall | 0.0958 |
| F1 Score | 0.1675 |

## 输出文件

- `data/10_Step10_结果分析/Run03_best-checkpoint-v2/`
  - metrics_summary.csv
  - error_analysis.csv
  - analysis_report.txt
  - figures/main_analysis_figure.png
