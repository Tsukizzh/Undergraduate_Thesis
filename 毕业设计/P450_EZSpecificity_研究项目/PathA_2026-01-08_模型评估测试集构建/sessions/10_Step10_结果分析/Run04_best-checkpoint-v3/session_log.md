# Run04 分析会话记录

## 基本信息
- **输入**: Run04_best-checkpoint-v3的推理结果
- **运行日期**: 2026-02-01
- **执行脚本**: step10_batch_full_analysis.py

## 核心指标

| 指标 | 值 |
|------|-----|
| AUC-ROC | 0.6624 |
| AUC-PR | 0.6440 |
| Accuracy | 0.5106 |
| Recall | 0.0805 |
| F1 Score | 0.1443 |

## 输出文件

- `data/10_Step10_结果分析/Run04_best-checkpoint-v3/`
  - metrics_summary.csv
  - error_analysis.csv
  - analysis_report.txt
  - figures/main_analysis_figure.png
