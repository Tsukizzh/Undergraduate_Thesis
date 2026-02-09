# Run03 推理会话记录

## 基本信息
- **检查点**: best-checkpoint-v2.ckpt
- **运行日期**: 2026-02-01
- **执行脚本**: step9_multi_checkpoint_inference.py --checkpoint best-checkpoint-v2

## 运行命令

```bash
python step9_multi_checkpoint_inference.py --checkpoint best-checkpoint-v2
```

## 输出

- 预测结果: `data/09_Step9_模型推理/Run03_best-checkpoint-v2/predictions.csv`
- 517条预测，Logit范围 [-22.27, 2.55]

## 后续分析

分析结果见: `data/10_Step10_结果分析/Run03_best-checkpoint-v2/`
