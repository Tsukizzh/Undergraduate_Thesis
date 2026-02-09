# Sessions - 09 Step9 模型推理

本目录记录使用不同检查点进行推理的会话记录。

## 目录结构

```
09_Step9_模型推理/
├── Run01_best-checkpoint/      # 使用best-checkpoint.ckpt的推理会话
│   ├── session_log.md           # 会话日志
│   └── file_index.md            # 文件索引
├── Run02_best-checkpoint-v1/   # 使用best-checkpoint-v1.ckpt的推理会话
├── Run03_best-checkpoint-v2/   # 使用best-checkpoint-v2.ckpt的推理会话
├── Run04_best-checkpoint-v3/   # 使用best-checkpoint-v3.ckpt的推理会话
└── Run05_best-checkpoint-v4/   # 使用best-checkpoint-v4.ckpt的推理会话
```

## 检查点说明

| Run | 检查点文件 | 说明 |
|-----|-----------|------|
| Run01 | best-checkpoint.ckpt | 训练最终保存的最佳模型（默认） |
| Run02 | best-checkpoint-v1.ckpt | 训练过程第1次刷新时保存 |
| Run03 | best-checkpoint-v2.ckpt | 第2次刷新时保存 |
| Run04 | best-checkpoint-v3.ckpt | 第3次刷新时保存 |
| Run05 | best-checkpoint-v4.ckpt | 第4次刷新时保存 |

## 执行脚本

- Run01: `scripts/09_Step9_模型推理/step9_inference.py`（原始脚本）
- Run02-05: `scripts/09_Step9_模型推理/step9_multi_checkpoint_inference.py --checkpoint <name>`

## 对应输出

各Run的推理结果位于：`data/09_Step9_模型推理/<RunID>/predictions.csv`
