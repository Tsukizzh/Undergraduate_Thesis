# Step 9: 模型推理 文件索引

## 脚本文件
| 文件 | 说明 |
|------|------|
| [step9_inference.py](../../scripts/09_Step9_模型推理/step9_inference.py) | Windows 安全推理脚本 |

## 数据文件
| 文件 | 说明 | 大小 |
|------|------|------|
| [predictions.csv](../../data/09_Step9_模型推理/predictions.csv) | 模型预测结果 | 517 条记录 |

## 输入依赖
| 文件 | 来源 |
|------|------|
| data.csv | Step 4 |
| enzyme_features.lmdb | Step 5 |
| reaction_features.lmdb | Step 6 |
| grover_fingerprint.lmdb | Step 7 |
| morgan_fingerprint.npy | Step 7 |
| structure_features.lmdb | Step 8 |
| high_quality_id.txt | Step 8 |

## 模型文件
| 文件 | 说明 |
|------|------|
| saved_model/model/run_0/models/best-checkpoint.ckpt | 预训练模型检查点 |
| saved_model/model/run_0/complete-full-random-all-0-complex.yml | 模型配置 |
