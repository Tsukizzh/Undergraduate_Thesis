# Chunk 03 配体三方深度验证 文件索引

## 输入文件
| 文件 | 路径 | 说明 |
|------|------|------|
| chunk_03.csv | `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_03.csv` | 57条待验证记录 (REC_0115-REC_0171) |

## 输出文件
| 文件 | 路径 | 说明 |
|------|------|------|
| verified_results.jsonl | `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_03_results/verified_results.jsonl` | 57条验证结果 (JSONL格式) |

## 会话文件
| 文件 | 路径 | 说明 |
|------|------|------|
| session_log.md | `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_03_session/session_log.md` | 操作日志 |
| file_index.md | `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_03_session/file_index.md` | 本文件 |

## 辅助脚本
| 文件 | 路径 | 说明 |
|------|------|------|
| batch_processor.py | `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_03_scripts/batch_processor.py` | 进度检查脚本 |

## 统计摘要
- **总记录**: 57条
- **任务02错误发现**: 5条 (REC_0115, REC_0118, REC_0127, REC_0132, REC_0154)
- **需人工审核**: 1条 (REC_0134)
- **分类分布**: SUBSTRATE(25) | INHIBITOR(19) | PRODUCT(3) | EXCLUDE(10)
