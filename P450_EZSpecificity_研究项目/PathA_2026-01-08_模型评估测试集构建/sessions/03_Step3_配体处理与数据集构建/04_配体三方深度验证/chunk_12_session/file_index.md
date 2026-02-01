# Chunk 12 文件索引

## 输入文件
- `chunk_12.csv` - 原始输入数据（56条P450酶-配体记录）

## 输出文件
- `verified_results.jsonl` - 验证结果（**已完成56条**）

## 辅助脚本
- `batch_process_remaining.py` - 统计剩余未处理记录的辅助脚本

## 中间数据
- 无

## 处理状态
| 状态 | 数量 |
|------|------|
| 已完成 | 56条 |
| 待处理 | 0条 |
| Task02错误 | 10条 |
| 需人工复审 | 11条 |

## 最终分类分布
| 分类 | 数量 |
|------|------|
| SUBSTRATE | 21条 |
| INHIBITOR | 21条 |
| EXCLUDE | 9条 |
| PRODUCT | 5条 |
| **总计** | **56条** |

## 发现的Task02错误汇总
| Record ID | 原分类 | 更正为 | 配体 |
|-----------|--------|--------|------|
| REC_0642 | SUBSTRATE | EXCLUDE | PALMITIC ACID |
| REC_0647 | SUBSTRATE | INHIBITOR | LAURIC ACID |
| REC_0657 | PRODUCT | SUBSTRATE | Filipin product |
| REC_0659 | SUBSTRATE | PRODUCT | Oligomycin A |
| REC_0660 | SUBSTRATE | PRODUCT | Amphotericin B |
| REC_0662 | SUBSTRATE | INHIBITOR | 13-hydroxy fatty acid |
| REC_0663 | SUBSTRATE | INHIBITOR | 13-hydroxy fatty acid |
| REC_0678 | SUBSTRATE | PRODUCT | Nystatin |
| REC_0680 | PRODUCT | SUBSTRATE | Hydroxy-farnesene |
| REC_0682 | INHIBITOR | SUBSTRATE | Compactin/Mevastatin |

## 备注
- 错误率：17.9% (10/56)
- 三方协作（Claude + Gemini + Codex）全程严格执行
- 详细日志见 `session_log.md`
