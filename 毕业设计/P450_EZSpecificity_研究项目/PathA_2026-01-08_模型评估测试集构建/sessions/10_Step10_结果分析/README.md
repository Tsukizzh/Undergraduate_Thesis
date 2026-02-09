# Sessions - 10 Step10 结果分析

本目录记录对不同Run推理结果的分析会话记录。

## 目录结构

```
10_Step10_结果分析/
├── Run01_best-checkpoint/      # Run01推理结果的分析会话
│   ├── session_log.md           # 会话日志（含论文背景、数据重叠查证等）
│   └── file_index.md            # 文件索引
├── Run02_best-checkpoint-v1/   # Run02分析会话
├── Run03_best-checkpoint-v2/   # Run03分析会话
├── Run04_best-checkpoint-v3/   # Run04分析会话
├── Run05_best-checkpoint-v4/   # Run05分析会话
├── 推理训练微调区别详解.md      # 技术文档：推理/训练/微调的区别
└── README.md                    # 本文件
```

## Run对应关系

| Run | 分析输入 | 分析输出 |
|-----|---------|---------|
| Run01 | data/09_Step9_模型推理/Run01_best-checkpoint/predictions.csv | data/10_Step10_结果分析/Run01_best-checkpoint/ |
| Run02 | data/09_Step9_模型推理/Run02_best-checkpoint-v1/predictions.csv | data/10_Step10_结果分析/Run02_best-checkpoint-v1/ |
| Run03 | data/09_Step9_模型推理/Run03_best-checkpoint-v2/predictions.csv | data/10_Step10_结果分析/Run03_best-checkpoint-v2/ |
| Run04 | data/09_Step9_模型推理/Run04_best-checkpoint-v3/predictions.csv | data/10_Step10_结果分析/Run04_best-checkpoint-v3/ |
| Run05 | data/09_Step9_模型推理/Run05_best-checkpoint-v4/predictions.csv | data/10_Step10_结果分析/Run05_best-checkpoint-v4/ |

## 执行脚本

- Run01: `scripts/10_Step10_结果分析/step10_analysis.py`（原始脚本）
- Run02-05: `scripts/10_Step10_结果分析/step10_batch_full_analysis.py`（批量分析）
- 多模型对比: `scripts/10_Step10_结果分析/step10_multi_run_analysis.py`

## 多模型对比汇总

所有Run的性能对比汇总见：`data/10_Step10_结果分析/多模型推理结果对比.md`

## 重要文档

### Run01/session_log.md
最详细的分析文档，包含：
- EZSpecificity论文背景（ESIBank、四种测试场景、负样本生成方式）
- P450测试集与论文的差异
- 数据重叠查证（酶0%重叠、底物10.8%重叠）
- 测试结果详细解读
- 生成图表说明
- 根因分析

### 推理训练微调区别详解.md
技术文档，包含：
- 推理/训练/微调的定义和区别
- 完整流程对比图
- 输入文件对比表
- EZSpecificity项目中的具体实现
- 常见问题解答
