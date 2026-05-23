# Q2 EXP002 actual-used 数据集整理记录

日期：2026-05-21

## 目录调整

服务器 PathD 内已将 Q2 数据层整理为：

- `data/q02_sequence_similarity_split/exp001_full_catalog_1622/`
  - 保存原先 Q2 根目录下的全量 1622 酶聚类审计产物。
- `data/actual_used_baseline/`
  - 保存 PathC 最佳 baseline 实际 cache 使用的数据集，作为后续 Q1/Q2/Q4 等实验的公共数据入口。
- `data/q02_sequence_similarity_split/exp002_actual_used_cache_valid/`
  - 保存 Q2 EXP002 的 actual-used FASTA、序列审计和后续聚类/划分产物。
- `experiments/q02_sequence_similarity_split/EXP002_actual_used_cache_valid/`
  - 保存 EXP002 实验说明和后续脚本/结果。

## actual_used_baseline 内容

关键文件：

- `tables/Enzymes_actual_used.csv`
- `tables/Substrates_actual_used.csv`
- `tables/samples_actual_used.csv`
- `tables/split_membership_actual_used.csv`
- `tables/missing_enzymes_from_full_catalog.csv`
- `tables/cache_filter_audit.csv`
- `manifests/actual_used_dataset_manifest.json`
- `manifests/source_hashes.tsv`
- `reports/actual_used_dataset_audit.md`

统计：

| 项目 | 数值 |
|---|---:|
| 全量酶 | 1622 |
| baseline cache 实际酶 | 1479 |
| 未进入实际 cache 的酶 | 143 |
| 实际样本 | 44090 |
| 正样本 | 3913 |
| 负样本 | 40177 |
| 实际底物 | 2111 |

按 split：

| split | 样本 | 正样本 | 负样本 |
|---|---:|---:|---:|
| train | 22083 | 1971 | 20112 |
| val | 11008 | 958 | 10050 |
| test | 10999 | 984 | 10015 |

cache 对齐检查：

- train/val/test 的 `enzyme_mismatch=0`
- train/val/test 的 `substrate_mismatch=0`

## EXP002 当前状态

已完成：

- 生成 `enzymes_actual_used.fasta`
- 生成 `enzyme_sequence_audit.csv`
- 生成 `duplicate_sequences.csv`
- 写入 `sequence_stats.json`
- 写入 EXP002 manifest 和 validation summary

序列审计：

| 项目 | 数值 |
|---|---:|
| 酶数 | 1479 |
| 空序列 | 0 |
| 非法序列 | 0 |
| 唯一精确序列 | 1469 |
| 重复序列组 | 10 |
| 位于重复序列组的酶 | 20 |
| 最短 / 最长序列 | 60 / 853 |
| 平均 / 中位长度 | 496.86 / 506 |

未完成：

- MMseqs2 聚类。曾尝试 EXP002 的 id80 `easy-cluster`，但进程在 prefilter 阶段被系统终止，失败日志保存在：
  `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q02_sequence_similarity_split/exp002_actual_used_cache_valid/clusters/mmseqs/id80/mmseqs_easy_cluster.log`

## 完整性检查结果

只读校验通过：

- `1479 + 143 = 1622`
- `samples_actual_used.csv` 为 44090 行
- 样本涉及 enzyme 数等于 `Enzymes_actual_used.csv` 行数
- 样本涉及 substrate 数等于 `Substrates_actual_used.csv` 行数
- manifest 与 validation summary 的样本数一致
- cache 与 PT CSV 的 enzyme/substrate 对齐 mismatch 均为 0

## 下一步

1. 让独立 Codex 只读审查目录整理、数据还原和验证逻辑。
2. 重新处理 EXP002 的 MMseqs2 聚类失败，优先尝试单独运行更保守的命令或拆分参数，而不是直接进入 split。
3. MMseqs2 聚类完成后，再做 id60/id40 的 cluster-held-out split 可行性模拟。
