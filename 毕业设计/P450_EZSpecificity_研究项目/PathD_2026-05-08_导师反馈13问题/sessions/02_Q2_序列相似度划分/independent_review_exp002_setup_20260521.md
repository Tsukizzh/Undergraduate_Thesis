# Q2 EXP002 actual-used 整理后独立复核记录

日期：2026-05-21

## 复核方式

优先尝试 Codex：

- `mcp__codex__codex`：超时，未返回审查结果。
- `collaborating-with-codex` 桥接：沙盒内网络被拦截；申请非沙盒后 5 分钟超时，未返回审查结果。

随后使用独立子代理做只读审查：

- agent id：`019e4a5f-5947-7e50-a5cd-7e77201ec186`
- 审查对象：本地下载的服务器产物副本、生成脚本、Q2 session 日志。

## 复核结论

独立审查未发现 actual-used 统计或 EXP002 当前状态的明显实质错误。

当前文件支持如下结论：

- EXP002 只是 actual-used 输入集和序列审计起点。
- MMseqs2 聚类未完成。
- 不能把 EXP002 当成正式 Q2 split。

## 核对通过的统计

| 项目 | 数值 |
|---|---:|
| actual-used enzyme | 1479 |
| missing enzyme | 143 |
| full catalog enzyme | 1622 |
| actual-used sample | 44090 |
| 正样本 | 3913 |
| 负样本 | 40177 |

按 split：

| split | 样本 | 正样本 | 负样本 |
|---|---:|---:|---:|
| train | 22083 | 1971 | 20112 |
| val | 11008 | 958 | 10050 |
| test | 10999 | 984 | 10015 |

审查还确认：

- `samples_actual_used.csv` 中 `cache_row` 在每个 split 内是 0 起连续编号。
- `sample_id_from_cache == pt_csv_row`。
- `enzyme_id_from_cache == enzyme_index`。
- `substrate_id_from_cache == substrate_index`。
- `split_membership_actual_used.csv` 与 `samples_actual_used.csv` 的核心列对应。
- 未见重复样本行。
- `exp002_manifest.json` 已明确 `mmseqs_status.status = not_completed`。
- MMseqs2 失败日志中存在 `Killed` 和 `Error: Prefilter step 0 died`。

## 复核指出的风险

1. `organize_q02_exp002_actual_used.py` 是一次性混合脚本，包含归档旧产物、物化 actual-used CSV/FASTA、运行 MMseqs2 三类动作，不适合以后直接重跑。
2. 当前本地审查基于导出 CSV 和 manifest；原始 `index.pt` 的复算已经在服务器端完成，但未把原始 `index.pt` 下载到本地复审目录。
3. `test_eval.json n_samples=11000` 与 `test/index.pt=10999` 的 1 条差异尚未关闭。

## 已采取/待采取措施

- 已在本记录中标注一次性混合脚本风险。
- 已准备在服务器 `scripts/` 下放置说明，提醒不要直接重跑该脚本。
- 正式 Q2 split 前仍需：
  - 重跑 EXP002 的 MMseqs2 id80/id60/id40/id30。
  - 解释 10999 vs 11000 的差异。
  - 确认 `enzyme_index` 等于 `Enzymes.csv` 行号。
  - 做重复序列同簇审计。
  - 做长度异常序列清单。
  - 做多随机种子 split 可行性模拟。
  - split 后审计同 enzyme、同精确序列、同 cluster 不跨 train/val/test。
