# Q2 EXP002 actual-used split 独立复核记录

日期：2026-05-21

## 复核方式

优先尝试 Codex MCP：

- `mcp__codex__codex`：120 秒超时，未返回结果。

随后使用独立子代理做只读审查：

- agent id：`019e4a93-f9f3-7f92-b973-26002a5f2446`
- 审查对象：最新 EXP002 manifests、reports、小型 split summary/validation 审查包、生成脚本、验证脚本和 session log。

## 复核结论

独立复核未发现阻塞性错误。

复核确认：

- `session_log.md` 已明确：EXP002 actual-used 是当前主线；EXP001 full-catalog 是归档。
- manifests/reports 已明确：当前 EXP002 cluster 来源为 `filtered_from_EXP001_full_catalog_1622_mmseqs_clusters`，不是 EXP002 actual-used MMseqs2 原地重跑。
- 四套候选 split 汇总一致：1479 enzyme、44090 samples、3913 positive、40177 negative。
- 四个阈值的 cluster/enzyme/重复序列跨 split 泄漏均为 0。
- id60 作为主候选、id40 作为严格对照的判断合理。

## 复核指出的问题

1. validation 脚本原先用 `set` 比较 sample key，理论上可能漏掉完全重复 key 的计数变化。
2. duplicate audit 原先对缺失 enzyme index 会静默跳过。
3. cluster leak 检查虽然通过，但它是按 cluster 赋值后的必然结果之一，仍需配合 sample key 和 cluster 覆盖检查。
4. 本地审查包没有包含所有大体量 `samples_with_*` 文件，因此复核方无法在本地逐行重算全部样本。

## 已修复/补强

已更新并重跑服务器脚本：

`/root/autodl-tmp/EZSpecificity/PathD/P450/scripts/validate_exp002_actual_used_splits_20260521.py`

补强内容：

- 使用多重集合（Counter）核对 sample key，而不是普通 set。
- 硬校验 cluster enzyme set 与 split sample enzyme set 一致。
- 硬校验 duplicate sequence group 中所有 enzyme 都存在于 actual-used split。

补强后验证结果：

| 阈值 | passes | sample_key_multiset_match | cluster_enzyme_set_match | duplicate_missing |
|---|---|---|---|---|
| id80 | true | true | true | [] |
| id60 | true | true | true | [] |
| id40 | true | true | true | [] |
| id30 | true | true | true | [] |

## 额外澄清

`test_eval.json n_samples=11000` 与 `test/index.pt=10999` 的差异已追溯。baseline 日志显示 DataModule 真实 test dataset 为 10999 samples；训练结束后的 4 卡 DDP auto-test 输出 11000 samples，原因是分布式测试采样会把 10999 补齐到可被 4 个进程整除的 11000。当前 actual-used 数据集仍以 cache index 的 10999 test samples 为准。

证据：

- baseline `nohub.out` 记录 `[DataModule] Test: 10999 samples`
- 同一日志稍后记录 `Samples      : 11000 (pos=984, neg=10016)`
- `10999` 在 4 卡 DDP 下会补齐到 `ceil(10999/4)*4 = 11000`

## 当前结论

EXP002 actual-used 的候选 cluster-held-out split 可以视为完成到“可选阈值并准备重建 pt_cache”的阶段。下一步不是继续 Q2 split 生成，而是选择主阈值，建议：

- id60：主方案。
- id40：严格对照。
