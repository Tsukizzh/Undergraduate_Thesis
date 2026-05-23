# Q2 聚类对象范围修正

日期：2026-05-21

## 结论

当前 Q2 已完成的 MMseqs2 聚类对象是 `Enzymes.csv` 全量 1622 条酶序列。这个结果可以作为全量 P450 酶序列聚类审计，但不能直接等同于 PathC/PathD baseline 实际训练 cache 的酶集合。

服务器只读核对发现，PathC 最佳 baseline 使用的是 `pt_cache_allfix_unified/random`，其 `train/val/test/index.pt` 源于 `splits/random/*_datas_0_pt.csv` 的一个过滤后子集。

## 已核对事实

全量表：

- `data/base_from_PathC/tables/Enzymes.csv`：1622 条酶序列。
- `data/base_from_PathC/tables/data.csv`：52254 条样本，覆盖 1622 个 enzyme 索引，范围 0-1621。
- `data/base_from_PathC/splits/random/data.csv`：52254 条样本，覆盖 1622 个 enzyme 索引，范围 0-1621。

baseline random PT split，fold0 原始 PT CSV：

- `training_datas_0_pt.csv`：26113 条样本，覆盖 1622 个 enzyme。
- `val_datas_0_pt.csv`：13054 条样本，覆盖 1616 个 enzyme。
- `testing_datas_0_pt.csv`：13056 条样本，覆盖 1617 个 enzyme。

baseline cache，fold0 过滤后实际进入 cache：

- train cache：22083 条样本，覆盖 1479 个 enzyme，范围 0-1576。
- val cache：11008 条样本，覆盖 1473 个 enzyme，范围 0-1576。
- test cache：10999 条样本，覆盖 1473 个 enzyme，范围 0-1576。
- train/val/test cache 合并后覆盖 1479 个 enzyme，缺失 143 个 Enzymes.csv 索引。
- cache 的 `enzymes/enzymes_index.pt` 中 `num_items=1577`，说明酶特征缓存本身只覆盖索引 0-1576，索引 1577-1621 未进入该酶特征缓存。

## 对 Q2 的影响

当前 1622 酶聚类没有白做，但后续正式生成 baseline-compatible 的序列相似度划分时，应该改成以 baseline 可训练样本宇宙为准：

1. 先追溯 `pt_cache_allfix_unified/random` 的过滤规则。
2. 构建 cache-valid 样本表，至少复现 fold0 cache 中的 train/val/test 样本集合。
3. 用 cache-valid enzyme 集合生成 Q2 的正式 cluster-held-out split。
4. 保留 1622 全量聚类作为上游审计，不直接用于训练结论。

## 仍需核对

- `test_eval.json` 记录 `n_samples=11000`，而当前 `test/index.pt` 统计为 10999 条。差 1 条样本，需要继续追溯评估脚本或 cache 加载逻辑。
- 143 个缺失 enzyme 的具体过滤原因：缺酶特征、缺结构、超过长度限制、缺对接 complex，还是样本过滤导致。
