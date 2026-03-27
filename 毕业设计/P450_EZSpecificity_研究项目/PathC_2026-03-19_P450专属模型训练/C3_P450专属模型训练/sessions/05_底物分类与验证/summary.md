# Step 5: 底物分类与多轮验证

> **日期**: 2026-03-27 ~ 2026-03-28
> **状态**: 分类完成 + 验证完成，~248 个错误已发现，修正待执行

## 做了什么

对 P450 数据集中 2,125 个底物进行化学类别分类，并通过 4 种工具多轮验证。

## 为什么做

导师方向：给定一个 P450 酶，预测它催化哪类底物（萜类？黄酮？生物碱？）。前提是底物分类必须准确。

## 分类方法

1. **NPClassifier**（神经网络）→ 2,125 全部分类为 15 类
2. **SMARTS 结构规则** → 7 类共 936 个确认
3. **ClassyFire**（规则分类器）→ 328 个有结果
4. **11 个 Opus 文献 agent** → ~400 个逐一查 PubChem/ChEBI（~500 次 web 调用）

## 15 个最终类别

| 类别 | 中文 | 数量 |
|------|------|------|
| Alkaloid | 生物碱 | 320 |
| Amino_acid | 氨基酸 | 292 |
| Fatty_acid | 脂肪酸 | 287 |
| Phenylpropanoid | 苯丙素 | 222 |
| Steroid | 甾体 | 211 |
| Diterpenoid | 二萜 | 144 |
| Unclassified | 未分类 | 104 |
| Triterpenoid | 三萜 | 97 |
| Polyketide | 聚酮 | 93 |
| Sesquiterpenoid | 倍半萜 | 93 |
| Monoterpenoid | 单萜 | 79 |
| Terpenoid_other | 其他萜类 | 67 |
| Flavonoid | 黄酮 | 53 |
| Macrolide | 大环内酯 | 38 |
| Coumarin | 香豆素 | 25 |

## 验证结果

| 类别 | 抽检正确率 | 主要问题 |
|------|-----------|---------|
| 甾体/二萜/三萜/单萜 | **100%** | 无 |
| 倍半萜 | 95% | 1 个 seco-甾体混入 |
| 其他萜类 | 85% | 3 个合成物 |
| 聚酮 | ~90% | 1 个 PAH |
| **苯丙素** | **64%** | **67 个 pathway fallback 错误**（PAH/卤代苯/合成物） |
| 脂肪酸 | ~52% | ~70 个烷烃/卤甲烷 |
| 生物碱 | ~71% | 合成药/无氮化合物 |
| 氨基酸 | ~74% | 环肽/合成物 |
| 未分类 | — | 39 个应重分类 |

**总计已发现 ~248 个明确错误，修正后预计可达 ~95%+。**

## 输出文件

- `data/05_底物分类/substrate_classifications_final.csv` — 当前分类结果
- `data/05_底物分类/npclassifier_cache.json` — NPClassifier API 缓存
- `data/05_底物分类/classyfire_cache.json` — ClassyFire API 缓存
- `scripts/05_底物分类/classify_substrates.py` — NPClassifier 批量分类
- `scripts/05_底物分类/cross_validate_smarts.py` — SMARTS 验证
- `scripts/05_底物分类/classify_classyfire_parallel.py` — ClassyFire 并行查询

## 下一步

与 Codex 多轮讨论修正方案 → 执行修正 → 达到 99%+ 准确率
