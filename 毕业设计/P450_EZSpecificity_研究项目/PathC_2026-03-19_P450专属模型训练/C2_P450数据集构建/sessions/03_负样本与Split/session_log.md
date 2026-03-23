# Phase 5: 负样本生成 + 4种Split Session Log

**日期**: 2026-03-24
**版本**: v6 (最终版)

---

## 1. 版本演进

| 版本 | 设计 | 问题 |
|------|------|------|
| v1 | 每种split独立生成负样本 | 4种split数据不同,不可比;对接54万对 |
| v2 | 全局负样本+连通分量all_split | 负样本比例仅7.28:1; all_split严重不平衡 |
| v3 | 全局负样本+软all_split | all_split酶86%/底物99%重叠,不是真正的unknown both |
| v4 | v3改为双向负样本 | all_split重叠率仍不可接受 |
| v5 | all_split独立data.csv | all_split和其他split不共享data.csv,不可接受 |
| **v6** | **共享data.csv + all_split对角线过滤** | **最终版** |

### 关键转折点
- 检查论文ESIBank数据发现: 4种split共享**完全相同的data.csv**(586,515行)
- 论文all_split实现了酶0%+底物0%重叠且不丢行(因25K酶+34K底物+多家族,互联度低)
- 我们P450数据(1.6K酶+2.1K底物+单家族)高度互联,无法不丢行实现0%重叠

## 2. 最终设计 (v6)

### 负样本生成
- **双向采样**(我们的改进): Direction A(固定底物换酶)×5 + Direction B(固定酶换底物)×5
- **全局采样**: 不限连通分量,从全部1,622酶/2,125底物中随机抽
- 排除已知正样本对,去重
- 合并成**一份共享data.csv**: 4,751正 + 47,503负 = 52,254行

### Split策略
**同一份data.csv(52,254行), 4种切法**:
- random_split: 行级随机分4折,**全部52,254行参与**
- reaction_split: 按底物分组分4折,**全部52,254行参与**,底物严格0%重叠
- enzyme_split: 按酶序列hash组分4折,**全部52,254行参与**,酶严格0%重叠
- all_split: 独立分配酶组和底物到4折,仅保留enzyme_fold==substrate_fold的行 → **18,159行参与**(34.8%),酶+底物严格0%重叠

### 列格式(匹配论文ESIBank)
`enzyme, reaction, label, ecnumber, difficulty, fake_ecnumber, structure_index`

## 3. 结果

### 全局统计
| 指标 | 值 |
|------|---|
| data.csv行数 | 52,254 (4种split共享) |
| 正样本 | 4,751 |
| 负样本 | 47,503 (A=23,755 + B=23,748) |
| 正负比 | 1:9.998 |
| 唯一对接对数 | 52,254 |

### 各Split泄露检查

| Split | fold文件行数 | 酶重叠 | 底物重叠 | 泄露检查 |
|-------|------------|--------|---------|---------|
| random_split | 52,254 | ~99%(允许) | ~99%(允许) | PASS |
| reaction_split | 52,254 | ~99%(允许) | **0%** | PASS |
| enzyme_split | 52,254 | **0%** | ~99%(允许) | PASS |
| all_split | **18,159** | **0%** | **0%** | PASS |

### 与论文ESIBank对比

| | 论文ESIBank | 我们P450 |
|---|---|---|
| 酶数量 | 25,225 | 1,622 |
| 底物数量 | 34,417 | 2,125 |
| data.csv | 586,515行 | 52,254行 |
| 负样本方向 | 单向(换酶) | **双向(换酶+换底物)** |
| 正负比 | 1:9.36 | 1:9.998 |
| all_split行数 | 586,515(100%) | 18,159(34.8%) |
| all_split酶重叠 | 0% | 0% |
| all_split底物重叠 | 0% | 0% |
| all_split丢行原因 | 不丢(多家族互联度低) | P450单家族高度互联 |

## 4. 输出文件

```
data/P450_Family/
├── Enzymes.csv              (1,622行, 共享)
├── Substrates.csv           (2,125行, 共享)
├── random_split/
│   ├── data.csv             (52,254行, 共享)
│   └── training/val/testing_datas_{0-3}.csv
├── reaction_split/
│   ├── data.csv             (52,254行, 共享)
│   └── ...
├── enzyme_split/
│   ├── data.csv             (52,254行, 共享)
│   └── ...
└── all_split/
    ├── data.csv             (52,254行, 共享, 和上面相同)
    └── training/val/testing_datas_{0-3}.csv  ← 只含18,159行(对角线子集)

data/combined/negatives/
├── biological.csv           (244条抑制剂)
└── split_stats.json         (完整统计)
```

## 5. Codex讨论记录 (8轮)

1. 负样本策略: 全局生成,先于split
2. Split策略: all_split连通分量→不可行
3. 代码实现v1: 独立负样本→发现与论文不一致
4. 对比论文ESIBank: 4种split共享同一data.csv
5. v2连通分量比例不足+不平衡 → v3全局+软all_split
6. v3/v4软all_split重叠86%/99% → 验证论文all_split是严格0%
7. v5独立data.csv → 用户不接受
8. v6: 共享data.csv + all_split对角线过滤 → 最终版

## 6. 下一步

Phase 6: 52,254个唯一对需要对接
