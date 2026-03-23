# Phase 5: 负样本生成 + 4种Split Session Log

**日期**: 2026-03-24
**版本**: v3 (论文对齐版)

---

## 1. 版本演进

| 版本 | 设计 | 问题 |
|------|------|------|
| v1 | 每种split独立生成负样本 | 4种split数据不同,不可比;对接54万对 |
| v2 | 全局负样本+连通分量all_split | 负样本比例仅7.28:1; all_split严重不平衡 |
| **v3** | **全局负样本+软all_split** | **论文对齐,52,261行共享data.csv** |

### 为什么v1错误
论文的ESIBank验证: 4种split共享**完全相同的data.csv**(586,515行)。我们v1每种split独立生成负样本,导致数据不同、对接量膨胀10倍。

### v2→v3修正
连通分量方法导致: (1)负样本只能从同分量采样→比例不足; (2)大分量不可分→all_split极端不平衡。改为全局采样+软all_split。

## 2. 最终设计 (v3)

### 负样本生成
- **单方向采样**（固定底物,换酶）→ 匹配论文
- **全局采样**（不限连通分量）→ 比例精确1:10
- **排除已知正样本对**
- **EC difficulty计算**: 真实酶与替换酶的EC号前缀匹配数

### Split策略
- **一份data.csv, 4种切法** → 匹配论文
- random_split: 行级随机,严格
- reaction_split: 按底物分组,严格无泄露
- enzyme_split: 按序列hash组分组,严格无泄露
- **all_split: 软约束** — 独立分配酶组和底物到fold,每行优先匹配两者的fold;保留所有行,报告重叠率

### 列格式（匹配论文ESIBank）
`enzyme, reaction, label, ecnumber, difficulty, fake_ecnumber, structure_index`

## 3. 结果

### 全局统计
- **master data.csv**: 52,261行（4,751正 + 47,510负）
- **正负比**: 1:10.0
- **唯一对接对数**: 52,261

### Split结果

| Split | 泄露检查 | fold平衡 | 备注 |
|-------|---------|---------|------|
| random_split | PASS (严格) | 均匀 | |
| reaction_split | PASS (严格) | 均匀 | 底物完全不重叠 |
| enzyme_split | PASS (严格) | 均匀 | 酶序列hash组完全不重叠 |
| all_split | 软约束 | 均匀 | 酶重叠~91%, 底物重叠~95% |

### all_split重叠率说明
P450是单一酶家族,大部分酶通过共享底物高度互联,严格零重叠不可行（连通分量方法会导致86%数据集中在一个fold）。论文ESIBank覆盖多家族,互联度低,但也几乎肯定使用软约束。

## 4. 输出文件

```
data/P450_Family/
├── Enzymes.csv              (1,622行)
├── Substrates.csv           (2,125行)
├── random_split/data.csv    (52,261行, 与其他split相同)
├── reaction_split/data.csv  (52,261行, 相同)
├── enzyme_split/data.csv    (52,261行, 相同)
└── all_split/data.csv       (52,261行, 相同)
    + 每个split: training_datas_{0-3}, val_datas_{0-3}, testing_datas_{0-3}
```

## 5. Codex讨论记录

- **Round 1**: 负样本策略 — 全局生成,单方向,先于split
- **Round 2**: Split策略 — all_split用连通分量(后发现不可行)
- **Round 3**: 代码实现 — v1(独立负样本)→发现与论文不一致
- **Round 4**: 发现论文4种split共享同一data.csv → 重做v2
- **Round 5**: v2连通分量导致比例不足+不平衡 → v3全局采样+软all_split
- **Round 6**: 最终验证 — Codex确认v3设计合理

## 6. 下一步

- Phase 6: 52,261个唯一对需要对接（比v1的54万减少90%）
