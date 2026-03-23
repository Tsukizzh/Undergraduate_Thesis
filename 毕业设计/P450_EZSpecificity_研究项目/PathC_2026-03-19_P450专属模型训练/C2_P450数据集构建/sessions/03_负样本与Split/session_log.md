# Phase 5: 负样本生成 + 4种Split Session Log

**日期**: 2026-03-24
**目标**: 为4,751条去重后的正样本生成双向负样本，构建4种split × 4折CV

---

## 1. 策略设计（Codex 4轮讨论）

### Round 1: 负样本生成策略
- **先split正样本，再在每个partition内独立生成负样本** — 防止跨fold泄露
- 双向采样: Direction A（固定底物换酶）5个 + Direction B（固定酶换底物）5个
- 正负比约1:10
- 排除已知正样本对

### Round 2: Split划分 + 边界情况
- **all_split**: 不用朴素25%对角线，用贪心优化最大化保留边
  - 先按交互数贪心分配酶序列hash组到4折
  - 再将底物分配到其交互最多的酶折，最大化同折保留
- **fold平衡**: 按交互数降序贪心bin-packing
- **序列hash分组**: 11个重复序列组必须在同一折
- **负样本域**: partition-local（只用该partition内的实体采样）
- **抑制剂**: 单独存放，不混入主pipeline
- **Dock Index=-1, positive_reactions=-1** 作占位

### Round 3: 代码实现
- Codex给出~500行完整原型
- 我基于原型重写，优化代码结构

### Round 4: 结果验证
- Codex确认所有结果正确，无重大问题
- shortfall是结构性的（partition内底物有限），可接受

## 2. 输入数据

| 项目 | 数量 |
|------|------|
| 正样本交互 | 4,751 |
| 唯一酶 | 1,622 |
| 唯一化合物 | 2,125 |
| 序列hash重复组 | 11（含24个酶） |
| 单交互酶 | 818 (50.4%) |
| 单交互化合物 | 1,352 (63.6%) |

## 3. Split结果

### 3.1 random_split
- 全部4,751正样本参与
- 正负比精确1:10（5A+5B全满）
- 泄露检查: PASS

### 3.2 reaction_split（按底物分组）
- 全部4,751正样本参与
- Direction B有少量shortfall（max -57, 底物partition内实体有限）
- 泄露检查: PASS（底物在train/val/test之间无重叠）

### 3.3 enzyme_split（按酶序列hash组分组）
- 全部4,751正样本参与
- Direction B有较大shortfall（max -1,238, 酶partition内底物多样性受限）
- 泄露检查: PASS（酶序列hash组在train/val/test之间无重叠）

### 3.4 all_split（酶+底物同时分组）
- **保留3,003条（63.2%）**，丢弃1,748条（酶/底物不在同一fold）
- 优化效果: 朴素方法仅保留~25%，优化后提升至63.2%
- Per-fold正样本: [759, 687, 788, 769]，平衡度=1.15
- 泄露检查: PASS（酶和底物同时在train/val/test之间无重叠）

### 汇总

| Split | 正样本 | 保留率 | 负样本(train/fold) | 泄露 |
|-------|--------|--------|-------------------|------|
| random_split | 4,751 | 100% | ~23,750 | PASS |
| reaction_split | 4,751 | 100% | ~23,750 | PASS |
| enzyme_split | 4,751 | 100% | ~23,750 | PASS |
| all_split | 3,003 | 63.2% | ~14,750 | PASS |

## 4. 输出文件

```
data/P450_Family/
├── Enzymes.csv              (1,622行, 共享)
├── Substrates.csv           (2,125行, 共享)
├── random_split/
│   ├── data.csv             (176,524行, fold-union)
│   ├── training_datas_{0-3}.csv
│   ├── val_datas_{0-3}.csv
│   └── testing_datas_{0-3}.csv
├── reaction_split/          (177,469行 data.csv)
├── enzyme_split/            (176,615行 data.csv)
└── all_split/               (111,481行 data.csv)

data/combined/negatives/
└── split_stats.json         (完整统计)
```

CSV格式: `Substrate Index, Enzyme Index, Label, Dock Index, positive_reactions`

**data.csv说明**: 各fold产生的所有行的去重并集（同一正样本在不同fold中以不同角色出现）。fold文件是权威的train/val/test输入。

## 5. 关键设计决策记录

1. **负样本在split后生成** — 防止train负样本使用test酶/底物
2. **all_split优化** — 基于贪心bin-packing + 底物跟随酶折，保留率从25%提升至63.2%
3. **sequence_hash分组** — 11个重复序列组在enzyme_split和all_split中始终在同一fold
4. **确定性** — 主种子20260324 + SHA256派生子种子，完全可复现
5. **抑制剂独立** — Path A的245条抑制剂不混入，后续可选混合实验

## 6. 下一步

- Phase 6: 结构获取与对接（需GPU服务器）
- Phase 7: 特征生成（ESM + GROVER + Morgan + 结构特征）
- Phase 8: 训练 + 评估
