# Session 10 — EXP004 论文基线外部评估

**日期**: 2026-04-15
**状态**: ✅ 完成
**核心数字**: 论文模型在过滤后 P450 测试集 AUC = **0.559**，我们 EXP001 = **0.921**，**差距 +0.362 AUC / +0.540 AUPR**。

---

## 一、为什么做这个实验

### 我们已有的数字

| 模型 | Test AUC | 说明 |
|---|---|---|
| EXP001_allfix_unified (bare 28 维) | **0.9320** | 我们在自建 P450 数据集上从头训练，allfix 双 bug 修复后的 baseline |

### 缺失的对照

**我们从未在同一测试集上评估过论文预训练模型的表现。** 论文自身报告的 Unknown enzyme+substrate 场景 AUC=0.7198，但那是论文的 BRENDA 测试集，不是 P450 专属。这导致我们无法回答一个关键问题：

> 如果直接拿论文 ckpt 跑我们的 P450 测试集，能拿多少 AUC？我们 0.932 的优势是"P450 专属数据集 + bug 修复"带来的真实提升，还是仅仅是测试分布差异？

### 公平性约束

直接把论文 ckpt 放到我们的 test 上会有**数据泄漏**：论文训练集里的 389 个 ESIBank P450 酶中，有 356 个也出现在我们的 P450 数据集里（占 91.5%）。如果不剔除这 356 个酶对应的样本，论文模型会被"它训练时见过的酶"拉抬，得出的 AUC 没有意义。

**目标**：在**排除所有论文训练见过的 P450 酶**后的测试集上，公平对比论文 ckpt 和我们 EXP001 ckpt，看论文模型对"未见过的 P450"到底有多大泛化能力。

---

## 二、最终结果（一页看懂）

### 5 路推理

| 模型 | 测试集 | 边模式 | **Test AUC** | Test AUPR | 样本数 |
|---|---|---|---|---|---|
| **论文 ckpt** | 过滤后 (主结果) | legacy_bug | **0.5586** | 0.1004 | 7963 |
| 论文 ckpt | 过滤后 | fixed | 0.5596 | 0.1007 | 7963 |
| **论文 ckpt** | **未过滤 (对照)** | legacy_bug | **0.5860** | 0.1124 | 10999 |
| 我们 EXP001 | 过滤后 | legacy_bug | 0.9154 | 0.6194 | 7963 |
| **我们 EXP001** | **过滤后 (主对比)** | **fixed** | **0.9205** | **0.6403** | 7963 |

### 三层证据闭环

```
┌─────────────────────────────────────────────────────────────────┐
│ 1. 论文未过滤 0.586  →  论文过滤后 0.559   ΔAUC = +0.027       │
│    证明：pipeline 无 bug，filter 机制工作正常                    │
│    证明：论文对它"见过的酶"也只有微弱记忆优势                    │
├─────────────────────────────────────────────────────────────────┤
│ 2. 论文过滤后 0.559  vs  我们过滤后 0.921   ΔAUC = +0.362      │
│    主结果：同一架构、同一数据、leakage-controlled 下            │
│    我们对未见 P450 泛化能力远超论文                              │
├─────────────────────────────────────────────────────────────────┤
│ 3. 我们未过滤 0.932  →  我们过滤后 0.921   ΔAUC = −0.012       │
│    证明：我们的高 AUC 不是靠记忆 ESIBank 酶                      │
│    过滤前后自身几乎不变 → 泛化真实                               │
└─────────────────────────────────────────────────────────────────┘
```

### 核心论断

1. **论文模型对 P450 整体就弱**。即使在它训练过的 P450 酶上（未过滤 test 含 27.6% 训练酶），AUC 也只有 0.586，离它自身的 0.72 相差甚远。
2. **我们的 +0.362 AUC 优势是真实的**，不是测试分布差异造成的假象。
3. **为什么论文对训练酶也弱**：论文训练的是完整的 (enzyme, substrate, complex) 三元组。我们 test 里酶熟悉，但**配对的底物来自 5 个新数据源**（RCSB / ESIBank / P450Rdb / PlantP450DB / PCPD），**对接复合物是用 Uni-Dock 重跑的**——三元组对论文整体仍是新样本。这正好解释了我们建立 P450 专属数据集的动因。
4. **Edge mode 在 inference 时不敏感**。论文 ckpt legacy vs fixed 只差 0.001，我们差 0.005。边排序 bug 主要影响训练收敛，不影响 inference 数值。

---

## 三、方法

### 黑名单

- **来源**：PathA 阶段已整理的 `389个P450的PDB映射_完整版.csv`，即 ESIBank 训练集全部 P450 酶的 UniProt ID 列表。这是一个**完备集**（ESIBank 里 P450 家族总共就 389 个，不存在隐藏的 P450）
- **匹配方式**：UniProt ID 精确匹配到我们 `P450/data/Enzymes.csv`（1622 行，每行 row idx = enzyme_global_id）
- **结果**：**356 个** `enzyme_global_id` 命中（91.5% of 389 / 21.9% of 1622）
- **残余限制**：我们 1622 个酶中有 165 个用合成 ID (`ENZ_G*`)，来自 PlantP450DB 等没有 UniProt 的来源。由于我们没有 ESIBank 389 酶的序列，无法做 sequence-hash 兜底，最多可能漏 33 个隐藏同源酶。写进限制声明即可。

### 过滤缓存（非破坏性 overlay）

守住"只加不减"原则，**原始 `pt_cache_allfix_unified/random/` 一个字节都不修改**：

```
pt_cache_allfix_unified_paperfilter/random/   ← 新建目录
├── enzymes     → symlink
├── substrates  → symlink
├── manifest.pt → symlink
├── train       → symlink
├── val         → symlink
└── test/
    ├── samples → symlink
    └── index.pt ← 唯一的真·新文件（~140 KB）
```

新 `test/index.pt` 的生成方式：读原 index.pt，对 5 个字段（sample_ids / enzyme_ids / substrate_ids / graph_shards / graph_rows）**按同一个 boolean mask 切一刀**——即 `~torch.isin(enzyme_ids, blacklist)`。位置对齐完全保留，sample_id 保持原值（本来就非连续）。

**过滤结果**：

| | 原始 | 过滤后 |
|---|---|---|
| 样本数 | 10999 | **7963** （丢 27.6%） |
| 唯一酶 | 1473 | 1125 （丢 348 个）|
| 正样本 | 984 | 680 |
| 负样本 | 10015 | 7283 |

### 推理前的预检（防止静默失败）

两个"看起来没必要但其实必须做"的检查：

1. **论文 ckpt strict load 预检**（本地 CPU, torch 2.3.0 + PL 1.9.0 与论文训练环境一致）
   - 实例化 `SS(config)` + `torch.load(paper_ckpt)` + `load_state_dict(strict=True)`
   - 结果：76/76 keys matched，0 missing / 0 unexpected / 0 shape mismatches
   - **为什么必须做**：`main_training_pt.py` 的 `test_evaluate()` 内部用的是 `strict=False`，如果 key 错位会静默加载空模型，AUC 看起来可能 ~0.5 但其实模型根本没载入任何权重

2. **pt_cache 数据完整性穿透验证**（服务器 CPU）
   - **Enzyme 侧**：12 spot check + 5 深度穿透 → csv_len == bin_len，enzymes_index.pt 稀疏保留 0..1621 原始行号，与 allfix 修复后的正确布局一致
   - **Substrate 侧**：全量 2124/2124 匹配 `GetNumAtoms() == bin_atom_len`，key=8 缺失（正好是 `*[H]` 这个 allfix GROVER fix 的已知例外），**验证 allfix GROVER rekey 持续有效**
   - **新 index.pt**：5 数组等长 + 0 blacklist leak + 所有 enzyme_id/substrate_id 在 lookup 里
   - **原缓存零改动**：所有 overlay 条目是 symlink，仅 test/index.pt 是新文件

3. **Pipeline smoke test**（开 GPU 后，正式跑前）
   - 1 batch 前向（88 样本）
   - 检查：logits 全 finite + std=3.55（非常数）+ 前 10 个 (label, logit, tag) 合理分布
   - 通过才允许跑 full grid

### GPU 推理

- 硬件：1× RTX 5090 (32GB) 无卡扩展模式
- 配置：`--test-only --batch-size 88 --num-workers 6`
- 5 次 inference 顺序跑，每次 ~43-52 秒，总 ~4 分钟
- 所有 inference 确定性（`eval()` + `torch.inference_mode()`），重复跑 logits 位对位一致

---

## 四、结果解读

### Q1：为什么 paper 过滤后只有 0.559？

两种可能解释：

| 假设 | 如果是真 |
|---|---|
| **A**：论文模型本来就弱于 P450，这是真实表现 | EXP004 主结果可信 |
| **B**：我们 preprocessing 坑了论文 ckpt，输出不能用 | EXP004 主结果作废 |

### Q2：如何区分 A 和 B？

**未过滤 test 对照**：10999 样本，含 27.6% 论文训练过的酶。

- 如果 B 是真：pipeline 对 paper 全崩，未过滤 AUC 也应该 ~0.56
- 如果 A 是真：paper 对训练酶有记忆优势，未过滤 AUC 应该明显高于 0.56

实际结果：**未过滤 AUC = 0.586**，比过滤后高 **+0.027**。

这说明：
- **Pipeline 正常**——paper 对见过的酶确实有微弱的记忆优势
- **但优势很小**——0.027 的 delta 远小于我们 vs paper 的 0.362 gap
- 所以 paper 即使在熟悉的 P450 酶上，也只比随机好一点点

**结论**：假设 A 成立。EXP004 主结果 +0.362 AUC 完全夯实。

### Q3：论文自身测试集报告的 0.72 AUC 去哪了？

论文的 Unknown enzyme+substrate AUC=0.72 是在 BRENDA 8 家族混合测试集上测的，P450 只是 8 个家族之一，且测试样本的负样本构造方式与论文训练一致。

我们的 test 对论文模型来说是**域外的**：
- 酶：356 个重合酶之外的 1117 个是论文没见过的
- 底物：来自 5 个新 P450 数据源，与论文 BRENDA P450 底物集合只 ~10% 重合
- 复合物：Uni-Dock 重新生成的（论文用 AutoDock Vina）
- 负样本生成方式：我们是双向生成（5A + 5B），论文是单向

这些 domain shift 叠加起来，让论文模型在我们 test 上的表现从 0.72 降到 0.56-0.59 区间。

### Q4：为什么我们的模型能到 0.92？

- 训练数据直接就是 P450 专属（5 库 merge, 4751 正样本, 47510 对）
- ESM + GROVER 双 LMDB bug 修复后的 allfix_unified baseline
- 我们的 filter 前 AUC 0.932 → filter 后 0.921，几乎不变，说明泛化是真的

---

## 五、限制声明（未来报告必须写的）

1. **黑名单仅 UniProt 匹配**：未做序列 hash 兜底。165 个合成 ID (`ENZ_G*`) 酶里可能隐藏和 33 个"不在我们数据的 ESIBank UniProt"对应的同一蛋白，最多漏 33 个。
2. **仅过滤 enzyme，未过滤 substrate**：论文 random_split 与我们 test 在底物层面仍有重合。
3. **过滤率 27.6%**：过滤后 test 应标注为 "ESIBank-P450-overlap-filtered subset"，不是"完整 external test"。
4. **单次 run, 无 seed variation**：推理确定性，但不同随机种子下模型训练结果可能略有差异；此处只能说"这一个 ckpt 对这一个 test 的表现"。
5. **论文 ckpt epoch=106**：论文公开的是 `best-checkpoint.ckpt`（Val 最高 epoch），如果 epoch 选择不同可能数字略有差异，但不会改变数量级结论。

---

## 六、文件位置索引

### 服务器

- **实验目录**：`/root/autodl-tmp/EZSpecificity/PathC/P450/experiments/EXP004_paper_baseline_unified/`
- **过滤缓存**：`/root/autodl-tmp/EZSpecificity/PathC/P450/data/pt_cache_allfix_unified_paperfilter/random/`
- **论文 ckpt**：`EXP004_paper_baseline_unified/results/checkpoints/paper_best-checkpoint.ckpt` (22MB, md5 `f4d87ea08fc64b62700aadef8bf151cf`)

### 本地

- **实验目录完整归档**：`results/10_论文基线外部评估_EXP004_paper_baseline_unified/`（含 src / scripts / configs / results / logs）
- **session log**：本文件
- **结果 JSON**（5 个）：同在上述实验目录归档下的 `results/` 子目录

### 关键脚本（本地 `scratch/`）

- `build_paper_blacklist.py` — 389 vs 1622 交集
- `build_paperfilter_cache.py` — 服务器 overlay + 切 mask
- `verify_paperfilter_cache.py` — 缓存完整性验证
- `preflight_ckpt_load.py` — 本地 ckpt strict 预检
- `exp004_smoke_test.py` — 开 GPU 后正式跑前的 pipeline 健康检查

---

## 七、后续方向

1. **近期**（1-2 天）
   - per-subfamily / per-substrate-class AUC 分解，看 paper 在哪类 P450 彻底崩、我们在哪类能保持
   - 画 ROC + PR 曲线 + logit 分布对比图，毕设直接可用
2. **中期**（毕设驱动）
   - 把 EXP004 + allfix 三对比写进毕设 Results 章节，这是最强的 external validation
   - 若有时间跑 3 个 natural 变体做最大数据量对比
3. **长期**
   - Step 14 双尺度结构编码——allfix 结果下已不推荐（原子级 EGNN 已吃完结构信号）
   - C3 Stage 2 多标签底物类别预测
   - Path D 区域选择性预测（S3 反应 SMILES + S9 反应图片）
