# Path C：P450 专属模型训练 — 计划与进度

> **创建日期**: 2026-03-19
> **当前状态**: C2 Phase 8 EXP001 ✅ (Test AUC=0.7730) | C3 底物分类 ✅ v6 FINAL (352 review全量Agent验证, 1870确认88.0%/255 other 12.0%/0 review) → 模型训练待启动
> **前置条件**: Path B Step 1-10 全部完成，legacy_bug 基线 Test AUC=0.7244 > 论文 0.7198

---

## 一、目标

在 Path B 完成的系统性诊断基础上，通过**模型架构优化**和**数据集扩建**两条线，进一步提升 P450 酶特异性预测性能。

**核心指标**: AUC-ROC（测试集 8841 样本）

**基线对比**:
| 基线 | Val AUC | Test AUC | 说明 |
|------|---------|----------|------|
| 论文 all_split | — | 0.7198 | 4×GPU, ~256ep, bs=16/GPU |
| Step 9 本地 (fixed) | 0.7667 (ep14) | 未测 | 单卡, bs=14, accum=2 |
| Step 10 Cloud-2 (legacy_bug) | 0.722 (ep22) | **0.7244** (ep27) | 2×4090, bs=56, 32ep |

---

## 二、总体架构

```
Path C
├── C1_论文基线训练与参数调整 ✅
│   ├── AllSplit 从头训练 (ESIBank) ✅ → Test AUC=0.7244 (超论文0.7198)
│   ├── .pt 缓存管线 + Cloud-2 DDP ✅ → 7.56 it/s (解决LMDB thrashing)
│   ├── fixed 基线训练 ✅ → Test AUC=0.7060 (边修复未提升)
│   └── dropout 消融 ✅ → Val↑但Test未迁移，不纳入后续
│
├── C2_P450全面数据集构建 ✅ Phase 1-8 全部完成
│   ├── Phase 1-3: 5个数据库下载+标准化 → 4,751交互
│   ├── Phase 4-5: 合并去重+4种Split → 52,254行
│   ├── Phase 6-7: 对接+特征生成 → 47,510可用对
│   └── Phase 8: EXP001 ✅ random_split Test AUC=0.7730
│
└── C3_P450专属模型训练 🔄 模型训练待启动
    ├── Step 1: 底物分类 ✅ (NPClassifier 15类→7+1类多标签)
    ├── Step 2: 多轮验证+修正 ✅ (v5, 5优先级管线, 150抽检~89%, Codex 8轮)
    ├── Step 3: 352 review/other 全量Agent验证 ✅ (v6 FINAL, 20批Sonnet Agent文献搜索+Codex审核, 97升级+255确认OTHER, review→0)
    ├── Step 4: 按类别聚合预测 ⏳ (酶序列→ESM→MLP→7-sigmoid→BCEWithLogitsLoss)
    └── Step 5: 其他split训练 + 评估 ⏳
```

---

## 三、C1 模型架构优化

### C1-Step 1: fixed 基线训练

**目的**: 量化边排序修复（fixed vs legacy_bug）对 AUC 的贡献

**方法**:
- 与 Step 10 legacy_bug 完全相同的配置，仅改 `--edge-mode fixed`
- Cloud-2, 2×4090 DDP, bs=56, max_epochs=50, EarlyStopping patience=15

**预期**:
- fixed 应该 ≥ legacy_bug（0.7244），因为边修复纠正了错误的边特征对齐
- Step 9 本地 fixed 的 val AUC=0.7667，但 effective batch 不同，不可直接比较

**结果**: Test AUC=0.7060, ΔAUC = 0.7060 - 0.7244 = -0.018（边修复未提升测试集性能）

**时间**: ~1 天（8-10 小时训练 + 评估）

---

### C1-Step 2: 消融实验

**目的**: 逐个测试每个改动的独立贡献，每次只改一个变量

**基线**: C1-Step 1 的 fixed 模型

| 实验 | 改动 | 预期效果 | 复杂度 | 需重新生成.pt? |
|------|------|---------|--------|---------------|
| **2a** | Dropout 0.9 → 0.1 | ✅ Val=0.7216(+0.007), **Test=0.6936(-0.012)** Val未迁移Test | 低（改 config） | 否 |
| **2a'** | Dropout 0.9 → 0.3 | ✅ Val=0.7397(+0.025), **Test=0.6959(-0.010)** Val未迁移Test | 低（改 config） | 否 |
| **2b** | lr 3e-4 → 4e-4 + warmup 400 | +0.01~0.02 | 低（改 config） | 否 |
| **2c** | weight_decay 0 → 1e-5 | +0.005~0.01 | 低（改代码） | 否 |
| **2d** | max_epochs 50 → 150 | 更充分训练 | 低（改参数） | 否 |
| **2e** | Fe 原子词汇表扩展 | 不确定（解决 Heme OOD） | 高（改 transforms.py + 重生成 .pt） | **是** |
| **2f** | LR scheduler 改为监控 auc/val | 修复 aupr/auc 不匹配 | 低（改代码） | 否 |

**原则**:
- 每个实验独立跑，与 fixed 基线对比
- 先跑低复杂度的（2a-2d, 2f），后跑高复杂度的（2e）
- 每个实验记录: config diff + val AUC 曲线 + test AUC + 训练时间

**时间**: 每个实验 ~8-10 小时训练，总计 3-5 天

---

### C1-Step 3: 组合最优改动

**目的**: 把 Step 2 中有效的改动叠加在一起

**方法**:
- 选择所有 ΔAUC > 0 的改动
- 组合训练一个完整模型
- 在测试集上评估最终 AUC

**时间**: 1-2 天

---

### C1-Step 4: 多 seed 验证

**目的**: 验证结果的稳定性，排除随机性

**方法**:
- 使用 3 个不同 seed（3407, 42, 12345）重复 Step 3 的最优组合
- 报告 mean ± std

**时间**: 2-3 天（3 轮训练）

---

## 四、C2 P450全面数据集构建（2026-03-22 重新规划）

**目标**: 为P450建立与`ESIBank/small_family/Phosphatase/`同等规格的专属数据集，支持论文4种benchmark场景的4折交叉验证。

**详细计划**: 见 `C2_P450数据集构建/C2_计划与进度.md`

**工作目录**: `C2_P450数据集构建/`

### 数据源（68个数据库调研，按优先级排序）

| 来源 | 预估规模 | 状态 |
|------|---------|------|
| ⭐ P450Rdb v2.0 | ~10,957反应, ~850酶, 400+物种 | 🔜 最高优先 |
| ⭐ Plant P450 DB (ERDA) | 874序列(有底物注释) | 🔜 最高优先 |
| ⭐ PCPD | 181植物P450酶 | 🔜 最高优先 |
| RCSB PDB (已有) | 682对, 153酶 | ✅ |
| ESIBank P450子集 (已有) | 12,329对, 367酶 | ✅ 已验证 |
| CYPED + CPK | 8,614序列, 3,257化合物 | 🔜 |
| EnzymeMap | 63K反应中P450子集 | 🔜 |
| 其他A类(BM3/BRENDA/SABIO-RK等) | 见C2计划 | 🔜 |

### 实施阶段

| Phase | 内容 | 资源需求 | 状态 |
|-------|------|---------|------|
| Phase 0 | 审计ESIBank已有数据 | 本地 | ✅ 12,329条已验证 |
| Phase 1-2 | 下载A类核心数据+B类化合物库 | 本地, ~1-2GB | 🔜 |
| Phase 3 | 统一格式标准化 | 本地 | ⏳ |
| Phase 4 | 合并去重+泄露检测 | 本地 | ⏳ |
| Phase 5 | 负样本生成+4种Split | 本地 | ⏳ |
| Phase 6 | 结构获取+对接 | Cloud-2 GPU | ⏳ |
| Phase 7 | 特征生成+.pt缓存 | Cloud-2 GPU | ⏳ |
| Phase 8 | 4种场景Benchmark+模型优化 | Cloud-2 GPU | ⏳ |

### 关键决策

- **负样本分三类**: positives.csv + biological.csv(抑制剂等) + generated.csv(随机配对)
- **反应信息**: 有就记录(reactions.csv)，没有留空
- **P450验证**: 必须用UniProt protein_families + InterPro域，不可用EC号
- **最终输出**: 与`small_family/`结构一致，直接接入训练管线

---

## 五、时间线

| 日期 | 任务 | 预期产出 |
|------|------|---------|
| 3/19-20 | C1-Step 1: fixed 基线训练 | fixed Test AUC, ΔAUC |
| 3/20-24 | C1-Step 2: 消融实验 (2a-2d, 2f) | 各实验 ΔAUC |
| 3/24-26 | C1-Step 2e: Fe 词汇表扩展 | Fe 扩展 ΔAUC |
| 3/26-27 | C1-Step 3: 组合最优改动 | 最优组合 Test AUC |
| 3/28-30 | C1-Step 4: 多 seed 验证 | mean ± std |
| 3/22-3/28 | C2 Phase 1-5: 数据收集+标准化+去重+Split | P450_Family数据集 |
| 3/28-4/3 | C2 Phase 6-7: 对接+特征生成+.pt缓存 | 模型输入就绪 |
| 4/3-5 | C2 Phase 8: 4种场景Benchmark+优化 | 最终 AUC |
| 4/6-10 | 论文写作 | 完整论文草稿 |

---

## 六、关键依赖

1. **Cloud-2 服务器**（2×RTX 4090）：C1 全部实验在此运行
2. **Path B 的 .pt 缓存**（57GB）：已部署在服务器上
3. **Path B 的诊断结论**：指导改动方向（5 项核心发现）
4. **边修复代码**：Step 9 已验证（monkey-patch 在 main_training_pt.py）

---

## 七、风险与应对

| 风险 | 概率 | 应对 |
|------|------|------|
| 消融实验无显著改进 | 中 | 聚焦数据扩建（C2）作为备选创新点 |
| Fe 词汇扩展需重生成全部 .pt | 确定 | 评估增量生成可行性，或只重生成结构特征部分 |
| 服务器临时不可用 | 低 | 本地单卡跑小规模验证 |
| 时间不足（截止 4 月中旬） | 中 | 优先 C1（架构优化），C2 可选 |

---

## 八、进度追踪

| Step | 内容 | 状态 | 结果 |
|------|------|------|------|
| C1-Step 1 | fixed 基线训练 | ✅ 已完成 | **Val AUC=0.7145 (ep16 best)**, **Test AUC=0.7060** (AUPR=0.2362, 8841样本). early stopped ep31, 32ep ~5.2h ~2.7 it/s. 边修复未提升AUC（legacy_bug Val=0.722/Test=0.7244 vs fixed Val=0.7145/Test=0.7060, Δ=-0.018）. |
| C1-Step 2a | Dropout 0.9→0.1 | ✅ 已完成 | Val AUC=0.7216 (+0.007, ep17), **Test AUC=0.6936** (-0.012). Val改善未迁移到Test. |
| C1-Step 2a' | Dropout 0.9→0.3 | ✅ 已完成 | Val AUC=0.7397 (+0.025, ep49), **Test AUC=0.6959** (-0.010). Val改善未迁移到Test. 结论: dropout改动不纳入Step 3组合. |
| C1-Step 2b | lr 4e-4 + warmup 400 | ⏳ | |
| C1-Step 2c | weight_decay 1e-5 | ⏳ | |
| C1-Step 2d | max_epochs 150 | ⏳ | |
| C1-Step 2e | Fe 原子词汇表扩展 | ⏳ | |
| C1-Step 2f | LR scheduler → auc/val | ⏳ | |
| C1-Step 3 | 组合最优 | ⏳ | |
| C1-Step 4 | 多 seed 验证 | ⏳ | |
| C2 Phase 0 | 审计ESIBank已有数据 | ✅ 已完成 | 12,329条确认可用, 367/389酶 |
| C2 Phase 1-2 | 数据下载(A类+B类) | 🔜 进行中 | 68个数据库调研完成, ⭐P450Rdb/PlantP450DB/PCPD最高优先 |
| C2 Phase 3 | 格式标准化 | ⏳ | |
| C2 Phase 4 | 合并去重 | ⏳ | |
| C2 Phase 5 | 负样本+4种Split | ⏳ | |
| C2 Phase 6 | 结构+对接 | ⏳ | |
| C2 Phase 7 | 特征+.pt缓存 | ⏳ | |
| C2 Phase 8 | EXP001 random_split 基线 | ✅ 已完成 | **Val AUC=0.7544, Test AUC=0.7730** (4×4090 DDP, 89ep, best=ep73) vs ESIBank P450 internal 0.638 → +0.135 |
| **C3-Step 1** | **底物分类** | ✅ 已完成 | NPClassifier → 15类, 2,125全部分类 |
| **C3-Step 2** | **多轮验证+修正** | ✅ 已完成 | v5: 5优先级管线(P1 Gold+P2 NPC Superclass+P4 AA SMARTS), confirmed 1,773(83.4%)/review 262(12.3%)/other 90(4.2%), 60多标签(2.8%), 150抽检~89%, Codex 8轮+3轮抽检 |
| **C3-Step 3** | **352 review/other 全量Agent验证** | ✅ 已完成 | v6 FINAL: 20批Sonnet Agent文献搜索+Codex审核, 97升级confirmed+255确认OTHER, review 262→0, confirmed 1,773→1,870(88.0%), other 90→255(12.0%), 多标签60→63 |
| C3-Step 4 | 按类别聚合预测 | ⏳ | 酶序列 → ESM → MLP → 7-sigmoid → BCEWithLogitsLoss |
| C3-Step 5 | 评估与分析 | ⏳ | |
