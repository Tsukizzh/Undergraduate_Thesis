# Step 9-10 会话完整总结 (2026-03-16)

## 核心成果（3个突破）

### 🎯 成果1: Step 9 AllSplit 训练完成 → ep14 AUC=0.7667

**训练执行（2026-03-10~16）**:
- ep0-27 全部完成，过拟合已确认
- **FINAL BEST = ep14 AUC-ROC = 0.7667** ★
  - 比 ep8 (0.7422) 提高 +0.0245 (+3.3%)
  - 比论文"未知酶+底物"(0.7198) 超出 +4.6%
- AUC 趋势: 0.654(ep2) → 0.742(ep8) → **0.7667(ep14)** → 0.766(ep19) → ~0.766(ep27)
- 过拟合指标：train_loss 0.352→0.088↓，val_loss 0.328→0.374↑
- EarlyStopping wait=13/15@ep27，若继续训练在ep29触发

**数据泄露审计（发现）**:
- 酶序列重叠：27个(0.15%) — **可忽略**
- 底物Morgan指纹重叠：585个(6.56%) — **不可忽略** ⚠️
  - 跨数据集污染（例：BRENDA test样本出现在Esterase train）
  - 含义：AUC 0.7667 严格说是"未知酶"不是"未知酶+底物"场景
  - **Phase 3a行动**：分离泄露/非泄露样本AUC以定量影响

**架构问题发现**:
- Edge ordering bug 已在 ep0 前通过 monkey-patch 修复
- ep14 基线包含 edge fix（非论文原始版本）
- **Phase 2行动**：创建"legacy_bug baseline"对照实验

### 🚀 成果2: .pt 缓存 v3 完成 → 7.56 it/s

**问题**:
- 60GB LMDB >> 32GB RAM 导致页面缓存污染
- 速度衰减：3.8→2.09 it/s（100+ min/epoch）
- 导师服务器 16GB RAM 会更严重

**三轮迭代**:

| 版本 | 方案 | 性能 | 状态 |
|------|------|------|------|
| v1 | 87×300MB shard | 不可用(SSD 100%) | ❌ |
| v2 | enzyme/substrate flatbin | 0.5 it/s | ⚠️ |
| **v3** | **176K×160KB per-sample** | **7.56 it/s** ★ | ✅ |

**存储布局（44GB 总）**:
```
ezspec_pt_v1/
├── train/samples/  26GB (176,843 个.pt)
├── val/samples/    0.1GB (1,482 个)
├── test/samples/   0.7GB (8,841 个)
├── enzymes/        17GB (26个 fp16 分片 + 索引)
└── substrates/     0.5GB (10个 GROVER 分片 + 索引)
```

**性能对比**:
- LMDB (衰减版)：2.09 it/s，100+ min/epoch
- .pt per-sample：7.56 it/s，27 min/epoch
- **改进：+262% 速度，-73% 时间** ★

**构建完成**:
- 执行时间：2 分钟 (10 workers)
- 验证：✅ 所有 186,166 样本通过
- 脚本：`scripts/09_Step9_AllSplit训练/build_pt_cache.py`

### ⚙️ 成果3: Dataset 类实现 (Step 10a) → Codex 3轮审核通过

**文件**: `scripts/10_Step10_pt训练管线/pt_dataset.py` (350 lines)

**功能**:
- 加载 per-sample .pt 图数据 (160KB)
- 动态加载 enzyme/substrate shard 特征
- Multi-worker safe，支持 fallback 到 flatbin.bin
- 单元测试全部通过

**Codex 三轮审核发现**:
1. manifest 结构不清 → 明确 split_id 映射
2. enzyme_shard 无缓存 → 改为内存缓存 16 个 shard
3. 多 worker 无法 pickle → 改用 collate_fn + shared memory

---

## 执行计划修正（5 个关键问题发现）

与 Codex 深度复审发现的问题及修正：

| 问题 | 原计划 | 改正方案 | 影响 |
|------|--------|---------|------|
| ❌ Edge fix 冗余 | 可跳过 | Phase 2：legacy_bug 对照实验 | Ablation 需重新训 |
| ❌ 消融可 resume | 从 ep14 继续 | 必须从头训 (random init) | 多轮实验 |
| ❌ Fe 词汇表直接有效 | 添加 Fe 即可 | 需重生成 structure (--include-heme) | B3 数据集扩展 |
| ❌ .pt schema 完整 | 按样本存 ESM | ESM 分片存，只样本存 graph | 已实现 v3 |
| ❌ 残基融合直接可行 | 直接融合 | 需完整映射链路 (atom→residue→ESM) | 架构修正 |

**改正后的 7 阶段计划仍可行，但需多轮并行**

---

## 性能优化方案评估（9 种方案排名）

根据代码质量、成本、时间进行评分：

| 排名 | 方案 | 性能 | 代价 | 时间 | 备注 |
|------|------|------|------|------|------|
| **1** | **F: readahead=False** | **+5-15%** | **1行改动** | 即时 | 已在源码 ✅ |
| 2 | A: 分片 LMDB | +5-10% | 小 | 1天 | 部分改善 |
| 3 | B: .pt 无 k-NN | +10-20% | 70-80GB | 2天 | 本方案已实施 |
| 4 | C: .pt 含 k-NN | +100-150% | 250-300GB | 5天 | 存储过大 |
| 5 | D: GPU k-NN | +50-80% | 中 | 3天 | 仍需改 k-NN |
| 备选 | G: 云服务器(100GB) | 6 it/s | ¥0.99/h | 1周 | 成本 ¥10-20 |
| 备选 | H: 升级服务器 RAM | +40-60% | ¥2000 | 2-3周 | 长期投资 |

**已选**：B (.pt per-sample v3)

---

## 性能基准测试（2026-03-16）

### 数据加载速度

```
配置                      速度         每 Epoch    内存占用    衰减
───────────────────────────────────────────────────────────────
LMDB baseline             1.71 it/s    2h          15GB        无
LMDB 优化                 4.3 it/s     49min       20GB        无
LMDB cache (final)        3.8→2.09     100+min     15-30GB     严重 ⚠️
.pt Per-sample (v3)       7.56 it/s    27 min      ~10GB       无 ★
```

### 内存占用详解

```
.pt Per-sample:           LMDB Flatbin:
PyTorch 参数：1.8MB      PyTorch 参数：1.8MB
Dataloader 缓冲：10GB    LMDB 污染：20-30GB
页缓存：0GB              Dataloader：5GB
总计：~10GB safe         总计：~30GB (16GB时严重)
```

---

## 下一步：7 阶段计划（2-3 周）

### Phase 1：服务器迁移与速度验证 (~3 天，2026-03-17~19)

**任务**:
1. 本地打包：44GB → 13GB (tar.gz, 压缩率 70%)
2. 上传：→ dc@jqlab4090:~/zhuangzeheng/EZSpecificity/data/
3. 解压验证：186,166 个文件
4. 首次速度基准：100 batches 测试
   - 期望：≥4 it/s (16GB RAM 下)
   - 如果 <3 it/s，调整 num_workers=2 或 1

**产出**:
- .pt 缓存在服务器就绪
- 速度基准数据（确认 I/O bound 已解除）

### Phase 2：Baseline 实验 (~3 天，2026-03-20~22)

**目标**：量化 edge fix 的贡献值

**实验**:
- Baseline A：edge_fix (当前 ep14)
- Baseline B：legacy_bug (重现论文原始 bug 版)
- 对比：edge fix 的 AUC 改进量

### Phase 3a：数据泄露量化 (~2 天，2026-03-23~24)

**目标**：分离泄露样本对 AUC 的影响

**方法**:
1. 标记 585 个泄露底物
2. 分离计算：AUC(泄露) vs AUC(非泄露)
3. 估计通胀值：0.01~0.04 范围

### Phase 3b：干净数据集（可选）

**触发条件**：若泄露影响 >0.02

**工作**：重生成数据集（expensive，2-3 周）

### Phase 4：Fe + Heme 扩展 (~3 天，2026-03-25~27)

**任务**:
1. 添加 Fe(26) 到蛋白原子词汇表
2. 重生成 structure_features (--include-heme)
3. 结合 B3 扩展数据集（含 Heme）
4. 训练对比

### Phase 5：融合有效性修复 (~3 天，2026-03-28~30)

**核心**：修复残基级融合映射

**实现**：
- atom → residue_instance → ESM_position 三级映射
- GNN output 聚合到残基
- 注意力激活模式分析

### Phase 6：P450 数据增强 (~3-5 天，2026-03-31~04-04)

**选项**:
1. 添加 ESIBank P450 (389 个) 到训练集
2. 或添加 PDBbind P450 数据
3. 重新训练验证改进

### Phase 7：最终多轮评估 (~2 天，2026-04-05~06)

**要求**:
- 3+ random seeds
- 固定测试集
- 论文级结果汇总

---

## 文件清单

### 新建脚本

| 文件 | 说明 | 状态 |
|------|------|------|
| `scripts/09_Step9_AllSplit训练/build_pt_cache.py` | .pt 缓存构建 | ✅ |
| `scripts/10_Step10_pt训练管线/pt_dataset.py` | Dataset 类 | ✅ |
| `scripts/10_Step10_pt训练管线/main_training_pt.py` | 训练脚本框架 | 🔄 |

### 新增数据

| 路径 | 说明 | 大小 | 状态 |
|------|------|------|------|
| `data/10_Step10_pt训练/ezspec_pt_v1/` | .pt 缓存 | 44GB | ✅ |

### 日志文件

| 文件 | 说明 | 状态 |
|------|------|------|
| `sessions/10_Step10_pt训练管线/session_log.md` | Session 记录 | ✅ |

---

## 关键数据表

### 性能基准数据（训练衰减分析）

| Epoch 阶段 | 速度 | 原因 |
|-----------|------|------|
| ep0-3 | 4.5-5.5 it/s | LMDB 缓存未满，页面缓存命中率高 |
| ep8-12 | 3.0 it/s | 数据访问模式多样化，缓存效率下降 |
| ep14-27 | 2.09 it/s | LMDB mmap 持续污染，99% RAM 使用，swap 活跃 |

### AUC 训练曲线完整数据

| Epoch | AUC-ROC | AUPR | Train Loss | 说明 |
|-------|---------|------|------------|------|
| 0 | 0.620 | 0.180 | 0.352 | 初始 |
| 2 | 0.654 | 0.291 | 0.277 | |
| 8 | 0.7422 | 0.2943 | — | 旧 best |
| **14** | **0.7667** | **0.2928** | **0.109** | **NEW BEST** ★ |
| 19 | 0.766 | 0.293 | 0.098 | AUC 小降 |
| 27 | ~0.766 | ~0.290 | ~0.088 | 过拟合：train↓ val↑ |

---

## 时间表（确认）

- **2026-03-17**：完成打包，开始迁移
- **2026-03-18~19**：服务器速度测试
- **2026-03-20~22**：Phase 2 baseline 实验
- **2026-03-23~24**：Phase 3a 数据泄露量化
- **2026-03-25~27**：Phase 4 Fe+Heme
- **2026-03-28~30**：Phase 5 融合修复
- **2026-03-31~04-04**：Phase 6 P450 增强
- **2026-04-05~06**：Phase 7 最终评估
- **2026-04-10 前**：论文结果交付

---

## 用户决策点

### 关键选择 1：数据泄露处理

**选项 A**（推荐）：保留数据，文档泄露量
- 成本：0，时间：0
- 如果 Phase 3a 显示影响 ≤0.02，直接文档+保留

**选项 B**：重生成干净数据集
- 成本：时间 2-3 周
- 仅在 Phase 3a 显示影响 >0.02 时执行

### 关键选择 2：Fe+Heme 优先级

**高优先级**（推荐）：Phase 4 执行
- P450 Heme OOD 问题（effect=-0.2322）
- 可能带来显著改进

**低优先级**：可选或后续
- 取决于 Phase 2 baseline 结果

---

## 已记录的内存文件

### Set A (Claude Code Auto Memory)
- ✅ MEMORY.md 更新
- ✅ research-progress.md 追加

### Set B (用户项目进度日志)
- ✅ 全局进度日志.md 大幅更新
- ✅ sessions/10_Step10_pt训练管线/session_log.md 新建

### Set C (Serena Project Memory)
- ✅ project_overview.md 更新

---

**文档版本**: v1.0
**完成时间**: 2026-03-16 23:00
**下一检查点**: 2026-03-18 (Phase 1 速度测试)
