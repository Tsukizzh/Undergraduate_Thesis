# Step 10: .pt 训练管线构建 — Session Log

> **开始日期**: 2026-03-16
> **目标**: 构建基于 .pt 分片缓存的完整训练管线，替代 LMDB 数据加载，消除 mmap 抖动和 CPU k-NN 瓶颈

---

## 背景

Step 9 训练（LMDB 模式）受到严重的性能瓶颈：
- 60GB LMDB 在 32GB RAM 上 mmap 抖动，稳态速度 2.09 it/s
- CPU k-NN 占数据加载时间 50-70%
- GPU 95% 时间在空等数据

Step 10 的解决方案：将 LMDB 数据预处理为 .pt 分片文件，预计算 k-NN 边索引，训练时直接加载。

---

## 子任务进度

| 编号 | 任务 | 状态 | 说明 |
|------|------|------|------|
| 10a | pt_dataset.py（Dataset 类） | ✅ 完成 | Codex 3 轮审核，所有修复已应用 |
| 10b | 训练脚本（PyTorch Lightning） | ⏳ 待做 | 接入 PtCacheDataset，支持 edge_mode 切换 |
| 10c | 正确性验证 | ⏳ 待做 | 跑几个 batch，确认输出和旧 LMDB 管线一致 |
| 10d | 数据传输到老师服务器 | ⏳ 待做 | 44GB .pt 文件 + 脚本 |

---

## 文件索引

| 文件 | 位置 | 说明 |
|------|------|------|
| pt_dataset.py | `scripts/10_Step10_pt训练管线/` | PyG Dataset 类 |
| build_pt_cache.py | `scripts/09_Step9_AllSplit训练/` | .pt 预处理脚本（共用基础设施） |
| validate_pt_cache.py | `scripts/09_Step9_AllSplit训练/` | 分片结构完整性验证（共用） |
| ezspec_pt_v1/ | `data/10_Step10_pt训练/` | .pt 缓存数据 |

### ezspec_pt_v1/ 目录结构
```
ezspec_pt_v1/
├── manifest.pt              # 全局配置（k=48, 原子词汇表等）
├── enzymes/                 # ESM fp16 分片（26 个，每个 ~900MB，总 26GB）
│   ├── index.pt             # enzyme_id → shard_id + row_id
│   ├── esm_0000.pt ~ esm_0025.pt
├── substrates/              # GROVER + Morgan 分片（10 个，总 0.5GB）
│   ├── index.pt
│   ├── grover_0000.pt ~ grover_0009.pt
├── train/                   # 87 个 graph 分片 + index.pt（总 26GB）
├── val/                     # 3 个 graph 分片 + index.pt
└── test/                    # 5 个 graph 分片 + index.pt
```

---

## 详细进度日志

### 2026-03-16 Session 1: 10a pt_dataset.py

#### 第一版编写
基于 Codex 技术规格（2 轮预审），编写了完整的 PtCacheDataset 类：
- 返回 StructureSequenceData（继承 PyG Data，有正确的 `__inc__`）
- 分片 LRU 缓存（ShardCache）
- 向量化 one-hot 重建
- edge_mode（fixed/legacy_bug）和 dist_noise 支持
- 兼容 `follow_batch=['ligand_index']`

#### Codex 审核第 1 轮：发现 7 个问题
1. **严重 - 分片字段名不匹配**: pt_dataset.py 读 `embeddings`，但 build_pt_cache.py 写的是 `embedding`（无 s）
2. **严重 - 酶/底物 ID 匹配问题**: Codex 担心 local vs global ID 不匹配，实测确认 graph 分片存的已经是 global ID，无需转换
3. **严重 - ligand_x 特征顺序错误**: atomic_number/100 放在了最后（位置 33），但原始代码放在位置 11（紧跟 element one-hot 之后）
4. **中等 - valid_idx 不兼容**: 返回 range(n) 而非原始 dataframe 行号
5. **中等 - mask_use_* 未提供**: 当前 config 不需要，但不完整
6. **中等 - 性能问题**: GaussianSmearing 每个样本重新实例化；vocab tensor 每次调用重新创建；无酶/底物解码缓存
7. **低 - src/ 路径引导脆弱**: 硬编码 4 层 `..`

#### 修复 + Codex 审核第 2 轮
全部 7 个问题修复：
- 字段名对齐 build_pt_cache.py 的实际输出
- ligand_x 顺序改为 [elem_oh, atomic/100, aromatic, degree_oh, numhs_oh, hybrid_oh]
- 从 ATOM_FEATS 导入动态获取 hybridization 类别数
- 添加 EntityTensorCache（酶/底物解码+padding 结果 LRU 缓存）
- 预构建 PROTEIN_ELEMENTS_T / LIGAND_ELEMENTS_T 全局 tensor
- 复用 self._smear（GaussianSmearing 单例）
- valid_idx 返回 sample_ids.tolist()

#### ESM 分片 >2GB Windows 崩溃问题
测试 Dataset 时发现 torch.load 无法加载 >2GB 的 ZIP 文件（PyTorch 2.3.0 Windows bug）。
原始 7 个 ESM 分片每个 2.1GB。解决方案：
- 删除旧分片，改用 enzyme_shard_size=1000（原来按 dataset 分 7 个）
- 重新生成 26 个分片，每个 ~900MB
- 验证通过

#### 最终验证
```
Dataset size: 1482
embedding:           [1450, 1280] float32 ✅
enzyme_padding_mask: [1, 1450] bool ✅
grover:              [280, 2400] float32 ✅
protein_x:           [295, 28] ✅
ligand_x:            [295, 34] ✅
complex_edge_attr:   [14160, 39] ✅
embedding nonzero:   True ✅（ESM 数据成功加载）
grover nonzero:      True ✅（GROVER 数据成功加载）
```

---

## 血泪教训（本 Step 及之前积累）

### 教训 1：第一版代码就必须包含性能优化
**不是"先写功能再优化"，而是"优化就是功能的一部分"。**

| 场景 | 优化前 | 优化后 | 加速倍数 |
|------|--------|--------|---------|
| build_pt_cache.py graph 分片 | 单线程 5 小时 | 10 worker 2 分钟 | **150x** |
| bond 边过滤 | Python for 循环 | torch.isin 向量化 | **~100x** |
| 酶分片生成 | 单线程 8 分钟 | 未优化（教训） | 本应更快 |

### 教训 2：文件大小必须检查
- ESM 分片 2.1GB → Windows PyTorch torch.load 崩溃（ZIP 格式 >2GB bug）
- 解决：分片大小限制 <1.5GB
- **规则：任何 torch.save 的文件，assert 最终大小 < 1.5GB**

### 教训 3：字段名必须从实际文件验证
- pt_dataset.py 假设字段名是 `embeddings`（复数），但实际是 `embedding`（单数）
- **规则：写 Dataset 前，先 torch.load 一个分片，打印 .keys()**

### 教训 4：特征顺序必须从源码逐行追踪
- ligand_x 的 atomic_number/100 应在位置 11，但错放到位置 33
- 形状（34 维）完全一样，不会报错，但训练结果会静默跑偏
- **规则：和 Codex 一起逐行对照 transforms.py 的特征拼接顺序**

### 教训 5：断点续跑必须在第一版就加
- build_pt_cache.py 最初没有断点续跑，中断后只能从头
- 加了 shard 存在检查后，重启自动跳过已完成部分
- **规则：任何耗时 >10 分钟的脚本，必须有中间状态持久化**

### 教训 6：移动文件后必须全局搜索引用
- 移动 ezspec_pt_v1/ 从 Step 9 到 Step 10 后，CLAUDE.md、全局进度日志等多处引用了旧路径
- **规则：移动文件后立刻 grep 旧路径，更新所有引用**

---

## 研究全局工作流（四个阶段）

### 阶段 1：复现论文模型（Step 10-11）
- **目的**：用 ESIBank 数据在 all_split（最难场景）上训练，得到可靠的 baseline
- **数据**：ESIBank all_split（176K train）
- **管线**：.pt 训练管线（ezspec_pt_v1/）
- **产出**：两个 checkpoint（legacy_bug + fixed），消融表前两行

### 阶段 2：架构优化（Step 12-13）
- **目的**：逐个测试架构改进，找到最优组合
- **数据**：同上（ESIBank all_split），不改数据只改模型
- **管线**：同上（.pt 管线），只改模型代码和超参
- **产出**：消融实验表（Dropout / Fe+Heme / 酶侧融合 / 组合）

### 阶段 3：Benchmark 对比（Step 14）
- **目的**：验证优化后的模型在 P450 上是否真正提升
- **数据**：Path A 的 682 条 P450 独立测试集（已有特征，不需要 .pt）
- **管线**：**推理管线**（main_testing.py），不是训练管线
- **操作**：
  - 原模型 checkpoint → 推理 → AUC = X（Path A 已有，AUC=0.6636）
  - 优化模型 checkpoint → 推理 → AUC = Y
  - 对比 X vs Y，量化改进

### 阶段 4：P450 数据扩充 + 重训（Step 15，可选）
- **目的**：加入更多 P450 训练数据，进一步提升
- **数据**：ESIBank + P450 扩充（Plant P450 DB、P450Rdb 等）
- **管线**：.pt 训练管线，但需要扩展：
  - 新酶 → 追加 enzymes/ 分片
  - 新底物 → 追加 substrates/ 分片
  - 新配对 → 生成新 graph 分片（需先做 Vina 对接）
  - pt_dataset.py 加 `shared_dir` 参数（到时候再改）
- **产出**：最终模型 + 再次 benchmark

### 最终评估（Step 16）
- 多 seed 重复（≥3 个种子）
- 在 P450 测试集上报告均值 ± 标准差

## .pt 缓存的适用范围

| 阶段 | 需要 .pt 训练管线？ | 需要新 .pt 文件？ |
|------|:--:|:--:|
| 阶段 1（复现） | ✅ | 当前 ezspec_pt_v1/ 够用 |
| 阶段 2（架构优化） | ✅ | 不需要，改模型不改数据 |
| 阶段 3（benchmark） | ❌ 用推理管线 | 不需要 |
| 阶段 4（数据扩充） | ✅ | 需要扩展（追加分片） |

**核心原则**：酶和底物分片是一次性投入，graph 分片按需重生成。当前 ezspec_pt_v1/ 服务阶段 1-2，阶段 3 走推理管线，阶段 4 再扩展。

---

## 后续具体步骤

### Step 10c-d：验证 + 数据传输（1-2 天）
- 本地跑几个 batch 验证 .pt 管线正确性
- 传 44GB .pt 文件到老师服务器

### Step 11：基线训练（3-4 天）
- legacy_bug 版本（从头训练）→ 基线 AUC
- fixed 版本（从头训练）→ 量化边修复贡献
- 在老师服务器上跑（RTX 4090 + .pt）

### Step 12：独立消融实验（5-8 天）
- Dropout 0.9→0.1
- Fe 词汇 + Heme 数据（需重新生成含 Heme 的结构 + graph 分片）
- 酶侧结构融合（需实现 atom→residue→ESM 映射）

### Step 13：组合有效修复（1-2 天）

### Step 14：Benchmark（1-2 天）
- 用阶段 1-2 的 checkpoint 在 P450 测试集上推理
- 对比原模型 vs 优化模型

### Step 15：P450 数据扩充 + 重训（5-7 天，可和 Step 12 并行收集数据）

### Step 16：多 seed 最终评估（4-7 天）

**总计约 3-4 周，论文截止 4 月中旬。**
