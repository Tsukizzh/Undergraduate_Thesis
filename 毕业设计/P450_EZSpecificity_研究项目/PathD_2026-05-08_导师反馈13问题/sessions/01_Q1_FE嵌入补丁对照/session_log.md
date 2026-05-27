# Q1 — FE 嵌入补丁对照实验

## 老师原话

> 对于原模型，图神经网络中其缺乏 FE 嵌入，那么我们重新训练一个有 FE 嵌入的模型，再用两个模型在同一 ESIBANK 测试集中进行对比。
> （完整的测试集 or 测试集中的 P450？）
> 假设我们发现有了铁后的模型，在测试集上表现更好，那么我们把这一类别单独拉出来（也就是 P450），无论其是否包含铁都没什么影响。

## 状态

🟡 **设计就绪，待启动训练**（Q7 已完成，实验设计已锁定）

## Q7 结论传入（2026-05-16）

✅ Q7 审计已完成（详见 `sessions/07_Q7_原论文Heme处理排查/session_log.md`），三方判定：
- 原模型蛋白侧 GNN **看不见** HEM/HETATM（解析器严格 ATOM 过滤）
- 蛋白原子词汇表**不含** Fe
- 任何代码路径都**未把** HEM-Fe 送入蛋白/结构通道

**Q1 因此立得住**：在 ESIBank 全集上做"加 Fe 嵌入 vs 不加"的对照实验是合理实验，可以填补原论文的实现 gap。

## 老师设计完整逻辑

| 阶段 | 内容 |
|---|---|
| 训练 | ESIBank **全集**（约 32 万对）训练两个模型：(A) 原版 / (B) 加 Fe+HEM 嵌入版 |
| 测试-步骤1 | 在**完整 ESIBank 测试集**上对比 A vs B 的 Test AUC（跨家族）|
| 测试-步骤2 | 把 **P450 类别从测试集拉出来单独评估**，比较 A vs B 在 P450 子集上的差距 |
| 判读 | 如果 P450 子集 A ≈ B（无差）→ 证明对 P450 来说 Fe 是常量、无信息量 → 整体提升来自跨家族判别 |

## 与已有工作的关系

### EXP002a (Fe/HEM 补丁，P450 数据)
- 在 P450 干净 AllFix 数据上加 Fe 元素 + HEM 残基类型 + is_hetero 标志
- feature_dim 28→31
- **结果**：Test AUC 0.9270 vs 基线 0.9320 = **-0.005**（Fe 反而掉点）
- 结论：**P450 干净数据上 Fe/HEM 不带来增益**

### 本问题 vs EXP002a 的关键差异

| 维度 | EXP002a | 老师本问题 |
|---|---|---|
| 训练集 | P450 子集（1622 酶 / 44090 样本） | ESIBank **全集**（323,783 对 / 8,124 酶） |
| 测试集 | P450 测试集（11,000 样本） | "完整测试集 or 测试集中的 P450"（待澄清） |
| 评估目的 | P450 内部 Fe/HEM 是否有用 | Fe 嵌入对**通用模型**是否有用 |

## 代码改造点（与 EXP002a 已做的改动一致，复用）

1. `src/Datasets/Structure/protein_ligand.py:67` — 解析器允许 HETATM
2. `src/Datasets/Structure/protein_ligand.py:20` — `AA_NAME_SYM` 加 HEM 及其变体（HEM/HEC/HEB/HEA）
3. `src/Datasets/Structure/transforms.py:23` — 蛋白原子词汇表 `[1,6,7,8,16,34]` → `[1,6,7,8,16,26,34]`（加 Fe）
4. 重建 LMDB → pt_cache → 训练（B 版与原论文 A 版唯一差异为上述三处）

## 资源需求与待解决问题

1. **ESIBank 32 万对数据是否在北京服务器**？需要现场确认（4×RTX5090 + 754Gi RAM 算力足够）
2. **ESIBank 的 ESM / GROVER cache 是否公开 / 已有**？如果都要重新生成是巨大工程
3. **原论文 ckpt** 是否可作为 A 版基线直接评估（避免重训 A 版的算力）
4. P450 在 ESIBank 测试集中如何识别 → 已有 PathC 的"389 个真 P450"列表

## 工作量估计

| 阶段 | 估时（4×RTX5090）| 备注 |
|---|---|---|
| 数据准备（下载 ESIBank + cache 检查）| 0.5 - 2 天 | 取决于 cache 公开性 |
| 改造代码 + 单元测试 | 0.5 天 | EXP002a 已有改造模板 |
| ESM / GROVER cache 重新生成（如需）| 1 - 2 天 | 看是否复用 |
| Pt_cache 构建 | 0.5 - 1 天 | |
| 训练 A 版 + B 版（或仅 B 版，A 版用原 ckpt）| 2 - 3 天 / 版 | 256 epoch 量级 |
| 评估（完整 + P450 子集）| 0.5 天 | |
| **合计** | **~5 - 10 天** | 取决于多少步可省 |

## 待办

- [x] 完成 Q7（原论文 Heme 处理排查）
- [x] 锁定改造代码点（3 处）
- [ ] 现场确认北京服务器是否有 ESIBank 数据 + cache
- [ ] 评估原论文 ckpt 能否直接当 A 版基线
- [ ] 跟老师确认是否走"原 ckpt 当 A 版 + 我们只训 B 版"的省算路线
- [ ] 重训 B 版
- [ ] 在 ESIBank 完整测试集评估 → 拆 P450 子集对比
- [ ] 报告 + 写入毕业论文

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |

## 2026-05-27 EXP002 Fe/HEM overlay 数据准备完成

本次目标是让 Q1-EXP002 在开 GPU 后可以直接开始训练。已完成的数据准备和验证如下。

### 结果位置

- EXP002 实验目录：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/`
- EXP002 训练缓存：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1/`
- 官方 ESIBank PDB：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/official_esibank_brenda_pdb_best1_20260527/pdb/`
- RCSB mmCIF：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/raw_mmcif_rcsb_heme_fe_v1/`
- 准备状态报告：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1/reports/full_v2/exp002_readiness_validation_20260527.json`
- 使用说明：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/README_EXP002_READY_20260527.md`
- 启动脚本：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/scripts/run_exp002_fe_overlay_train.sh`

### 数据生成结果

EXP002 使用 EXP001 的官方 A-only PT cache 作为底座，通过硬链接复用未改变样本，只替换能加入 HEM/Fe 的样本文件。这样不会修改 EXP001 原缓存，也避免再次复制整套大缓存。

四个分片全部完成：

| 分片 | 计划样本 | ok | failed | missing_base_sample |
|---|---:|---:|---:|---:|
| shard_0 | 797 | 797 | 0 | 0 |
| shard_1 | 796 | 796 | 0 | 0 |
| shard_2 | 796 | 796 | 0 | 0 |
| shard_3 | 796 | 796 | 0 | 0 |
| 合计 | 3185 | 3185 | 0 | 0 |

按 split 统计：

| split | 总样本数 | 增强样本数 |
|---|---:|---:|
| train | 229855 | 2362 |
| val | 22952 | 257 |
| test | 53588 | 566 |

### 验证结果

已抽样读取 train、val、test 的增强样本。增强样本均新增 HEM/Fe 异质原子；典型样本新增 `43` 个异质原子，其中 Fe 元素数为 `1`。未增强样本仍与 EXP001 base cache 保持硬链接，不带 `protein_is_hetero` 字段。

训练读取脚本已改为 31 维 protein 输入：元素 7 维、氨基酸 22 维、backbone 1 维、hetero 1 维。无卡模式下做了 CPU 前向验证：

- 单个增强样本 `sample_id=130`：`protein_x=(344, 31)`，hetero 特征列合计 `43`，模型前向通过。
- 普通样本 `sample_id=0` 和增强样本 `sample_id=130` 混合成 batch：`protein_x=(655, 31)`，模型输出 `logits shape=[2,1]`。

服务器空间检查：`/root/autodl-tmp` 剩余约 `96G`。

### 明天启动命令

3 卡正式训练建议先从这个命令开始：

```bash
cd /root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay
DEVICES=3 BATCH_SIZE=40 NUM_WORKERS=8 MAX_EPOCHS=200 scripts/run_exp002_fe_overlay_train.sh
```

如果先做短试跑：

```bash
cd /root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay
DEVICES=1 BATCH_SIZE=40 NUM_WORKERS=8 MAX_EPOCHS=1 LIMIT_TRAIN_BATCHES=400 LIMIT_VAL_BATCHES=20 scripts/run_exp002_fe_overlay_train.sh
```

默认不会自动关机。需要自动关机时再加 `EXTRA_ARGS="--shutdown"`。

### 补充对齐审计

子智能体复核提醒后，又补做了 EXP002 与 EXP001 base cache 的全量对齐审计，结果写入：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1/reports/full_v2/exp002_alignment_validation_20260527.json`

审计结果：

| 项目 | 结果 |
|---|---|
| train 样本数和顺序 | 与 EXP001 一致，`229855` |
| val 样本数和顺序 | 与 EXP001 一致，`22952` |
| test 样本数和顺序 | 与 EXP001 一致，`53588` |
| 增强样本 label 是否变化 | `0` 个变化 |
| 增强样本 enzyme_id 是否变化 | `0` 个变化 |
| 增强样本 substrate_id 是否变化 | `0` 个变化 |

增强样本标签分布：

| split | 负样本 | 正样本 |
|---|---:|---:|
| train | 2196 | 166 |
| val | 243 | 14 |
| test | 533 | 33 |

这个验证说明 EXP002 cache 是在 EXP001 完整 train/val/test 上做增量增强，不是只剩下一个 Fe/HEM 子集。

## 2026-05-27 EXP001 结果目录重命名

服务器上 EXP001 的正式结果目录已改成更短、更靠前的名字：

- `/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/results/00_EXP001_final/`

同时新增中文说明文件：

- `/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/results/00_README_EXP001_FINAL.md`

原目录内容没有删除，只是从较长目录名移动到新的短目录名下。后续查看 EXP001 最终测试结果时，优先看 `00_EXP001_final/`。

## 2026-05-27 EXP001 停止继续训练并做测试集评估

### 当前决定

EXP001 不继续训练到 200 epoch。按当前已经得到的最佳 checkpoint，直接在官方 A-only 测试集上做最终测试评估。

这次操作只停止了正在运行的 EXP001 训练进程，没有删除训练日志、checkpoint 或结果文件。

### 停止训练前的状态

正式训练 run：

`Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040`

停止时间：`2026-05-27 04:32:15`

停止前最新完整训练记录到 `epoch 23`。最新几轮验证结果如下：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 19 | 0.209507 | 0.234218 | 0.890818 | 0.545088 |
| 20 | 0.203991 | 0.232862 | 0.891122 | 0.545077 |
| 21 | 0.202312 | 0.233275 | 0.893959 | 0.539946 |
| 22 | 0.197529 | 0.233592 | 0.893906 | 0.548956 |
| 23 | 0.196293 | 0.232295 | 0.894143 | 0.551813 |

用于测试的 checkpoint：

`results/Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040/checkpoints/pt-Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040-ep23-auc0.8965.ckpt`

说明：checkpoint 文件名里的 `auc0.8965` 来自训练脚本保存 checkpoint 时的命名；`metrics.csv` 里 epoch 23 对应的完整行记录为 `val_auc = 0.894143`。后续汇报时优先引用 `metrics.csv` 和测试 JSON 里的实际数值。

### 测试设置

测试命令使用同一个 PT cache 和同一个数据划分，只读取 `test` split：

| 项目 | 设置 |
|---|---|
| 脚本 | `scripts/main_training_pt_speed_v2_nograd_cache256.py` |
| 模式 | `--test-only` |
| 配置 | `configs/config_accum1_perf.yml` |
| cache | `data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1` |
| edge mode | `legacy_bug` |
| batch size | `40` |
| num workers | `8` |
| GPU | 单卡 `cuda:0` |

测试输出：

`results/Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040/test_eval_best_checkpoint_20260527_043410.json`

测试日志：

`logs/Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040_test_eval_20260527_043410.log`

### 测试结果

| 指标 | 数值 |
|---|---:|
| test AUC-ROC | 0.894871 |
| test AUPR | 0.569394 |
| 测试样本数 | 53588 |
| 正样本数 | 5941 |
| 负样本数 | 47647 |
| 推理耗时 | 264.0 秒 |
| 推理速度 | 约 203 samples/s |

分家族结果：

| family | AUC-ROC | AUPR | 样本数 | 正样本数 |
|---|---:|---:|---:|---:|
| brenda | 0.889530 | 0.484944 | 43235 | 4095 |
| Duf | 0.834040 | 0.514861 | 473 | 54 |
| Esterase | 0.945296 | 0.851341 | 2451 | 532 |
| Gt_acceptor | 0.895373 | 0.736088 | 834 | 177 |
| Nitrilase | 0.936893 | 0.645473 | 117 | 14 |
| Phosphatase | 0.876881 | 0.645052 | 6291 | 978 |
| Thiolase | 0.886676 | 0.888675 | 187 | 91 |

### 当前结论

EXP001 已经得到一个可汇报的官方 A-only 测试集结果。它可以作为 Q1 后续 Fe/HEM 嵌入对照实验的原模型基线。

下一步再处理 EXP002：准备带 Fe/HEM 信息的数据，并训练“加入 Fe/HEM 嵌入”的对照模型。

## 2026-05-27 EXP002 官方 PDB 子集下载与上传完成

为了准备 EXP002 的 Fe/HEM 对照实验，已从 Google Drive 元数据定位官方 ESIBank `brenda/structure/complex` 中与 RCSB best1 候选相关的 PDB 文件。

本次只使用官方 brenda complex 目录下的 PDB。对于只出现在 `small_family` 的同名 PDB，没有拿来替代 brenda 缺失项。

### 本地位置

`D:\ESIBank\EXP002_official_brenda_pdb_best1_20260527`

压缩包：

`D:\ESIBank\archives\EXP002_official_brenda_pdb_best1_20260527_package_v2.tar.gz`

### 服务器位置

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/official_esibank_brenda_pdb_best1_20260527`

服务器压缩包：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/archives/EXP002_official_brenda_pdb_best1_20260527_package_v2.tar.gz`

### 审计结果

| 项目 | 数值 |
|---|---:|
| RCSB best1 涉及的 structure_index | 7184 |
| 在 Google Drive 官方 brenda/structure/complex 找到的 PDB | 6596 |
| 未在官方 brenda/structure/complex 找到的 PDB | 588 |
| 已下载并上传的 PDB | 6596 |
| PDB 总字节数 | 4188644489 |
| 本地 PDB 文件头异常 | 0 |
| 服务器剩余空间 | 约 97G |

### 关键文件

- `manifests/exp002_google_drive_pdb_manifest_for_rcsb_best1_20260527.csv`
- `official_esibank_brenda_pdb_best1_20260527/local_audit_20260527.json`
- `official_esibank_brenda_pdb_best1_20260527/server_audit_20260527.json`
- `archives/EXP002_official_brenda_pdb_best1_20260527_package_v2.tar.gz`

### 当前状态

EXP002 的结构下载层已经准备好：服务器上同时有 RCSB HEM/Fe mmCIF 和官方 brenda PDB 子集。下一步应先做小样本闭环，确认 Fe/HEM 特征如何并入现有 PT cache，再批量生成 EXP002 训练 cache。

## 2026-05-27 三卡正式训练加速记录

### 当前结论

Q1-EXP001 正式训练已经从旧脚本切到新脚本继续跑。新脚本保留模型、数据、batch size、三卡 DDP、checkpoint 和结果目录，只去掉每个训练 step 的梯度范数计算，并改用已有的 enzyme cache 256 数据集版本。

这次切换是从 `epoch 14` 的 checkpoint 之后进行的，没有丢掉已完成 epoch 的结果。

### 为什么要切换

旧脚本实际仍在每个优化器 step 里计算 `grad_norm`：

- 遍历模型参数；
- 对 GPU 上的梯度做范数计算；
- 用 `.item()` 把结果拉回 CPU。

这个过程会造成 CPU 和 GPU 同步等待。三卡 DDP 下，一个 rank 等待会拖慢整体 step。旧日志里 `grad_norm=inf`，说明这个记录项已经出现无效值，继续保留没有实验价值。

### 新增脚本

服务器新增脚本：

- `scripts/main_training_pt_speed_v2_nograd_cache256.py`

它调用已有数据集脚本：

- `scripts/pt_dataset_speed_cache_v2.py`

原脚本没有覆盖：

- `scripts/main_training_pt_speed_v1.py`

### 正式训练命令核心参数

| 项目 | 当前值 |
|---|---|
| run name | `Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040` |
| 训练脚本 | `scripts/main_training_pt_speed_v2_nograd_cache256.py` |
| cache | `full_legacy_cuda_20260526_161909_per_sample_graph_v1` |
| edge mode | `legacy_bug` |
| batch size | `40` |
| GPU | `3` |
| workers | `8` |
| prefetch factor | `1` |
| train_in_order | `false` |
| resume | `results/<run>/checkpoints/last.ckpt` |
| skip auto test | 开启 |

### 速度对比

| 阶段 | 脚本 | 训练段速度 | 完整 epoch 耗时 |
|---|---|---:|---:|
| epoch 14 | `main_training_pt_speed_v1.py` | 训练段约 `3.6-3.7 it/s` | `9:38` |
| epoch 15 | `main_training_pt_speed_v2_nograd_cache256.py` | `4.06 it/s` | `8:40` |

这里的“完整 epoch”包含训练后验证。新脚本完整一轮节省约 58 秒。

### 指标状态

| epoch | val AUC | val AUPR | 备注 |
|---:|---:|---:|---|
| 13 | `0.871588` | `0.502659` | 旧脚本 |
| 14 | `0.879724` | `0.521623` | 旧脚本，切换前最后完整 epoch |
| 15 | checkpoint 显示 best 约 `0.890` | `metrics.csv` 中该行 AUC/AUPR 为空 | 新脚本，需后续从 checkpoint 或 Lightning 日志复核 |

epoch 15 的 `metrics.csv` 行里 AUC/AUPR 为空，但 checkpoint 日志显示 `auc/val` improved 到约 `0.890`。这属于记录侧差异，后续整理结果时要回到 checkpoint 和 Lightning 日志核对。

### GPU 利用率观察

新脚本运行时做了 60 秒采样：

| GPU | 平均利用率 | 平均显存 | 平均功耗 |
|---|---:|---:|---:|
| GPU0 | `34.5%` | `27156 MiB` | `230.1 W` |
| GPU1 | `98.5%` | `24753 MiB` | `225.3 W` |
| GPU2 | `98.0%` | `26993 MiB` | `213.6 W` |

GPU1 和 GPU2 基本满载。GPU0 的利用率数值偏低，但功耗和显存都很高，说明它不是完全空闲。后续若继续压榨三卡，需要重点处理 DDP rank 负载不均和 CPU 侧构图开销。

### 下一步可加速方向

优先级最高的是两个方向：

1. 做按图规模分桶的 batch sampler，让每张卡拿到的图复杂度更接近，减少某个 rank 先算完后等待。
2. 把 CPU 侧 edge 特征中不随训练变化的部分预计算进缓存，只保留训练噪声相关部分在线计算。

这两个方向会动训练代码和缓存结构，适合作为后续新的加速实验，不应悄悄混入当前 Q1-EXP001 baseline。

### 复核意见

子智能体只读复核后给出的建议是：当前不应继续中断正式训练去试新方案。`main_training_pt_speed_v2_nograd_cache256.py` 已经带来约 10% 的完整 epoch 提速，继续打断会消耗确定的服务器时间，还会引入新的正确性风险。

后续更深的提速方案，例如按图规模分桶的 batch sampler 或 static edge cache，适合作为新的 speed_v3 实验单独验证，不混入当前 Q1-EXP001 正式 baseline。

## 2026-05-27 EXP002 结构数据准备

### 当前结论

服务器上的 ESIBank A-only 官方数据包只有处理后的 CSV、LMDB、特征缓存和 split 文件，没有找到原始 PDB、mmCIF、SDF 或 docking 输入文件。因此，无法直接从官方 A-only 数据包重建“原论文作者当时使用的原始结构输入”。

为了先推进 EXP002 的 Fe/HEM 结构准备，当前采用 RCSB 候选路线：

1. 从 ESIBank A-only 处理后数据里提取 `structure_index -> UniProt -> substrate` 映射。
2. 用 UniProt 查询 RCSB 可用结构。
3. 从 RCSB 命中里筛出带 HEM/Fe 标记的条目。
4. 下载这些条目的 mmCIF 文件到服务器。

这条路线能用于后续构建“有 Fe/HEM 信息的候选结构数据”，但不等同于原论文作者当时使用的官方原始 docking 结构。后续汇报时要保守说明。

### 数据位置

服务器结构文件目录：

- `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/raw_mmcif_rcsb_heme_fe_v1`

服务器 manifest 目录：

- `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/manifests`

本地已同步 manifest：

- `D:\EZSpecificity_Project\.q1_exp002_manifests\manifests_20260527`

本地说明文件：

- `D:\EZSpecificity_Project\.q1_exp002_manifests\EXP002_data_status_20260527.md`

### 下载和审计结果

| 项目 | 数值 |
|---|---:|
| 计划下载 HEM/Fe RCSB 条目 | `1120` |
| 服务器已下载 `.cif` 文件 | `1120` |
| 缺失文件 | `0` |
| 空文件 | `0` |
| 总大小 | `1.733 GiB` |
| 快速 HEM/Fe token 检查通过 | `1120` |

审计文件：

- `exp002_rcsb_heme_fe_download_audit_v3_20260527.csv`
- `exp002_rcsb_heme_fe_download_audit_v3_20260527.json`

审计表记录每个 RCSB 条目的文件路径、是否存在、字节数、SHA256 和 HEM/Fe token 快速检查结果。

### 后续处理建议

1. EXP002 训练前，需要决定 RCSB 候选结构如何映射回 ESIBank 样本。
2. 如果后续脚本需要 PDB 格式，可以从已下载的 mmCIF 转换，不必重新下载。
3. 如果要做严格复现实验，还需要继续寻找原论文作者使用的官方原始结构包或 docking 输入文件。

### 样本到结构的候选映射

已生成“ESIBank 结构记录 -> RCSB HEM/Fe 结构”的候选映射：

| 项目 | 数值 |
|---|---:|
| ESIBank A-only 结构记录 | `317577` |
| 涉及 UniProt | `7741` |
| 能映射到 HEM/Fe RCSB 结构的记录 | `7184` |
| 能映射到 HEM/Fe RCSB 结构的 UniProt | `145` |
| 候选映射行数 | `77613` |
| 每条结构记录 best1 映射行数 | `7184` |
| 样本覆盖比例 | `2.2621%` |

新增文件：

- `exp002_structure_records_to_rcsb_heme_fe_candidates_v3_20260527.csv`
- `exp002_structure_records_to_rcsb_heme_fe_best1_v3_20260527.csv`
- `exp002_structure_records_to_rcsb_heme_fe_summary_v3_20260527.json`

`candidates` 文件保留一个 UniProt 对应多个 RCSB 结构的全部候选。`best1` 文件按参考序列覆盖度、实体序列覆盖度和 entry id 选择每条 ESIBank 结构记录的一个优先候选，方便后续先跑通流程。

## 2026-05-26 进展：EXP001 第一次正式训练已停止，原因是数据供给瓶颈明显

### 结论

`Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453` 已手动停止。它只跑到 epoch 0 的训练阶段，未完成第一个 epoch，未进入验证阶段，因此没有有效的 val AUC/AUPR，也没有 checkpoint。

这次停止不是因为模型报错或数据校验失败，而是因为训练效率不达标。GPU 利用率出现大量深谷，服务器证据显示 GPU 在等 CPU 端准备数据。继续跑 200 epoch 会浪费资源，也不适合作为 Q1-EXP001 的正式结果。

### 停止前状态

| 项目 | 状态 |
|---|---|
| 运行名 | `Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453` |
| 停止时间 | 2026-05-26 约 19:03 |
| 训练阶段 | epoch 0，约 76% |
| 验证阶段 | 尚未进入 |
| metrics.csv | 只有表头，无有效指标 |
| checkpoint | 未生成 |
| GPU 停止后状态 | 利用率 0%，显存 0 MiB |

保留文件：

- 训练日志：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/logs/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453.log`
- GPU 监控：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/logs/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453_gpu.csv`
- 结果目录：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/results/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453`

### 为什么判断是数据供给瓶颈

只读检查得到的证据：

| 证据 | 说明 |
|---|---|
| GPU 利用率采样有 `96% -> 37% -> 0% -> 0% -> 96%` 这类波动 | GPU 有时能吃满，但经常等下一批数据 |
| 显存稳定在约 `25737 MiB / 32607 MiB` | 显存占用稳定不代表 GPU 一直在计算 |
| 主训练进程约 `899% CPU`，4 个 DataLoader worker 各约 `100% CPU` | CPU 端正在持续处理数据 |
| `vmstat` 显示 `wa=0` | 不是磁盘等待导致 |
| 系统可用内存约 700 GiB | 不是内存不足导致 |
| 日志显示仍在 epoch 0 训练，约 `3.96 batch/s` | 不是验证阶段低利用率 |

代码证据：

- `scripts/pt_dataset.py` 第 410-416 行会按样本读取 `sample_xxxxxx.pt`。
- `scripts/pt_dataset.py` 第 545-597 行每次 `__getitem__` 会重新构造图对象、重建边特征、加载 ESM 和底物特征。
- `scripts/main_training_pt.py` 第 167-178 行使用 PyG 的 DataLoader 拼 batch。

因此，目前的瓶颈主要在 Python/PyG 数据准备链路：每个样本单独读取、张量类型整理、图对象构造、边特征重建、再拼成 batch。这个流程能保证正确性，但无法持续喂满 RTX 5090。

### 下一步建议

1. 不再把这次 `per_sample_b40_w4_pf2_full` 当作正式训练结果。
2. 保留 per-sample cache，因为它证明了数据可用、索引可对齐，但它还不是高效训练形态。
3. 下一步应做 Q1-EXP001 的高效缓存版本：把训练时反复做的图对象重建和小文件读取尽量前置，生成更适合训练的 batch-level 或 shard-level cache。
4. 在新缓存上先跑短测：确认 GPU 利用率、batch/s、显存和第一个验证指标，再决定是否正式跑满。
5. 4 卡训练应等数据供给瓶颈解决后再开。当前形态上 4 卡可能只是让多个 GPU 一起等数据，速度不一定线性提升。
| 2026-05-16 | Q7 完成，实验设计锁定，改造点列出（3 行代码） |

## 2026-05-26 进展：EXP001 原版基线训练已启动

### 当前结论

Q1 先跑 A 版，也就是原模型、原特征维度、原 ESIBank A-only 数据。目标是得到一个可复现的原版基线，后面再和 B 版 Fe/HEM 补丁模型做同数据、同划分、同训练设置的对照。

由于直接从大 PT shard 随机读取时 GPU 利用率波动很大，本次新增了一个 per-sample graph PT cache。这个 cache 只改变读取形态，把原来 shard 里的 graph 样本拆成单样本 `.pt` 文件；实体特征、样本顺序、标签、索引和训练逻辑不改变。

### 已完成检查

| 检查项 | 结果 |
|---|---|
| per-sample cache 生成 | 已完成 |
| train 样本数 | 229855 |
| val 样本数 | 22952 |
| test 样本数 | 53588 |
| train/val/test 抽样校验 | 均通过 |
| GPU 启动前状态 | 空闲，显存 0 MiB |
| 早停策略 | 训练脚本自带 EarlyStopping，监控 `auc/val`，patience=15 |

服务器关键路径：

- 原始 PT cache：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909`
- per-sample PT cache：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1`
- 转换报告：`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1/reports/per_sample_graph_conversion_report.json`
- 校验结果：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/results/pt_cache_validation/per_sample_graph_v1_{train,val,test}.json`

### 正式训练参数

| 参数 | 值 |
|---|---|
| 实验 | Q1-EXP001，ESIBank A-only 原版基线 |
| 运行名 | `Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453` |
| cache | per-sample graph PT cache |
| edge mode | `legacy_bug` |
| batch size | 40 |
| gradient accumulation | 2 |
| effective batch size | 80 |
| num workers | 4 |
| prefetch factor | 2 |
| max epochs | 200 |
| early stopping | `auc/val`，patience=15 |
| GPU | 1 张 RTX 5090 |
| shutdown | 未启用 |

输出路径：

- 启动脚本：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/scripts/run_Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453.sh`
- 训练日志：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/logs/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453.log`
- GPU 监控：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/logs/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453_gpu.csv`
- 结果目录：`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP001_esibank_aonly_original_baseline/results/Q1_EXP001_per_sample_b40_w4_pf2_full_20260526_184453`

### 启动后观察

训练已正常进入 epoch 0。约 5 分钟后检查，进度到 `1162/5747` batch，约 20%，速度约 `3.89 batch/s`。GPU 监控采样 62 次，平均利用率 `37.65%`，峰值 `95%`，显存峰值 `24767 MiB / 32607 MiB`。

这个结果说明 per-sample cache 已经明显改善大 shard 随机读取造成的训练阻塞。GPU 曲线仍会有低谷，后续需要继续看完整 epoch 的验证阶段和 checkpoint 阶段表现。

### EXP002 准备判断

当前服务器的 A-only 官方数据源位于：

`/root/autodl-tmp/EZSpecificity/PathD/data_sources/ESIBank_official_A_only_20260526`

已检查该目录，当前没有 `.pdb`、`.cif`、`.sdf` 或结构压缩包。它主要包含 CSV、LMDB、NPY/NPZ 等已加工特征。

这对 EXP002 很关键。Fe/HEM 补丁需要让 HETATM/Fe 进入结构图。如果只有已经解析好的 `structure_features.lmdb`，原解析阶段丢掉的 HETATM/Fe 通常无法再补回来。后续准备 EXP002 时，应优先确认官方是否有原始结构文件包，或者是否需要重新下载包含 PDB/pocket 的数据源。

## 2026-05-26 晚：Q1-EXP001 GPU 利用率诊断与参数短测

### 当前结论

Q1-EXP001 现在还没有达到“GPU 长时间稳定 95% 以上”的状态。已经验证的现象是：GPU 峰值可以到 95%-97%，显存也能吃到 25-29 GiB，但训练过程中仍会频繁掉到低利用率。当前瓶颈主要出现在数据准备到 PyG batch 拼接这条链路，单纯继续增大 batch size 或 worker 数已经没有明显空间。

截至 20:13，服务器上没有 Q1 训练进程，GPU 利用率 0%，显存 0 MiB。下面这些都是短测或诊断运行，不作为正式实验结果。

### 已确认的事实

| 检查项 | 结果 | 判断 |
|---|---:|---|
| `batch_size=48, workers=4, prefetch=4` | OOM | 32 GiB 5090 放不下 |
| `batch_size=44, workers=4, prefetch=4` | OOM | batch 继续增大空间很小 |
| `batch_size=42, workers=4, prefetch=4` | 完成 400 train batch + 20 val batch | 速度不优于当前最佳设置 |
| `batch_size=40, workers=4, prefetch=4` | 完成 400 train batch + 20 val batch | 稳定，但 GPU 仍有明显低谷 |
| `batch_size=40, workers=8, prefetch=4` | 训练中可到约 4.67 batch/s，但验证阶段报 `received 0 items of ancdata` | 当前服务器 `ulimit -n=1024` 限制文件描述符，worker/prefetch 过高会不稳 |
| `batch_size=40, workers=6, prefetch=2` | 完成 | 没有超过 `workers=8, prefetch=1` |
| `batch_size=40, workers=8, prefetch=1` | 完成 | 当前已知最快且安全的组合 |

当前最好的安全短测：

| 项目 | 数值 |
|---|---:|
| 运行名 | `Q1_EXP001_perf_pathB_like_b40_w8_pf1_acc1_20260526_200739` |
| train batch | 400 |
| val batch | 20 |
| batch size | 40 |
| workers | 8 |
| prefetch | 1 |
| gradient accumulation | 1 |
| 训练段末尾速度 | 约 `4.42 batch/s` |
| 全流程平均 | `0.262 s/batch` |
| 短测 val AUC | `0.5919` |
| GPU 采样峰值 | 93%左右 |
| 显存峰值 | 约 28.3 GiB |

### 和 PathB 的对照

PathB 里确实有较好的单卡记录：`batch_size=48, workers=4` 在单张 4090 上约 `3.74-3.79 it/s`，日志称 GPU 可到 100%。但 Q1-EXP001 不能直接照搬这个设置，原因有三个：

1. 当前 Q1 A-only 数据规模更大，train 有 `229855` 个样本。
2. 当前 5090 上 `batch_size=48` 和 `batch_size=44` 都已经 OOM。
3. PathB 启动脚本里会尝试 `ulimit -n 65536`，当前 AutoDL 容器里实际是 `1024`，高 worker + 高 prefetch 会触发文件描述符相关错误。

PathB 也记录过多卡训练的锯齿问题，原因和变长图样本、batch 难度不均、进程等待有关。多卡可以缩短总训练时间，但不能保证每张卡都稳定 95%。

### 下一步判断

当前可用于正式跑的保守组合是：

| 参数 | 建议值 |
|---|---:|
| cache | per-sample graph PT cache |
| batch size | 40 |
| workers | 8 |
| prefetch | 1 |
| train_in_order | false |
| gradient accumulation | 2 用于保持之前 effective batch 80；若优先测速度可用 1 |

如果继续追求 GPU 利用率，应新建一个小范围性能优化版本，优先减少 `scripts/pt_dataset.py` 里每次 `__getitem__` 的 CPU 工作，例如把可确定的 edge 静态部分、dtype 整理、重复构造步骤提前固化到 cache。这个方向需要做等价性验证，尤其要保证 `edge_mode=legacy_bug`、样本键、标签和训练划分完全一致；训练时的 `dist_noise` 也要保留或明确记录是否关闭。

### `--preload` 额外短测

随后测试了现有脚本自带的 `--preload` 参数：

| 项目 | 数值 |
|---|---:|
| 运行名 | `Q1_EXP001_perf_pathB_like_b40_w8_pf1_acc1_preload_20260526_201857` |
| batch size | 40 |
| workers | 8 |
| prefetch | 1 |
| train batch | 400 |
| val batch | 20 |
| train preload | `229855` 个样本，约 `33260 MB`，耗时 `226.9s` |
| val preload | `22952` 个样本，约 `3323 MB`，耗时 `31.0s` |
| 训练段末尾速度 | 约 `4.57 batch/s` |
| 包含 preload 后的总耗时 | `373s` |
| 短测 val AUC | `0.5925` |

这个结果说明：预加载可以减少一部分小文件读取开销，但启动代价很大，GPU 利用率锯齿仍然存在。现在更明确的瓶颈是每个样本进入模型前的 CPU 构造流程，包括 dtype 转换、mask 构造、边距离计算、GaussianSmearing、PyG DataLoader 拼 batch，以及不同样本图大小差异造成的 batch 计算量波动。

## 2026-05-26 单卡提速定位更新

### 当前结论

Q1-EXP001 单卡慢，主要由两件事叠加造成：

1. Q1 A-only 的训练样本数是 `229855`，在 `batch_size=40` 时每个 epoch 约 `5747` 个 batch。之前 PathD/Q2 常用训练集约 `2.2万` 个样本，batch 数量只有几百级。因此 Q1 一个 epoch 变长，首先是 batch 数量变多。
2. 单 batch 内确实存在 CPU 等待。手写训练循环拆分后，`batch_size=40, workers=8, prefetch=1` 下，每个 batch 平均约：

| 阶段 | 平均耗时 |
|---|---:|
| DataLoader 等待 batch | `104.96 ms` |
| batch 搬到 GPU | `9.52 ms` |
| GPU 前向 + 反向 + 优化器 | `91.65 ms` |
| 总 step | `206.13 ms` |

这说明 GPU 会等 DataLoader，但 GPU 计算本身也不是零成本。单卡想长期稳定 100% 并不现实，除非能继续增大 batch；目前 32 GiB 显存下 `batch_size=44/48` 已经 OOM。

### 完整数据管线 profile

对 `PtCacheDataset.__getitem__` 抽样 200 个 train 样本，结果如下：

| 环节 | 平均耗时 |
|---|---:|
| `__getitem__` 总耗时 | `15.414 ms/样本` |
| 读取 per-sample `.pt` | `0.863 ms/样本` |
| 重建 protein_x | `0.102 ms/样本` |
| 重建 ligand_x | `0.104 ms/样本` |
| 重建 edge feature | `10.893 ms/样本` |
| 读取 ESM enzyme embedding | `2.727 ms/样本` |
| 读取 substrate features | `0.583 ms/样本` |

因此，单样本 CPU 侧最大开销是 edge feature 重建，不是磁盘读取。

### `speed_v1` 轻量改写测试

新增了不覆盖原脚本的临时 speed 版本：

- `scripts/pt_dataset_speed_v1.py`
- `scripts/main_training_pt_speed_v1.py`
- `scripts/q1_speed_v1_smoke.sh`

这个版本只优化 edge feature 的张量构造和距离噪声生成方式。微基准里，带噪声 edge 重建约从 `11.684 ms/样本` 降到 `9.342 ms/样本`。但放回完整训练短测后，收益很小：

| 测试 | 平均 step | 估算 train epoch |
|---|---:|---:|
| 原始 `pt_dataset`，workers=8 | `206.13 ms/batch` | `1184.6 s` |
| `speed_v1`，workers=8 | `199.04 ms/batch` | `1143.9 s` |
| 原始 `pt_dataset`，workers=12 | `198.55 ms/batch` | `1141.1 s` |

正式 Lightning 短测中，`speed_v1` 400 train batch + 20 val batch 仍约 `103s`，与原始版本基本一致。这个版本不足以作为主要提速方案。

### static edge cache 可行性

抽样 1000 个真实样本估算 static edge cache 空间：

| 方案 | 估算空间 |
|---|---:|
| 完整 static edge，float32 | `236.1 GiB` |
| 紧凑版，int64/float32 | `118.1 GiB` |
| 紧凑版，int32/float16/uint8 | `54.1 GiB` |

当前服务器 `/root/autodl-tmp` 剩余约 `64G`。因此完整 static edge 不现实；紧凑全量非常贴边；更现实的是先只对 train 做紧凑 static edge cache，预计约 `41 GiB`，并且必须做成分片缓存，不能做成几十万个额外小文件。

### 下一步建议

单卡下真正值得继续做的提速方案是：

1. 新建 train-only 紧凑 static edge cache，保留原数据不动。
2. 新建 `pt_dataset_static_edge_v1.py`，训练时读取静态 edge_index、基础距离、cross 标记和 bond 类别，只在训练时继续生成 `dist_noise` 并计算 GaussianSmearing。
3. 先用 400 batch 短测比较 GPU 利用率和 `s/batch`，再决定是否用于正式训练。

如果目标是显著缩短 wall-clock，一个更现实的方向是 4 卡 DDP。多卡不能从根上消除每张卡的利用率波动，但可以把每个 epoch 的 batch 分摊到多张卡上；当前脚本已经有 `--devices` 和 DDP 分支，后续需要在 4 卡资源上单独做 DDP 冒烟测试。

## 2026-05-26 单卡提速补充排查

这轮继续排查了几个没有改变数据划分、标签和模型结构的单卡方案。

### 1. 关闭 grad norm 日志

新增脚本：

- `scripts/main_training_pt_nogradnorm_v1.py`
- `scripts/q1_nogradnorm_v1_smoke.sh`

目的：原脚本每次 optimizer step 会计算梯度范数，并通过 `.item()` 读回 CPU。这个操作可能导致 GPU 同步。

结果：

| 项目 | 数值 |
|---|---:|
| 400 train batch + 20 val batch | `103s` |
| 平均 | `0.259s/batch` |
| GPU 平均利用率 | `35.20%` |
| GPU 峰值 | `97%` |

判断：关闭 grad norm 日志没有明显提速。

### 2. `torch.compile`

目的：测试模型侧编译能否减少 GPU 前向/反向耗时。

结果：

| 项目 | 原始手写循环 | `torch.compile` |
|---|---:|---:|
| DataLoader 等待 | `104.96ms` | `114.38ms` |
| GPU 前向+反向+优化器 | `91.65ms` | `86.13ms` |
| 总 step | `206.13ms` | `210.09ms` |
| 估算 train epoch | `1184.6s` | `1207.4s` |

编译过程中出现 PyG/`torch_scatter` 相关 graph break 和反复重编译提示。GPU 计算略降，但总耗时没有下降。

判断：`torch.compile` 不适合当前动态图 PyG 训练路径。

### 3. 手写最小 PyTorch 训练循环

目的：绕开 Lightning 的框架层开销，测试训练主循环本身的上限。

400 train batch 结果：

| 项目 | 数值 |
|---|---:|
| DataLoader 等待 | `109.42ms` |
| GPU 前向+反向+优化器 | `97.62ms` |
| 总 step | `216.57ms` |
| 估算 train epoch | `1244.6s` |
| GPU 平均利用率 | `41.62%` |
| GPU 峰值 | `94%` |

判断：手写循环能减少一部分验证和回调开销，但训练主路径没有数量级提升。要把它改成正式训练脚本，还需要复刻验证、checkpoint、学习率调度、早停和结果记录，收益不够大。

### 4. 增大 enzyme embedding cache

新增脚本：

- `scripts/pt_dataset_enzcache256_v1.py`

目的：当前每个 worker 只缓存 64 个 enzyme embedding，Q1 样本多，可能反复读 ESM。测试把 enzyme cache 提到 256。

结果：

| 项目 | 原始 | enzyme cache 256 |
|---|---:|---:|
| DataLoader 等待 | `104.96ms` | `96.48ms` |
| GPU 前向+反向+优化器 | `91.65ms` | `92.51ms` |
| 总 step | `206.13ms` | `198.51ms` |
| 估算 train epoch | `1184.6s` | `1140.8s` |

判断：这个方案有小幅改善，但每个 worker 会占用更多 CPU 内存。它可以作为正式脚本的小优化，但不能解决 GPU 平均利用率偏低的问题。

### 5. edge 轻改 + enzyme cache 256 合并

新增脚本：

- `scripts/pt_dataset_speed_cache_v2.py`

结果：

| 项目 | 数值 |
|---|---:|
| 总 step | `206.88ms` |
| 估算 train epoch | `1188.9s` |

判断：两个小优化没有稳定叠加收益，不建议作为正式方案。

### 补充结论

目前已排除的单卡小方案：

| 方案 | 结论 |
|---|---|
| 继续调 batch/worker/prefetch | 已接近边界，收益很小 |
| `--preload` | 启动耗时大，GPU 锯齿仍在 |
| edge 重建轻改 | 微基准有效，端到端无明显提升 |
| 关闭 grad norm 日志 | 无明显提升 |
| `torch.compile` | graph break 多，总耗时不降 |
| 手写训练循环 | 只减少少量框架开销 |
| enzyme cache 256 | 小幅改善，不能根治 |

单卡可用的小优化是 `pt_dataset_enzcache256_v1.py` 这一类更大的实体缓存，预期只带来约 3% 到 4% 的 train step 改善。它不能把 GPU 平均利用率拉到长期 90% 以上。

现在更可靠的判断是：Q1-EXP001 如果要明显缩短训练总时间，需要做 4 卡 DDP 冒烟测试。单卡路径已经没有发现值得继续投入的大幅提速方案。

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |
