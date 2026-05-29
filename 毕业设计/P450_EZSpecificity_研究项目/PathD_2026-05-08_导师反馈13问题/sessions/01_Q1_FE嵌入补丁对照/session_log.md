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

## 2026-05-27 EXP003 全 P450 版本数据准备启动

### 当前任务

EXP002 正在服务器上继续正式训练。训练期间不再改 EXP002 的脚本、参数和结果目录，只做只读监控。

同时启动 EXP003 的数据准备。EXP003 暂定为“全 P450 版本”：入口是本地 389 个 P450 UniProt 清单，先找出这些 P450 在当前 ESIBank A 版训练数据中对应哪些样本，再判断哪些样本能沿用 EXP002 的 Fe/HEM 结构增强规则。

### 数据入口

| 类型 | 路径 | 作用 |
|---|---|---|
| 389 个 P450 清单 | `D:\EZSpecificity_Project\毕业设计\提取P450过程日志\2026-01-02_01-46_P450精确验证\数据\P450酶列表_最终版389个.csv` | EXP003 的 P450 定义入口 |
| ESIBank 结构记录索引 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/manifests/exp002_structure_index_to_uniprot_substrate_20260527_014800.csv` | 把 UniProt 映射到 ESIBank structure_index / substrate |
| EXP002 HEM/Fe 审计表 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/manifests/exp002_best1_with_mmcif_atom_audit_atom_audit_v1_20260527.csv` | 判断哪些结构有可用 RCSB HEM/Fe 证据 |
| 官方 ESIBank PDB manifest | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/manifests/exp002_google_drive_pdb_manifest_for_rcsb_best1_20260527.csv` | 判断对应官方 PDB 是否已准备好 |
| EXP001 A 版 split | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/normalized_splits/` | 找到 P450 样本属于 train / val / test 的位置 |
| EXP001 PT cache | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp001_esibank_aonly_baseline/pt_cache/full_legacy_cuda_20260526_161909_per_sample_graph_v1/` | 检查样本是否已有可训练缓存 |

### 已生成脚本

服务器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/scripts/q01_exp003_prepare_full_p450_manifest_20260527.py`

本地：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\scripts\q01_exp003_prepare_full_p450_manifest_20260527.py`

脚本只生成 manifest，不修改 PT cache、训练脚本、检查点或原始数据。

### 输出目录

服务器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/manifests/`

本地：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\manifests\`

最新一版 tag：`20260527_214331`。

### 关键结果

| 项目 | 数量 | 解释 |
|---|---:|---|
| P450 清单行数 | 389 | 来自 `P450酶列表_最终版389个.csv` |
| 唯一 UniProt | 389 | 没有重复 UniProt |
| 当前 ESIBank 结构记录总数 | 317577 | 用于查 P450 是否在现有训练数据里出现 |
| 命中 ESIBank 结构记录的 P450 样本 | 8448 | 这是 EXP003 全 P450 样本的初始集合 |
| 命中 ESIBank 结构记录的 P450 UniProt | 216 | 389 个 P450 中有 216 个在当前数据里出现 |
| P450 正样本 | 586 | 在 8448 条 P450 样本里，label 为正的样本 |
| 没有命中 ESIBank 结构记录的 P450 | 173 | 389 个 P450 中有 173 个暂时没有进入当前训练数据 |
| 有可用 Fe/HEM overlay 证据的样本 | 1152 | 同时满足 RCSB HEM/Fe 审计可用、官方 PDB 已准备好 |
| 有可用 Fe/HEM overlay 证据的 UniProt | 22 | 说明直接按 EXP002 方式补结构级 Fe/HEM 只能覆盖一小部分 P450 |
| 可训练缓存存在的 P450 样本 | 8229 | 8448 条中有 219 条在当前 PT cache 里缺失 |

### train / val / test 分布

| split | 全 P450 样本 | 正样本 | 负样本 | 正样本比例 | 唯一 UniProt | cache 缺失 |
|---|---:|---:|---:|---:|---:|---:|
| train | 6428 | 433 | 5995 | 0.067362 | 216 | 165 |
| val | 621 | 47 | 574 | 0.075684 | 171 | 14 |
| test | 1399 | 106 | 1293 | 0.075768 | 203 | 40 |

去掉 cache 缺失后，可直接训练的样本数：

| split | 可训练样本 |
|---|---:|
| train | 6263 |
| val | 607 |
| test | 1359 |

### 最新 manifest 文件

| 文件 | 作用 |
|---|---|
| `exp003_full_p450_manifest_summary_20260527_214331.json` | 总结输入、输出和核心数量 |
| `exp003_p450389_split_sample_manifest_20260527_214331.csv` | 8448 条全 P450 样本，带 split、sample_id、Dock Index、UniProt 和 cache 检查 |
| `exp003_p450389_trainable_split_sample_manifest_20260527_214331.csv` | 去掉 cache 缺失后的 8229 条可训练样本 |
| `exp003_p450389_split_summary_20260527_214331.csv` | train / val / test 数量表 |
| `exp003_fe_overlay_target_candidates_20260527_214331.csv` | 1152 条能直接沿用 EXP002 Fe/HEM overlay 的样本 |
| `exp003_p450389_uniprot_summary_20260527_214331.csv` | 389 个 P450 每个 UniProt 的命中情况 |
| `exp003_missing_from_esibank_structure_records_20260527_214331.csv` | 173 个没有出现在当前 ESIBank 结构记录中的 P450 |
| `exp003_missing_from_fe_overlay_targets_20260527_214331.csv` | 367 个暂时没有直接 Fe/HEM overlay 证据的 P450 |

### 当前判断

已验证事实：

1. EXP003 的全 P450 初始样本集合是 8448 条，来自 216 个 P450 UniProt。
2. 当前 PT cache 可直接训练的 P450 样本是 8229 条。
3. 直接沿用 EXP002 的结构级 Fe/HEM overlay 规则，只能覆盖 1152 条样本、22 个 UniProt。

基于证据的判断：

1. 若 EXP003 目标是“全 P450 样本训练”，下一步可以先用 8229 条可训练样本构建 EXP003 子集 cache。
2. 若 EXP003 目标是“所有 P450 都补上结构级 Fe/HEM 坐标”，当前证据不足。原因是 367 个 P450 还没有可直接使用的 RCSB HEM/Fe 结构映射。
3. `has_heme` 字段只能作为 UniProt 注释证据，不能直接替代三维坐标；训练图里要加 Fe/HEM 原子，仍然需要结构坐标或明确的模板规则。

子智能体复核：

Franklin 做了只读复核，建议把三个集合分清楚：

1. `p450_uniprot389`：EXP003 的入口，来自 389 个 UniProt。
2. `exp002_fe_overlay_enhanced`：EXP002 中能做结构级 Fe/HEM overlay 的子集。
3. `q2_actual_used`：Q2 的 actual-used 样本集合，只能作为参考或交集审计。

主助手已经按这个建议生成 split 级 manifest 和可训练样本 manifest。

### EXP002 同步监控点

EXP003 manifest 生成完成后，只读检查 EXP002。服务器时间 `2026-05-27 21:46:31`：

| 项目 | 状态 |
|---|---|
| 当前进度 | `epoch 5`，约 `581 / 2874` |
| 当前速度 | 约 `3.25 it/s` |
| GPU 状态 | GPU0 瞬时 0%，GPU1 瞬时 99%；仍有 DDP 等待波动 |
| 最新完整验证轮次 | `epoch 4` |
| epoch 4 train loss | 0.282091 |
| epoch 4 val loss | 0.271798 |
| epoch 4 val AUC | 0.790452 |
| epoch 4 val AUPR | 0.397927 |

EXP002 训练仍在继续，没有因为 EXP003 数据准备中断。

## 2026-05-27 EXP003 P450 子集 cache 与 Fe/HEM overlay 准备

### 本轮目标

在不打断 EXP002 正式训练的情况下，继续准备 EXP003。EXP003 目前按下面边界推进：

1. 入口是 389 个 P450 UniProt 清单。
2. 训练样本先继承 EXP001 A-only baseline 的 train / val / test。
3. 只保留当前 PT cache 里实际存在的 P450 样本。
4. 对其中有可靠结构证据的样本，沿用 EXP002 的规则补入 HEM/Fe 原子。

### 新建 cache

服务器 cache：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1`

本地记录：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\pt_cache_reports\exp003_p450_subset_cache_summary_20260527_2150.json`

构建脚本：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\scripts\q01_exp003_build_p450_subset_cache_20260527.py`

服务器对应脚本：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/scripts/q01_exp003_build_p450_subset_cache_20260527.py`

### cache 构建结果

| split | 样本数 | 正样本 | 负样本 | 样本文件链接 |
|---|---:|---:|---:|---|
| train | 6263 | 431 | 5832 | 6263 个 hardlink |
| val | 607 | 47 | 560 | 607 个 hardlink |
| test | 1359 | 106 | 1253 | 1359 个 hardlink |
| 合计 | 8229 | 584 | 7645 | 8229 个 hardlink |

共享特征：

| 目录 | 文件处理 |
|---|---:|
| `enzymes/` | 3 个 hardlink |
| `substrates/` | 4 个 hardlink |

说明：这些 hardlink 只作为读取来源。后续如果改样本内容，必须用替换文件的方式，不能原地写 hardlink 指向的文件。

### Fe/HEM overlay 结果

已在 EXP003 新 cache 上应用 EXP002 的 Fe/HEM overlay 规则。

服务器报告目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/overlay_reports/fe_overlay_from_exp002_rules_v1`

本地报告目录：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\overlay_reports\fe_overlay_from_exp002_rules_v1`

关键文件：

| 文件 | 作用 |
|---|---|
| `fe_overlay_sample_summary.json` | overlay 脚本输出总结 |
| `fe_overlay_sample_audit_dedup_20260527_2216.csv` | 去重后的 overlay 明细 |
| `fe_overlay_verification_20260527_2216.json` | 实际 `.pt` 文件内容验证 |
| `fe_overlay_audit_dedup_summary_20260527_2216.json` | audit 去重说明 |

overlay 验证结果：

| 项目 | 数量 |
|---|---:|
| Fe/HEM target candidates | 1152 |
| 其中 cache 存在、可处理的样本 | 1140 |
| 实际写入 Fe/HEM 的样本 | 1140 |
| 失败样本 | 0 |
| 缺失 `.pt` 文件 | 0 |
| 写入后缺少 `protein_is_hetero` 的样本 | 0 |
| 写入后找不到 Fe 元素的样本 | 0 |

中途 SSH 守护连接断过一次，原始 `fe_overlay_sample_audit.csv` 里留下 100 条重复记录。已额外生成去重版：

`fe_overlay_sample_audit_dedup_20260527_2216.csv`

去重后是 1140 条唯一 `(split, sample_id, dock_index)` 记录。

### cache 校验

已生成并同步校验报告：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\pt_cache_reports\exp003_cache_index_manifest_verification_20260527_2218.json`

校验结果：

| 检查项 | 结果 |
|---|---|
| train index 行数等于 manifest 行数 | 通过 |
| val index 行数等于 manifest 行数 | 通过 |
| test index 行数等于 manifest 行数 | 通过 |
| `sample_id` 顺序与 manifest 一致 | 通过 |
| `enzyme_id` 与 manifest 一致 | 通过 |
| `substrate_id` 与 manifest 一致 | 通过 |
| 样本文件存在 | 通过 |
| 抽样 `.pt` 内部 label/enzyme/substrate 与 manifest 一致 | 通过 |
| Dataset 读取未增强样本 | 通过 |
| Dataset 读取已增强 Fe/HEM 样本 | 通过 |

复用 EXP002 overlay 脚本后，`manifest.pt` 曾短暂带有 `exp002_fe_overlay=True` 标记。该标记容易误导，已在 EXP003 新 cache 内修正为：

| manifest 字段 | 当前值 |
|---|---|
| `exp003_full_p450_subset` | `True` |
| `exp003_fe_overlay_from_exp002_rules` | `True` |
| `exp003_fe_overlay_sample_count` | `1140` |
| `protein_elements` | `[1, 6, 7, 8, 16, 34, 26]` |
| `protein_has_is_hetero` | `True` |

修正记录：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\data\q01_fe_embedding_patch\exp003_full_p450_fe_heme_overlay\pt_cache_reports\exp003_manifest_flag_correction_20260527_2220.json`

### 当前边界

已验证事实：

1. EXP003 已有独立可训练 cache：8229 条 P450 样本。
2. 其中 1140 条样本已经按 EXP002 规则补入 HEM/Fe 结构原子。
3. 训练脚本可以用现有 `PtCacheDataset` 读取这个 cache。
4. EXP001 原始 cache 未被覆盖；EXP003 目录里被增强的样本是替换后的新文件，未增强样本仍是硬链接读取。

仍需注意：

1. 8229 条样本只覆盖 216 个 P450 UniProt，不等于 389 个 P450 全覆盖。
2. EXP003 当前继承 EXP001 random split，不是酶级隔离 split。
3. 只有 1140 条样本带结构级 Fe/HEM 原子，其余 P450 样本仍是 A-only 原图结构。

### EXP002 同步监控点

服务器时间 `2026-05-27 22:21:53`：

| 项目 | 状态 |
|---|---|
| 当前进度 | `epoch 7`，约 `1199 / 2874` |
| 当前速度 | 约 `3.28 it/s` |
| 最新完整验证轮次 | `epoch 6` |
| epoch 6 train loss | 0.264587 |
| epoch 6 val loss | 0.262130 |
| epoch 6 val AUC | 0.827763 |
| epoch 6 val AUPR | 0.430451 |

EXP002 没有因为 EXP003 cache 构建和 overlay 准备中断。

## 2026-05-27 22:32 EXP003 训练入口已补齐

### 当前状态

EXP003 没有放弃。它现在已经完成了三类准备：

1. 数据清单：从 389 个本地 P450 列表出发，匹配到 ESIBank A 版中可训练的 P450 子集。
2. PT cache：从 EXP001 A-only baseline cache 派生出 EXP003 独立 cache，并对可覆盖样本加入 Fe/HEM 异质原子信息。
3. 训练入口：新增 EXP003 专用配置和启动脚本，后续可以直接开训。

EXP003 训练尚未启动，当前 GPU 仍主要给 EXP002 正式训练使用。

### EXP003 新增训练文件

服务器目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/`

| 文件 | 用途 |
|---|---|
| `configs/config_exp003_p450_subset_fe_overlay.yml` | EXP003 专用配置，`data.tag` 为 `Q1-EXP003-p450389-subset-fe-heme-overlay` |
| `scripts/main_training_pt_fe_overlay.py` | 从 EXP002 复制的 Fe/HEM 训练脚本副本 |
| `scripts/pt_dataset_fe_overlay.py` | 从 EXP002 复制的 PT cache Dataset 副本，可读取 `protein_is_hetero` |
| `scripts/run_exp003_p450_subset_fe_overlay_train.sh` | EXP003 专用启动脚本 |

本地对应文件：

| 文件 | 用途 |
|---|---|
| `configs/q01_exp003_p450_subset_fe_overlay.yml` | 本地保存的 EXP003 配置副本 |
| `scripts/q01_exp003_run_p450_subset_fe_overlay_train.sh` | 本地保存的 EXP003 启动脚本副本 |

### EXP003 启动脚本指向的数据

启动脚本使用的 cache：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1`

当前 cache 的样本数：

| split | 样本数 |
|---|---:|
| train | 6263 |
| val | 607 |
| test | 1359 |

Fe/HEM 覆盖信息：

| 检查项 | 结果 |
|---|---|
| `protein_elements` | `[1, 6, 7, 8, 16, 34, 26]`，包含 Fe |
| `protein_has_is_hetero` | `true` |
| `exp003_fe_overlay_from_exp002_rules` | `true` |
| `exp003_fe_overlay_sample_count` | `1140` |
| 误导性 `exp002_fe_overlay` 标记 | 已移除 |

### 静态验证

已在服务器完成：

```bash
bash -n scripts/run_exp003_p450_subset_fe_overlay_train.sh
/root/miniconda3/bin/python -m py_compile scripts/main_training_pt_fe_overlay.py scripts/pt_dataset_fe_overlay.py
```

已检查 EXP003 配置和启动脚本中没有残留这些容易跑错实验的关键词：

```text
EXP002
exp002_fe_heme_structures
00_EXP002
Q1_EXP002
```

启动脚本中的关键路径已确认：

```text
CACHE_DIR=/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1
CONFIG=configs/config_exp003_p450_subset_fe_overlay.yml
DEVICES=${DEVICES:-1}
RESULTS_DIR=${RESULTS_DIR:-results/00_EXP003_P450_SUBSET_FE_OVERLAY_${STAMP}}
```

### 后续如何启动

等 EXP002 结束或用户确认分配 GPU 后，可以在服务器 EXP003 目录执行：

```bash
DEVICES=1 BATCH_SIZE=40 NUM_WORKERS=8 MAX_EPOCHS=200 bash scripts/run_exp003_p450_subset_fe_overlay_train.sh
```

如果使用多卡，需要先重新确认 GPU 空闲情况，再设置 `CUDA_VISIBLE_DEVICES` 和 `DEVICES`，避免和 EXP002 或其他进程抢卡。

### 子智能体复核意见

Franklin 做了只读复核。主要风险提醒：

1. EXP003 不能残留 EXP002 的 cache 路径、结果目录和实验名。
2. 第一次正式跑不要误续接 EXP002 的检查点。
3. 如果 EXP002还在跑，EXP003不要默认占用多卡。
4. 这次复制的是 EXP002 的 Fe/HEM 训练入口，Dataset 仍应使用 `pt_dataset_fe_overlay.py`，不能换回普通 PT Dataset。

主助手已按这些风险点复核了实际文件。当前结论是：EXP003 的训练入口已经具备，尚未启动训练。

## 2026-05-27 22:36 EXP002 仍在训练，EXP003 暂不抢 GPU

### 服务器当前状态

检查时间：`2026-05-27 22:36:55`

EXP002 仍在运行，命令核心参数如下：

```text
--cache-dir /root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1
--batch-size 40
--max-epochs 200
--devices 2
--results-dir results/00_EXP002_FE_OVERLAY_20260527_184505
--resume last
```

当时 GPU 占用：

| GPU | 利用率 | 显存 |
|---|---:|---:|
| 0 | 75% | 14366 / 32607 MiB |
| 1 | 97% | 18940 / 32607 MiB |

两张卡都被 EXP002 占用。此时启动 EXP003 会和 EXP002 抢计算资源，所以没有启动 EXP003 训练。

### EXP002 最新完整指标

`metrics.csv` 最新完整行到 epoch 7：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 4 | 0.282091 | 0.271798 | 0.790452 | 0.397927 |
| 5 | 0.271018 | 0.271012 | 0.814643 | 0.420826 |
| 6 | 0.264587 | 0.262130 | 0.827763 | 0.430451 |
| 7 | 0.260572 | 0.265369 | 0.840238 | 0.438965 |

当前观察：

1. EXP002 的验证 AUC 和 AUPR 仍在上升。
2. epoch 7 的 val loss 比 epoch 6 高一点，但 AUC/AUPR 继续提高，暂时不能只看 loss 判定变差。
3. 当前最佳验证 AUC 记录为 `0.840238`。

### EXP003 当前动作

EXP003 已具备训练入口，但未启动。后续启动前需要确认：

1. EXP002 已结束，或用户明确要求暂停 EXP002。
2. 至少有一张 GPU 空闲，或者用户确认用多卡重新分配资源。
3. 第一次 EXP003 启动先用新结果目录，避免续接错误检查点。

## 2026-05-27 22:43 EXP003 开训前自检脚本完成

### 为什么加这个脚本

EXP003 是从 EXP002 的 Fe/HEM 训练入口复制出来的，最容易出现两类错误：

1. 启动脚本仍然指向 EXP002 的 cache 或结果目录。
2. cache 里 Fe/HEM 标记、样本数、split 索引和 manifest 不一致。

因此新增一个开训前自检脚本。后续每次正式启动 EXP003 前，先跑它，确认关键路径和数据状态。

### 新增文件

本地：

`scripts/q01_exp003_preflight_check_20260527.py`

服务器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/scripts/q01_exp003_preflight_check_20260527.py`

### 自检覆盖内容

| 检查项 | 结果 |
|---|---|
| EXP003 实验目录存在 | 通过 |
| EXP003 cache 目录存在 | 通过 |
| EXP003 config 存在 | 通过 |
| EXP003 run script 存在 | 通过 |
| EXP003 训练脚本和 Dataset 脚本存在 | 通过 |
| config tag 为 `Q1-EXP003-p450389-subset-fe-heme-overlay` | 通过 |
| run script 无 EXP002 路径或实验名残留 | 通过 |
| run script 指向 EXP003 cache | 通过 |
| run script 写入 `00_EXP003_P450_SUBSET_FE_OVERLAY_*` | 通过 |
| 训练脚本和 Dataset 脚本 Python 编译 | 通过 |
| run script bash 语法 | 通过 |
| manifest 包含 Fe 元素 26 | 通过 |
| manifest 有 `protein_has_is_hetero=true` | 通过 |
| manifest 有 EXP003 Fe/HEM overlay 标记 | 通过 |
| manifest 中 overlay 样本数为 1140 | 通过 |
| manifest 无误导性 `exp002_fe_overlay` 键 | 通过 |
| train/val/test index 数量为 6263 / 607 / 1359 | 通过 |
| 当前没有 EXP003 训练进程 | 通过 |

### 自检报告

服务器报告：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/preflight_reports/exp003_preflight_20260527_224315.json`

本地报告：

`data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/preflight_reports/exp003_preflight_20260527_224315.json`

第一次运行时，自检脚本把自检进程自身误判成 EXP003 进程，已修正为只识别 EXP003 训练进程。第二次运行通过。

### 当前结论

EXP003 训练入口和数据自检已经准备好。当前没有启动 EXP003，因为 EXP002 仍占用两张 GPU。

### 2026-05-27 22:51 EXP002 最新训练状态

检查时间：`2026-05-27 22:51:39`

GPU 状态：

| GPU | 利用率 | 显存 |
|---|---:|---:|
| 0 | 96% | 14366 / 32607 MiB |
| 1 | 55% | 18940 / 32607 MiB |

EXP002 仍在运行，命令仍然是：

```text
--cache-dir /root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1
--batch-size 40
--max-epochs 200
--devices 2
--results-dir results/00_EXP002_FE_OVERLAY_20260527_184505
--resume last
```

`metrics.csv` 最新行到 epoch 8：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 5 | 0.271018 | 0.271012 | 0.814643 | 0.420826 |
| 6 | 0.264587 | 0.262130 | 0.827763 | 0.430451 |
| 7 | 0.260572 | 0.265369 | 0.840238 | 0.438965 |
| 8 | 0.256100 | 0.252495 | 0.834843 | 0.437313 |

同时，checkpoint 目录出现了新的最好模型：

`pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep08-auc0.8554.ckpt`

直接读取 checkpoint 元数据得到：

| checkpoint | ckpt epoch | global step | current `auc/val` | best `auc/val` so far |
|---|---:|---:|---:|---:|
| `ep04-auc0.8146.ckpt` | 4 | 14370 | 0.814643 | 0.814643 |
| `ep05-auc0.8278.ckpt` | 5 | 17244 | 0.827763 | 0.827763 |
| `ep06-auc0.8402.ckpt` | 6 | 20118 | 0.840238 | 0.840238 |
| `ep07-auc0.8348.ckpt` | 7 | 22992 | 0.834843 | 0.840238 |
| `ep08-auc0.8554.ckpt` | 8 | 25866 | 0.855377 | 0.855377 |

当前判断：

1. EXP002 仍在训练，没有完成。
2. 当前 best checkpoint 是 `ep08-auc0.8554.ckpt`。
3. `metrics.csv` 和 checkpoint 文件名/元数据存在一轮左右的写入时机差异；后续报告 best 模型时，以 checkpoint 元数据为准。
4. 训练脚本已包含早停：`EarlyStopping(monitor="auc/val", patience=15)`。
5. 训练脚本训练结束后默认会自动用 best checkpoint 跑 test evaluation，除非启动时使用 `--skip-auto-test`。当前启动命令没有这个参数。

### 2026-05-27 22:56 EXP002 仍在运行，尚未产生 test_eval

检查时间：`2026-05-27 22:56:25`

结果：

| 项目 | 状态 |
|---|---|
| EXP002 训练进程 | 仍在运行 |
| `test_eval.json` | 尚不存在 |
| GPU0 | 利用率 0%，显存 14366 / 32607 MiB |
| GPU1 | 利用率 6%，显存 18940 / 32607 MiB |
| `metrics.csv` 最新更新时间 | `2026-05-27 22:48:07` |
| `last.ckpt` 最新更新时间 | `2026-05-27 22:48:20` |

进一步查看进程状态后，训练没有结束，也没有证据显示已经卡死：

| 进程类型 | 观察 |
|---|---|
| 主训练进程 | CPU 仍有占用 |
| DDP 子进程 | CPU 仍有占用 |
| DataLoader worker | 多个 worker 仍在运行 |
| GPU 显存 | 仍被两个训练进程占用 |

当前判断：

1. EXP002 还没有完成，不能开始 EXP003 抢卡。
2. 当前低 GPU 利用率只是这个检查时刻的瞬时状态，结合 CPU 占用和进程状态看，不能判定为训练结束。
3. 继续等待下一轮 metrics 或自动 test evaluation。

### checkpoint 与 metrics 的记录差异说明

Franklin 只读复核后确认：当前记录“最佳 checkpoint 为 `ep08-auc0.8554`”是合理的。

原因：

1. 训练脚本 `ModelCheckpoint` 明确使用 `monitor="auc/val"`。
2. checkpoint 文件名中的 `auc` 来自 Lightning 监控值。
3. 直接读取 checkpoint 元数据得到 `current_score=0.8553769588`，`best_model_score=0.8553769588`。
4. 自定义 `metrics.csv` 与 checkpoint 元数据存在约一轮验证的错位。

因此，后续记录 best 模型时，以 checkpoint 元数据和最终 `test_eval.json` 为准；`metrics.csv` 主要用于训练曲线参考。

### 2026-05-27 23:03 EXP002 当前结果汇总脚本

新增脚本：

本地：

`scripts/q01_exp002_collect_result_summary_20260527.py`

服务器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/scripts/q01_exp002_collect_result_summary_20260527.py`

用途：

1. 读取 EXP002 的 `metrics.csv`。
2. 读取 checkpoint 元数据，找出当前 best checkpoint。
3. 如果 `test_eval.json` 已存在，同时纳入测试集结果。
4. 写出新增摘要文件，不修改训练进程和原始结果文件。

当前已生成摘要：

服务器：

| 文件 | 用途 |
|---|---|
| `results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_20260527_230320.json` | 当前机器可读摘要 |
| `results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_20260527_230320.md` | 当前人读摘要 |
| `results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_latest.json` | 最新摘要指针 |
| `results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_latest.md` | 最新摘要指针 |

本地：

`results/q01_fe_embedding_patch/EXP002_fe_heme_overlay/result_summary_20260527_230320.json`

`results/q01_fe_embedding_patch/EXP002_fe_heme_overlay/result_summary_20260527_230320.md`

当前摘要内容：

| 项目 | 值 |
|---|---|
| training process still running | `true` |
| test_eval exists | `false` |
| 当前 best checkpoint | `pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep08-auc0.8554.ckpt` |
| checkpoint `auc/val` | `0.8553769588470459` |
| checkpoint epoch | `8` |
| global step | `25866` |

备注：第一次运行脚本时，摘要把 `last.ckpt` 显示为 best 文件。已修正为优先使用 checkpoint 元数据中的 `best_model_path_recorded`，所以当前摘要指向真实 best 文件 `ep08-auc0.8554.ckpt`。

### 2026-05-27 23:12 EXP002 最新状态

检查时间：`2026-05-27 23:11:12`

EXP002 仍在训练，`test_eval.json` 仍未生成。GPU 仍被 EXP002 训练进程占用：

| GPU | 利用率 | 显存 |
|---|---:|---:|
| 0 | 100% | 14370 / 32607 MiB |
| 1 | 45% | 18940 / 32607 MiB |

`metrics.csv` 已更新到 epoch 9：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 9 | 0.249226 | 0.249149 | 0.855377 | 0.474261 |

重新运行结果汇总脚本后，当前 best checkpoint 更新为：

`pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep09-auc0.8599.ckpt`

checkpoint 元数据：

| 项目 | 值 |
|---|---|
| checkpoint epoch | 9 |
| global step | 28740 |
| checkpoint `auc/val` | 0.8599215149879456 |
| best `auc/val` so far | 0.8599215149879456 |
| test_eval exists | false |

当前摘要文件：

服务器：

`results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_20260527_231208.json`

`results/00_EXP002_FE_OVERLAY_20260527_184505/result_summary_20260527_231208.md`

本地：

`results/q01_fe_embedding_patch/EXP002_fe_heme_overlay/result_summary_20260527_231208.json`

`results/q01_fe_embedding_patch/EXP002_fe_heme_overlay/result_summary_20260527_231208.md`

当前判断：EXP002 仍在变好，继续等待早停或训练结束自动测试；暂不启动 EXP003，避免抢占 GPU。

### 2026-05-27 22:48 自检脚本补强

Franklin 只读复核后指出，自检脚本已经能检查路径和 manifest，但还可以补两个点：

1. 报告名用更精细的时间戳，避免同一秒内重复运行时覆盖报告。
2. 加入 Fe/HEM overlay audit 文件存在性和少量样本级抽查。

已采纳并重新运行，新的通过报告为：

服务器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/preflight_reports/exp003_preflight_20260527_224858_134618.json`

本地：

`data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/preflight_reports/exp003_preflight_20260527_224858_134618.json`

新增检查结果：

| 检查项 | 结果 |
|---|---|
| overlay audit CSV 存在 | 通过 |
| overlay audit 行数 | 1140 |
| overlay audit `ok` 行数 | 1140 |
| overlay audit 唯一样本键数 | 1140 |
| 抽查 PT 文件存在 | 通过 |
| 抽查 PT 文件含 `protein_is_hetero` | 通过 |
| 抽查 PT 文件 `hetero_sum` | 5 个样本均为 43 |

抽查样本：

| split | sample_id | hetero_sum |
|---|---:|---:|
| train | 142 | 43 |
| train | 546 | 43 |
| train | 551 | 43 |
| train | 777 | 43 |
| train | 1077 | 43 |

这一步仍然没有启动训练，也没有修改原始数据。脚本只读取 cache、manifest、index、audit CSV 和少量 PT 文件，新增写入的是 EXP003 `preflight_reports/` 下的 JSON 报告。

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |

## 2026-05-28 最新结果索引

详细记录见本文件上方同名章节：`2026-05-28：EXP002 停止训练后四组测试与 3185 条 overlay 解释`。

四组测试结果：

| 模型 | 测试集 | 样本数 | 正样本 | 测试 AUC | 测试 AUPR |
|---|---|---:|---:|---:|---:|
| EXP001 原始基线 | 全量 test | 53588 | 5941 | 0.894870 | 0.569394 |
| EXP001 原始基线 | 389 P450 子 test | 1359 | 106 | 0.890075 | 0.457316 |
| EXP002 HEM/Fe overlay | 全量 test | 53588 | 5941 | 0.893670 | 0.559092 |
| EXP002 HEM/Fe overlay | 389 P450 子 test | 1359 | 106 | 0.890369 | 0.403287 |

当前判断：EXP002 没有超过 EXP001；P450 子测试集上 AUC 接近，AUPR 下降明显。EXP002 当前只能作为全库 HEM/Fe overlay 对照，不能作为严格 P450 专项补丁实验。

`3185` 条 overlay 的含义：

| 分类 | 样本条数 | UniProt 数 |
|---|---:|---:|
| 在 389 P450 清单内 | 1140 | 22 |
| 不在 389 P450 清单内 | 2045 | 32 |
| 合计 | 3185 | 54 |

原因：EXP002 的筛选条件是“ESIBank brenda 样本能找到 HEM/Fe 结构证据并能写回 PT cache”，没有加 389 P450 清单过滤。后续 P450 专项应以 EXP003 为准。

## 2026-05-28：EXP003 当前执行状态与自动收尾设置

### 当前状态

EXP003 已按严格 P450 口径进入正式全量重训练。训练使用全量 ESIBank A 版 brenda cache，模型使用包含 Fe 的原子嵌入表；其中有 raw PDB 的 7770 条 trainable P450 样本已通过 AlphaFill / 同源模板方式补入 HEM/Fe，并完成缓存完整性检查。

截至服务器时间 `2026-05-28 07:06:58`，训练仍在运行，尚未完成 35 epoch，因此还不能写最终结果。

已确认的运行设置：

| 项目 | 当前值 |
|---|---|
| 实验 | `Q1_EXP003_full_af_heme_fe_ddp2_b40_w4_20260528_055549` |
| 结果目录 | `/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/results/00_EXP003_FULL_AF_HEME_FE_20260528_055549` |
| full cache | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/full_exp003_alphafill_heme_fe_v1` |
| max epoch | `35` |
| 早停 | 已禁用，命令包含 `--disable-early-stopping` |
| batch size | `40` |
| GPU 数 | `2` |

当前已写入的阶段性训练记录：

| epoch | train loss | val loss | 阶段性 checkpoint |
|---:|---:|---:|---|
| 0 | 0.351663 | 0.340517 | `ep00-auc0.5942.ckpt` |
| 1 | 0.325542 | 0.313724 | `ep01-auc0.6778.ckpt` |
| 2 | 0.310306 | 0.299656 | `ep02-auc0.7386.ckpt` |

这些只是中途结果，最终必须等 epoch 34 完成后再选最佳 checkpoint。

### 389 P450 子 test cache

EXP003 的 389 P450 子 test cache 已准备完成，来源是 EXP003 full cache，不使用 EXP001 或 EXP002 的缓存。

| 项目 | 数值 |
|---|---:|
| test 样本数 | 1359 |
| 正样本 | 106 |
| 负样本 | 1253 |

路径：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/eval_subsets/p450389_test_subset_20260528_033240/exp003_p450389_test_cache`

### 已新增的自动收尾脚本

为了避免训练结束后无人值守浪费服务器费用，已新增两个脚本。它们只在 EXP003 目录下工作，不修改原始数据。

| 脚本 | 作用 |
|---|---|
| `scripts/q01_exp003_eval_full_and_p450389_20260528.sh` | 训练完成后选择最终最佳 checkpoint，分别跑 full test 和 389 P450 子 test，生成汇总 JSON |
| `scripts/q01_exp003_finish_watch_eval_shutdown_20260528.sh` | 后台监控训练；确认 35 epoch 完成、两个测试 JSON 成功且样本数正确后，自动关机 |

保护条件：

1. 如果 `metrics.csv` 最大 epoch 小于 `34`，不测试。
2. 确认训练已经跑到 epoch `34` 后，再从所有已保存 checkpoint 中按验证 AUC 选择最佳 checkpoint；最佳 checkpoint 不要求一定来自 epoch `34`。
3. 如果 full test 或 389 P450 子 test 的 JSON 缺失，不关机。
4. 如果 389 P450 子 test 样本数不是 `1359`，或正负样本不是 `106 / 1253`，不关机。
5. 只有上述检查全部通过，才执行关机。
6. 测试成功后，watcher 会先把最终测试表追加到服务器 `sessions/q01_fe_embedding_patch/session_log.md`，再延迟 2 分钟关机。

自动收尾日志：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/logs/exp003_finish_watch_eval_shutdown_20260528_070915.log`

最终测试结果预计写入：

`/root/autodl-tmp/EZSpecificity/PathD/P450/results/q01_fe_embedding_patch/exp003_eval_<时间戳>/EXP003_eval_summary.json`

同时会写入指针文件：

`/root/autodl-tmp/EZSpecificity/PathD/P450/results/q01_fe_embedding_patch/EXP003_FINAL_EVAL_POINTER.txt`

## 2026-05-28：EXP002 停止训练后四组测试与 3185 条 overlay 解释

### 当前结论

已停止 EXP002 当前训练进程，并用 EXP001、EXP002 各自最佳检查点完成四组测试：

| 模型 | 测试集 | 样本数 | 正样本 | 测试 AUC | 测试 AUPR |
|---|---|---:|---:|---:|---:|
| EXP001 原始基线 | 全量 test | 53588 | 5941 | 0.894870 | 0.569394 |
| EXP001 原始基线 | 389 P450 子 test | 1359 | 106 | 0.890075 | 0.457316 |
| EXP002 HEM/Fe overlay | 全量 test | 53588 | 5941 | 0.893670 | 0.559092 |
| EXP002 HEM/Fe overlay | 389 P450 子 test | 1359 | 106 | 0.890369 | 0.403287 |

结果说明：EXP002 在全量 test 和 P450 子 test 上都没有超过 EXP001。P450 子 test 上 AUC 基本接近，AUPR 明显下降，因此当前 EXP002 不能作为“加 HEM/Fe 后效果更好”的证据。

### 使用的检查点

| 实验 | 检查点 | epoch |
|---|---|---:|
| EXP001 | `pt-Q1_EXP001_ddp3_nogradnorm_b40_w8_pf1_inorderfalse_resumeep00_full_20260527_0040-ep23-auc0.8965.ckpt` | 23 |
| EXP002 | `pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep24-auc0.8932.ckpt` | 24 |

四组测试输出目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/results/q01_fe_embedding_patch/fourway_eval_20260528_033240/`

汇总文件：

`fourway_eval_summary.json`

### 389 P450 子测试集

子测试集来自 389 个 P450 UniProt 清单与 EXP001/EXP002 test cache 的交集：

| 项目 | 数值 |
|---|---:|
| 样本数 | 1359 |
| 正样本 | 106 |
| 负样本 | 1253 |
| UniProt 数 | 203 |

子测试集目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/eval_subsets/p450389_test_subset_20260528_033240/`

### 3185 条 overlay 的真实含义

`3185` 指 EXP002 中实际成功写入 HEM/Fe overlay 的样本级条目数，也就是 3185 个 `(split, sample_id, dock_index)` 对应的 `.pt` 样本被增强。它不是 3185 个酶，也不是 3185 个 P450。

按 389 个 P450 UniProt 清单重新核对后：

| 分类 | 样本条数 | UniProt 数 |
|---|---:|---:|
| 在 389 P450 清单内 | 1140 | 22 |
| 不在 389 P450 清单内 | 2045 | 32 |
| 合计 | 3185 | 54 |

这说明 EXP002 的 3185 条不是 P450 专项集合。它来自全 ESIBank brenda 中“能找到 HEM/Fe 结构证据并能写回 PT cache”的样本。筛选条件看的是 HEM/Fe 结构可用性，没有加 389 P450 清单过滤。

来源链条：

1. 从 ESIBank A-only 可训练结构记录中提取 `structure_index -> UniProt -> substrate`，共 317577 条结构记录、7741 个 UniProt。
2. 用 UniProt 查 RCSB HEM/Fe 候选结构，得到 7184 条 best1 映射。
3. 审计目标链是否同时包含 HEM 和 Fe，得到 3496 条可用目标。
4. 再要求 Google Drive 官方 `ESIBank/brenda/structure/complex` 中有对应 PDB，并且能成功写回 PT cache，最终成功 overlay 3185 条。

2045 条不在 389 P450 清单内的样本全部来自 `brenda`。按当前 ESIBank 行记录的 EC 一级类统计，多数属于 EC 1 氧化还原相关条目：

| EC 一级类 | 条数 |
|---|---:|
| EC 1 | 1614 |
| EC 2 | 163 |
| EC 3 | 159 |
| EC 4 | 57 |
| EC 6 | 21 |
| EC 7 | 18 |
| EC 5 | 13 |

高频非 389-P450 UniProt 示例：

| UniProt | 条数 | RCSB 结构 | 结构配体标记 |
|---|---:|---|---|
| P04963 | 402 | 7RST | HEM; IOD; MAN; MN; NAG |
| B9W4V6 | 334 | 9JLH | HEM; MG; NAG |
| Q46444 | 290 | 1KB0 | CA; GOL; HEC; PQQ; TFB |
| Q8GR64 | 201 | 1KV9 | ACN; CA; EPE; GOL; HEC; PQQ |
| P05164 | 171 | 6BMT | CA; CL; HEM; NAG |

### 对 EXP002 和 EXP003 的影响

EXP002 当前更适合作为“全库 HEM/Fe overlay 对照”，不适合作为严格的 P450 专项补丁实验。

后续如果继续做老师提到的 P450 补 Fe/HEM 对照，EXP003 应从 389 P450 清单出发，只在 P450 可训练样本中做增强。当前 EXP003 汇总显示，389 P450 可训练样本为 8229 条；其中满足 EXP002 这套 HEM/Fe overlay 条件的候选约 1152 条，实际 overlay 成功为 1140 条。

## 2026-05-27 23:18 EXP003 当前状态核对

本次只读查看服务器，没有改动训练文件、数据文件或结果目录。

### 当前结论

EXP003 已完成数据准备、缓存准备、训练入口和预检查；还没有开始正式训练。

它暂时没有启动的原因是 EXP002 仍在占用两张 GPU。服务器当前只看到 EXP002 的训练进程，没有看到 EXP003 的训练进程，EXP003 的 `results/` 下也还没有正式训练结果目录。

### 服务器现场证据

检查时间：`2026-05-27 23:18:54`

GPU 当前占用：

| GPU | 利用率 | 显存 |
|---|---:|---:|
| GPU0 | 76% | 14370 / 32607 MiB |
| GPU1 | 32% | 18940 / 32607 MiB |

仍在运行的训练命令属于 EXP002：

```bash
scripts/run_exp002_fe_overlay_train.sh
scripts/main_training_pt_fe_overlay.py ... --devices 2 ... --results-dir results/00_EXP002_FE_OVERLAY_20260527_184505 --resume last
```

EXP002 还没有生成最终测试结果：

```text
EXP002_TEST = NO_TEST_EVAL
```

EXP002 最新写入到 `metrics.csv` 的完整 epoch 是 `epoch 9`：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 7 | 0.260572 | 0.265369 | 0.840238 | 0.438965 |
| 8 | 0.256100 | 0.252495 | 0.834843 | 0.437313 |
| 9 | 0.249226 | 0.249149 | 0.855377 | 0.474261 |

当前检查点目录里已有：

```text
last.ckpt
pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep09-auc0.8599.ckpt
pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep08-auc0.8554.ckpt
pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep07-auc0.8348.ckpt
pt-Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505-ep06-auc0.8402.ckpt
```

注意：`metrics.csv` 的 epoch 行和 Lightning 保存的 checkpoint 名称存在一轮左右的记录差异；最终看 EXP002 结果时，应以 checkpoint 元数据和训练结束后的 `test_eval.json` 为准。

### EXP003 已准备好的内容

EXP003 路径：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP003_full_p450_fe_heme_overlay/
```

EXP003 数据缓存：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/pt_cache/p450389_trainable_from_exp001_v1/
```

EXP003 预检查报告：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay/preflight_reports/exp003_preflight_20260527_224858_134618.json
```

预检查已确认：

| 项目 | 数值 |
|---|---:|
| P450 清单中的 UniProt | 389 |
| ESIBank 中能匹配到结构记录的 P450 | 216 |
| EXP003 可训练样本 | 8229 |
| train | 6263 |
| val | 607 |
| test | 1359 |
| Fe/HEM 叠加样本 | 1140 |

EXP003 当前还没有正式结果目录，说明它没有误启动。

### 下一步

等 EXP002 释放 GPU 后，先重新跑一遍 EXP003 预检查；预检查仍通过后，再启动 EXP003 正式训练。

## 2026-05-27 23:30 EXP003 定义纠正

### 纠正结论

用户指出：EXP003 的目标不应是只拿 P450 子集训练。正式 EXP003 应该保留 ESIBank A 版全量训练集，在全量样本中把 389 个 P450 对应的样本尽量补上 HEM/Fe，再从头训练全量模型。

因此，前面已准备的 `8229` 条 P450 子集 cache 不能作为正式 EXP003 开训数据。它只能作为“P450 命中情况和 Fe/HEM 覆盖不足的审计草稿”，后续正式 EXP003 不应启动这个子集训练。

### 已确认的 EXP002 数据事实

服务器只读核对结果：

| 数据层级 | 数量 | 说明 |
|---|---:|---|
| EXP001 A-only baseline PT cache | 306395 | train 229855，val 22952，test 53588 |
| EXP002 Fe/HEM PT cache | 306395 | 仍是全量训练集，样本顺序与 EXP001 一致 |
| EXP002 实际写入 Fe/HEM 的样本 | 3185 | train 2362，val 257，test 566 |
| EXP002 覆盖的 UniProt | 54 | 来自严格 RCSB HEM/Fe 匹配和官方 PDB 可用条件 |
| EXP002 中 P450 相关可覆盖样本 | 约 1140 | 来自 389 P450 清单和 EXP002 strict overlay 规则的交集 |

EXP002 的 `3185` 条不是全 P450 覆盖，而是“全 ESIBank 中满足严格结构证据条件的样本”。它适合做第一版 Fe/HEM 补丁实验，但覆盖太窄。

关键证据路径：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1/reports/full_v2/full_v2_combined_audit_summary.json
/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1/reports/full_v2/exp002_readiness_validation_20260527.json
```

### 已确认的 P450 覆盖事实

389 个 P450 清单核对结果：

| 数据层级 | 数量 | 说明 |
|---|---:|---|
| P450 清单 UniProt | 389 | 来自 `P450酶列表_最终版389个.csv` |
| 在当前 ESIBank 结构记录中命中的 P450 UniProt | 216 | 有样本可映射到当前训练数据 |
| 这些 P450 对应的 ESIBank 样本 | 8448 | 这是“当前训练数据中能看到的 P450 样本” |
| 当前 PT cache 中可直接找到的 P450 样本 | 8229 | 有 219 条缺 PT cache |
| 用 EXP002 严格结构规则能补 Fe/HEM 的 P450 样本 | 1152 候选，1140 已写入 | 覆盖太少 |

因此，正式 EXP003 不能只训练 `8229` 条；它应该训练全量 `306395` 条，只是在其中对 P450 样本做更严格、更完整的 Fe/HEM 处理。

### 正式 EXP003 新定义

正式 EXP003 应该这样定义：

1. 训练集、验证集、测试集仍使用 EXP001/EXP002 的全量 split。
2. 样本总数仍应是 `306395`，不是 `8229`。
3. 非 P450 样本保持原始 A-only 图，不参与 Fe/HEM 修改。
4. P450 样本先识别出所有可训练样本，当前已知约 `8229` 条。
5. 对这些 P450 样本重新设计 Fe/HEM 添加方案，目标是显著超过 EXP002 的 `1140` 条 P450 覆盖。
6. 只有当新的全量 cache、覆盖 manifest、抽样验证、前向验证全部通过后，才启动正式 EXP003 训练。

### 当前 EXP003 状态

正式 EXP003 数据准备状态：未完成。

已经完成且仍有价值的是：

1. 389 个 P450 清单与当前 ESIBank 训练数据的映射审计。
2. P450 样本在 train / val / test 中的分布统计。
3. EXP002 严格 Fe/HEM 规则对 P450 覆盖不足的证据。

需要重新做的是：

1. 不再用 `p450389_trainable_from_exp001_v1` 作为正式训练 cache。
2. 新建一个全量 EXP003 cache，样本数应保持 `306395`。
3. 在全量 cache 内只替换或增强 P450 样本，非 P450 样本保持与 baseline 一致。
4. 新建正式 EXP003 启动脚本和预检查，避免沿用 `p450_subset` 命名和路径。

## 2026-05-27 EXP002 二卡利用率两次尝试记录

### 背景

EXP002 正式训练在二卡模式下仍出现一张 GPU 利用率较高、另一张 GPU 利用率明显偏低的问题。用户限定最多尝试两种方案；若不能明显改善，就恢复正式训练继续跑完。

本次没有修改原始数据，也没有删除正式训练结果。新增的诊断脚本和短跑结果均放在 EXP002 自己的 `scripts/`、`logs/` 和 `results/speed_diag/` 下。

### 先确认的事实

1. 硬件不是直接瓶颈：之前的 synthetic 双卡测试能让两张 5090 同时接近满载。
2. 日志、进度条、`torchrun` 启动方式不是主因：去日志、换 `torchrun` 后，二卡不均衡仍存在。
3. Q2 的 GPU dense table 方案不能直接用于 EXP002 全量训练：预估 train split 的表显存约 `65.5 GB`，单张 5090 只有约 `32 GB`。

GPU dense table 显存预检结果：

| split | samples | enzyme 数 | substrate 数 | 估算表显存 |
|---|---:|---:|---:|---:|
| train | 229855 | 7205 | 30931 | 约 65518 MB |
| val | 22952 | 5205 | 11436 | 约 33223 MB |
| test | 53588 | 6449 | 19197 | 约 47666 MB |

子智能体 Franklin 做了只读复核：Q2 的 GPU dense table 只适合在显存预检通过时尝试；如果表放不下，不应把短跑机会浪费在这个方向。主助手已用服务器当前 cache 做了显存估算，确认 train 表超出单卡显存。

### 尝试 1：SizeBalancedSampler

新增脚本：

| 文件 | 说明 |
|---|---|
| `scripts/main_training_pt_fe_overlay_balanced_try1.py` | 从已有 balanced 训练脚本复制而来，改为使用 `pt_dataset_fe_overlay.py`，保持 Fe/HEM 输入格式 |

运行：

```bash
batch_size=40
devices=2
num_workers=8
limit_train_batches=500
limit_val_batches=1
```

结果：

| 项目 | 结果 |
|---|---|
| run | `Q1_EXP002_try1_balanced_b40_w8_500b_20260527_205041` |
| 训练速度 | 500 个训练 batch 约 3 分钟，和旧短跑基线接近 |
| GPU 状态 | 没有明显改善二卡不均衡 |
| 后续处理 | 训练段完成后进入完整 test，为避免继续占卡，手动停止 |

判断：不采用。

### 尝试 2：balanced_cost batch sampler

新增和生成：

| 文件 | 说明 |
|---|---|
| `scripts/main_training_pt_fe_overlay_speed_v3_balanced_try2.py` | 从已有 speed_v3 balanced 脚本复制而来，改为使用 `pt_dataset_fe_overlay.py` |
| `results/speed_diag/cost_manifests/train_sample_costs_file_size_v1.pt` | 训练样本 cost manifest，用每个 sample `.pt` 文件大小近似图复杂度 |
| `results/speed_diag/cost_manifests/train_sample_costs_file_size_v1.json` | cost manifest 摘要 |

cost manifest 扫描结果：

| 项目 | 数值 |
|---|---:|
| train samples | 229855 |
| missing sample | 0 |
| 文件大小最小值 | 7323 bytes |
| 文件大小最大值 | 580763 bytes |
| 文件大小平均值 | 149247.33 bytes |
| 扫描耗时 | 5.57 秒 |

运行：

```bash
batch_size=40
devices=2
num_workers=8
train_sampler=balanced_cost
limit_train_batches=500
limit_val_batches=1
skip_auto_test=true
```

结果：

| 项目 | 结果 |
|---|---|
| run | `Q1_EXP002_try2_costbalanced_b40_w8_500b_20260527_210405` |
| 训练速度 | rank0 181 秒，rank1 189 秒 |
| 平均速度 | rank0 0.362 s/batch，rank1 0.378 s/batch |
| GPU 状态 | 仍然不均衡，采样时 GPU0 多在 20% 到 30%，GPU1 多在 94% 到 98% |
| 验证 | 只跑 1 个 val batch，不作为模型指标 |

判断：相对旧短跑基线只有很小变化，不能解决 GPU 利用率问题，不采用为正式方案。

### 当前正式训练状态

两次尝试结束后，已恢复 EXP002 正式训练，并从旧正式结果的 `last.ckpt` 续跑。

正式训练目录：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/results/00_EXP002_FE_OVERLAY_20260527_184505
```

恢复日志：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay/logs/Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505_resume_20260527_211005.log
```

续跑确认：

| 项目 | 状态 |
|---|---|
| 恢复 checkpoint | `results/00_EXP002_FE_OVERLAY_20260527_184505/checkpoints/last.ckpt` |
| 恢复前最新 epoch | epoch 2 |
| 恢复后当前 epoch | epoch 3 |
| 正式参数 | `devices=2`，`batch_size=40`，`num_workers=8`，`prefetch_factor=1`，`train_in_order=false`，`max_epochs=200` |

截至恢复前，正式训练已有指标：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 0 | 0.351604 | 0.340488 | 空 | 空 |
| 1 | 0.325845 | 0.313744 | 0.594163 | 0.180194 |
| 2 | 0.316298 | 0.308017 | 0.675017 | 0.309695 |

### 本次判断

两次短跑都没有把二卡利用率拉平。当前最稳妥的做法是继续跑原本 EXP002 正式训练，等待完整曲线和 test 结果。后续若还要根治性能问题，需要单独做更大的代码改造，例如静态 edge cache、真正按每一步 rank 计算量配对的 batch plan，或把模型改成更适合多卡均衡的图 batch 形式。这些已经超出本次“两次机会”的范围。

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

## 2026-05-27 Q1-EXP002 二卡训练利用率诊断

### 当前正式训练

EXP002 正式训练仍在运行，未被中断。

| 项目 | 当前值 |
|---|---|
| 实验目录 | `/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q01_fe_embedding_patch/EXP002_fe_heme_overlay` |
| 运行名 | `Q1_EXP002_fe_overlay_ddp2_b40_w8_full_20260527_184505` |
| 结果目录 | `results/00_EXP002_FE_OVERLAY_20260527_184505` |
| 主脚本 | `scripts/main_training_pt_fe_overlay.py` |
| 启动脚本 | `scripts/run_exp002_fe_overlay_train.sh` |
| 数据缓存 | `data/q01_fe_embedding_patch/exp002_fe_heme_structures/pt_cache/exp002_fe_overlay_v1` |
| 关键参数 | `devices=2`, `batch_size=40`, `num_workers=8`, `prefetch_factor=1`, `train_in_order=false`, `edge_mode=legacy_bug` |

截至检查时，`epoch 1` 已完成并写入 `metrics.csv`：

| epoch | train loss | val loss | val AUC | val AUPR |
|---:|---:|---:|---:|---:|
| 0 | 0.351604 | 0.340488 | 空 | 空 |
| 1 | 0.325845 | 0.313744 | 0.594163 | 0.180194 |

`epoch 0` 的 AUC 在日志中存在，但 CSV 第一行为空，属于记录逻辑问题；`epoch 1` 已正常写入。

### 现场看到的 GPU 情况

20 秒连续采样：

| GPU | 平均利用率 | 范围 |
|---|---:|---|
| GPU0 | 约 30.3% | 0% 到 93% |
| GPU1 | 约 97.0% | 84% 到 100% |

此前也采到过相反方向：GPU0 接近满载、GPU1 较低。低利用率会换卡，说明问题不像某张物理 GPU 损坏，更像两个 DDP rank 之间有等待。

### 已核对的脚本证据

1. `scripts/main_training_pt_fe_overlay.py` 的 DataLoader 只按样本数组织 batch，没有按图大小或边数均衡：

   - `batch_size`
   - `shuffle`
   - `num_workers`
   - `persistent_workers`
   - `prefetch_factor`
   - `in_order`

   这能让两张卡各自拿到相同数量的样本，但不能保证两张卡拿到的总边数接近。

2. `scripts/pt_dataset_fe_overlay.py` 在每个样本读取时仍会做 CPU 侧工作：

   - `torch.load(sample_xxxxxx.pt)`
   - 重建 `protein_x`、`ligand_x`
   - 调用 `rebuild_edge_features`
   - 读取并 padding ESM enzyme embedding
   - 读取并 padding GROVER / Morgan substrate features

3. EGNN 的计算量和边数强相关。当前 Q1 图样本的边数差异很大，抽样中单样本边数大约从 4,702 到 48,984，差异接近 10 倍。

4. 训练回调中每个 batch 都有 `loss.item()`：

   - 位置：`scripts/main_training_pt_fe_overlay.py:253-262`
   - 影响：GPU 计算默认异步，`.item()` 会把 GPU tensor 拉回 CPU，造成同步等待。
   - 脚本自己已经在 `on_before_optimizer_step` 注释里写明 grad norm 的 `.item()` 会强制同步，但 train loss 的 `.item()` 仍保留着。

5. 验证集也用了 `persistent_workers=True`：

   - 位置：`scripts/main_training_pt_fe_overlay.py:174-175`
   - 影响：验证完成后，验证 DataLoader 的 worker 也会留着；每个 worker 有 Dataset 实例和缓存，会增加内存、文件句柄和调度压力。

6. 当前训练开着进度条：

   - 位置：`scripts/main_training_pt_fe_overlay.py:746`
   - 现场日志里同一个 step 会出现两次，说明两个 rank 都在向日志写进度。

### 官方文档对照

查阅 PyTorch / Lightning 官方文档后，目前判断如下：

| 项目 | 官方含义 | 当前脚本情况 | 判断 |
|---|---|---|---|
| DDP 输入分配 | DDP 同步梯度，不会自动按输入计算量切分数据 | 当前按样本数分给 rank | 对图数据不够，无法保证边数均衡 |
| DistributedSampler | 负责让每个进程读不同样本，并让样本数可整除 | Lightning 默认会处理分布式采样 | 只能均分样本数，不知道图边数 |
| `find_unused_parameters=False` | Lightning 建议能关就关，避免额外开销 | 已经是 False | 这里做对了 |
| `gradient_as_bucket_view=True` | 可减少梯度拷贝和显存峰值 | 已经开启 | 这里做对了 |
| `persistent_workers=True` | 可减少每轮重启 worker 的开销，但会让 worker 常驻 | train/val 都开启 | 训练集可以保留，验证集不建议常驻 |
| `in_order=False` | 不强制先进先出，官方提示可能影响可复现性和数据分布 | 当前训练开启 | 需要短跑对照，不应直接当正式默认 |

### 目前结论

这次二卡锯齿不是单一 batch size 问题。当前最可能的主因是：

1. DDP 每个 rank 拿到的样本数一样，但图边数不一样。
2. 边数不同导致每个 rank 的 EGNN 计算时间不同。
3. DDP 每一步需要同步梯度，先算完的卡会等待慢的卡。
4. `loss.item()`、多 rank 进度条、验证 worker 常驻会进一步增加等待和调度开销。

因此，继续只调 `batch_size`、`num_workers`、`prefetch_factor`，很难把两张卡都长期稳定拉到 90% 以上。

### 已新增的短跑诊断脚本

只新增文件，未改动正式训练脚本，未中断当前 EXP002 正式训练。

| 文件 | 用途 |
|---|---|
| `scripts/main_training_pt_fe_overlay_speeddiag.py` | EXP002 诊断版训练脚本 |
| `scripts/run_exp002_speeddiag.sh` | EXP002 诊断版启动脚本 |

诊断版改动：

1. 关闭 Lightning 进度条，减少两个 rank 同时刷日志。
2. `loss.item()` 改成只在 rank0 每 50 个 batch 记录一次，避免每个 batch 都 CPU-GPU 同步。
3. 训练集保留 `persistent_workers=True`，验证集不再常驻 worker。
4. 默认 `TRAIN_IN_ORDER=true`，便于先测试更保守的数据返回顺序。
5. 默认短跑 `500` 个训练 batch 和 `1` 个验证 batch，结果写入 `results/speed_diag/`。

已完成语法检查：

```bash
/root/miniconda3/bin/python -m py_compile scripts/main_training_pt_fe_overlay_speeddiag.py
bash -n scripts/run_exp002_speeddiag.sh
```

### 下一步建议

当前正式 EXP002 可以继续跑。要验证提速方案时，不建议和正式训练并发抢 GPU。合适的做法是：

1. 等当前正式训练到一个可接受检查点，或用户确认暂停。
2. 跑短诊断版：`TRAIN_IN_ORDER=true`，500 个 batch。
3. 再跑对照：`TRAIN_IN_ORDER=false`，500 个 batch。
4. 对比两组的 it/s、20 秒 GPU 平均利用率和日志重复情况。
5. 如果短诊断版仍然一高一低，下一步做 `num_edges` / `num_nodes` 统计表和按图规模分桶的 batch sampler。

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |

## 2026-05-28 最新结果索引

详细记录见本文件上方同名章节：`2026-05-28：EXP002 停止训练后四组测试与 3185 条 overlay 解释`。

四组测试结果：

| 模型 | 测试集 | 样本数 | 正样本 | 测试 AUC | 测试 AUPR |
|---|---|---:|---:|---:|---:|
| EXP001 原始基线 | 全量 test | 53588 | 5941 | 0.894870 | 0.569394 |
| EXP001 原始基线 | 389 P450 子 test | 1359 | 106 | 0.890075 | 0.457316 |
| EXP002 HEM/Fe overlay | 全量 test | 53588 | 5941 | 0.893670 | 0.559092 |
| EXP002 HEM/Fe overlay | 389 P450 子 test | 1359 | 106 | 0.890369 | 0.403287 |

当前判断：EXP002 没有超过 EXP001；P450 子测试集上 AUC 接近，AUPR 下降明显。EXP002 当前只能作为全库 HEM/Fe overlay 对照，不能作为严格 P450 专项补丁实验。

`3185` 条 overlay 的含义：

| 分类 | 样本条数 | UniProt 数 |
|---|---:|---:|
| 在 389 P450 清单内 | 1140 | 22 |
| 不在 389 P450 清单内 | 2045 | 32 |
| 合计 | 3185 | 54 |

原因：EXP002 的筛选条件是“ESIBank brenda 样本能找到 HEM/Fe 结构证据并能写回 PT cache”，没有加 389 P450 清单过滤。后续 P450 专项应以 EXP003 为准。
