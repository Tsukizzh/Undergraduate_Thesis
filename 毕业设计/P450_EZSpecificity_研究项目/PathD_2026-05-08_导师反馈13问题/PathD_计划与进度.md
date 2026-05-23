# PathD — 导师反馈 13 问题推进计划

> **创建日期**: 2026-05-08
> **依据**: 2026-04-22 组会反馈（PPT: 组会4.22_庄泽恒.pdf）
> **指导老师**: 王家祺
> **学生**: 庄泽恒

---

## 一、PathD 定位

PathD 取代原先的 "区域选择性预测" 占位主题（旧目录 `PathD_待定_区域选择性预测/` 已删除），作为 2026-04-22 组会反馈的 **13 个待办问题** 的统一推进容器。

每个问题独立成 session，详见 `sessions/01_Q1_*` ~ `sessions/13_Q13_*`，每个目录下有一份 `session_log.md`，包含问题原文、状态、与已有工作的关系、待澄清点、待办列表、工作量估计。

2026-05-21 起，PathD 的服务器端执行区已经独立建立：

`/root/autodl-tmp/EZSpecificity/PathD/P450`

后续新实验原则上只在该服务器目录下推进。PathC 只作为历史证据区和只读来源，不再作为新实验运行位置。

## 二、目录结构

### 2.1 本地规划目录

```
PathD_2026-05-08_导师反馈13问题/
├── PathD_计划与进度.md              ← 本文件
├── configs/                         ← 共享配置
├── data/00_shared/                  ← 跨问题共享数据
├── scripts/utils/                   ← 共享工具
├── results/
├── logs/
└── sessions/
    ├── 01_Q1_FE嵌入补丁对照/
    ├── 02_Q2_序列相似度划分/        ← 本地记录；服务器工作在 PathD/P450
    ├── 03_Q3_PubChem_3D重对接/
    ├── 04_Q4_EGNN换GVP/
    ├── 05_Q5_数据源版本号与长期网站/
    ├── 06_Q6_正样本数据验证/
    ├── 07_Q7_原论文Heme处理排查/
    ├── 08_Q8_糖基转移酶等家族扩展/
    ├── 09_Q9_底物分类下游预测/      ← 分类部分已在 PathC/C3 完成
    ├── 10_Q10_打分系统/
    ├── 11_Q11_6a15模板重建结构/
    ├── 12_Q12_负样本杂泛性核验/
    └── 13_Q13_FE催化原子距离筛选/
```

本地目录主要承担计划、日志、讨论记录和报告整理功能；真正的大数据、cache、结构文件和训练实验在服务器端 PathD 执行区中维护。

### 2.2 服务器执行目录（2026-05-21 已建立）

```
/root/autodl-tmp/EZSpecificity/PathD/P450/
├── README.md
├── manifests/
├── data/
│   ├── base_from_PathC/
│   │   ├── tables/
│   │   ├── splits/
│   │   ├── features/
│   │   ├── structure/
│   │   ├── cache_best_baseline/
│   │   └── manifests/
│   ├── q01_fe_embedding_patch/
│   ├── q02_sequence_similarity_split/
│   ├── q03_pubchem_3d_redock/
│   ├── q06_q12_sample_validation/
│   ├── q09_substrate_class_prediction/
│   └── q13_fe_catalytic_distance/
├── baselines/
│   └── EXP001_allfix_unified_best/
├── experiments/
├── scripts/
├── sessions/
├── logs/
└── results/
```

服务器端目录规则：

1. `data/` 只放可复用数据资产。
2. `experiments/` 放某次实验的代码快照、配置、日志、结果和临时工作目录。
3. `data/base_from_PathC/` 是从 PathC 复制来的 PathD 起步数据，复制后默认只读使用。
4. 新实验产生的可复用数据进入对应 `data/qXX_*` 目录；一次性中间产物留在对应 `experiments/qXX_*/EXP*/work/`。
5. 不删除、不覆盖 PathC；PathD 内不保留指向 PathC 的符号链接。

## 三、13 问题状态总表

| # | 简称 | 状态 | 与已有工作 | 优先级 | 工作量 |
|---|---|---|---|---|---|
| Q1 | FE 嵌入补丁对照 | 🟡 **设计就绪**（Q7 已解锁，待 PathD/GPU 上启动对照）| 与 EXP002a 部分重合；需在 PathD baseline 口径下重新组织 | P0 | 大（~5-10 天）|
| Q2 | 序列相似度划分 | 🟢 **PathD 进行中**（FASTA、四阈值聚类、样本平衡审计已完成） | 旧 PathC/EXP006 不继续；只作为历史参考 | P0 | 中 |
| Q3 | PubChem 3D 重对接 | 🔴 待启动 | 已有 `pubchem_3d_audit_2026-05-01.md` | TBD | 大 |
| Q4 | EGNN 换 GVP | 🟡 待启动 | EXP005 是双图（不是替换）→ -0.007；老师可能要"完全替换" | TBD | 中 |
| Q5 | 数据源版本号 + 长期网站 | 🔴 待启动 | — | TBD | 5A 小 / 5B 大 |
| Q6 | 正样本数据验证 | 🔴 待启动 | RCSB 682 条已通过多智能体配体分类 | TBD | 大 |
| Q7 | 原论文 Heme 处理排查 | ✅ **已完成**（2026-05-16，Claude + Codex 双方审计）| — | P0 | 极小（已耗时 ~2h） |
| Q8 | 糖基转移酶等家族扩展 | 🔴 待澄清 | — | TBD | TBD |
| Q9 | 底物分类下游预测（Stage 1）| 🟡 部分完成 | C3 v6 FINAL 8 类分类已完成 | TBD | 中 |
| Q10 | 打分系统 | 🔴 待澄清 | 可与 Q5B 合并 | TBD | TBD |
| Q11 | 6a15 模板重建结构 | 🔴 待澄清 | — | TBD | TBD |
| Q12 | 负样本杂泛性核验 | 🔴 待启动 | 可与 Q6 合并为"全数据集真实性核验" | TBD | 大 |
| Q13 | FE 催化原子距离筛选 | 🔴 待澄清 | "能量筛选"指代不明 | TBD | 中 |

**图例**：🔴 待启动 / 🟡 部分进行 / 🟢 进行中 / ✅ 完成 / ⏸️ 搁置

## PathD 服务器初始化状态（2026-05-21）

已完成：

- 新建服务器执行目录：`/root/autodl-tmp/EZSpecificity/PathD/P450`。
- 从 PathC 复制起步数据到 `data/base_from_PathC/`。
- 复制方式为解引用复制，PathD 内不保留指向 PathC 的符号链接。
- 复制最佳 baseline 快照到 `baselines/EXP001_allfix_unified_best/`。
- 修复 baseline 运行脚本：当前 `run_train.sh` 只指向 PathD 路径，并去掉原始 `--shutdown`；原始 PathC 脚本保留为 `run_train.original_from_PathC.sh`。

已复制数据：

| 数据块 | PathD 位置 | 大小 | 文件数 | 说明 |
|---|---|---:|---:|---|
| tables | `data/base_from_PathC/tables/` | 2.5M | 4 | `Enzymes.csv`、`Substrates.csv`、`data.csv`、`build_random_config.yml` |
| splits | `data/base_from_PathC/splits/` | 38M | 100 | PathC 原有划分文件 |
| features | `data/base_from_PathC/features/` | 13G | 17 | ESM、GROVER、Morgan、reaction features 等 |
| structure | `data/base_from_PathC/structure/` | 24G | 150545 | 结构文件、复合物、结构特征 LMDB 等 |
| cache_best_baseline | `data/base_from_PathC/cache_best_baseline/pt_cache_allfix_unified/random/` | 9.3G | 47422 | EXP001_allfix_unified 实际使用的 baseline cache，已解引用 |
| EXP001 baseline | `baselines/EXP001_allfix_unified_best/` | 130M | 60 | PathC 最佳 baseline 的代码、配置、脚本、结果、检查点 |

复查结论：

- `find PathD -type l` 结果为 0，PathD 内没有符号链接。
- 关键 CSV、baseline `ss.py`、`egnn.py`、`test_eval.json`、`metrics.csv` 哈希均与 PathC 源文件一致。
- `splits`、`features`、`structure`、`cache_best_baseline` 与 PathC 源目录文件数一致。
- 当前 `/root/autodl-tmp` 约 180G，总已用约 134G，剩余约 47G。

服务器端记录文件：

- `/root/autodl-tmp/EZSpecificity/PathD/P450/manifests/pathc_to_pathd_copy_manifest.md`
- `/root/autodl-tmp/EZSpecificity/PathD/P450/manifests/data_governance.md`
- `/root/autodl-tmp/EZSpecificity/PathD/P450/manifests/pathd_integrity_recheck_20260521_1751.md`
- `/root/autodl-tmp/EZSpecificity/PathD/P450/sessions/00_pathd_setup/session_log.md`

## Q7 关键发现（2026-05-16）

原模型 **完全看不见 HEM/HETATM 和 Fe** —— 蛋白原子词汇表无 Fe / 氨基酸词汇表无 HEM / PDB 解析器严格 ATOM 过滤丢弃 HETATM。原作者文档主观上想保留辅因子（example.ipynb 文档明确要求 cofactor 用 HETATM 标识），但下游 `PDBProtein` 解析器**静默丢弃** —— 这是**实现层 intent vs implementation gap**，是原论文的隐性缺陷。

→ 老师 Q1 的直觉对：加 Fe 嵌入是合理实验。改造点最小（3 行代码），与 EXP002a 在 P450 子集上做过的改造一致。详见 `sessions/07_Q7_原论文Heme处理排查/session_log.md`。

## Q2 当前进展（2026-05-21）

服务器 PathD 已重新启动 Q2：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q02_sequence_similarity_split/EXP001_mmseqs_cluster_split`

已完成：

- 从 `data/base_from_PathC/tables/Enzymes.csv` 导出 1622 条酶序列 FASTA。
- 完成序列审计：无空序列、无非法字符；1609 条唯一精确序列；11 组完全重复序列，涉及 24 个酶。
- 使用 `/root/miniconda3/bin/mmseqs` 跑完 id80 / id60 / id40 / id30 四个阈值聚类。
- 将 MMseqs2 TSV 解析为标准 CSV。
- 将聚类结果映射回 52254 条样本，完成正负样本平衡审计。

聚类摘要：

| 阈值 | clusters | 酶数 | 最大序列簇 | singletons |
|---|---:|---:|---:|---:|
| id80 | 1139 | 1622 | 19 | 900 |
| id60 | 758 | 1622 | 38 | 496 |
| id40 | 328 | 1622 | 148 | 162 |
| id30 | 140 | 1622 | 481 | 65 |

样本平衡摘要：

| 阈值 | 单簇最大样本数 | 单簇最大正样本数 | p95 样本/簇 |
|---|---:|---:|---:|
| id80 | 2085 | 348 | 134 |
| id60 | 2111 | 350 | 239 |
| id40 | 4027 | 377 | 596 |
| id30 | 12990 | 982 | 1199 |

当前判断：

- id30 过粗，最大簇太大，后续划分很难平衡。
- id80 过细，singleton 太多，可能更接近随机划分。
- id40 和 id60 是当前最值得继续比较的候选阈值。

下一步：选择主阈值或候选阈值组合，生成 cluster-held-out train/val/test split，并审计是否有簇泄漏、样本不平衡或正样本过少。

已补做独立 Codex 只读复核（SESSION_ID=`019e4a27-dc49-7e41-9524-0b16b20c67e1`）。复核未发现当前 Q2 导出、聚类、解析、样本回填和硬校验的实质错误；正式生成 split 前，需要把索引语义对齐、精确重复序列同簇、长度异常序列、多随机种子可行性模拟、跨划分最近邻相似度审计作为硬门槛。

2026-05-21 进一步核对 baseline cache 后，Q2 范围需要修正：当前聚类对象是 `Enzymes.csv` 全量 1622 条酶序列；PathC 最佳 baseline 实际使用的 `pt_cache_allfix_unified/random` 是过滤后的 cache 宇宙，fold0 train/val/test cache 合并后覆盖 1479 个 enzyme。后续正式 baseline-compatible 的 Q2 split 应先构建 cache-valid 样本和 enzyme 集合；1622 酶聚类保留为全量上游审计。

已按此修正整理数据层：

- `data/q02_sequence_similarity_split/exp001_full_catalog_1622/`：归档全量 1622 酶聚类审计。
- `data/actual_used_baseline/`：新增公共数据入口，物化 baseline cache 实际使用数据。
- `data/q02_sequence_similarity_split/exp002_actual_used_cache_valid/`：新增 Q2 EXP002 起点，基于实际使用数据生成 FASTA 和序列审计。

`actual_used_baseline` 当前统计：1479 个 enzyme、2111 个 substrate、44090 条样本，其中正样本 3913、负样本 40177。cache 与 PT CSV 的 enzyme/substrate 对齐 mismatch 均为 0。

EXP002 当前已完成 actual-used FASTA 和序列审计。MMseqs2 id80 原地重跑曾在 prefilter 阶段被系统终止；为避免继续卡在该步骤，当前候选 split 使用已成功完成的 EXP001 全量 1622 酶 MMseqs2 聚类结果，过滤到 actual-used 1479 酶集合，再生成 EXP002 cluster-held-out split。这个来源已经写入 manifest，不把它伪装成 actual-used MMseqs2 原地重跑结果。

EXP002 候选 split 已完成：

| 阈值 | clusters | 最大簇 enzyme | train 样本/正样本 | val 样本/正样本 | test 样本/正样本 | 泄漏审计 |
|---|---:|---:|---:|---:|---:|---|
| id80 | 1044 | 18 | 22031 / 1776 | 11021 / 1084 | 11038 / 1053 | 通过 |
| id60 | 691 | 37 | 21975 / 1861 | 11008 / 1046 | 11107 / 1006 | 通过 |
| id40 | 306 | 138 | 20254 / 1846 | 12515 / 1044 | 11321 / 1023 | 通过 |
| id30 | 127 | 455 | 20110 / 1887 | 12807 / 983 | 11173 / 1043 | 通过 |

四个阈值均覆盖 1479 enzyme / 44090 samples，cluster、enzyme、精确重复序列跨 train/val/test 泄漏均为 0。当前建议 id60 作为主方案，id40 作为更严格的对照。

EXP002 split 后已做独立复核，未发现阻塞性错误。根据复核建议，已补强并重跑验证脚本：sample key 改为多重集合核对；重复序列组中的 enzyme 必须显式存在；cluster enzyme set 必须等于 split sample enzyme set。四个阈值均通过补强验证。

`test_eval.json n_samples=11000` 与 `test/index.pt=10999` 的差异已追溯：baseline 日志中 DataModule 真实 test dataset 为 10999 samples；训练结束的 4 卡 DDP auto-test 会把 10999 补齐到 11000。PathD actual-used 数据集仍以 cache index 的 10999 test samples 为准。

整理后已做独立复核。Codex MCP 与 Codex 桥接均因超时或网络限制未能返回审查结果，随后使用独立子代理只读审查本地下载的服务器产物副本、生成脚本和 Q2 日志。复核未发现 actual-used 统计或 EXP002 当前状态的明显实质错误；指出 `organize_q02_exp002_actual_used_20260521.py` 是一次性混合脚本，不适合后续直接重跑。

## 四、推荐执行顺序（2026-05-21 更新）

不建议按 Q1 到 Q13 机械执行，应按依赖关系、实验成本和论文产出价值排序。

### 第 0 步：PathD 基线确认

PathD 服务器目录、基础数据和最佳 baseline 已经建立。正式训练前还需要做一次轻量 smoke test，确认 PathD 路径下能读取 cache、加载配置和跑最小 batch。

### 第 1 步：Q7 Heme/Fe 审计

Q7 已基本完成，是 Q1 和 Q13 的理论前置：原模型实现看不见 HEM/HETATM/Fe，因此 Fe 嵌入补丁和 Fe-催化原子距离筛选具有实验动机。

### 第 2 步：Q2 序列相似度划分

优先在 PathD 重新启动，不继续旧 PathC/EXP006。原因是 Q2 会影响后续所有模型对照的可信度，尤其是是否存在序列相似性导致的评估乐观。

无卡模式下可先做：FASTA 导出、重复序列审计、MMseqs2 脚本准备、划分策略设计。

### 第 3 步：Q1 Fe 嵌入补丁对照

建议等 Q2 划分方案明确后再启动正式训练。否则只在随机划分上证明 Fe 有用，结论容易被质疑。

对照关系应至少包含：

- PathD baseline：`baselines/EXP001_allfix_unified_best/`
- Fe embedding patch：基于同一数据、同一划分、同一训练设置，仅改变 Fe/HEM 相关表示。

### 第 4 步：Q3 + Q13

Q3 和 Q13 建议合并推进，因为都围绕结构与对接质量：

- Q3：用 PubChem 3D SDF 优化配体三维结构来源，并重新对接。
- Q13：把原先不清楚的“能量筛选”推进为 Fe 与催化原子距离筛选。

### 第 5 步：Q6 + Q12

Q6 和 Q12 本质上都是标签质量审计，建议合并为正负样本真实性核验：

- Q6：正样本是否真的有证据支持。
- Q12：负样本是否可能是假负样本，尤其考虑 P450 杂泛性。

### 第 6 步：Q4 EGNN 换 GVP

Q4 是模型结构大改，成本较高。建议放在 Q1/Q2 之后，避免在划分和 baseline 还没稳定时消耗大量时间。

### 第 7 步：Q9 + Q10

Q9 先做底物类别下游预测，Q10 再考虑包装成打分系统。Q10 不应先于 Q9，否则容易变成没有模型证据支撑的界面工程。

### 第 8 步：Q5 数据源版本号与长期网站

Q5 可以穿插做小版本记录，但长期网站适合等主要数据和验证规则稳定后再做。可作为论文未来工作或项目工程化方向。

### 第 9 步：Q8 + Q11

Q8 家族扩展和 Q11 6a15 模板重建结构容易扩大项目范围，建议后置。除非老师明确要求，否则不作为当前主线。

近期最推荐的实际顺序：

1. 选择 Q2 的主阈值或候选阈值组合，优先比较 id60 与 id40。
2. 生成 Q2 的 cluster-held-out train/val/test split，并做泄漏和平衡审计。
3. 准备 Q1 Fe 嵌入补丁的实验目录和代码差异。
4. 开 GPU 后，先跑 Q2 新划分下 baseline，再跑 Q1 Fe 嵌入对照。

## 五、与其他 Path 的关联

| Path 资源 | 对应问题 | 关系 |
|---|---|---|
| `/root/autodl-tmp/EZSpecificity/PathD/P450` | 全部 | 后续服务器主执行区 |
| `/root/autodl-tmp/EZSpecificity/PathD/P450/data/base_from_PathC/` | 全部 | PathD 起步数据，已从 PathC 解引用复制 |
| `/root/autodl-tmp/EZSpecificity/PathD/P450/baselines/EXP001_allfix_unified_best/` | Q1/Q2/Q4 | PathC 最佳 baseline 的可复现快照 |
| `PathC/.../C3_P450专属模型训练/EXP006_cluster_split_allfix_unified` | Q2 | 旧预实验记录；后续不继续，只作历史参考 |
| `PathC/.../C3_P450专属模型训练/.../C3 底物分类` | Q9 | 上游（分类已做完） |
| `PathC/.../pubchem_3d_audit_2026-05-01.md` | Q3 | 已有审计报告 |
| `PathC/.../C2 P450数据集构建/.../RCSB配体分类` | Q6 | 已对 682 条 RCSB 完成多智能体分类 |

## 六、待跟老师对齐的关键问题（先于动手）

1. **Q1 / Q4 矛盾**：EXP002a (Fe/HEM 在 P450 干净数据 -0.005) 和 EXP005 (双图 GVP -0.007) 老师是否知晓？是否影响他对 Q1/Q4 的判断？
   - 注意 PPT 第 15 页 "训练结果"页里 EXP005 反而比 EXP001 高 0.0008（vs 记忆里 -0.0067），数字差异需要核实
2. **Q1 训练集范围**：ESIBank 全集 retrain（大工程）or P450 子集（已做完）?
3. **Q13 "能量筛选" 指代**：对接 vina score？hard negative 阈值？还是别的？
4. **整体优先级**：13 条不可能全做到能写论文的深度。毕业论文核心贡献定哪 2-3 条？

详见各 session 内的 "待澄清点" 章节。

## 七、变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | PathD 创建，13 问题 session 骨架建立 |
| 2026-05-21 | 服务器端 PathD/P450 正式建立；复制 PathC 起步数据和 EXP001 最佳 baseline；完成解引用、哈希、符号链接和脚本路径复查；更新推荐执行顺序 |
| 2026-05-21 | Q2 在服务器 PathD 重新启动；完成 FASTA、序列审计、MMseqs2 四阈值聚类和样本平衡审计 |
| 2026-05-21 | 补做 Q2 独立 Codex 只读复核，SESSION_ID=`019e4a27-dc49-7e41-9524-0b16b20c67e1`；当前结果通过复核，但正式 split 前新增若干硬审计门槛 |
| 2026-05-21 | 修正 Q2 聚类对象口径：当前 1622 酶聚类是全量审计；baseline-compatible split 应改以 1479 个 cache-valid enzyme 为起点 |
| 2026-05-21 | 整理 Q2 数据层并生成 `actual_used_baseline`；创建 EXP002 actual-used 起点，完成实际使用数据 CSV、FASTA 和序列审计；MMseqs2 聚类待重跑 |
| 2026-05-21 | EXP002 整理后复核完成：actual-used 数据自洽性通过；记录 Codex 工具不可用、独立子代理复核结果和一次性混合脚本风险 |
| 2026-05-21 | EXP002 actual-used 候选 split 已生成并通过完整性验证；id60 建议为主方案，id40 建议为严格对照 |
| 2026-05-21 | EXP002 split 独立复核完成并补强验证；11000 vs 10999 差异已解释为 4 卡 DDP auto-test padding |
