# P450-EZSpecificity 研究项目

本仓库用于整理和推进 P450 酶-底物特异性预测研究。项目基于 EZSpecificity 模型，围绕 P450 家族的数据构建、模型诊断、实验复现和后续改进任务展开。

当前主线已经从 PathC 进入 PathD。PathC 作为历史证据区保留；PathD 用于后续 13 个导师反馈问题的整理、实验设计、日志记录和结果汇报。

## 当前状态

截至 2026-05-29 本地日志，重点进展如下：

| 模块 | 状态 | 说明 |
|---|---|---|
| PathC 历史基线 | 已完成 | 修复 ESM / GROVER 对齐问题后，随机划分基线达到约 0.93 Test AUC，作为后续对照来源 |
| PathD 目录整理 | 已完成 | 已建立独立 PathD 规划目录和服务器执行区，PathC 数据只作为只读来源 |
| Q7 原论文 Heme 处理排查 | 已完成 | 原模型蛋白侧看不见 HEM/Fe，可作为 Q1 和 Q13 的实验动机 |
| Q2 序列相似度划分 | 已完成一轮可汇报实验 | 已完成 PathD random 对照、id60 聚类簇划分、strict NN40/60/80 和 strict NN60 分布敏感性实验 |
| Q1 Fe/HEM 嵌入补丁 | 阶段性完成 | 主线保留 EXP001 / EXP003 对比；EXP003 是 full-data experiment，但未训满 35 epoch |
| Q9 底物分类下游预测 | 部分完成 | 2125 个化合物已分为 8 类；下游酶到类别预测模型尚未启动 |
| Q10 湿实验打分系统 | 第一版完整排序已完成 | Q10_EXP001 已对 diosgenin + 204 条候选 P450 完成结构、对接、PT cache、模型打分和排序；只能作为 preliminary ranking |
| Q3/Q4/Q5/Q6/Q8/Q11/Q12/Q13 | 待启动或待澄清 | 主要还停留在计划或已有材料阶段，后续按依赖关系推进 |

Q2 当前最重要的结果如下。这里优先看 Test AUC；Test AUPR 受测试集正样本数影响，汇报时要一起说明测试集规模。

| 实验 | 划分方式 | Test AUC | Test AUPR | Test 样本 / 正样本 | 说明 |
|---|---|---:|---:|---:|---|
| EXP008 | PathD random split 对照 | 0.934206 | 0.686618 | 10999 / 984 | PathD 管线能复现随机划分高分 |
| EXP004 | id60 聚类簇划分 | 0.807786 | 0.321236 | 11107 / 1006 | 同一 id60 聚类簇不跨 train / val / test |
| EXP007 | strict NN80 | 0.819789 | 0.383128 | 11123 / 1160 | 80% 阈值较宽松，保留更多近缘信息 |
| EXP005 | strict NN60 主实验 | 0.670801 | 0.218569 | 11091 / 1006 | 回答老师提出的“test 与 train 序列相似度 <60%” |
| EXP006 | strict NN40 | 0.638403 | 0.102411 | 11755 / 899 | 更严格的远缘泛化设置 |
| EXP010 | strict NN60，约 7:1.5:1.5 | 0.733258 | 0.225581 | 6621 / 632 | strict NN60 下提高训练集比例后的补充结果 |
| EXP011 | strict NN60，约 8:1:1 | 0.602512 | 0.169727 | 4466 / 435 | 测试集更小且更集中，作为分布敏感性证据 |

这组结果说明：随机划分下的高分不能直接代表远缘新酶泛化能力。把 test/val 到 train 的近缘酶限制住之后，模型评估难度明显上升。EXP010 和 EXP011 说明同为 strict NN60 时，train/val/test 比例和测试集组成也会影响最终指标，因此它们作为敏感性证据，不替代 EXP005 主结果。

Q1 当前进度：

| 实验 | Full test AUC | Full test AUPR | 389 P450 子 test AUC | 389 P450 子 test AUPR | 当前含义 |
|---|---:|---:|---:|---:|---|
| EXP001 原始 A-only baseline | 0.894870 | 0.569394 | 0.890075 | 0.457316 | 原模型对照 |
| EXP003 P450 专项 HEM/Fe overlay | 0.863196 | 0.493476 | 0.908318 | 0.617855 | 更贴近老师问题的 P450 专项补丁 |

当前解释：

1. EXP003 是 full-data experiment，不能写成只训练 P450 子集；P450 子集只用于 subgroup analysis 和 Fe/HEM target accounting。
2. EXP003 的 P450 子 test 表现高于 EXP001，说明补充 HEM/Fe 对 P450 子集有阶段性收益。
3. EXP003 的 full test 表现低于 EXP001，说明当前补丁会牺牲全量 ESIBank 泛化表现。
4. EXP003 未训满 35 epoch，当前结果适合组会阶段性汇报，不应写成完整训练闭环。
5. EXP002 不放入主线。它是早期宽口径尝试，`3185` 条 overlay 中只有 `1140` 条属于 389 P450 清单，另外 `2045` 条来自非 P450 清单样本。

Q10 当前进度：

| 内容 | 当前结论 |
|---|---|
| 实验 | Q10_EXP001_diosgenin_20260528 |
| 输入 | diosgenin + 204 条候选 P450，其中 MDP450 43 条、ARATH 161 条 |
| 完成环节 | ColabFold 结构、HEM/Fe receptor、Uni-Dock 对接、PT cache、模型打分和 candidate_id 一致性审计 |
| 样本数 | PT cache `204` 条，模型打分 `204` 条，全流程审计 `ok=true` |
| 模型 | Q2 EXP008 random best checkpoint，`ep79-auc0.9316.ckpt` |
| 分数范围 | `score_min = 1.93e-16`，`score_max = 0.088214606`，`score_mean = 0.002178010` |
| 本地表格 | `results/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/model_scores/Q10_EXP001_两个酶列表完整打分表.md` |

Q10_EXP001 可以作为第一版候选排序和后续严格版对照，不能解释成湿实验概率。已知风险包括：ColabFold 曾依赖在线 MSA，HEM/Fe 来自少数模板 overlay，Uni-Dock 未按官方高吞吐入口重做，且 docking score 存在异常正分。下一版建议新建 Q10_EXP002，按官方文档和本机 `--help` 重新设计结构、对接和小批性能基准。

## 推荐阅读顺序

新打开项目时，建议按下面顺序看：

1. [AGENTS.md](AGENTS.md)
   项目协作规则、证据优先级、服务器安全边界和子智能体复核要求。

2. [docs/PROJECT_CONTEXT_EZSpecificity.md](docs/PROJECT_CONTEXT_EZSpecificity.md)
   长期项目背景归档，适合了解 PathA / PathB / PathC / PathD 的整体关系。

3. [PathD_计划与进度.md](毕业设计/P450_EZSpecificity_研究项目/PathD_2026-05-08_导师反馈13问题/PathD_计划与进度.md)
   PathD 的 13 个导师反馈问题、目录规则和当前推荐执行顺序。

4. [Q2序列相似度划分实验汇总.md](毕业设计/P450_EZSpecificity_研究项目/PathD_2026-05-08_导师反馈13问题/sessions/02_Q2_序列相似度划分/Q2序列相似度划分实验汇总.md)
   Q2 的整理版说明，适合快速了解实验逻辑和 PPT 可用结论。

5. [Q2 session_log.md](毕业设计/P450_EZSpecificity_研究项目/PathD_2026-05-08_导师反馈13问题/sessions/02_Q2_序列相似度划分/session_log.md)
   Q2 的详细过程记录，包含服务器路径、审计结果、训练参数、监控和复核说明。

6. [Q1 session_log.md](毕业设计/P450_EZSpecificity_研究项目/PathD_2026-05-08_导师反馈13问题/sessions/01_Q1_FE嵌入补丁对照/session_log.md)
   Q1 原模型基线、EXP002 宽口径尝试、EXP003 P450 专项 HEM/Fe overlay 和结果边界。

7. [Q10 session_log.md](毕业设计/P450_EZSpecificity_研究项目/PathD_2026-05-08_导师反馈13问题/sessions/10_Q10_打分系统/session_log.md)
   Q10 湿实验打分系统的第一版流程、事故复盘、完整排序表和下一版硬规则。

## 目录入口

```text
EZSpecificity_Project/
├── src/                         原始 EZSpecificity 代码
├── data/                        原始示例数据，大文件默认不入库
├── saved_model/                 原始模型权重，大文件默认不入库
├── docs/                        长期项目背景与归档文档
├── AGENTS.md                    当前项目协作规则
├── README.md                    当前入口说明
└── 毕业设计/
    ├── 毕业论文/                论文正文、图表和写作材料
    └── P450_EZSpecificity_研究项目/
        ├── PathA_*/             模型评估测试集构建
        ├── PathB_*/             P450 数据集构建与结构优化
        ├── PathC_*/             P450 专属模型训练历史阶段
        └── PathD_*/             导师反馈 13 问题与后续实验主线
```

PathD 的服务器执行区为：

```text
/root/autodl-tmp/EZSpecificity/PathD/P450
```

服务器端文件默认按只读资产处理。需要新结果时优先新建目录、日志或带说明的版本文件，避免覆盖旧产物。

## PathD 当前任务地图

| 问题 | 主题 | 当前建议 |
|---|---|---|
| Q1 | Fe/HEM 嵌入补丁对照 | EXP001 / EXP003 已有阶段性对比；EXP003 的 P450 子 test 提升、full test 下降，且未训满 35 epoch |
| Q2 | 序列相似度划分 | 已有可汇报结果；EXP005 是 strict NN60 主实验，EXP010/EXP011 是分布敏感性补充 |
| Q3 | PubChem 3D 重对接 | 后续单独规划，需要处理结构与对接数据 |
| Q4 | EGNN 换 GVP | 暂不优先，历史 GVP 旁路结果没有稳定增益 |
| Q5 | 数据源版本号与长期网站 | 适合作为工程化整理任务 |
| Q6 / Q12 | 正负样本真实性核验 | 建议合并规划为标签质量审计 |
| Q7 | 原论文 Heme 处理排查 | 已作为 Q1 / Q13 的理论前置 |
| Q8 | 糖基转移酶等家族扩展 | 工作量大，后置 |
| Q9 | 底物分类与下游预测 | 8 类底物分类已完成；下游酶到类别预测模型尚未启动 |
| Q10 | 湿实验完整打分系统 | Q10_EXP001 已完成 preliminary ranking；Q10_EXP002 应按严格结构和对接流程重做 |
| Q11 | 6a15 模板重建结构 | 后置，需老师明确优先级 |
| Q13 | Fe 与催化原子距离筛选 | 可与 Q3 的结构更新合并考虑 |

## 当前汇报重点

当前最适合向老师汇报的是 Q2，Q1 可以作为第二条阶段性结果，Q10 作为应用落地的 preliminary ranking：

1. 原随机划分基线表现高，但可能受到近缘酶影响。
2. Q2 先还原 baseline 真正使用的 1479 个 enzyme 和 44090 条样本。
3. EXP008 在 PathD 内复现 random split 高分，Test AUC 为 0.934206。
4. EXP005 做 strict NN60 划分，Test AUC 为 0.670801。
5. EXP006/EXP007 给出 strict NN40/80 阈值梯度，EXP010/EXP011 给出 strict NN60 分布敏感性。
6. 性能下降说明严格新酶泛化比随机划分更难，也更接近老师提出的问题。
7. Q1-EXP003 在 389 P450 子 test 上高于 EXP001，但 full test 下降，且 EXP003 未训满 35 epoch。
8. Q10-EXP001 已跑通 diosgenin + 204 条候选 P450 的第一版完整排序，但必须带 preliminary 和质量风险标签使用。

详细数字、划分机制和图表建议见 Q2 汇总文档。

## 协作与安全规则

- 回答实验事实时先查日志、配置、结果文件和检查点。
- 不把旧日志中的“当前阶段”直接当作今天的状态。
- 不删除、不覆盖服务器原有文件；需要改动远程资产前先确认。
- 重要实验、目录调整、结果解释和下一步安排应写入对应 `session_log`。
- 临时目录、本地缓存、训练产物和大文件默认不入库。
