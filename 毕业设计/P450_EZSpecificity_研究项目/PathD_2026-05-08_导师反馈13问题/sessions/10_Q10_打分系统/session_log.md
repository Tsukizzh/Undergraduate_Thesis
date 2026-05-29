# Q10 湿实验完整打分系统

## 老师原话与最新定义

> 做成打分系统

2026-05-28 最新定义：Q10 不再只理解为“输出一个模型分数”或“做一个网页界面”，而定义为面向真实湿实验的完整打分系统。

目标输入是真实湿实验给出的底物 SMILES 和候选 P450 蛋白序列列表；目标输出是候选 P450 的排序、分数、结构与对接依据、湿实验优先级建议，以及后续实验结果回填表。

第二版按完整流程执行，不把“只用序列和 SMILES 的快速打分”作为目标结果。结构准备、底物 3D、重新对接、模型输入缓存、打分排序、湿实验回填都属于 Q10 第二版的必要组成部分。

## 状态

🟡 **已重新定义，可启动设计**

## 当前真实输入场景

| 输入 | 当前材料 | 用途 |
|---|---|---|
| 候选酶 | 去除无跨膜域后的 MDP450 序列，约 128 条 | 构建待排序 P450 列表 |
| 补充酶库 | 拟南芥 P450 FASTA | 后续可扩展为第二批候选 |
| 底物 | diosgenin SMILES | 构建 1 个底物对 128 个 P450 的候选配对 |
| 任务 | 用最新模型给候选 P450 排名 | 指导湿实验优先验证 |

## 第二版系统的强制范围

第二版先做成可复现的完整湿实验排序流程，暂时不优先做网页界面。这个流程必须稳定接收真实湿实验输入，并产出排序表、结构证据表、对接证据表和湿实验回填表。网页化可以留给 Q5，完整结构流程不能省略。

建议拆成 6 层：

| 层级 | 做什么 | 输出 |
|---|---|---|
| 0. 输入审计 | 读取 P450 序列和底物 SMILES，检查 ID、重复序列、长度、非法字符、SMILES 是否能解析 | `01_input_audit.csv` |
| 1. P450 结构重建 | 为 128 条 P450 获取或预测结构，记录结构来源、质量、是否需要 HEM/Fe 补全 | `02_structure_status.csv` |
| 2. 底物 3D 准备 | 为 diosgenin 获取 PubChem 3D SDF，缺失时用 RDKit 生成并记录原因 | `03_substrate_3d_status.csv` |
| 3. 重新对接 | 对 128 条 P450 与 diosgenin 重新对接，保留 pose、能量、口袋和失败原因 | `04_docking_results.csv` |
| 4. 模型缓存与打分 | 把新结构和新对接结果转换成模型可读缓存，再用最新模型输出分数 | `05_model_scores.csv` |
| 5. 湿实验优先级 | 融合模型分数、对接分数、Fe/催化原子距离、结构质量等信息，给出 Top 候选 | `06_final_priority_list.csv`、`07_wetlab_feedback_template.csv` |

## 和现有模型的关系

当前 Q2 random split 最优模型可以作为排序核心。它的输出本质是酶-底物配对的模型分数，经过 sigmoid 后可以得到 0 到 1 的分数，适合用于排序。

需要注意：这个分数不能直接当成湿实验成功率。更稳妥的用法是给 128 条 P450 排序，并把 Top 10、Top 20、Top 30 作为不同优先级给湿实验验证。

## 还缺哪些东西

| 缺口 | 为什么需要 |
|---|---|
| 原始文件落盘 | 需要把 128 条 P450 序列、拟南芥 FASTA、diosgenin SMILES 放入 Q10 输入目录 |
| 序列 ID 规范 | 后续结果表要能对应回湿实验编号 |
| 结构来源 | 完整版必须为每条 P450 准备结构；没有现成结构时需要预测结构并记录来源 |
| HEM/Fe 处理 | P450 的结构证据需要检查 HEM/Fe 是否存在；缺失时要记录是否补全以及补全方法 |
| diosgenin 3D 构象 | 重新对接必须使用 3D 底物构象，优先 PubChem 3D SDF，缺失时再用 RDKit |
| 对接参数 | 需要确定对接盒子、口袋中心、失败重跑规则和结果保留格式 |
| 排序阈值 | 要和老师确认湿实验先做 Top 10、Top 20 还是全部 128 条 |

## 建议的第一个实验

实验名：`Q10-EXP001_diosgenin_128P450_wetlab_scoring`

目标：把这次真实湿实验输入整理成标准输入，并完成结构重建、底物 3D 准备、重新对接、模型缓存重建、模型打分和最终湿实验优先级排序。

必须输出：

1. `input_manifest.csv`
2. `enzyme_sequence_audit.csv`
3. `substrate_smiles_audit.csv`
4. `structure_status.csv`
5. `substrate_3d_status.csv`
6. `docking_results.csv`
7. `model_cache_manifest.csv`
8. `model_scores.csv`
9. `final_priority_list.csv`
10. `wetlab_feedback_template.csv`

## 与其他问题的关联

- Q3：Q10 第二版会直接用到 PubChem 3D 或 RDKit 3D，并要求重新对接。
- Q9：底物类别预测可以作为粗筛先验，帮助解释 diosgenin 属于哪类底物。
- Q13：结构增强版可以加入 Fe 与催化原子距离，作为最终排序的物理约束。
- Q5：等流程稳定后，再考虑做成长更新网站或在线打分入口。

## 变更日志

| 日期 | 变更 |
|---|---|
| 2026-05-08 | session 创建 |
| 2026-05-28 | Q10 重新定义为“湿实验完整打分系统”，补充真实输入场景和第二版流程 |
| 2026-05-28 | 明确 Q10 第二版必须走完整结构流程：结构重建、底物 3D、重新对接、模型缓存和最终排序都不能省略 |

## 2026-05-28 晚间服务器实际推进记录

### 1. 当前输入数量已经重新核对

这次不再按截图里的“128 条”估算。按服务器上当前真实上传文件解析结果，Q10-EXP001 的候选酶一共有 204 条：

| 来源列表 | 当前解析条数 | 说明 |
|---|---:|---|
| MDP450_no_TM_removed | 43 | 来自 `MDP450-400-600.no-TM-removed.txt`，原始文件保留不改，处理版只去掉 FASTA 终止符 `*` |
| ARATH_uniprot_p450 | 161 | 来自拟南芥 P450 FASTA |
| 合计 | 204 | 都和 diosgenin 组成候选酶-底物配对 |

服务器输入目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528`

### 2. ColabFold 已经部署并通过试跑

已按“优先 ColabFold”的方案处理。服务器直接部署和权重下载容易受网络影响，所以采用了本地克隆/上传补齐：

| 项目 | 当前状态 |
|---|---|
| ColabFold 源码 | 已从本地克隆上传到服务器 `deps/ColabFold` |
| 运行环境 | `/root/autodl-tmp/envs/colabfold_q10` |
| ColabFold 版本 | 1.6.1 |
| JAX/GPU | JAX 0.6.2，可见 GPU0 和 GPU1 |
| AlphaFold 参数 | 已放到 `/root/.cache/colabfold/params`，存在 `download_finished.txt` |

先做了单条 MDP450 的 MSA 模式试跑，输出路径：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/structures/colabfold_smoke/results_msa_v1`

试跑结果：

| 指标 | 数值 |
|---|---:|
| 平均 pLDDT | 86.774 |
| pTM | 0.88 |
| 输出 PDB | 已生成 unrelaxed PDB |

这说明 ColabFold 主程序、在线 MSA、参数缓存和 GPU 推理链路都能跑通。

### 3. 正式批量结构预测已经启动

已生成 ColabFold 批量输入：

| 文件 | 用途 |
|---|---|
| `q10_candidates_204.fasta` | 204 条候选酶总表 |
| `MDP450_no_TM_removed.fasta` | 43 条 MDP450 |
| `ARATH_uniprot_p450.fasta` | 161 条拟南芥 P450 |
| `q10_colabfold_input_manifest.csv` | 每条序列的 ID、来源、长度和描述 |

正式批量预测分两组跑：

| 分组 | GPU | 输出目录 | 参数 |
|---|---:|---|---|
| ARATH_uniprot_p450，161 条 | GPU0 | `structures/colabfold_batch/CF_ARATH_msa_m1r1_v1` | MSA 模式，1 model，1 recycle |
| MDP450_no_TM_removed，43 条 | GPU1 | `structures/colabfold_batch/CF_MDP450_msa_m1r1_v1` | MSA 模式，1 model，1 recycle |

当前这版是全量结构的第一版，目标是先让 204 条候选都进入完整 Q10 流程。后续如果需要更高质量结构，可以对 Top 候选、低 pLDDT 条目或失败条目单独用更多 model/recycle 重跑，并记录为独立 rerun，不覆盖当前结果。

### 4. 仍未完成的部分

| 部分 | 当前状态 |
|---|---|
| 204 条结构预测 | 正在服务器后台运行 |
| 结构 manifest | 等批量结果完成后生成 |
| HEM/Fe 与对接盒子规则 | 需要结合预测结构和 P450 保守位点继续确定 |
| Uni-Dock 重新对接 | 等结构输出后执行 |
| 模型可读 PT 样本 | 等对接结构后构建 |
| Q2 random 最优模型打分 | 等 PT 样本完成后执行 |
| 最终湿实验排序表 | 等模型分数和对接结果合并后生成 |

### 5. 2026-05-28 22:05 当前进度复查

服务器当前仍在跑 ColabFold 批量预测，两个进程都存在，没有看到 fatal error。

| 分组 | 总数 | 已完成结构 | 仍待完成 |
|---|---:|---:|---:|
| MDP450_no_TM_removed | 43 | 16 | 27 |
| ARATH_uniprot_p450 | 161 | 18 | 143 |
| 合计 | 204 | 34 | 170 |

当前已完成结构的平均 pLDDT 为 88.427。已同步更新：

- `structures/manifests/q10_colabfold_structure_manifest_msa_m1r1_v1.csv`
- `structures/manifests/q10_colabfold_structure_summary_msa_m1r1_v1.json`
- `docking/manifests/q10_unidock_box_manifest_msa_m1r1_v1.csv`
- `docking/manifests/q10_unidock_box_summary_msa_m1r1_v1.json`

对接盒子目前按三类来源记录：

| box 来源 | 当前数量 | 含义 |
|---|---:|---|
| P450 heme-binding motif | 24 | 找到 P450 常见保守 Cys 附近位点，用该位点附近坐标作为对接中心 |
| C-terminal CxG fallback | 2 | 没找到完整 motif，但找到末端附近 CxG 线索 |
| protein centroid fallback | 8 | 当前没有可靠 motif，只能先用蛋白质中心，后续需要单独标为低置信对接 |

ColabFold 输出是 protein-only 结构，不自带 HEM/Fe。当前第一版 Q10 用 motif/几何中心作为对接盒子来源，manifest 里已经记录每个候选的 box 来源。后续如果要做更严格版本，需要再做 HEM/Fe overlay 或 AlphaFill 之类的补全，并把它作为独立实验版本记录。

Uni-Dock 链路已经做过单条 smoke test。重要结论：

1. 需要把 `/root/autodl-tmp/envs/ambertools_q10/bin` 和 `/root/autodl-tmp/envs/unidock_env/bin` 同时加入 PATH。
2. ColabFold PDB 没有氢，`unidock_pipeline` 需要加 `-ph`，触发 reduce/tleap 补氢流程。
3. 单条 smoke 输出了 docking score，说明 ligandprep、proteinprep、补氢和 Uni-Dock 主程序能跑通。

已上传批量 Uni-Dock 脚本：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/scripts/q10_run_unidock_batch_20260528.py`

脚本已通过服务器端 `py_compile` 检查。等 ColabFold 结构完成后，可以直接用这个脚本批量跑 204 条候选结构的 diosgenin 对接。

### 6. 2026-05-28 22:38 链路冒烟测试进展

这一轮没有把 Q10 当作已完成，只是先把最小链路打通，避免等 204 条结构全部完成后才发现脚本或数据格式不通。

#### 6.1 ColabFold 当前状态

服务器时间 `2026-05-28 22:38:55 CST` 复查时，两个 ColabFold 批量进程仍在运行。

| 分组 | 总数 | 当前已输出 PDB | 当前已输出 scores JSON |
|---|---:|---:|---:|
| MDP450_no_TM_removed | 43 | 24 | 24 |
| ARATH_uniprot_p450 | 161 | 27 | 27 |
| 合计 | 204 | 51 | 51 |

GPU 在复查瞬间显示 0% 利用率，但进程仍存在；ColabFold 在 MSA、写文件、模型间隙会出现这种阶段性低占用。后续是否完成以结构 manifest 的 complete 数为准。

#### 6.2 Uni-Dock 批处理 smoke 已通过

已用当前完成结构中的 2 条 MDP450 样本跑通 Uni-Dock 批处理脚本。输出结果位于：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/docking/manifests/q10_unidock_results_msa_m1r1_v1.csv`

| candidate_id | docking score | box 来源 | 状态 |
|---|---:|---|---|
| MDP450_no_TM_removed_0028 | -4.494 | P450 heme-binding motif | complete |
| MDP450_no_TM_removed_0029 | -6.722 | P450 heme-binding motif | complete |

这只能说明批处理对接链路可运行，不能代表完整 Q10 排序结果。

#### 6.3 Q10 模型 PT 缓存构建 smoke 已通过

已把上面 2 条对接结果转换成 Q2 random 模型可读取的 PT 缓存：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/model_cache/q10_unidock_pt_msa_m1r1_v1`

构建过程中发现一个重要格式问题：Uni-Dock 输出的 SDF 带显式氢，第一次读取时 diosgenin 被读成 72 个原子，和 Q2 训练缓存中的底物原子规模不一致。已在 Q10 PT 构建脚本中改为读取后去掉显式氢，修正后 diosgenin 为 30 个重原子，更接近训练缓存里的底物表示。

当前 2 条 smoke 样本统计：

| candidate_id | ligand atoms | pocket protein atoms | bonds | KNN edges | pocket 方法 |
|---|---:|---:|---:|---:|---|
| MDP450_no_TM_removed_0028 | 30 | 452 | 70 | 23066 | radius_10A |
| MDP450_no_TM_removed_0029 | 30 | 289 | 70 | 15242 | radius_10A |

验证脚本已确认 EXP008 的 `PtCacheDataset` 能读取这个 Q10 缓存，batch 也能正常拼接。

#### 6.4 Q2 random 最优模型 smoke 推理已通过

已使用 Q2 random split 最优 checkpoint 做 2 条样本的推理：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q02_sequence_similarity_split/EXP008_random_gdtable_b80_retry_20260525/results/checkpoints/pt-Q2_EXP008_random_gdtable_b80_retry_full_20260525_174150-ep79-auc0.9316.ckpt`

输出目录：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/model_scores`

当前 2 条 smoke 分数：

| rank | candidate_id | model_score_sigmoid | model_logit | docking score |
|---:|---|---:|---:|---:|
| 1 | MDP450_no_TM_removed_0029 | 0.009170 | -4.682597 | -6.722 |
| 2 | MDP450_no_TM_removed_0028 | 0.000220 | -8.422556 | -4.494 |

这些分数只能当作模型排序分数，不能解释为湿实验成功概率。当前只有 2 条 smoke 样本，不写入正式结论。

本轮也补了 checkpoint 权重键审计，`load_state_dict_missing_count = 0`，`load_state_dict_unexpected_count = 0`，说明当前推理脚本加载的 EXP008 模型定义和 checkpoint 权重键完整对上。

#### 6.5 当前还不能写成完成的内容

| 内容 | 当前状态 |
|---|---|
| 204 条结构预测 | 还在跑，当前已见 51 条 PDB |
| 204 条 Uni-Dock 对接 | 还没开始全量，只完成 2 条 smoke |
| Q10 全量 PT 缓存 | 还没构建，只构建 2 条 smoke |
| Q10 全量模型打分 | 还没完成，只验证 2 条 smoke |
| HEM/Fe 参与结构或对接 | 当前 ColabFold 是 protein-only，第一版盒子来自 motif/centroid 规则，不等于真实 HEM/Fe 中心 |

下一步等待 ColabFold 结构完成后，重新生成结构 manifest 和盒子 manifest，再跑完整 Uni-Dock、重建完整 PT 缓存，最后用 Q2 random 最优 checkpoint 对两组酶分别输出正式排序表。

### 7. 2026-05-28 22:48 全量执行脚本准备

为了减少 ColabFold 完成后手动输入长命令造成的错误，已新增并上传一个全量流程脚本：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/scripts/q10_run_full_after_colabfold_20260528.sh`

本地对应文件：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\scripts\q10_run_full_after_colabfold_20260528.sh`

脚本已在服务器上通过 `bash -n` 语法检查，并已加执行权限。

这个脚本的执行顺序：

1. 重新收集 ColabFold 结构 manifest，并链接 selected PDB。
2. 检查 `total=204` 且 `complete=204`；如果结构没齐，脚本直接停止。
3. 生成 Uni-Dock 对接盒子 manifest。
4. 用 GPU0 和 GPU1 跑全量 Uni-Dock。
5. 检查对接结果是否 204 条全部 complete。
6. 构建完整 Q10 PT 缓存。
7. 检查 `samples_written=204` 且 `failures=0`。
8. 用 EXP008 的数据集封装验证缓存可读。
9. 用 Q2 random best checkpoint 推理打分。
10. 检查输出样本数为 204，且 checkpoint 权重键完全对齐。

当前没有运行这个全量脚本，因为 ColabFold 结构还没全部完成。它只是作为后续的安全执行入口。

### 8. 2026-05-28 22:58 启动等待器

由于 ColabFold 当前慢在 MSA/序列搜索阶段，GPU 推理本身只需数秒，暂时没有中断当前两个 ColabFold 进程。为了避免结构完成后无人接续，已启动一个后台等待器：

`/root/autodl-tmp/EZSpecificity/PathD/P450/experiments/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/scripts/q10_wait_and_run_full_20260528.sh`

后台 PID：

`112587`

等待器日志：

`/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/logs/q10_full_pipeline/q10_wait_and_run_full_20260528.nohup.log`

等待器行为：

1. 每隔约 5 分钟重新收集 ColabFold manifest。
2. 如果 `total=204` 且 `complete=204`，自动执行 `q10_run_full_after_colabfold_20260528.sh`。
3. 如果结构没齐，继续等待，不启动 Uni-Dock。

启动时第一次刷新结果：

| 分组 | total | complete | pending |
|---|---:|---:|---:|
| MDP450_no_TM_removed | 43 | 29 | 14 |
| ARATH_uniprot_p450 | 161 | 31 | 130 |
| 合计 | 204 | 60 | 144 |

当前 60 条完成结构的平均 pLDDT 为 `88.723`。

### 9. 2026-05-28 23:02 修正 Q10 对接前提

用户指出一个关键问题：ColabFold 预测出来的是 protein-only 结构，没有 HEM/Fe。这个判断是正确的。当前用 motif 或蛋白几何中心做对接盒子，只能算 smoke test，不能作为正式 Q10 对接证据。

已立即处理：

1. 停止后台等待器 PID `112587`，避免它在 204 条结构完成后自动进入 protein-only Uni-Dock。
2. 确认没有误伤 ColabFold 结构预测进程，两个原始 ColabFold 进程仍在运行。
3. 更新 `q10_run_full_after_colabfold_20260528.sh`，让它在收集完 ColabFold manifest 后直接停止，并提示“正式对接需要 HEM/Fe-aware receptor preparation”。
4. 更新 `q10_wait_and_run_full_20260528.sh` 的提示，避免误以为它会直接跑正式对接。

当前判断：

| 项目 | 状态 |
|---|---|
| ColabFold 结构 | 可继续作为 protein backbone 来源 |
| HEM/Fe | 还没有进入预测结构，正式对接前必须补或采用明确替代策略 |
| protein-only Uni-Dock | 只保留为 smoke，不作为正式 Q10 结果 |
| Q10 正式对接 | 暂停，等 HEM/Fe receptor 准备方案完成后再跑 |

### 10. 2026-05-28 23:07 ColabFold 加速批次

慢速原因已由日志确认：主要慢在在线 MMseqs2 MSA 搜索，每条序列常见需要数分钟；真正 AlphaFold GPU 推理常常只有数秒。因此单纯有 2 张 5090、50 核 CPU、180G 内存，也不能让在线 MSA 阶段线性加速。

为减少等待，已采用不破坏原结果的并行加速策略：

1. 不停止原来的两个 ColabFold 进程。
2. 重新收集当前 pending 序列。
3. 把 pending 序列拆成 6 个新批次。
4. 新批次输出到 `CF_ACCEL_*_msa_m1r1_v1` 目录。
5. 修改 manifest 收集脚本，使它可以同时识别原始批次和 `CF_ACCEL_*` 批次，哪个目录先出完整 PDB/score，就采用哪个结果。

新增/更新脚本：

| 脚本 | 用途 |
|---|---|
| `q10_collect_colabfold_manifest_20260528.py` | 已支持扫描 `CF_ACCEL_*` 加速批次 |
| `q10_launch_colabfold_accel_pending_20260528.sh` | 拆分 pending FASTA 并启动 6 个加速批次 |

加速批次已启动：

| shard | 条数 | GPU |
|---|---:|---:|
| CF_ACCEL_00 | 24 | 0 |
| CF_ACCEL_01 | 23 | 1 |
| CF_ACCEL_02 | 23 | 0 |
| CF_ACCEL_03 | 23 | 1 |
| CF_ACCEL_04 | 23 | 0 |
| CF_ACCEL_05 | 23 | 1 |

启动时 ColabFold manifest 状态：

| 分组 | total | complete | pending |
|---|---:|---:|---:|
| MDP450_no_TM_removed | 43 | 31 | 12 |
| ARATH_uniprot_p450 | 161 | 34 | 127 |
| 合计 | 204 | 65 | 139 |

### 11. ligand 和 protein 准备现状

底物 diosgenin 当前已经按用户要求使用 PubChem 3D SDF，而不是只用 SMILES 或 RDKit 随机构象：

| 文件 | 状态 |
|---|---|
| `ligands/pubchem_3d/diosgenin_CID99474_pubchem_3d.sdf` | 已存在，来源 PubChem3D，CID 99474 |
| `docking/ligands_pdbqt/diosgenin_CID99474_pubchem_3d.pdbqt` | 已存在 |
| `ligands/rdkit_fallback/diosgenin_rdkit_etkdg.sdf` | 备用，不是当前首选 |

已确认 `substrate_3d_status_v2.json` 记录当前 selected SDF 来源是 `PubChem3D`。

蛋白准备方面，之前只做过 2 条 protein-only smoke test。Uni-Dock 的 `-ph` 补氢/proteinprep 链路能跑通，但正式 Q10 还缺 HEM/Fe-aware receptor 准备，所以不能直接进入全量 docking。

### 12. 2026-05-29 Q10 流程事故复盘与硬性教训

用户指出 Q10 执行中存在严重流程问题：ColabFold 预测慢且没有充分利用服务器资源；P450 结构预测没有把 HEM/Fe 缺失作为前置风险；HEM/Fe 补充没有充分查询 AlphaFill 数据库；Uni-Dock 没有先查官方高吞吐模式。经重新查日志、脚本、服务器 `unidock --help` 和官方资料，这些批评大体成立，Q10_EXP001 不能当成理想正式流程。

#### 已确认的问题

| 问题 | 已查证据 | 影响 |
|---|---|---|
| ColabFold 依赖在线 MSA，出现限流 | `logs/colabfold_accel/*.log` 多次出现 `RATELIMIT` / `RUNNING` 等等待；部分 AlphaFold 推理本身只需数秒到百秒 | 慢点主要在 MSA 和网络等待，不是 5090 算力不足；正式批量前应先评估本地 MSA 或替代结构来源 |
| 存在残留 ColabFold 进程 | 2026-05-29 复查时仍看到旧 `colabfold_batch` ARATH 进程，虽然 manifest 已显示 204 条完成 | 说明流程收尾管理不合格，继续占用资源 |
| HEM/Fe 风险发现太晚 | Q10 初始 smoke 曾按 protein-only receptor 尝试；用户提醒后才停止等待器并改为 HEM/Fe-aware receptor | P450 对接和 Fe 距离相关任务必须把 HEM/Fe 作为设计前提 |
| AlphaFill 使用不充分 | Q10 receptor 脚本使用固定 `6VBY,5YLW,7X2Q,6A15,8E83,3RUK,2VE3`，每个最多 2 个模板；没有逐条候选查询完整 AlphaFill entry | 少数模板 overlay 只能作为临时应急，不能作为正式高置信度补 HEM/Fe 方案 |
| Uni-Dock 没有先用官方一对一批处理模式 | 服务器 `unidock --help` 显示支持 `--paired_batch_size`；当前脚本是每个 receptor 单独启动一次 `unidock --gpu_batch` | 浪费 Uni-Dock 高吞吐优势，启动开销和 proteinprep 开销被放大 |

#### Uni-Dock 关键纠正

服务器安装的 Uni-Dock v1.1.3 支持：

```bash
--paired_batch_size
```

但不能把 Uni-Dock 优化理解成只有这一个参数。官方 README 和服务器 `unidock --help` 显示，Uni-Dock 的入口要按任务形态选择：

| 任务形态 | 应优先评估的官方入口 |
|---|---|
| 一个 receptor / 一个 pocket 对大量 ligand | `--gpu_batch`；ligand 太多导致命令行过长时用 `--ligand_index` |
| 少量 ligand 对一个 pocket | 官方 FAQ 明确说明这种场景开销比例大，速度可能变慢；需要考虑合并批次或换任务组织方式 |
| 一条 protein 配一条 ligand 的多组一对一任务 | `--paired_batch_size` + JSON `--ligand_index`，符合 `paired_batching.schema.json` |
| 从 PDB/SDF 等常见格式出发，需要准备蛋白、配体、补氢或统一工作目录 | `unidocktools proteinprep`、`ligandprep`、`unidock_pipeline` |
| 只做打分或局部搜索 | `--score_only` / `--local_only`，并核对是否需要 maps 或 autobox |
| AD4 / Vina / Vinardo 不同打分函数 | 分别核对 `--scoring`、`--maps` 或 `--receptor` 的输入要求 |

同时必须核对 `--search_mode`、`--exhaustiveness`、`--max_step`、`--num_modes`、`--refine_step`、`--max_gpu_memory`、box 参数、输出目录和输入格式。还要区分两种 `--ligand_index` 语义：常规模式下通常是 ligand 路径列表；`paired_batch` 模式下是 protein-ligand-box JSON 配对表。不能只盯 `paired_batch`，也不能默认每个 receptor 单独启动一个低吞吐命令。

Q10 的 204 条 P450 对同一个 diosgenin，本质上可以组织成 204 个 receptor-ligand pair，也可以视作多 receptor 同底物的批量任务。后续应先做小批模式对比，而不是直接押注单一参数。至少对比：

1. 运行时间。
2. GPU 利用率。
3. docking score 是否和旧脚本一致或接近。
4. 异常正分是否仍然出现。
5. 输出文件能否被后续 PT cache 和模型推理稳定读取。

#### 后续硬规则

1. 不再把 Q10_EXP001 包装成正式理想流程。它保留为临时尝试和事故复盘证据。
2. 下一版建议新建 `Q10_EXP002`，不要覆盖 Q10_EXP001。
3. 正式运行前必须先查官方文档、本机 `--help`、已有最佳实践和当前日志。
4. ColabFold 大批量预测前必须评估公共 MSA 限流、本地 MSA、序列去重、分片并行、复用 AlphaFold / AlphaFill 结构等方案。
5. P450 结构预测默认没有 HEM/Fe。凡是对接、Fe-催化原子距离、口袋定位或湿实验打分，必须先明确 HEM/Fe 来源、补充方法和质量审计。
6. AlphaFill 补 HEM/Fe 要优先查完整数据库或候选 entry，记录模板来源、identity、local RMSD、clash、Fe-Cys 距离和质量标签。少数模板 overlay 只能标记为低置信度应急方案。
7. Uni-Dock 必须先按任务形态查官方入口，不能只盯 `paired_batch`。一个 pocket 对大量 ligand、一对一批量、格式准备、只打分、局部优化、不同 scoring function、不同 box/memory 设置，都要分别核对最合适参数。常规 `--ligand_index` 和 paired JSON `--ligand_index` 不是同一种输入语义，必须分清。
8. 昂贵服务器任务必须先做真实性能基准，报告 GPU 利用率、显存、每样本耗时、瓶颈和预计总耗时。发现限流、低利用率、残留进程或异常慢日志时先停下排查。

#### 记录位置

本次教训已经同步写入：

| 位置 | 用途 |
|---|---|
| `C:\Users\Administrator\.codex\memories\extensions\ad_hoc\notes\2026-05-29-q10-colabfold-alphafill-unidock-lessons.md` | 长期记忆 |
| `C:\Users\Administrator\.codex\AGENTS.md` | 全局协作硬规则 |
| `D:\EZSpecificity_Project\AGENTS.md` | 本项目硬规则 |
| 本 session_log | Q10 任务内事故复盘 |

### 13. 2026-05-29 Q10_EXP001 恢复并跑完第一版打分

用户决定先按原来的 Q10_EXP001 继续推进，拿到第一版完整排序；后续再用更严格流程做 EXP002，并考虑 EXP003 把构象选择从“能量优先”改成“Fe 与催化原子距离优先”。

#### 当前完成状态

Q10_EXP001 父流程从 Step 4 后的暂停状态恢复，继续执行：

1. Step 5：构建 204 条模型 PT cache。
2. Step 6：校验 PT cache。
3. Step 7：加载 Q2 random 最优 checkpoint 对 204 条候选打分。
4. Step 8：做 candidate_id 一致性审计。

最终流程显示：

```text
[Q10] Full pipeline finished.
```

#### 样本数和审计

| 环节 | 样本数 | 结果 |
|---|---:|---|
| 输入 enzyme candidates | 204 | 43 条 MDP450，161 条 ARATH |
| ColabFold 结构 | 204 | complete |
| HEM/Fe receptor | 204 | 每条均有 Fe |
| Uni-Dock 对接 | 204 | complete |
| PT cache | 204 | samples_written = 204 |
| 模型打分 | 204 | n_samples = 204 |
| 全流程一致性审计 | 204 | `ok=true` |

HEM/Fe receptor 质量标签：

| 标签 | 数量 |
|---|---:|
| good_fe_cys | 184 |
| bad_or_missing_fe_cys | 17 |
| poor_fe_cys | 3 |

底物图样本中 ligand heavy atom 数为 30，204 条均一致。

#### 输出路径

| 输出 | 服务器路径 |
|---|---|
| 全部模型打分 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/model_scores/q10_model_scores_all_test_msa_m1r1_v1.csv` |
| MDP450 单独排序 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/model_scores/q10_model_scores_MDP450_no_TM_removed_test_msa_m1r1_v1.csv` |
| ARATH 单独排序 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/model_scores/q10_model_scores_ARATH_uniprot_p450_test_msa_m1r1_v1.csv` |
| 模型打分 summary | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/model_scores/q10_model_scores_summary_test_msa_m1r1_v1.json` |
| 全流程审计 | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/results/audits/q10_full_output_audit_test_msa_m1r1_v1.json` |
| PT cache summary | `/root/autodl-tmp/EZSpecificity/PathD/P450/data/q10_wetlab_scoring/Q10_EXP001_diosgenin_20260528/model_cache/q10_unidock_pt_msa_m1r1_v1/manifests/q10_model_cache_summary_test_msa_m1r1_v1.json` |

模型使用：

| 项目 | 值 |
|---|---|
| checkpoint | Q2 EXP008 random best，`ep79-auc0.9316.ckpt` |
| device | `cuda:0` |
| use_gdtable | false |
| score_min | `1.93e-16` |
| score_max | `0.088214606` |
| score_mean | `0.002178010` |

注意：这里的 `model_score_sigmoid` 是排序分数，不是湿实验概率。Q10 没有标签，不能当作已校准概率解释。

#### Top 排名快照

全体 Top 10：

| rank | candidate_id | 来源 | model_score | docking_score | HEM/Fe 质量 |
|---:|---|---|---:|---:|---|
| 1 | ARATH_uniprot_p450_0014 | ARATH | 0.08821461 | -7.181 | good_fe_cys |
| 2 | ARATH_uniprot_p450_0109 | ARATH | 0.05679045 | -10.879 | good_fe_cys |
| 3 | MDP450_no_TM_removed_0021 | MDP450 | 0.04018906 | -9.714 | good_fe_cys |
| 4 | MDP450_no_TM_removed_0003 | MDP450 | 0.03694136 | 8.573 | good_fe_cys |
| 5 | MDP450_no_TM_removed_0026 | MDP450 | 0.02905585 | -9.422 | good_fe_cys |
| 6 | ARATH_uniprot_p450_0149 | ARATH | 0.02303494 | 5.207 | good_fe_cys |
| 7 | ARATH_uniprot_p450_0148 | ARATH | 0.01653275 | -6.715 | good_fe_cys |
| 8 | ARATH_uniprot_p450_0016 | ARATH | 0.01587803 | -11.210 | bad_or_missing_fe_cys |
| 9 | ARATH_uniprot_p450_0129 | ARATH | 0.01537628 | -6.990 | good_fe_cys |
| 10 | ARATH_uniprot_p450_0138 | ARATH | 0.01487731 | -7.757 | good_fe_cys |

MDP450 Top 5：

| rank | candidate_id | model_score | docking_score | HEM/Fe 质量 |
|---:|---|---:|---:|---|
| 1 | MDP450_no_TM_removed_0021 | 0.04018906 | -9.714 | good_fe_cys |
| 2 | MDP450_no_TM_removed_0003 | 0.03694136 | 8.573 | good_fe_cys |
| 3 | MDP450_no_TM_removed_0026 | 0.02905585 | -9.422 | good_fe_cys |
| 4 | MDP450_no_TM_removed_0029 | 0.01166793 | 6.547 | good_fe_cys |
| 5 | MDP450_no_TM_removed_0035 | 0.00124703 | 88.971 | good_fe_cys |

ARATH Top 5：

| rank | candidate_id | model_score | docking_score | HEM/Fe 质量 |
|---:|---|---:|---:|---|
| 1 | ARATH_uniprot_p450_0014 | 0.08821461 | -7.181 | good_fe_cys |
| 2 | ARATH_uniprot_p450_0109 | 0.05679045 | -10.879 | good_fe_cys |
| 3 | ARATH_uniprot_p450_0149 | 0.02303494 | 5.207 | good_fe_cys |
| 4 | ARATH_uniprot_p450_0148 | 0.01653275 | -6.715 | good_fe_cys |
| 5 | ARATH_uniprot_p450_0016 | 0.01587803 | -11.210 | bad_or_missing_fe_cys |

#### 结果解释边界

Q10_EXP001 已经形成完整第一版排序，但必须带风险标签使用：

1. ColabFold 批量阶段依赖在线 MSA，曾出现限流，效率流程不合格。
2. HEM/Fe receptor 是少数模板 overlay，AlphaFill 查询不够完整。
3. Uni-Dock 使用逐 receptor 单独命令，没有按任务形态评估官方高吞吐模式。
4. docking score 存在异常正分，例如 MDP450 Top 5 里有 `88.971`，说明对接构象或 receptor/box 需要进一步异常筛查。
5. 模型分数整体偏低，适合作为候选排序线索，不适合作为湿实验概率或最终推荐。

因此：Q10_EXP001 可以作为“第一版可交付排序”和后续严格版对照；正式汇报或交给老师时，应明确写成 preliminary ranking。

#### 后续安排

| 实验 | 目标 |
|---|---|
| Q10_EXP002 | 按新硬规则重做：先查官方文档和 `--help`，更完整查询 AlphaFill，按任务形态评估 Uni-Dock 入口，小批基准后再全量 |
| Q10_EXP003 | 在 EXP002 的基础上，尝试把对接构象选择从能量优先改为 Fe 与催化原子距离优先 |

### 2026-05-29 追加：两个酶列表完整打分表已整理

为方便直接查看两个输入酶列表的排序，已把服务器上的 Q10_EXP001 打分 CSV 下载到本地，并整理成中文 Markdown：

`D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathD_2026-05-08_导师反馈13问题\results\q10_wetlab_scoring\Q10_EXP001_diosgenin_20260528\model_scores\Q10_EXP001_两个酶列表完整打分表.md`

这份表包含：

1. `MDP450_no_TM_removed` 的 43 条完整排名。
2. `ARATH_uniprot_p450` 的 161 条完整排名。
3. 当前打分系统的字段解释。
4. 服务器原始 CSV 和摘要 JSON 路径。

复核提醒：Uni-Dock 的 pose 已进入模型输入，所以对接步骤不是无效步骤；当前没有把 `docking_score` 数值作为最终排名，是因为 EXP001 中存在异常正分，需要后续严格版再处理。
