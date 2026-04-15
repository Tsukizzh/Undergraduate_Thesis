# Session 10 — EXP004 论文基线外部评估（准备阶段）

**日期**: 2026-04-15
**状态**: 无 GPU 准备阶段全部完成 ✅，等用户开 GPU 跑推理

## 一、目标

拿论文训练好的 checkpoint（`saved_model/model/run_0/models/best-checkpoint.ckpt`）在我们 P450 allfix_unified 测试集上推理，得到 **外部基线 AUC**。为了公平，测试集中必须剔除论文训练时见过的酶（389 个 ESIBank P450）。

- **对比参照**：我们自训的 EXP001_allfix_unified (bare 28 维) Test AUC = **0.9320**
- **论文模型架构**：与我们 EXP001 字节级一致（bare 28 维、hidden_dim=128、k=48、3 EGNN、8 头 cross-attention）
- **唯一差别**：权重 + edge ordering（论文训的是 legacy_bug，我们训的是 fixed）

## 二、准备阶段 5 步（都在无 GPU 模式下完成）

### Step 1 — 建黑名单（本地 CPU）

**脚本**: `scratch/build_paper_blacklist.py`

**输入**:
- 本地 `389个P450的PDB映射_完整版.csv`（ESIBank 训练集全部 P450 UniProt）
- 服务器 `P450/data/Enzymes.csv`（我们 1622 个 P450 酶）

**结果**:

| 项目 | 数量 | 说明 |
|---|---|---|
| ESIBank P450 全集 | 389 | 完备黑名单，不存在隐藏 P450 |
| 我们 P450 enzymes | 1622 | 1457 有真 UniProt + 165 合成 ID (`ENZ_G*`) |
| **交集（需拉黑）** | **356** | 91.5% of 389, 21.9% of 1622 |
| ESIBank 不在我们数据里 | 33 | 论文训过但我们没有，无泄漏风险 |

详细报告：`paper_blacklist_report.json`

### Step 2 — 建过滤缓存 overlay（服务器）

**脚本**: `scratch/build_paperfilter_cache.py`

**新目录**: `/root/autodl-tmp/.../data/pt_cache_allfix_unified_paperfilter/random/`

```
paperfilter/random/
├── enzymes      → symlink → ../../pt_cache_allfix_unified/random/enzymes
├── substrates   → symlink → ../../pt_cache_allfix_unified/random/substrates
├── manifest.pt  → symlink → ../../pt_cache_allfix_unified/random/manifest.pt
├── train        → symlink → ../../pt_cache_allfix_unified/random/train
├── val          → symlink → ../../pt_cache_allfix_unified/random/val
└── test/
    ├── samples  → symlink → ../../../pt_cache_allfix_unified/random/test/samples
    └── index.pt ← 唯一新文件（~几百 KB），对 enzyme_ids ∈ blacklist 的行做 boolean mask
```

**过滤结果**:

| 项目 | 原始 | 过滤后 | 丢失 |
|---|---|---|---|
| 样本数 | 10999 | **7963** | 3036 (**27.6%**) |
| 唯一酶 | 1473 | 1125 | 348 |

原 cache **一个字节都没动**，守住"加法不减法"原则。详细报告：`filter_report.json`。

### Step 3 — 缓存端到端验证（服务器）

**脚本**: `scratch/verify_paperfilter_cache.py`（4 组检查）

| 检查 | 内容 | 结果 |
|---|---|---|
| **A** enzyme_id 映射 | 12 个 enzyme_global_id 的 Enzymes.csv 序列长度 vs enzymes.bin 存储长度 | 12/12 OK |
| **B** index.pt 自洽 | 5 个数组等长 (7963)，0 黑名单泄漏，所有 substrate_id/enzyme_id 在 lookup 里存在 | PASS |
| **C** PtCacheDataset 加载 | 17 个混合位置（边界+分位+随机）通过真实 loader 加载，edge_mode=legacy_bug，无 NaN | 17/17 PASS |
| **D** 原缓存未改动 | 所有 overlay 条目是 symlink（test/index.pt 除外） | PASS |
| **OVERALL** | — | **PASS ✓** |

详细报告：`verify_report.json`

### Step 4 — 论文 ckpt 预检（本地 CPU）

**脚本**: `scratch/preflight_ckpt_load.py`

**环境**: `D:\anaconda3\envs\torch` (Python 3.12.5, torch 2.3.0, pytorch_lightning **1.9.0**，与论文训练时完全一致)

**做法**:
1. 实例化 `SS(config)` 使用论文的 `complete-full-random-all-0-complex.yml`
2. `torch.load(paper_ckpt)` → 提取 `state_dict`
3. key diff: `model.state_dict() ∩ ckpt.state_dict()`
4. shape check: 所有匹配 key 的 tensor shape 逐一比对
5. `model.load_state_dict(ckpt_sd, strict=True)` 必须不抛异常

**结果**:

| 项目 | 值 |
|---|---|
| ckpt state_dict keys | 76 |
| model state_dict keys | 76 |
| **Matched** | **76/76** |
| Missing | 0 |
| Unexpected | 0 |
| Shape mismatches | 0 |
| **`strict=True` load** | **PASSED** |

ckpt metadata: `epoch`, `global_step`, `pytorch-lightning_version`, `state_dict`, `loops`, `callbacks`, `optimizer_states`, `lr_schedulers`（标准 PL checkpoint 结构）。

样本 keys 检查通过：`graph_net.protein_atom_emb.weight`, `graph_net.encoder_complex.net.0.edge_mlp.net.0.weight`, ... —— 命名完全一致，无需 prefix rewriting。

详细报告：`ckpt_preflight_report.json`

**⚠️ 重要背景**：`main_training_pt.py` 的 `test_evaluate` 内部用 `strict=False` 加载 ckpt。这次预检的意义是在 strict=True 下 **主动** 验证一次 key/shape 完全匹配，确认没有任何权重静默丢失——消除"AUC 看起来合理但模型实际是随机初始化"的隐形失败。

### Step 5 — 建 EXP004 实验目录（服务器）

**位置**: `/root/autodl-tmp/.../experiments/EXP004_paper_baseline_unified/`

**建立方式**:
```bash
cp -a EXP001_allfix_unified EXP004_paper_baseline_unified
rm -rf results/checkpoints/* results/test_eval.json results/metrics.csv results/fig_*.png logs/*
cp paper_best-checkpoint.ckpt results/checkpoints/paper_best-checkpoint.ckpt      # 22MB, md5 f4d87ea...
ln -sf ../../../EXP001_allfix_unified/.../pt-EXP001_...-ep43-auc0.9250.ckpt \
       results/checkpoints/ours_EXP001_ep43_auc0.9250.ckpt
```

`src/` 从 EXP001 完整拷贝，`scripts/main_training_pt.py` + `scripts/pt_dataset.py` 未修改（preflight 已确认兼容）。

**两个新脚本**:

`scripts/run_train.sh`（主结果，单次）:
```bash
python ... --test-only --checkpoint paper_best-checkpoint.ckpt \
  --edge-mode legacy_bug --cache-dir .../paperfilter/random \
  --output-json results/test_eval_paper_legacy_filtered.json
```

`scripts/run_test_grid.sh`（4 路诊断网格）:

| # | ckpt | edge_mode | 意义 | 输出 |
|---|---|---|---|---|
| 1 | paper | legacy_bug | **主结果**：论文模型在我们数据上的表现 | `test_eval_paper_legacy_filtered.json` |
| 2 | paper | fixed | paper 权重对 edge 的敏感度 | `test_eval_paper_fixed_filtered.json` |
| 3 | ours EXP001 | legacy_bug | 对称控制：我们的权重应偏好 fixed | `test_eval_ours_legacy_filtered.json` |
| 4 | ours EXP001 | fixed | 我们 native 模式在过滤子集上（对比原始 0.9320） | `test_eval_ours_fixed_filtered.json` |

## 三、最终目录结构

```
/root/autodl-tmp/EZSpecificity/PathC/P450/
├── data/
│   ├── pt_cache_allfix_unified/random/           [原始，零改动]
│   └── pt_cache_allfix_unified_paperfilter/random/
│       ├── enzymes → ../../pt_cache_allfix_unified/random/enzymes
│       ├── substrates → ...
│       ├── manifest.pt → ...
│       ├── train → ...
│       ├── val → ...
│       └── test/
│           ├── samples → ...
│           └── index.pt ← 新文件（7963 行 × 5 数组）
│
└── experiments/
    ├── EXP001_allfix_unified/                     [原始，零改动]
    ├── ...
    └── EXP004_paper_baseline_unified/             [新]
        ├── README.md
        ├── src/                                   [EXP001 拷贝]
        ├── configs/config.yml                    [未改，tag 留 EXP001]
        ├── scripts/
        │   ├── main_training_pt.py               [EXP001 拷贝，未改]
        │   ├── pt_dataset.py                     [EXP001 拷贝，未改]
        │   ├── run_train.sh                      [新：test-only, paper + legacy]
        │   └── run_test_grid.sh                  [新：4 路诊断]
        ├── results/
        │   └── checkpoints/
        │       ├── paper_best-checkpoint.ckpt    [22MB, 论文 ckpt 上传]
        │       └── ours_EXP001_ep43_auc0.9250.ckpt → symlink
        └── logs/                                  [空]
```

## 四、代码产物清单

| 路径 | 类型 | 说明 |
|---|---|---|
| `scratch/build_paper_blacklist.py` | 本地 CPU | 389 vs 1622 交集，生成 blacklist JSON |
| `scratch/build_paperfilter_cache.py` | 服务器 CPU | overlay cache + filtered index.pt |
| `scratch/verify_paperfilter_cache.py` | 服务器 CPU | 4 组端到端验证（A/B/C/D） |
| `scratch/preflight_ckpt_load.py` | 本地 CPU | strict=True key+shape 检查 |
| `scratch/exp004_run_train.sh` | 服务器模板 | EXP004 单次主运行 |
| `scratch/exp004_run_test_grid.sh` | 服务器模板 | EXP004 4 路诊断 |
| `scratch/exp004_README.md` | 文档 | EXP004 实验说明 |

## 五、所有报告本地归档（本 session 目录下）

- `paper_blacklist_report.json` — 389/1622 交集详情
- `filter_report.json` — 过滤前后计数
- `verify_report.json` — 缓存 4 组检查结果
- `ckpt_preflight_report.json` — ckpt 预检结果

## 六、待执行（等 GPU）

1. 用户开 GPU（4×RTX4090）
2. `cd /root/autodl-tmp/.../experiments/EXP004_paper_baseline_unified && bash scripts/run_test_grid.sh`
3. 4 个 `test_eval_*.json` 产出后，回填到本 session log 并做 4 路对比分析

**预计耗时**: 每次 ~3-5 min，全部 4 次 ~15-20 min

## 七、限制声明（未来报告时必须提的）

1. **黑名单 UniProt-only 匹配**：我们只有 389 个 ESIBank P450 的 UniProt ID，没有序列，无法做序列 hash 兜底。165 个合成 ID (`ENZ_G*`) 的酶里可能隐藏着和 33 个"不在我们数据的 ESIBank UniProt"对应的同一蛋白——最多漏 33 个
2. **只过滤 enzyme，不过滤 substrate**：论文 random_split 与我们 test 在底物层面有重叠，本实验不控制
3. **过滤率 27.6%**：略低于 30% 门槛，可作为 leakage-controlled 外部基线报告，但必须明确标注 "filtered subset"，不能省略为 "the external test"

## 八、多轮 codex 协作记录

每一步都和 codex 做了独立审查：

| 步骤 | codex 主要输入 |
|---|---|
| 方案整体可行性 | 确认 plan 成立，非破坏性 overlay 正确 |
| filter 结果审查 | 27.6% drop 可接受，8 个黑名单酶不在 test 是正常现象 |
| verify 脚本设计 | 加入 substrate lookup + enzyme lookup + 随机位置 + provenance 检查（D 组） |
| preflight 脚本设计 | 强调 shape check 必做、避免静默 prefix 重写、分离 SS 构造失败与 ckpt 载入失败 |
| EXP004 dir 设计 | 移除 test-only 下无用的 `--devices 4` `--max-epochs`，用 symlink 代替 copy |

所有审查结果都让我主动改进了脚本或方案，没有一步是盲走。
