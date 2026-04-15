# Session 11：EXP005 双图架构 Dualgraph 2+ 训练执行

**实验编号**：EXP005
**方案代号**：dualgraph 2+（residue backfill + g_res bypass）
**基线**：EXP001_allfix_unified Test AUC = **0.9320**
**起始日期**：2026-04-16
**完整计划文档**：[../../EXP005_双图架构完整计划.md](../../EXP005_双图架构完整计划.md)

---

## 一、总体架构

在基线 EXP001_allfix_unified 的原子级 SE(3)-不变 EGNN + 双向交叉注意力基础上，**新增一条残基级 GVP-GNN 通路**，双出口融合：

1. **出口 1（residue backfill）**：GVP 残基嵌入 `h_res` 按 `pocket_residue_idx` 注入回 `x_pro` 对应 UniProt 位置，再进交叉注意力 — **深融合**
2. **出口 2（g_res bypass）**：`scatter_mean(h_res)` 得到图级 `g_res`，作为第 8 个向量拼接到末端预测头，header 输入维度 896→1024 — **浅旁路**

**关键设计**：两个新加路径都**零初始化**（`h_res_proj` 最后一层 Linear + `specificity_header` 第一层 weight 的 new 128 列块），训练 step 0 时 SSDualgraph 严格等价基线 SS；step 1 之后 GVP 分支通过零权重更新被激活，开始学习。

---

## 二、Phase 1：enzyme_resid_map 构建 ✅

### 问题溯源

单 pocket smoke test 扩展到 10 个后发现 **5/10 pocket_residue_idx 和 UniProt 序列位置对不上**。追本地 + 服务器数据生成链全部字节级 copy，没改过 resid。codex 独立 review 确认：bug 不在代码里，而在**源头 PDB 文件本身的编号约定异质性**：

- `alphafold_heme_transplant_v2` (677) / `alphafill` (509) / `pcpd_predicted` (大部分, 246)：delta = 0，UniProt 编号
- `experimental_pdb` (103)：作者编号，N 端 Met 常被切除 → delta ≠ 0
- 少数 alphafill / pcpd 也存在异常（如 dock=30789 AlphaFill delta=+68，dock=985 PCPD delta=+1）

### 方案 D（codex 最终推荐）：Pairwise Alignment + SIFTS 交叉验证

- 对 **1479 个酶**（pt_cache_allfix_unified 的 sidecar 集合）逐个做 BioPython 全局比对（match=2, mismatch=0, open=-10, extend=-0.5）
- `experimental_pdb` 的 96 个还原 SIFTS REST API (https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/) 做独立交叉验证
- 多链处理：同源二聚体（如 1JPZ）同时映射 A+B 链

### 结果（2026-04-16 早）

| tier | count | 条件 |
|---|---|---|
| gold | 1471 | SIFTS 100% 一致 OR identity≥0.99 AND coverage≥0.95 |
| trusted | 6 | identity≥0.90 AND coverage≥0.90 |
| partial | 2 | identity≥0.75 AND coverage≥0.75 |
| suspect | 0 | — |
| **total** | **1479 / 1479 = 100%** | |

**扰动稳定性测试**（codex round 3 补充要求）：对 8 个非 gold 酶分别测 5 组不同 gap penalty 配置，**全部 0 位置变动** — 全局对齐在 79% identity 下依然局部唯一。

### 六轮审计验证

| 审计 | 范围 | 结果 |
|---|---|---|
| A | `split("\|")` 回环保 icode 空格 | ✓ |
| B | BioPython `aln.aligned` 半开 0-indexed | ✓ codex 确认 |
| C | parser-consumer key 一致性 | ✓ 100% 覆盖隐含 |
| E | sidecar / index.pt / CSV 三方对齐 | **44090/44090 × 3 splits 0 mismatch** |
| F | Enzymes.csv vs ESM cache 源漂移 | ✓ 时间戳无漂移 |
| G | **全量 44090 pocket sanity** | **2,991,586 / 2,991,586 残基 100% 覆盖, 99.97% aa match** |
| H | 5 组 gap penalty 扰动稳定性 | **0 位移** |

产出：
- `pt_cache_dualgraph_allfix_unified/enzyme_resid_map.pt` (14.5 MB)
- `pt_cache_dualgraph_allfix_unified/sifts_cache/*.json` × 96
- `pt_cache_dualgraph_allfix_unified/enzyme_alignment_report.md`
- `pt_cache_dualgraph_allfix_unified/random/{train,val,test}/dock_sidecar.pt`

---

## 三、Phase 2：GVP 特征抽取 ✅

### 核心脚本

`scratch/exp005_build_gvp_cache.py`，对 **44,090 个 training dock** 抽取：

```python
{
    "node_s":             [N, 6]     float16   # phi/psi/omega cos/sin
    "node_v":             [N, 3, 3]  float32   # forward/backward/sidechain 单位向量
    "edge_index":         [2, E]     int32     # kNN k=30, src=neighbor, dst=self
    "edge_s":             [E, 32]    float16   # 16 RBF 距离 + 16 位置编码
    "edge_v":             [E, 1, 3]  float32   # CA-CA 单位向量
    "pocket_residue_idx": [N]        int64     # 0-indexed UniProt 位置 (Phase 1 映射)
    "aa_type":            [N]        uint8     # AA 类型 0-19
}
```

### 关键 bug 发现 + 修复（codex 4 轮审核）

1. **Round 1** — BioPython `aligned` 0-indexed / tuple_cat dim 语义：codex 最初质疑 tuple_cat 错误，我 push back（对 2D 标量 + 3D 向量的情况 dim=1 自动对齐），codex 自我纠正
2. **Round 2** — dtype 安全：fp16 `node_s/edge_s` 进 fp32 Linear 会报错，在 encoder 入口统一 `.float()`
3. **Round 3** — **edge 方向 bug**：我原先 `src=arange, dst=nbr`，和 torch_cluster 默认 `flow='source_to_target'` 相反 → 每个中心变成"reverse-kNN 聚合"而非"kNN 聚合"。修正为 `src=neighbor, dst=self`，和 EnzymeCAGE 原版一致
4. **Round 4** — 性能：`OMP_NUM_THREADS=1` + `torch.set_num_threads(1)` 避免 20 worker 抢核

### 速度对比

| 配置 | 44090 samples 耗时 |
|---|---|
| OMP 默认（多线程争用）| ~10 分钟（1778 files / 3 min 推算）|
| OMP=1（每 worker 单线程）| **8.6 秒** ← 70 倍提速 |

### 结果

| 状态 | count |
|---|---|
| ok | **44,026 (99.85%)** |
| empty_pocket | 60 |
| too_few_residues=1 | 4 |
| 总残基数 | 1,627,439 |
| 总边数 | 47,033,436 |
| 映射覆盖 | **1,627,439 / 1,627,439 = 100.00%** |

### 全量 sanity（11 类检查，codex 推荐）

| 检查 | 结果 |
|---|---|
| schema (9 keys) | 44026/44026 ✓ |
| dtypes (float16/float32/int32/int64/uint8) | ✓ |
| count consistency (N × 5 places, E × 3 places) | ✓ |
| graph semantics (no self-loops, E=N*k, in-degree uniform) | ✓ |
| n_mapped == N | ✓ |
| N ≥ 2 | ✓ |
| aa_type range [0,19] | ✓ |
| pocket_residue_idx within [0, seq_len) | ✓ |
| NaN/Inf free | ✓ |
| PDB re-parse cross-check | **2,991,586 / 2,991,586 100% 位置正确** |

### 64 失败样本的 Option 3 fallback（codex 推荐）

不生成占位 .pt 文件，**单独建 `gvp_invalid_docks.pt` manifest** 含 per-dock 元数据 `{reason, enzyme_id, split, uniprot}`。Dataset `__getitem__` 查 set 发现 invalid dock 就合成 1-node-0-edge 占位 + `gvp_valid=False`，模型 forward 时通过 mask 跳过。

失败分布：train 37 / val 16 / test 11，主要来自 enzyme 1334 (Q9SB21) 这类结构有问题的酶。

---

## 四、Phase 3：模型代码集成 ✅

### 目录结构（严格只在 EXP005 内新增 / 修改，不动其他 src）

```
experiments/EXP005_dualgraph_2plus_allfix_unified/
├── configs/config.yml                          [改: +model.gvp + tag]
├── scripts/
│   ├── main_training_pt.py                     [改: import + DDP find_unused + cache 验证]
│   ├── pt_dataset_dualgraph.py                 [新]
│   └── run_train.sh                            [新/改]
└── src/Models/
    ├── ss_dualgraph.py                         [新] 继承 SS, forward 里插入 GVP
    └── Structure/gvp.py                        [新] 移植 EnzymeCAGE GVP encoder
```

### ss_dualgraph.py 关键设计

```python
class SSDualgraph(SS):
    def __init__(self, config):
        super().__init__(config)
        self.gvp_encoder = GVP_embedding(
            node_in_dim=(6, 3), node_h_dim=(128, 16),
            edge_in_dim=(32, 1), edge_h_dim=(32, 1),
            seq_in=True, num_layers=3, drop_rate=0.1)
        # GVP_embedding 输出 [N, 2*128=256]
        self.h_res_proj = Sequential(LayerNorm(256), Linear(256, 128))
        nn.init.zeros_(self.h_res_proj[-1].weight)   # ★ 零初始化
        nn.init.zeros_(self.h_res_proj[-1].bias)
        self.g_res_proj = Linear(256, 128)
        self.h_res_gate_logit = Parameter(tensor(-2.2))  # sigmoid ≈ 0.10
        # 扩展 specificity_header.mlp.net[0] 从 896 → 1024
        new_first[:, :896] = old_first.weight.clone()
        new_first[:, 896:] = 0  # ★ 零初始化 g_res 列块
    
    def forward(self, G):
        x_pro = self.protein_mlp(G.embedding).view(-1, 1450, hidden_dim)
        # GVP encode
        h_res = self.gvp_encoder((G.gvp_node_s, G.gvp_node_v), G.gvp_edge_index,
                                  (G.gvp_edge_s, G.gvp_edge_v), seq=G.gvp_aa_type)
        # h_res → x_pro injection with scatter_mean duplicate aggregation
        valid = (G.gvp_pocket_residue_idx >= 0) & G.gvp_valid[G.gvp_node_s_batch]
        key = G.gvp_node_s_batch[valid] * 1450 + G.gvp_pocket_residue_idx[valid]
        unique_keys, inverse = key.unique(return_inverse=True)
        h_agg = scatter_mean(h_res_proj_valid, inverse, dim=0, dim_size=unique_keys.numel())
        delta = (torch.sigmoid(h_res_gate_logit) * h_agg).to(x_pro.dtype)
        x_pro[sid_unique, rix_unique] += delta
        # ... 跑 baseline EGNN + cross-attention + pooling
        g_res = scatter_mean(h_res, G.gvp_node_s_batch, dim=0, dim_size=B)
        g_res = g_res * G.gvp_valid.float().unsqueeze(-1)  # 屏蔽占位样本
        embeddings.append(self.g_res_proj(g_res))  # 8th vector
        return self.specificity_header(embeddings)  # 1024 → 128 → 1
```

### pt_dataset_dualgraph.py 设计（codex 2 轮审核）

- `DualgraphData(StructureSequenceData)` 重写 `__inc__` / `__cat_dim__`，`gvp_edge_index` offset by node count，`gvp_pocket_residue_idx` **NO offset**（绝对 UniProt 位置）
- `PtCacheDualgraphDataset(PtCacheDataset)` 重写 `__getitem__`，基线 sample + GVP 字段合并到 DualgraphData
- `follow_batch=['ligand_index', 'gvp_node_s']` → PyG 自动生成 `batch.gvp_node_s_batch`
- 失败样本合成占位 dict 而非加载文件
- `preload_gvp=True` 将 350 MB GVP cache 常驻内存（每 rank）

### codex 多轮审核（共 7 轮）

| 轮次 | 主题 | 发现 |
|---|---|---|
| 1 | 架构选型 Q1/Q2/Q3 | Option 1A (subclass) + Option 3 (manifest) + 部分注入掩码 |
| 2 | PyG idioms | scatter_mean + unique + inverse 正确写法 |
| 3 | forward 具体位置 | specificity_header.mlp.net[0] 替换 + AMP dtype 安全 |
| 4 | tuple_cat push back | 我对，codex 自我纠正 |
| 5 | 7 个 ss_dualgraph 检查 | delta = `.to(x_pro.dtype)` 安全补丁 |
| 6 | 路径 bug | 我自己发现 + 用户提醒双重触发 |
| 7 | 最终 review + --test-only 安全检查 | 拒绝加载 legacy ckpt 到 SSDualgraph |

---

## 五、L1-L7 渐进式 smoke test ✅

| Level | 内容 | 结果 |
|---|---|---|
| **L1** | import + syntax compile | ✓ 所有新模块 + main_training_pt.py 编译通过 |
| **L2** | SSDualgraph(config) 构造 + 参数量 | **2,684,654 params** (+838K vs baseline 1,846,660), header 1024, 零初始化验证 |
| **L3** | dataset `__getitem__` 单样本 (valid + placeholder) | 两种 schema 都正确 |
| **L4** | batch collate + gvp_node_s_batch | 4 样本 (3 valid + 1 placeholder) 完美 offset：node 0-83 range 按 [0,40,65,84] 切分，pocket_residue_idx 保持绝对位置 [-1, 489] |
| **L5** | CPU forward + baseline equivalence | **max abs diff 4e-8** 和基线 SS 等价（只有 FP 浮点噪声）— 证明零初始化严格等价 |
| **L6** | loss + backward + GVP 梯度流 | **延迟解冻**：step 0 `gvp_encoder/g_res_proj/alpha_logit` 梯度为 0（零权重切断路径），step 1 全部激活（0.011 / 0.003 / -5e-6）。`h_res_proj` 和 `header[896:]` 从输出端接收梯度一直活跃 |
| **L7** | 完整 pl.Trainer mini-fit (1 GPU, 3 train + 4 val batch) | fit/test 全路径通过，`auc/val` logging 正确，checkpoint save/load 含 GVP keys，header 第一层 shape (128, 1024) ✓ |

---

## 六、Phase 4：正式训练（进行中）

### 第一次尝试 — `--preload` + 4 卡 → OOM

配置：bs=88 × devices=4 × workers=6 × `--preload`

**崩溃点**：Epoch 0 结束后 val phase，`rank: 3 Child process with PID 2053 terminated with code -9`（SIGKILL）

**根因分析**（和 codex 对话 round 8）：
1. `--preload` 让每 rank 主进程加载 ~7 GB dict-of-tensors
2. DataLoader 6 workers 通过 fork 继承这 ~7 GB，Python refcount COW 每次访问写 object header → COW 页面逐步被复制
3. 24 workers × 7 GB = 168 GB worker RSS
4. 同时 DataLoader prefetch_factor=4 × bs=88 × PyG Batch 入 /dev/shm 队列 → 225 GB shared tmpfs 超过 180 GB /dev/shm 限
5. cgroup 360 GB 限被击穿 → OOM killer

### 第二次尝试 — 去掉 `--preload`，严格对标 EXP001_allfix_unified → 🟢 运行中

```bash
python main_training_pt.py \
  --config configs/config.yml \
  --cache-dir pt_cache_dualgraph_allfix_unified/random \
  --edge-mode fixed \
  --batch-size 88 \
  --max-epochs 200 \
  --devices 4 \
  --num-workers 6 \
  --run-name EXP005_dualgraph_2plus_allfix_unified \
  --shutdown
```

**配置和 EXP001_allfix_unified 完全一致**（唯一区别：cache 路径 + 模型类，没有 `--preload`）。

### DDP `find_unused_parameters` 必须为 True

L6 发现的延迟解冻效应在 DDP 下会触发 `RuntimeError: It looks like your LightningModule has parameters that were not used in producing the loss`。修复：`DDPStrategy(find_unused_parameters=True, gradient_as_bucket_view=True)`。代价约 5-10% 速度，但不可避免。

### 运行状态（2026-04-16 06:05 UTC，训练中）

```
Epoch 14 完成:
  train_loss: 0.171, val_loss: 0.209
  val_auc: 0.874, val_aupr: 0.456
  
Loss / AUC 曲线:
  ep1:  tl=0.324  vl=0.293  
  ep5:  tl=0.302  vl=0.285  auc=0.591
  ep8:  tl=0.261  vl=0.254  auc=0.711
  ep10: tl=0.218  vl=0.233  auc=0.800
  ep12: tl=0.199  vl=0.215  auc=0.842
  ep14: tl=0.171  vl=0.209  auc=0.874

速度:       ~62 s / epoch (warmup 后 ~50 s / epoch)
  Train:     63 batch × ~1.85 it/s = 34 s
  Val:       29 batch × ~2.00 it/s = 14 s
  overhead:  checkpoint + metric + init ~14 s
GPU util:   99-100% × 4 卡
GPU mem:    26-27 GB / 32 GB per 5090
RAM:        165 GB used, 144 shared, 479 available (健康)
```

### 预计

- 早停 epoch：对标 EXP001 约 ~43-58 epoch
- 总训练时长：**~45-60 min**
- 加 auto-test + shutdown：**~50-70 min 总**

---

## 七、关键决策与教训

### 决策记录

1. **Option 2+ 融合方案**：residue backfill + g_res bypass，对称 EZSpecificity 原架构的 `x_str` + `str_mean` 双出口哲学
2. **Dataset 改动范围**：Option 1A 子类化 `PtCacheDataset`，不改基线数据流
3. **失败样本处理**：Option 3 manifest + 运行时合成，不生成占位 .pt 污染 cache
4. **零初始化保守启动**：`h_res_proj[-1]` + `specificity_header[896:]` 双零初始化，step 0 严格等价基线 SS
5. **preload 放弃**：按 allfix 系列惯例不开 `--preload`，避开 fork COW + /dev/shm 坑
6. **find_unused_parameters=True**：零初始化导致的 step 0 梯度断流必须 DDP 层容忍

### 教训

1. **单 pocket smoke test 错觉**：早期 dock 3/13/14/15/17 恰好全是 alphafill（delta=0），扩大到 10 pocket 才暴露 5/10 错位。sanity 样本必须跨所有数据源。
2. **DDP 的不可见代价**：单 GPU smoke test 不能测出 DDP `find_unused_parameters=False` 的严格检查。必须在真实 DDP 下跑一小段才能暴露。
3. **fork COW 成本**：dict-of-tensors 大对象 preload + fork 是 Python 的反模式，refcount 写入触发 COW 会逐步复制到每个 worker，累积击穿 cgroup。长期方案是 mmap/flatbin 式共享存储（codex Option 4）。
4. **日志 buffer**：DDP 下 Python stdout 有 buffer 延迟，第一次 batch 出来前看不到进度。要结合 `free -g` / `nvidia-smi` / `ps` 判断状态。

### 文件索引

| 文件 | 路径 | 作用 |
|---|---|---|
| gvp.py | `src/Models/Structure/gvp.py` | GVP encoder 移植 (335 行) |
| ss_dualgraph.py | `src/Models/ss_dualgraph.py` | SSDualgraph 模型 (180 行) |
| pt_dataset_dualgraph.py | `scripts/pt_dataset_dualgraph.py` | Dataset + DualgraphData batching (290 行) |
| main_training_pt.py | `scripts/main_training_pt.py` | 训练入口（改 9 处）|
| config.yml | `configs/config.yml` | +model.gvp 子节 + tag |
| run_train.sh | `scripts/run_train.sh` | 4 卡 DDP 启动脚本 |
| enzyme_resid_map.pt | `data/pt_cache_dualgraph_allfix_unified/enzyme_resid_map.pt` | 1479 酶 resid→UniProt 位置表 |
| gvp_cache/samples/ | `data/pt_cache_dualgraph_allfix_unified/gvp_cache/samples/` | 44026 个 gvp_{dock}.pt |
| gvp_invalid_docks.pt | `data/pt_cache_dualgraph_allfix_unified/gvp_cache/gvp_invalid_docks.pt` | 64 失败样本 manifest |

---

## 八、下一步

等训练完成后：
1. 下载 `results/checkpoints/` + `metrics.csv` + `test_eval.json` + `train_live.log`
2. 对比 EXP005 Test AUC vs 基线 0.9320
3. 画图：训练曲线 + val AUC 比对
4. 分析 GVP 分支是否真的贡献了信号（可选消融：关掉 h_res 注入或 g_res 旁路）
5. 写入 MEMORY.md 和项目进度日志

Auto-shutdown 会在 test 完成后 30 秒关机，实例停止计费。
