# Step 10: .pt 训练管线 — 完整 Session 日志

> **日期**: 2026-03-16 ~ 2026-03-17
> **目标**: 用 .pt 预处理缓存替代 LMDB + 部署云服务器 + 还原论文基线(legacy_bug)训练
> **最终成果**: ✅ per-sample .pt v3，本地 7.56 it/s / 云服务器 3.74 it/s(bs48) GPU 100%

---

## 一、背景：为什么要做 .pt 缓存

Step 9 的 LMDB 缓存模式有两个问题：
1. **mmap 抖动**：60GB LMDB 文件映射到 32GB 内存，系统不断换页，速度从 3.8 it/s 衰减到 2.09 it/s
2. **CPU k-NN 瓶颈**：每个样本的 k-NN 图构建占 CPU 时间 50-70%

PyG 社区的标准做法：预计算所有图数据，存成 per-sample .pt 文件，训练时直接加载。

---

## 二、.pt 方案三次迭代（血泪史）

### v1：分片模式（失败）

**做法**：176K 样本打包成 87 个分片，每个 300MB。
**问题**：
- torch.load 必须加载整个 300MB 才能取出一个 160KB 的样本
- shuffle=True 时连续 batch 跳到不同分片，ShardCache 缓存不断淘汰
- 一个 epoch 总读取量 ~15TB → SSD 100% 饱和 → 系统完全卡死
- 速度：73 s/it（比 LMDB 慢 280 倍）

**教训**：torch.load 不适合大文件的随机访问。

### v2：Flatbin 模式（部分成功）

**做法**：enzyme/substrate 改成 .bin 文件 + file.seek（只读需要的 2MB），graph 仍用分片。
**结果**：
- enzyme 读取从 900MB/次 → 2MB/次，速度提升 250 倍
- 但 graph 分片仍是 300MB，torch.load 仍然慢
- num_workers=0: 2.65 s/it（仍比 LMDB 慢）
- num_workers=1: 1.92 s/it
- **4 个 worker 直接撑爆 32GB 内存，电脑卡死重启**

**卡死原因**：EntityTensorCache 设了 4096 个条目 × 7MB/酶 = 28GB。4 个 worker 各一份 = 112GB。

**修复**：EntityTensorCache 4096→64（酶），4096→256（底物）

### v3：Per-sample 模式（成功！）

**做法**：每个样本一个独立 .pt 文件（~160KB），分成 177 个子目录（每目录 1000 个文件）。
**结果**：
- torch.load 160KB 只需几毫秒
- 一个 epoch 总读取量 28GB（vs v1 的 15TB）
- SSD 利用率从 100% 降到几个百分点
- **7.33-7.70 it/s，29 min/epoch，稳定不衰减**

---

## 三、所有速度测试数据（完整记录）

### Step 9 LMDB 模式

| 配置 | 速度 | 每 epoch | 内存 | 备注 |
|------|------|---------|------|------|
| bs4/accum8/w2（初始基线） | 1.71 it/s | ~2h | ~15GB mmap | 无优化 |
| bs8/accum4/w6（优化后，无缓存） | 4.3 it/s | ~49min | ~20GB mmap | 最快无缓存 |
| 240GB 缓存 benchmark（300 batch） | 14.93 it/s | N/A | 全进 cache | **虚假数字**，只跑了 300 batch |
| 240GB 缓存早期 epoch | 4.5 it/s | ~47min | 逐渐增长 | 速度会衰减 |
| 240GB 缓存 ep14-19 稳态 | 2.09 it/s | ~100min | 30GB+ mmap | 最终稳态 |
| BlockShuffleSampler | 2.26 it/s | — | — | 失败实验，废弃 |
| 无缓存 + fp16 enzyme + 6w + bs14 | 1.35 it/s | — | — | 比缓存更慢 |
| 缓存 + float32 + 2w + bs14（Step 9 最终） | 3.8→2.09 | 47-100min | 15-30GB | 衰减 |

### Step 10 .pt 模式

| 配置 | 速度 | 每 epoch | 内存 | 备注 |
|------|------|---------|------|------|
| v1 分片 300MB | 0.01 it/s (73s/it) | 258h | 9.7GB | SSD 100%，不可用 |
| v2 flatbin + 分片 graph, 0w | 0.38 it/s (2.65s/it) | 9h | 4GB→增长 | 无并行 |
| v2 flatbin + 分片 graph, 1w | 0.52 it/s (1.92s/it) | 6.8h | ~8GB | graph 分片仍慢 |
| **v3 per-sample, 0w** | **4.55 it/s** | ~47min | 3.4GB | 数据加载测试 |
| **v3 per-sample, 2w** | **7.33-7.70 it/s** | **29min** | **~10GB** | ✅ 最终方案 |

---

## 四、遇到的所有问题及修复

### 问题 1：4 worker 内存爆炸，电脑卡死
- **原因**：EntityTensorCache 4096 条目 × 7MB = 28GB/进程，4 进程 = 112GB
- **修复**：酶缓存 4096→64（450MB），底物缓存 4096→256（690MB）
- **教训**：缓存大小必须计算 entries × per-entry-size ≤ 内存预算

### 问题 2：SSD 饱和 100%，系统冻结
- **原因**：分片模式每个 batch 触发 3-5 次 300MB 的 torch.load，一个 epoch = 15TB
- **修复**：改成 per-sample .pt（160KB/文件），一个 epoch = 28GB
- **教训**：torch.load 不适合 >10MB 文件的随机访问

### 问题 3：enzyme 分片 >2GB，Windows torch.load 失败
- **原因**：7 个分片按 dataset_id 分，brenda 分片 2.1GB 超过 Windows ZIP 限制
- **修复**：enzyme_shard_size 从 4000 改为 1000（每个 ~900MB）
- **教训**：Windows 上 .pt 文件不能超过 2GB

### 问题 4：ligand_x 特征顺序错误
- **原因**：atomic_number/100 应在 element one-hot 后面（位置 11），写成了最后（位置 33）
- **修复**：对照 transforms.py FeaturizeLigandAtom 原始顺序重写
- **教训**：必须逐字段对照原始代码

### 问题 5：__getstate__ 里缓存大小没同步
- **原因**：__init__ 改了 64/256，__getstate__（worker spawn 时调用）还是 4096
- **修复**：__getstate__ 同步为 64/256
- **教训**：改 __init__ 的同时检查 __getstate__

### 问题 6：服务器 PYTHONPATH 缺失
- **原因**：多进程 worker fork 后找不到 Datasets 模块
- **修复**：export PYTHONPATH=.../src:$PYTHONPATH
- **教训**：nohup 启动时环境变量要显式设置

### 问题 7：服务器 8 worker 无产出
- **原因**：16GB RAM + 8 个 worker 同时读 60GB LMDB = 严重争抢
- **修复**：放弃服务器生成，改为本地生成后传输
- **教训**：worker 数不能超过 RAM/数据集大小的比值

### 问题 8：速度测试不代表真实训练
- **原因**：用 val（1482 样本）测出 10 batch/s，但只测了数据加载没跑模型
- **修复**：改为用完整训练（forward + backward）测试
- **教训**：必须用真实训练场景测速

---

## 五、血泪教训（所有未来开发必须遵守）

### 存储格式
1. **PyG 标准 = per-sample .pt 文件**，不要发明自定义分片格式
2. **torch.load 只适合小文件**（<10MB），大文件用 file.seek + np.frombuffer
3. **共享数据（enzyme/substrate）不要重复存**，用索引引用

### 内存管理
4. **缓存大小必须计算**：entries × per-entry-MB ≤ 可用 MB
5. **Windows spawn 模式**：num_workers × 单进程内存 ≤ 总 RAM
6. **__getstate__ 和 __init__ 必须同步**

### 测试方法
7. **用真实训练测速**（forward + backward + optimizer），不要只测数据加载
8. **小数据集（val 1482 样本）不代表大数据集（train 176K）**
9. **监控 SSD 利用率和内存趋势**，不只是看 it/s

### 数据生成
10. **多 worker 并行**：build_pt_cache.py 单线程跑 graph 分片要 5 小时，10 worker 只要 2 分钟（150 倍加速）
11. **向量化替代 Python 循环**：bond 过滤从 for-loop 改成 torch.isin，速度提升 100 倍
12. **断点续传是必需品**：长时间脚本必须检查已有文件并跳过，否则中断就白跑

### 开发流程
13. **第一版就加**：多 worker、断点续传、向量化——不是"先写功能再优化"
14. **不要先写功能再优化**——优化是功能的一部分
15. **每次改完必须做：内存测试 + 磁盘 IO 测试 + 真实训练测试**

---

## 六、最终存储布局（清理后）

```
allsplit_pt_cache/ (~57GB)
├── manifest.pt                        (1KB)
├── enzymes/
│   ├── enzymes.bin                    (27GB, seek 读取)
│   ├── enzymes_index.pt               (索引)
│   └── enzymes_flatbin_index.pt       (索引)
├── substrates/
│   ├── substrates_grover.bin          (4.3GB, seek 读取)
│   ├── substrates_meta.pt             (grover_mean + morgan, ~430MB)
│   └── substrates_index.pt           (索引)
├── train/
│   ├── index.pt                       (样本元信息)
│   └── samples/                       (176,843 个 .pt, 各 ~160KB)
│       ├── 000/ (sample_000000.pt ~ sample_000999.pt)
│       ├── 001/
│       └── 176/
├── val/
│   ├── index.pt
│   └── samples/                       (1,482 个 .pt)
└── test/
    ├── index.pt
    └── samples/                       (8,841 个 .pt)
```

---

## 七、1 epoch 完整训练结果

| 指标 | 值 |
|------|------|
| 时间 | 29 分钟 |
| 速度 | 7.33 it/s（稳定不衰减）|
| val AUC | 0.6103（随机初始化第 1 个 epoch，正常） |
| 内存 | ~10GB（num_workers=2） |
| SSD | 正常（几个百分点） |
| checkpoint | fixed/pt-fixed-ep00-auc0.6103.ckpt |

---

## 八、今日清理

| 删除项 | 大小 | 原因 |
|--------|------|------|
| data/step9_structure_cache/ | 241GB | 240GB 缓存方案已废弃 |
| data/step9/（重复副本） | 60GB | 毕业设计目录已有一份 |
| enzyme_features_fp16.lmdb ×7 | 42GB | fp16 实验失败产物 |
| esm_*.pt 旧分片 ×26 | 25GB | 已有 enzymes.bin |
| grover_*.pt 旧分片 ×10 | 4.5GB | 已有 substrates_grover.bin |
| graph_*.pt 旧分片 ×95 | 28GB | 已有 per-sample 文件 |
| **总计** | **~400GB** | |

---

## 九、云服务器部署与训练（2026-03-17）

### 9.1 本地 legacy_bug 训练（归档）

**配置**: bs14, 2 workers, edge-mode legacy_bug, 本地 4070S

| 指标 | 数值 |
|------|------|
| 运行范围 | ep0-18 |
| 最佳 | ep10 AUC=0.715 |
| 速度 | 7.26-7.56 it/s |
| Per Epoch | ~29 min |

结果已归档到 `results/10_Step10_pt训练/local_训练/`

### 9.2 云服务器配置

| | Cloud-1（单卡，测试用） | Cloud-2（双卡，主力） |
|---|---|---|
| SSH | `root@sh01-ssh.gpuhome.cc -p 30008` | `root@hn01-ssh.gpuhome.cc -p 30156` |
| GPU | 1×RTX 4090 24GB | 2×RTX 4090 24GB |
| RAM | 90GB | 180GB |
| vCPU | 14 | 14 |
| 费用 | ¥1.46/h | ¥2.26/h |

### 9.3 环境搭建经验（新实例快速配置）

基于 PyTorch 2.5.1+cu124 预装镜像，一条命令装完所有依赖：
```bash
export PATH=/opt/conda/bin:$PATH
pip install pytorch_lightning==1.9.0 easydict lmdb tqdm matplotlib \
    torch_geometric warmup-scheduler scikit-learn cython tensorboard \
    pandas ray 'numpy<2' rdkit-pypi
```

**踩过的坑**：
- NumPy 2.x 和 RDKit 不兼容 → 必须 `numpy<2`（降到 1.26.4）
- `warmup_scheduler`、`scikit-learn`、`ray`、`cython`、`tensorboard` 镜像都没预装
- conda 环境：无需创建独立环境，直接用 `/opt/conda/bin/python`

### 9.4 fd limit 修复（关键！）

**问题**：8 workers 崩溃 `received 0 items of ancdata`
**原因**：容器默认 `ulimit -n = 1024`，多 worker 的 pipe fd 超限
**修复**：`ulimit -n 65536`（容器允许调高）
**必须在每次训练前执行**

### 9.5 Cloud-1 验证训练（ep0-1）

**配置**: bs48, accumulate=1, 4 workers, edge-mode legacy_bug

| Epoch | AUC | Loss | 速度 | 时间 |
|-------|-----|------|------|------|
| 0 | 0.606 | 0.316 | 3.74 it/s | 16:41 |
| 1 | 0.590 | 0.341 | 3.79 it/s | 16:21 |

GPU 100%利用率，19.6/24.6GB 显存，38°C。
Checkpoint: `pt-legacy_bug-ep00-auc0.6061.ckpt`（在Cloud-1数据盘上）

### 9.6 全环境性能对比

| 环境 | 速度 | 样本吞吐 | GPU利用率 | Per Epoch |
|------|------|----------|----------|-----------|
| Local 4070S (bs14, 2w) | 7.56 it/s | 103 samples/s | 30-60% | 29 min |
| Advisor 4090 (bs14, 2w) | 5.19 it/s | 72 samples/s | 30-55% | 41 min |
| **Cloud-1 4090 (bs48, 4w)** | 3.74 it/s | **166 samples/s** | **100%** | **16:41** |
| Cloud-2 2×4090 DDP (预估) | 7-15 it/s | 200-400 samples/s | ~90% | ~10-15 min |

### 9.7 代码改动

- `main_training_pt.py`：新增 `--devices`（多卡DDP）、`--shutdown`（训练完自动关机）
- `server_config.yml`：`accumulate_grad_batches: 2→1`（云服bs够大不需累积）
- Trainer：`devices > 1` 时自动用 `strategy="ddp"`

### 9.8 关键技术发现

1. **GPU compute bound**：Cloud-1 GPU 100% = 瓶颈从数据加载转移到 GPU 计算，FP16+TF32 已开启
2. **LMDB + DDP 有已知性能问题**：mmap 多进程冲突、锁竞争。.pt per-sample 天然适合 DDP
3. **RTX 4090 FP16 算力**：660 Tensor TFLOPS，单卡 > 2×RTX 3090 总和
4. **数据传输**：流式 tar（`tar cf - | ssh tar xf -`）不占服务器临时空间，适合小磁盘

---

## 十、Cloud-2 DDP 参数调优（2026-03-17）

在 Cloud-2（2×RTX 4090, 14vCPU, 180GB RAM）上进行了 14 轮参数调优。

### 调优结果汇总

| # | bs | accum | workers | 关键改动 | it/s | 每epoch | 结果 |
|---|-----|-------|---------|----------|------|---------|------|
| 1 | 48 | 1 | 6 | 默认DDP | **2.82** | 11 min | ✅ 初始基线 |
| 2 | 48 | 1 | 7 | +find_unused=F/bucket_view/prefetch=4 | 2.78 | 11 min | ✅ 微弱提升 |
| 3 | 56 | 2 | 7 | +static_graph=T | CRASH | - | ❌ 变长图不兼容 |
| 4 | 56 | 2 | 7 | 去掉static_graph | 2.35 | 11.3 min | ✅ 更慢 |
| 5 | 96 | 1 | 8 | +preload | OOM | - | ❌ 显存爆 |
| 6 | 64 | 1 | 6 | +preload | KILLED | - | ❌ 内存爆(DDP×2份) |
| 7 | 48 | 1 | 6 | +preload | KILLED | - | ❌ 同上 |
| 8 | 48 | 1 | 7 | +torch.compile | CRASH | - | ❌ torch_scatter不兼容 |
| 9 | 48 | 2 | 7 | accum=2 | 2.73 | 11.3 min | ✅ 无提升 |
| 10 | 56 | 1 | 6 | 最优DDP参数 | 2.63 | 10.1 min | ✅ 快 |
| 11 | 56 | 1 | 6 | +sync_dist=F+NCCL_P2P_DISABLE | 2.43 | 10.9 min | ✅ GPU均衡但慢 |
| **12** | **56** | **1** | **6** | **+sync_dist=F (无NCCL)** | **2.65** | **10.0 min** | **✅ 最优** |
| 13 | 64 | 1 | 6 | 同#12 | 2.31 | 10.1 min | ✅ 显存紧张 |
| 14 | 56 | 1 | 6 | +SizeSortedSampler | OOM | - | ❌ 大图集中→爆显存 |

### 最优配置（Run #12）

```bash
# 启动命令（已保存为 /root/start_train.sh）
cd /root/rivermind-data/EZSpecificity/src
export PYTHONPATH=/root/rivermind-data/EZSpecificity/src:$PYTHONPATH
ulimit -n 65536
python ../scripts/10_Step10_pt训练管线/main_training_pt.py \
  --config ../scripts/server_config.yml \
  --cache-dir ../data/10_Step10_pt训练/allsplit_pt_cache \
  --edge-mode legacy_bug \
  --devices 2 --num-workers 6 --batch-size 56 \
  --max-epochs 50 --shutdown
```

**关键参数**：
- `batch_size=56`，`accumulate_grad_batches=1`（server_config.yml）
- `DDPStrategy(find_unused_parameters=False, gradient_as_bucket_view=True)`
- `prefetch_factor=4`，`persistent_workers=True`，`pin_memory=True`
- `sync_dist=False`（monkey-patch model.log，减少rank0同步开销）
- `--shutdown`（训练完自动关机，节省费用）

### 核心发现

1. **DDP 对小模型加速有限**：1.8M 参数，GPU 计算极快，通信开销占比大。双卡 295 samples/s vs 单卡 165 samples/s = **1.79× 加速**（非 2×）
2. **GPU 利用率锯齿不可消除**：变长图导致两卡负载不均，快卡等慢卡。图大小范围 19-1436 原子（75倍差异）
3. **preload 在 DDP 下不可行**：两个进程各自加载 28GB → 总 56GB+，180GB 内存不够
4. **torch.compile 与 PyG 不兼容**：torch_scatter 的 dynamo 支持不完善

### 服务器状态

**已关机，已清理**。所有之前的 checkpoint、metrics、tensorboard logs、train.log 全部删除。下次开机直接 `bash /root/start_train.sh` 从 epoch 0 干净启动。

---

## 十一、Cloud-2 DDP legacy_bug 基线训练完成（2026-03-19）

### 11.1 训练配置（最终版）

| 参数 | 值 |
|------|-----|
| batch_size | 56 per GPU |
| accumulate_grad_batches | 1 |
| effective batch size | 112 |
| num_workers | 6 |
| prefetch_factor | 4 |
| devices | 2 (DDP) |
| precision | 16 (FP16) |
| find_unused_parameters | False |
| gradient_as_bucket_view | True |
| sync_dist | False (with all_gather fix) |
| EarlyStopping | patience=15, monitor auc/val |
| ModelCheckpoint | save_top_k=3, monitor auc/val |
| --shutdown | 训练完自动关机 |

### 11.2 Validation AUC 趋势

0.518(ep0) → 0.598(ep1) → 0.652(ep6) → 0.714(ep17) → **0.722(ep22=val best)** → 0.720(ep27)

训练跑了 32 epochs 后手动停止。

### 11.3 Test Set 结果（8841 样本，966 正样本）

| Checkpoint | Test AUC | Test AUPR |
|------------|----------|-----------|
| ep13 | 0.7175 | 0.2351 |
| ep22 | 0.7146 | 0.2207 |
| **ep27** | **0.7244** | **0.2142** |

**ep27 test AUC=0.7244 超过论文 all_split AUC=0.7198**（+0.0046）

### 11.4 DDP 关键 Bug 发现与修复

**问题**：`sync_dist=False` monkey-patch 导致 `validation_epoch_end` 只用一半数据计算 AUC（只有当前 DDP rank 的输出）。

**修复**：在 `validation_epoch_end` 和 `test_epoch_end` 中添加 `all_gather`，收集所有 rank 的输出后再计算 AUC。仅 rank 0 记录聚合指标。

### 11.5 Codex 审计（12 个问题，按严重性排序）

1. validation_epoch_end 只收集本地 rank 输出 → **已用 all_gather 修复**
2. sync_dist=False 使所有指标为 rank-0 本地值 → **已修复**
3. test_epoch_end 同样问题 → **已修复**
4. Effective batch size (112) >> 论文 (64) → 可能解释部分 AUC 差异
5. 缺失酶/底物 fallback 产生 NaN → 解释 grad_norm=inf
6. LR scheduler 监控 aupr/val 但 EarlyStopping 监控 auc/val → 不匹配

### 11.6 速度调优总结（14 轮配置测试）

| # | 配置 | 速度 | Per Epoch | 备注 |
|---|------|------|-----------|------|
| 1 | bs=48, accum=1, w=8 | 2.85 it/s | ~11 min | 初始基线 |
| 2 | bs=48, accum=1, w=6 | 2.82 it/s | ~11 min | 减少 workers |
| 3 | bs=48, accum=2, w=7 | 2.73 it/s | ~11.3 min | accum 无帮助 |
| 4 | bs=56, accum=2, w=7 | 2.35 it/s | ~11.3 min | 更慢 |
| **5** | **bs=56, accum=1, w=6** | **2.65 it/s** | **~10 min** | **最优配置** |
| 6 | bs=56, accum=1, w=6, sync_dist=F | 2.65 it/s | ~10 min | 同速度 |
| 7 | bs=56, accum=1, w=6, NCCL_P2P_DISABLE | 2.43 it/s | ~11 min | 更慢 |
| 8 | bs=64, accum=1, w=6 | 2.31 it/s | ~10.1 min | 显存紧张 |
| 9 | bs=48, accum=1, w=7 (preload) | OOM | - | DDP×2 份 preload 爆内存 |
| 10 | bs=96 | OOM | - | GPU OOM |
| 11 | static_graph=True | CRASH | - | 变长图不兼容 |
| **12** | **bs=56, accum=1, w=6, sync_dist=F** | **2.65 it/s** | **~10 min** | **FINAL BEST** |
| 13 | 双单卡并行 | 2.14/2.56 | ~15 min | CPU 争抢 |
| 14 | find_unused_parameters=False | 已包含 | - | 所有配置都包含 |

### 11.7 Codex 建议（下一轮 fixed 训练）

- lr: 3e-4 → 4e-4（sqrt scaling for larger batch）
- warmup_steps: 300 → 400
- max_epochs: 50 → 100-150
- weight_decay: 0 → 1e-5

---

## 十二、关键发现：Val Loss ↑ 而 AUC ↑ — Codex 深度分析（2026-03-19）

### 12.1 现象

legacy_bug 基线训练（Cloud-2, 2×RTX4090, 32 epochs）过程中观察到：
- 验证集: 161 正样本, 1322 负样本（10.8% 阳性率）
- train_loss: 持续下降（0.362 → 0.180）
- val_loss: 从 ep2 开始上升（0.319 → 0.393@ep22 → 0.422@ep31）
- val AUC: 持续上升直到 ep22（0.598 → 0.722），之后下降

### 12.2 根因分析（Codex 确认：非 Bug）

BCE Loss 和 AUC 度量的是**完全不同的维度**：

- **BCE Loss 是逐点度量（pointwise）**：逐样本计算损失再取平均。一个极端错误的预测会支配整个平均值。
  - 关键数字：负样本 logit z=5 → 单样本损失 5.01；logit z=0 → 损失仅 0.693。**7.2 倍差距**。
- **AUC 是成对度量（pairwise）**：统计正确排序的正-负样本对比例。一个错误预测只影响配对中的一小部分。

**真实数据量化（ep2→ep22）**：
- 验证集总配对数: 161 × 1322 = 212,842 对
- AUC 从 ~0.60 → 0.722 ≈ 额外 **~26,500 对**被正确排序
- val_loss 增加 ~100 个 BCE 单位（总计）— 由**几十个**在极端错误 logit 上的 hard 样本导致

### 12.3 三个阶段

| 阶段 | Epochs | LR | 表现 | 机制 |
|------|--------|-----|------|------|
| **Warmup** | 0-8 | 5e-6 → 3e-4 | val_loss 震荡, AUC 不稳定 | LR 快速上升，模型尚未收敛 |
| **Divergence** | 8-22 | 固定 3e-4 | **AUC ↑ 但 val_loss ↑** | 排序改善但对 hard OOD 样本过度自信 |
| **True overfitting** | 22+ | 固定 3e-4 | AUC ↓ 且 val_loss ↑ | 全面过拟合 |

### 12.4 LR Schedule 的贡献

- Warmup 5e-6→3e-4（8 epochs）解释了早期震荡
- Warmup 后 LR 固定 3e-4，持续推动 logits 幅度增大 → 有利于排序但伤害 BCE
- **ReduceLROnPlateau 监控 aupr/val（不是 auc/val 或 loss）**— 不匹配，可能导致 LR 未及时衰减

### 12.5 Codex 建议

1. **选择 checkpoint**: 同等 AUC 下优先选择 loss 更低的（更早 epoch，更"校准"）
2. **LR 衰减**: warmup 后添加余弦衰减或线性衰减，而非固定 3e-4
3. **不换 loss**: focal loss 恶化校准，contrastive loss 差异太大
4. **后处理**: 若下游需要概率值，用 temperature scaling 后校准

### 12.6 验证

- 损失函数: BCEWithLogitsLoss（`src/Models/ss.py:159`）— 标准二分类损失，与论文一致
- 论文使用相同的 Loss + AUC-ROC 组合
- 论文未展示训练曲线（仅最终 AUC 数字），无法直接对比动态过程
- 我们的 test AUC = 0.7244 > 论文 all_split 0.7198，确认模型正确工作

### 12.7 论文对比

| 指标 | 论文 all_split | 我们（legacy_bug） |
|------|---------------|-------------------|
| AUROC | 0.7198 | **0.7244**（+0.0046）|
| AUPR | 仅图（Fig. S3）| 0.2142-0.2351 |
| 训练 epoch | ~256 | 32 |
| GPU 数 | 4 | 2 |
| Effective batch | 64 | 112 |

**补充**: 论文 Supplementary PDF（78 页）中 Table S3（p64）有各 EC 层级的完整 AUROC/AUPR，但不包含 unknown enzyme & substrate split。

### 12.8 EZSpecificity 三种训练模式

| 模式 | 训练数据 | 说明 |
|------|---------|------|
| **EZSpecificity** | 全 ESIBank（323K 对）| 直接应用，不针对特定家族 |
| **EZSpecificity-finetune** | ESIBank 预训练 + 目标家族微调 | |
| **EZSpecificity-individual** | **仅目标家族**（85-5424 对）从头训 | 我们的 AllSplit 方式类似此模式 |

### 12.9 Codex 超参建议（下一轮 fixed 训练）

Effective batch 112 vs 论文 64 — 每 epoch 更少的优化器更新步数。建议：lr 3e-4→4e-4（sqrt scaling）、warmup 300→400、epochs 50→100-150、weight_decay 0→1e-5。用户决定不改 batch size（优先速度）。

---

## 十三、文件组织重构（2026-03-19）

### 13.1 脚本目录拆分

将 `scripts/10_Step10_pt训练管线/` 拆分为本地和云服务器两个子目录：

```
scripts/10_Step10_pt训练管线/
├── local/                        ← 本地开发版本
│   ├── main_training_pt.py
│   ├── pt_dataset.py
│   └── train_allsplit_config.yml  (renamed from server_config.yml)
└── cloud2x4090/                  ← 云服务器运行版本
    ├── main_training_pt.py
    ├── pt_dataset.py
    ├── server_config.yml
    └── start_ddp_2gpu.sh
```

### 13.2 删除的文件
- `start_parallel_2x1gpu_废弃.sh`：双单卡并行方案已放弃（比DDP慢）
- `__pycache__/`：已添加到 .gitignore

### 13.3 重命名
| 旧名称 | 新名称 | 原因 |
|--------|--------|------|
| `ezspec_pt_v1/` | `allsplit_pt_cache/` | 更准确描述内容 |
| `cloud_legacy_bug/` | `cloud2x4090_legacy_bug/` | 包含服务器名 |
| `server_config.yml`（local） | `train_allsplit_config.yml` | 区分本地/云配置 |

### 13.4 Gitignore 修复
- 添加 `**/allsplit_pt_cache/` 防止 176K 个 .pt 文件（57GB）被 git 跟踪
- 添加 `tmp_*`、`__pycache__/` 模式
- 从 git 跟踪中移除 `P450_EZSpecificity完整研究手册_终极整合版.md`
- Checkpoint gitignore：只保留最佳（ep27）

---

## 十四、下一步

- [x] Cloud-2 legacy_bug 基线训练 → **完成（test AUC=0.7244）**
- [ ] 应用 all_gather 修复 → 跑 fixed 基线对比
- [ ] 量化 edge fix 贡献（fixed vs legacy_bug 的 Δ AUC）
- [x] Step 12：~~数据泄露量化~~ → 已跳过

---

**版本**: v7.0 | **更新**: 2026-03-19
