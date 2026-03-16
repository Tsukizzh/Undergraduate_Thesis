# Step 10: .pt 训练管线 — 完整 Session 日志

> **日期**: 2026-03-16
> **目标**: 用 .pt 预处理缓存替代 LMDB，消除 mmap 抖动和 CPU k-NN 瓶颈
> **最终成果**: ✅ per-sample .pt v3 方案，7.33-7.70 it/s，29 min/epoch，内存 ~10GB

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
ezspec_pt_v1/ (~57GB)
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

## 九、下一步

- [ ] 传 per-sample 文件到老师服务器（train 26GB + val + test + enzymes.bin + substrates）
- [ ] Step 11：跑两轮基线训练（legacy_bug + fixed）
- [ ] Step 12：数据泄露量化
- [ ] Step 13：消融实验（Dropout / Fe+Heme / 融合）

---

**版本**: v2.0 | **更新**: 2026-03-16 08:00
