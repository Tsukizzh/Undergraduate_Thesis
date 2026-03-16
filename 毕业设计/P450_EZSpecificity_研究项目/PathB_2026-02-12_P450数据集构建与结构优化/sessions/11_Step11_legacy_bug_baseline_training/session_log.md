# Step 11 Session Log: Legacy_bug Baseline Training + Server Optimization

**Date**: 2026-03-16 evening ~ 2026-03-17 02:00
**Status**: 🔄 IN PROGRESS (legacy_bug done, fixed pending)
**Key Achievement**: Local training shows 7.56 it/s stable (vs LMDB 2.09 it/s decay)

---

## 背景与目标

### 问题
- Step 10 完成了 .pt 缓存v3，但 ep14 基线包含了"边排序修复"(monkey-patch)
- 需要创建 "legacy_bug 基线版本"，来隔离边修复的贡献
- 同时验证不同硬件上的训练速度

### 目标
1. 运行 legacy_bug 版本从头训（与 fixed ep14 对标）
2. 量化边排序修复对 AUC 的贡献 (expect ±0.01-0.05)
3. 优化不同硬件的训练参数（本地/导师服务器/云服务器）

---

## 执行内容

### 1. 本地 Legacy_bug 训练（4070S GPU）✅

**配置**:
```bash
cd D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化
python ../../scripts/10_Step10_pt训练管线/main_training_pt.py \
  --config ../../scripts/server_config.yml \
  --cache-dir data/10_Step10_pt训练/ezspec_pt_v1 \
  --edge-mode legacy_bug \
  --num-workers 2
```

**结果**:
| 指标 | 数值 |
|------|------|
| 运行范围 | ep0-18 |
| 最佳 | ep10 AUC=0.715 |
| 速度 | 7.26-7.56 it/s (稳定，无衰减) |
| Per Epoch | ~29 min |
| RAM | ~10GB (稳定) |
| 状态 | 暂停，可用 --resume last 续跑 |

**关键发现**:
- 速度与 ep14 fixed 版本相同 (7.56 it/s vs LMDB 的 2.09 it/s)
- 无论 --edge-mode legacy_bug 还是 fixed，基础速度不变
- 修复仅影响数据加载逻辑，不影响 I/O 性能

---

### 2. 导师服务器参数调优（jqlab4090 RTX 4090）✅

**测试矩阵** (16GB 单通道 RAM, i7-10700 CPU):

| Config | Speed | iters/epoch | Per Epoch | RAM | Swap |
|--------|-------|-------------|-----------|-----|------|
| **bs14, 2w** | **5.19 it/s** | 12738 | **41 min** | 9.2GB | 6.6GB |
| bs24, 2w | 3.07 it/s | 7431 | 40 min | 10.3GB | 7.8GB |
| bs24, 1w | 2.30 it/s | 7431 | 54 min | 8.9GB | 4.6GB |
| bs32, 4w | 1.39 it/s | 5574 | 67 min | 13.8GB | 13.5GB |
| bs32, 2w | 2.31 it/s | 5574 | 40 min | 10.7GB | 8.4GB |

**核心发现**:
- **最优配置**: bs14 + 2 workers = 5.19 it/s, 41 min/epoch
- **瓶颈确认**: CPU (i7-10700 单线程) + 16GB 单通道RAM，非 GPU
- **GPU 利用率**: 30-55% (远低于本地 95%)
- **40 min/epoch 是硬极限** (12738 iter × 单线程速率 ~3000 iter/min)
- **增加 workers 无益**: 4w 反而变慢 (CPU过载)

**.pt 方案的价值**:
- 不是绝对速度差异，而是**避免 LMDB mmap 污染**
- 16GB < 60GB LMDB → mmap thrashing(pageout/pagein 循环)
- .pt per-sample = 160KB/sample × 并发 = ~10GB 工作集
- 避免了无休止的页面置换，保持稳定性

**结论**: 导师服务器适合 backup 训练，不是主力。用.pt缓存绕过 mmap 污染是关键。

---

### 3. 云服务器租赁（智川云 Rivermind）✅

**服务规格** (新增):

| 指标 | 配置 |
|------|------|
| 平台 | 智川云 (https://gpuhome.cc) |
| GPU | RTX 4090 24GB |
| vCPU | **14 核心** (vs i7-10700 的 8 核) |
| RAM | **90GB** (关键) |
| SSD | 75GB (50 free + 25 paid) |
| 网络 | 1Gbps |
| CUDA | 12.4, PyTorch 2.5.1 预装 |
| 成本 | ¥1.46/h + ¥0.25/day storage |
| SSH | root@sh01-ssh.gpuhome.cc -p 30006 |

**期望性能** (基于 RAM + CPU规格):
- 14 vCPU vs i7-10700 的 8 核 → 应该能充分喂饱 GPU
- 90GB RAM vs 16GB → 完全避免 swap，避免 mmap 污染
- 预期速度: **7-15 it/s** (should be GPU-bound like local)
- 若达到 7.56 it/s，per epoch = 29 min，20 epochs = 570 min ≈ 10h = ¥14.6

**状态**:
- ✅ 已创建(无 GPU 模式，准备环境)
- 📋 待做: conda 激活, .pt 数据传输, 速度基准测试

---

### 4. 五种硬件配置完整对比

**所有测试汇总** (2026-03-16 之前 + 新增):

| 配置 | 硬件 | 速度 | Per Epoch | 内存 | 用途 |
|------|------|------|-----------|------|------|
| LMDB bs4/a8 | Local 4070S | 1.71 it/s | 2h | 15GB | baseline |
| LMDB bs8/a4/6w | Local 4070S | 4.3 it/s | 49min | 20GB | 优化无缓存 |
| LMDB cache | Local 4070S | 3.8→2.09 衰减 | 47-100min | 30GB+ | Step 9 |
| **.pt v3 2w** | **Local 4070S** | **7.56 it/s** | **29 min** | **10GB** | **Step 10/11** |
| .pt v3 2w | Advisor 4090 | 5.19 it/s | 41 min | 9.2+6.6 | backup |
| (预期) | Cloud 4090 | 7-15 it/s | 15-25 min | ~15GB | 优先 |

---

## 计算环境角色定位（最终确定）

### Local (Windows, 4070S 12GB VRAM, 32GB RAM)
- **用途**: 代码开发、测试、.pt 构建、小规模训练(<20 epochs)
- **优势**: 交互快，迭代快，本地调试方便
- **限制**: 显存 12GB (batch_size=14 天花板), RAM 32GB 与系统共享
- **性能**: 7.56 it/s, 29 min/epoch

### Advisor 服务器 (jqlab4090, 4090 24GB, i7-10700 8核, 16GB RAM)
- **用途**: 备用训练（不是主力）
- **优势**: 硬件质量好，GPU 驱动稳定，局部化（无网络延迟）
- **限制**: 16GB 单通道 RAM 极其有限，CPU 低端（无法喂饱 GPU）
- **性能**: 5.19 it/s, 41 min/epoch (最优配置)
- **状态**: 有 legacy_bug 训练在后台(nohup)，可以备选

### Cloud 云服务器 (智川云, 4090 24GB, 14核, 90GB RAM)
- **用途**: 优先训练（支撑 Phase 2-7 主要实验）
- **优势**: 大内存避免 swap，足够 CPU 喂 GPU，按小时计费便宜
- **限制**: 网络延迟，初始数据传输成本
- **性能**: 预期 7-15 it/s (GPU-bound), 15-25 min/epoch
- **成本**: ¥1.46/h, 20 epochs ≈ ¥15-20
- **状态**: 准备中，待速度基准

---

## Step 12 状态：已跳过 ✅

**用户决定** (2026-03-16):
- 底物 6.56% 泄露已完整审计，不单独量化影响
- 跳过 Phase 3a 三轮训练（clean/leaked/all）
- **节省**: 3-5 天训练时间

---

## 7 阶段计划时间表（修订）

| Phase | 内容 | 完成期 | 时长 | 状态 |
|-------|------|--------|------|------|
| **1** | **.pt 迁移到云 + 速度基准** | **2026-03-20** | 3天 | 🔄 准备中 |
| **2** | **legacy_bug vs fixed 对标** | **2026-03-23** | 3天 | ⏳ 当前 |
| 3 | ~~数据泄露量化~~ **跳过** | — | — | ✅ |
| 4 | Fe+Heme 扩展实验 | 2026-03-30 | 3天 | ⏳ |
| 5 | 通道融合修复 | 2026-04-02 | 3天 | ⏳ |
| 6 | P450 数据增强 | 2026-04-05 | 3-5天 | ⏳ |
| 7 | 最终多轮评估+论文结果 | 2026-04-10前 | 2天 | ⏳ |

**关键里程碑**:
- 2026-03-20: 云服基准速度确定 (≥7 it/s 则 GPU bound，计划可行)
- 2026-03-25: edge fix 贡献量化 (预期 Δ AUC ±0.01-0.05)
- 2026-04-10: 论文级最终结果交付

---

## 下一步（Phase 1, 2026-03-17~20）

### 立即任务
1. **云服务器环境准备**:
   - SSH 连接验证
   - conda 激活 ezspecificity 环境
   - pip 依赖确认

2. **.pt 数据迁移**:
   - 打包本地 data/10_Step10_pt训练/ezspec_pt_v1/ (57GB)
   - 传输到云服 (初估 2-3 小时@100Mbps)
   - 验证数据完整性

3. **云服速度基准测试**:
   - 运行 1 epoch with legacy_bug + fixed
   - 目标: ≥7 it/s (表示 GPU bound，plan viable)
   - 若 <5 it/s: 网络瓶颈，需优化传输

4. **本地 fixed 版本启动**:
   - 启动 --edge-mode fixed 版本基线训
   - 与 legacy_bug (best=ep10) 对标

### 预期完成日期
- 2026-03-20: Phase 1 完成，云服基准定下
- 2026-03-25: Phase 2 完成，edge fix 贡献量化

---

## 文件与产出

### 新增脚本
- `scripts/10_Step10_pt训练管线/main_training_pt.py` — .pt 训练入口 (已支持 --edge-mode legacy_bug/fixed)
- `scripts/10_Step10_pt训练管线/pt_dataset.py` — Dataset 类

### 新增配置
- `scripts/server_config.yml` — 云服务器配置（待编写）

### 数据
- `data/10_Step10_pt训练/ezspec_pt_v1/` — .pt 缓存 v3 (57GB, 完成)

### 会话日志
- 本文件: `sessions/11_Step11_legacy_bug_baseline_training/session_log.md`

---

## 关键发现与经验教训 💡

### 1. LMDB vs .pt 性能真相
- **非绝对速度差异**，而是**内存污染避免**
- LMDB 60GB mmap > 32GB RAM → 页面缓存颠簸(thrashing) → 2.09 it/s
- .pt per-sample 160KB/iter → ~10GB 工作集 → 稳定 7.56 it/s
- **结论**: 大数据场景，小文件随机读 > 单个大文件 mmap

### 2. 硬件选择关键点
- **本地**: GPU 充足，RAM 充足，但显存限制(batch_size)
- **导师服务器**: GPU 好但 CPU 弱 + RAM 少 → CPU bottleneck (GPU idle 50%)
- **云服务器**: 关键是**足够 RAM + 足够 CPU 喂 GPU** → GPU bound (95%+)
- **CPU 不是余数**: 14 vCPU can feed GPU much better than i7-10700 8核

### 3. 成本效益分析
- 本地 free (电费+硬件折旧)，per epoch 30 min
- 导师服务器 free，per epoch 41 min (backup ok)
- 云服务器 ¥1.46/h ≈ ¥30/20epoch，per epoch 15-25 min (3x faster!)
- **对于 7 阶段计划(20-40 epochs 实验)**: 云服总成本 ¥100-200，节省时间 5-10 天 → 值得

### 4. Edge Fix 的预期贡献
- 当前假设: Δ AUC ±0.01-0.05
- 若 <0.01: 可放后，非关键优化
- 若 >0.05: 非常重要，需深入分析原因

---

**Session 记录人**: Claude Code
**最后更新**: 2026-03-17 02:00
