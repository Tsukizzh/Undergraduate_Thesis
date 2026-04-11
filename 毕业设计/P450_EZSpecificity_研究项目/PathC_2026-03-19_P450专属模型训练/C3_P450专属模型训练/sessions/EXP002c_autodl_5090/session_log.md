# EXP002c: AutoDL 5090 调参实验

> 2026-04-11 服务器迁移至AutoDL 2×RTX5090 + 性能调优

## 服务器信息

- **平台**: AutoDL
- **GPU**: 2×RTX 5090 32GB
- **CPU**: 50 vCPU
- **RAM**: 180GB（cgroup硬限制）
- **数据盘**: 150GB SSD
- **PyTorch**: 2.8.0+cu128, PL 2.6.1, Python 3.12, CUDA 12.8
- **费用**: ¥5.57/时

## PL 2.6 迁移改动

从PL 1.9迁移到2.6需要的代码修改（用autodl_migration.py自动化脚本完成）：

| 文件 | 改动 |
|------|------|
| ss.py | validation_epoch_end→on_validation_epoch_end + 手动缓存outputs + optimizer_step 8→4参数 + __init__加output lists |
| main_training_pt.py | precision=16→"16-mixed" + DDP patch用all_gather_object + on_before_optimizer_step去掉opt_idx |
| transforms.py | 注释掉pyximport（Python 3.12不兼容） |
| cpi.py | 同ss.py的epoch_end改动 |
| run_train.sh | 路径rivermind→autodl-tmp + conda路径 + NCCL环境变量 |
| config.yml | root_dir路径 |

## 关键发现

### NCCL eager connect问题（致命）
- PyTorch 2.8新增`eager_connect_single_device`在init_process_group时预分配NCCL缓冲区
- 在5090上触发`Cuda failure 'out of memory'`
- **修复**: `export TORCH_NCCL_ENABLE_EAGER_CONNECT=0`
- 同时需要: `NCCL_P2P_DISABLE=1`, `NCCL_IB_DISABLE=1`, `NCCL_NVLS_ENABLE=0`

### 180GB cgroup硬限制
- AutoDL严格执行cgroup内存限制（gpuhome不严格执行）
- preload + DDP双rank + workers在180GB下空间很紧
- 多次OOM kill的根因都是内存超限

### static_graph不可用
- `DDPStrategy(static_graph=True)` 在PyTorch 2.8触发内部assert错误
- `expect_autograd_hooks_ INTERNAL ASSERT FAILED at reducer.cpp:1660`
- 是PyTorch bug，不是我们的代码问题

### 不兼容的优化
- **bf16-mixed**: BF16 tensor不能直接.numpy()，代码里有cpu().numpy()操作
- **fused AdamW**: 和PL的gradient_clip_val冲突
- **torch.compile**: 模型太小(1.85M)，compile overhead > 收益，反而从2.1降到1.16 it/s

## 性能调优实验记录

所有实验基于: devices=2, preload=✅, lr=4e-4, warmup=12, wd=1e-5, dropout=0.9

| # | bs | workers | accum | val freq | 其他 | 速度(it/s) | epoch时间 | 结果 |
|---|-----|---------|-------|----------|------|-----------|----------|------|
| 1 | 56 | 6 | 2 | 1 | 基线 | 1.84 | ~2.4min | ✅ GPU利用率波动大 |
| 2 | 56 | 12 | 2 | 1 | — | — | — | ❌ OOM killed |
| 3 | 72 | 6 | 1 | 1 | bf16+compile+fused | — | — | ❌ bf16 numpy错误 + fused冲突 |
| 4 | 72 | 6 | 1 | 1 | compile+static_graph | 1.16 | — | ❌ compile反而慢 + OOM |
| 5 | 80 | 12 | 1 | 1 | 无preload+static_graph | — | — | ❌ shm bus error |
| 6 | 80 | 8 | 1 | 1 | 无preload | — | — | ❌ OOM killed |
| 7 | 80 | 6 | 1 | 1 | 无preload | 1.11 | — | ❌ 验证阶段OOM |
| 8 | 56 | 6 | 2 | 1 | 确认基线 | 1.84 | ~2.0min | ✅ |
| **9** | **64** | **8** | **2** | **2** | **—** | **2.27** | **~1.6min** | **✅ 最快(val=2)** |
| 10 | 64 | 10 | 2 | 1 | — | 2.09 | ~2.0min | ✅ GPU 98%稳定 |
| 11 | 72 | 8 | 2 | 1 | — | 1.74 | ~2.1min | ✅ GPU1显存96%满 |
| 12 | 64 | 8 | 2 | 1 | val batch=128 | — | — | ❌ OOM killed |
| 13 | 64 | 8 | 2 | 1 | val batch=96 | — | — | ❌ OOM killed |

## 最优配置

**bs=64, workers=8, preload=✅, accumulate=2, val_frequency=1**
- effective batch = 64×2×2 = 256
- 训练速度: ~2.1-2.27 it/s
- 每epoch: ~2.0min（含验证）
- GPU利用率: 改善但仍有波动

## EXP002c当前代码状态

相比EXP002b（原始PL 1.9代码）的改动：
1. PL 2.6兼容性修改（ss.py, main_training_pt.py, transforms.py, cpi.py）
2. ss.py __init__增加validation_step_outputs/test_step_outputs初始化
3. NCCL环境变量（EAGER_CONNECT=0等）
4. 路径autodl-tmp + shutdown命令适配
5. val batch恢复为统一batch_size（val加倍方案OOM）
6. static_graph已去掉（PyTorch bug）

## 结论

180GB cgroup硬限制是性能天花板。5090的32GB显存和50核CPU无法被充分利用。
最优配置下~2.0min/epoch，vs 旧gpuhome 4×4090的~1.4min/epoch，实际更慢。
如果能拿到cgroup更宽松的服务器（或360GB RAM），性能会好很多。
