# Step 10 .pt 预处理缓存与训练管线 — Session日志

> **日期**: 2026-03-16
> **目标**: 构建高效.pt预处理缓存解决LMDB mmap污染
> **成果**: ✅ .pt v3方案 7.56 it/s 完成

---

## 一、性能优化三迭代

### v1 Shard模式 ❌
- 87个300MB shard = 15TB/epoch读取
- SSD饱和100% → 不可用

### v2 Flatbin模式 ⚠️  
- enzyme/substrate加速，但graph仍为300MB
- 速度 0.38-0.52 it/s → 瓶颈未解

### v3 Per-sample模式 ✅
- 176K个160KB .pt文件
- 速度 7.56 it/s（稳定，无衰减）
- **per epoch: 27min** (vs 100min LMDB衰减版)

---

## 二、存储布局（44GB）

```
ezspec_pt_v1/
├── train/samples/  26GB (176,843 per-sample .pt)
├── val/samples/    0.1GB (1,482)
├── test/samples/   0.7GB (8,841)
├── enzymes/        17GB (26 fp16 shards + index)
└── substrates/     0.5GB (10 GROVER shards + index)
```

Per-sample .pt包含：
- enzyme_id, substrate_id
- graph (nodes, edges, 3D coords, bond types, labels)

---

## 三、构建完成统计

| 指标 | 值 |
|------|---|
| 总样本 | 186,166 |
| 构建时间 | 2分钟 (10 workers) |
| 验证 | ✅ 所有样本通过 |
| 总大小 | 44GB |

---

## 四、Dataset类（Step 10a）

**文件**: `scripts/10_Step10_pt训练管线/pt_dataset.py`

核心：
- 加载per-sample .pt图数据 (160KB)
- 从enzyme/substrate shard动态读特征
- Multi-worker safe
- Fallback to flatbin.bin

**Codex审核**: ✅ 3轮通过

---

## 五、性能对比

```
方案              速度        Per Epoch
LMDB (衰减后)    2.09 it/s    100+ min
.pt Per-sample   7.56 it/s    27 min ★

改进：+262% 速度，-73% 时间
```

---

## 六、服务器迁移计划

1. **本地**：tar.gz打包 44GB → 13GB (压缩率70%)
2. **上传**：→ dc@jqlab4090:~/zhuangzeheng/EZSpecificity/data/
3. **解压**：验证186,166文件
4. **测试**：Phase 1速度基准 (期望≥4 it/s in 16GB RAM)

**时间**: 
- 2026-03-17：打包+迁移
- 2026-03-18~19：服务器速度测试
- 2026-03-20+：Phase 2 baseline实验

---

## 七、Next Steps

- ✅ Step 10a Dataset完成
- 🔄 Step 10b 训练脚本编写中
- ⏳ Phase 1 服务器就绪 (2026-03-18)

---

**版本**: v1.0 | **更新**: 2026-03-16 22:45
