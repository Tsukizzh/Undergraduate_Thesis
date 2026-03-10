# Step 9: Unknown Enzyme+Substrate (all_split) 模型训练

## 目标
从头训练 EZSpecificity 模型，使用 all_split（unknown enzyme+substrate）分割方式，这是论文中最难的场景。

## 时间线

| 时间 | 事件 | 备注 |
|------|------|------|
| 2026-03-10 开始 | Step 9 启动 | 目标：先跑 1 epoch 看时间 |
| 2026-03-10 03:48 | 冒烟测试通过 | config/CSV/model init 全部正确 |
| 2026-03-10 04:15 | Pipeline测试启动 | 20 train + 5 val batches (G: Drive) |
| 2026-03-10 ~04:30 | Pipeline测试完成 | checkpoint saved, loss ~0.677 |
| 2026-03-10 04:43 | 等待数据拷贝 + 自动启动训练 | enzyme 56% (28/50GB) |

## 硬件配置
- GPU: NVIDIA 4070 Super, 12GB VRAM
- 磁盘: D: 盘 874GB 可用
- OS: Windows 11
- Python: 3.10.12, CUDA 12.1, PyTorch 2.1.0

## 训练配置
```yaml
split: all_split (unknown enzyme+substrate)
fold: 0
dataset: BRENDA only
batch_size: 4
accumulate_grad_batches: 8
effective_batch: 32
num_workers: 2  # with LMDB close + lazy reconnect fix
persistent_workers: True
precision: fp16
optimizer: AdamW
lr: 3e-4
warmup_epochs: 8
```

## 数据来源
- Google Drive: `G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\`
- 分割: `brenda/all_split/` fold 0
- 训练集: 340,787 样本（有效: ~143,355, 42%）
- 验证集: 2,876 样本（有效: ~1,176, 41%）
- 测试集: 16,316 样本

## 数据文件大小
| 文件 | 大小 | 状态 |
|------|------|------|
| enzyme_features.lmdb | 50.2GB | 拷贝中 |
| grover_fingerprint.lmdb | 8.9GB | ✅ 完成 |
| reaction_features.lmdb | 197MB | ✅ 完成 |
| morgan_fingerprint.npy | 310MB | ✅ 完成 |
| structure/structure_features.lmdb | ~GB | ✅ 完成 |
| structure/sequence_features.lmdb | ~GB | ✅ 完成 |

## 关键代码修复

### 1. CSV 列名映射
Google Drive CSVs 使用小写列名 (`enzyme, reaction, label, structure_index`)，代码期望大写 (`Enzyme Index, Substrate Index, Label, Dock Index`)。解决方案：预处理 CSVs 并 rename 列名，保存为 `*_renamed.csv`。

### 2. LMDB map_size
代码硬编码 `map_size=10GB`，但 enzyme_features.lmdb 有 51GB。解决方案：monkey-patch `lmdb.open()`，对 readonly 打开使用 128GB map_size。

### 3. Windows 多进程 + LMDB (num_workers > 0)
Windows spawn 模式无法 pickle LMDB handles。解决方案：
- 在 `valid_idx` 构建完成后关闭所有 LMDB handles（设为 None）
- Worker 进程通过 `_connect_db()` 延迟重连
- `persistent_workers=True` 避免重复创建/销毁 worker

### 4. 其他修复
- `high_quality_id_path`: 使用字符串 `'None'` 而非 YAML null
- `sample_weight`: `[1.0, 1.0]`（complex_weight, ligand_weight）
- `precision=16`（PL 1.9.0 不支持 `'16-mixed'`）

## Codex 审核记录

### 审核 #1: 可行性评估 (2026-03-10)
- **结论**: 可行，但有多个坑需要填
- **关键风险**:
  1. CSV 列名不匹配（小写 vs 大写）→ 已修复
  2. LMDB map_size 硬编码 10GB → 已修复（monkey-patch）
  3. batch_size=16 可能炸显存 → 使用 4
  4. 必须从头训练，不能用 random_split 检查点
  5. Windows num_workers 问题 → 已修复（LMDB close + lazy reconnect）
  6. 训练时间估计: 2-7 天/折

### 审核 #2: 训练脚本代码审核 (2026-03-10)
- **修复**: sample_weight 从 `[1.0]` 改为 `[1.0, 1.0]`
- **修复**: high_quality_id_path 从 YAML null 改为字符串 `'None'`

### 审核 #3: num_workers 优化方案 (2026-03-10)
- SESSION_ID: `019cd44e-09e7-7a60-855d-f45e1eba7d0a`
- **结论**: LMDB close + lazy reconnect 方案可行
- **建议**:
  1. batch_size=4 更安全（模型有 1450 长度蛋白+cross-attention）
  2. num_workers=2 + persistent_workers=True 是好的起点
  3. morgan_dbs (numpy) 保留不关闭（可 pickle，避免重复加载）
  4. 建议未来: 拆分 `_connect_db()` 为 "open" 和 "build_key_dict"
  5. 建议未来: 添加 `__getstate__()` 自动剥离 LMDB handles

## 文件结构
```
PathB_.../
├── scripts/09_Step9_AllSplit训練/
│   ├── train_allsplit_config.yml    # 训练配置
│   ├── main_training.py             # 主训练脚本
│   ├── smoke_test.py                # 冒烟测试
│   ├── pipeline_test.py             # Pipeline测试 (20 batch)
│   └── wait_and_train.py            # 等待数据+自动启动
├── data/09_Step9_AllSplit训练/brenda/      # 本地数据副本
│   ├── all_split/                   # CSV splits (renamed)
│   ├── enzyme_features.lmdb         # 蛋白特征 (50GB)
│   ├── grover_fingerprint.lmdb      # GROVER指纹 (9GB)
│   ├── reaction_features.lmdb       # 反应特征 (197MB)
│   ├── morgan_fingerprint.npy       # Morgan指纹 (310MB)
│   └── structure/                   # 结构特征
├── logs/09_Step9_AllSplit训练/      # TensorBoard日志
├── results/09_Step9_AllSplit训练/    # checkpoints + 训练日志
└── sessions/09_Step9_.../
    └── session_log.md               # 本文件
```

## 进度
- [x] 数据拷贝到本地（enzyme 拷贝中，其余完成）
- [x] 代码分析 + 修复（CSV、LMDB、num_workers）
- [x] 训练脚本编写
- [x] 冒烟测试 + Pipeline测试
- [ ] 1 epoch 训练计时
- [ ] 结果分析
