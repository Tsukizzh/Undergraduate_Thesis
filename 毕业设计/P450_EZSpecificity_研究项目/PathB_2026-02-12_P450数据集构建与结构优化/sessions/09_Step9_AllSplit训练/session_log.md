# Step 9: Unknown Enzyme+Substrate (all_split) 模型训练

## 目标
从头训练 EZSpecificity 模型，使用 all_split（unknown enzyme+substrate）分割方式，这是论文中最难的场景。

## 时间线

| 时间 | 事件 | 备注 |
|------|------|------|
| 2026-03-10 03:48 | 冒烟测试通过 | config/CSV/model init 全部正确 |
| 2026-03-10 04:15 | Pipeline测试启动 | 20 train + 5 val batches (G: Drive) |
| 2026-03-10 ~04:30 | Pipeline测试完成 | checkpoint saved, loss ~0.677 |
| 2026-03-10 04:43 | 等待数据拷贝 | enzyme 56% (28/50GB) |
| 2026-03-10 ~20:00 | 基线性能测量 | 5分钟监控: 1.71 it/s, GPU 8% |
| 2026-03-10 ~20:30 | 方法1参数调优实施 | batch=8, workers=6, prefetch=4 |
| 2026-03-10 21:20 | 方法1基准测试启动 | 44分钟运行 + GPU监控 |
| 2026-03-10 22:05 | 方法1基准测试完成 | 11488/22292 batches, 4.3 it/s |
| 2026-03-10 22:10 | Codex方法2辩论(3轮) | 结构缓存方案设计完成 |
| 2026-03-10 22:30 | num_workers=10尝试 | 32GB RAM不够，回退到6 |
| 2026-03-10 23:00 | 方法2验证咨询(4轮Codex) | 发现原始代码边排序Bug，确定修复方案 |
| 2026-03-11 ~01:00 | 方法2实施完成 | cache_utils.py + build_structure_cache.py + main_training_cached.py |
| 2026-03-11 ~02:00 | 冒烟测试全部通过 | 导入/端到端构建/运行时重建(fixed+legacy)/缓存往返 |
| 2026-03-11 ~02:30 | Codex审核#5-7(构建脚本) | resume支持 + 计数验证 + 性能建议 |
| 2026-03-11 ~03:00 | 缓存构建启动 | `--split all`, 176K样本, 预计2-4小时 |

## 硬件配置
- GPU: NVIDIA RTX 4070 Super, 12GB VRAM
- CPU: Intel i5-14600KF, 14 cores (6P+8E), 20 threads
- RAM: 32GB（6个DataLoader worker是极限）
- 磁盘: D: 盘 SSD, 874GB 可用
- OS: Windows 11
- Python: 3.10.12, CUDA 12.1, PyTorch 2.1.0

## 训练配置
```yaml
split: all_split (unknown enzyme+substrate)
fold: 0
dataset: BRENDA + 6 small families (halogenase excluded)
batch_size: 8          # 方法1优化: 4→8
accumulate_grad_batches: 4  # 方法1优化: 8→4 (effective batch仍=32)
effective_batch: 32
num_workers: 6         # 方法1优化: 2→6 (32GB RAM极限)
prefetch_factor: 4     # 方法1优化: 2→4
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

## 性能优化

### 基线性能（优化前）
- GPU利用率: 8% avg, max 40%
- 速度: 1.71 it/s (0.58s/batch)
- 44,582 batches/epoch, ~7.2h/epoch, 50 epochs ≈ 15天
- 瓶颈: CPU数据加载（EdgeConnection k-NN占50-70%）

### 方法1: 参数调优（已完成 + 基准测试）
| 指标 | 优化前 | 方法1 | 变化 |
|------|--------|-------|------|
| 速度 | 1.71 it/s | **4.3 it/s** | +2.5x |
| GPU利用率(avg) | 8% | **29.0%** | +3.6x |
| GPU利用率(max) | 40% | **94%** | 可满载 |
| VRAM | ~8GB | **11.4GB (93%)** | 接近极限 |
| 每epoch | ~7.2h | **~1.44h** | -5x |
| 50 epochs | ~15天 | **~3天** | 可行 |

GPU脉冲模式: 646次活跃段(avg 2.0s, ~60% util) + 646次空闲段(avg 2.1s)
瓶颈仍是CPU数据加载 — GPU每次~2s算完，等~2s等下一个batch

### 方法2: 结构缓存（设计完成 + 4轮Codex验证完成，待实施）

#### 初始设计（Codex 3轮辩论）
- SESSION_ID: `019cd7e8-8259-7b93-a2bf-362d9cce45d2`
- 仅缓存结构分支（k-NN图+原始距离），序列分支保持实时读取
- 不修改src/，新增3个文件
- 磁盘: ~90-140 GB, 构建时间: 2-4小时（一次性）
- 预期: 4.3→8-12 it/s, 50 epochs ≈ 1-1.5天

#### 4轮Codex验证咨询（SESSION_ID: 019cd830-4107-7621-9360-97676a1dcec0）

##### 第1轮: 边排序Bug发现与确认
- **重大发现**: `src/Datasets/Structure/transforms.py` 中 `EdgeConnection.__call__` 存在边排序不匹配Bug
  - `complex_edge_attr` (line 132/137/140/142) = `[knn边属性, 真实键边属性]` 顺序
  - `complex_edge_index` (line 147) = `[真实键边索引, knn边索引]` 顺序
  - **边i的属性不对应边i的索引！**
- `_get_dist` (line 109) 不改变edge_index顺序，确认无隐藏重排
- EGNN消息传递 (egnn.py:34-52) 将 `edge_attr[i]` 与 `edge_index[:, i]` 一一配对使用
- **结论**: 这是真实Bug，不是顺序不变的

##### 第2轮: Bug严重性评估与修复方案
- **Bug严重程度: 高**
  - 真实化学键（单键/双键/芳香键）丢失键类型身份，被分配knn的type-5属性
  - 真实键距离(~1-2Å)被替换为随机knn距离
  - 交叉边标志(cross_bond)被污染
  - 典型样本中 Ek >> Er（knn_graph k=48生成的边远多于化学键），大多数边都被错配
  - attention门控 (egnn.py:47-48) 使用了被污染的消息，进一步放大错误
- **模型为何仍能工作**: 通过抑制/下调结构分支来"补偿"bug，非图通道（ESM、GROVER、Morgan）提供主要信号
  - 这可能**解释了Step 8 E8v2中结构通道变化不显著的现象** — 结构通道本就在被错误数据训练
- **修复方案决策（由于从头训练）**:
  - **Option A（复现bug）**: 安全，与论文可比
  - **Option B（修复bug）**: 结构模块可正确学习，潜在更好 ← **推荐**
  - **Option C（两者都训）**: 最充分但费时
- **Codex建议**: 以Option B为主训练，保留legacy_bug模式作为可选基线
- **最佳缓存设计**: 分区存储real/knn边（不存打包后的final tensors），运行时按模式组装
  - 一次缓存构建，同时支持fixed和legacy_bug两种模式

##### 第3轮: 实现细节验证（6项）
1. **str_tag/ligand-only**: 当前配置 `full_data: False`，不会出现ligand-only样本。`str_tag`始终为`'complex'`。但代码应支持Ek=0以兼容未来配置
2. **字段名冲突**: 无冲突。结构侧字段（`protein_x`, `ligand_x`, `complex_edge_*`等）通过前缀机制与序列侧字段（`element`, `edge_index`, `embedding`等）完全隔离
3. **sample_weight/mask_use**: 来自`StructureDataset.getitem_with_real_idx`（非transform）。当前配置下均为常量（sample_weight=1.0，mask_use_complex=全1，mask_use_ligand=全0）。mask_use在`mode: complex`下未使用，但保留以兼容
4. **LMDB多进程安全**: 多进程并发只读安全。推荐设置:
   - writer: `map_size=256GB+`, `readonly=False`
   - reader: `map_size=256GB+`, `readonly=True`, `create=False`, `lock=False`（缓存不可变）
   - 每个spawn worker是独立进程，各自开LMDB环境，安全
5. **str_tag来源**: `StructureDataset.getitem_with_real_idx` (structure.py:324/332)设置，`EdgeConnection`仅消费。模型也使用(ss.py:154)，必须缓存
6. **y vs label**: `y`来自结构侧(structure.py:39)，`label`来自序列侧(data_representer.py:195)。训练用`label`(ss.py:159)。`y`仅用于merge后断言，可选缓存

##### 第4轮: 最终实现规格

**LMDB缓存Schema（18个字段）**:

| 字段 | dtype | Shape | 说明 |
|------|-------|-------|------|
| `str_tag` | str | scalar | 'complex' 或 'ligand' |
| `sample_weight` | float32 | [] | 样本权重（当前=1.0） |
| `mask_use_complex` | float32 | [N_total] | complex分支掩码 |
| `mask_use_ligand` | float32 | [N_total] | ligand分支掩码 |
| `y` | float32 | [] | 可选，仅用于断言 |
| `protein_x` | float32 | [N_total, 28] | 蛋白原子特征（含零填充） |
| `ligand_x` | float32 | [N_total, F_lig] | 配体原子特征（含零填充） |
| `ligand_mask` | float32 | [N_total] | 配体原子掩码 |
| `protein_mask` | float32 | [N_total] | 蛋白原子掩码 |
| `ligand_index` | int64 | [N_total] | 配体索引（蛋白位填280） |
| `real_edge_index` | int64 | [2, Er] | 真实化学键边索引 |
| `knn_edge_index` | int64 | [2, Ek] | k-NN邻近边索引 |
| `real_dist` | float32 | [Er] | 真实键原始距离（无噪声/无smearing） |
| `knn_dist` | float32 | [Ek] | knn边原始距离（无噪声/无smearing） |
| `real_bond_onehot` | float32 | [Er, 6] | 真实键类型one-hot（含12→5映射） |
| `knn_bond_onehot` | float32 | [Ek, 6] | knn键类型one-hot（全为class-5） |
| `real_cross` | float32 | [Er, 1] | 真实键交叉标志（全0） |
| `knn_cross` | float32 | [Ek, 1] | knn边交叉标志 |

**运行时组装逻辑**:
- `complex_edge_index = cat([real_edge_index, knn_edge_index])` — 始终 [real, knn]
- **fixed模式**: `complex_edge_attr = cat([real_attr, knn_attr])` — 与index对齐
- **legacy_bug模式**: `complex_edge_attr = cat([knn_attr, real_attr])` — 复现原始bug
- dist_noise增强: 运行时对raw距离加Laplace噪声 → Gaussian smearing → 拼接

**文件结构**:
```
scripts/step9/
├── cache_utils.py              # ~400行
│   ├── BuildStructureCacheData  # 离线构建transform（替代EdgeConnection，无噪声/无smearing）
│   ├── RebuildComplexEdgeAttr   # 运行时重建transform（支持fixed/legacy_bug模式）
│   ├── CachedStructureSequenceDataset  # 序列live + 结构cached数据集
│   └── ordered_intersection()   # 修复set()排序bug
├── build_structure_cache.py    # ~200行，离线LMDB缓存构建器
└── main_training_cached.py     # ~200行，缓存训练脚本（替代main_training.py）
```

**关键设计决策**:
1. 分区存储（real/knn分开） → 一次构建支持两种运行模式
2. 不缓存`complex_edge_attr` → 运行时重建保留dist_noise随机增强
3. 断言 `max_substrate_length == 280`（多处硬编码，暂不改）
4. Windows lazy LMDB: `__getstate__`剥离handles，worker内`_ensure_*_connections()`延迟重连
5. valid_idx用ordered_intersection修复set()非确定性排序bug
6. 缓存元数据: `cache_schema_version`, `cache_builder_k`, `cache_builder_max_substrate_length`

#### 冒烟测试结果 (2026-03-11)
| 测试项 | 结果 | 详情 |
|--------|------|------|
| 基础导入 | ✅ | BuildStructureCacheData, RebuildComplexEdgeAttr, ordered_intersection |
| 维度验证 | ✅ | protein_x=[N,28], ligand_x=[N,34], ATOM_FEATS.Hybridization=9 |
| 端到端构建 | ✅ | 5个真实样本, 平均~960KB/样本 |
| 运行时重建(fixed) | ✅ | edge_attr=[N_edges, 39], edge_index=[2, N_edges] |
| 运行时重建(legacy_bug) | ✅ | attr排序与fixed正确交换 |
| 缓存往返(pickle) | ✅ | 序列化→反序列化→重建, 数据完全一致 |

#### 构建脚本改进 (Codex Round 5-7)
- **断点续传**: LMDB存在时扫描已有keys, 跳过已缓存条目
- **计数器分离**: `dup`(本次运行重复) vs `skipped_resume`(上次已缓存)
- **后置验证**: 构建完成后重新扫描LMDB, 比对expected vs actual entry counts
- **LMDB参数**: `readahead=False, meminit=False` for scan, `map_size=512GB` for write

#### 数据规模估算
| dataset_id | 名称 | 总行数 | 唯一Dock | 估计大小 |
|------------|------|--------|---------|---------|
| 0 | brenda | 359,979 | 157,929 | ~150GB |
| 1 | Duf | 1,669 | 1,669 | ~1.6GB |
| 2 | Esterase | 8,432 | 8,432 | ~8GB |
| 3 | Gt_acceptor | 3,042 | 3,042 | ~2.9GB |
| 4 | Nitrilase | 429 | 429 | ~0.4GB |
| 5 | Phosphatase | 21,630 | 21,630 | ~20GB |
| 6 | Thiolase | 683 | 683 | ~0.6GB |
| **合计** | | **395,864** | **193,814** | **~184GB** |

注: 实际缓存大小取决于StructureDataset过滤后的有效样本数(176,843)和去重后的唯一dock数。

#### 创新点价值
- **发现并修复论文原始代码中的边特征对齐Bug**
  - `complex_edge_attr`与`complex_edge_index`边排序不匹配
  - 真实化学键丢失键类型身份和正确距离
  - 位置: `src/Datasets/Structure/transforms.py:130-147`
- 修复后结构模块（SE(3)-EGNN）可正确学习边-属性对应关系
- 可解释Step 8 E8v2中结构通道变化不显著的现象（结构通道在被错误数据训练时已被抑制）
- 可作为毕业设计创新贡献点

## 进度
- [x] 数据拷贝到本地
- [x] 代码分析 + 修复（CSV、LMDB、num_workers）
- [x] 训练脚本编写 + 冒烟测试 + Pipeline测试
- [x] 方法1参数调优 + 基准测试（4.3 it/s, GPU 29%）
- [x] 方法2设计（Codex 3轮辩论）
- [x] 方法2验证（Codex 4轮验证咨询，发现边排序Bug）
- [x] 方法2实施（3文件: cache_utils.py, build_structure_cache.py, main_training_cached.py）
- [x] 方法2冒烟测试（导入/端到端/运行时重建/缓存往返 全部通过）
- [x] 方法2 Codex审核（累计7轮: 4轮设计 + 3轮构建脚本）
- [ ] 缓存构建（~2-4小时一次性）← **进行中**
- [ ] 正式训练（50 epochs, fixed模式）
- [ ] 可选: legacy_bug基线训练（对比用）
- [ ] 结果分析
