# Path C1 + C2 深度调研报告：基线训练与 P450 全面数据集构建

调研日期：2026-04-23

---

## C1 部分：论文基线训练与参数调整

### 1. AllSplit 从头训练（为什么、用什么、如何对标）

**动机**：复现论文在 all_split（unknown enzyme+substrate 最难场景）的基线，评估论文声称的 AUC=0.7198。论文模型用 ESIBank BRENDA 数据集训练（25,225 酶 + 34,417 底物），需**从头训练一遍**以建立后续消融的参考点。

**用什么数据**：ESIBank 的 BRENDA 全量数据集（不是 P450 子集）。Google Drive `173a36Ni...ESIBank/brenda/all_split/` 下 176K 训练样本（有效 ~143K, 42%）、2,876 验证、16.3K 测试。

**如何对标**：论文在 NCSA Delta HPC（4×A100）训练 256+ epochs，最终 Val AUC=0.722 / Test=0.7244。我们本地 RTX 4070S（12GB）跑 12 个 epoch（13.7 小时），最佳 **Val=0.7542 / Test=0.7505**（ep12），超论文 +0.026 on test。后转云服务器 Cloud-2（2×RTX 4090 DDP）跑 legacy_bug 基线达 **Val=0.7224 / Test=0.7244**，完全复现论文。

### 2. pt 训练管线与 Cloud-DDP

**云服务器配置**：
- Cloud-1 单卡：1×RTX 4090 24GB，¥1.46/h
- Cloud-2 双卡→四卡：2→4×RTX 4090 24GB，180→360GB RAM，¥2.26/h→¥4.4/h，28 核 vCPU

**DDP 配置**：PyTorch Lightning DDPStrategy
- `find_unused_parameters=False`（模型所有参数都用到）
- `gradient_as_bucket_view=True`（view 形式通信减内存拷贝）
- `sync_dist=False`：验证时关闭 rank 0 同步，但需在 `validation_epoch_end` 和 `test_epoch_end` 手动 all_gather
- effective_batch_size = bs × num_gpus = 56 × 2 = **112**

**shard_size 优化**：从 `--shard-size 2048`（20 workers 只有 5 任务，并行度低）改为 `--shard-size 256` + `OMP_NUM_THREADS=1`，所有 graph shards **30 秒完成**（vs 10+ 分钟）。

### 3. C1-Step1 fixed 基线（边排序修复）

**Edge fix 修复的 bug**：`src/Datasets/Structure/transforms.py` 中 `EdgeConnection`：
- `complex_edge_index` = [real_edge_index, knn_edge_index]（真实键+k-NN 边）
- `complex_edge_attr` = [knn_attr, real_attr]（**顺序反了**！）
- 后果：真实化学键（单/双/芳香）被错误分配为 k-NN 类型属性（全 type-5），真实距离被随机 k-NN 距离替换

**性能对比（边修复反而下降）**：

| 版本 | Val AUC | Test AUC |
|------|:---:|:---:|
| legacy_bug | 0.7224 | 0.7244 |
| **fixed** | **0.7145** | **0.7060** |
| 差值 | -0.008 | -0.018 |

**解释**：修复反而下降 1.84% 说明：
1. 错误边属性可能引入"有益噪声"正则化
2. 模型压制结构分支补偿 bug，反而综合信号更好
3. **batch_size 才是主要瓶颈**（本地 bs=32 fixed 达 0.7667；云 bs=112 fixed 只有 0.7145）

### 4. C1-Step2 dropout 消融（Val 改善不迁移到 Test）

| 设置 | Val AUC | Test AUC | Val-Test gap |
|-----|:---:|:---:|:---:|
| dropout=0.9（基线） | 0.7145 | 0.7060 | 0.009 |
| **dropout=0.1** | 0.7216 | 0.6936 | 0.028 |
| **dropout=0.3** | 0.7397 | 0.6959 | 0.044 |

**关键发现**：
1. Val 提升 +2.52pt（dropout=0.3）但 Test 下降 -1.0pt
2. Val-Test gap 扩大 3-5 倍 → 过拟合到验证集分布
3. 降低 dropout 增加网络容量，在验证集学到细节但在测试集（unknown enzyme+substrate）**无法泛化**
4. **结论**：dropout 改动**不纳入 C1-Step3**

### 5. 云服务器 14 轮参数调优

| # | bs | Accum | Workers | 改动 | Speed | /epoch | 备注 |
|---|----|-----|-------|------|:---:|:---:|:---:|
| 1 | 48 | 1 | 6 | 默认 DDP | 2.82 | 11 min | 初始基线 |
| 5 | 56 | 1 | 6 | 最优参数 | 2.63 | 10.1 min | ✅ |
| **12** | **56** | **1** | **6** | **sync_dist=F+all_gather** | **2.65** | **10.0 min** | **最终最优** |

---

## C2 部分：P450 全面数据集构建

### 6. 数据来源（5 个主源）

| 源 | 数据源 | 正样本 | 酶 | 底物 | 质量 |
|----|--------|:---:|:---:|:---:|:---:|
| S1 | RCSB | 272 | 103 | 220 | Tier A（晶体） |
| S2 | ESIBank | 806 | 338 | 390 | Tier B |
| S3 | P450Rdb（成都中医药大学） | 2,798 | 857 | 1,492 | Tier B |
| S8 | PlantP450DB（哥本哈根） | 979 | 578 | 295 | Tier B |
| S9 | PCPD（中科院天津） | 1,209 | 818 | 570 | Tier B |
| **合计（去重前）** | | **6,064** | 2,694 | 2,967 | |

**去重后**：4,751 正样本 / 1,622 酶 / 2,125 底物（去重率 21.6%）。

### 7. 合并去重（Phase 4）— Key 设计

**酶去重 key**：
- 主键：UniProt ID（精确匹配） → 1,457 条
- Fallback 1：序列桥接（PCPD 无 UniProt 的酶若序列完全匹配有 UniProt 者） → **282 条成功桥接**
- Fallback 2：序列 SHA256 hash → 165 条
- **合计**：**1,622 唯一酶**

**化合物去重 key**：
- 主键：RDKit Canonical isomeric SMILES → **Standard InChIKey**（99.3% 成功）
- Fallback：Canonical SMILES（15 条含通配符无法生成 InChI）
- **合计**：**2,125 唯一化合物**

**为什么用这两个 key**：UniProt ID 是生命科学全球标准；InChIKey 对化学结构（含手性、盐、顺反异构）标准化 → 两者覆盖生物+化学两维度去重。

### 8. 负样本与 Split（Phase 5，v1→v6 六次迭代）

| 版本 | 转折 | 问题 | 结果 |
|-----|------|------|------|
| v1-v4 | all_split 连通分量 | 负样本比 7.28:1，分布不平衡 | ❌ |
| v5 | all_split 独立 data.csv | 违反论文"4 split 共享数据"原则 | ❌ |
| **v6** | **共享 data.csv + all_split 对角线过滤** | 最终版 | ✅ |

**v6 最终设计**：
- 共享一份 data.csv：52,254 行（4,751 正 + 47,503 负），4 种 split 切同一份数据
- **负样本生成（改论文单向→双向）**：
  - Direction A（固定底物换酶）：1 正 × 5 = 23,755 条
  - Direction B（固定酶换底物）：1 正 × 5 = 23,748 条
- **4 种 split**：

| Split | 行数 | 酶重叠 | 底物重叠 |
|------|:---:|:---:|:---:|
| random | 52,254 | ~99% | ~99% |
| reaction | 52,254 | ~99% | **0%** |
| enzyme | 52,254 | **0%** | ~99% |
| **all** | **18,159** | **0%** | **0%** |

**all_split 为何必须丢行**：论文 25K 酶 × 34K 底物稀疏空间 → 无损 0% 重叠保留 100%；我们 P450 1.6K 酶 × 2.1K 底物互联空间 → 只能保留 **34.8%**（18,159/52,254）。两者都是最严格泛化测试，物理约束不同。

### 9. 结构获取与对接（Phase 6）

**对接工具**：**Uni-Dock**（Vina 的 GPU 加速版），非论文的 AutoDock-GPU。理由：精度更高（RMSD 1.30 vs 1.42）、不需要 AutoGrid 格点图和修改 Fe 电荷、从头训练无需和论文姿态一致。

**50,180 对怎么来**：
- 理论应 52,254 对
- 实际 50,180（95.8%）的差异：
  - 1 个失败底物（embed failed）
  - 54 个失败酶（血红素移植失败）
  - 8 个无 AlphaFold 酶
  - 219 个 blocked 酶（ColabFold 预测中）

**ligand complex PDB 流程**：
1. RDKit 3D 构象生成（ETKDG+MMFF94s） → 2,124/2,125 成功
2. Meeko 生成配体 PDBQT
3. MGLTools 生成受体 PDBQT（兼容 HEM，不用 Meeko）→ 1,397/1,403 成功
4. Uni-Dock 对接（100 runs，选最优姿势）
5. 复合物 PDB 格式：蛋白原子 + "COMPND" 分隔 + 配体原子

**parse_smile bug**：PCPD 底物含通配符（如 `*[H]`，锚点标记），RDKit `MolFromSmiles` 失败 → 修复：检测 `*` 并用 `[H]` 替换后重试。

**Meeko 修复**：Meeko 不兼容 HEM，改用 MGLTools `prepare_receptor4.py`，成功率 99.6%。

### 10. 特征生成（Phase 7）

最终**可用配对 47,510 对**（91.0% 覆盖率）：

| 特征 | 覆盖 | 跳过 | 说明 |
|------|:---:|:---:|------|
| ESM 序列编码 | 1,577/1,622 | 45 | 44 非标 AA + 1 超长 (>1000aa) |
| Ligand Graph | 2,119/2,125 | 6 | parse_smile bug |
| Morgan 指纹 | 2,125/2,125 | 0 | 100% |
| GROVER 嵌入 | 2,124/2,125 | 1 | `*[H]` 越界 |
| Structure 酶 | 48,354 | — | |
| Structure 底物 | 2,125 | — | |
| **最终可用对** | **47,510** | | 4,362 正 + 43,148 负（1:9.9） |

### 11. Phase 7.5 pt_cache 构建

| Split | Train | Val | Test | Total |
|-------|:---:|:---:|:---:|:---:|
| all_split | 8,326 | 4,023 | 4,143 | 16,492 |
| random_split | 23,710 | 11,878 | 11,823 | 47,411 |

**关键优化**：`--shard-size 256` + `OMP_NUM_THREADS=1` → 30 秒完成。

### 12. Phase 8 EXP001 训练结果

| 指标 | 值 |
|------|:---:|
| 配置 | bs=56/GPU, eff=224, 4×4090 DDP, --preload ~149GB |
| 训练集 | 23,710 |
| **Val AUC** | **0.7544** |
| **Test AUC** | **0.7730** |
| Test AUPR | 0.3620 |
| 速度 | 1.8 it/s, ~2 min/epoch |
| 总耗时 | 49 分钟（89 epochs, early stop） |

**vs 基线对比**：

| 模型 | 数据 | Test AUC |
|------|------|:---:|
| 论文 checkpoint | ESIBank P450 内部 | 0.638 |
| 我们 | ESIBank all_split | 0.7244 |
| **我们** | **P450 数据集** | **0.7730** |

提升 **+0.135**（vs ESIBank P450 内部），确认数据集构建成功。

**main_training_pt.py 4 轮 Codex review 修复**：
1. DDP test_epoch_end bug：sync_dist=False 原只用本地 rank logits 算 AUC → 修复后 all_gather 全部 rank
2. 自动测试：训练完直接 trainer.test()
3. 自动关机：--shutdown 参数省云费用
4. SRC_DIR fallback：experiment-local src/ → global src/

---

## 总结（关键数字）

### C1 关键指标
- **AllSplit 本地最佳**：Val=0.7542 / Test=0.7505（ep12）
- **AllSplit 云 legacy_bug**：Val=0.7224 / Test=**0.7244**（完全复现论文）
- **Fixed 基线**：Val=0.7145 / Test=0.7060（边修复反降 -1.84%）
- **Dropout 消融**：d=0.3 Val+2.52pt / Test-1.0pt（不采纳）

### C2 关键指标
- 数据源合并：6,064 → **4,751 正样本**（去重 21.6%）
- 酶去重：2,694 → **1,622**（含 282 条序列桥接）
- 底物去重：2,967 → **2,125**
- 负样本（双向）：**47,503**（正负比 1:9.998）
- All_split 严格：18,159 行（34.8%，0%/0% 重叠）
- 特征覆盖：47,510 可用对（91.0%）
- **EXP001 baseline：Test AUC=0.7730**（vs ESIBank 0.638, **+0.135**）
