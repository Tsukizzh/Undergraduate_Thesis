# Session 06: Fe/血红素编码 — EXP002a

> 日期: 2026-04-01
> 实验: EXP002a_fe_heme_10A
> 服务器: Cloud-2 (4×RTX4090, 360GB RAM, 28核)
> 环境: conda ezspec (Python 3.10)

---

## 一、实验目标

P450酶的催化完全依赖血红素中心的Fe原子，但EZSpecificity原始模型：
1. PDB解析器跳过所有HETATM行 → 血红素从第一步就被丢弃
2. 蛋白质原子词汇表没有Fe → 即使进来也编码为全零
3. 氨基酸表没有HEM → 被标为UNK（未知）

本实验将Fe和血红素纳入模型的结构通道，评估对预测性能的影响。

## 二、实验前数据验证

### Fe-配体距离分析（20个含HEM的对接结构抽样）

| 指标 | 值 |
|------|-----|
| Fe到底物最近原子平均距离 | 4.9Å |
| 中位距离 | 4.5Å |
| 10Å内覆盖率 | 90% (18/20) |
| PDB中HEM存在率 | 95.5% (191/200) |

### 正负样本Fe距离对比（各50个样本）

| | 正样本 | 负样本 |
|---|---|---|
| 平均距离 | 5.0Å | 5.3Å |
| p值 | 0.637（不显著） |

结论：Fe距离本身不能区分正负样本（Uni-Dock都把底物对接到Fe上方），Fe的价值在于提供催化中心的完整3D上下文。

### 关键发现：旧pocket文件已含HEM

旧的口袋提取脚本(phase7_step1_pocket.py)对HETATM行不做距离筛选直接写入pocket。验证8/10个抽样pocket文件包含HEM行。因此**不需要重新提取口袋**，只需重生成结构特征LMDB。

## 三、代码改动

### 3.1 protein_ligand.py（7处改动）

Codex审核通过（session 019d4534）。

1. 新增类属性: `HEM_AA_TYPE=21`, `HEM_RES_NAMES={'HEM','HEC'}`, `AA_NAME_NUMBER['HEM']=21`, `AA_NAME_NUMBER['HEC']=21`
2. `__init__`: 新增 `self.is_hetero=[]`
3. `_enum_formatted_atom_lines()`: 读HETATM中的HEM/HEC残基，**过滤在yield之前**（Codex发现的bug修正）
4. `_parse()`: HEM原子 aa_type=21, is_backbone=False, is_hetero=True
5. `to_dict_atom()`: 输出字典加 is_hetero
6. `get_zero_protein_feature()`: 加 is_hetero
7. `query_residues_ligand()`: 默认criterion从 center_of_mass → min_atom_distance

### 3.2 transforms.py（3处改动）

1. `atomic_numbers`: [1,6,7,8,16,34] → [1,6,7,8,16,34,26]（+Fe）
2. `max_num_aa`: 21 → 22; `feature_dim`: 28 → 31（+Fe+HEM+is_hetero）
3. `__call__()`: 拼接时加入 is_hetero

### 3.3 build_pt_cache.py（额外发现的改动）

- SRC_DIR: 改为优先使用实验本地src/
- PROTEIN_ATOMIC_NUMBERS: 加Fe
- protein_is_hetero: 在提取、chunk收集、shard存储、shard读取、per-sample .pt共7处添加
- self-test: num_classes 21→22, 加is_hetero

### 3.4 pt_dataset.py（Codex发现的关键遗漏）

pt_dataset.py有自己的特征重建逻辑（不走transforms.py），硬编码了旧参数。
- PROTEIN_ELEMENTS: 加Fe
- PROTEIN_NUM_AA: 21→22
- PROTEIN_FEAT_DIM: 28→31
- rebuild_protein_x(): 加protein_is_hetero参数
- 3处shard loading: 全部加protein_is_hetero

### 3.5 验证结果

端到端单条数据测试通过：
- pocket/0.pdb: 423个蛋白原子（含74个HEM原子，1个Fe）
- Fe编码: element=[0,0,0,0,0,0,1], aa=HEM(idx21), is_backbone=0, is_hetero=1
- 特征重建: [423, 31]维 ✓

## 四、数据管线

| Step | 内容 | 耗时 | 输出 |
|------|------|------|------|
| 1 | 重生成structure_features_heme.lmdb | 20秒 | 1.1GB, 48,750条 |
| 2 | 构建pt_cache_heme/random | 21秒 | 7.6GB, train/val/test |
| 3 | 训练EXP002a | ~2.5小时 | 114 epochs |

不需要重新生成的: ESM嵌入、Morgan指纹、GROVER指纹、分子图特征、口袋PDB。

## 五、训练结果

### 训练配置（与EXP001完全一致，仅cache-dir不同）

| 参数 | 值 |
|------|-----|
| 环境 | conda ezspec, Python 3.10 |
| batch_size | 56/GPU × 4 = 224 |
| lr | 3e-4, warmup 8 epochs |
| dropout | 0.9 |
| edge-mode | fixed |
| max-epochs | 200 (early stopped at 114) |

### 训练过程

```
阶段1 (ep0-8):   warmup, AUC 0.46→0.60
阶段2 (ep8-73):  lr=3e-4, AUC 0.60→0.774
阶段3 (ep75-95): lr=1.5e-4 (ReduceLR触发), AUC 0.774→0.778
阶段4 (ep96-109): lr=7.5e-5, AUC平稳
阶段5 (ep110-114): lr=3.75e-5, early stop
Best: ep99, Val AUC=0.7784
```

### 最终对比

| | EXP001 (无Fe) | EXP002a (含Fe) | 差异 |
|---|---|---|---|
| Best Val AUC | 0.754 (ep73) | 0.778 (ep99) | **+0.024** |
| **Test AUC-ROC** | **0.773** | **0.782** | **+0.009** |
| Test AUPR | 0.362 | 0.352 | -0.010 |
| Best Epoch | 73 | 99 | +26 |
| Total Epochs | 89 | 114 | +25 |
| 训练速度 | ~2 min/epoch | ~1.4 min/epoch | 更快 |

### 结论

加入Fe/HEM后Test AUC从0.773提升到0.782（+0.009），Val AUC提升更明显（+0.024）。提升幅度有限，可能原因：
1. 正负样本的Fe-配体距离无显著差异（对接都把底物放到Fe附近）
2. dropout=0.9可能淹没了血红素信号（~43节点 / ~400+总节点）
3. 口袋中已有的蛋白质残基已隐式编码了催化位点信息

## 六、经验教训

### 启动错误（浪费了约30分钟GPU时间）

1. **用错Python环境**: 用了conda base(3.11)而不是ezspec(3.10) → 训练能跑但速度慢
2. **旧进程没杀干净**: killall没杀到所有DDP子进程
3. **metrics.csv没清**: 多次启动导致数据叠加，绘图出现多条线
4. **TensorBoard events没清**: 同上
5. **脚本模板vs实际版本差异**: P450/scripts/模板和EXP001实际运行的版本不同

### 代码遗漏（被Codex发现）

1. **pt_dataset.py忘了改**: 它有自己的特征重建逻辑，硬编码了旧的28维参数
2. **build_pt_cache.py SRC_DIR指向全局src**: 导致第一次build用了旧的28维特征
3. **build_pt_cache.py第三处shard loading缺少protein_is_hetero**: grep所有is_backbone出现位置逐一核对才发现

### 教训总结

- 改特征维度时，必须grep所有相关硬编码（PROTEIN_ELEMENTS, PROTEIN_NUM_AA, PROTEIN_FEAT_DIM等），逐个文件核对
- 启动训练前：killall python → nvidia-smi确认GPU干净 → 清空旧metrics/events → 确认conda env
- 每改一处代码都用Codex审查 + 端到端单条数据测试验证

## 七、文件位置

```
服务器:
  实验目录: ~/rivermind-data/EZSpecificity/PathC/P450/experiments/EXP002a_fe_heme_10A/
  结构特征: ~/rivermind-data/EZSpecificity/PathC/P450/data/structure/structure_features_heme.lmdb
  训练缓存: ~/rivermind-data/EZSpecificity/PathC/P450/data/pt_cache_heme/random/

本地:
  结果: C3_P450专属模型训练/results/06_Fe血红素_EXP002a/
  Session: C3_P450专属模型训练/sessions/06_Fe血红素编码_EXP002a/
```
