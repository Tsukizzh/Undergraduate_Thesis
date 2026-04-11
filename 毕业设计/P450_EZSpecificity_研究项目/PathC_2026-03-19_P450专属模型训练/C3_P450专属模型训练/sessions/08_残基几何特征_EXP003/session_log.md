# EXP003: 残基级几何特征 (φ/ψ/χ1)

## 实验概述

| 项目 | 值 |
|------|-----|
| 实验名 | EXP003_residue_geometry |
| 基线 | EXP002b (Test AUC=0.7889, lr=4e-4, warmup=12, wd=1e-5) |
| 创新点 | 残基级二面角特征 (φ/ψ/χ1) → EGNN节点 |
| feature_dim | 31 → **37** (+6: sin/cos of phi, psi, chi1) |
| 灵感来源 | EnzymeCAGE (Nature Chem Eng 2025) GVP-GNN |
| 日期 | 2026-04-12 |

## 架构改动

### 新增特征 (6维)

| 特征 | 物理含义 | 编码 | 维度 |
|------|---------|------|------|
| φ (phi) | 骨架弯折角1: C(i-1)─N─Cα─C | sin(φ), cos(φ) | 2 |
| ψ (psi) | 骨架弯折角2: N─Cα─C─N(i+1) | sin(ψ), cos(ψ) | 2 |
| χ1 (chi1) | 侧链朝向: N─Cα─Cβ─Xγ | sin(χ1), cos(χ1) | 2 |

- sin/cos编码处理角度周期性 (-180°≈+180°)
- 残基→原子广播: 同残基所有原子共享角度值
- 配体原子和HETATM (Heme/Fe) 填0
- GLY/ALA 无χ1 → 填0
- 链末端残基缺φ或ψ → 填0
- 肽键断裂检查 (>2.0Å) 防止跨链误算

### χ1残基特异性映射

18种有侧链的氨基酸各自的第4原子:
- 多数: CG (ARG, ASN, ASP, GLN, GLU, HIS, LEU, LYS, MET, PHE, PRO, TRP, TYR)
- SER: OG, THR: OG1, CYS: SG, VAL: CG1, ILE: CG1

### 代码改动

| 文件 | 改动数 | 内容 |
|------|--------|------|
| protein_ligand.py | 8处 | `_calc_dihedral()` + `CHI1_ATOMS` + `_annotate_residue_angle_features()` + `to_dict_atom()` + `get_zero_protein_feature()` |
| transforms.py | 3处 | feature_dim 31→37, 拼接角度特征, None fallback |
| pt_dataset.py | 6处 | PROTEIN_GEOM_DIM=6, `rebuild_protein_x()` + 3处加载点 |
| build_pt_cache.py | 12处 | 提取/打包(float16)/存储/self_test(37维)/per-sample复制 |

### 不改的文件

- egnn.py: nn.Linear(37, 128) 自动适配
- ss.py: 不涉及结构特征维度
- structure.py / utils.py: 自动透传新字段

## Codex验证 (6轮)

| 轮次 | 验证内容 | 结果 |
|------|---------|------|
| R1 | 完整数据管道追踪 + 初版diff | 确认5个文件改动点 |
| R2 | 调用位置修正 + 命名链验证 | 修正`_annotate`位置在LOOP 2之后 |
| R3 | 4个具体issue验证 | 修正`get_zero_protein_feature()`遗漏 |
| R4 | 二面角数学正确性 + self_test | Praxeolitic公式正确, 37维确认 |
| R5 | 设计层面5点验证 | 命名链/维度/兼容性全部一致 |
| R6 | 读已打patch文件实物验证 | 4文件验证通过, 1个stale注释已修 |

## 数据管道

```
PDB pocket → protein_ligand.py._parse() + _annotate_residue_angle_features()
  → to_dict_atom() {residue_angle_feature: (N_atoms, 6)}
  → structure.py: auto-prefix → data.protein_residue_angle_feature
  → pickle into structure_features_geom.lmdb
  → build_pt_cache.py: extract as protein_residue_geom (float16)
  → .pt shards/per-sample files
  → pt_dataset.py: load → rebuild_protein_x() → 37-dim protein features
  → EGNN: nn.Linear(37, 128) → message passing
```

## 实验进度

- [x] Structure LMDB重建 (structure_features_geom.lmdb) ✅ 48,750条, 1.5GB
- [x] pt_cache重建 (pt_cache_geom/) ✅ 47,807 samples, self-test dim=37 PASSED
- [x] 训练 EXP003 (random_split, 4×RTX5090) ✅ 73 epochs, ~1:04/ep, ~78min
- [x] 结果对比 EXP002b (Test AUC=0.7889) ✅

## 最终结果

| 指标 | EXP003 | EXP002b (基线) | 变化 |
|------|--------|----------------|------|
| **Test AUC-ROC** | **0.7914** | 0.7889 | **+0.0025** |
| **Test AUPR** | **0.3814** | 0.3542 | **+0.0272** |
| Val AUC (best) | 0.7850 (ep58) | — | — |
| Early stop | ep73 (patience=15) | — | — |
| Test samples | 11,924 (pos=1,109, neg=10,815) | — | — |
| Feature dim | 37 | 31 | +6 |
| Trainable params | 1,847,812 | 1,847,044 | +768 (128×6) |
| Training | 4×RTX5090, ~1:04/ep, 73ep, ~78min | — | — |

### AUC提升链

EXP001(0.7730) → EXP002a(0.782, +Fe/HEM) → EXP002b(0.7889, +调参) → **EXP003(0.7914, +φ/ψ/χ1)**

### 结论

- φ/ψ/χ1残基级几何特征有效提升性能
- Test AUC +0.0025, Test AUPR +0.0272（AUPR提升显著）
- Step 14（双尺度结构编码）值得推进——角度信息在独立GNN通道中可能效果更好

**日期**: 2026-04-12
