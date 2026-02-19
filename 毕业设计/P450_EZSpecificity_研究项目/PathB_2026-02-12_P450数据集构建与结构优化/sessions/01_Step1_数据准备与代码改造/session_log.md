# Step 1: 数据准备与代码改造 — Session Log

## 执行日期
2026-02-12

## 概述
完成 PathB Step 1 的全部子任务：B6 数据集迁移、共享特征迁移、代码改造、实验运行器编写。

---

## 子任务执行记录

### 1.1 B6 数据集迁移

**来源**: `PathA_.../data/04_Step4_格式修正后数据/`
- `B6_底物正抑制剂负_271pos245neg/data.csv` → B6 子目录
- `Enzymes.csv` → Step4 根目录（非 B6 子目录内，Codex 审核发现）
- `Substrates.csv` → Step4 根目录

**目标**: `data/00_shared/datasets/B6_v1/`

**验证结果**:
| 文件 | 行数（含表头） | 内容 |
|------|:---:|------|
| data.csv | 517 | 516 条记录（272 正 + 244 负） |
| Enzymes.csv | 293 | 292 个酶序列 |
| Substrates.csv | 437 | 436 个底物 SMILES |

**交叉验证**: PASSED — data.csv 引用的所有 Enzyme Index / Substrate Index 均存在于对应 CSV 中。

**注意**: 实际标签分布为 272/244（非计划文档中的 271/245），以实际数据为准。

### 1.2 共享特征迁移

**来源/目标**: PathA Steps 5-7 → `data/00_shared/features/`

| 特征文件 | 大小 | 来源 | 验证 |
|---------|:---:|------|:---:|
| enzyme_features.lmdb | 4.29 GB | Step 5 ESM 嵌入 | 大小匹配 |
| reaction_features.lmdb | 4.29 GB | Step 6 反应图 | 大小匹配 |
| grover_fingerprint.lmdb | 10.74 GB | Step 7 GROVER | 大小匹配 |
| morgan_fingerprint.npy | 447 KB | Step 7 Morgan | shape=(436,1024) int8 |

**LMDB 完整性**: Windows 上因大文件 mmap 限制，无法直接 lmdb.open 读取 key 数量。改用文件大小完全匹配验证（源/目标字节数一致）。

### 1.3 代码改造

**产出**: `scripts/02_Step2_因子实验/extract_pocket_ligand.py`

**与 PathA 原版的关键差异**:

| 特性 | PathA | PathB |
|------|-------|-------|
| 口袋半径 | 硬编码 10.0Å | CLI `--pocket_radius`（默认 10.0） |
| Heme 处理 | 排除所有 HETATM | `--include_heme` 启用时纳入 |
| 路径配置 | 全部硬编码 | CLI 参数（`--mapping_csv`, `--pdb_dir`, `--output_dir`） |
| Heme 残基集 | N/A | HEM, HEC, HEA, HEO, DHE, HEB |

**Heme 设计决策**:
- Heme 原子始终全量加入 pocket，不受 radius 过滤
- 仅排除氢原子
- 原因：保证 2×2 因子实验的正交性

**Codex 审核要点**:
1. 下游 `protein_ligand.py` 只解析 ATOM 记录 → heme 的 HETATM 可能被忽略
2. 需在 Step 2 执行时同步修改 `generate_structure_lmdb.py`
3. `_dist()` 函数已修复空残基 edge case

### 1.4 实验运行器

**产出**: `scripts/utils/run_experiment.py`

**功能**:
- 自动解析项目路径（PathB root → PathA → Project root）
- 单实验模式：`--name EXP01 --pocket_radius 10.0 [--include_heme]`
- 批量模式：`--run_all --config base.yaml` 运行 4 组因子实验
- Dry-run 模式：`--dry_run` 仅打印命令不执行
- Preflight 检查：验证所有必需文件存在

**当前实现范围**: Step 8.1（提取）已完成，Steps 8.2-10 标记 TODO。

**Codex 审核修复**:
- 修复 `resolve_paths()` off-by-one 错误（`scripts/utils/` 需要 `.parent.parent`）
- 移除未使用的 `shutil` 导入
- PathA 选择改为排序后取第一个（确定性）

### 1.5 基础配置

**产出**: `configs/base.yaml`

包含 2×2 因子矩阵定义：
- EXP01: 10Å, 无 Heme（基线，与 PathA 条件一致）
- EXP02: 10Å, 有 Heme
- EXP03: 6Å, 无 Heme
- EXP04: 6Å, 有 Heme

---

## 产出文件清单

```
data/00_shared/
├── datasets/B6_v1/
│   ├── data.csv          (516 条)
│   ├── Enzymes.csv       (292 酶)
│   └── Substrates.csv    (436 底物)
└── features/
    ├── enzyme_features.lmdb    (4.29 GB)
    ├── reaction_features.lmdb  (4.29 GB)
    ├── grover_fingerprint.lmdb (10.74 GB)
    └── morgan_fingerprint.npy  (447 KB)

scripts/
├── 02_Step2_因子实验/
│   └── extract_pocket_ligand.py
└── utils/
    └── run_experiment.py

configs/
└── base.yaml
```

## 已知问题与 Step 2 待处理事项

1. **下游 HETATM 兼容性**: `step8_generate_structure_lmdb.py` 的 PDBProtein 解析器只读 ATOM 记录。Step 2 执行前需创建 PathB 版本，增加 HETATM 支持。
2. **Step 8.2 align_ligand.py**: 需要 PathB 版本（可配置路径），但对齐逻辑本身不受 heme/radius 影响。
3. **Step 9/10**: 需要 PathB 版本以支持可配置的特征路径和输出目录。
4. **B6 标签分布**: 实际为 272/244（非 271/245），需更新计划文档。
