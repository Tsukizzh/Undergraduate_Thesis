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

### 1. [阻断级] 下游 HETATM 兼容性（已验证）

**2026-02-20 Codex 四轮深度验证结论**：

**风险确认**：`extract_pocket_ligand.py` 的 `--include_heme` 在当前链路中**对模型输入无效**，2×2 实验会退化为只测半径。

**根因链路**：
1. `extract_pocket_ligand.py` 用 BioPython PDBIO 写出 pocket PDB
2. BioPython PDBIO 根据残基的 `hetfield` 决定写 `ATOM` 还是 `HETATM`（已查阅 PDBIO.py 源码确认）
3. Heme 残基的 `hetfield != " "`，因此写出为 `HETATM` 行
4. 下游 `PDBProtein._parse()` 只匹配 `line[0:6].strip() == 'ATOM'`，完全跳过 `HETATM`
5. 结果：Heme 原子被静默丢弃，`include_heme=True/False` 产生的 LMDB 数据相同

**受影响代码**：
- PathA 的 `step8_generate_structure_lmdb.py:127` — 独立 `PDBProtein` 副本，同样只读 ATOM
- `src/Datasets/Structure/protein_ligand.py:67` — 原始代码（gitignored，仅供参考）
- 两处实现完全一致，都只解析 ATOM 记录

**额外发现**：
- PDB 格式中 ATOM 与 HETATM 列布局完全一致（wwPDB 标准确认），可直接复用解析逻辑
- PathA 的 539 个 pocket PDB 中 HETATM 行数为 0（全量扫描已验证）
- Heme 原子量级：单个 HEM 残基约 43 个非氢重原子，占 pocket 平均原子数（~411）的 ~10.5%
- 元素编码：Fe(26) 不在 `FeaturizeProteinAtom.atomic_numbers` 中（会得到全零向量），但 Heme 中的 C/N/O 可正确编码
- `AA_NAME_NUMBER` 中无 "HEM"，Heme 原子会被标记为 `UNK`（index=20），与非标准残基处理一致
- residue key 冲突风险极低（实际 P450 数据中 ATOM 与 HETATM 的 chain+resSeq+iCode 无冲突）
- `_parse()` 不依赖 `atom_id` 连续性（用 `next_ptr = len(self.element)` 作内部索引）

**推荐修复方案**：方案 A — 在 PathB 创建 fork 版 `PDBProtein`，支持 `ATOM` + `HETATM`，不修改 `src/` 原始代码
- 改动点：`line[0:6].strip() == 'ATOM'` → `line[0:6].strip() in ('ATOM', 'HETATM')`
- 建议 residue key 增加 record_type 字段（防未来数据冲突）
- 增加日志：每样本的 `n_atom`, `n_hetatm`, `n_heme_atoms`, `n_fe_atoms`

**实验科学意义评估**：
- 修复后 Step 2 仍有价值，但定位为**推理敏感性/结构扰动实验**
- 不能直接宣称模型"学会了 Heme 生化机理"（训练分布未包含 HETATM，且 Fe 特征缺失）
- 若要机理层面结论，需在 Path C 重训时纳入 Heme 和 Fe 特征通道

### 2. Step 8.2 align_ligand.py — ✅ 已完成
PathB 版已创建（可配置路径 CLI），对齐逻辑与 PathA 一致。

### 3. Step 8.3 generate_structure_lmdb.py — ✅ 已完成
PathB 版已创建，核心改动：PDBProtein 支持 HETATM。

### 4. run_experiment.py — ✅ 已更新
已集成 Step 8.2 和 8.3。Steps 9-10 仍 TODO。

### 5. Step 9/10
需要 PathB 版本以支持可配置的特征路径和输出目录。

### 6. B6 标签分布
实际为 272/244（非 271/245），需更新计划文档。

---

## Codex 三轮代码审核记录（2026-02-21）

### 审核范围
3 个 PathB 脚本：`step8_generate_structure_lmdb.py`、`step8_align_ligand.py`、`run_experiment.py`

### 已修复问题（7 项）

| # | 严重度 | 文件 | 问题 | 修复 |
|---|--------|------|------|------|
| 1 | MEDIUM | step8_generate_structure_lmdb.py | `element_symb` 回退 `line[13:14]` 对 Fe 等双字母元素只取 1 字符 | → `line[12:14].strip().capitalize()` |
| 2 | MEDIUM | step8_generate_structure_lmdb.py | `r["success"] == "True"` 大小写敏感 | → `.strip().lower() == "true"` |
| 3 | LOW | step8_generate_structure_lmdb.py | AtomMapNum 全零检查误拒单原子配体 | → 加 `mol.GetNumAtoms() > 1` 守卫 |
| 4 | HIGH | run_experiment.py | 共享对齐实验顺序依赖 | → 排序 radius DESC + heme DESC |
| 5 | HIGH | run_experiment.py | `alignment_summary.exists()` 不检测空/损坏文件 | → 行数验证 + 重跑逻辑 |
| 6 | MEDIUM | step8_align_ligand.py | 缺失 raw_ligand 静默成功 | → 成功率 < 50% 返回错误码 1 |
| 7 | LOW | step8_generate_structure_lmdb.py | 空 target_docks 仍生成空 LMDB | → fail fast 返回 2 |

### 确认无问题
- `is_backbone` 只对 ATOM 记录为 True
- PDB 列布局 ATOM/HETATM 一致
- pickle 兼容性保持
- `project_root` 路径解析正确

### 接受的残余风险
- alignment_summary 验证用 `row_count > 0`（完美验证需循环依赖，当前已足够实用）
