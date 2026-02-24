# Step 2 待办清单

> 上次更新：2026-02-23
> 上次状态：代码+环境全部就绪，**卡在缺少大文件**

---

## 立即要做（文件迁移）

从 Windows 机器复制以下文件到 Mac（或直接在 Windows 上运行）：

### 必需文件 — Step 8.1 需要

| 文件 | 大小 | Windows 路径 | Mac 目标路径 |
|------|------|-------------|-------------|
| 627 个 PDB 对接文件 | ~数 GB | `PathA_.../data/01_Step1_PDB文件/` | 同路径（PathA 下） |

### 必需文件 — Step 9 需要

| 文件 | 大小 | Windows 路径 | Mac 目标路径 |
|------|------|-------------|-------------|
| enzyme_features.lmdb | 4.29 GB | PathA Step 5 产出 | `PathB/data/00_shared/features/` |
| reaction_features.lmdb | 4.29 GB | PathA Step 6 产出 | `PathB/data/00_shared/features/` |
| grover_fingerprint.lmdb | 10.74 GB | PathA Step 7 产出 | `PathB/data/00_shared/features/` |
| morgan_fingerprint.npy | 447 KB | PathA Step 7 产出 | `PathB/data/00_shared/features/` |

**总大小约 22 GB**。如果用外接硬盘/U盘传输，预计需要一些时间。

**替代方案**：如果不方便传文件，也可以直接在 Windows 机器上运行（需先装依赖：`pip install rdkit lmdb torch_geometric`）。

---

## 文件迁移后执行

### 1. 先跑单个实验验证

```bash
conda activate torch  # 或你的环境名
cd "PathB_2026-02-12_P450数据集构建与结构优化/"

# 跑 EXP02（10Å + Heme）— 最能验证 HETATM 修复是否生效
python scripts/utils/run_experiment.py \
  --name EXP02_B6_10A_Heme \
  --pocket_radius 10.0 \
  --include_heme
```

**验证要点**：
- Step 8.1 输出的 pocket PDB 中应包含 HETATM 行
- Step 8.3 日志中 `n_hetatm > 0`（表示 Heme 原子被成功读取）
- 生成的 LMDB 文件不为空

### 2. 跑全部 4 个实验

```bash
python scripts/utils/run_experiment.py --run_all
```

### 3. 编写 Step 9/10 脚本

Step 9（推理）和 Step 10（分析）的 PathB 版脚本还没写。需要：
- Step 9：可配置的特征路径 + 结构 LMDB 路径 + 输出路径
- Step 10：4 组实验的 AUC-ROC 对比分析

### 4. 让 Codex 审核运行结果

用户要求"Codex 也帮你验证，你们一起审核代码，确保正确"。

---

## 已完成的工作回顾

| 完成项 | 说明 |
|--------|------|
| Step 1 数据迁移 | B6 数据集 + 共享特征（Windows 上） |
| extract_pocket_ligand.py | 支持 --pocket_radius 和 --include_heme |
| step8_align_ligand.py | PathB 版，CLI 可配置 |
| step8_generate_structure_lmdb.py | PathB 版，HETATM 修复 |
| run_experiment.py | 集成 8.1→8.2→8.3，智能排序 |
| Codex 三轮审核 | 修复 7 个 bug |
| MacBook 环境 | RDKit + lmdb + PyG 已安装 |
| dry-run | 4 个实验命令验证通过 |
