# CLAUDE.md

本文件为 Claude Code (claude.ai/code) 在此代码库中工作时提供指导。

## 项目概述

EZSpecificity 是一个基于交叉注意力机制的 SE(3)-等变图神经网络，用于预测酶-底物特异性。它基于以下特征预测给定酶是否能催化特定底物的反应：
- 蛋白质序列嵌入（ESM）
- 分子指纹（Morgan、GROVER）
- 酶-底物对接复合物的 3D 结构特征

**Web 界面**：https://ezspecificity.platform.moleculemaker.org/home

**主要工作流程**：本项目专注于推理（inference）。使用 `src/example.ipynb` 获取完整的端到端流程。

## 快速开始

```bash
# 1. 创建环境
conda env create -f src/environment.yml
conda activate ezspecificity

# 2. 安装 PyTorch Geometric（激活 conda 环境后）
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html

# 3. 运行示例
cd src
jupyter notebook example.ipynb
```

## 系统要求

**测试环境**：
- 操作系统：Ubuntu 22.04（主要）、Windows 10/11（需调整）
- Python：3.10.12
- CUDA：12.1
- GPU：推荐 NVIDIA GPU ≥8GB 显存
- 磁盘：≥20GB 可用空间用于模型和特征

**Windows 注意事项**：
- 路径使用反斜杠 `\` 或原始字符串
- 某些 shell 命令可能需要 PowerShell/CMD 等效命令
- GROVER PYTHONPATH 设置：`set PYTHONPATH=<project_root>\src\other_softwares\grover_software;%PYTHONPATH%`

## 常用命令

### 环境设置

```bash
# 创建 conda 环境
conda env create -f src/environment.yml
conda activate ezspecificity

# 安装 PyTorch Geometric（激活 conda 环境后）
pip install torch_geometric
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-2.1.0+cu121.html
```

**注意**：PyTorch Geometric 包版本必须与您的 CUDA 版本匹配。对于 CUDA 12.1，使用上述 URL。对于其他 CUDA 版本，请查看 https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html

### 运行推理

```bash
# 主要方法：Jupyter notebook（推荐用于交互式探索）
cd src
jupyter notebook example.ipynb

# 替代方法：Python 脚本（如果在 Jupyter 中遇到 DDP 策略错误，请使用此方法）
cd src
python main_testing.py
```

### GROVER 特征生成

GROVER 生成分子指纹，需要多步骤流程：

```bash
# 导航到 GROVER 目录
cd src/other_softwares/grover_software

# 设置 PYTHONPATH（Linux/Mac）
export PYTHONPATH=<project_root>/src/other_softwares/grover_software:$PYTHONPATH

# 设置 PYTHONPATH（Windows）
set PYTHONPATH=<project_root>\src\other_softwares\grover_software;%PYTHONPATH%

# 步骤 1：生成特征
python scripts/save_features.py --data_path <input.csv> --save_path <output.npz> --features_generator fgtasklabel --restart

# 步骤 2：构建词汇表
python scripts/build_vocab.py --data_path <input.csv> --vocab_save_folder <vocab_dir> --dataset_name test

# 步骤 3：生成指纹
python main.py fingerprint --data_path <input.csv> --features_path <features.npz> --checkpoint_path <grover_large.pt> --fingerprint_source both --output <output.npz> --save_lmdb_path <output.lmdb>
```

**重要**：GROVER 检查点 `grover_large.pt` 应位于 `data/pretrain_model/grover_large.pt`。如果缺失，请从原始 GROVER 仓库下载。

### 验证步骤

安装后，验证您的设置：

```bash
# 1. 检查 conda 环境
conda activate ezspecificity
python -c "import torch; print(f'PyTorch: {torch.__version__}, CUDA: {torch.cuda.is_available()}')"

# 2. 检查 PyTorch Geometric
python -c "import torch_geometric; print(f'PyG: {torch_geometric.__version__}')"

# 3. 检查关键依赖项
python -c "import esm, rdkit, Bio; print('ESM, RDKit, BioPython: OK')"

# 4. 验证示例数据存在
ls data/toy_example/
# 应显示：Enzymes.csv, Substrates.csv, data.csv, Structure/
```

## 架构

### 项目结构概览

```
EZSpecificity_Project/
├── src/                           # 源代码
│   ├── Models/                    # 模型实现
│   │   ├── ss.py                  # 主要 SS 模型（交叉注意力）
│   │   └── Structure/             # 基于结构的模型（EGNN）
│   ├── Datasets/                  # 数据处理
│   │   ├── brenda.py              # 数据集加载和批处理
│   │   ├── create_features.py     # 特征生成（ESM、图）
│   │   └── Structure/             # 结构数据处理
│   ├── other_softwares/
│   │   └── grover_software/       # GROVER 分子指纹
│   ├── example.ipynb              # 完整工作流程演示
│   └── main_testing.py            # 推理脚本
├── data/                          # 数据文件
│   ├── pretrain_model/            # 预训练权重
│   │   └── grover_large.pt
│   └── toy_example/               # 示例输入/输出
└── saved_model/                   # 训练好的模型检查点
    └── model/run_0/
        └── complete-full-random-all-0-complex.yml
```

### 核心模型流程

```
输入：酶序列 + 底物 SMILES + 对接复合物 PDB
         ↓
┌────────────────────────────────────────────────────────────┐
│  特征提取                                                    │
│  ├── ESM 嵌入（protein_mlp → 1280 → hidden_dim）            │
│  ├── GROVER 嵌入（原子级 + 分子级）                          │
│  ├── Morgan 指纹（1024-bit）                                │
│  └── 结构特征（口袋 + 配体坐标）                             │
└────────────────────────────────────────────────────────────┘
         ↓
┌────────────────────────────────────────────────────────────┐
│  SE(3)-等变 GNN（Models/Structure/egnn.py）                 │
│  └── EnBaseLayer × num_layers（图上消息传递）               │
└────────────────────────────────────────────────────────────┘
         ↓
┌────────────────────────────────────────────────────────────┐
│  交叉注意力（Models/ss.py）                                  │
│  ├── enzyme_attention：酶查询底物                            │
│  └── reaction_attention：底物查询酶                          │
└────────────────────────────────────────────────────────────┘
         ↓
┌────────────────────────────────────────────────────────────┐
│  特异性预测头（MLP）                                         │
│  └── 聚合嵌入 → 二分类                                       │
└────────────────────────────────────────────────────────────┘
```

### 关键源文件

| 文件 | 用途 |
|------|------|
| [src/Models/ss.py](src/Models/ss.py) | 主模型（SS 类）- 交叉注意力 + 预测头 |
| [src/Models/Structure/egnn.py](src/Models/Structure/egnn.py) | SE(3)-等变 GNN 实现 |
| [src/Models/Structure/structure.py](src/Models/Structure/structure.py) | 从结构构建图 |
| [src/Datasets/brenda.py](src/Datasets/brenda.py) | 数据集加载和批处理 |
| [src/Datasets/create_features.py](src/Datasets/create_features.py) | ESM 嵌入、反应特征 |
| [src/Datasets/Structure/protein_ligand.py](src/Datasets/Structure/protein_ligand.py) | 口袋/配体解析 |
| [src/main_testing.py](src/main_testing.py) | 推理脚本 |
| [src/example.ipynb](src/example.ipynb) | 端到端工作流程示例 |

### 推理数据流程

**两个主要阶段**：1）特征准备，2）预测生成

#### 阶段 1：所需输入文件

- `Enzymes.csv`：必须包含列 "Protein sequence"
- `Substrates.csv`：必须包含列 "Substrate_SMILES"
- `data.csv`：必须包含列 "Dock Index"、"Enzyme Index"、"Substrate Index"、"Label"
  - "Dock Index" d、"Enzyme Index" e、"Substrate Index" s 表示对接 d 是酶 e 和底物 s 的复合物
  - 索引必须与词汇文件（Enzymes.csv、Substrates.csv）对齐
  - "Label" 列可以随机设置为 0 或 1（仅用于流程执行，对预测无意义）
- `Structure/complex/*.pdb`：对接结果，需符合格式要求（见下方 PDB 部分）

#### 阶段 2：特征生成顺序

遵循 `example.ipynb` 工作流程：

**阶段 1：酶特征**
- ESM 嵌入 → `enzyme_features.lmdb`（蛋白质 < 1000 个氨基酸，仅标准氨基酸）

**阶段 2：底物特征**
- 图特征 → `reaction_features.lmdb`（底物 < 275 个原子）
- Morgan 指纹 → `morgan_fingerprint.npy`（1024-bit）
- GROVER 嵌入 → `grover_fingerprint.lmdb`（需要多个子步骤，见 GROVER 部分）

**阶段 3：结构特征**
- 从 PDB 文件提取口袋和配体 → `str_tmp_data/pocket/`、`str_tmp_data/ligand/`
- 对齐对接和 SMILES 表示之间的底物编号
- 处理结构特征 → `structure_features.lmdb`、`str_features.lmdb`
- 生成高质量 ID → `high_quality_id.txt`（有效口袋和配体的交集）

#### 阶段 3：预测输出

包含 "score" 列的 CSV（表示特异性的 logits）

### PDB 对接文件格式要求

对接 PDB 文件必须遵循以下严格格式要求：

- **辅因子处理**：任何辅因子相关原子必须以 "HETATM" 标识符开头（辅因子是可选的）
- **原子排序**：蛋白质/辅因子原子必须在文件中位于配体原子之前
- **分隔符要求**：蛋白质/辅因子和配体部分必须通过以下内容分隔：`COMPND    ***.pdbqt`
- **标准合规性**：
  - 蛋白质部分必须符合标准 PDB 格式
  - 配体部分必须可被 RDKit 包解析
- **坐标提取**：口袋定义为配体原子 10Å 范围内的蛋白质原子（不包括氢）

## 配置

模型配置基于 YAML（参见 `saved_model/model/run_0/complete-full-random-all-0-complex.yml`）：

```yaml
model:
  hidden_dim: 128
  use_gnn: True
  graph:
    num_layers: 3
    attention: True
  cross_attention:
    n_head: 8
    dropout: 0.9
```

**关键超参数**：
- `max_enzyme_length`：1450 个氨基酸
- `max_substrate_length`：280 个原子
- `batch_size`：16
- 特征：Morgan（1024-bit）+ GROVER mean（4885-dim）+ GROVER atom（2400-dim）

## 依赖项

**核心依赖项**：
- Python 3.10.12
- PyTorch 2.1.0 with CUDA 12.1
- PyTorch Geometric（torch_scatter、torch_sparse、torch_cluster）
- PyTorch Lightning 1.9.0
- RDKit 2023.9.2（分子处理）
- BioPython 1.81（PDB 解析）
- fair-esm 2.0.0（蛋白质嵌入）
- LMDB 1.4.1（特征存储）

**完整依赖列表**：参见 `src/environment.yml` 获取完整的 conda 环境规范。

## 研究背景

本项目基于 EZSpecificity（Nature 2025）原始代码库，当前正在进行 P450 酶家族的特异性研究。

### EZSpecificity 原始数据集
- 训练数据：来自 BRENDA 的 323,783 个酶-底物对
- 酶种类：8,124 种不同的酶
- 负样本生成：单向随机采样（仅采样负酶），正负比例 1:9.36

### 当前 P450 研究方向

本项目正在进行 P450 细胞色素酶家族的专项研究：

#### 研究文档索引

| 文档 | 内容 | 用途 |
|------|------|------|
| [P450_EZSpecificity完整研究手册_终极整合版.md](毕业设计/P450_EZSpecificity完整研究手册_终极整合版.md) | ESIBank 构建原理、研究路径分析、核心概念详解 | 总纲性文档，理解全局 |
| [提取P450过程日志/项目进度日志.md](毕业设计/提取P450过程日志/项目进度日志.md) | 动态更新的工作日志、版本化的工作节点 | 追踪进度、回溯历史 |

#### 关键数据资产

```
P450 数据层级结构：

ESIBank 训练集中的 P450
├── 389 个已验证的 P450 酶（从 25,225 个酶中筛选）
│   ├── 来源：UniProt + InterPro 双重验证
│   ├── 结构：大多使用 AlphaFold 预测结构
│   └── 用途：EZSpecificity 原始模型的训练数据子集
│
独立测试集（本研究重点）
└── 682 条酶-配体对（153个唯一P450，627个PDB）
    ├── 来源：RCSB PDB 数据库（实验结构）
    ├── 质量：X-ray/Cryo-EM 实验解析的高质量结构
    ├── 特点：模型从未见过的独立数据
    ├── 数据净化：已排除58条Photosystem I污染（740→682条）
    ├── 配体分类：已完成四类分类（SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE）
    └── 用途：验证模型在P450上的预测能力
```

#### 当前研究阶段

**阶段五十**：C3 底物分类 ✅ 完成 — 模型训练待启动

| 路径 | 名称 | 状态 | 说明 |
|------|------|------|------|
| **A** | 模型评估测试集构建 | ✅ 已完成 | AUC-ROC 0.6636 |
| **B** | P450数据集构建与结构优化 | ✅ 已完成 | Step 1-10, legacy_bug test AUC=0.7244 |
| **C** | P450专属模型训练 | 🔄 C3模型训练待启动 | |
| │    | └─ C1 论文基线训练与参数调整 | ✅ | AllSplit Test=0.7244(超论文) + fixed/dropout消融 |
| │    | └─ C2 P450全面数据集构建 | ✅ | **EXP001 Test AUC=0.7730** (4×4090, 89ep) |
| │    | └─ **C3 P450专属模型训练** | **底物分类✅** | 7+1类多标签分类(96%准确率, 2125化合物, ~40 agents验证) → 模型训练待启动 |
| D | 区域选择性预测 | ⏳ 待定 | S3反应SMILES + S9 PCPD反应图片 |

**工作目录**: `毕业设计/P450_EZSpecificity_研究项目/PathC_2026-03-19_P450专属模型训练/`

**Step 9 最终结果（AllSplit训练，已完成）**:
- **核心指标: AUC-ROC**（用户明确要求，AUPR仅辅助参考）
- **FINAL BEST: ep14 AUC=0.7667**（超论文 0.7198，+0.047）
- 训练到 ep27，过拟合确认：train_loss 0.107→0.088↓，val_loss 0.328→0.374↑
- AUC趋势: ep2(0.654) → ep8(0.742) → **ep14(0.7667 BEST)** → ep19(0.766) → ep27(0.751)
- EarlyStopping wait=13/15 at ep27, ep14 is definitive best
- 边排序Bug已在训练前修复（monkey-patch fixed模式），**baseline 已包含边修复**
- **数据泄露**: 酶0.15%（可忽略），底物6.56%（跨数据集Morgan FP重叠）
- 原始论文run_0: 4块GPU训了~256 epochs → AUC=0.893

**Step 10: .pt训练管线（已完成，2026-03-16）**:
- **缓存v3 (per-sample)**：176K个独立.pt文件(~160KB/个)，PyG标准模式
- 三层存储: enzymes.bin(27GB,seek读取) + substrates_grover.bin(4.3GB) + per-sample .pt(28GB)
- **速度: 7.56 it/s**（比LMDB缓存3.8→2.09 it/s衰减快2倍，不衰减）
- **内存: ~10GB**（vs LMDB 15-30GB mmap thrashing）
- 训练脚本: `scripts/10_Step10_pt训练管线/local/main_training_pt.py`（本地）/ `cloud2x4090/main_training_pt.py`（云服务器DDP）
- Dataset类: `scripts/10_Step10_pt训练管线/local/pt_dataset.py` / `cloud2x4090/pt_dataset.py`
- 支持 --edge-mode fixed/legacy_bug 切换（消融实验用）
- **目录组织（2026-03-19重构）**: scripts拆分为 `local/` 和 `cloud2x4090/` 子目录，数据缓存从 `ezspec_pt_v1` 重命名为 `allsplit_pt_cache`
- 迭代历史: v1(分片300MB,SSD爆)→v2(flatbin,慢)→**v3(per-sample,最优)**
- **为什么LMDB失败（关键洞察）**:
  - 论文作者: 256GB+ RAM服务器, LMDB 60GB mmap完全驻留内存 = 随机读=内存速度
  - 我们: 32GB RAM < 60GB LMDB → mmap thrashing(pageout/pagein循环) → 2.09 it/s decay, 99% page cache换页
  - .pt方案: 仅160KB/样本随机读 → ~10GB工作集 → no thrashing → 7.56 it/s stable, 3.6倍提速
- **P450数据集已完成（C2 Phase 7, 2026-03-26）**: 47,510可用对(1,514酶×2,117底物), 全部LMDB特征已生成
- **Phase 7.5 pt_cache构建完成（2026-03-26）**: LMDB→per-sample .pt转换
  - all_split: 8,326 train + 4,023 val + 4,143 test = 16,492 samples
  - random_split: 23,710 train + 11,878 val + 11,823 test = 47,411 samples
  - 三步流程: main(graph shards) → convert-per-sample → convert-flatbin(enzymes.bin+substrates)
  - 性能修复: --shard-size 2048→256 + OMP_NUM_THREADS=1, 从10+min→**30秒**
  - 目录: P450/data/pt_cache/{shared/,random/,all/}, shared通过symlink共享到各split
- **Phase 8 EXP001已完成（2026-03-26）**: P450 random_split fold0基线
  - **Val AUC=0.7544, Test AUC=0.7730**（89 epochs, early stopped ep88, best=ep73）
  - vs ESIBank P450 internal AUC=0.638 → **+0.135提升**，确认P450专属数据集构建成功
  - 服务器: Cloud-2 **4×RTX4090** (升级自2×4090), 360GB RAM, 28核(64 vCPU)
  - 配置: bs=56/GPU, effective=224, 4 GPU DDP, num-workers=6/GPU(24 total), --preload(~149GB RAM), edge-mode=fixed
  - 速度: ~1.8 it/s, ~2 min/epoch, GPU利用率92-100%, **总耗时49分钟**
  - main_training_pt.py模板重写（Codex 4轮review）: 输出路径扁平化、max-epochs=200、save_top_k=5、训练后自动测试、DDP test_epoch_end修复、自动关机、SRC_DIR fallback
  - 实验管理: 每个实验含完整src/副本实现架构隔离
  - 实验目录: P450/experiments/EXP001_random_baseline/
- **Step 10 legacy_bug 已完成**: Cloud-2 DDP 32ep, test AUC=0.7244 > paper 0.7198
- **Val Loss ↑ while AUC ↑（Codex深度分析，非Bug）**: BCE=逐点（outlier主导），AUC=成对排序（outlier鲁棒）。ep2→ep22: +26,500对正确排序但few dozen hard samples极端logit(z=5→loss=5.01)主导loss均值。三阶段: warmup(ep0-8)→divergence(ep8-22,排序↑过度自信↑)→true overfitting(ep22+)。ReduceLROnPlateau监控aupr（非auc/loss）不匹配。建议: 同AUC选低loss ckpt, warmup后加LR衰减, temperature scaling后校准。
- **EZSpecificity-individual**: 论文中"仅目标家族数据从头训练"模式（85-5424对），我们的AllSplit方式类似此模式
- **C1-Step 1 已完成**: Cloud-2 DDP 2×4090, --edge-mode fixed, bs=56, 32ep ~5.2h ~2.7 it/s. **Val AUC=0.7145 (ep16 best)**, early stopped ep31. 边修复未提升AUC（legacy_bug Val=0.722 vs fixed Val=0.7145, Δ=-0.008）. **Test AUC=0.7060**（AUPR=0.2362, 8841样本, vs legacy_bug test=0.7244, Δ=-0.018）. Log: PathC/logs/train_fixed.log
- **C1-Step 2 dropout消融已完成 (2026-03-21)**:
  - dropout=0.1: Val AUC=0.7216(+0.007), **Test AUC=0.6936(-0.012)**, AUPR=0.2012
  - dropout=0.3: Val AUC=0.7397(+0.025), **Test AUC=0.6959(-0.010)**, AUPR=0.2044
  - **关键发现: Val改善未迁移到Test**。Val-Test gap随dropout降低而增大(基线0.009→d=0.1: 0.028→d=0.3: 0.044)
  - 结论: dropout改动不纳入C1-Step 3组合，需关注能改善Test泛化的改动
  - 数据泄露验证: all_split模式0%酶/底物重叠（per-family），跨家族ID重用: 酶3.2%/底物1.6%
  - 代码整合: run_test_eval.py合并到main_training_pt.py --test-only, start_ablation.sh已删除
- **Cloud-2 服务器重构(2026-03-19→03-26)**: PathB/(归档) + PathC/P450/(自包含工作区) + PathC/ESIBank_baseline/(冻结) + allsplit_pt_cache/(ESIBank共享)
- **Cloud-2 服务器升级(2026-03-26)**: 2×RTX4090 → **4×RTX4090**, 360GB RAM, 28核(64 vCPU)

**Path B 诊断阶段总结（Step 7-8）**:
- **Step 7 核心发现**: 底物身份驱动评分（我们的P450）vs 配对交互信号（ESIBank P450）
- **Step 8 核心发现**: E8v1（EXP01）两个家族的通道重要性排序反转，但E8v2（Step 5）所有通道的AUC变化均不显著——瓶颈在任务/数据设计
- **假说验证（2026-03-09）**:
  - 假说1 ✅: ESIBank用EC层级混合难度0-5，我们≈难度3。ESIBank P450在难度≥3时AUC=0.549（vs我们0.517，差距仅0.032）
  - 假说2 ✅: 仅有random_split检查点，从头训不现实。微调方向确定
- Git 分支: `pathb-ablation`
- 下一步: 等待服务器（用户2026-03-12找导师要），在服务器上继续训练或Path C（P450微调）

**Path B Step 7 关键成果**:
- E4: R²(底物)=0.37, R²(酶)=-0.06, 残差AUC=0.509 → **底物身份驱动评分**
- E6 扩展: 7家族对比 — Esterase=0.934 ... Duf=0.796 ... **P450=0.517（最低）**
- E7(ESIBank内部基准): AUC=0.638，方向B占61%，R²(底物)=0.08，残差AUC=0.612
- **两种失败模式**: 我们=底物身份捷径；ESIBank=配对信号存在但P450家族本质困难
- 详见: `sessions/07_Step7_Tier1_诊断实验/session_log.md` 第四/七节

**Path B Step 8 关键成果**:
- **E8v1**（EXP01，晶体+抑制剂，495样本，基线AUC=0.712）：P450 关闭结构通道→AUC暴跌0.275，两个家族的通道重要性排序完全反过来
- **E8v2**（Step 5，Vina+随机负样本，2766样本，基线AUC=0.517）：**P450 关闭或保留任何通道，AUC变化均不显著！**
  - 关结构=-0.028, 关序列=+0.011, 关分子=+0.023，所有置信区间包含0
  - 磷酸酶对照组与E8v1完全一致（验证方法论正确）
  - **关键洞察**：E8v1的"结构通道最关键"结论只适用于EXP01条件，在真正崩溃的Step 5数据上，瓶颈在任务/数据设计层面
- E9 废弃：分子通道对P450仅贡献1.4%
- **修正后的Path C指导**: (1)优先修复任务/数据设计 (2)小规模试验对比"仅序列"/"序列+分子"/"序列+结构"三种配置 (3)不在分子通道上花力气
- 详见: `sessions/08_Step8_补充诊断实验/session_log.md`

**Path A 核心结论**:
- AUC-ROC: 0.6636（比论文Unknown enzyme+substrate场景低5.6%）
- 主因: 任务不匹配（抑制剂负样本 vs 随机配对负样本）+ 0%酶重叠
- 结论: 结果在预期范围内，不代表模型失败

#### 重要概念澄清

**配体（Ligand）vs 底物（Substrate）**：
- **PDB 中的配体**：晶体结构中与蛋白结合的任何小分子（可能是底物、产物、抑制剂、辅因子、溶剂等）
- **BRENDA 中的底物**：文献证实的酶催化反应的反应物
- **关键区别**：并非所有 PDB 配体都是酶的催化底物！需要查阅文献验证

**负样本生成机制**（EZSpecificity）：
- 方法：单向随机采样（只随机配对负酶，保持底物不变）
- 原因：避免采样到真实的正样本对（某个底物可能被多个酶催化）
- 比例：实际数据显示正负样本比为 1:9.36

详见各研究文档了解完整背景和技术细节。

## 常见问题与故障排除

### 特征生成问题

**蛋白质序列过长**：
- ESM 嵌入生成会跳过 > 1000 个氨基酸的蛋白质
- 解决方案：拆分长序列或使用替代嵌入方法

**非标准氨基酸**：
- 特征生成过程中会跳过含有非标准氨基酸的序列
- 解决方案：将非标准氨基酸替换为标准等效物或删除受影响的序列

**底物过大**：
- 图特征生成会跳过 > 275 个原子的底物
- 解决方案：分割大分子或调整 `max_substrate_length` 参数

**RDKit 解析失败**：
- PDB 文件中的配体部分必须与 RDKit 兼容
- 解决方案：检查键序和原子类型；确保正确的 PDB 格式

**parse_smile ValueError（已修复）**：
- `create_features.py` 中 `_smilesAtomOutputOrder` 解析使用 `[1:-2]` 会保留尾逗号导致 `int('')` 报错
- 已修复为 `[1:-1]` + 空字符串过滤（Cloud-2 `src/Datasets/create_features.py` 已打补丁）

### 结构处理问题

**对齐失败**：
- 如果对接 PDB 和 SMILES 之间的底物对齐失败，工具会尝试使用 SMILES 作为模板重新分配键序
- 解决方案：验证 SMILES 与对接配体结构匹配

**缺失原子**：
- 确保配体中的所有重原子都在蛋白质原子的合理距离内
- 解决方案：检查对接质量；必要时重新对接

**高质量 ID 过滤**：
- 最终预测仅使用口袋提取和配体对齐都成功的结构
- 解决方案：检查 `high_quality_id.txt` 查看哪些结构通过了过滤

### 环境问题

**Jupyter 中的 DDP 策略错误**：
- PyTorch Lightning 的 DDP 策略可能在交互式环境中失败
- 解决方案：使用 `main_testing.py` 脚本而不是 notebook 进行最终预测生成

**CUDA 版本不匹配**：
- PyTorch Geometric 包必须与您的 CUDA 版本匹配
- 解决方案：为您的 CUDA 安装正确的 PyG 版本（见安装命令）

**Windows 路径问题**：
- Unix 风格路径（/）可能在 Windows 上不起作用
- 解决方案：使用 Windows 风格路径（\）或原始字符串（r"path\to\file"）

### GROVER 特征生成

**PYTHONPATH 要求**：
- 运行脚本前必须将 PYTHONPATH 设置为 GROVER 目录
- 解决方案：`export PYTHONPATH=<project_root>/src/other_softwares/grover_software:$PYTHONPATH`（Linux/Mac）或 `set PYTHONPATH=<project_root>\src\other_softwares\grover_software;%PYTHONPATH%`（Windows）

**检查点位置**：
- 验证 `grover_large.pt` 路径与您的设置匹配
- 解决方案：如果缺失，从原始 GROVER 仓库下载；放置在 `data/pretrain_model/`

**多步骤流程**：
- GROVER 需要三个连续步骤（save_features → build_vocab → fingerprint）
- 解决方案：按 GROVER 部分中的确切顺序执行；不要跳过步骤

### 常见错误消息

**"CUDA out of memory"**：
- 在配置文件中减少批量大小
- 使用更小的模型或更少的层
- 分块处理数据

**"LMDB file not found"**：
- 特征文件尚未生成
- 解决方案：按正确顺序运行特征生成步骤（见推理数据流程）

**"Invalid PDB format"**：
- PDB 文件不符合格式要求
- 解决方案：检查 PDB 格式要求部分；确保正确的分隔符行

**"Sequence contains invalid characters"**：
- 蛋白质序列中含有非标准氨基酸
- 解决方案：仅使用标准 20 种氨基酸


