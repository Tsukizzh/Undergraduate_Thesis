# Phase 6: 结构获取与对接 Session Log

**日期**: 2026-03-24~25
**状态**: 结构获取大部分完成，对接未开始

---

## 1. 目标

EZSpecificity模型需要三种输入：酶序列、底物SMILES、酶-底物复合物3D结构。Phase 5产出了52,254行的data.csv，但`structure_index`全部为-1（没有3D结构）。Phase 6的目标是：
- 为1,622个酶获取含血红素的3D结构
- 为2,125个底物生成3D坐标
- 最终完成52,254个酶-底物对的分子对接，填充`structure_index`

---

## 2. 酶结构获取（1,622个酶）

### 2.1 总体数据流

```
1,622个P450酶
    │
    ├─ 349个已有PDB文件 ─────────────────────── 直接可用 ✅
    │   (103 S1实验 + 246 PCPD预测)
    │
    └─ 1,273个没有PDB文件
        │
        ├─ 1,165个有UniProt（后修正为1,141个真UniProt）
        │   │
        │   ├─ AlphaFill下载 ─── 602个成功
        │   │                      │
        │   │                      ├─ 493个 RMSD<5Å ─── 直接可用 ✅
        │   │                      │
        │   │                      └─ 109个 RMSD≥5Å ─── v2重新移植
        │   │                                              ├─ 89个成功 ✅
        │   │                                              └─ 20个失败 ❌
        │   │
        │   └─ AlphaFill失败 ─── 563个
        │       │
        │       ├─ AlphaFold下载裸结构 ─── 532个成功
        │       │   │
        │       │   └─ v2移植血红素
        │       │       ├─ 472个成功 ✅
        │       │       └─ 60个失败 ❌
        │       │
        │       └─ AlphaFold也没有 ─── 8个 ❌
        │
        └─ 131个无真正UniProt ──── 需ColabFold预测 ❌
            (108原始无UniProt + 24假UniProt修正)
```

### 2.2 最终结果

| 来源 | 数量 | 占比 | 说明 |
|------|------|------|------|
| S1实验PDB | 103 | 6.4% | RCSB晶体结构，质量最高 |
| PCPD预测PDB | 246 | 15.2% | Maestro生成，全部含HEM |
| AlphaFill下载（过滤后） | 493 | 30.4% | 从AlphaFill网站下载，全部RMSD<5Å |
| v2血红素移植 | 561 | 34.6% | AlphaFold裸结构 + 序列比对移植 |
| **有结构总计** | **1,403** | **86.5%** | |
| 无结构 | 219 | 13.5% | 131 ColabFold + 54移植失败 + 26无模板 + 8无AF |

### 2.3 详细过程

#### 步骤1：盘点已有资源

在Phase 6开始前，我们手上已有两批PDB文件：
- **103个S1_RCSB实验结构**：Path A阶段从RCSB获取的晶体结构，有真实PDB编号
- **246个PCPD预测结构**：从PCPD数据库下载的Maestro生成结构，检查确认全部423个文件都含血红素(HEM)，其中246个匹配到全局酶表
- PCPD文件分类：290个Maestro标注 + 133个无头信息，0个实验结构

#### 步骤2：AlphaFill批量下载

对剩余1,273个酶中有UniProt ID的1,165个，从AlphaFill网站批量下载含血红素的预测结构。

- 脚本：`scripts/04_结构获取/download_alphafill.py`
- API：`https://alphafill.eu/v1/aff/{UniProt_ID}-F1/json`（元数据）+ `https://alphafill.eu/v1/aff/{UniProt_ID}-F1`（结构）
- 结果：602个成功含HEM，548个404，15个无HEM

**遇到的问题**：
1. `data.get('hits', [])` 返回None导致TypeError → 修复为 `data.get('hits') or []`
2. AlphaFill只覆盖了53%的UniProt候选（AlphaFill是模板驱动的，不是所有蛋白都有同源模板）

#### 步骤3：发现并修复SEQHASH假UniProt Bug

审查发现24个酶的`canonical_uniprot_id`是内部标识符（23个`SEQHASH_...` + 1个GenBank ID `KAJ5738141`），不是真正的UniProt。

- **根因**：Phase 4 merge_dedup.py缺少UniProt格式验证，序列hash fallback合并时未清理
- **影响**：这24个被错误发送到AlphaFill API（当然404）
- **修复**：重新分类为needs_colabfold，真正UniProt从1,457→1,433
- **验证**：用正则 `^[A-Z][A-Z0-9]{4,9}$` 验证所有UniProt ID

#### 步骤4：AlphaFold裸结构下载

对540个AlphaFill失败的酶（525真404 + 15无HEM），从AlphaFold DB下载裸结构（不含任何配体/辅因子）。

**关于AlphaFold DB版本体系**（容易混淆，这里澄清）：

AlphaFold DB有**两套版本号**，含义完全不同：

| 维度 | 版本 | 含义 |
|------|------|------|
| **AI预测模型** | AlphaFold 2 (Monomer v2.0) | 实际生成3D坐标的深度学习模型。AlphaFold DB自2021年建库至今（2025年v6），所有结构**始终由AlphaFold 2预测**，从未更换为AlphaFold 3 |
| **数据库文件版本** | model_v6 (2025-09-15) | 数据库的发布版本号。v4→v6对已有蛋白的坐标**完全相同**（原文："coordinates carried over from v4"），仅新增蛋白/更新元数据 |

数据库版本简史：v1(2021,365K) → v2(2021,995K) → v3(2022,214M大扩展) → v4(2022,修复4.4%坐标bug) → v5(跳过,未公开) → **v6(2025,241M,当前版)**

PDB文件中可验证：`TITLE ALPHAFOLD MONOMER V2.0 PREDICTION FOR ...`（所有641个文件一致）。文件HEADER日期（如`01-JUN-22`或`01-AUG-25`）是蛋白被收录进数据库的日期，不是预测日期。

**遇到的问题**：
1. URL版本错误：脚本用`model_v4.pdb`，但AlphaFold DB的v4文件已下线（仅保留最新版v6），v4全部404
2. 验证：`v4→404, v5→404, v6→200`
3. 修复URL为`model_v6.pdb`后532个下载成功，8个AlphaFold DB也没有（UniProt ID未被收录）

#### 步骤5：血红素移植v1（失败）

用脚本从951个含血红素结构中选模板，对齐后复制血红素坐标到532个裸结构中。

**v1结果极差**：519完成但78.6% RMSD≥5Å（不可用），平均RMSD=12.03Å

**Codex审查发现三个根本Bug**：
1. **残基编号对齐**：按residue number配对Cα原子，但不同P450编号系统不同 → 配对了完全不对应的原子
2. **k-mer模板选择**：k=3 overlap太粗糙，不是真正的序列一致性
3. **实验PDB模板被忽略**：路径解析只处理了AlphaFill CIF和PCPD PDB，103个S1实验结构没被加载

#### 步骤6：血红素移植v2（成功）

完全重写移植脚本，修复三个Bug：

- 脚本：`scripts/04_结构获取/transplant_heme_v2.py`
- **序列比对驱动**：BioPython PairwiseAligner做全局比对，找到真正对应的残基对
- **真正的序列一致性**：计算match/aligned_positions，最低25%阈值
- **三个模板来源全启用**：AlphaFill CIF + PCPD PDB + S1实验PDB
- **质量分级**：Tier1(<2Å) / Tier2(2-3Å) / Tier3(3-5Å) / Tier4(≥5Å拒绝换下一个模板)

**v2结果**：472/532成功，全部Tier1-3，0个Tier4
- Mean RMSD: 2.49Å（v1=12.03Å）
- Median: 2.27Å（v1=9.59Å）

#### 步骤7：修复AlphaFill中质量差的结构

用户发现AlphaFill下载的602个中有109个RMSD≥5Å。对这109个也用v2重新移植：
1. 下载这109个的AlphaFold裸结构
2. v2移植替换AlphaFill的坏结果
3. 89个成功（RMSD<5Å），20个失败

修复后AlphaFill剩余493个全部RMSD<5Å。

### 2.4 v2移植质量统计

| 指标 | AlphaFill过滤后(493) | v2移植(561) |
|------|---------------------|------------|
| Mean RMSD | 2.52Å | 2.55Å |
| Median RMSD | — | 2.35Å |
| Tier1 (<2Å) | — | 188 (33.5%) |
| Tier2 (2-3Å) | — | 203 (36.2%) |
| Tier3 (3-5Å) | — | 170 (30.3%) |
| Tier4 (≥5Å) | 0 | 0 |

---

## 3. 底物3D坐标生成（2,125个底物）

### 3.1 流程

SMILES → RDKit解析 → 加氢 → ETKDG生成50个3D构象 → MMFF94s能量最小化（失败则UFF fallback）→ 选能量最低构象 → 保存SDF

### 3.2 遇到的问题

| 问题 | 原因 | 解决方法 |
|------|------|---------|
| 第1次运行：全部2,125个失败 | 默认Python没有RDKit | 换用torch conda环境的Python |
| 第2次运行：SDF写入失败 | RDKit SDWriter不支持中文路径 | 先写到临时文件（英文路径），再移动到目标 |
| 第3次运行：913个没被重处理 | 第1次的error记录残留在manifest中 | 清理manifest中的error行，重新运行 |

### 3.3 结果

| 指标 | 值 |
|------|---|
| 成功 | 2,124 (99.95%) |
| 失败 | 1 (embed_failed: 分子太复杂无法生成3D) |
| MMFF94s力场 | 2,107个 |
| UFF fallback | 17个 |

---

## 4. 产出文件

```
data/structures/
├── manifests/
│   ├── receptor_manifest.csv          ← 1,622个酶的结构状态+质量指标
│   ├── ligand_manifest.csv            ← 2,125个底物的3D生成状态
│   ├── structure_coverage_final.png   ← 覆盖率+质量分析图
│   └── phase6_funnel.png              ← 数据流漏斗图
├── alphafill/
│   ├── cif/                           ← 602个AlphaFill下载的CIF
│   └── json/                          ← 602个AlphaFill元数据JSON
├── alphafold/pdb/                     ← 641个AlphaFold裸结构PDB(532+109)
├── heme_transplant/pdb/               ← 561个v2移植后的PDB
└── ligands/sdf/                       ← 2,124个底物SDF文件

scripts/04_结构获取/
├── download_alphafill.py              ← AlphaFill批量下载
├── download_alphafold_and_transplant_heme.py ← v1(有Bug,已废弃)
├── transplant_heme_v2.py              ← v2血红素移植(序列比对驱动)
└── generate_substrate_3d.py           ← 底物3D生成
```

---

## 5. AlphaFold DB 版本澄清

下载的641个AlphaFold结构涉及两套容易混淆的版本号：

### 5.1 AI预测模型 vs 数据库文件版本

| 维度 | 版本 | 含义 |
|------|------|------|
| **AI预测模型** | AlphaFold 2 (Monomer v2.0) | 实际生成3D坐标的深度学习模型。AlphaFold DB自2021年建库至今，所有结构**始终由AlphaFold 2预测**，从未更换为AlphaFold 3 |
| **数据库文件版本** | model_v6 (2025-09-15) | 数据库的发布版本号。v4→v6对已有蛋白的坐标**完全相同**（CHANGELOG原文："coordinates carried over from v4"），仅新增蛋白/更新元数据 |

### 5.2 数据库版本历史

| DB版本 | 发布日期 | 结构数 | 主要变化 |
|--------|---------|--------|---------|
| v1 | 2021-07-22 | 365K | 初始发布，21个模式生物 |
| v2 | 2021-12-09 | 995K | 加入Swiss-Prot，TITLE改名为"Monomer v2.0" |
| v3 | 2022-07-28 | 214M | 大扩展：收录UniProt全部序列 |
| v4 | 2022-11-01 | 214M | Bug修复：4.4%结构因数值bug重新计算坐标 |
| v5 | — | — | 跳过（内部版本，从未公开发布） |
| v6 | 2025-09-15 | 241M | 同步UniProt 2025_03，新增26M+，首次加入isoform |

### 5.3 验证证据

每个PDB文件头部可确认：`TITLE ALPHAFOLD MONOMER V2.0 PREDICTION FOR ...`（全部641个文件一致）。文件HEADER日期（如`01-JUN-22`或`01-AUG-25`）是该蛋白被收录进数据库的日期，不是预测日期。

AlphaFold 3虽然2024年5月已发布、2025年2月开源，但AlphaFold DB至今没有用AF3重新预测任何结构。AF3仅通过AlphaFold Server（alphafoldserver.com）或本地安装使用。

---

## 6. 导师服务器与结构预测方案分析（2026-03-25）

### 6.1 导师服务器（wanglab jqlab4090）环境

**访问**: `dc@10.45.246.249`（Zerotier VPN），port 22，password: dachuang111

**硬件**:

| 项目 | 详情 |
|------|------|
| 主机名 | jqlab4090 |
| OS | Ubuntu 22.04, Kernel 6.5.0 |
| CPU | Intel i7-10700 @ 2.90GHz, 8核16线程 |
| RAM | **16GB**（较紧张，AF2预测500残基约需10-16GB） |
| GPU | RTX 4090 24GB, CUDA 12.2, 驱动535.154.05 |
| 系统盘 | 442GB（剩余257GB） |
| Home盘 | 1.8TB（剩余466GB） |

**已安装的结构预测工具**:

| 工具 | 位置 | 说明 |
|------|------|------|
| **LocalColabFold** (AF2) | `~/localcolabfold/` (8.5GB) | 已测试通过（P61823蛋白，5个模型成功生成） |
| **ESMFold** | `~/miniconda3/envs/esmfold/` | fair-esm 2.0.1 + torch 1.13.1+cu117 |
| **OpenFold** | `~/esm-repository/openfold-main/` | 完整安装 |
| Boltz-1/2 | ❌ 未安装 | |
| AlphaFold 3 | ❌ 未安装 | |

用户工作目录: `~/zhuangzeheng/EZSpecificity/`（已创建data/scripts/src/logs/results子目录）

### 6.2 导师提供的两个方案

1. **用AlphaFold 3预测结构**
2. **用导师服务器上本地部署的AlphaFold 2 + 导师提供的CYP81刚性结构模板**

### 6.3 导师的CYP81刚性结构模板——详细说明

#### 背景

导师课题组之前研究**CYP81**（一种植物P450酶，478个氨基酸残基），催化**薯蓣皂苷元（diosgenin，一种甾体化合物）**的反应。这就是导师说的"甾体研究"。

#### 什么是MD模拟

MD = 分子动力学模拟（Molecular Dynamics Simulation）。AlphaFold预测的结构是一张静态照片，但真实蛋白在细胞中泡在水里，会不停抖动。MD模拟就是用电脑模拟蛋白质在水中抖动的全过程：

1. 把蛋白放进虚拟水盒子（盒子大小89.3×96.7×86.8 Å）
2. 用物理公式（牛顿力学+分子力场）计算每个原子的受力和运动
3. 每0.2纳秒保存一帧，连续模拟500纳秒
4. 总共得到**2,502帧**快照

使用软件：**GROMACS 2023.5**

#### CYP81 WT（Root, 500ns）模拟的量化指标

**蛋白整体稳定性（RMSD）**:

| 指标 | 值 | 含义 |
|------|-----|------|
| 蛋白骨架RMSD平均值 | **4.2 Å** | 蛋白在模拟中和起始结构的平均偏移 |
| 蛋白骨架RMSD最大值 | 4.75 Å | 偏移最大的一帧 |
| 底物(diosgenin)RMSD平均值 | **2.8 Å** | 底物在口袋中波动较小→结合稳定 |
| 底物RMSD最大值 | 5.0 Å | 底物偶尔稍有晃动 |

蛋白RMSD在4.2Å左右趋于稳定→**已达平衡态**。

**每个残基的柔性（RMSF）**:

RMSF衡量每个残基在整个模拟过程中抖动幅度。477个残基分布：

| 柔性类别 | RMSF范围 | 残基数 | 占比 | 含义 |
|---------|-----------|-------|------|------|
| **刚性** | < 1.0 Å | 118 | 24.7% | 几乎不动，活性口袋核心、α螺旋内部 |
| **中等** | 1.0 - 2.0 Å | 258 | 54.1% | 有一定柔性但不剧烈 |
| **柔性** | > 2.0 Å | 101 | 21.2% | 大幅晃动，蛋白表面的环和末端 |

- RMSF最小值: 0.48Å（最刚性的残基）
- RMSF最大值: 5.05Å（最柔性的残基）
- RMSF平均值: 1.52Å

**导师说的"刚性结构"，就是指那24.7%（118个）RMSF < 1.0Å的残基区域**——活性口袋核心残基通常在此类别。

**铁原子与底物的距离（Fe-DIO C2）**:

P450催化核心是血红素铁原子（Fe），导师监测了Fe到底物diosgenin的C2碳的距离变化（2,502帧）：

| 指标 | 值 | 含义 |
|------|-----|------|
| 平均距离 | **5.0 Å** | 底物稳定待在铁原子附近 |
| 最小距离 | 4.0 Å | 最近时接近催化距离(P450催化约3-5Å) |
| 最大距离 | 7.1 Å | 底物偶尔稍远离 |

底物始终在活性口袋内→**结合模式稳定**。

**氢键分析（关键残基）**:

| 供体 | 受体 | 存在概率 |
|------|------|---------|
| THR333 N → HIE329 O | 主链氢键 | **96.4%** |
| THR333 OG1 → HIE329 O | 侧链氢键 | **95.2%** |
| HIE329 NE2 → LEU322 O | His-Leu氢键 | **80.8%** |

这些高占有率氢键构成活性口袋的骨架结构，是"刚性"的物理基础。

**聚类分析**: 2,502帧聚为**6个聚类**。clusters1.pdb = 出现频率最高的稳定构象。

#### MD结构 vs AlphaFold结构

| 对比 | AlphaFold 2 预测 | MD模拟后的结构 |
|------|-----------------|---------------|
| 来源 | 算法从序列预测，一次性输出一个静态构象 | 虚拟水环境中模拟500ns后蛋白自然平衡的构象 |
| 活性口袋 | 可能不够准确（未考虑底物和血红素影响） | 在底物+血红素存在条件下口袋形状更真实 |
| 柔性信息 | 无（只有冻结姿势） | 有（知道哪里刚性、哪里柔性） |
| 底物结合 | 不含底物 | 含底物，且验证了结合稳定性(RMSD=2.8Å) |

#### 模板的局限性

CYP81 模板只是**一个P450亚家族**的结构。我们219个酶横跨几十个不同亚家族（CYP1A2、CYP2D6、CYP71D等），CYP81与这些远亲的序列一致性可能只有20-30%。

#### ColabFold不支持"只用活性口袋片段当模板"

ColabFold的`--custom-template-path`只接受完整蛋白PDB文件。AF2内部的模板处理流程是：
1. 模板整条链（478残基）与目标序列做全局序列比对
2. 提取匹配残基的3D坐标作为"结构提示"
3. 输入到neural network的template track

如果只给50个口袋残基，和目标500残基比对时只有10%能匹配，效果反而变差。AF2代码也不支持片段输入。

### 6.4 219个酶的三种结构预测方案对比

| 方案 | 命令 | 模板来源 | 血红素 | 耗时估算 |
|------|------|---------|--------|---------|
| AF2无模板 | `colabfold_batch input.fasta output/` | 无 | ❌需要移植 | ~11-36h |
| **AF2自动搜索模板** | `colabfold_batch input.fasta output/ --templates` | PDB数据库自动匹配（每个酶找最相似的已知结构） | ❌需要移植 | ~15-40h |
| AF2+CYP81模板 | `colabfold_batch ... --templates --custom-template-path ~/test/DIO_5_12/` | 固定用CYP81 | ❌需要移植 | ~15-40h |
| Boltz-1 (AF3替代) | `boltz predict input_dir/` | 内置 | ✅直接预测 | ~7-11h |
| AlphaFold Server (AF3) | 网页手动提交 | 内置 | ✅直接预测 | ~8天(30/天限制) |

**当前推荐**: AF2自动搜索模板（`--templates`），在导师服务器上用LocalColabFold运行。理由：
1. 每个酶自动匹配最佳模板（比固定CYP81更好）
2. 符合导师"用模板"的意图
3. 导师服务器已装好LocalColabFold，零安装成本
4. 只比不用模板多一个参数

**决策待定**: 用户尚未最终确认方案

### 6.5 可用结构预测工具速度对比

| 工具 | 每个P450(~500aa)耗时 | 含血红素 | 可批量 | 硬件要求 |
|------|---------------------|---------|--------|---------|
| AlphaFold Server (AF3) | 5-15分钟 | ✅ | ❌手动 | 无(Google服务器) |
| Boltz-1 on RTX 4090 Linux | 2-3分钟 | ✅ | ✅ | 24GB GPU |
| Boltz-1 on RTX 4070S 12GB Win | 5-15分钟(有triton问题) | ✅ | ✅ | 12GB GPU(Windows需workaround) |
| LocalColabFold on RTX 4090 | 3-10分钟 | ❌ | ✅ | 24GB GPU + 16GB+ RAM |
| ESMFold | 秒级 | ❌ | ✅ | 较低 |

---

## 7. ColabFold 实测（2026-03-25 02:29-02:41）

### 7.1 选定方案

**AF2 + 自动 PDB 模板搜索**（`colabfold_batch input.fasta output/ --templates`），在导师服务器（wanglab）上运行。

选择理由：
1. 每个酶自动从 PDB 数据库匹配最相似的已知结构当模板（比固定 CYP81 更好）
2. 导师服务器已装好 LocalColabFold，零安装成本
3. 符合导师"用模板"的意图

### 7.2 测速结果（3 个蛋白）

| 蛋白 | 残基数 | AF2推理耗时 | pLDDT | 自动找到的PDB模板 |
|------|--------|-----------|-------|-----------------|
| ENZ_G000044 | 495 | 68s | **80.3** | 6c93_A, 6ma6_A, 6ma8_A 等20个 |
| ENZ_G000059 | 640 | 75s | **78.5** | 5te8_A, 6ma6_A, 6oo9_A 等20个 |
| ENZ_G000014 | 653 | 81s | **87.7** | 5nik_J, 5ws4_B, 5lil_A 等20个 |

全部成功，质量良好（pLDDT 78-88）。

### 7.3 遇到的问题

1. **hhsearch 不在 PATH 中**：`--templates` 模式需要 hhsearch，但默认 PATH 没加 ColabFold 的 bin 目录
   - 修复：`export PATH=$HOME/localcolabfold/localcolabfold/colabfold-conda/bin:$PATH`

### 7.4 全量运行（219 个，02:41 启动）

- 命令：`colabfold_batch input.fasta output/ --templates --num-models 1 --num-seeds 1`
- 输入：`~/zhuangzeheng/EZSpecificity/data/colabfold_predict/input.fasta`（219 条序列）
- 输出：`~/zhuangzeheng/EZSpecificity/data/colabfold_predict/output/`
- 预计耗时：219 × ~2.3 分钟 ≈ **8-9 小时**（预计 ~11:00 完成）
- 预测的结构**不含血红素**，跑完后需要用 v2 脚本移植

---

## 8. 分批对接管线计划（Codex 审核，2026-03-25）

### 8.1 核心原则

- **`Dock_Index` 是唯一主键**，从头到尾不改编号。分批只影响执行顺序，不影响 ID
- **永远从登记表（CSV）查任务**，不从文件夹内容推断
- **底物准备只做一次**（底物不分批，所有酶共享同一套底物）

### 8.2 分批定义

52,254 个酶-底物对根据**酶的结构就绪状态**分为三类：

| 分类 | 条件 | 酶数 | 对数（估算） |
|------|------|------|------------|
| **Batch 1** | 酶已有结构+血红素（1,403 个酶） | 1,403 | ~47,000+ |
| **Batch 2** | 酶在 ColabFold 预测中（219 个酶） | 219 | ~5,000+ |
| **Blocked** | 酶/底物结构失败 | 极少 | 极少 |

### 8.3 详细执行步骤

```
阶段 A：建登记表（本地，立即）
├── A1. 建 pair_registry_master.csv
│    ├── 导入 data.csv 的 52,254 行
│    ├── 根据 receptor_manifest.csv 标记每个酶 → batch_1 / batch_2 / blocked
│    ├── 根据 ligand_manifest.csv 标记每个底物 → ready / failed
│    └── 验证：总行数=52,254，Dock_Index 唯一，无遗漏
│
├── A2. 拆分为执行队列
│    ├── pair_registry_batch1.csv（酶有结构 且 底物有3D）
│    ├── pair_registry_batch2.csv（酶待预测 且 底物有3D）
│    └── pair_registry_blocked.csv（酶失败 或 底物失败）

阶段 B：底物准备（本地，只做一次，与 Batch 1 并行）
├── B1. 2,124 个 SDF → PDBQT（Meeko mk_prepare_ligand.py）
├── B2. 记录到 ligand_prep_registry.csv
└── B3. 1 个失败底物标记，对应的所有对标记 blocked

阶段 C：Batch 1 酶准备（本地，现在开始）
├── C1. 1,403 个酶 PDB 清理
│    ├── 去水分子（HOH）
│    ├── 去盐离子（Na, Cl 等）
│    ├── 去非血红素配体
│    ├── 保留：蛋白链 + HEM
│    └── 加氢
├── C2. PDB → PDBQT（Meeko mk_prepare_receptor.py）
├── C3. 定义对接盒子：以 HEM 的 Fe 原子为中心，24×24×24 Å
├── C4. 记录到 receptor_batch_assignment.csv
└── C5. 生成 docking_jobs_batch1.csv
     └── 列：Dock_Index, Enzyme_Index, Substrate_Index, Receptor_PDBQT, Ligand_PDBQT, Box_Center

阶段 D：Batch 1 对接（Cloud-2，2×RTX4090）
├── D1. 安装 AutoDock-GPU
├── D2. 上传 Batch 1 的 PDBQT + job list
├── D3. 跑 AutoDock-GPU
│    ├── 输出文件名包含 Dock_Index：dock_{Dock_Index}_e{Enzyme_Index}_s{Substrate_Index}
│    └── 每对 100 runs, 2.5M evals/run
├── D4. 后处理：提取最优 pose → 保存为 complex PDB
└── D5. 记录到 docking_results_batch1.csv

阶段 E：ColabFold 完成后，Batch 2 酶准备（本地）
├── E1. 下载 219 个预测结构到本地
├── E2. v2 脚本移植血红素
├── E3. PDB 清理 + PDBQT 转换（同 C1-C3）
├── E4. 更新 receptor_batch_assignment.csv
└── E5. 生成 docking_jobs_batch2.csv

阶段 F：Batch 2 对接（Cloud-2）
├── F1. 上传 Batch 2 的 PDBQT + job list
├── F2. 跑 AutoDock-GPU
└── F3. 记录到 docking_results_batch2.csv

阶段 G：合并与验证
├── G1. 合并 docking_results_batch1 + batch2 → docking_results_merged.csv
├── G2. 按 Dock_Index 关联回 pair_registry_master.csv
├── G3. 最终检查：
│    ├── completed + failed + blocked = 52,254
│    ├── 无重复 Dock_Index
│    ├── 无遗漏 Dock_Index
│    └── 每个 completed 行有且仅有一个 complex PDB 路径
└── G4. 更新 data.csv 的 structure_index 列（填入 Dock_Index → complex PDB 映射）
```

### 8.4 文件组织

```
data/structures/
├── receptors_clean/              ← 清理后的 PDB（1,622个）
├── receptors_pdbqt/              ← 对接输入-受体（1,622个）
├── ligands_pdbqt/                ← 对接输入-底物（2,124个）
├── complexes/                    ← 对接输出-复合物 PDB（52,254个）
├── colabfold_predict/            ← ColabFold 预测（导师服务器上运行）
│   ├── input.fasta               ← 219 条序列
│   └── output/                   ← ColabFold 输出
├── manifests/                    ← 已有的 manifest 文件
└── ...（已有的 alphafill/, alphafold/, heme_transplant/, ligands/ 目录）

data/registries/                  ← 新建：对接管线登记表
├── pair_registry_master.csv      ← 主表（52,254行，标记 batch/状态）
├── pair_registry_batch1.csv      ← Batch 1 执行队列
├── pair_registry_batch2.csv      ← Batch 2 执行队列
├── pair_registry_blocked.csv     ← 被阻塞的对
├── receptor_batch_assignment.csv ← 酶的 batch 分配 + 准备状态
├── ligand_prep_registry.csv      ← 底物准备状态
├── docking_jobs_batch1.csv       ← Batch 1 对接任务清单
├── docking_jobs_batch2.csv       ← Batch 2 对接任务清单
├── docking_results_batch1.csv    ← Batch 1 对接结果
├── docking_results_batch2.csv    ← Batch 2 对接结果
└── docking_results_merged.csv    ← 最终合并结果
```

### 8.5 Codex 提醒的关键陷阱

1. **不要重新编号 Dock_Index** — Batch 1 和 Batch 2 共用 data.csv 的原始编号
2. **不要从文件夹内容推断任务** — 永远从登记表（CSV）查
3. **不要忘记那 1 个失败底物对应的所有对** — 标记为 blocked
4. **输出文件名必须包含 Dock_Index** — 否则合并时无法对应
5. **不要假设 4 个 split 的 data.csv 可以独立处理** — 必须用合并后的主登记表
6. **batch 2 不能覆盖 batch 1 的输出** — 因为文件名包含 Dock_Index，不会冲突

### 8.6 当前执行状态（截至 2026-03-25 04:00）

| 阶段 | 状态 | 说明 |
|------|------|------|
| A. 建登记表 | ✅ 完成 | 52,254行，Batch1=46,728 / Batch2=5,495 / Blocked=31 |
| B. 底物准备 | ✅ 完成 | 2,118成功（含11个救回），7个无法救回 |
| C. Batch 1 酶准备 | 🔄 pilot待跑 | 文献调研+Codex审核完成，待写pilot脚本 |
| D. Batch 1 对接 | ⏳ 待 Cloud-2 | |
| E. Batch 2 酶准备 | ⏳ 待 ColabFold | ColabFold 02:41启动，~15/219完成 |
| F. Batch 2 对接 | ⏳ 待 E 完成 | |
| G. 合并验证 | ⏳ 待 D+F 完成 | |

---

## 9. 阶段 A 执行记录：建登记表 ✅（2026-03-25 03:00）

### 9.1 脚本

`scripts/05_对接管线/build_pair_registry.py`，经 Codex 两轮审核。

### 9.2 数据映射过程

**酶映射**（Enzymes.csv → receptor_manifest.csv）需要两步 join：
- Enzymes.csv 的 `uniprots` 列包含**混合内容**：1,457 行是真 UniProt ID，165 行是 `ENZ_G` 开头的全局酶 ID
- Step 1：uniprots → manifest.canonical_uniprot_id（匹配 1,457 个）
- Step 2：uniprots(ENZ_G*) → manifest.global_enzyme_id（匹配 165 个）
- 结果：1,622/1,622 全部匹配 ✅

**底物映射**（Substrates.csv → ligand_manifest.csv）：
- ligand_manifest 有 2,126 行（比 Substrates.csv 多 1 行，因为有 1 行 NaN compound_id）
- 去掉 NaN 行后，通过 SMILES 精确匹配：2,125/2,125 全部匹配 ✅

### 9.3 Codex 审核修正

**第一轮**：Codex 发现分批定义有误——88 个酶（54 heme_transplant_failed + 26 no_template + 8 alphafill_not_found）被标记为 `blocked`，但实际上这 88 个也在 ColabFold 上跑着，应该归入 Batch 2。

修正：`BATCH2_STATUSES` 从仅 `needs_colabfold`(131) 扩展为全部 219 个无结构酶。

**修正前后对比**：

| 分类 | 修正前 | 修正后 |
|------|--------|--------|
| Batch 1 | 46,728 | 46,728（不变） |
| Batch 2 | 3,041 | **5,495** |
| Blocked | 2,485 | **31** |

### 9.4 验证检查（全部通过）

- 总行数 = 52,254 ✅
- dock_index 唯一且为 0-52,253 ✅
- enzyme_index 在 0-1,621 范围内 ✅
- substrate_index 在 0-2,124 范围内 ✅
- join 后无行数变化 ✅
- 所有对分配到且仅分配到一个 batch ✅
- 46,728 + 5,495 + 31 = 52,254 ✅

### 9.5 产出文件

```
data/registries/
├── enzyme_registry.csv              ← 1,622 行（酶 → batch 分配 + 结构状态）
├── substrate_registry.csv           ← 2,125 行（底物 → 3D 状态 + PDBQT 状态）
├── pair_registry_master.csv         ← 52,254 行（主表，含 batch/status 标记）
├── pair_registry_batch_1.csv        ← 46,728 行（Batch 1 执行队列）
├── pair_registry_batch_2.csv        ← 5,495 行（Batch 2 执行队列）
└── pair_registry_blocked.csv        ← 31 行（底物 3D 失败导致的阻塞对）
```

---

## 10. 阶段 B 执行记录：底物 SDF → PDBQT ✅（2026-03-25 03:15-03:30）

### 10.1 脚本

`scripts/05_对接管线/prepare_ligands_pdbqt.py`，8 workers 并行，ProcessPoolExecutor。

### 10.2 遇到的问题与修复

| 问题 | 原因 | 修复 |
|------|------|------|
| RDKit SDMolSupplier 无法读中文路径 | Windows 中文路径编码 | 先 copy 到英文临时目录，在临时目录中验证和转换 |
| Meeko mk_prepare_ligand 无法写中文路径 | 同上 | 在临时目录生成 PDBQT，再 copy 回目标目录 |
| Path.resolve() 损坏中文路径 | Python pathlib 的编码问题 | 用 Path.parent 代替 resolve() |

### 10.3 第一轮结果

| 状态 | 数量 | 说明 |
|------|------|------|
| 成功 | 2,107 | 直接转换成功 |
| 通配符失败 | 15 | 含 `*` dummy atom |
| 含金属失败 | 1 | CMP_G000474 = Fe(III)-heme b（辅因子，不是底物） |
| 电荷计算失败 | 1 | CMP_G001812 = SMILES 录入错误（`I[I-]I`） |

### 10.4 底物抢救（Codex 两轮审核）

对 17 个失败底物逐个分析来源和化学合理性：

**可救回（11 个）**：

| ID | 原始 SMILES | 修复方法 | 修复后 | Codex 审核 |
|---|---|---|---|---|
| CMP_G000001-4 | `*C1CCC2...`（甾体骨架） | 去掉 `*`，RDKit 自动加 H 封端 | 甾体母核 | ✅ 标记为 capped_surrogate |
| CMP_G000005 | `*CC(=O)[O-]` | 去 `*` → 乙酸根 | `CC(=O)[O-]` | ✅ 保留带电形式（Codex修正：不要中和为-OH） |
| CMP_G000006 | `*CCC(=O)[O-]` | 去 `*` → 丙酸根 | `CCC(=O)[O-]` | ✅ 同上 |
| CMP_G000007-8 | `*[C@H](NC(=O)...)` | 去 `*`，保留立体化学 | 氨基酸衍生物 | ✅ Codex修正：保留 `/C=C\` E/Z标记 |
| CMP_G000013-14 | `C*C(=O)[O-]`, `CC*C(=O)[O-]` | 去 `*` → 乙/丙酸根 | 同 G000005/6 | ✅ 与 G000005/6 去重但保留独立 registry 行 |
| CMP_G001812 | `I[I-]I`（录入错误） | PubChem CID 656506 查到正确 SMILES | glucobrassicin（葡萄糖苷） | ✅ 标记为 exact_corrected |

**无法救回（6 个）**：

| ID | 原因 | 处理 |
|---|---|---|
| CMP_G000009 | `*[H]` = 占位符，不是真实底物 | 排除（excluded_placeholder） |
| CMP_G000010-12 | 多通配符黄酮（6-8 个 `*`），无法确定取代基 | 排除（excluded_generic_multiwildcard） |
| CMP_G000474 | Fe(III)-heme b = 辅因子，不是底物 | 排除（excluded_cofactor） |
| CMP_G001211 | embed_failed（原始 3D 生成就失败的） | 排除（embed_failed） |

### 10.5 抢救执行

对 11 个可救回底物：
1. 用修正后的 SMILES 重新生成 3D 坐标（ETKDG + MMFF94s）
2. 保存 SDF（覆盖原来的坏 SDF）
3. Meeko 转 PDBQT
4. 11/11 全部成功 ✅

### 10.6 最终结果

| 状态 | 数量 |
|------|------|
| PDBQT 成功 | **2,118**（2,107 原始 + 11 救回） |
| 无法转换 | 7（6 排除 + 1 embed_failed） |
| **总计** | 2,125 |

### 10.7 产出文件

```
data/structures/ligands_pdbqt/     ← 2,118 个 .pdbqt 文件 + 2,118 个 .ok 标记文件
```

耗时：4.6 分钟（8 workers 并行）+ 抢救约 1 分钟

---

## 11. 阶段 C 调研记录：受体 PDB 清理 + PDBQT 准备（2026-03-25 03:30-04:00）

### 11.1 文献调研

搜索了 P450 分子对接受体准备的标准做法（ResearchGate 讨论、PMC 文献、Meeko/AutoDock 官方文档）。

**P450 对接中血红素的处理（文献共识）**：
- **School A**：去掉血红素，不参与对接 → 不推荐用于 P450（会扭曲结合口袋形状）
- **School B**：保留血红素作为受体的一部分 → **P450 推荐做法**
- Fe 原子需要**手动赋电荷**（AutoDock/Meeko 不自动处理金属离子），文献中通常赋 +2.000（Fe²⁺）

**参考文献**：
- [How to dock heme proteins with AutoDock (ResearchGate)](https://www.researchgate.net/post/How_to_do_docking_heme_proteins_with_AutoDock_and_MGL_Tools)
- [AutoDock suite tutorial (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4868550/)
- [Meeko receptor CLI docs](https://meeko.readthedocs.io/en/develop/rec_cli_options.html)
- [AutoDock Vina zinc metalloprotein tutorial](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html)

### 11.2 关键发现：EZSpecificity 需要 PDB 格式而非 PDBQT

检查 toy_example 中的对接复合物文件后发现，EZSpecificity 期望的是 **PDB 格式**的复合物文件，格式为：
1. 蛋白原子（ATOM 行，含 HETATM 的辅因子）
2. 分隔行：`COMPND    filename.pdbqt`
3. 配体原子（ATOM 行 + CONECT 键连接）

因此每个酶需要输出**两个文件**：
- `clean_receptor.pdb` → 最终与配体合并给 EZSpecificity
- `receptor.pdbqt` → 给 AutoDock-GPU 对接用

### 11.3 Codex 审核后的受体准备协议

**每个酶的处理流程**：
```
原始文件(PDB/CIF) → 格式统一(CIF先转PDB) → 清理(去水/盐/杂质，保留蛋白+HEM)
  → 验证(恰好1个HEM + 1个Fe) → 提取Fe坐标
  → 输出1: clean_receptor.pdb（给EZSpecificity）
  → 输出2: receptor.pdbqt（给AutoDock，需检查/修补Fe电荷）
  → 输出3: box参数（Fe为中心，24×24×24Å）
```

**清理规则**：
- 保留：蛋白原子 + 恰好 1 个 HEM 残基（含 Fe）
- 去除：水（HOH/WAT）、离子（Na/Cl/Mg 等）、非 HEM 配体、缓冲液分子、altloc 重复
- 多链处理：保留含血红素的链 + 与活性位点接触的链，去掉结晶伴侣/抗体/远端链

**4 种来源的特殊处理**：

| 来源 | 数量 | 格式 | 文件位置 | 特殊注意 |
|------|------|------|---------|---------|
| S1 实验 PDB | 103（47个PDB文件） | PDB | `downloads/RCSB/PDB/` | 清理负担最重：多条链、水、离子、共结晶配体、altloc |
| PCPD 预测 PDB | 246（423个文件） | PDB | `downloads/PCPD/PDB/` | 较简单，验证 HEM 完整性 |
| AlphaFill CIF | 493 | **mmCIF** | `data/structures/alphafill/cif/` | 需先用 BioPython MMCIFParser 转 PDB（Meeko 不支持 CIF） |
| v2 移植 PDB | 561 | PDB | `data/structures/heme_transplant/pdb/` | 最简单，仅含蛋白+HEM |

### 11.4 Codex 建议的 pilot 测试

**在全量跑之前，先测 20 个（每种来源 5 个）**，验证：
1. CIF → PDB 转换是否保留 HEM 和 Fe
2. 清理后是否恰好 1 个 HEM
3. Meeko 能否正确处理含 HEM 的受体
4. PDBQT 中 Fe 是否保留且有电荷
5. 对接盒子中心是否正确

### 11.5 Pilot 测试（5 轮迭代）

先尝试用 **Meeko mk_prepare_receptor** 做 PDBQT 转换，4 种来源各 1 个酶。

**第 1 轮**：直接跑
- PDB 清理：4/4 通过（修了 AlphaFill CIF 76 条链问题 + PCPD Fe 元素识别）
- Meeko PDBQT：仅 v2 通过（1/4），其他 3 种失败（HEM 模板匹配、元素 X、Se 原子）

**第 2-4 轮**：尝试修复 Meeko（normalize_pdb 函数、HEM 原子排序、split-and-merge 方案）
- 修了 MSE→MET、元素列修复、HEM 排序等问题
- 最终 Meeko split-and-merge 方案：v2 + PCPD + AlphaFill 通过（3/4），实验 PDB 仍失败

**第 5 轮（关键转折）**：查 Path B 文档发现**之前用的就是 MGLTools `prepare_receptor4.py`**，不是 Meeko！
- MGLTools 就装在 `D:\autodock\MGLTools-1.5.7\`
- 命令：`prepare_receptor4.py -r input.pdb -o output.pdbqt -A hydrogens -U nphs_lps_waters`
- 注意：不能用 `-U nonstdres`（会删掉 HEM）
- **4/4 全部通过！** HEM 和 Fe 全部保留

### 11.6 全量执行 ✅（2026-03-25 04:15-04:25）

脚本：`scripts/05_对接管线/batch_receptor_pdbqt.py`（6 workers 并行）

处理流程（每个酶）：
1. PDB 清理（BioPython 解析 → 保留蛋白+1个HEM → normalize_pdb 修复格式）
2. 复制到英文临时目录（MGLTools = Python 2.7，不支持中文路径）
3. 修复 altloc 列（col 17 → 空格，Path B 验证的 fix）
4. MGLTools `prepare_receptor4.py` 转换
5. 验证 Fe 在 PDBQT 中存在

**结果**：

| 状态 | 数量 | 说明 |
|------|------|------|
| 成功 | **1,397** | 99.6% |
| 失败 | 6 | 5 个 MGLTools 报错 + 1 个无 HEM |
| 总计 | 1,403 | |

耗时：**6.3 分钟**（6 workers 并行）

Fe 电荷在 PDBQT 中为 0.000（Gasteiger 不处理金属离子）。但后续确认使用 Vina 打分函数（忽略电荷），所以 Fe 电荷不影响对接结果。

产出：`data/structures/receptors_pdbqt/`（1,397 个 .pdbqt + .ok 标记）
产出：`data/structures/receptors_clean/`（1,397 个干净 PDB）

### 11.7 对接工具选型（Codex + 文献调研）

**决策：使用 Uni-Dock（Vina GPU 加速版）而非 AutoDock-GPU（AD4）**

| 对比 | AutoDock-GPU (AD4) | Uni-Dock (Vina) |
|------|-------------------|-----------------|
| 精度（姿态 RMSD） | 1.42 Å | **1.30 Å**（更好） |
| 成功率（RMSD<2Å） | 77% | **81%** |
| 速度 | 快 | **更快（~1000x vs CPU Vina）** |
| Fe 电荷 | 需要手动修补 | **不需要（Vina 忽略电荷）** |
| 额外准备 | 需要 AutoGrid 格点图 | **不需要** |
| 输入 | receptor.pdbqt + ligand.pdbqt + .maps.fld | **receptor.pdbqt + ligand.pdbqt** |

选择 Vina 的理由：
1. **姿态精度更高**（[JCIM 2020](https://pubs.acs.org/doi/10.1021/acs.jcim.9b00778)）
2. **操作更简单**：不需要 AutoGrid 格点图，不需要修补 Fe 电荷
3. **速度更快**：Uni-Dock GPU 比 CPU Vina 快 1000-2000 倍
4. **我们从头训练自己的模型**：不需要和论文的 AD4 姿态分布一致
5. **P450 柔性口袋问题**：[2024 ChemRxiv benchmark](https://chemrxiv.org/engage/chemrxiv/article-details/67470b577be152b1d00cfc8f) 表明 Vina 和 AD4 在 P450 上表现类似，Fe 电荷不是主要瓶颈

### 11.8 Cloud-2 设置（2026-03-25 04:30）

**Uni-Dock v1.1.3 已安装**（之前 conda env 已存在）

**目录结构**：
```
/root/rivermind-data/EZSpecificity/PathC/docking/
├── receptors_pdbqt/    ← 1,397 个受体（待传）
├── receptors_clean/    ← 1,397 个干净 PDB（待传）
├── ligands_pdbqt/      ← 2,118 个底物（待传）
├── registries/         ← 登记表（待传）
├── scripts/            ← 对接脚本
├── results/            ← 对接原始输出
├── complexes/          ← 最终复合物 PDB
└── logs/               ← 日志
```

**磁盘**：删除 enzymes.bin（25.5GB）后，可用空间 43GB，足够对接使用。

### 11.9 当前执行状态（2026-03-25 05:00）

| 阶段 | 状态 | 说明 |
|------|------|------|
| A. 建登记表 | ✅ | 52,254行，Batch1=46,728 / Batch2=5,495 / Blocked=31 |
| B. 底物 PDBQT | ✅ | 2,118/2,125 成功（含 11 个救回） |
| C. 受体 PDBQT | ✅ | 1,397/1,403 成功（MGLTools，6.3 分钟） |
| Cloud-2 设置 | ✅ | Uni-Dock 已装，目录已建，43GB 可用 |
| ColabFold | 🔄 | ~53/219 完成（老师服务器运行中） |
| D. 数据传输 | ⏳ | ~1GB 待传到 Cloud-2 |
| E. Batch 1 对接 | ⏳ | 46,728 对，Uni-Dock GPU |
| F. Batch 2 准备 | ⏳ | 待 ColabFold 完成 |
| G. Batch 2 对接 | ⏳ | |
| H. 复合物组装 | ⏳ | 对接结果 → EZSpecificity PDB 格式 |
| I. 合并验证 | ⏳ | completed + failed + blocked = 52,254 |
