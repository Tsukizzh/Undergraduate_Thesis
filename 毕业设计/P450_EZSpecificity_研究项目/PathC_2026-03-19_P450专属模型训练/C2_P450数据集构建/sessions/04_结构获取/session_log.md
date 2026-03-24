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

对540个AlphaFill失败的酶（525真404 + 15无HEM），从AlphaFold DB下载裸结构。

**遇到的问题**：
1. URL版本错误：脚本用`model_v4.pdb`，但AlphaFold DB已升级到v6，v4全部404
2. 验证：`v4→404, v5→404, v6→200`
3. 修复后532个下载成功，8个AlphaFold DB也没有

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

## 5. 待处理

### 5.1 219个无结构酶

| 类型 | 数量 | 解决方案 |
|------|------|---------|
| needs_colabfold | 131 | ColabFold/ESMFold从序列预测结构 → 再移植血红素 |
| heme_transplant_failed | 54 | 上传到AlphaFill API处理，或用保守Cys锚定放置 |
| no_template | 26 | 同上 |
| alphafill_not_found | 8 | ColabFold预测 + 移植 |

### 5.2 后续Phase 6步骤

| 步骤 | 内容 | 状态 |
|------|------|------|
| 清理所有酶PDB | 去水/盐/非血红素配体，加氢 | 未开始 |
| 受体准备 | PDB → PDBQT，定义对接盒子(以HEM Fe为中心) | 未开始 |
| 底物准备 | SDF → PDBQT (Meeko) | 未开始 |
| Pilot对接 | 100-500对，CPU Vina验证流程 | 未开始 |
| 全量对接 | 52,254对，AutoDock-GPU，Cloud-2 | 未开始 |
