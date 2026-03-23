# Session Log: Phase 1 数据收集

**日期**: 2026-03-22
**目标**: C2 P450全面数据集构建 — Phase 1 数据源下载与提取

---

## 1. 会话前状态

- Phase 0 已完成: ESIBank P450 12,329条确认可用, RCSB 682条确认完整
- C2_计划与进度.md 已写入8阶段计划
- 用户确认计划, 新增需求: **双向负样本生成** (固定底物换酶 + 固定酶换底物, 各5个/正样本)

---

## 2. S3_P450Rdb 下载与提取

### 2.1 P450Rdb 是什么

P450Rdb（P450 Reaction Database）是成都中医药大学团队维护的数据库，专门收集已发表文献中 P450 酶催化的化学反应。里面全部都是 P450 酶，没有其他类型的酶。

### 2.2 下载内容 (2026-03-22)

从 https://www.cellknowledge.com.cn/p450rdb_v2/download.html 下载4个文件，保存到 `downloads/P450Rdb/`：

| 文件 | 大小 | 内容 |
|------|------|------|
| Reactions.csv | 84MB | 反应记录，每行是一个P450催化的化学反应 |
| P450s.csv | 536KB | P450酶的基本信息（名称、物种、EC号等） |
| Compounds.csv | 352KB | 参与反应的化合物信息（底物和产物） |
| Sequences.fasta | 445KB | 858个P450酶的蛋白质序列 |

### 2.3 Reactions.csv 的结构（85列）

每一行代表一个P450催化的反应，信息非常丰富：

**酶的信息：**
- [0] ID — 反应编号（如 P450R0001）
- [1] P450 symbol — CYP命名（如 CYP726A20、CYP3A4）
- [2] Protein name, [3] Gene ID, [4] Uniprot ID
- [5] EC number, [6] Species（物种）, [7] Species classification（植物/动物/微生物）
- [8] Txid（NCBI分类ID）
- [80] blast — 完整蛋白序列（这是为什么文件有84MB）

**反应的信息：**
- [9-14] Substrate1（底物1：名称、PubChem CID/SID、分子式、CAS号、SMILES）
- [15-29] Substrate2-4（同样结构，最多4个底物）
- [30-66] Product1-7（最多7个产物，同样字段）
- [69] Solvents, [70] Conditions, [71] Transformations（反应类型：氧化/还原等）
- [81] reaction, [82] bond（键变化信息）

**文献来源：**
- [72] PMID, [74] DOI, [75] Title, [76] Journal, [77] Time
- [73] Patent ID（专利号）

**验证状态：**
- [84] Reaction existence — "Experimentally validated" 或 空

### 2.4 数据实际情况：两张表拼在一起

Reactions.csv 有 132,764 行，但实际是两张完全不同的表拼在一起：

**第一部分（3,753 行）—— 真正的反应数据：**
- 每行ID以 "P450R" 开头（P450R0001, P450R0002, ...）
- 标记为 "Experimentally validated"
- 有完整的酶-底物-产物-文献信息
- 858个P450酶，来自400+个物种
- 物种分布：微生物 46.5%、植物 27.1%、动物 26.4%
- 79.6% 有PubMed文献引用

**第二部分（第3866行起，129,006 行）—— UniProt蛋白+Rhea反应映射表：**

这是一张独立的表，有自己的表头（15列，后面填空到85列）：

| 列 | 含义 | 示例 |
|-----|------|------|
| [0] Uniprot ID | P450蛋白的UniProt ID | A0A061I7U6 |
| [1] Reviewed | UniProt审核状态 | reviewed / unreviewed |
| [2] Protein names | 蛋白名称 | Cytochrome P450 2E1 (EC 1.14.13.n7) |
| [3] Gene Names | 基因名 | H671_3g10614 |
| [4] GeneID | NCBI Gene ID | (空或数字) |
| [5] Organism | 物种 | Cricetulus griseus (Chinese hamster) |
| [6] Organism (ID) | NCBI Taxonomy ID | 10029 |
| [7] PubMed ID | 文献PMID | 23929341 |
| [8] Length | 蛋白长度（氨基酸数） | 580 |
| [9] Protein existence | 蛋白存在证据 | Inferred from homology |
| [10] **Rhea ID** | 反应数据库编号 | 50088 |
| [11] EC number | EC分类 | 1.14.13.n7; 1.14.14.1 |
| [12] Subcellular location | 亚细胞定位 | Endoplasmic reticulum membrane |
| [13] **Equation** | 反应方程式（文字） | docosahexaenoate + NADPH → ... |
| [14] Reaction existence | 反应验证状态 | 全部为 "Inferred" |

关键统计：
- 42,394 个唯一 UniProt ID（1,548 reviewed + 127,458 unreviewed）
- 128,758 行有 Rhea ID（但只有 **614 个唯一 Rhea ID**）
- 127,759 个唯一 (UniProt, RheaID) 配对
- 129,005 行有文字反应方程式
- **全部 129,006 行都标记为 "Inferred"**（计算推断，零实验验证）

**为什么不作为正样本使用：**
- 全部是 "Inferred"（基于序列同源性推断的酶-反应关联），不是实验验证
- 只有 614 个唯一 Rhea ID，意味着绝大多数行共享同一个反应模板（同一 EC 类的不同物种同源酶被分配了相同的反应）
- 没有底物 SMILES（只有文字方程式和 Rhea ID）
- 如果当作正样本训练，模型可能学到 EC/反应模板的规律，而不是酶-底物特异性

**已提取为独立文件：**
- `inferred_enzymes.csv`：42,394 个唯一 P450 UniProt ID（384 reviewed + 42,010 unreviewed），含蛋白名称、物种、EC号等
  - 用途：以后扩充负样本的酶池（需先验证确实是 P450）
- `rhea_templates.csv`：614 个唯一 Rhea 反应模板，含文字反应方程式和 EC 号
  - 用途：化学元数据参考，以后可通过 Rhea→ChEBI→SMILES 获取化合物结构

### 2.5 验证反应统计（3,753行）

| 指标 | 数值 |
|------|------|
| 验证行总数 | 3,753 |
| 有UniProt ID | 3,303 (88.0%) |
| 无UniProt ID | 450 (12.0%) |
| 有Substrate1 SMILES | 3,706 (98.7%) |
| 无Substrate1 SMILES | 47 (1.3%) |
| 有PMID文献引用 | 2,988 (79.6%) |
| 有Substrate2 SMILES | 1,733 (46.2%) |
| 有Substrate3 SMILES | 1,627 (43.4%) |
| 物种分布(有UniProt) | 微生物 1,538 / 植物 894 / 动物 871 |

### 2.6 提取过程中的5个问题及修复

#### 问题1: 129K行非验证数据
- **现象**: 132K行中只有3,753条标记为 "Experimentally validated"
- **原因**: 是UniProt蛋白目录（42K个ID）被拼接在反应表后面，不是损坏的反应数据
- **处理**: 只使用3,753条验证行，129K行保留在 `downloads/P450Rdb/` 原始文件中不动
- **以后可用**: 42K个UniProt ID可用于扩充负样本的P450酶池（需先验证是否确实都是P450）

#### 问题2: 450条验证反应没有UniProt ID
- **现象**: 有CYP符号（如 CYP81Q58）和物种名（如 Cucumis sativus），但UniProt ID为 "/"
- **原因**: CYP命名系统（Nelson系统）和UniProt的基因名不一样。UniProt存的是基因名如 "CYP3A4"，但P450Rdb可能用变体格式
- **处理**: 用UniProt REST API按 (gene_exact:{CYP符号} AND organism_name:"{物种}") 查询
- **结果**: 只匹配上2条（cyp102A1 → P14779，即著名的P450 BM3）。大部分CYP符号在UniProt中不是标准基因名

#### 问题3: 75行无UniProt但有蛋白序列
- **现象**: 问题2中解析失败的448行中，75行在Reactions.csv的blast列有完整蛋白序列
- **处理**: 用序列的SHA256哈希值（取前16位，加 "SEQHASH_" 前缀）作为临时酶ID。我们的模型需要蛋白序列做ESM编码，有序列就够了，不一定需要UniProt ID
- **结果**: 恢复75行。剩余373行既没有UniProt也没有序列，无法使用

#### 问题4: 9行有复合UniProt ID
- **现象**: 如 `P18326/P18327/P26911`（斜杠分隔）或 `A0A5B8ND22;P0DXH4;/`（分号分隔）
- **原因**: 多个P450同工酶都能催化同一个底物。如 cyp105A1/cyp105B1/cyp105D1 都能催化 Alpha-Ionone
- **处理**: 拆分为独立的正样本对。如 (P18326, Alpha-Ionone, +)、(P18327, Alpha-Ionone, +)、(P26911, Alpha-Ionone, +)
- **结果**: 9条拆成27条（新增18行）

#### 问题5: 7个UniProt ID不在FASTA文件中
- **现象**: 验证反应引用了7个UniProt ID，但下载的Sequences.fasta（858条）里没有
- **具体**: A0A2V1NWK9(TxtE), O24782(CYP152), P9WPP0/P9WPP1(cyp125A1, 结核分枝杆菌), Q4WAW5/Q4WAW8(ftmP450, 烟曲霉), Q9XHE8(CYP71D18, 留兰香)
- **处理**: 从UniProt REST API获取FASTA序列；A0A2V1NWK9从blast列恢复
- **结果**: 全部8个（含1个blast恢复）成功获取序列

#### 问题6: 47行没有底物SMILES
- **现象**: 只有底物名称（如 "5-Keto-Casbene"、"13-Hydroxyl-8-Abietene"），没有SMILES也没有PubChem CID
- **原因**: 这些是生物合成途径的中间产物，太专业，PubChem里没有标准条目
- **处理**:
  1. PubChem名称查询 → 0/47成功
  2. 按PMID聚类分析：47行分布在26篇论文中，最大的组只有6行（PMID 26936895, CYP720B2/CYP720B12的松香烷类中间产物）。太分散，手动查26篇论文恢复几十条SMILES的收益太低
  3. Pathway chaining检查：部分底物（如"5-Keto-Casbene"）是同一酶前一步反应的产物，但产物也没有SMILES，无法链式恢复
- **结果**: 这47行保留在reactions.csv元数据中（有底物名称、有产物名称、有PMID），不进入interactions.csv（因为没有化学结构，模型无法编码）

#### 问题7: 编码问题
- **现象**: Reactions.csv包含非UTF-8字符（0xa6位置）
- **处理**: Python读取时使用 `errors='replace'` 容错

### 2.7 最终提取结果

**修复前 vs 修复后对比：**

| 项目 | 首次提取（只看UniProt+SMILES） | 修复后（加入序列hash等） | 增量 |
|------|------|------|------|
| 酶 | 833 | **857** | +24 |
| 化合物 | 1,437 | **1,492** | +55 |
| 正样本对 | 2,710 | **2,798** | +88 |
| 反应记录 | 3,248 | **3,352** | +104 |

**输出文件（保存在 `data/sources/Source_P450Rdb/`）：**

来自验证反应（3,753行）：

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 857个酶 | 全部有蛋白序列。列: enzyme_id, uniprot_id, p450_symbol, species, species_class, ec_number, sequence, sequence_length |
| compounds.csv | 1,492个化合物 | 全部有SMILES。列: compound_id, smiles, name |
| interactions.csv | 2,798条正样本对 | 列: interaction_id, enzyme_id, compound_id, label, source, quality_tier, num_reactions, has_pmid, has_products |
| reactions.csv | 3,352条反应记录 | 含产物SMILES和PMID。列: reaction_id, enzyme_id, compound_id, substrate_smiles, product_smiles, pmids |
| unresolved.csv | 145个CYP | 未能解析UniProt的CYP符号，供以后参考 |

来自推断酶目录（129K行）：

| 文件 | 数量 | 说明 |
|------|------|------|
| inferred_enzymes.csv | 42,394个酶 | 384 reviewed + 42,010 unreviewed。列: uniprot_id, reviewed, protein_name, gene_names, organism, taxid, ec_number, ... |
| rhea_templates.csv | 614个反应模板 | Rhea反应数据库编号+文字方程式。列: rhea_id, equation, ec_number, reaction_existence |

**未能使用的数据：**
- 373行（约10%的验证反应）因既无UniProt也无蛋白序列而无法使用
- 47行因无底物SMILES而不进入interactions.csv（保留在reactions.csv元数据中）
- 145个CYP符号未能解析为UniProt ID

**提取脚本**: `scripts/01_数据下载/S3_P450Rdb_extract.py`

### 2.8 和ESIBank P450的对比

| | ESIBank P450 | P450Rdb |
|--|-------------|---------|
| 正样本 | 884 | 2,798 |
| 酶数量 | 367 | 857 |
| 底物数量 | 393 | 1,492 |
| 有产物信息 | 无 | 有（大部分有产物SMILES） |
| 有文献引用 | 无 | 有（79.6%有PMID） |
| 数据来源 | BRENDA数据库 | 文献手动整理 |
| 数据质量 | 数据库自动提取 | 人工审核+实验验证标记 |

两个数据源有一定重叠（都涉及P450酶-底物关系），具体重叠率需要合并时测量。

---

## 3. S2_ESIBank P450 子集提取

### 3.1 数据来源

已有提取文件: `PathA_.../source_data/02_底物数据/P450酶底物反应详表_完整版.csv`
提取过程记录: `毕业设计/提取P450过程日志/` (2025-12-31 ~ 2026-01-04)

原始数据: 12,329条（884正样本 + 11,445负样本）
- 正样本的 difficulty = -1，负样本 difficulty = 0-5
- enzyme/reaction 列是 ESIBank 主数据集的行索引
- 10列: uniprot_id, protein_families, ecnumber, substrate_smiles, label, has_heme, enzyme, reaction, difficulty, sequence

### 3.2 索引验证

P450文件中的 `enzyme` 和 `reaction` 列是 0-based 行索引，分别指向 G 盘 ESIBank 原始文件：
- `enzymes.csv`（25,225行）：enzyme 索引 → 该行的 UniProt ID 和序列
- `reaction.csv`（39,635行）：reaction 索引 → 该行的完整反应 SMILES 和底物 SMILES

**全量验证结果（880条有效正样本）：**
- 酶索引不匹配: 0
- 反应索引不匹配: 0
- 结论: 索引完全正确，无 off-by-one 问题

### 3.3 辅因子过滤（Codex审查发现的问题）

首次提取时只排除了 4 条（SMILES = "O" 和 "OO"），但 Codex 审查发现 `O=O`（分子氧）也是辅因子，被遗漏了。

**分子氧 (O=O) 的问题：**
- 所有 P450 反应都需要分子氧作为辅因子（P450 = 单加氧酶）
- ESIBank 把 O=O 记录为底物是因为它参与了反应方程式
- 但 O=O 不是我们要预测特异性的有机底物，应该排除

**修复后排除列表：** O, OO, O=O, [O][O], [H]O[H], [OH2], [H][H], N=N, [H]O, O=[O]

**排除结果：**
- 总排除: 78条（原来只排除了4条，多排除了74条 O=O）
- 884 → **806 条有效正样本**
- 28个酶的唯一正样本就是 O=O，这些酶在 ESIBank 中丢失（但它们在 P450Rdb 中可能有其他底物记录）

### 3.4 反应 SMILES 提取

从 G 盘 ESIBank `reaction.csv` 中通过 reaction 索引获取完整反应 SMILES。
格式: `substrate.cofactor.O=O>>product.cofactor.water`（包含辅因子，未清洗）

806/806 条全部匹配成功。

### 3.5 提取结果

**输出文件（保存在 `data/sources/Source_ESIBank/`）：**

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 338个酶 | 全部有序列 |
| compounds.csv | 390个化合物 | 全部有SMILES |
| interactions.csv | 806条正样本对 | Tier B |
| reactions.csv | 806条反应记录 | 含完整反应SMILES（substrate>>product） |

**与 S3_P450Rdb 重叠统计：**

| 比较维度 | S2 数量 | S3 数量 | 重叠 | S2 重叠率 |
|----------|---------|---------|------|-----------|
| 酶 (UniProt) | 338 | 834 | 175 | 51.8% |
| 化合物 (SMILES) | 390 | 1,492 | 8 | 2.1% |
| 配对 (酶+底物) | 806 | 2,798 | 2 | 0.2% |

**解读：** 两个数据库收录了一半相同的 P450 酶，但它们从文献中找到的底物几乎完全不同。合并后能大幅增加酶-底物对的覆盖。

**提取脚本**: `scripts/01_数据下载/S2_ESIBank_extract.py`

---

## 4. S1_RCSB P450 子集提取

### 4.1 数据来源

PathA 已处理的数据（2026-01-08 ~ 01-21）：
- `独立测试集_682条.csv`：682条酶-配体记录，来自 RCSB PDB 晶体结构
  - 分类分布：substrate=84, inhibitor=164, product=36, unknown=398（原始PDB标注）
  - PathA 中经过人工审核重分类为：SUBSTRATE=271, INHIBITOR=245, PRODUCT=23, EXCLUDE=143
- `B1_仅底物_271pos/data.csv`：272条正样本（PathA审核后的底物）
- `Enzymes.csv`：292个酶（有UniProt、序列、PDB_ID、物种）
- `Substrates.csv`：436个底物（有SMILES）

### 4.2 提取方案

**使用 B1（审核后底物）作为正样本来源**，不用原始682条的 classification 字段。原因：
- 原始682条中的 "substrate" 只有84条，是PDB原始标注，很不完整
- PathA 中经过三方验证（Claude+Codex+Gemini逐条审核），最终确认了271条底物
- B1 是这个审核结果的最终输出

**质量等级：Tier A**（最高级），因为有实验解析的晶体结构直接证据。

### 4.3 提取结果

**输出文件（保存在 `data/sources/Source_RCSB/`）：**

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 103个酶 | 全部有序列，额外有 pdb_id 字段 |
| compounds.csv | 220个化合物 | 全部有SMILES |
| interactions.csv | 272条正样本对 | Tier A，额外有 pdb_id 字段 |

无辅因子排除（0条被过滤）。无 reactions.csv（RCSB 没有反应 SMILES 信息）。

### 4.4 与 S2/S3 重叠统计

| 比较 | 酶重叠 | 化合物重叠 | 配对重叠 |
|------|--------|----------|---------|
| S1 vs S2_ESIBank | **0/103 (0.0%)** | 16/220 (7.3%) | **0/272 (0.0%)** |
| S1 vs S3_P450Rdb | 54/103 (52.4%) | 9/220 (4.1%) | 5/272 (1.8%) |

**S1 和 S2 零酶重叠**：RCSB 的 P450 和 ESIBank 的 P450 完全是不同的酶。这是因为 ESIBank 主要收录 BRENDA 文献数据，而 RCSB 收录有晶体结构的酶——两个来源覆盖的酶种群不同。

**S1 和 S3 有 52.4% 酶重叠**：P450Rdb 覆盖面更广，包含了部分 RCSB 中有结构的酶。

**提取脚本**: `scripts/01_数据下载/S1_RCSB_extract.py`

---

## 5. 其他数据源调查结果

### S4 CYPED (cyped.biocatnet.de) — 无底物数据（2026-03-23重新调查）
- 斯图加特大学维护的P450工程数据库，之前认为服务器故障
- 实际可访问，有REST API（`/API/sequences`, `/API/proteins`）
- 发现41,514条P450蛋白记录、595个结构
- **但反应表为空**（`/API/reactions`返回空数组），无substrates/ligands关系
- Codex建议：对接阶段可用作P450结构索引和模板库（预组织的单体化PDB）
- 决策：**现阶段跳过**，对接阶段再用

### 旧S4 Plant P450 DB (p450.kvl.dk) — 不可用
- 主要是拟南芥 P450 的基因目录和序列数据
- **没有底物/化合物数据**
- 决策：**跳过**

### S9 PCPD (p450.biodesign.ac.cn) — 数据已下载（2026-03-23）

#### 9.1 数据库概述

PCPD（Plant Cytochrome P450 Database），中国科学院天津工业生物技术研究所维护。原始论文描述181个植物P450，但当前网站已扩展到**1,427个P450**，覆盖所有生物界（不仅植物）。

网站是React SPA，数据通过AWS API Gateway（需Cognito认证）和CloudFront CDN提供。

#### 9.2 数据获取过程

**找到数据入口**：
- 通过浏览器DevTools的Network标签，发现网站加载 `resource_new.json`
- 直接访问 `https://p450.biodesign.ac.cn/resource_new.json` 获取全部1,427条数据
- 通过分析JS源码找到CDN路径：`https://d1en57qlwrlmqu.cloudfront.net/P450/20240612/`

**Codex审核补充**：
- 确认 `resource_new.json` 是主要公开数据文件
- 找到PDB结构文件路径：`CDN/PDB/{CYP_ID}.pdb`
- 确认不存在结构化的反应SMILES数据（SDF/MOL/RXN/JSON全部404）
- 反应信息仅以图片形式存在

**可下载的全部资源**：

| 资源 | CDN路径 | 说明 |
|------|---------|------|
| 列表数据 | `p450.biodesign.ac.cn/resource_new.json` | 1,427条，含ID/Family/Kingdom/Function |
| FASTA序列 | `CDN/FASTA/{ID}.fasta` | 含物种名、UniProt ID、PDB ID |
| PDB结构 | `CDN/PDB/{ID}.pdb` | 预测结构（AlphaFold/Rosetta） |
| 反应图片 | `CDN/img_reaction/{ID}.png` | 底物→产物反应方程式 |

#### 9.3 下载结果

| 类型 | 成功 | 失败 | 说明 |
|------|------|------|------|
| FASTA序列 | 1,425 | 2 | 几乎全部有序列 |
| PDB结构 | 423 | 1,004 | 30%有结构（很多P450无预测结构） |
| 反应图片 | 857 | 570 | 60%有反应图片（Function=none的无图） |

#### 9.4 数据分布

| Kingdom | 数量 | 占比 |
|---------|------|------|
| Viridiplantae（植物） | 578 | 40.5% |
| Fungi（真菌） | 480 | 33.6% |
| Bacteria（细菌） | 182 | 12.7% |
| Metazoa（动物） | 175 | 12.3% |
| OtherEukaryotes | 10 | 0.7% |
| Archaea | 1 | 0.1% |
| Virus | 1 | 0.1% |

有Function描述（非none）：1,317个（92.3%）

#### 9.5 FASTA header结构

每个FASTA文件的header包含丰富信息（格式不完全统一）：
```
>Bacteria_CYP101A1 Bacteria P450cam;Camphor 5-monooxygenase pdb_id=2zwu_A [Pseudomonas putida]
>Fungi_CYP652D1 sdnE;A0A1B4XBH0;sdnE like [Sordaria araneosa]
>Viridiplantae_CYP71BN1 311_3_CYP71BN1;K4CI56;lycosantalene oxygenase [Solanum lycopersicum]
```
可解析出：Kingdom、CYP ID、UniProt ID（部分）、PDB ID（部分）、功能描述、物种名（[方括号]内）

#### 9.6 反应图片的价值（Path D 区域选择性预测）

857张反应图片包含**底物→产物**的化学结构转变图，是区域选择性预测的珍贵数据源。

**区域选择性预测**需要知道底物分子上**哪个位置**被酶修饰（如C-5位羟基化）。通过对比底物SMILES和产物SMILES，可以用RDKit做原子映射（atom mapping）找到反应位点。

**图片→结构化数据的转换方案（Path D启动时执行）**：
1. 用 **RxnScribe**（MIT, 2023）识别反应图片 → 拆分底物图和产物图
2. 用 **MolScribe/DECIMER** 把分子图转SMILES
3. 用大模型做质量审核

**与其他数据源的反应数据对比**：

| 数据源 | 有产物SMILES？ | 物种覆盖 |
|--------|---------------|---------|
| S3_P450Rdb | ✅ 结构化SMILES | 跨物种（植物57%+微生物+动物） |
| S9_PCPD | ✅ 图片形式 | 全覆盖（真菌34%是S3没有的） |

PCPD的反应图片补充了S3没有的真菌P450反应数据。

#### 9.7 待办

1. **解析FASTA header** → 提取物种名、UniProt ID、PDB ID
2. **从Function文本提取底物名** → PubChem查SMILES
3. **生成 Source_PCPD/ 标准文件**
4. **可选（Path D）**：反应图片→结构化反应SMILES

**关键文件（保存在 `downloads/PCPD/`）：**
- `resource_new.json` — 1,427条P450列表数据
- `FASTA/` — 1,425个FASTA序列文件
- `PDB/` — 423个PDB结构文件
- `img_reaction/` — 857张反应方程式图片

**提取脚本**: `scripts/01_数据下载/S9_PCPD_download.py`

### S7 P450 BM3 Variants DB — 不适用
- >1,500个BM3突变体的反应数据
- 只涉及单一酶 CYP102A1 的变体
- 适合变体级别分析，不适合酶家族多样性研究
- 决策：**跳过**

---

## 6. S6_Figshare CYP450 — 辅助数据源

### 6.1 数据来源

Figshare "Curated CYP450 Interaction Dataset"（article 26630515, Scientific Data 2025）

专门针对负责人类 90% I相药物代谢的6个CYP同工酶，收集了约2000个化合物/酶的底物/非底物分类数据。

### 6.2 为什么不放入主benchmark

- 只有6个人类CYP（CYP1A2, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4）
- 加入3,610条正样本会让6个酶占总数据60%，严重扭曲enzyme_split和all_split
- 这6个CYP已在S2和S3中存在
- 会增加不成比例的对接工作量（6个酶的上千条对接）

### 6.3 辅助用途

1. **化合物池扩展**：3,258个唯一化合物（含药物分子），其中1,954个是S2/S3中没有的新化合物 → 极大丰富底物换向负样本的化合物池
2. **确认生物学负样本**：11,395条非底物标签（label=0）是文献验证的"这个化合物不是这个酶的底物" → 比随机生成的负样本质量高得多
3. **注意**：非底物标签是酶-化合物对级别的，不是化合物全局的。一个化合物对CYP1A2是非底物，但可能对CYP3A4是底物

### 6.4 提取结果

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 6个酶 | 人类CYP，无序列（需从UniProt获取） |
| compounds.csv | 3,258个化合物 | 全部有SMILES和化合物名称 |
| interactions_positive.csv | 3,610条正样本 | label=1, usage=auxiliary |
| biological_negatives.csv | 11,395条负样本 | label=0, usage=biological_negative |

正样本和负样本**分开存储**，避免被意外合并到主数据集。

**提取脚本**: `scripts/01_数据下载/S6_Figshare_CYP450_extract.py`

---

## 7. 所有数据源最终汇总

### 主 benchmark 数据源（用于训练和4场景评估）

| 数据源 | 正样本对 | 唯一酶 | 唯一底物 | Tier | 状态 |
|--------|---------|--------|---------|------|------|
| S1_RCSB | **272** | **103** | **220** | A (晶体结构) | ✅ 已提取 |
| S2_ESIBank | **806** | **338** | **390** | B (文献验证) | ✅ 已提取 |
| S3_P450Rdb | **2,798** | **857** | **1,492** | B (实验验证) | ✅ 已提取 |

**去重前总正样本**: 272 + 806 + 2,798 = 3,876条
**配对重叠**: S1∩S3=5, S2∩S3=2, S1∩S2=0 → 去重后约 **3,869条**
**去重后预估唯一酶**: ~1,069个

### 辅助数据源（化合物池 + 生物学负样本）

| 数据源 | 正样本 | 负样本 | 唯一化合物 | 用途 |
|--------|--------|--------|-----------|------|
| S6_Figshare | 3,610 | 11,395 | 3,258 | 化合物池扩展 + 确认负样本 |

### 不可用/跳过

| 数据源 | 调查结果 | 决策 |
|--------|---------|------|
| S4_CYPED | 服务器SSL证书错误 + API返回PHP Fatal Error，网站处于维护不善状态，无法访问数据 | **跳过（服务器故障）** |
| S5_EnzymeMap | 已下载（349K反应）。EC 1.14有8,704条、515个UniProt。但与S2+S3重叠检查：477个"新"ID中只有2个在P450目录中，说明绝大部分是非P450单加氧酶（酪氨酸羟化酶等） | **跳过（P450含量极低）** |
| S7_BM3_Variants | >1,500个CYP102A1突变体，单一酶的变体数据 | **跳过（不适合家族研究）** |
| S9_PCPD | 网站是JS SPA难以爬取，GitHub只有配置文件，181个酶规模小 | **跳过（提取困难）** |

---

## 8. S8_PlantP450DB — 植物 P450 数据库

### 8.1 数据来源

Plant Cytochrome P450 Database，哥本哈根大学维护，托管在 ERDA（`erda.dk/public/vgrid/PlantP450/`）。收录"已知有内源底物的植物 P450"，是植物 P450 领域最权威的分类数据库。

引用论文：Hansen et al. (2021) Molecular Plant 14:1244-1265

数据库特点：
- 913 个条目，每个条目 = 一个 CYP 基因（CYP名+物种+Clan/Family）
- 详情页包含 Function、Compound class、Pathway、Accession、References、DOIs
- 底物信息分散在文献中，不是结构化字段，需要从论文中提取

### 8.2 提取过程（多轮 agent 协作）

整个提取过程分 8 步，历时多个会话，涉及 20+ 个子 agent 并行工作。

#### Step 1：爬取列表页
- 主页 table.html 有 913 个条目（CYP名+物种+Clan/Family）
- 每个条目有详情页（`sub/{CYP名}.html`），包含 Function、Compound class、Pathway、Accession、References、DOIs

#### Step 2：爬取全部 910 个详情页
- 3 个条目页面不存在，实际爬取 910 个
- 数据完整度：72.2% 有 Compound class，27.7% 有 Function，99.7% 有 DOI
- Function 字段含有底物信息的只有 252 条，其余需从论文确认

#### Step 3：获取 567 篇论文摘要
- 910 条条目引用了 567 个唯一 DOI
- PubMed 获取摘要成功 546 篇，CrossRef 补充 20 篇，1 篇失败
- 保存在 `downloads/PlantP450DB/abstracts.jsonl`

#### Step 4：20 个子 agent 并行分析摘要
- 910 条分 20 组，每组 ~46 条，并行调用 agent 分析
- 每个 agent 的任务：读取条目的 Function/Compound class 字段 + 对应论文摘要，判断能否确定底物
- 规则：只提取摘要中明确写了的底物名称，不确定的标记 needs_fulltext
- 结果统计：
  - 589 条 found (64.7%) — 从摘要中找到了具体底物
  - 310 条 needs_fulltext — 摘要信息不足，需要看全文
  - 11 条 no_info — 完全没有底物相关信息

#### Step 5：获取开放获取全文
- 310 条 needs_fulltext 涉及 125 篇唯一论文
- Unpaywall API 查询：76 篇有开放获取版本
- 下载方式：PMC XML + 期刊 HTML 两个通道
- 最终成功下载 52 篇全文，覆盖 130 条 needs_fulltext 记录
- 保存在 `downloads/PlantP450DB/track2_fulltext/fulltext.jsonl`

#### Step 6：9 个子 agent 分析全文
- 52 篇全文分 9 组，并行分析
- 每个 agent 读取全文内容，寻找 CYP 对应的底物名称和催化反应描述
- 结果：额外找到 88 条底物信息（从 needs_fulltext → found）
- 全文分析比摘要分析多找到约 68% 的底物信息

#### Step 7：用户手动下载 PDF + PDF 分析
- 剩余 70 篇论文需要用户手动下载（付费期刊，无 OA 版本）
- 用户实际下载了 61 篇 PDF（11 篇 Science/Nature 论文无法获取）
- 10 个子 agent 并行分析 61 篇 PDF
- 结果：额外找到 132 条底物信息
- 部分论文虽有下载但内容与目标 CYP 不直接相关

#### Step 8：质量验证（QC）
- 对 44 条疑似 agent 推断（非文献直接证据）的条目进行复核
- 2 个子 agent 重新审查，对照原始文献判断可靠性
- 结果：
  - 4 条 verified（确认正确）
  - 33 条 plausible（合理推断但非直接证据，保留使用）
  - 7 条 rejected（证据不足，改回 needs_fulltext）
- 最终 plausible 条目保留为正样本（来源标记为 plausible）

#### Step 9：底物名称转 SMILES（多轮查询）

**第一阶段（Session 1, 2026-03-22）：4 轮查询**
- 从 found 条目中提取 432 个唯一底物名称
- 查询策略（4轮递进）：
  1. 原始名称 → PubChem 精确匹配
  2. 清洗名称（去除括号内修饰语、统一大小写）→ PubChem
  3. 括号内简称单独查（如 "flavanone (naringenin)" → 查 naringenin）
  4. 缩写词展开（如 GA12 → gibberellin A12）→ PubChem + KEGG
- 结果：336 个成功 (77.8%)，96 个失败
- 但只覆盖了 432/480 个底物名（PDF/全文提取的 127 个新名字未查询）

**第二阶段（Session 2, 2026-03-23）：完整覆盖 + 3 轮重试**

发现 smiles_cache.jsonl 只有 432 条，但 final_merged_results.jsonl 中有 480 个唯一底物名。缺 127 个（全文/PDF提取的新名字从未查过PubChem），另有 79 个 orphan（cache中有但results不再引用，名字格式变了）。

完整重建 smiles_cache.jsonl：
1. **精确匹配**：353 个名字在旧 cache 中精确匹配
2. **模糊匹配**：11 个名字通过去除括号注释匹配到 orphan 的 SMILES（如 `(+)-abscisic acid` ← `(+)-abscisic acid (ABA)`）
3. **PubChem Round 1**：116 个名字查 PubChem，找到 55 个
4. **Codex 审核**：发现 3 个错误 SMILES（8-oxogeraniol 返回染料、compound 3 返回含卤素药物、13-HPOT 位置错误）+ 1 个需修正（zealexin A3 阴离子形式）
5. **PubChem Round 2**：用替代名字重试（如 cheilanthifoline → (S)-cheilanthifoline），找到 29 个
6. **PubChem Round 3**：最后一轮（如 furanocoumarins → psoralen 代表性化合物），找到 12 个
7. **[O-] 质子化**：9 个 SMILES 中的 `[O-]` 改为 `O`（统一为中性形式）

最终结果：**480 个底物名全部覆盖**
- 有 SMILES：**422 个 (87.9%)**
- 无 SMILES：**58 个**
  - 46 个类别描述（如 "C12-C18 fatty acids"、"diterpene alcohol intermediates"），无法对应单一化合物
  - 12 个冷门化合物（如 "beta-macrocarpene"、"cysteine-indole-3-acetonitrile"），PubChem/KEGG 均无收录

### 8.3 数据来源分布

按底物信息获取来源统计（806 条 found）：

| 来源 | 条数 | 说明 |
|------|------|------|
| 摘要分析（Step 4） | 556 | 20 个 agent 从 567 篇摘要中提取 |
| PDF 全文分析（Step 7） | 132 | 10 个 agent 从 61 篇用户下载 PDF 中提取 |
| HTML 全文分析（Step 6） | 88 | 9 个 agent 从 52 篇 OA 全文中提取 |
| plausible 推断（Step 8） | 27 | QC 复核后保留的合理推断 |
| QC verified（Step 8） | 3 | QC 确认的条目 |

### 8.4 最终结果（2026-03-23 更新）

**910条原始记录的完整去向：**

| 类别 | 数量 | 说明 |
|------|------|------|
| found | 818 条 | 从摘要/全文/PDF中找到了底物名称 |
| needs_fulltext | 58 条 | 摘要无法确定，论文下不到 |
| partial | 23 条 | 部分信息（19条有底物，4条无） |
| no_info | 11 条 | 数据库条目本身无底物相关内容 |

**837条有底物的记录展开为 1,357 条酶-底物对（一个酶可能催化多个底物）**

**SMILES 查询结果（480个唯一底物名，全部已查询）：**
- 有 SMILES：422 个底物名 → **1,277 条酶-底物对可用**
- 无 SMILES：58 个底物名 → 80 条酶-底物对不可用
  - 46 个类别描述（如 "C12-C18 fatty acids"）
  - 12 个冷门化合物（PubChem/KEGG 无收录）

| 指标 | 数值 |
|------|------|
| **最终可用酶-底物对** | **1,277 条** |
| 唯一 CYP 名称 | ~656 个 |
| 唯一底物（有 SMILES） | 422 个 |
| 缺 UniProt 序列 | 全部（需查询） |

**关键决策和问题：**

1. **plausible 条目的处理**：33 条 plausible 标记的条目保留为正样本。理由是这些底物来自合理的文献推断（同一 CYP 家族的已知功能），且 P450 数据本身稀缺
2. **化合物类别 vs 具体分子**：57 个 "底物" 实际是类别名（如 "terpenoids"），无法转为 SMILES。需要后续从全文中找到具体分子名
3. **Science/Nature 论文**：11 篇高影响力论文用户无法下载，影响约 30-40 条 CYP 的底物信息。这些论文覆盖的 CYP 可能在其他来源中有重叠
4. **植物 P450 的独特价值**：S8 带来 557 个植物 P450 酶，与 S1-S3 几乎不重叠（S1-S3 主要是微生物/动物/人类 P450）

**关键文件（保存在 `downloads/PlantP450DB/`）：**
- `all_entries.json` — 910 条原始爬取数据
- `final_merged_results.jsonl` — 最终分析结果（含 QC 标记和底物名称）
- `smiles_cache.jsonl` — SMILES 查询缓存（480 条全覆盖：422 有 SMILES + 58 无 SMILES）
- `pdfs/` — 用户手动下载的 72 篇 PDF

### 8.5 UniProt 蛋白序列查询（2026-03-23）

#### 为什么需要查 UniProt

Plant P450 DB 只记录了 CYP 名称和物种，没有蛋白序列。但我们的模型需要蛋白序列做 ESM-2 编码，所以必须获取每个酶的氨基酸序列。

#### 植物 P450 的发现流程（为什么很多酶在 UniProt 里查不到）

一个植物P450酶被发现和命名的典型流程如下：

1. **克隆基因+测序**：研究者从植物的基因组或 mRNA 中找到一段 DNA 序列（核苷酸序列，ATGCTTGAC...），提交到 **GenBank**（NCBI 维护的核酸序列数据库），得到 accession 号（如 AK250744）
2. **DNA 翻译成蛋白序列**：DNA 序列可以按遗传密码表直接翻译成蛋白质的氨基酸序列（每3个碱基对应1个氨基酸，如 ATG→M, CTT→L, GAC→D），所以 DNA 序列 = 蛋白序列，只是存储形式不同
3. **CYP 命名**：Nelson 命名委员会拿到蛋白序列后，和所有已知 P450 做序列比对（逐个氨基酸位置对齐，计算相同位置的比例）。>40% 相同 → 同一家族（如 CYP71），>55% 相同 → 同一亚族（如 CYP71C），然后按发现顺序编号 → CYP71C113
4. **功能验证**：在酵母/大肠杆菌中表达该蛋白，测试它能催化哪些底物 → 确定酶-底物关系
5. **数据入库**：GenBank 有 DNA 序列 ✅，论文报道了功能 ✅，Plant P450 DB 收录了 CYP 名+底物 ✅，但 **UniProt 可能还没收录**（需要人工审核或批量导入，很多新发现的植物 P450 尤其是非模式植物的还没被导入）

**这就是为什么 205 个酶在 UniProt 里查不到**——序列存在于 GenBank（DNA 形式），但 UniProt（蛋白形式数据库）还没收录。

#### 查询过程

**Round 1：按 CYP 名+物种搜索 UniProt（783 个酶）**
- 3种策略：gene name 搜索、protein name 搜索、全文搜索
- 结果：413/783 找到（52.7%）

**Round 2：用 GenBank accession 交叉引用（370 个 misses）**
- 从 all_entries.json 中获取 GenBank accession 号
- 224 个 miss 有 accession，通过 UniProt 的 EMBL 交叉引用搜索
- 另外尝试只用物种属名（genus）搜索
- 结果：新找到 165 个（159 通过 xref + 6 通过 genus 搜索）

**最终：578/783 有蛋白序列（73.8%），205 个没有**

205 个查不到的原因分布：
- 137 个：Plant P450 DB 中没有 GenBank accession 记录，无法进一步追溯
- 65 个：有 GenBank accession 但 UniProt 未收录（序列只存在于 NCBI GenBank 的 DNA 形式）
- 3 个：非标准命名（如 "Ai BX4"），无法匹配

**可选的后续恢复方案**（暂不执行）：
- 65 个有 GenBank accession 的可以从 NCBI Nucleotide API 下载 DNA 序列 → 翻译成蛋白序列
- 137 个无 accession 的可以从论文 DOI 中提取 accession 号（需要 agent 读论文）

#### 查询脚本
- `scripts/01_数据下载/S8_uniprot_query.py` — Round 1
- `scripts/01_数据下载/S8_uniprot_retry.py` — Round 2
- 结果文件：`downloads/PlantP450DB/uniprot_results.jsonl`

### 8.6 S8 标准化文件生成（2026-03-23）

只保留**同时满足两个条件**的酶-底物对：酶有蛋白序列 + 底物有 SMILES。

**完整损耗链：**
```
910 条原始记录
  ↓ -73 条没有底物信息
837 条有底物 → 展开为 1,357 条酶-底物对
  ↓ -80 条底物没有 SMILES（46 类别描述 + 12 冷门 + 其他）
1,277 条有 SMILES 的酶-底物对
  ↓ -296 条酶没有蛋白序列（205 个酶在 UniProt 查不到）
  ↓ -2 条内部去重（同酶同底物）
979 条最终可用的酶-底物对
```

**输出文件（保存在 `data/sources/Source_PlantP450DB/`）：**

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 578 个酶 | 全部有蛋白序列，无 UniProt 的用 SEQHASH_ 前缀 |
| compounds.csv | 295 个化合物 | 全部有 SMILES，按 SMILES 去重 |
| interactions.csv | 979 条正样本对 | label=1, source=S8_PlantP450DB, quality_tier=B |
| unresolved.csv | 279 条丢弃记录 | 205 个无序列 + 73 个无底物 |

**列格式与 S1-S3 完全一致**，可直接接入合并管线。

**提取脚本**: `scripts/01_数据下载/S8_generate_source_files.py`

### 8.7 待办（时间允许时）

1. 从 NCBI GenBank 恢复 65 个有 accession 的酶序列（翻译 DNA→蛋白）
2. 从论文中提取 137 个无 accession 酶的序列信息
3. 从 KEGG/MetaCyc 恢复部分缺失 SMILES

---

## 9. S9 PCPD 提取（2026-03-23）

### 9.1 数据源发现

PCPD（Plant Cytochrome P450 Database，http://p450.biodesign.ac.cn/）是天津工业生物技术研究所维护的数据库。虽然名字叫"Plant"，实际包含植物、真菌、细菌、动物等多个物种类群的P450。

网站是React SPA，通过浏览器DevTools的Network面板发现数据加载自 `resource_new.json`。进一步在JS代码中找到CloudFront CDN地址，试出了FASTA/PDB/反应图片的下载路径。

### 9.2 下载内容

| 资源 | CDN路径 | 成功数 | 说明 |
|------|---------|--------|------|
| JSON列表 | `p450.biodesign.ac.cn/resource_new.json` | 1,427条 | ID, Family, Kingdom, Function |
| FASTA序列 | `cloudfront.net/.../FASTA/{ID}.fasta` | 1,425 | 含物种名+序列+UniProt(部分)+PDB(部分) |
| PDB结构 | `cloudfront.net/.../PDB/{ID}.pdb` | 423 | 很多P450没有预测结构 |
| 反应图片 | `cloudfront.net/.../img_reaction/{ID}.png` | 857 | 底物→产物反应方程式图 |

2个FASTA缺失：CYP113E1、CYP3A38（CDN上不存在）。

反应图片（底物→产物化学结构图）**没有结构化数据**，只有图片格式。Codex确认CDN上无SDF/MOL/RXN文件，JS中SMILES/SDF仅用于用户提交查询。图片保存供Path D区域选择性预测使用（可用RxnScribe+MolScribe转换为SMILES）。

### 9.3 FASTA解析

FASTA header格式：`>Kingdom_ID description pdb_id=xxx [Species]`

从1,425个header中提取：
- 物种名（1,425个全有）
- UniProt ID（664个）
- PDB ID（152个）
- 功能描述文本（比JSON的Function有时更详细）
- 蛋白序列

**因为FASTA已包含序列，不需要像S8那样单独查UniProt。**

### 9.4 底物提取（三阶段）

**第一阶段：规则提取**（7种正则模式）
- "X monooxygenase/hydroxylase" → 底物X
- "hydroxylation of X" → 底物X
- "from X to Y" → 底物X
- 管道符分隔的多底物
- 等等
- 结果：654 specific, 57 class_only, 582 unclear

**第二阶段：10个opus agent并行处理582个unclear**
- 每批~59条，10个agent同时跑
- 结果：450 specific, 61 class_only, 71 no_substrate

**第三阶段：Codex审核 + 清洗**
Codex抽查发现系统性错误：
- FASTA占位符污染（"An Organic Molecule"被当成底物）
- 反应描述当底物（"regio- and stereoselective hydroxylation"）
- 标识符混入（PDB ID、UniProt ID）
- 底物名带多余文本（"the natural product compactin"→"compactin"）

清洗脚本修复：移除110个错误底物，修正57个名字，修正61个状态标签。

### 9.5 SMILES查询（三轮）

**Round 1：PubChem查询**（866个唯一底物名）
- 首次运行因字段名bug全部失败（PubChem返回"SMILES"不是"IsomericSMILES"）
- 修复后：615/866找到（71.0%）

**Round 2：清洗+重试**（Codex审核后）
- 拆分多底物名（"nicotine and cotinine"→分别查）
- 修复拼写（"aplpha"→"alpha"）
- 去掉酶名后缀（"Mevastatin hydroxylase"→"Mevastatin"）
- 展开缩写（VD3→Vitamin D3, OMST→O-methylsterigmatocystin）
- 结果：恢复26个

**Round 3：KEGG + PubChem变体查询**（252个仍缺SMILES的）
- PubChem变体查询（去立体化学前缀、大小写变化等）：9个
- KEGG化合物搜索：19个
- 结果：恢复28个

最终：**643个底物有SMILES**（866个中的74.2%），223个确实查不到（冷门天然产物中间体）。

### 9.6 筛选漏斗

```
1,427 JSON条目
  ↓ -2 无FASTA序列
1,425 有序列
  ↓ -108 Function="none"
1,317 有功能描述
  ↓ -276 (28基因名 + 131无法确定底物 + 119仅有类别名)
1,041 提取到具体底物名
  ↓ -223 底物无SMILES (PubChem+KEGG都查不到)
  818 有底物+有SMILES+有序列
  ↓ 展开为酶-底物对，去重
1,209 条最终可用酶-底物对（818酶 × 570化合物）
```

### 9.7 最终标准文件

输出到 `data/sources/Source_PCPD/`：

| 文件 | 数量 | 说明 |
|------|------|------|
| enzymes.csv | 818个酶 | 全部有蛋白序列（来自FASTA） |
| compounds.csv | 570个化合物 | 全部有SMILES，按SMILES去重 |
| interactions.csv | 1,209条正样本对 | label=1, source=S9_PCPD, quality_tier=B |
| unresolved.csv | 609条丢弃记录 | 原因分布见漏斗 |

818个酶的物种分布：
- 植物 (Viridiplantae): 380
- 真菌 (Fungi): 216
- 动物 (Metazoa): 117
- 细菌 (Bacteria): 84
- 其他: 2

### 9.8 附加资源（Path D区域选择性预测）

PCPD的857张反应图片（底物→产物化学结构图）是区域选择性预测的珍贵数据源。转换方案：
1. **RxnScribe**：识别反应图片→拆分底物图和产物图
2. **MolScribe/DECIMER**：单个分子图→SMILES
3. **Claude/GPT-4o审核**：验证SMILES正确性

这些专用化学图像识别工具比通用多模态LLM更精确（对计算机生成的清晰结构图准确率95%+）。

### 9.9 关键脚本

- `scripts/01_数据下载/S9_PCPD_download.py` — 批量下载FASTA+PDB+图片
- `scripts/01_数据下载/S9_extract_substrates.py` — 规则提取底物名
- `scripts/01_数据下载/S9_cleanup_substrates.py` — 清洗（Codex审核后修复）
- `scripts/01_数据下载/S9_smiles_lookup.py` — PubChem SMILES查询
- `scripts/01_数据下载/S9_generate_source_files.py` — 生成标准文件

---

## 10. 所有数据源最终汇总

### 主 benchmark 数据源（用于训练和4场景评估）

| 数据源 | 正样本对 | 唯一酶 | 唯一底物 | Tier | 状态 |
|--------|---------|--------|---------|------|------|
| S1_RCSB | **272** | **103** | **220** | A (晶体结构) | ✅ 已提取 |
| S2_ESIBank | **806** | **338** | **390** | B (文献验证) | ✅ 已提取 |
| S3_P450Rdb | **2,798** | **857** | **1,492** | B (实验验证) | ✅ 已提取 |
| S8_PlantP450DB | **979** | **578** | **295** | B (文献验证) | ✅ 标准文件已生成 |
| S9_PCPD | **1,209** | **818** | **570** | B (跨物种) | ✅ 标准文件已生成 |

**去重前总正样本**: 272 + 806 + 2,798 + 979 + 1,209 = **6,064条**
**预估唯一酶**: ~2,200+（S9带来818个跨物种P450，含真菌/动物/细菌）

### 辅助数据源（化合物池 + 生物学负样本）

| 数据源 | 正样本 | 负样本 | 唯一化合物 | 用途 |
|--------|--------|--------|-----------|------|
| S6_Figshare | 3,610 | 11,395 | 3,258 | 化合物池扩展 + 确认负样本 |

### 附加资源（Path D区域选择性预测）

| 资源 | 来源 | 数量 | 说明 |
|------|------|------|------|
| 反应SMILES | S3_P450Rdb | 3,352条 | 结构化底物→产物 |
| 反应图片 | S9_PCPD | 857张 | 需RxnScribe+MolScribe转换 |

---

## 11. 下一步

1. **继续搜集其他数据库**（用户尚未完成所有数据源调查）
2. **合并去重**：将 S1+S2+S3+S8+S9（+其他）合并为主数据集
3. **双向负样本生成**
4. **4种Split生成**
5. **对接与特征生成**（需要服务器）
6. **可选**：从 NCBI GenBank 恢复 S8 中 65 个缺失酶的序列
