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

### S4 Plant P450 DB (p450.kvl.dk) — 不可用
- 主要是拟南芥 P450 的基因目录和序列数据
- 有 FASTA 序列可下载，有染色体定位和系统发育树
- **没有底物/化合物数据**，无法提供酶-底物配对
- 决策：**跳过**（只有序列价值，无配对价值）

### S5 PCPD (p450.biodesign.ac.cn) — 暂不可用
- 181个植物P450，有序列、结构、功能信息
- 网站是JavaScript SPA，难以自动提取
- GitHub (JiangLab2020/PCPD) 只有配置文件
- 181个酶规模小，提取成本高，收益不确定
- 决策：**暂时跳过**，时间允许再考虑

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

Plant Cytochrome P450 Database，哥本哈根大学维护，托管在 ERDA（`erda.dk/public/vgrid/PlantP450/`）。收录"已知有内源底物的植物 P450"。

引用论文：Hansen et al. (2021) Molecular Plant 14:1244-1265

### 8.2 提取过程

**Step 1：爬取列表页**
- 主页 table.html 有 913 个条目（CYP名+物种+Clan/Family）
- 每个条目有详情页（`sub/{CYP名}.html`），包含 Function、Compound class、Pathway、Accession、References、DOIs

**Step 2：爬取全部 910 个详情页**
- 数据完整度：72.2% 有 Compound class，27.7% 有 Function，99.7% 有 DOI
- Function 字段含有底物信息的只有 252 条，其余需从论文确认

**Step 3：获取 567 篇论文摘要**
- 567 个唯一 DOI → PubMed 获取摘要
- 成功 546 篇（PubMed），补充 20 篇（CrossRef），1 篇失败
- 保存在 `downloads/PlantP450DB/abstracts.jsonl`

**Step 4：20 个子 agent 并行分析摘要**
- 910 条分 20 组，每组 ~46 条，并行分析
- 规则：只提取摘要中明确写了的底物，不确定的标记 needs_fulltext
- 结果：589 条 found (64.7%)，310 条 needs_fulltext，11 条 no_info

**Step 5：获取 OA 全文**
- Unpaywall 查询：125 篇论文中 76 篇有开放获取版本
- PMC + HTML 方式下载成功 52 篇全文

**Step 6：9 个子 agent 分析全文**
- 52 篇全文覆盖 130 条 needs_fulltext
- 额外找到 92 条底物信息

**Step 7：质量验证**
- 对 44 条疑似 agent 推断（非摘要直接证据）的条目，用 2 个子 agent 重新验证
- 结果：4 条 verified，33 条 plausible（合理但非直接证据），7 条 rejected（改回 needs_fulltext）

**Step 8：底物名称转 SMILES**
- 432 个唯一底物名称 → PubChem + KEGG 查询
- 4 轮查询（原始名 → 清洗名 → 括号内简称 → 缩写词）
- 成功 336 个 (77.8%)，失败 96 个（57 个化合物类别 + 39 个专业中间产物）

### 8.3 当前结果

| 类别 | 数量 | 说明 |
|------|------|------|
| **可直接使用** | **898 条酶-底物对** | 557 个植物 P450 + 336 个底物（有 SMILES） |
| 有底物名但缺 SMILES | 115 条 | 57 个是化合物类别需全文确认，47 个是专业中间产物 |
| 缺底物信息 | 202 条 | 需下载 70 篇论文后用子 agent 分析 |
| 部分信息 | 23 条 | 有线索但不完整 |
| 无信息 | 11 条 | |

**关键文件（保存在 `downloads/PlantP450DB/`）：**
- `all_entries.json` — 910 条原始爬取数据
- `abstracts.jsonl` — 567 篇论文摘要
- `final_merged_results.jsonl` — 最终分析结果（含 QC 标记）
- `track1_output/smiles_cache.jsonl` — SMILES 查询缓存
- `track2_fulltext/需要手动下载的论文_全部70篇.csv` — 待手动下载的论文列表
- `track2_fulltext/fulltext.jsonl` — 已下载的 52 篇全文

### 8.4 待办

1. **用户手动下载 70 篇论文**（优先前 10 篇，覆盖 80+ CYP）→ 放到 `track2_fulltext/pdfs/`
2. 下载后用子 agent 分析全文提取底物
3. 跑 UniProt 查询获取蛋白序列
4. 生成 S8 最终提取文件（`data/sources/Source_PlantP450DB/`）

---

## 9. 所有数据源最终汇总

### 主 benchmark 数据源（用于训练和4场景评估）

| 数据源 | 正样本对 | 唯一酶 | 唯一底物 | Tier | 状态 |
|--------|---------|--------|---------|------|------|
| S1_RCSB | **272** | **103** | **220** | A (晶体结构) | ✅ 已提取 |
| S2_ESIBank | **806** | **338** | **390** | B (文献验证) | ✅ 已提取 |
| S3_P450Rdb | **2,798** | **857** | **1,492** | B (实验验证) | ✅ 已提取 |
| S8_PlantP450DB | **898** | **557** | **336** | B (文献验证) | ✅ 部分可用（559条有SMILES） |

**去重前总正样本**: 272 + 806 + 2,798 + 898 = 4,774条
**预估唯一酶**: ~1,600+（S8 带来大量植物 P450，与 S1-S3 几乎不重叠）

### 辅助数据源（化合物池 + 生物学负样本）

| 数据源 | 正样本 | 负样本 | 唯一化合物 | 用途 |
|--------|--------|--------|-----------|------|
| S6_Figshare | 3,610 | 11,395 | 3,258 | 化合物池扩展 + 确认负样本 |

---

## 10. 下一步

1. **用户下载 70 篇论文** → 用子 agent 分析，可额外解锁 ~200 条数据
2. **UniProt 查询**：获取所有酶的蛋白序列
3. **S8 最终提取**：生成 `data/sources/Source_PlantP450DB/` 的标准文件
4. **合并去重**：将 S1+S2+S3+S8 合并为主数据集
5. **双向负样本生成**
6. **4种Split生成**
7. **对接与特征生成**（需要服务器）
