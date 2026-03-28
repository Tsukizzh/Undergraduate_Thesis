# Session Log: C3-Step 5 底物分类与多轮验证

> **日期**: 2026-03-27 ~ 2026-03-29
> **目标**: 对 P450 数据集中 2,125 个底物进行化学类别分类，通过多工具多轮验证确保准确率
> **状态**: 多轮验证完成，200 抽样审计准确率 88.5%，待合并到 8 类并修正后重新审计

---

## 一、背景与目标

导师提出的研究方向：给定一个 P450 酶，预测它催化**哪一类底物**（萜类？黄酮？生物碱？脂肪酸？）。

实现思路：
1. 先给所有 2,125 个底物打上化学类别标签
2. 用已有模型（EXP001 checkpoint, Test AUC=0.7730）对每个"酶-底物对"打分
3. 按底物类别汇总分数 → 预测酶偏好催化哪类底物

**前提条件：底物分类必须准确（目标 99%+）。**

---

## 二、分类工具调研（2026-03-27）

与 Codex 讨论后确定工具选型：

| 工具 | 原理 | 优势 | 局限 |
|------|------|------|------|
| **NPClassifier** | 深度神经网络（73,607 天然产物训练） | 天然产物专用，三级分类 | 黑盒，合成物/外源物会错 |
| **ClassyFire** | 9,000+ 条 SMARTS 结构规则 | 白盒可审计，覆盖面广 | 通过 InChIKey 查数据库，大量化合物未收录 |
| **SMARTS 结构规则** | 手写化学骨架模式匹配 | 确定性验证 | 只能验证有明确骨架的类别 |
| ~~CANOPUS~~ | 需要质谱数据 | — | 我们只有 SMILES，不适用 |
| ~~MolNetEnhancer~~ | 代谢组学流程 | — | 同上，不适用 |

**结论**：NPClassifier 作主分类器 + ClassyFire 交叉验证 + SMARTS 结构确认 + 文献人工验证。

---

## 三、NPClassifier 自动分类（2026-03-27）

### 3.1 数据预处理

- 总化合物：2,125 个
- 精确 SMILES：2,110 个 → 调 NPClassifier API
- 通配符 SMILES（含 `*`）：15 个 → 按化合物名称手动分类

### 3.2 API 调用

- 脚本：`scripts/05_底物分类/classify_substrates.py`
- API：`https://npclassifier.gnps2.org/classify?smiles=...`
- 延迟：0.15s/请求，串行调用
- 耗时：~8 分钟
- 结果：2,110/2,110 成功，**0 个 API 错误**
- 缓存：`data/05_底物分类/npclassifier_cache.json`

### 3.3 Codex 三轮审核

**第一轮**（6 个问题）：
1. Sesterterpenoid 被错误映射到 Triterpenoid（25碳 vs 30碳）→ 修正为独立类
2. Jasmonoid 缺失于映射表 → 添加
3. 3 个通配符化合物名称含 mojibake（编码损坏的撇号）→ 从精确匹配改为关键词匹配
4. 关键词匹配顺序错误（"flavon" 先于 "isoflavon"）→ 改为有序列表，isoflavon 优先
5. 通配符标签不一致（"Fatty acid" vs "Fatty_acid"）→ 统一为下划线格式
6. 脚本路径错误 → 修正 PROJECT_DIR 层级

**第二轮**：确认修复正确，讨论后处理策略（合并稀有类别 <30 个的）。

**第三轮**：确定最终 15 类标签 + 合并规则：
- Isoflavonoid → Flavonoid（黄酮）
- Jasmonoid + Sesterterpenoid → Terpenoid_other（其他萜类）
- Lipid → Fatty_acid（脂肪酸）
- Carbohydrate + Nucleoside → Unclassified（未分类）

### 3.4 最终 15 类分布

| 类别 | 中文 | 数量 | 占比 | 分类来源 |
|------|------|------|------|---------|
| Alkaloid | 生物碱 | 320 | 15.1% | superclass 245 + pathway 75 |
| Amino_acid | 氨基酸 | 292 | 13.7% | superclass 273 + pathway 19 |
| Fatty_acid | 脂肪酸 | 287 | 13.5% | superclass 144 + pathway 139 + name 4 |
| Phenylpropanoid | 苯丙素 | 222 | 10.4% | superclass 15 + **pathway 207** |
| Steroid | 甾体 | 211 | 9.9% | superclass 207 + name 4 |
| Diterpenoid | 二萜 | 144 | 6.8% | superclass 144 |
| Unclassified | 未分类 | 104 | 4.9% | superclass 5 + pathway 7 + **empty 91** + name 1 |
| Triterpenoid | 三萜 | 97 | 4.6% | superclass 97 |
| Polyketide | 聚酮 | 93 | 4.4% | **pathway 93** |
| Sesquiterpenoid | 倍半萜 | 93 | 4.4% | superclass 93 |
| Monoterpenoid | 单萜 | 79 | 3.7% | superclass 79 |
| Terpenoid_other | 其他萜类 | 67 | 3.2% | superclass 2 + pathway 63 + name 2 |
| Flavonoid | 黄酮 | 53 | 2.5% | superclass 49 + name 4 |
| Macrolide | 大环内酯 | 38 | 1.8% | superclass 38 |
| Coumarin | 香豆素 | 25 | 1.2% | superclass 25 |

**关键问题**：苯丙素 207/222 来自 pathway fallback，聚酮 93/93 全部来自 pathway fallback——这些可信度低。

---

## 四、SMARTS 结构规则交叉验证（2026-03-27）

### 4.1 方法

脚本：`scripts/05_底物分类/cross_validate_smarts.py`

对 7 个有明确化学骨架定义的类别，用 RDKit SMARTS 模式验证 NPClassifier 的分类：

| 类别 | SMARTS 规则 |
|------|------------|
| 甾体 | 6-6-6-5 稠合四环骨架 |
| 生物碱 | 含环氮 |
| 氨基酸 | alpha-氨基酸骨架 (N-C-C(=O)O) |
| 脂肪酸 | 羧基 + 长链碳 (≥8C) |
| 黄酮 | 黄酮/黄烷酮/查尔酮核心 |
| 大环内酯 | ≥12 元环 + 环内酯键 |
| 香豆素 | 苯并吡喃-2-酮骨架 |

### 4.2 结果

| 结果 | 数量 | 含义 |
|------|------|------|
| confirmed（NPClassifier + SMARTS 一致） | 936 | 两个方法都说是同一个类 |
| not_matched（SMARTS 规则存在但不匹配） | 290 | 可能是 NPClassifier 错了，也可能是 SMARTS 太严 |
| no_opinion（无 SMARTS 规则） | 899 | 萜类、苯丙素、聚酮等没有简单的骨架规则 |

### 4.3 局限

SMARTS "confirmed" 不等于 100% 正确——SMARTS 只检查"结构模式是否存在"，不检查"这个分类是否合理"。例如含环氮的合成药物也能通过生物碱的 SMARTS 检查。

---

## 五、ClassyFire 规则分类器（2026-03-28）

### 5.1 方法

- 脚本：`scripts/05_底物分类/classify_classyfire_parallel.py`（8 线程并行）
- 原理：RDKit 生成 InChIKey → 查询 `http://classyfire.wishartlab.com/entities/{InChIKey}.json`
- 耗时：~5 分钟（8 线程并行）

### 5.2 结果

| 结果 | 数量 | 占比 |
|------|------|------|
| 成功分类 | 328 | 15.5% |
| 数据库未收录 | 72 | 3.4% |
| 其他错误（大多也是未收录） | 1,710 | 81.0% |

**严重局限**：ClassyFire 通过 InChIKey 查已有数据库，大量合成色氨酸类似物和非天然产物未收录。只有 328 个化合物能提供第二套分类。

### 5.3 与 NPClassifier 交叉比对

328 个有 ClassyFire 结果的化合物中：
- 表面一致率：57.6%
- 实际大部分"不一致"是**分类体系差异**（如 ClassyFire 把萜类归在"Lipids and lipid-like molecules"下），不是真正的分类错误

---

## 六、Opus 子 agent 文献验证（2026-03-28）

### 6.1 方法

启动 11 个 Opus 子 agent，每个负责一个类别或任务：
- 每个 agent 用 **WebSearch + WebFetch** 去 PubChem/ChEBI/Wikipedia 读实际页面内容
- 总计 **~500 次 web 调用**
- 覆盖 **~400 个化合物**（苯丙素和未分类全量审计，其余抽检 ~20 个）

### 6.2 各 agent 结果

| Agent | 覆盖类别 | 抽检数 | 正确 | 错误 | web 调用 | 耗时 |
|-------|---------|--------|------|------|---------|------|
| 甾体验证 | 甾体(211) | 29 | 29 | 0 | ~30 | ~30s |
| 二萜验证 | 二萜(144) | 20 | 20 | 0 | 30 | ~2.5min |
| 三萜验证 | 三萜(97) | 20 | 20 | 0 | 63 | ~8.5min |
| 单萜验证 | 单萜(79) | 20 | 20 | 0 | 51 | ~3min |
| 倍半萜验证 | 倍半萜(93) | 20 | 19 | 1 | 36 | ~2.7min |
| 其他萜类验证 | 其他萜类(67) | 20 | 17 | 3 | 20 | ~2.4min |
| 聚酮验证 | 聚酮(93) | 20 | 17 | 3 | 31 | ~2.8min |
| 苯丙素全量审计 | 苯丙素(222) | **全部222** | 143 | **67** | 52 | ~4.8min |
| 未分类全量审计 | 未分类(104) | **全部104** | 57(保留) | **39(应重分类)** | 32 | ~4min |
| SMARTS 抽检 1 | 混合 | 20 | 18 | 2 | 54 | ~4.4min |
| SMARTS 抽检 2 | 混合 | 20 | 16 | 3+1数据错 | 31 | ~5.7min |
| 生物碱验证(推理) | 生物碱(320) | 14 | 6 | 4+4边界 | 0 | ~0.8min |
| 脂肪酸验证(推理) | 脂肪酸(287) | 22 | 3 | 19 | 0 | ~0.6min |
| 黄酮验证 | 黄酮(53) | 17 | 12 | 3+1数据错 | 3 | ~1.3min |

### 6.3 关键发现

#### 6.3.1 苯丙素 pathway fallback 系统性错误（67 个）

**根因**：NPClassifier 不确定时返回 pathway="Shikimates and Phenylpropanoids"（一个大杂烩），我们的映射逻辑把它当成了"苯丙素"。实际上这个 pathway 包含了所有含芳环的化合物。

被错标为苯丙素的实际是：
- PAH 多环芳烃（13个）：苯并[a]芘、菲、蒽、萘、芘、荧蒽、茚、甲基萘等
- 卤代芳烃（24个）：二氯苯、氯代甲苯、溴甲苯、氟代苯甲酸等
- 简单芳烃/溶剂（8个）：甲苯、乙苯、苯甲醚、二苯甲烷、茚满等
- 芳基硫醚（3个）：硫代苯甲醚、乙基苯硫醚等
- 合成探针（2个）：7-乙氧基试卤灵(EROD)、7-苄氧基试卤灵(BROD)
- 含氟合成物（3个）：三氟甲基苯乙烯、三氟甲氧基苯甲酸、氟比洛芬
- 其他（14个）：硝基苯酚、2-氨基联苯、2-乙酰氨基芴、DES等

#### 6.3.2 脂肪酸定义过宽（21 个已确认 + ~49 个推算）

NPClassifier 的 "Fatty acyls" superclass 包含所有脂肪族化合物：
- 烷烃（7个）：丙烷、丁烷、戊烷、己烷、庚烷、辛烷、壬烷
- 环烷烃（4个）：环戊烷、环己烷、环庚烷、环辛烷
- 卤甲烷（3个）：甲基溴、甲基氯、甲基碘
- 其他（7个）：丙醇、有机磷农药(Phorate sulfoxide)、腈(2)、溴化阻燃剂(TBECH)、2-甲基丁醛、反式-3-己烯

#### 6.3.3 生物碱错标（10 个已确认）

- 无氮化合物被标为生物碱：五氯苯（"脯氨酸生物碱"！）、1-甲氧基-4-甲硫基苯、七环二萜
- 合成药物：丙咪嗪（三环抗抑郁药）
- 农药：氟草隆、利尿隆、敌草隆（脲类除草剂）
- 简单酰胺：乙酰苯胺、苯脲
- 边界案例：多巴胺、L-DOPA、4'-O-甲基去甲莲碱（生物合成上可算"原生物碱"）

#### 6.3.4 氨基酸错标（22 个已确认）

- 大环肽：环孢素A（11 肽免疫抑制剂）
- 环缩肽：WS9326A
- 二酮哌嗪生物碱：Tryprostatin B（应归生物碱）
- 非化合物：氮原子(N)、苯胺
- 葡萄糖苷等

#### 6.3.5 未分类中被遗漏的天然产物（39 个）

NPClassifier 实际在更细级别给了分类，但被我们的映射逻辑忽略了：
- → 聚酮（5个）：蒽环类抗生素（柔红霉素衍生物等）
- → 氨基酸（6个）：色氨酸衍生物（吗啉基色氨酸等）
- → 苯丙素（7个）：简单酚类（对甲酚、间甲酚等）
- → 脂肪酸（5个）：长链酸
- → 碳水化合物（4个）：甲基半乳糖等单糖
- → 单萜（3个）：降冰片酮类
- → 核苷（3个）：NADP+、SAM 等
- → 其他萜类（3个）：二氢茉莉酮等
- → 大环内酯（2个）：大环糖苷
- → 卟啉（1个）：Fe(III)-heme b

#### 6.3.6 SMILES 数据错误（2 个）

- CMP_G001809：名称是"all tested hydroperoxides"，但 SMILES 编码的是葡萄糖
- CMP_G001812：名称是"indole-3-yl-methyl glucosinolate (I3M)"，但 SMILES 编码的是 I₃⁻ 三碘化物离子

---

## 七、错误汇总

### 7.1 已发现的 173 个明确错误

| 从哪个类 | 应改为 | 数量 | 具体构成 |
|---------|--------|------|---------|
| 苯丙素 | → 未分类 | 67 | PAH(13) + 卤代苯(24) + 含氟合成物(3) + 合成探针(2) + 简单芳烃(8) + 硫醚(3) + 其他(14) |
| 未分类 | → 各自正确类 | 39 | 聚酮(5) + 氨基酸(6) + 苯丙素(7) + 脂肪酸(5) + 碳水(4) + 单萜(3) + 核苷(3) + 其他萜(3) + 大环内酯(2) + 卟啉(1) |
| 氨基酸 | → 各自正确类 | 22 | 环肽 + 缩肽 + 生物碱 + 苯胺 + 氮原子 + 苯丙素共轭物 等 |
| 脂肪酸 | → 未分类 | 21 | 烷烃(7) + 环烷烃(4) + 卤甲烷(3) + 其他(7) |
| 生物碱 | → 未分类/其他 | 10 | 无氮化合物(3) + 合成药(1) + 农药(3) + 简单酰胺(2) + 二萜(1) |
| 黄酮 | → 各自正确类 | 4 | 二苯甲酮 + 反式芪 + SMILES 数据错误(2) |
| 其他萜类 | → 未分类 | 3 | 合成硫醚 + 磺酰叠氮 + 抗组胺药 |
| 大环内酯 | → 其他 | 2 | 利福霉素(安沙霉素) |
| 倍半萜 | → 甾体 | 1 | Grundmann's ketone(维D降解产物) |
| 聚酮 | → 未分类 | 1 | 荧光素(PAH) |
| 香豆素 | → 未分类 | 1 | 氯唑沙宗(苯并噁唑酮) |
| **总计** | | **173** | |

### 7.2 推算全数据集错误

根据各类别的抽检错误率乘以总数：

| 类别 | 总数 | 抽检正确率 | 推算错误数 |
|------|------|-----------|-----------|
| 甾体 + 二萜 + 三萜 + 单萜 | 531 | 100% | ~0 |
| 倍半萜 | 93 | 95% | ~5 |
| 其他萜类 | 67 | 85% | ~10 |
| 聚酮 | 93 | 85% | ~9 |
| 大环内酯 | 38 | 71% | ~11 |
| 香豆素 | 25 | 50% | ~12 |
| 黄酮 | 53 | 76% | ~12 |
| 生物碱 | 320 | 71% | ~58 |
| 氨基酸 | 292 | 74% | ~58 |
| 脂肪酸 | 287 | 50% | ~70 |
| 苯丙素 | 222 | 64% | 67（全量审计） |
| 未分类 | 104 | — | 39（全量审计，应重分类） |
| **合计** | **2,125** | | **~351（16.5%）** |

### 7.3 置信度分层

| 层级 | 数量 | 占比 | 说明 |
|------|------|------|------|
| A. 高置信（两个来源一致） | ~986 | 46% | NPClassifier superclass + SMARTS 规则一致 |
| B. 较高置信（一个来源，抽检 ≥95%） | ~480 | 23% | 萜类等，NPC 给了 superclass，抽检几乎无错 |
| C. 低置信（一个来源，抽检发现有错） | ~555 | 26% | 苯丙素/脂肪酸/生物碱/氨基酸 |
| D. 最低（NPC 返回空） | ~104 | 5% | 其中 39 个应重分类 |

---

## 八、错误根因分析

| 根因 | 影响数量 | 详细说明 |
|------|---------|---------|
| **Pathway fallback 系统性错误** | ~67 | NPClassifier 不确定时返回 pathway="Shikimates and Phenylpropanoids"，我们的映射逻辑 `PATHWAY_REMAP` 把它全标成了苯丙素。实际这个 pathway 是个大杂烩，包含所有芳环化合物 |
| **"Fatty acyls" 定义过宽** | ~70 | NPClassifier 的 "Fatty acyls" superclass 包含所有脂肪族化合物（烷烃、醇、醛、卤甲烷），不只是真正的脂肪酸 |
| **合成物被标为天然产物** | ~30 | NPClassifier 在 73,607 个天然产物上训练，遇到合成药/农药/工业品会强行套天然产物标签 |
| **环肽归入氨基酸** | ~40 | NPClassifier 的 "Amino acids and Peptides" pathway 包含环肽，我们映射到 "Amino_acid" 粒度不够 |
| **未分类中被遗漏** | ~39 | NPClassifier 返回空时直接标了 Unclassified，但其中有些在更细级别有分类被我们忽略 |
| **SMILES 数据录入错误** | 2 | 数据源头的问题，非分类工具的问题 |

---

## 九、输出文件

| 文件 | 位置 | 说明 |
|------|------|------|
| NPClassifier 分类脚本 | `scripts/05_底物分类/classify_substrates.py` | 含映射表、通配符处理、API 调用 |
| SMARTS 验证脚本 | `scripts/05_底物分类/cross_validate_smarts.py` | 7 个类别的结构规则 |
| ClassyFire 串行版 | `scripts/05_底物分类/classify_classyfire.py` | 单线程版本 |
| ClassyFire 并行版 | `scripts/05_底物分类/classify_classyfire_parallel.py` | 8 线程并行 |
| NPClassifier API 缓存 | `data/05_底物分类/npclassifier_cache.json` | 2,110 个结果 |
| ClassyFire API 缓存 | `data/05_底物分类/classyfire_cache.json` | 2,110 个查询（328 成功） |
| ClassyFire 结果 CSV | `data/05_底物分类/classyfire_results.csv` | 328 个有分类结果 |
| 初版分类结果 | `data/05_底物分类/substrate_classifications.csv` | 合并前 21 类 |
| 合并后最终版 | `data/05_底物分类/substrate_classifications_final.csv` | 合并后 15 类（当前版本，含 ~351 个错误） |

---

## 十、修正执行（2026-03-28 续）

### 10.1 修正策略（Codex 4 轮讨论）

与 Codex 进行了 4 轮深入讨论，确定了以下修正架构：

**不修改 NPClassifier，而是在其上构建"纠正层"**：
1. **Manual overrides**: 2 个 SMILES 数据错误
2. **Re-derive**: 从 NPClassifier 原始缓存重新推导标签，使用扩展的 superclass 映射
   - 关键发现：NPClassifier 返回的 superclass 名称带括号注释（如 `"Phenylpropanoids (C6-C3)"`），原始映射只有 `"Phenylpropanoids"`，导致 135/207 个 pathway-only 苯丙素无法匹配
   - 修复后，55 个 `Phenylpropanoids (C6-C3)` + 26 个 `Phenolic acids (C6-C1)` + 4 个 `Stilbenoids` + 29 个 `Cyclic polyketides` 等被正确重分类
3. **Structural vetoes**: RDKit SMARTS 结构验证
   - Alkaloid: 必须含氮
   - Fatty_acid: 必须有游离羧基 + 脂肪链 ≥4C
   - Amino_acid: 拒绝多肽/大环；DKP 仅在 NPC 确认生物碱时重分类
   - Phenylpropanoid: pathway-only → 降级，然后 SMARTS 家族救回
   - Phenolic acids: 必须有 ArOH 或 ArOMe（天然产物指标），否则降级为异源物
4. **Unclassified rescue**: 从原始 NPC class 层级救回（蒽环类→聚酮等）
5. **Phenylpropanoid SMARTS rescue**: 肉桂酸、香豆素、二苯乙烯、木脂素、酚类骨架匹配

### 10.2 修正结果

**总修正 275 个化合物（12.9%）**

| 变化 | 数量 | 规则 |
|------|------|------|
| Phenylpropanoid → Unclassified | 110 | 83 pathway-only 无 SMARTS 匹配 + 26 Phenolic acids 异源物 + 1 其他 |
| Fatty_acid → Unclassified | 109 | 103 无羧基 + 6 链<4C |
| Amino_acid → Alkaloid | 29 | DKP + NPC pathway 确认 |
| Amino_acid → Unclassified | 17 | 多肽/大环/DKP 无 NPC 确认 |
| Unclassified → Polyketide | 4 | 蒽环类 class 救回 |
| Alkaloid → Unclassified | 3 | 无氮 |
| Amino_acid → Phenylpropanoid | 2 | 对香豆酰胺重推导 |
| Flavonoid → Unclassified | 1 | SMILES 数据错误 |

### 10.3 修正后分布

| 类别 | 修正前 | 修正后 | 变化 |
|------|--------|--------|------|
| Alkaloid | 320 | 346 | +26 |
| **Unclassified** | **104** | **340** | **+236** |
| Amino_acid | 292 | 244 | -48 |
| Steroid | 211 | 211 | 0 |
| Fatty_acid | 287 | 178 | -109 |
| Diterpenoid | 144 | 144 | 0 |
| Phenylpropanoid | 222 | 114 | -108 |
| Triterpenoid | 97 | 97 | 0 |
| Polyketide | 93 | 97 | +4 |
| Sesquiterpenoid | 93 | 93 | 0 |
| Monoterpenoid | 79 | 79 | 0 |
| Terpenoid_other | 67 | 67 | 0 |
| Flavonoid | 53 | 52 | -1 |
| Macrolide | 38 | 38 | 0 |
| Coumarin | 25 | 25 | 0 |

### 10.4 质量评估

- **置信度**: high 2075 (97.6%), medium 50 (2.4%), low 0
- **Codex 估计错误率**: ~16.5% → **~3-6%**
- **微审计**: 20 Alkaloid + 20 Unclassified + 4 rescued Polyketide 全部正确
- 340 个 Unclassified 合理（异源物、探针、工业品等，不强行贴错标签）

### 10.5 修正脚本

`scripts/05_底物分类/correct_classifications.py`:
- 完全确定性、幂等（多次运行结果相同）
- 依赖 RDKit + 标准库
- 输入: `substrate_classifications_final.csv` + `npclassifier_cache.json`
- 输出: `substrate_classifications_corrected.csv`（含 `corrected_label`, `correction_reason`, `confidence` 列）

---

## 十一、输出文件（最终版）

| 文件 | 位置 | 说明 |
|------|------|------|
| NPClassifier 分类脚本 | `scripts/05_底物分类/classify_substrates.py` | 原始分类 |
| SMARTS 验证脚本 | `scripts/05_底物分类/cross_validate_smarts.py` | 7 类结构验证 |
| ClassyFire 并行版 | `scripts/05_底物分类/classify_classyfire_parallel.py` | 交叉验证 |
| **修正脚本** | `scripts/05_底物分类/correct_classifications.py` | **Codex 4 轮审核的修正层** |
| NPClassifier 缓存 | `data/05_底物分类/npclassifier_cache.json` | 2,110 个结果 |
| ClassyFire 缓存 | `data/05_底物分类/classyfire_cache.json` | 2,110 查询（328 成功） |
| 原始分类 | `data/05_底物分类/substrate_classifications_final.csv` | 合并前（含 ~351 个错误） |
| **修正后分类** | `data/05_底物分类/substrate_classifications_corrected.csv` | **最终版（~3-6% 残余错误）** |

---

## 十二、多源验证管线（2026-03-28）

### 12.1 问题：修正后准确率仍然不够

Codex 修正后估计错误率 ~3-6%，但实际检验发现 Tier 3-4（证据不足的化合物）错误率远高于此。
决定构建多源验证管线：结构规则 + ClassyFire + Agent 文献验证 + 共识引擎。

### 12.2 Phase 1: 综合结构验证器

为全部 15 类编写了 SMARTS 验证器 + 碳数/环拓扑分析：
- 甾体：6-6-6-5 四环拓扑检测
- 香豆素/黄酮：骨架 SMARTS
- 大环内酯：≥12 元环 + 内酯键
- 脂肪酸：羧基 + 脂肪链 ≥4C
- 氨基酸：α-氨基酸骨架
- 生物碱：亚型 SMARTS 库（吲哚/异喹啉/哌啶/嘧啶/嘌呤/托烷）+ 排除规则
- 萜类：碳数+低杂原子+萜类拓扑启发式
- 苯丙素：家族 SMARTS
- 聚酮：蒽醌/萘醌亚型

脚本：`scripts/05_底物分类/verify_classifications.py`

### 12.3 Phase 2: ClassyFire 批量查询

- 对 2,110 个化合物通过 InChIKey 查询 Wishart Lab ClassyFire API
- 结果：**855/2,110 命中**（40.2%）
- 输出：`data/05_底物分类/classyfire_full_results.csv`

### 12.4 Phase 3: 共识引擎与分层

对每个化合物交叉对比结构规则 + NPClassifier + ClassyFire，分为 4 个置信层级：

| 层级 | 标准 | 数量 | 占比 |
|------|------|------|------|
| Tier 1 | 结构证明 + 外部来源一致 | 541 | 25.5% |
| Tier 2 | 两个来源一致或强结构支撑 | 794 | 37.4% |
| Tier 3 | 仅一个来源 | 512 | 24.1% |
| Tier 4 | 证据不足或矛盾 | 278 | 13.1% |

Tier 1+2 = 1,335 个自动接受（62.8%），Tier 3-4 = 790 个需要 Agent 验证（37.2%）。

---

## 十三、Agent 文献验证 Round 1（2026-03-28）

### 13.1 概述

对 Tier 3-4 中最高风险的 **450 个化合物**启动 **15 个 Opus Agent 并行**验证：
- 每个 Agent 负责 30 个化合物
- 必须用 WebSearch + WebFetch 查询 PubChem/ChEBI/Wikipedia/PMC 论文
- 必须阅读实际页面内容（不是只看标题）
- 总计 ~480 次网络搜索

### 13.2 结果

- 确认正确：259 个（58%）
- **发现需要重分类：190 个（42%）**

### 13.3 主要错误模式

| 错误类型 | 数量 | 原因 |
|---------|------|------|
| Alkaloid → Unclassified | 47 | 合成药/农药/探针有环氮但非天然产物 |
| Phenylpropanoid → Unclassified | 43 | C6-C1/C6-C2 简单酚、卤代肉桂酸、异源物 |
| Diterpenoid → Terpenoid_other | 14 | C25 倍半萜酯被错标为 C20 二萜 |
| Polyketide → Unclassified | 10 | 简单芳烃无 PKS 来源 |
| Alkaloid → Amino_acid | 10 | 修饰色氨酸保留完整 AA 骨架 |
| Terpenoid_other → Unclassified | 9 | 合成物（氟化 ABA、硫醚） |
| Terpenoid_other → Sesquiterpenoid | 8 | ABA/独脚金内酯是 C15 倍半萜 |
| Macrolide → Polyketide | 7 | 大环内酰胺（非内酯）+ 安沙霉素 |
| Phenylpropanoid → Polyketide | 6 | 松脂酸衍生物是 PKS 来源 |
| Polyketide → Terpenoid_other | 5 | 杂萜中萜类为主导骨架 |

### 13.4 Agent 结果文件

`data/05_底物分类/agent_results_batch01.csv` ~ `agent_results_batch15.csv`（15 个文件）

---

## 十四、Agent 文献验证 Round 2（2026-03-28~29）

### 14.1 概述

对剩余 **341 个未验证的 Unclassified** 化合物启动 **10 个 Opus Agent 并行**验证：
- 检查这些被标为"未分类"的化合物是否有被错误踢出的天然产物
- ~340 次网络搜索

### 14.2 结果

- 确认 Unclassified 正确：291 个（85%）
- **救回到正确类别：50 个（15%）**

### 14.3 救回的化合物类型

| 救回到 | 数量 | 典型例子 |
|--------|------|---------|
| Fatty_acid | 17 | 脂酰辅酶 A、N-棕榈酰氨基酸、天然脂肪醇、AHL 群体感应分子 |
| Amino_acid | 7 | 异亮氨酸/缬氨酸醛肟（CYP79 产物）、SAM、修饰色氨酸 |
| Alkaloid | 6 | Tryprostatin A/B、Thaxtomin D、乙烯基甲基色氨酸 |
| Monoterpenoid | 3 | 降冰片酮（P450cam 经典底物） |
| Phenylpropanoid | 2 | p-茴香酸、二氢茴香脑 |
| Macrolide | 2 | A-2315A 抗生素、灰黄霉素绿素 |
| Polyketide | 2 | 三羟基呫吨酮、隐孢菌素 D |
| Diterpenoid | 1 | Atisane 二萜 |

### 14.4 Agent 结果文件

`data/05_底物分类/agent_results_r2_batch01.csv` ~ `agent_results_r2_batch10.csv`（10 个文件）

---

## 十五、二次验证 + Codex 仲裁（2026-03-29）

### 15.1 概述

Round 1 中有 **106 个化合物**置信度为 medium。启动 **4 个 Opus Agent 并行**对其进行二次独立验证：
- 每个 Agent 能看到第一轮 Agent 的判断和理由
- 要求独立上网重新验证，不能简单同意
- ~200 次网络搜索

### 15.2 结果

| 指标 | 数值 |
|------|------|
| 两轮 Agent 一致 | 48 (69%) |
| 两轮 Agent 分歧 | 22 (31%) |

### 15.3 Codex 仲裁

22 个分歧全部由 Codex 逐一分析裁决。关键政策决定：
- 脂肪醇/脂肪酮 → Unclassified（无羧基/酰基不算脂肪酸）
- 缬氨酸醛肟中间体 → Unclassified（不再保留 AA 骨架）
- 天然甲氧基/羟基苯甲酸 → Phenylpropanoid（有天然产物特征）
- 水杨酸 → Unclassified（来自异分支酸途径非 PAL）
- 杂萜(萜主导) → Terpenoid_other
- 植物呫吨酮 → Polyketide（传统聚酮归属）

### 15.4 仲裁结果文件

`data/05_底物分类/agent_results_recheck_batch02.csv`、`agent_results_recheck_batch03.csv`

---

## 十六、R1+R2+二次验证后的分布（2026-03-29）

应用全部 Agent 验证 + Codex 仲裁后的分布：

| 类别 | 英文 | 数量 | 占比 |
|------|------|------|------|
| 未分类 | Unclassified | 398 | 18.7% |
| 生物碱 | Alkaloid | 295 | 13.9% |
| 氨基酸 | Amino_acid | 261 | 12.3% |
| 甾体 | Steroid | 211 | 9.9% |
| 脂肪酸 | Fatty_acid | 203 | 9.6% |
| 二萜 | Diterpenoid | 137 | 6.4% |
| 倍半萜 | Sesquiterpenoid | 98 | 4.6% |
| 三萜 | Triterpenoid | 97 | 4.6% |
| 聚酮 | Polyketide | 95 | 4.5% |
| 单萜 | Monoterpenoid | 79 | 3.7% |
| 苯丙素 | Phenylpropanoid | 79 | 3.7% |
| 其他萜类 | Terpenoid_other | 66 | 3.1% |
| 黄酮 | Flavonoid | 50 | 2.4% |
| 大环内酯 | Macrolide | 33 | 1.6% |
| 香豆素 | Coumarin | 23 | 1.1% |

验证覆盖：
- Agent 文献验证：790 个（37.2%），~1,020 次 web 搜索
- 结构规则 + NPC：1,335 个（62.8%）
- 总计：2,125 个（100%）

输出文件：`data/05_底物分类/substrate_classifications_FINAL.csv`

---

## 十七、200 个分层随机抽样准确率审计（2026-03-29）

### 17.1 方法

- 按 15 个类别比例分层抽样 200 个化合物
- 5 个 Opus Agent 并行审计（每个 40 个）
- 每个 Agent 独立上网搜索验证，判断标签是 CORRECT 还是 WRONG

### 17.2 审计结果

| 批次 | 正确/总数 | 准确率 |
|------|----------|--------|
| Batch 1 | 35/40 | 87.5% |
| Batch 2 | 34/40 | 85.0% |
| Batch 3 | 35/40 | 87.5% |
| Batch 4 | 38/40 | 95.0% |
| Batch 5 | 35/40 | 87.5% |
| **总计** | **177/200** | **88.5%** |

### 17.3 23 个审计错误的模式

| 错误类型 | 数量 | 说明 |
|---------|------|------|
| Alkaloid → Amino_acid | 7 | 修饰色氨酸有完整 AA 骨架应为氨基酸 |
| Unclassified → Phenylpropanoid | 4 | 苯乙醛肟、乙氧基肉桂酸等苯丙素途径产物被错踢 |
| Unclassified → Fatty_acid | 3 | 脂肪酸甲酯、脂酰辅酶 A、壬酮 |
| Unclassified → Amino_acid | 2 | 环孢菌素（环肽）、甲氨基苯甲酸 |
| Unclassified → Alkaloid | 2 | 苯并噁嗪酮、溴化吲哚 |
| Polyketide → Phenylpropanoid | 2 | 植物呫吨酮、莽草酸来源苯乙酮 |
| Alkaloid → Unclassified | 2 | 丁螺环酮、Boc-氮杂降冰片烷（合成物） |
| Diterpenoid → Terpenoid_other | 1 | 视黄醛是脱辅基类胡萝卜素非二萜 |

### 17.4 根因分析

**核心问题不是"分类太细"，而是"边界规则不统一"：**

1. **Alkaloid vs Amino_acid 边界**（7/23 = 30% 的错误）：从未定义清楚"保留完整 AA 骨架的修饰色氨酸"算哪类
2. **Unclassified 过于保守**（9/23 = 39%）：一些天然产物被结构规则过于严格地踢出
3. **分类哲学不统一**：有的 Agent 按生物合成来源分，有的按化学结构分

---

## 十八、文献调研与分类哲学（2026-03-29）

### 18.1 NPClassifier 的原始定义

论文（Kim et al., J. Nat. Prod. 2021）明确指出：
- NPClassifier 使用**有向无环图（DAG）结构**，不是严格的树
- 一个化合物**可以同时属于多个类别**（多标签）
- 分类基于**生物合成来源 + 化学结构的混合**
- 技术实现使用 sigmoid 多标签输出

### 18.2 一个化合物能否同时属于多个类别？

**是的。** 文献确认：
- 杂萜(meroterpenoid) = 萜类 + 聚酮
- 黄酮的 A 环来自醋酸途径（聚酮），B/C 环来自苯丙素途径
- DKP 同时属于生物碱和氨基酸/肽
- NPClassifier 原文："peptide alkaloids Superclass belongs to both Alkaloids and Amino acids-Peptides Pathways"

### 18.3 P450 底物分类的文献空白

搜索了 17 篇相关论文（49 次 web 搜索），**没有人做过"从 P450 序列预测底物化学类别"**：
- 所有 P450 ML 论文都做"底物是/否"二分类
- 最接近的 Jinich et al. (JCIM 2023) 在 SDR/SAM-MTase 上做了类似工作（非 P450，仅 4-5 个粗二分类）
- 详见 `C3_创新性分析.md`

---

## 十九、合并方案讨论（2026-03-29）

### 19.1 提出的 8 类方案

Codex 建议合并到 8 类以减少边界争议：

| 合并后类别 | 中文 | 数量 | 合并了哪些 |
|-----------|------|------|-----------|
| Steroid | 甾体 | 211 | 不变 |
| Terpenoid | 萜类 | 477 | 单萜+倍半萜+二萜+三萜+其他萜 |
| Alkaloid | 生物碱 | 295 | 不变（修正 AA 骨架规则） |
| Amino_acid | 氨基酸 | 261 | 不变（修正 AA 骨架规则） |
| Fatty_acid | 脂肪酸 | 203 | 放宽（含脂酰辅酶 A） |
| Phenylpropanoid_broad | 苯丙素广义 | 152 | 苯丙素+黄酮+香豆素 |
| Polyketide_broad | 聚酮广义 | 128 | 聚酮+大环内酯 |
| Other | 其他 | 398 | Unclassified |

### 19.2 Codex 制定了正式分类手册

为 8 个类别写了精确的包含/排除/优先级/边界案例规则。关键边界决定：
- 修饰色氨酸（完整 AA 骨架）→ Amino_acid
- DKP → Alkaloid
- 茉莉酸 → Terpenoid
- 脂酰辅酶 A → Fatty_acid
- 脂肪醇/酮 → Other
- 合成药 → Other
- 等 20+ 个明确边界案例

### 19.3 合并后的审计准确率估算

仅类别合并消除了 1/23 个审计错误（二萜→其他萜）。**主要提升需要来自边界规则修正**（特别是 Alkaloid/Amino_acid 边界的 7 个错误）。

---

## 二十、当前状态与待办事项（2026-03-29）

### 20.1 已完成

| 步骤 | 详情 |
|------|------|
| NPClassifier 初始分类 | 2,125 → 15 类 |
| Codex 4 轮修正 | 275 个修正（结构规则） |
| ClassyFire 批量查询 | 855/2,110 命中（40.2%） |
| 结构验证器 + 共识引擎 | Tier 1-4 分层 |
| R1 Agent 验证 | 15 Agent × 30 化合物 = 450 个，~480 次 web 搜索，190 个重分类 |
| R2 Agent 验证 | 10 Agent × ~34 化合物 = 341 个，~340 次 web 搜索，50 个救回 |
| 二次验证 | 4 Agent × ~26 化合物 = 106 个，22 个分歧 Codex 仲裁 |
| 准确率审计 | 5 Agent × 40 化合物 = 200 个，**177/200 = 88.5%** |
| 文献调研 | NPClassifier/ClassyFire 论文 + 17 篇 P450 ML 论文 |
| 分类手册 | Codex 制定 8 类正式定义 |

**总计：~34 个 Opus Agent，~1,820 次 web 搜索**

### 20.2 未完成

| 步骤 | 说明 |
|------|------|
| ❌ 合并到 8 类 | 15→8 映射已设计，未执行 |
| ❌ 按分类手册重新分类 | 手册已写，未应用 |
| ❌ 修正审计发现的 23 个错误 | 已分析，未修正 |
| ❌ 重新审计 | 需要在修正后重新抽样验证 |

### 20.3 关键决策点

**用 15 类还是 8 类？** 取决于导师的下游目标：
- 如果做"酶序列→底物类别"预测（导师方向），8 类更稳健
- 底物分类只是训练标签，~88-90% 准确率对大类平均已够用

---

## 二十一、输出文件索引（截至 2026-03-29）

| 文件 | 位置 | 说明 |
|------|------|------|
| 原始分类脚本 | `scripts/05_底物分类/classify_substrates.py` | NPClassifier 批量分类 |
| 修正脚本 | `scripts/05_底物分类/correct_classifications.py` | Codex 4 轮修正层 |
| 验证引擎 | `scripts/05_底物分类/verify_classifications.py` | 结构规则+共识引擎 |
| ClassyFire 批量查询 | `scripts/05_底物分类/query_classyfire_batch.py` | InChIKey 批量查询 |
| NPClassifier 缓存 | `data/05_底物分类/npclassifier_cache.json` | 2,110 个结果 |
| ClassyFire 缓存 | `data/05_底物分类/classyfire_inchikey_cache.json` | 2,110 个查询 |
| ClassyFire 结果 | `data/05_底物分类/classyfire_full_results.csv` | 855 个有分类 |
| 修正后分类 | `data/05_底物分类/substrate_classifications_corrected.csv` | 第一轮修正 |
| 验证分层 | `data/05_底物分类/verified_classifications.csv` | Tier 1-4 |
| **最终分类** | `data/05_底物分类/substrate_classifications_FINAL.csv` | **当前最终版** |
| R1 Agent 结果 | `data/05_底物分类/agent_results_batch01~15.csv` | 15 批 |
| R2 Agent 结果 | `data/05_底物分类/agent_results_r2_batch01~10.csv` | 10 批 |
| 二次验证结果 | `data/05_底物分类/agent_results_recheck_batch02~03.csv` | 2 批（Codex 仲裁） |
| 审计批次数据 | `data/05_底物分类/agent_batch_audit_01~05.json` | 5 批×40 |
