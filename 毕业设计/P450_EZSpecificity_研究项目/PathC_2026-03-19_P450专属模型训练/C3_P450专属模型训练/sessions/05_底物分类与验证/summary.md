# Step 5 底物分类与验证：详细数据记录

> **用途**: 存放 session_log 中省略的详细数据
> **日期**: 2026-03-27 ~ 2026-03-29
> **最终状态**: ✅ 完成。2,125 化合物 → 7+1 类多标签分类，50 抽检 96% 准确率
> **最终文件**: `data/05_底物分类/substrate_multilabel_FINAL.csv`

---

## 一、NPClassifier 初始 15 类分布（修正前）

| 类别 | 中文 | 数量 | 占比 | 分类来源 |
|------|------|------|------|---------|
| Alkaloid | 生物碱 | 320 | 15.1% | superclass 245 + pathway 75 |
| Amino_acid | 氨基酸 | 292 | 13.7% | superclass 273 + pathway 19 |
| Fatty_acid | 脂肪酸 | 287 | 13.5% | superclass 144 + pathway 139 + name 4 |
| Phenylpropanoid | 苯丙素 | 222 | 10.4% | superclass 15 + **pathway 207**（大量 fallback） |
| Steroid | 甾体 | 211 | 9.9% | superclass 207 + name 4 |
| Diterpenoid | 二萜 | 144 | 6.8% | superclass 144 |
| Unclassified | 未分类 | 104 | 4.9% | superclass 5 + pathway 7 + empty 91 + name 1 |
| Triterpenoid | 三萜 | 97 | 4.6% | superclass 97 |
| Polyketide | 聚酮 | 93 | 4.4% | **pathway 93**（全部来自 fallback） |
| Sesquiterpenoid | 倍半萜 | 93 | 4.4% | superclass 93 |
| Monoterpenoid | 单萜 | 79 | 3.7% | superclass 79 |
| Terpenoid_other | 其他萜类 | 67 | 3.2% | superclass 2 + pathway 63 + name 2 |
| Flavonoid | 黄酮 | 53 | 2.5% | superclass 49 + name 4 |
| Macrolide | 大环内酯 | 38 | 1.8% | superclass 38 |
| Coumarin | 香豆素 | 25 | 1.2% | superclass 25 |

**关键问题**：苯丙素 207/222 来自 pathway fallback，聚酮 93/93 全部来自 pathway fallback。

---

## 二、Codex 修正详情（275 个化合物）

### 2.1 修正规则

1. **Manual overrides**：2 个 SMILES 数据错误
2. **Re-derive**：从 NPClassifier 原始缓存重新推导标签
   - 关键 bug：NPClassifier 返回 `"Phenylpropanoids (C6-C3)"`，映射表只有 `"Phenylpropanoids"` → 135 个匹配失败
   - 修复 `normalize_npc_name()` 函数：去掉括号注释
3. **Structural vetoes**：
   - 生物碱必须含氮
   - 脂肪酸必须有羧基 + 脂肪链≥4C
   - 苯丙素 pathway-only → 降级后 SMARTS 救回
   - Phenolic acids 需有 ArOH/ArOMe（天然产物指标）
4. **Unclassified rescue**：从 NPC class 层级救回（蒽环→聚酮等）
5. **Phenylpropanoid SMARTS rescue**：肉桂酸/香豆素/二苯乙烯/木脂素骨架

### 2.2 修正数量

| 变化 | 数量 | 规则 |
|------|------|------|
| Phenylpropanoid → Unclassified | 110 | pathway-only 无 SMARTS + Phenolic acids 异源物 |
| Fatty_acid → Unclassified | 109 | 无羧基(103) + 链<4C(6) |
| Amino_acid → Alkaloid | 29 | DKP + NPC pathway 确认 |
| Amino_acid → Unclassified | 17 | 多肽/大环/DKP 无 NPC 确认 |
| Unclassified → Polyketide | 4 | 蒽环类 class 救回 |
| Alkaloid → Unclassified | 3 | 无氮 |
| Amino_acid → Phenylpropanoid | 2 | 对香豆酰胺重推导 |
| Flavonoid → Unclassified | 1 | SMILES 数据错误 |

---

## 三、Agent 文献验证详细结果

### 3.1 阶段1 Opus Agent 审计（修正前，3-27~28）

| Agent | 覆盖类别 | 抽检数 | 正确 | 错误 | web调用 |
|-------|---------|--------|------|------|---------|
| 甾体验证 | Steroid(211) | 29 | 29 | 0 | ~30 |
| 二萜验证 | Diterpenoid(144) | 20 | 20 | 0 | 30 |
| 三萜验证 | Triterpenoid(97) | 20 | 20 | 0 | 63 |
| 单萜验证 | Monoterpenoid(79) | 20 | 20 | 0 | 51 |
| 倍半萜验证 | Sesquiterpenoid(93) | 20 | 19 | 1 | 36 |
| 其他萜类验证 | Terpenoid_other(67) | 20 | 17 | 3 | 20 |
| 聚酮验证 | Polyketide(93) | 20 | 17 | 3 | 31 |
| 苯丙素全量审计 | Phenylpropanoid(222) | 222 | 143 | 67 | 52 |
| 未分类全量审计 | Unclassified(104) | 104 | 57 | 39 | 32 |

### 3.2 R1 Agent 验证（Tier 3-4，15 Agent × 30 化合物）

190 个重分类的主要模式：

| 错误类型 | 数量 | 原因 |
|---------|------|------|
| Alkaloid → Unclassified | 47 | 合成药/农药/探针有环氮但非天然产物 |
| Phenylpropanoid → Unclassified | 43 | C6-C1/C2 简单酚、卤代肉桂酸 |
| Diterpenoid → Terpenoid_other | 14 | C25 倍半萜酯被错标为 C20 二萜 |
| Polyketide → Unclassified | 10 | 简单芳烃无 PKS 来源 |
| Alkaloid → Amino_acid | 10 | 修饰色氨酸保留完整 AA 骨架 |
| Terpenoid_other → Unclassified | 9 | 合成物 |
| Terpenoid_other → Sesquiterpenoid | 8 | ABA/独脚金内酯是 C15 |
| Macrolide → Polyketide | 7 | 大环内酰胺 + 安沙霉素 |
| Phenylpropanoid → Polyketide | 6 | 松脂酸衍生物是 PKS 来源 |
| Polyketide → Terpenoid_other | 5 | 杂萜中萜类为主导 |

### 3.3 R2 Agent 验证（Unclassified，10 Agent × ~34 化合物）

50 个救回到正确类别：

| 救回到 | 数量 | 典型例子 |
|--------|------|---------|
| Fatty_acid | 17 | 脂酰辅酶A、AHL群体感应分子 |
| Amino_acid | 7 | 缬氨酸醛肟（CYP79产物）、SAM |
| Alkaloid | 6 | Tryprostatin A/B |
| Monoterpenoid | 3 | 降冰片酮（P450cam经典底物） |
| Phenylpropanoid | 2 | p-茴香酸 |
| Macrolide | 2 | A-2315A抗生素 |
| Polyketide | 2 | 隐孢菌素D |
| Diterpenoid | 1 | Atisane二萜 |

### 3.4 Codex 仲裁（22 个分歧案例）

关键政策决定：
- 脂肪醇/脂肪酮 → Unclassified（无羧基/酰基不算脂肪酸）
- 缬氨酸醛肟 → Unclassified（不再保留AA骨架）
- 天然甲氧基/羟基苯甲酸 → Phenylpropanoid
- 水杨酸 → Unclassified（来自异分支酸途径非PAL）
- 杂萜(萜主导) → Terpenoid_other
- 植物呫吨酮 → Polyketide

---

## 四、已发现的 173 个明确错误（修正前全量统计）

| 从哪个类 | 应改为 | 数量 | 构成 |
|---------|--------|------|------|
| 苯丙素 | → 未分类 | 67 | PAH(13)+卤代苯(24)+含氟合成物(3)+合成探针(2)+简单芳烃(8)+硫醚(3)+其他(14) |
| 未分类 | → 正确类 | 39 | 聚酮(5)+氨基酸(6)+苯丙素(7)+脂肪酸(5)+碳水(4)+单萜(3)+核苷(3)+其他萜(3)+大环内酯(2)+卟啉(1) |
| 氨基酸 | → 正确类 | 22 | 环肽+缩肽+DKP生物碱+苯胺+N原子+共轭物 |
| 脂肪酸 | → 未分类 | 21 | 烷烃(7)+环烷烃(4)+卤甲烷(3)+其他(7) |
| 生物碱 | → 未分类/其他 | 10 | 无氮化合物(3)+合成药(1)+农药(3)+简单酰胺(2)+二萜(1) |
| 黄酮 | → 正确类 | 4 | 二苯甲酮+反式芪+SMILES数据错误(2) |
| 其他萜类 | → 未分类 | 3 | 合成硫醚+磺酰叠氮+抗组胺药 |
| 大环内酯 | → 其他 | 2 | 利福霉素(安沙霉素) |
| 倍半萜 | → 甾体 | 1 | Grundmann's ketone |
| 聚酮 | → 未分类 | 1 | 荧光素(PAH) |
| 香豆素 | → 未分类 | 1 | 氯唑沙宗 |

---

## 五、错误根因分析

| 根因 | 影响数量 | 说明 |
|------|---------|------|
| Pathway fallback 系统性错误 | ~67 | NPC 不确定时返回 pathway="Shikimates and Phenylpropanoids"，映射逻辑全标成苯丙素 |
| "Fatty acyls" 定义过宽 | ~70 | NPC superclass 含所有脂肪族（烷烃/醇/醛/卤甲烷） |
| 合成物被标为天然产物 | ~30 | NPC 在天然产物上训练，遇到合成物强行套标签 |
| 环肽归入氨基酸 | ~40 | NPC "Amino acids and Peptides" pathway 含环肽 |
| 未分类中被遗漏 | ~39 | NPC 返回空时直接标 Unclassified，但细级别有分类被忽略 |
| SMILES 数据录入错误 | 2 | CMP_G001809(SMILES是葡萄糖非对应名称) + CMP_G001812(SMILES是I₃⁻) |

---

## 六、200 样本审计分批结果

| 批次 | 正确/总数 | 准确率 | Agent |
|------|----------|--------|-------|
| Batch 1 | 35/40 | 87.5% | Opus Agent 1 |
| Batch 2 | 34/40 | 85.0% | Opus Agent 2 |
| Batch 3 | 35/40 | 87.5% | Opus Agent 3 |
| Batch 4 | 38/40 | 95.0% | Opus Agent 4 |
| Batch 5 | 35/40 | 87.5% | Opus Agent 5 |
| **总计** | **177/200** | **88.5%** | |

---

## 七、每类的验证来源分布

| 类别 | structural_npc | agent_literature | double_checked | codex_arbitrated | high | medium |
|------|---------------|-----------------|----------------|-----------------|------|--------|
| Unclassified(398) | 0 | 395 | 3 | 0 | 362 | 35 |
| Alkaloid(295) | 265 | 29 | 1 | 0 | 261 | 34 |
| Amino_acid(261) | 231 | 30 | 0 | 0 | 253 | 8 |
| Steroid(211) | 207 | 4 | 0 | 0 | 211 | 0 |
| Fatty_acid(203) | 174 | 29 | 0 | 0 | 188 | 15 |
| Diterpenoid(137) | 124 | 13 | 0 | 0 | 135 | 2 |
| Sesquiterpenoid(98) | 88 | 10 | 0 | 0 | 98 | 0 |
| Triterpenoid(97) | 95 | 2 | 0 | 0 | 97 | 0 |
| Polyketide(95) | 0 | 93 | 2 | 0 | 81 | 14 |
| Monoterpenoid(79) | 62 | 17 | 0 | 0 | 77 | 2 |
| Phenylpropanoid(79) | 0 | 69 | 6 | 4 | 65 | 14 |
| Terpenoid_other(66) | 2 | 61 | 2 | 1 | 56 | 10 |
| Flavonoid(50) | 33 | 17 | 0 | 0 | 50 | 0 |
| Macrolide(33) | 31 | 2 | 0 | 0 | 33 | 0 |
| Coumarin(23) | 23 | 0 | 0 | 0 | 23 | 0 |

---

## 八、NPClassifier 多标签化合物统计

207/2,110 个化合物（9.8%）返回了多个 pathway。按组合分布：

通过 `npclassifier_cache.json` 检查可得具体组合（待统计完整分布）。

---

## 九、工作量统计

| 项目 | 数量 |
|------|------|
| Opus Agent | ~34 个 |
| Web 搜索 | ~1,820 次 |
| Codex 讨论 | ~8 轮 |
| 覆盖率 | 2,125/2,125 = 100% |

---

## 十、文件索引

### 脚本

| 文件 | 说明 |
|------|------|
| `scripts/05_底物分类/classify_substrates.py` | NPClassifier 批量分类（原始） |
| `scripts/05_底物分类/correct_classifications.py` | Codex 4 轮修正层 |
| `scripts/05_底物分类/verify_classifications.py` | 结构规则+共识引擎 |
| `scripts/05_底物分类/query_classyfire_batch.py` | ClassyFire InChIKey 批量查询 |

### 数据

| 文件 | 说明 |
|------|------|
| `data/05_底物分类/npclassifier_cache.json` | NPClassifier 2,110 个结果缓存 |
| `data/05_底物分类/classyfire_inchikey_cache.json` | ClassyFire 2,110 查询缓存 |
| `data/05_底物分类/classyfire_full_results.csv` | ClassyFire 855 个有分类结果 |
| `data/05_底物分类/substrate_classifications.csv` | 初版分类（合并前 21 类） |
| `data/05_底物分类/substrate_classifications_corrected.csv` | Codex 修正后（15 类） |
| `data/05_底物分类/verified_classifications.csv` | 共识引擎 Tier 1-4 |
| **`data/05_底物分类/substrate_classifications_final.csv`** | **当前最终版（15 类，含全部 Agent 验证结果）** |
| `data/05_底物分类/agent_results_batch01~15.csv` | R1 Agent 结果（15 批） |
| `data/05_底物分类/agent_results_r2_batch01~10.csv` | R2 Agent 结果（10 批） |
| `data/05_底物分类/agent_results_recheck_batch02~03.csv` | 二次验证结果（Codex 仲裁） |
| `data/05_底物分类/agent_batch_audit_01~05.json` | 200 样本审计批次数据 |
