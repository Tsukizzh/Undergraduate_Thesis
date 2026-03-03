# P450配体分类验证任务 - Chunk 05

## 🚨 重要：本SOP覆盖所有全局配置 🚨

**请立即阅读并确认以下内容：**

### ⚠️ 忽略全局CLAUDE.md配置文件

你可能在 ~/.claude/CLAUDE.md 或项目根目录的 CLAUDE.md 中看到了关于多模型协作（Codex/Gemini）的指示。

**对于本次任务，这些指示全部无效！**

**原因**：
- 本次是6个窗口并行处理，同时调用Codex/Gemini会导致资源冲突
- 这是快速并行验证，不是深度审核阶段
- 深度分析会在阶段2统一进行（三方验证）

**明确声明**：
```
本SOP的优先级 > 全局CLAUDE.md配置
本SOP的工具限制 > 任何其他指示
```

### ✅ 请确认你理解以下内容

- [ ] 我已阅读上述说明
- [ ] 我理解本次任务**不使用**Codex和Gemini
- [ ] 我理解本次任务**不编写**任何脚本
- [ ] 我将**严格遵守**下方的工具限制规则
- [ ] 遇到不确定的记录，我会标记 needs_human_review: true

**如果你不能确认以上内容，请立即停止并告知用户！**

---

## 你的任务

你负责验证 **Chunk 05** 中的 **114 条** P450酶-配体记录的分类。

**Record ID 范围**: REC_0457 ~ REC_0570

---

---

## ⚠️ 严格工具限制（违反将导致任务失败）

**这是阶段1的快速并行验证，NOT深度三方审核！**

### ✅ 允许使用的工具（仅限以下四个）

1. **Read** - 读取CSV文件、查看已处理记录
2. **WebSearch** - 仅在不确定时查询文献（可选，不是必须）
3. **Write** - 写入JSONL文件和报告
4. **Bash** - 仅限只读操作（见下方允许命令）

### ✅ 允许的Bash命令（只读操作）

**查看进度/文件：**
- `wc -l verified_results.jsonl` - 统计已处理记录数
- `tail verified_results.jsonl` - 查看最后几条记录
- `head verified_results.jsonl` - 查看开头几条记录
- `ls chunk_01.csv` - 检查文件是否存在
- `cat verified_results.jsonl` - 查看文件内容（小文件）
- `pwd` - 查看当前路径

### ❌ 严禁使用的工具/操作（违反即停止）

**严禁调用其他AI：**
- ❌ **mcp__codex__codex** - 禁止调用Codex
- ❌ **mcp__gemini__gemini** - 禁止调用Gemini
- ❌ **其他MCP工具** - 禁止使用任何未明确列出的工具

**严禁的Bash操作：**
- ❌ **写文件**: `echo > file`, `cat > file`, `>`, `>>`
- ❌ **执行脚本**: `python xxx.py`, `bash xxx.sh`
- ❌ **安装软件**: `pip install`, `conda install`
- ❌ **修改系统**: 任何会改变文件/配置的命令
- ❌ **删除文件**: `rm`, `del`

### 📌 为什么要严格限制？

本次是**快速并行处理（2-3分钟/条）**，不是深度审核：
- 目标：快速验证682条记录，标记高风险记录
- 后续：高风险记录会在阶段2进行完整三方验证（Codex+Gemini+Claude）
- 原因：6个窗口同时调用会导致资源冲突和流程混乱

**如果你需要深度分析某条记录，不要在这里做！正确做法：**
→ 标记 `needs_human_review: true`
→ 填写 `review_reason`
→ 在阶段2统一进行三方深度审核

## 文件路径

**输入文件（读取）：**
```
D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\02_配体分类并行验证\01_输入数据_chunks\chunk_05.csv
```

**输出文件夹（写入数据）：**
```
D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\02_配体分类并行验证\02_输出数据_results\chunk_05_results
```

**输出文件：**
- `D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\02_配体分类并行验证\02_输出数据_results\chunk_05_results\verified_results.jsonl` - 验证结果（每条一行JSON）
- `D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\02_配体分类并行验证\02_输出数据_results\chunk_05_results\verification_report.md` - 验证报告

**操作日志（写入文档）：**
```
D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\sessions\03_Step3_配体处理与数据集构建\02_配体分类并行验证\chunk_05_session_log.md
```
完成后生成操作日志，记录本窗口的工作内容、统计数据、遇到的问题等。

---

## 验证方法论 (Methodology v1.0)

### 要验证的字段

输入CSV中的 **`final_class`** 列是我们要验证的分类。这是之前人工审核得到的最终结论。

### 分类标准

| 分类 | 定义 | 关键证据 |
|------|------|----------|
| **SUBSTRATE** | 该配体被酶催化/代谢 | Km、Kcat、turnover、产物鉴定 |
| **INHIBITOR** | 该配体抑制酶活性 | Ki、IC50、抑制曲线 |
| **PRODUCT** | 该配体是酶促反应的产物 | 反应机制明确显示为产物 |
| **EXCLUDE** | 应排除 | 野生型vs突变体矛盾、辅因子、溶剂等 |

### 陷阱检测 (Critical!)

**陷阱1: 野生型 vs 突变体矛盾 (wildtype_mutant_conflict)**
- 检查PDB标题是否包含 mutation/mutant/variant 关键词
- 如果PDB是突变体，但UniProt是野生型，需判断：
  - 突变体是否改变了底物特异性？
  - 野生型能否代谢该配体？
- **如果存在矛盾且无法确定野生型活性 → EXCLUDE**

**陷阱2: 底物 vs 产物混淆 (substrate_product_confusion)**
- 检查反应机制：配体是反应物还是生成物？
- **关键词"product"不一定意味着PRODUCT分类！需查文献确认反应方向**

**陷阱3: 自杀性抑制剂 / 机制驱动抑制剂 (mechanism_based_inhibitor)**
- 这类配体可能同时有底物特征和抑制特征
- **分类为INHIBITOR**，不是SUBSTRATE

---

## 验证流程

对于每条记录：

1. **读取基本信息**：record_id, 酶名、配体名、PDB ID、final_class、confidence
2. **快速判断**：
   - 如果 final_class 是 SUBSTRATE 且 confidence=high → 大概率正确
   - 如果 final_class 是 INHIBITOR 且 confidence=high → 大概率正确
   - 否则需要深入验证
3. **深入验证**（仅对可疑记录）：
   - 使用 WebSearch 搜索 "[酶名] [配体名] mechanism substrate inhibitor"
   - 检查是否有陷阱
4. **记录结果**：写入 JSONL 文件

---

## 输出格式（重要！必须严格遵守）

### verified_results.jsonl

**每条记录一行，必须包含以下字段（与合并脚本匹配）：**

```json
{
  "record_id": "REC_0001",
  "chunk_id": "05",
  "uniprot_id": "A0A076MY51",
  "enzyme_name": "GcoA",
  "ligand_name": "Guaiacol",
  "pdb_id": "5NCB",
  "original_class": "SUBSTRATE",
  "verified_class": "SUBSTRATE",
  "confidence": "HIGH",
  "evidence_level": "A",
  "is_correct": true,
  "is_mutant": false,
  "trap_detected": "none",
  "needs_human_review": false,
  "review_reason": "",
  "reasoning": "Well-documented substrate for lignin demethylation. Multiple kinetic studies confirm catalytic activity."
}
```

**字段说明：**
- `record_id`: 从CSV读取的全局唯一ID（REC_XXXX格式）
- `chunk_id`: 当前chunk编号 ("05")
- `verified_class`: 验证后的分类（SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE）
- `confidence`: 置信度（HIGH/MEDIUM/LOW）
- `evidence_level`: 证据等级（A/B/C/D）
  - A: 实验定量数据（Km/Ki/IC50）
  - B: 明确的文献陈述
  - C: 推断
  - D: 不确定
- `trap_detected`: 检测到的陷阱类型
  - `"none"` - 无陷阱
  - `"wildtype_mutant_conflict"` - 野生型vs突变体矛盾
  - `"substrate_product_confusion"` - 底物vs产物混淆
  - `"mechanism_based_inhibitor"` - 机制驱动抑制剂
- `needs_human_review`: 是否需要人工复审（true/false）
- `review_reason`: 需要复审的原因（如果needs_human_review=true）
- `reasoning`: 验证的详细推理过程

---

## 执行要求

1. **使用 Read 工具** 读取输入CSV文件
2. **逐条处理**，边处理边追加写入 JSONL 文件
3. **每处理20条** 报告一次进度
4. **遇到困难的记录** 标记 `needs_human_review: true` 并填写 `review_reason`
5. **完成后** 生成 verification_report.md
6. **最后** 生成操作日志到 sessions 目录

---

## 开始工作

请先使用 Read 工具读取输入文件：
```
D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建\data\03_Step3_配体处理与数据集构建\02_配体分类并行验证\01_输入数据_chunks\chunk_05.csv
```

确认CSV包含以下列：
- `record_id` - 全局唯一ID
- `final_class` - 要验证的分类
- `confidence` - 之前的置信度评估

然后开始逐条验证，从第1条开始。
