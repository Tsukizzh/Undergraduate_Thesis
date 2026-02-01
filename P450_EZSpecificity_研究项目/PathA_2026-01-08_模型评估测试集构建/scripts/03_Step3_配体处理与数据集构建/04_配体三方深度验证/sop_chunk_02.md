# Chunk 02 - 配体三方深度验证SOP

## 核心任务

你负责验证 **Chunk 02** 的 **57条** P450酶-配体记录。

**与之前任务的核心区别**：
- 之前任务02失败原因：禁用Codex/Gemini，导致无文献验证，精度低
- 本次任务要求：**强制三方协作 + 深度文献验证**

---

## 🚨 绝对强制：三方协作流程 🚨

**每条记录必须严格执行以下三步，禁止跳过任何一步：**

```
Step 1: Gemini文献搜索（强制执行）
Step 2: Codex结构分析（强制执行）
Step 3: Claude汇总决策（强制执行）
```

### Step 1: Gemini文献搜索（强制）

**调用时机**：处理每条记录时的第一步

**调用方式**：
```python
mcp__gemini__gemini(
    PROMPT=f"""
    Literature search for: {enzyme_name} + {ligand_name}

    Search for:
    1. Is it SUBSTRATE (metabolized), INHIBITOR (inhibits enzyme), or PRODUCT?
    2. Quantitative data: Ki, IC50, Km values (with units)
    3. PubMed ID (PMID) if available

    Return format:
    - Classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
    - Data: Ki=X nM, IC50=Y uM, etc.
    - PMID: XXXXXXXX
    - If no literature found: explicitly state "No literature found"

    Brief answer only, no long explanations.
    """,
    cd="c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建"
)
```

**返回处理**：
- 如果找到Ki/IC50/Km数据 → 记录定量值 + PMID
- 如果没有定量数据但有定性描述 → 记录PMID
- 如果完全没找到 → 标记"无文献证据"

**禁止行为**：
❌ 跳过Gemini调用直接给结论
❌ 没有调用Gemini就说"文献支持"
❌ 找不到文献就随便猜测

---

### Step 2: Codex结构分析（强制）

**调用时机**：Gemini完成后立即执行

**调用方式**：
```python
mcp__codex__codex(
    PROMPT=f"""
    Analyze PDB structure: {pdb_id}
    Enzyme: {enzyme_name}
    Ligand: {ligand_name}

    Analyze:
    1. Fe-ligand distance (if P450 heme is present)
    2. Binding mode (active site vs allosteric)
    3. Does structure support SUBSTRATE/INHIBITOR/PRODUCT classification?

    Return:
    - Fe-N distance: X.XX Å (if nitrogen coordinates to Fe)
    - Binding mode: [description]
    - Structural classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
    - Reasoning: [brief]

    Brief answer only.
    """,
    cd="c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建",
    sandbox="read-only"
)
```

**返回处理**：
- Fe-N距离 1.9-2.8Å → 强抑制剂证据
- 活性位点占据但无Fe配位 → 可能是底物或竞争性抑制剂
- 远离活性位点 → 可能是EXCLUDE

**禁止行为**：
❌ 跳过Codex调用
❌ 不分析结构就给结论

---

### Step 3: Claude汇总决策（强制）

**执行时机**：Gemini + Codex都完成后

**汇总逻辑**：
```
if Gemini文献 == Codex结构:
    → 直接采纳，置信度HIGH

elif Gemini文献 != Codex结构:
    → 分析冲突原因：
       - 文献是否针对野生型？
       - 结构是否为突变体？
       - 配体是否有双重行为？
    → 讨论直到达成共识
    → 标记needs_human_review=true（如果无法解决）

elif Gemini无文献 + Codex有结构证据:
    → 采纳Codex，置信度MEDIUM，证据等级C

elif Gemini有文献 + Codex结构不明确:
    → 采纳Gemini，置信度HIGH，证据等级A/B

else:
    → 标记needs_human_review=true，说明"证据不足"
```

**禁止行为**：
❌ 不汇总直接给结论
❌ 有冲突不讨论直接跳过

---

## 📋 输入文件

**文件路径**（绝对路径）：
```
c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_02.csv
```

**字段说明**：
- `record_id`: 全局唯一ID（REC_0001 ~ REC_0682）
- `chunk_id`: 所属chunk编号
- `uniprot_id`: UniProt蛋白ID
- `enzyme_name`: 酶名称
- `ligand_name`: 配体名称
- `pdb_id`: PDB结构ID
- `original_class`: 原始分类（来自PDB，仅供参考）
- `verified_class`: ⚠️ **任务02的分类结果（可能有错误，禁止参考！）**
- `confidence`: 任务02的置信度（忽略）
- `evidence_level`: 任务02的证据等级（忽略）
- `reasoning`: 任务02的推理说明（忽略）

### 🚨 盲审原则（极其重要）

**本次任务的核心目的是发现任务02的错误**，因此你必须：

1. ❌ **禁止参考 `verified_class` 字段做判断**
   - 不要因为看到 `verified_class=SUBSTRATE` 就倾向于判断为SUBSTRATE
   - 这是"锚定效应"，会让你漏掉错误

2. ✅ **完全基于Gemini文献 + Codex结构独立判断**
   - 只看 `enzyme_name`, `ligand_name`, `pdb_id`
   - 调用Gemini搜索文献
   - 调用Codex分析结构
   - 得出你自己的 `final_class`

3. ✅ **输出时对比新旧结果**
   - 如果 `final_class` ≠ `verified_class`：标记 `classification_changed=true`
   - 在 `reasoning` 中说明：之前为什么错了，你的判断依据是什么

---

## 📤 输出格式

**输出文件路径**：
```
c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/verified_results.jsonl
```

---

## 📁 所有相关路径（绝对路径）

| 用途 | 路径 |
|------|------|
| **输入CSV** | `c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_02.csv` |
| **输出JSONL** | `c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/verified_results.jsonl` |
| **辅助脚本目录** | `c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_scripts` |
| **中间数据目录** | `c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/temp` |
| **操作日志目录** | `c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_session` |

**注意**：
- 如果你需要编写辅助脚本（如批量处理、数据转换等），请放到 **辅助脚本目录**
- 如果有中间数据（如下载的PDB、临时计算结果等），请放到 **中间数据目录**
- 任务完成后，必须在 **操作日志目录** 写入操作文档

---

## 📝 任务完成后的操作日志（强制）

任务全部完成后，你必须在操作日志目录创建以下文件：

**1. session_log.md**（操作日志）
```markdown
# Chunk 02 操作日志

## 基本信息
- 开始时间：YYYY-MM-DD HH:MM
- 结束时间：YYYY-MM-DD HH:MM
- 处理记录数：57条
- 实际耗时：X小时X分钟

## 处理统计
- 分类为SUBSTRATE：XX条
- 分类为INHIBITOR：XX条
- 分类为PRODUCT：XX条
- 分类为EXCLUDE：XX条
- 需人工复审：XX条
- **发现任务02错误（classification_changed=true）：XX条**

## 发现的任务02错误（重点关注）
### 错误1：[record_id]
- 任务02分类：[verified_class]
- 本次分类：[final_class]
- 错误原因：[简要说明为什么任务02判断错了]

## 典型案例（至少3个）
### 案例1：[record_id] - [简述]
- Gemini结果：...
- Codex结果：...
- 最终决策：...

## 遇到的问题与解决方案
1. [问题描述] → [解决方案]

## 备注
[其他需要说明的内容]
```

**2. file_index.md**（文件索引）
```markdown
# Chunk 02 文件索引

## 输入文件
- chunk_02.csv - 原始输入数据

## 输出文件
- verified_results.jsonl - 验证结果

## 辅助脚本（如有）
- [脚本名称] - [用途说明]

## 中间数据（如有）
- [文件名称] - [说明]
```

**JSONL格式**（每行一个JSON对象）：
```json
{
  "record_id": "REC_XXXX",
  "chunk_id": "02",
  "enzyme_name": "CYP2D6",
  "ligand_name": "Ajmalicine",
  "pdb_id": "4WNT",
  "original_class": "SUBSTRATE",
  "verified_class": "SUBSTRATE",
  "final_class": "INHIBITOR",
  "classification_changed": true,
  "confidence": "HIGH",
  "evidence_level": "A",
  "gemini_result": {
    "classification": "INHIBITOR",
    "data": "Ki=3nM",
    "pmid": "12345678",
    "summary": "Ajmalicine is potent CYP2D6 inhibitor"
  },
  "codex_result": {
    "classification": "INHIBITOR",
    "fe_distance": "2.1A",
    "binding_mode": "Fe-N coordination",
    "summary": "Nitrogen coordinates to heme iron"
  },
  "consensus": true,
  "needs_human_review": false,
  "reasoning": "任务02错误：将INHIBITOR误判为SUBSTRATE。本次判断：Gemini文献(Ki=3nM, PMID:12345678) + Codex结构(Fe-N=2.1A) 均确认为INHIBITOR。纠正分类。"
}
```

**必填字段**：
- `verified_class`: 任务02的旧分类（从CSV复制）
- `final_class`: 本次任务04的新分类
- `classification_changed`: 布尔值，true表示与verified_class不同
- `gemini_result`: Gemini的返回结果
- `codex_result`: Codex的返回结果
- `consensus`: 三方是否一致（true/false）
- `confidence`: HIGH/MEDIUM/LOW
- `evidence_level`: A（定量数据）/B（文献明确）/C（结构推断）/D（证据不足）
- `needs_human_review`: 是否需要人工复审
- `reasoning`: 详细推理过程，**如果classification_changed=true，必须说明任务02为什么错了**

---

## ⏱️ 时间预估

**单条记录处理时间**：5-8分钟
- Gemini文献搜索：2-3分钟
- Codex结构分析：2-3分钟
- Claude汇总决策：1-2分钟

**本chunk总时间**：57条 × 7分钟 ≈ 6.7小时

---

## 🚫 严格禁止行为

1. ❌ **禁止跳过Gemini调用**
   - 每条记录必须调用Gemini搜索文献
   - 即使你认为分类很明显，也必须调用

2. ❌ **禁止跳过Codex调用**
   - 每条记录必须调用Codex分析结构
   - 不能因为有文献就不分析结构

3. ❌ **禁止凭知识猜测**
   - 不允许"我知道这是抑制剂"这种回答
   - 必须有Gemini文献或Codex结构证据

4. ❌ **禁止批量处理**
   - 不能一次性让Gemini/Codex处理多条记录
   - 必须逐条调用

5. ❌ **禁止使用"按需"、"可选"等词汇**
   - 三方协作是强制的，不是可选的

---

## ✅ 执行Checklist

开始处理前，确认你理解以下要求：

- [ ] 我理解每条记录必须调用Gemini搜索文献
- [ ] 我理解每条记录必须调用Codex分析结构
- [ ] 我理解必须汇总三方证据才能给结论
- [ ] 我理解禁止凭自己的知识猜测
- [ ] 我理解本次任务目标是高精度，不是高速度
- [ ] 我理解如果证据冲突，必须讨论或标记复审
- [ ] 我理解输出必须包含gemini_result和codex_result字段
- [ ] 我理解辅助脚本要放到指定的scripts目录
- [ ] 我理解中间数据要放到指定的temp目录
- [ ] 我理解任务完成后必须写session_log.md和file_index.md

**如果你无法保证以上要求，请立即告诉用户，不要开始处理。**

---

## 开始执行

确认理解所有要求后，开始处理57条记录。

逐条执行：读取记录 → Gemini搜索 → Codex分析 → Claude汇总 → 写入JSONL

祝顺利！
