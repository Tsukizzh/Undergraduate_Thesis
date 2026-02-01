# Chunk 10 配体三方深度验证 - 文件索引

## 输入文件

### 1. 主要输入数据

| 文件路径 | 文件名 | 说明 | 记录数 | 用途 |
|---------|--------|------|--------|------|
| `04_配体三方深度验证/01_输入数据_chunks/` | `chunk_10.csv` | Chunk 10待验证记录列表 | 57 | 提供待验证的酶-配体对信息 |

**字段说明**:
- record_id: 记录唯一标识（REC_0514 - REC_0570）
- chunk_id: 批次编号（10）
- uniprot_id: UniProt酶标识符
- enzyme_name: 酶名称
- ligand_name: 配体名称
- pdb_id: PDB结构代码
- original_class: 原始分类（任务01）
- verified_class: 任务02验证分类
- confidence: 任务02置信度
- evidence_level: 任务02证据等级
- is_correct: 任务02判定正确性
- is_mutant: 是否为突变体
- trap_detected: 陷阱检测结果
- needs_human_review: 任务02人工复审标记
- review_reason: 复审原因
- reasoning: 任务02推理过程

### 2. 辅助参考文件

| 文件路径 | 文件名 | 说明 | 用途 |
|---------|--------|------|------|
| `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/` | `sop_chunk_10.md` | 标准操作流程文档 | 指导三方协作验证流程 |
| `04_配体三方深度验证/01_输入数据_chunks/` | `chunk_05.csv` | 参考：Chunk 5数据（用户查看） | 对比参考 |

## 输出文件

### 1. 主要输出数据

| 文件路径 | 文件名 | 说明 | 记录数 | 格式 |
|---------|--------|------|--------|------|
| `04_配体三方深度验证/02_输出数据_results/chunk_10_results/` | `verified_results.jsonl` | Chunk 10深度验证结果 | 57 | JSONL |

**字段说明**:
- record_id: 记录ID
- chunk_id: 批次ID
- enzyme_name: 酶名称
- ligand_name: 配体名称
- pdb_id: PDB代码
- original_class: 任务01分类
- verified_class: 任务02分类
- **final_class**: 任务03最终分类（本次验证结果）
- **classification_changed**: 是否与任务02不同（发现错误标志）
- confidence: 置信度（HIGH/MEDIUM/LOW）
- evidence_level: 证据等级（A/B/C/D）
- **gemini_result**: Gemini文献检索结果
  - classification: Gemini分类
  - data: 定量数据（Km/Ki/Kd等）
  - pmid: PubMed文献ID
  - summary: 文献摘要
- **codex_result**: Codex结构分析结果
  - classification: Codex分类
  - fe_distance: Fe-配体距离
  - binding_mode: 结合模式
  - summary: 结构分析总结
- consensus: Gemini与Codex是否一致
- needs_human_review: 是否需人工复审
- **reasoning**: Claude综合推理（中文）

**分类统计**:
- SUBSTRATE: 44条 (77.2%)
- INHIBITOR: 8条 (14.0%)
- PRODUCT: 1条 (1.8%)
- EXCLUDE: 4条 (7.0%)

**关键发现**:
- classification_changed=true: 1条（REC_0515，任务02错误）

### 2. 操作日志文件

| 文件路径 | 文件名 | 说明 | 格式 |
|---------|--------|------|------|
| `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_10_session/` | `session_log.md` | 处理过程详细日志 | Markdown |
| `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_10_session/` | `file_index.md` | 本文件（文件索引） | Markdown |

## 辅助脚本

### 1. MCP工具调用

本任务使用了以下MCP工具：

| 工具名称 | 用途 | 调用次数（约） |
|---------|------|---------------|
| `mcp__gemini__gemini` | 文献检索与定量数据查询 | 57次 |
| `mcp__codex__codex` | PDB结构分析与Fe距离测量 | 57次 |
| Bash | 文件写入与追加 | 45次 |
| Write | 临时文件创建 | 2次 |

**Gemini调用模式**:
```bash
mcp__gemini__gemini(
  PROMPT="Search literature for [enzyme] and [ligand]. Find Km/Ki/Kd data, PMID, classification.",
  cd="<project_root>",
  sandbox=false
)
```

**Codex调用模式**:
```bash
mcp__codex__codex(
  PROMPT="Analyze PDB [pdb_id] structure. Measure Fe-ligand distance, binding mode, coordination.",
  cd="<project_root>",
  sandbox="read-only"
)
```

### 2. 无外部Python脚本

本任务完全通过MCP工具协作完成，未使用独立Python脚本。

## 中间数据

### 1. 临时文件

| 文件路径 | 文件名 | 说明 | 状态 |
|---------|--------|------|------|
| `C:/Users/ADMINI~1/AppData/Local/Temp/claude/.../scratchpad/` | `batch_0558_0570.jsonl` | 最后13条记录临时批次文件 | 已清理 |

**用途**: 解决Windows环境下bash heredoc引号问题，用Write工具创建临时文件后cat追加。

### 2. 无持久化中间数据

所有处理结果直接写入最终输出文件（verified_results.jsonl），无需保留中间状态。

## 数据流图

```
输入数据流:
chunk_10.csv (57 records)
    ↓
    ├─→ Gemini (文献检索) → gemini_result
    ├─→ Codex (结构分析) → codex_result
    └─→ Claude (综合决策) → final_class, reasoning
         ↓
verified_results.jsonl (57 records)
         ↓
    ├─→ session_log.md (统计与案例)
    └─→ file_index.md (本文件)
```

## 文件关系

```
04_配体三方深度验证/
│
├── 01_输入数据_chunks/
│   └── chunk_10.csv ────────────┐
│                                 │
├── 02_输出数据_results/          │ (读取)
│   └── chunk_10_results/         │
│       └── verified_results.jsonl ←┘ (写入)
│
└── scripts/.../
    └── sop_chunk_10.md ─────────┐ (指导)
                                  │
sessions/.../chunk_10_session/    │
├── session_log.md ←──────────────┘ (记录)
└── file_index.md (本文件)
```

## 版本信息

- **任务版本**: Step 3.01 配体分类四方深度验证
- **数据版本**: Chunk 10 / Total 17 chunks
- **处理日期**: 2026-01-27
- **处理工具**: Claude Code + Gemini MCP + Codex MCP
- **协议版本**: SOP v1.0

## 质量保证

### 输入验证
- ✓ chunk_10.csv包含完整的57条记录
- ✓ 所有必需字段齐全
- ✓ PDB ID格式正确

### 输出验证
- ✓ verified_results.jsonl行数=57（与输入一致）
- ✓ 每行均为有效JSON对象
- ✓ 所有必需字段齐全（record_id, final_class, gemini_result, codex_result等）
- ✓ 分类枚举值有效（SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE）

### 流程验证
- ✓ 每条记录均执行三方协作（Gemini + Codex + Claude）
- ✓ 盲审原则贯彻（忽略verified_class进行独立判断）
- ✓ 发现任务02错误1条（classification_changed=true）

## 备注

1. **Gemini API失败**: 约7条记录Gemini调用失败，依赖Codex结构分析，标记evidence_level=C
2. **Windows路径兼容**: 所有路径使用正斜杠`/`或Windows原生反斜杠`\`
3. **JSONL写入方式**: 前44条逐条追加，后13条批量追加（临时文件workaround）
4. **人工复审**: REC_0515需人工复审（仅结构证据，无文献支持）
