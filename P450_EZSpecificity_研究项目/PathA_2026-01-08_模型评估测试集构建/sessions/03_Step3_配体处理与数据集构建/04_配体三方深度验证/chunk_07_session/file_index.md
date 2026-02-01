# Chunk 07 File Index

## 输入文件

### 1. 核心输入数据
| 文件路径 | 说明 | 记录数 |
|---------|------|--------|
| `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_07.csv` | Task 02初步分类结果 | 57条 |

**字段说明**:
- record_id: REC_0343 - REC_0399
- uniprot_id: P14779, P15128, P15538, P18326, P20813, P20815
- enzyme_name: CYP102A1, CYP11B1, CYP105A1, CYP2B6, CYP3A5等
- ligand_name: 配体名称
- pdb_id: PDB结构ID
- original_class: 原始分类
- verified_class: Task 02分类
- reasoning: 初步reasoning

### 2. SOP指导文件
| 文件路径 | 说明 |
|---------|------|
| `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/sop_chunk_07.md` | 详细操作规范 (366行) |

**SOP核心内容**:
- 三方验证流程 (Gemini → Codex → Claude)
- 盲审原则 (blind review)
- 输出格式规范
- 证据等级定义 (A/B/C/D)

### 3. 参考脚本
| 文件路径 | 说明 |
|---------|------|
| `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_07_scripts/batch_process_remaining.py` | 剩余记录参考模板 |

## 输出文件

### 1. 主要输出
| 文件路径 | 说明 | 格式 |
|---------|------|------|
| `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_07_results/verified_results.jsonl` | 57条完整验证记录 | JSONL |

**记录结构**:
```json
{
  "record_id": "REC_XXXX",
  "chunk_id": "07",
  "enzyme_name": "...",
  "ligand_name": "...",
  "pdb_id": "...",
  "original_class": "...",
  "verified_class": "...",
  "final_class": "...",
  "classification_changed": true/false,
  "confidence": "HIGH/MEDIUM/LOW",
  "evidence_level": "A/B/C/D",
  "gemini_result": {...},
  "codex_result": {...},
  "consensus": true/false,
  "needs_human_review": true/false,
  "trap_detected": "...",
  "reasoning": "..."
}
```

### 2. Session文档
| 文件路径 | 说明 |
|---------|------|
| `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_07_session/session_log.md` | 本次验证总结 |
| `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_07_session/file_index.md` | 文件索引 (本文件) |

## PDB结构文件访问

所有PDB文件位于: `data/01_Step1_PDB文件/`

**涉及的PDB ID** (57个):
- CYP102A1 (27): 7D0U, 7D1F, 7E46, 7W97, 7W9J, 7WG0, 7WY1-3, 7XZK, 7Y0P-R, 7Y9J-L, 7YD9-L, 7YFT, 8DSG, 8JC3-4, 8YB2, 9ISU
- CYP2B6 (14): 3IBD, 3QOA, 3QU8, 3UA5, 4I91, 4RQL, 4RRT, 5UAP, 5UDA, 5UFG, 5WBG
- CYP105A1 (10): 3CV9, 5X7E, 9JHW, 9JI1, 9JKW, 9JKZ, 9KW2-4
- CYP11B1 (3): 6M7X, 7E7F
- CYP3A5 (3): 5VEU, 7LAD, 7SV2

## 数据统计

| 统计项 | 数值 |
|--------|------|
| 总记录数 | 57 |
| SUBSTRATE | 27 (47.4%) |
| INHIBITOR | 20 (35.1%) |
| PRODUCT | 3 (5.3%) |
| EXCLUDE | 7 (12.3%) |
| Task 02错误 | 5 (8.8%) |
| 需人工复审 | 1 (1.8%) |
