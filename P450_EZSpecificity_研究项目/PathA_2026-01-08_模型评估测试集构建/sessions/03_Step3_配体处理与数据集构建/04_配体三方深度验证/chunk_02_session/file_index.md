# Chunk 02三方深度验证文件索引

## 会话基本信息

| 属性 | 值 |
|------|------|
| **Chunk ID** | 02 |
| **记录范围** | REC_0058 - REC_0114 (57条) |
| **处理日期** | 2026-01-27 |
| **状态** | ✅ 已完成 |

---

## 输入文件

### 主输入数据

| 文件名 | 路径 | 说明 |
|--------|------|------|
| `chunk_02.csv` | `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/` | Chunk 02原始数据,包含57条待验证记录 |

**文件格式**: CSV
**字段列表**:
- `record_id`: 记录唯一标识 (REC_0058-REC_0114)
- `chunk_id`: Chunk编号 (02)
- `uniprot_id`: UniProt ID
- `enzyme_name`: 酶名称
- `ligand_name`: 配体名称
- `pdb_id`: PDB结构ID
- `original_class`: PDB中的原始分类
- `verified_class`: 任务02验证分类
- `confidence`: 任务02置信度
- `evidence_level`: 任务02证据级别
- `is_correct`: 任务02正确性标记
- `is_mutant`: 是否突变体
- `trap_detected`: 陷阱检测
- `needs_human_review`: 是否需人工审核
- `review_reason`: 审核原因
- `reasoning`: 任务02推理

---

## 输出文件

### 主输出数据

| 文件名 | 路径 | 说明 |
|--------|------|------|
| `verified_results.jsonl` | `data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_02_results/` | 三方验证完整结果 |

**文件格式**: JSONL (每行一个JSON对象)
**记录数**: 57条

**字段说明**:
```json
{
  "record_id": "记录ID (REC_XXXX)",
  "chunk_id": "Chunk编号",
  "uniprot_id": "UniProt ID",
  "enzyme_name": "酶名称",
  "ligand_name": "配体名称",
  "pdb_id": "PDB ID",
  "original_class": "PDB原始分类",
  "verified_class": "任务02验证分类",
  "final_class": "三方验证最终分类",
  "classification_changed": "是否与任务02不同 (true/false)",
  "confidence": "置信度 (HIGH/MEDIUM/LOW)",
  "evidence_level": "证据级别 (A/B/C/D)",
  "gemini_result": {
    "classification": "Gemini分类结果",
    "data": "定量数据 (Ki, IC50等)",
    "pmid": "PubMed ID",
    "summary": "Gemini总结"
  },
  "codex_result": {
    "classification": "Codex分类结果",
    "fe_distance": "Fe-配体距离",
    "binding_mode": "结合模式",
    "summary": "Codex总结"
  },
  "consensus": "Gemini+Codex是否达成共识 (true/false)",
  "needs_human_review": "是否需人工审核 (true/false)",
  "reasoning": "Claude最终推理"
}
```

---

## 会话日志

| 文件名 | 路径 | 说明 |
|--------|------|------|
| `session_log.md` | `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_session/` | 操作日志 (本次任务的完整记录) |
| `file_index.md` | `sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_session/` | 文件索引 (本文件) |

---

## 辅助脚本

| 文件名 | 路径 | 说明 | 状态 |
|--------|------|------|------|
| `batch_verify.py` | `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_scripts/` | 批量验证框架脚本 | 未执行 (仅模板) |
| `append_records_0094_0103.py` | `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_scripts/` | 追加REC_0094-0103记录 | ✅ 已执行 |
| `append_records_0104_0114.py` | `scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_02_scripts/` | 追加REC_0104-0114记录 | ✅ 已执行 |

**脚本说明**:
- `batch_verify.py`: 展示批量处理逻辑的框架脚本,包含Gemini/Codex调用逻辑和共识决策机制
- `append_records_*.py`: 由于Bash heredoc语法问题,改用Python脚本批量追加记录

---

## 目录结构

```
PathA_2026-01-08_模型评估测试集构建/
├── data/
│   └── 03_Step3_配体处理与数据集构建/
│       └── 04_配体三方深度验证/
│           ├── 01_输入数据_chunks/
│           │   └── chunk_02.csv                    # 输入: 57条待验证记录
│           └── 02_输出数据_results/
│               └── chunk_02_results/
│                   └── verified_results.jsonl       # 输出: 验证结果
├── scripts/
│   └── 03_Step3_配体处理与数据集构建/
│       └── 04_配体三方深度验证/
│           └── chunk_02_scripts/
│               ├── batch_verify.py                 # 批量验证框架
│               ├── append_records_0094_0103.py     # 追加脚本1
│               └── append_records_0104_0114.py     # 追加脚本2
└── sessions/
    └── 03_Step3_配体处理与数据集构建/
        └── 04_配体三方深度验证/
            └── chunk_02_session/
                ├── session_log.md                  # 操作日志
                └── file_index.md                   # 文件索引(本文件)
```

---

## 统计摘要

### 输出文件统计

| 分类 | 数量 | 占比 |
|------|------|------|
| SUBSTRATE | 30 | 52.6% |
| INHIBITOR | 18 | 31.6% |
| EXCLUDE | 8 | 14.0% |
| PRODUCT | 1 | 1.8% |

### 质量指标

| 指标 | 数值 |
|------|------|
| 总记录数 | 57 |
| 三方共识率 | 91.2% (52/57) |
| 发现任务02错误 | 4条 |
| 需人工审核 | 2条 |
| A级证据占比 | 86.0% (49/57) |

---

## 重要记录索引

### 任务02错误记录 (classification_changed=true)

| record_id | 任务02分类 | 三方验证分类 | 说明 |
|-----------|-----------|-------------|------|
| REC_0094 | PRODUCT | INHIBITOR | VT-1161 tetrazole是CYP51抑制剂 |
| REC_0100 | PRODUCT | EXCLUDE | Bis-Tris是结晶缓冲剂 |
| REC_0102 | PRODUCT | SUBSTRATE | Farnesol是CYP124A1底物 |
| **REC_0112** | **PRODUCT** | **SUBSTRATE** | **Tirandamycin C是TamI的底物(被迭代氧化)** |

### 需人工审核记录 (needs_human_review=true)

| record_id | 最终分类 | 审核原因 |
|-----------|---------|---------|
| REC_0066 | SUBSTRATE | Gemini与历史记录冲突 |
| REC_0077 | EXCLUDE | 文献与结构分析冲突 |

---

**创建时间**: 2026-01-27
**创建者**: Claude Code
**版本**: v1.0
