# Chunk 06 验证结果 - 文件索引

## 目录结构

```
chunk_06_results/
├── verified_results.jsonl   # 主验证结果文件
├── session_log.md           # 操作日志
└── file_index.md            # 本索引文件
```

## 文件说明

### 1. verified_results.jsonl

**文件类型**: JSON Lines格式 (每行一个JSON对象)
**记录数量**: 57条
**编码**: UTF-8
**生成时间**: 2026-01-27

#### 字段说明

| 字段名 | 类型 | 说明 |
|-------|------|------|
| record_id | string | 唯一记录标识 (REC_0286-REC_0342) |
| chunk_id | string | 分块编号 ("06") |
| uniprot_id | string | UniProt蛋白ID |
| enzyme_name | string | 酶名称 |
| ligand_name | string | 配体化学名称 |
| pdb_id | string | PDB结构ID |
| original_class | string | Task 02原始分类 |
| verified_class | string | Task 02验证后分类 |
| final_class | string | 本次三方验证最终分类 |
| classification_changed | boolean | 分类是否发生变更 |
| confidence | string | 置信度 (HIGH/MEDIUM/LOW) |
| evidence_level | string | 证据等级 (A/B/C) |
| gemini_result | object | Gemini文献检索结果 |
| codex_result | object | Codex结构分析结果 |
| consensus | boolean | 三方是否达成共识 |
| needs_human_review | boolean | 是否需要人工复核 |
| reasoning | string | 最终判定理由 |

#### 分类结果分布

| 分类 | 数量 | 说明 |
|-----|------|------|
| EXCLUDE | 29 | 非催化配体(辅因子/溶剂/decoy等) |
| INHIBITOR | 18 | 抑制剂(Type I/II) |
| SUBSTRATE | 10 | 底物 |
| PRODUCT | 0 | 产物 |

#### 使用示例

```python
import json

results = []
with open('verified_results.jsonl', 'r', encoding='utf-8') as f:
    for line in f:
        results.append(json.loads(line.strip()))

# 统计分类变更
changed = [r for r in results if r['classification_changed']]
print(f"Task 02 errors found: {len(changed)}")  # 输出: 15
```

### 2. session_log.md

**文件类型**: Markdown文档
**主要内容**:
- 任务概览和三方协作分工
- 验证统计和分类分布
- **15条Task 02错误的详细分析**
- Decoy分子识别规则总结
- Type II抑制剂判定标准
- 特殊情况记录
- 数据质量评估
- 下一步建议

### 3. file_index.md

**文件类型**: Markdown文档
**主要内容**: 本索引文件，描述目录结构和文件用途

## 关键发现速查

### Task 02错误记录 (15条)

| Record ID | 错误分类 | 正确分类 | 原因 |
|-----------|---------|---------|------|
| REC_0297 | SUBSTRATE | EXCLUDE | Decoy (perfluorinated) |
| REC_0303 | SUBSTRATE | EXCLUDE | Decoy (perfluorinated) |
| REC_0304 | SUBSTRATE | EXCLUDE | Decoy (perfluorinated) |
| REC_0305 | SUBSTRATE | EXCLUDE | Decoy (perfluorinated) |
| REC_0306 | SUBSTRATE | EXCLUDE | Decoy (perfluorinated) |
| REC_0314 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0315 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0319 | SUBSTRATE | EXCLUDE | Wrong protein (FABP4) |
| REC_0322 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0327 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0328 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0330 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0331 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0335 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |
| REC_0337 | SUBSTRATE | EXCLUDE | Decoy (Stanfield) |

### 需要人工复核

- **REC_0319**: PDB 6JLV可能是FABP4而非BM3，需确认源数据

## 数据版本控制

| 版本 | 日期 | 变更说明 |
|-----|------|---------|
| 1.0 | 2026-01-27 | 初始版本，完成57条记录三方验证 |

---

**索引生成**: 2026-01-27
