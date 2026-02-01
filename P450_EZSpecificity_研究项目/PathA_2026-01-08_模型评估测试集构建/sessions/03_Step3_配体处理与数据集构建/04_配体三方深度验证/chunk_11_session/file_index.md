# Chunk 11 验证文件索引

## 输入文件

### 主数据文件
- **chunk_11.csv**
  - 路径: `01_输入数据_chunks/chunk_11.csv`
  - 描述: 包含56条待验证的P450酶-配体记录（REC_0571-REC_0626）
  - 字段: record_id, chunk_id, uniprot_id, enzyme_name, ligand_name, pdb_id, original_class, verified_class等

### SOP文档
- **chunk_11_sop.md**
  - 路径: `01_输入数据_chunks/chunk_11_sop.md`
  - 描述: Chunk 11三方验证标准操作流程
  - 内容: Gemini→Codex→Claude三方协作验证步骤

## 输出文件

### 主输出
- **verified_results.jsonl**
  - 路径: `02_输出数据_results/chunk_11_results/verified_results.jsonl`
  - 格式: JSONL（每行一个JSON对象）
  - 记录数: 56条
  - 字段说明:
    - `record_id`: 记录唯一标识符
    - `enzyme_name`: 酶名称
    - `ligand_name`: 配体名称
    - `pdb_id`: PDB结构编号
    - `original_class`: 原始分类
    - `verified_class`: 任务02验证分类
    - `final_class`: 本次三方验证最终分类
    - `classification_changed`: 是否改变分类（true=发现错误）
    - `confidence`: 置信度（HIGH/MEDIUM/LOW）
    - `evidence_level`: 证据等级（A/B/C/D）
    - `gemini_result`: Gemini文献检索结果
    - `codex_result`: Codex结构分析结果
    - `consensus`: 三方是否达成共识
    - `needs_human_review`: 是否需要人工审核
    - `reasoning`: 综合推理与决策过程

### 文档
- **session_log.md**
  - 路径: `02_输出数据_results/chunk_11_results/session_log.md`
  - 描述: 本次验证的详细执行日志
  - 内容: 统计摘要、错误发现、特殊案例、执行时间线

- **file_index.md**
  - 路径: `02_输出数据_results/chunk_11_results/file_index.md`
  - 描述: 本文件，索引所有输入输出文件

## 错误发现详情

### REC_0572
- **文件**: verified_results.jsonl (第2行)
- **错误**: SUBSTRATE → INHIBITOR
- **酶**: CYP199A4
- **配体**: 4-benzoylbenzoic acid
- **证据**: PMID:26053303，诱导Type II光谱位移

### REC_0583
- **文件**: verified_results.jsonl (第13行)
- **错误**: SUBSTRATE → PRODUCT
- **酶**: CYP199A4
- **配体**: P-HYDROXYBENZOIC ACID
- **证据**: Kd=458µM，是4-methoxybenzoic acid的O-去甲基化产物

### REC_0623
- **文件**: verified_results.jsonl (第53行)
- **错误**: SUBSTRATE → EXCLUDE
- **酶**: P450_Q70AS3
- **配体**: Bisphenol A
- **证据**: Fe-O=11.28Å，远离活性位点

### REC_0624
- **文件**: verified_results.jsonl (第54行)
- **错误**: SUBSTRATE → INHIBITOR
- **酶**: P450_Q70KH6
- **配体**: Cyclopropyl pyrimidine
- **证据**: Fe-N=2.40Å，Type II抑制剂

## 统计数据

### 按酶分类
| 酶类型 | 记录数 | 文件行号范围 |
|--------|--------|--------------|
| CYP199A4 | 16 | 1-16 |
| CYP51 (T. brucei) | 8 | 17-24 |
| CYP74, CYP154E1, CYP119 | 7 | 25-31 |
| MycG | 5 | 31-36 |
| OleP | 6 | 36-42 |
| CYP120A1 | 1 | 42 |
| CYP154C5 | 7 | 43-49 |
| CYP51 (M. capsulatus) | 3 | 49-51 |
| 其他 | 3 | 52-56 |

### 按分类结果
| 最终分类 | 记录数 | 示例记录ID |
|----------|--------|------------|
| SUBSTRATE | 31 | REC_0571, REC_0573, REC_0579等 |
| INHIBITOR | 16 | REC_0572, REC_0587, REC_0596等 |
| PRODUCT | 4 | REC_0575, REC_0577, REC_0583, REC_0584 |
| EXCLUDE | 5 | REC_0574, REC_0576, REC_0601, REC_0604, REC_0625 |

## 质量指标文件

### 高置信度记录 (47条)
- 占比: 83.9%
- 特征: 文献+结构双重支持，或经典案例

### 中等置信度记录 (9条)
- 占比: 16.1%
- 记录ID: REC_0571, REC_0575, REC_0576, REC_0582, REC_0590, REC_0602, REC_0605, REC_0609-0611
- 原因: 文献缺失、结构模糊或突变体特异性

### 需人工审核记录 (1条)
- REC_0590: Lanosterol analog（底物vs产物争议）

## 文件使用说明

### 查询特定记录
```bash
# 查询REC_0572
grep '"record_id":"REC_0572"' verified_results.jsonl

# 查询所有分类改变的记录
grep '"classification_changed":true' verified_results.jsonl
```

### 统计分析
```bash
# 统计各分类数量
grep -o '"final_class":"[^"]*"' verified_results.jsonl | sort | uniq -c

# 统计置信度分布
grep -o '"confidence":"[^"]*"' verified_results.jsonl | sort | uniq -c
```

### 导入Python分析
```python
import json

records = []
with open('verified_results.jsonl', 'r') as f:
    for line in f:
        records.append(json.loads(line))

# 分析错误率
changed = sum(1 for r in records if r['classification_changed'])
print(f"错误率: {changed}/{len(records)} = {changed/len(records)*100:.1f}%")
```

## 版本信息
- **验证日期**: 2026-01-27
- **执行者**: Claude Code (Opus 4.5)
- **协作模型**: Gemini 2.5 Flash, Codex
- **SOP版本**: v1.0
