# Chunk 05 文件索引

**任务**: Chunk 05 配体三方深度验证
**日期**: 2026-01-27
**记录范围**: REC_0229 - REC_0285 (共57条)

---

## 文件结构

```
04_配体三方深度验证/
├── 01_输入数据_chunks/
│   └── chunk_05.csv                    # 输入数据 (57条记录)
├── 02_输出数据_results/
│   └── chunk_05_results/
│       └── verified_results.jsonl     # 验证结果 (57条，JSONL格式)
└── 03_操作日志/
    ├── chunk_05_session_log.md        # 会话日志 (本次验证)
    └── chunk_05_file_index.md         # 本文件
```

---

## 文件详细说明

### 1. 输入文件

#### chunk_05.csv
- **路径**: `01_输入数据_chunks/chunk_05.csv`
- **格式**: CSV
- **记录数**: 57条
- **记录范围**: REC_0229 - REC_0285
- **字段**:
  - `record_id`: 记录ID
  - `enzyme_name`: UniProt ID
  - `ligand_name`: 配体名称
  - `pdb_id`: PDB编号
  - `verified_class`: 任务02分类结果
  - 其他元数据字段

### 2. 输出文件

#### verified_results.jsonl
- **路径**: `02_输出数据_results/chunk_05_results/verified_results.jsonl`
- **格式**: JSONL (每行一个JSON对象)
- **记录数**: 57条
- **文件大小**: ~45 KB
- **编码**: UTF-8

**JSONL结构**:
```json
{
  "record_id": "REC_0XXX",
  "chunk_id": "05",
  "enzyme_name": "P12345",
  "ligand_name": "compound_name",
  "pdb_id": "1ABC",
  "original_class": "UNKNOWN",
  "verified_class": "SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE",
  "final_class": "SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE",
  "classification_changed": true/false,
  "confidence": "HIGH/MEDIUM/LOW",
  "evidence_level": "A/B/C/D",
  "gemini_result": {
    "classification": "...",
    "data": "...",
    "pmid": "...",
    "summary": "..."
  },
  "codex_result": {
    "classification": "...",
    "fe_distance": "...",
    "binding_mode": "...",
    "summary": "..."
  },
  "consensus": true/false,
  "needs_human_review": true/false,
  "reasoning": "详细决策理由..."
}
```

**字段说明**:
- `classification_changed`: 是否纠正了任务02的分类
- `confidence`: 分类置信度 (HIGH/MEDIUM/LOW)
- `evidence_level`: 证据等级
  - A: 量化数据 (Ki, IC50, Km等)
  - B: 文献确认
  - C: 结构推断
  - D: 证据不足
- `consensus`: Gemini和Codex是否达成共识
- `needs_human_review`: 是否需要人工复核

### 3. 操作日志

#### chunk_05_session_log.md
- **路径**: `03_操作日志/chunk_05_session_log.md`
- **类型**: Markdown文档
- **内容**:
  - 处理统计（总数、分类分布、质量指标）
  - 发现的任务02错误详情
  - 典型案例展示（6个代表性案例）
  - 遇到的问题与解决方案
  - 数据质量评估
  - SOP执行情况
  - 后续改进建议

#### chunk_05_file_index.md
- **路径**: `03_操作日志/chunk_05_file_index.md`
- **类型**: Markdown文档
- **内容**: 本文件（文件索引和使用说明）

---

## 数据统计摘要

### 分类分布
| 分类 | 数量 | 占比 |
|------|------|------|
| SUBSTRATE | 27 | 47.4% |
| INHIBITOR | 23 | 40.4% |
| PRODUCT | 1 | 1.8% |
| EXCLUDE | 6 | 10.5% |
| **总计** | **57** | **100%** |

### 质量指标
- **任务02错误纠正**: 2条 (3.5%)
- **需人工审核**: 1条 (1.8%)
- **高置信度分类**: 54条 (94.7%)
- **证据等级A**: 42条 (73.7%)

### 涉及酶种类
- **CYP11A1**: 6条 (胆固醇侧链裂解酶)
- **CYP17A1**: 15条 (17α-羟化酶/17,20裂解酶)
- **CYP2D6**: 18条 (药物代谢酶)
- **CYP2A6**: 7条 (烟碱代谢酶)
- **CYP2C9**: 11条 (华法林代谢酶)

---

## 如何使用验证结果

### Python读取JSONL
```python
import json

# 读取所有记录
with open('verified_results.jsonl', 'r', encoding='utf-8') as f:
    data = [json.loads(line) for line in f]

# 筛选SUBSTRATE
substrates = [d for d in data if d['final_class'] == 'SUBSTRATE']

# 筛选高置信度记录
high_conf = [d for d in data if d['confidence'] == 'HIGH']

# 查找分类更改的记录
errors = [d for d in data if d['classification_changed']]
print(f"发现 {len(errors)} 条任务02错误")

# 需人工审核的记录
review = [d for d in data if d['needs_human_review']]
```

### 提取特定信息
```python
# 提取所有SUBSTRATE的PDB ID
substrate_pdbs = [d['pdb_id'] for d in data if d['final_class'] == 'SUBSTRATE']

# 统计各酶的底物数量
from collections import Counter
enzyme_counts = Counter(d['enzyme_name'] for d in data if d['final_class'] == 'SUBSTRATE')

# 查找特定记录
rec_0232 = next(d for d in data if d['record_id'] == 'REC_0232')
print(f"REC_0232分类: {rec_0232['final_class']}")
print(f"决策理由: {rec_0232['reasoning']}")
```

### 生成汇总报告
```python
# 统计分类分布
from collections import Counter
class_dist = Counter(d['final_class'] for d in data)

# 统计证据等级
evidence_dist = Counter(d['evidence_level'] for d in data)

# 统计置信度
conf_dist = Counter(d['confidence'] for d in data)

print("分类分布:", class_dist)
print("证据等级:", evidence_dist)
print("置信度:", conf_dist)
```

---

## 质量保证

### 数据完整性验证
```bash
# 检查记录数
wc -l verified_results.jsonl
# 预期输出: 57

# 验证JSON格式
python -c "import json; [json.loads(l) for l in open('verified_results.jsonl', encoding='utf-8')]"
# 无错误输出表示格式正确
```

### 必需字段验证
```python
required_fields = [
    'record_id', 'chunk_id', 'enzyme_name', 'ligand_name', 'pdb_id',
    'original_class', 'verified_class', 'final_class',
    'classification_changed', 'confidence', 'evidence_level',
    'gemini_result', 'codex_result', 'consensus',
    'needs_human_review', 'reasoning'
]

with open('verified_results.jsonl', 'r', encoding='utf-8') as f:
    data = [json.loads(line) for line in f]

for record in data:
    missing = [f for f in required_fields if f not in record]
    if missing:
        print(f"{record['record_id']} 缺少字段: {missing}")
```

---

## 版本历史

### v1.0 (2026-01-27)
- ✅ 完成57条记录的三方深度验证
- ✅ 发现并纠正2条任务02分类错误
- ✅ 标记1条需人工审核记录
- ✅ 生成完整JSONL输出和会话日志

---

## 相关文件引用

### 上游文件
- 任务02输出: `../../02_任务02_配体分类审核/02_输出/配体分类审核结果.xlsx`
- 原始PDB数据: `../../../01_Step1_下载PDB文件/data/pdb_files/`

### 下游文件
- 最终数据集: `../05_最终数据集/P450_test_dataset.csv` (待整合)
- 统计报告: `../06_统计分析/chunk_05_statistics.xlsx` (待生成)

---

## 技术支持

### 验证方法
- **文献检索**: Google Gemini (session-based)
- **结构分析**: OpenAI Codex (read-only mode)
- **决策汇总**: Claude Opus 4.5

### SOP文档
- `sop_chunk_05.md`: 详细验证流程
- `../00_SOP/三方验证标准操作流程.md`: 通用SOP

### 问题反馈
如发现数据问题或需要复核，请参考:
- `chunk_05_session_log.md` - 第六节"需人工复核的记录"
- `chunk_05_session_log.md` - 第四节"遇到的问题与解决方案"

---

**文档创建日期**: 2026-01-27
**最后更新**: 2026-01-27
**维护者**: Claude Opus 4.5 (三方协作验证系统)
