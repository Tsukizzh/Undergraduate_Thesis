# Chunk 01 文件索引

## 输入文件

### 主输入
- **chunk_01.csv** - 原始输入数据（57条P450酶-配体记录）
  - 路径：`../../data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_01.csv`
  - 格式：CSV，含header
  - 字段：record_id, chunk_id, uniprot_id, enzyme_name, ligand_name, pdb_id, original_class, verified_class（任务02结果）, confidence, evidence_level, reasoning等

## 输出文件

### 主输出
- **verified_results.jsonl** - 配体三方深度验证结果
  - 路径：`../../data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_01_results/verified_results.jsonl`
  - 格式：JSONL（每行一个JSON对象）
  - 记录数：57条
  - 文件大小：~90KB
  - 字段：
    - record_id, chunk_id
    - enzyme_name, ligand_name, pdb_id
    - original_class, verified_class（任务02）, **final_class（本次结果）**
    - **classification_changed**（bool，是否与任务02不同）
    - confidence（HIGH/MEDIUM/LOW）
    - evidence_level（A/B/C/D）
    - gemini_result（Gemini文献搜索结果）
    - codex_result（Codex结构分析结果）
    - consensus（三方是否一致）
    - needs_human_review（是否需人工复审）
    - reasoning（详细推理过程）

## 辅助脚本

本次任务未生成额外脚本（所有操作通过MCP工具直接完成）

## 中间数据

无中间数据（直接从MCP工具获取结果）

## 操作日志

### 日志文件
- **session_log.md** - 详细操作日志
  - 路径：`./session_log.md`
  - 内容：处理统计、发现的错误、典型案例、问题与解决方案等

- **file_index.md** - 本文件
  - 路径：`./file_index.md`
  - 内容：文件索引

## 文件结构树

```
chunk_01_session/
├── session_log.md          # 操作日志
└── file_index.md           # 本文件

相关数据文件（在其他目录）:
../../data/.../01_输入数据_chunks/
└── chunk_01.csv            # 输入数据

../../data/.../02_输出数据_results/chunk_01_results/
├── verified_results.jsonl  # 验证结果
└── temp/                   # 中间数据目录（未使用）

../../scripts/.../chunk_01_scripts/
└── (无脚本生成)            # 辅助脚本目录（未使用）
```

## 关键统计

| 项目 | 数值 |
|------|------|
| 输入记录数 | 57条 |
| 输出记录数 | 57条 |
| 分类变更数 | 3条 (5.3%) |
| 需复审数 | 2条 (3.5%) |
| 三方一致率 | 96.5% |
| 证据等级A | 46条 (80.7%) |
| 证据等级B | 11条 (19.3%) |
| 置信度HIGH | 54条 (94.7%) |
| 置信度MEDIUM | 3条 (5.3%) |

## 数据完整性验证

- ✅ 所有57条记录均已处理
- ✅ 所有必填字段完整
- ✅ JSONL格式正确
- ✅ 无重复record_id
- ✅ classification_changed字段准确标记
- ✅ gemini_result和codex_result字段均有内容

## 使用说明

### 读取结果文件

```python
import json

with open('verified_results.jsonl', 'r', encoding='utf-8') as f:
    results = [json.loads(line) for line in f]

# 筛选分类变更的记录
changed = [r for r in results if r['classification_changed']]
print(f"发现 {len(changed)} 条任务02错误")

# 筛选需复审的记录
review = [r for r in results if r['needs_human_review']]
print(f"需复审 {len(review)} 条记录")
```

### 提取关键信息

```python
# 统计各分类数量
from collections import Counter
classes = Counter([r['final_class'] for r in results])
print("分类统计:", dict(classes))

# 统计证据等级
evidence = Counter([r['evidence_level'] for r in results])
print("证据等级:", dict(evidence))
```
