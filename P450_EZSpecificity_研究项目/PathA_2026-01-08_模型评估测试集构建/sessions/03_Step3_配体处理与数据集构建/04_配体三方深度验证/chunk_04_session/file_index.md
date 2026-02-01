# Chunk 04 文件索引

## 文件目录结构

```
04_配体三方深度验证/
├── 01_输入数据_chunks/
│   └── chunk_04.csv                    # 输入：57条待验证记录
│
├── 02_输出数据_results/
│   └── chunk_04_results/
│       ├── verified_results.jsonl     # ✅ 主要输出：验证结果
│       ├── session_log.md             # ✅ 会话日志
│       ├── file_index.md              # ✅ 本文件
│       ├── rebuild_jsonl.py           # 辅助脚本：重建JSONL
│       ├── fix_jsonl.py               # 辅助脚本：修复编码（废弃）
│       ├── fix_jsonl_v2.py            # 辅助脚本：修复v2（废弃）
│       └── temp/                      # 临时文件目录
│
└── 03_工具脚本/
    └── chunk_04_scripts/              # （规划但未使用）
```

---

## 核心文件详情

### 1. 输入文件

#### `chunk_04.csv`
- **位置**: `01_输入数据_chunks/chunk_04.csv`
- **来源**: Task 02配体初步分类的chunk分割结果
- **格式**: CSV with header
- **记录数**: 57条
- **列字段**:
  - `record_id`: REC_0172 ~ REC_0228
  - `chunk_id`: 04
  - `uniprot_id`, `enzyme_name`, `ligand_name`, `pdb_id`
  - `original_class`, `verified_class` (Task 02分类)
  - `confidence`, `evidence_level`
  - `is_correct`, `is_mutant`, `trap_detected`, `needs_human_review`, `review_reason`, `reasoning`
- **用途**: 三方深度验证的输入源

---

### 2. 主要输出文件

#### `verified_results.jsonl`
- **位置**: `02_输出数据_results/chunk_04_results/verified_results.jsonl`
- **格式**: JSONL (每行一个JSON对象)
- **记录数**: 57条
- **编码**: UTF-8
- **大小**: ~45KB
- **生成方式**: rebuild_jsonl.py从CSV重建
- **数据结构**:
```json
{
  "record_id": "REC_0172",
  "chunk_id": "04",
  "enzyme_name": "CYP2B4",
  "ligand_name": "phenanthren-9-ylacetaldehyde",
  "pdb_id": "3UAS",
  "original_class": "unknown",
  "verified_class": "PRODUCT",          // Task 02分类
  "final_class": "PRODUCT",             // 三方验证最终分类
  "classification_changed": false,      // 是否纠正Task 02错误
  "confidence": "HIGH",
  "evidence_level": "B",
  "gemini_result": {                    // Gemini文献搜索结果
    "classification": "PRODUCT",
    "pmid": "23276288",
    "summary": "..."
  },
  "codex_result": {                     // Codex结构分析结果
    "classification": "PRODUCT",
    "fe_distance": "4.64A",
    "summary": "..."
  },
  "consensus": true,                    // Gemini和Codex是否一致
  "needs_human_review": false,
  "reasoning": "三方验证确认Task 02分类正确。..."
}
```

- **关键字段说明**:
  - `verified_class`: Task 02的初步分类（盲审基准）
  - `final_class`: 三方验证的最终分类（权威分类）
  - `classification_changed`: true表示纠正了Task 02错误
  - `gemini_result`: 包含PMID、文献摘要、代谢参数等
  - `codex_result`: 包含Fe距离、结合模式、结构特征等
  - `consensus`: 三方意见一致性标志
  - `reasoning`: 中文综合判断理由

- **使用场景**:
  1. EZSpecificity模型评估的测试集构建
  2. Task 02质量评估与错误分析
  3. P450配体分类的黄金标准参考
  4. 合成探针识别的案例库

---

#### `session_log.md`
- **位置**: `02_输出数据_results/chunk_04_results/session_log.md`
- **格式**: Markdown
- **大小**: ~30KB
- **内容**:
  - 执行摘要（记录数、完成情况）
  - 分类统计（分布、错误率）
  - 酶种类分析（CYP2B4/2C5/101A1/1A1）
  - 关键发现（探针识别、Type-II抑制剂、特殊案例）
  - 文献支持质量（PMID列表、覆盖率）
  - Gemini与Codex一致性分析
  - 技术挑战与解决方案
  - 质量保证措施
  - 结论与评估
- **用途**:
  - 验证过程的完整记录
  - 数据质量报告
  - Task 02性能评估依据
  - 其他chunk验证的参考模板

---

#### `file_index.md`
- **位置**: `02_输出数据_results/chunk_04_results/file_index.md`
- **格式**: Markdown
- **内容**: 本文件，提供文件结构和索引

---

### 3. 辅助脚本

#### `rebuild_jsonl.py`
- **位置**: `02_输出数据_results/chunk_04_results/rebuild_jsonl.py`
- **语言**: Python 3
- **功能**: 从CSV重建JSONL，整合三方验证结果
- **输入**: `chunk_04.csv`
- **输出**: `verified_results.jsonl`
- **核心逻辑**:
  1. 读取CSV中的Task 02分类
  2. 应用三方验证决策（`three_way_decisions`字典）
  3. 填充Gemini/Codex结果摘要
  4. 计算consensus标志
  5. 排序并写入JSONL
- **执行**:
```bash
python rebuild_jsonl.py
```
- **输出示例**:
```
=== Chunk 04 Three-Way Verification Results ===
Total records: 57

=== Classification Distribution ===
  EXCLUDE: 29
  INHIBITOR: 14
  PRODUCT: 2
  SUBSTRATE: 12

=== Task 02 Errors Found ===
  classification_changed=true: 1
    REC_0178: INHIBITOR -> SUBSTRATE
```

#### `fix_jsonl.py` *(已废弃)*
- **位置**: `02_输出数据_results/chunk_04_results/fix_jsonl.py`
- **状态**: 废弃（编码修复失败）
- **原因**: 原始JSONL存在严重UTF-8损坏，无法通过编码转换修复
- **替代方案**: 使用rebuild_jsonl.py从CSV重建

#### `fix_jsonl_v2.py` *(已废弃)*
- **位置**: `02_输出数据_results/chunk_04_results/fix_jsonl_v2.py`
- **状态**: 废弃（正则修复失败）
- **原因**: 替换字符导致JSON结构破坏，无法解析

---

## 文件关系图

```
chunk_04.csv (输入)
     |
     | (rebuild_jsonl.py读取)
     v
three_way_decisions (代码中)
     |
     | (整合)
     v
gemini_codex_results (代码中)
     |
     | (写入)
     v
verified_results.jsonl (输出)
     |
     +---> 用于EZSpecificity模型评估
     +---> 用于Task 02质量分析
     +---> 用于P450研究参考
```

---

## 数据完整性检查清单

- [x] 记录数：57条（与输入CSV一致）
- [x] record_id连续性：REC_0172 ~ REC_0228（无遗漏）
- [x] 无重复记录（已修复REC_0180重复问题）
- [x] 所有记录包含gemini_result和codex_result
- [x] 所有记录包含final_class和reasoning
- [x] classification_changed标志正确（1条true，56条false）
- [x] JSON格式有效（无损坏字符）
- [x] UTF-8编码正确（中文显示正常）

---

## 使用指南

### 读取JSONL文件
```python
import json

records = []
with open('verified_results.jsonl', 'r', encoding='utf-8') as f:
    for line in f:
        record = json.loads(line)
        records.append(record)

print(f"Total records: {len(records)}")
```

### 筛选特定分类
```python
substrates = [r for r in records if r['final_class'] == 'SUBSTRATE']
inhibitors = [r for r in records if r['final_class'] == 'INHIBITOR']
print(f"Substrates: {len(substrates)}, Inhibitors: {len(inhibitors)}")
```

### 查找纠正记录
```python
corrections = [r for r in records if r['classification_changed']]
for c in corrections:
    print(f"{c['record_id']}: {c['verified_class']} -> {c['final_class']}")
```

### 统计Gemini-Codex一致性
```python
consensus = sum(1 for r in records if r['consensus'])
print(f"Consensus rate: {consensus}/{len(records)} ({consensus/len(records)*100:.1f}%)")
```

---

## 版本历史

| 版本 | 日期 | 变更说明 |
|------|------|----------|
| v1.0 | 2026-01-28 | 初始版本，完成57条记录的三方验证 |
| v1.1 | 2026-01-28 | 修复JSONL编码问题，从CSV重建 |
| v1.2 | 2026-01-28 | 添加session_log.md和file_index.md |

---

## 相关文档

- [../../../SOP_chunk_04.md](SOP文件位置待确认) - 三方验证SOP
- [session_log.md](./session_log.md) - 详细验证日志
- [../../01_输入数据_chunks/chunk_04.csv](../../01_输入数据_chunks/chunk_04.csv) - 输入数据

---

## 联系与支持

如发现数据问题或需要额外信息，请参考：
- 会话日志：session_log.md
- 原始输入：chunk_04.csv
- 重建脚本：rebuild_jsonl.py

---

**文档生成时间**: 2026-01-28
**数据版本**: v1.2
**验证状态**: ✅ 完成
