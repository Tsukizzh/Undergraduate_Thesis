"""
一键准备并行验证 (v2.0 - 修复Schema匹配问题)
================

修复内容：
1. 添加全局唯一 record_id
2. 统一输出Schema（与合并脚本匹配）
3. 修复路径推导
4. 标准化 trap_detected 值
5. 明确验证 final_class
"""

import csv
import json
import math
from pathlib import Path
from datetime import datetime


def main():
    import sys
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    # 配置
    NUM_CHUNKS = 6

    # 路径（使用相对路径，更健壮）
    SCRIPT_DIR = Path(__file__).parent
    BASE_DIR = SCRIPT_DIR.parent.parent.parent
    INPUT_CSV = BASE_DIR / "data" / "03_Step3_配体处理与数据集构建" / "01_配体分类审核" / "03_最终输出" / "主审核文件_682条P450配体分类.csv"

    # 数据文件放在 data 目录下
    DATA_DIR = BASE_DIR / "data" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证"
    CHUNKS_DIR = DATA_DIR / "01_输入数据_chunks"
    RESULTS_DIR = DATA_DIR / "02_输出数据_results"

    # Sessions 目录
    SESSIONS_DIR = BASE_DIR / "sessions" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证"

    print("=" * 70)
    print("一键准备并行验证 v2.0")
    print("=" * 70)

    # 检查输入文件
    if not INPUT_CSV.exists():
        print(f"[ERROR] 未找到输入文件: {INPUT_CSV}")
        return

    # 读取数据并添加全局唯一 record_id
    records = []
    with open(INPUT_CSV, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        original_fieldnames = list(reader.fieldnames)
        for idx, row in enumerate(reader, 1):
            # 添加全局唯一ID: REC_{全局索引:04d}
            row['record_id'] = f"REC_{idx:04d}"
            records.append(row)

    print(f"\n读取了 {len(records)} 条记录，添加了全局唯一 record_id")

    # 新的字段列表（record_id 放在第一列）
    fieldnames = ['record_id'] + original_fieldnames

    # 创建目录
    CHUNKS_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # 分块
    chunk_size = math.ceil(len(records) / NUM_CHUNKS)
    chunks = []
    for i in range(0, len(records), chunk_size):
        chunks.append(records[i:i + chunk_size])

    print(f"\n分成 {len(chunks)} 个chunk:")

    # 保存每个chunk并生成SOP
    for i, chunk in enumerate(chunks, 1):
        chunk_id = f"{i:02d}"

        # 记录每个chunk包含的record_id范围
        record_ids = [r['record_id'] for r in chunk]
        start_id = record_ids[0]
        end_id = record_ids[-1]

        # 1. 保存chunk CSV (输入数据，包含record_id)
        chunk_file = CHUNKS_DIR / f"chunk_{chunk_id}.csv"
        with open(chunk_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(chunk)

        # 2. 创建结果目录 (输出数据)
        result_dir = RESULTS_DIR / f"chunk_{chunk_id}_results"
        result_dir.mkdir(exist_ok=True)

        # 3. 生成SOP Prompt
        sop_content = generate_sop(
            chunk_id=chunk_id,
            record_count=len(chunk),
            record_id_range=f"{start_id} ~ {end_id}",
            input_file=str(chunk_file),
            output_dir=str(result_dir),
            session_log_file=str(SESSIONS_DIR / f"chunk_{chunk_id}_session_log.md")
        )

        sop_file = SCRIPT_DIR / f"sop_prompt_chunk_{chunk_id}.md"
        with open(sop_file, 'w', encoding='utf-8') as f:
            f.write(sop_content)

        print(f"  Chunk {chunk_id}: {len(chunk)} 条记录 ({start_id} ~ {end_id})")
        print(f"    - 数据文件: {chunk_file.name}")
        print(f"    - SOP文件: {sop_file.name}")

    # 保存元数据
    manifest = {
        'created_at': datetime.now().isoformat(),
        'total_records': len(records),
        'num_chunks': len(chunks),
        'chunks': [
            {
                'chunk_id': f'{i+1:02d}',
                'count': len(chunks[i]),
                'record_id_range': f"{chunks[i][0]['record_id']} ~ {chunks[i][-1]['record_id']}",
                'input_file': f'chunk_{i+1:02d}.csv',
                'output_dir': f'chunk_{i+1:02d}_results'
            }
            for i in range(len(chunks))
        ]
    }

    manifest_file = CHUNKS_DIR / "manifest.json"
    with open(manifest_file, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    print(f"\n元数据: {manifest_file}")

    # 打印操作指南
    print("\n" + "=" * 70)
    print("准备完成！接下来的操作步骤：")
    print("=" * 70)
    print(f"""
【步骤1】打开6个新的Claude Code窗口
  - 在VS Code中按 Ctrl+Shift+P
  - 输入 "Claude" 并选择 "Claude Code: New Session"
  - 重复6次

【步骤2】在每个窗口粘贴对应的SOP Prompt
  - 窗口1: 打开 {SCRIPT_DIR / 'sop_prompt_chunk_01.md'}
           全选复制(Ctrl+A, Ctrl+C) → 粘贴到窗口1 → 回车
  - 窗口2: 打开 sop_prompt_chunk_02.md → 粘贴到窗口2
  - 窗口3: 打开 sop_prompt_chunk_03.md → 粘贴到窗口3
  - 窗口4: 打开 sop_prompt_chunk_04.md → 粘贴到窗口4
  - 窗口5: 打开 sop_prompt_chunk_05.md → 粘贴到窗口5
  - 窗口6: 打开 sop_prompt_chunk_06.md → 粘贴到窗口6

【步骤3】等待所有窗口完成
  - 每个窗口会自动处理约 {chunk_size} 条记录
  - 每个窗口会把结果保存到对应的 chunk_XX_results 文件夹
  - 预计每个窗口需要 6-8 小时

【步骤4】所有窗口完成后
  - 回到原来的窗口（这个窗口）
  - 告诉Claude "都完成了"
  - Claude会运行合并脚本，生成最终报告
""")


def generate_sop(chunk_id: str, record_count: int, record_id_range: str,
                 input_file: str, output_dir: str, session_log_file: str) -> str:
    """生成SOP Prompt内容（v2.0 - Schema与合并脚本匹配）"""

    return f'''# P450配体分类验证任务 - Chunk {chunk_id}

## 你的任务

你负责验证 **Chunk {chunk_id}** 中的 **{record_count} 条** P450酶-配体记录的分类。

**Record ID 范围**: {record_id_range}

---

## 文件路径

**输入文件（读取）：**
```
{input_file}
```

**输出文件夹（写入数据）：**
```
{output_dir}
```

**输出文件：**
- `{output_dir}\\verified_results.jsonl` - 验证结果（每条一行JSON）
- `{output_dir}\\verification_report.md` - 验证报告

**操作日志（写入文档）：**
```
{session_log_file}
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
{{
  "record_id": "REC_0001",
  "chunk_id": "{chunk_id}",
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
}}
```

**字段说明：**
- `record_id`: 从CSV读取的全局唯一ID（REC_XXXX格式）
- `chunk_id`: 当前chunk编号 ("{chunk_id}")
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
{input_file}
```

确认CSV包含以下列：
- `record_id` - 全局唯一ID
- `final_class` - 要验证的分类
- `confidence` - 之前的置信度评估

然后开始逐条验证，从第1条开始。
'''


if __name__ == "__main__":
    main()
