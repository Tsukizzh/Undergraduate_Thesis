"""
SOP Prompt 生成器
=================

为每个并行窗口生成完整的上下文Prompt
"""

from pathlib import Path
from datetime import datetime

# SOP模板 - 这是发送给每个新窗口的完整指令
SOP_TEMPLATE = '''# P450配体分类验证任务 - Chunk {chunk_id}

## 你的任务

你负责验证 **Chunk {chunk_id}** 中的 {record_count} 条P450酶-配体记录的分类。

### 工作目录
```
{work_dir}
```

### 输入文件
```
{input_file}
```

### 输出文件
```
{output_file}
```

---

## 验证方法论 (Methodology v1.0)

### 分类标准

| 分类 | 定义 | 关键证据 |
|------|------|----------|
| **SUBSTRATE** | 该配体被酶催化/代谢 | Km、Kcat、turnover、产物鉴定 |
| **INHIBITOR** | 该配体抑制酶活性 | Ki、IC50、抑制曲线 |
| **PRODUCT** | 该配体是酶促反应的产物 | 反应机制明确显示为产物 |
| **EXCLUDE** | 应排除（辅因子、溶剂、野生型vs突变体矛盾等） | 见下方陷阱检测 |

### 陷阱检测 (Critical!)

**陷阱1: 野生型 vs 突变体矛盾**
- 检查PDB标题是否包含mutation/mutant/variant关键词
- 如果PDB是突变体，但UniProt是野生型，需判断：
  - 突变体是否改变了底物特异性？
  - 野生型能否代谢该配体？
- **如果存在矛盾且无法确定野生型活性 → EXCLUDE**

**陷阱2: 底物 vs 产物混淆**
- 检查反应机制：配体是反应物还是生成物？
- 检查PDB是否为"底物复合物"还是"产物复合物"
- **关键词"product"不一定意味着PRODUCT分类！需查文献确认**

**陷阱3: 自杀性抑制剂 / 机制驱动抑制剂**
- 这类配体可能同时有底物特征和抑制特征
- **分类为INHIBITOR**，不是SUBSTRATE

**陷阱4: PDB配体 ≠ 催化底物**
- PDB中的配体只是"结合"在酶上
- 结合不等于催化
- **需要文献证据证明存在催化反应**

### 验证流程

对于每条记录：

1. **读取基本信息**
   - 酶名、UniProt ID
   - 配体名、CCD代码
   - PDB ID、PDB标题
   - 当前分类和置信度

2. **检查PDB信息**
   - 是否突变体？
   - 标题中的关键词

3. **搜索文献** (使用WebSearch)
   - 搜索 "[酶名] [配体名] mechanism"
   - 搜索 "[酶名] [配体名] substrate/inhibitor"
   - 提取定量数据（Ki、IC50、Km）

4. **做出判断**
   - 基于证据判断分类
   - 识别任何陷阱
   - 给出置信度

5. **记录结果**
   - 写入输出文件
   - 包含验证理由和引用

---

## 输出格式

对于每条记录，输出以下JSONL格式（一行一条）：

```json
{{
  "record_id": "原始记录ID",
  "enzyme_name": "酶名",
  "ligand_name": "配体名",
  "pdb_id": "PDB ID",
  "original_class": "原始分类",
  "verified_class": "SUBSTRATE|INHIBITOR|PRODUCT|EXCLUDE",
  "confidence": "HIGH|MEDIUM|LOW",
  "evidence_level": "A|B|C|D",
  "is_mutant": true/false,
  "trap_detected": "none|wt_mutant|substrate_product|suicide_inhibitor|other",
  "quant_data": {{"Ki_uM": null, "IC50_uM": null, "Km_uM": null}},
  "key_evidence": ["证据1", "证据2"],
  "citations": ["PMID:xxx", "DOI:xxx"],
  "reasoning": "简要推理过程",
  "needs_human_review": false,
  "review_reason": ""
}}
```

---

## 执行流程

1. **开始前**：读取输入CSV，确认记录数量
2. **逐条处理**：按顺序处理每条记录
3. **实时保存**：每完成一条立即追加写入输出文件
4. **遇到困难**：如果某条记录无法确定，标记`needs_human_review: true`并继续
5. **完成后**：报告处理统计

---

## 开始工作

请确认你已理解上述方法论，然后：

1. 读取输入文件: `{input_file}`
2. 创建输出文件: `{output_file}`
3. 开始逐条验证

**注意**：
- 使用 WebSearch 工具搜索文献
- 使用 Read 工具读取CSV文件
- 使用 Write 工具写入结果
- 每完成10条记录报告一次进度

请开始处理第1条记录。
'''


def generate_sop_prompts(num_chunks: int = 6):
    """生成所有chunk的SOP prompt"""

    BASE_DIR = Path(__file__).parent.parent.parent.parent
    CHUNKS_DIR = BASE_DIR / "sessions" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证" / "chunks"
    OUTPUT_DIR = Path(__file__).parent

    # 确保目录存在
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    CHUNKS_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("生成SOP Prompts")
    print("=" * 60)

    all_prompts = []

    for i in range(1, num_chunks + 1):
        chunk_id = f"{i:02d}"
        input_file = CHUNKS_DIR / f"chunk_{chunk_id}.csv"
        output_file = CHUNKS_DIR / f"chunk_{chunk_id}_results" / "verified_results.jsonl"

        # 尝试获取记录数
        record_count = "待确认"
        if input_file.exists():
            with open(input_file, 'r', encoding='utf-8') as f:
                record_count = sum(1 for _ in f) - 1  # 减去header

        prompt = SOP_TEMPLATE.format(
            chunk_id=chunk_id,
            record_count=record_count,
            work_dir=str(CHUNKS_DIR),
            input_file=str(input_file),
            output_file=str(output_file)
        )

        # 保存单独的prompt文件
        prompt_file = OUTPUT_DIR / f"sop_prompt_chunk_{chunk_id}.md"
        with open(prompt_file, 'w', encoding='utf-8') as f:
            f.write(prompt)

        all_prompts.append({
            'chunk_id': chunk_id,
            'prompt_file': str(prompt_file),
            'input_file': str(input_file),
            'output_file': str(output_file)
        })

        print(f"  Chunk {chunk_id}: {prompt_file.name}")

    # 保存汇总
    summary_file = OUTPUT_DIR / "all_prompts_summary.json"
    with open(summary_file, 'w', encoding='utf-8') as f:
        import json
        json.dump({
            'created_at': datetime.now().isoformat(),
            'num_chunks': num_chunks,
            'prompts': all_prompts
        }, f, indent=2, ensure_ascii=False)

    print(f"\n汇总文件: {summary_file}")

    print("\n" + "=" * 60)
    print("使用说明")
    print("=" * 60)
    print(f"""
1. 打开 {num_chunks} 个新的Claude Code窗口

2. 在每个窗口中：
   - 复制对应的 sop_prompt_chunk_XX.md 内容
   - 粘贴到新窗口作为第一条消息
   - 等待Claude开始处理

3. 窗口分配：
""")

    for i in range(1, num_chunks + 1):
        print(f"   窗口 {i}: 复制 sop_prompt_chunk_{i:02d}.md")

    print(f"""
4. 所有窗口完成后，运行合并脚本

5. 预计时间：
   - 每条记录约2-5分钟
   - 每个窗口约 {record_count if isinstance(record_count, int) else '?'} 条
   - 并行处理总时间取决于最慢的窗口
""")


if __name__ == "__main__":
    import sys
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    generate_sop_prompts(6)
