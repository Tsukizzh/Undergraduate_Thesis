"""
Chunk 04 - 配体三方深度验证自动化脚本

该脚本用于批量处理chunk_04的57条记录，执行三方协作验证流程。
注意：此脚本仅用于辅助处理，实际的Gemini/Codex调用仍需手动执行。
"""

import csv
import json
from pathlib import Path
from typing import Dict, List
from datetime import datetime

# 路径配置
BASE_DIR = Path("c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建")
INPUT_CSV = BASE_DIR / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_04.csv"
OUTPUT_JSONL = BASE_DIR / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/verified_results.jsonl"
TEMP_DIR = BASE_DIR / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/temp"

def load_records() -> List[Dict]:
    """加载输入CSV记录"""
    records = []
    with open(INPUT_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records

def save_result(result: Dict):
    """保存单条验证结果到JSONL"""
    with open(OUTPUT_JSONL, 'a', encoding='utf-8') as f:
        f.write(json.dumps(result, ensure_ascii=False) + '\n')

def create_result_template(record: Dict) -> Dict:
    """创建结果模板"""
    return {
        "record_id": record["record_id"],
        "chunk_id": record["chunk_id"],
        "enzyme_name": record["enzyme_name"],
        "ligand_name": record["ligand_name"],
        "pdb_id": record["pdb_id"],
        "original_class": record["original_class"],
        "verified_class": record["verified_class"],  # 任务02的分类
        "final_class": "",  # 待填写：本次任务04的分类
        "classification_changed": False,  # 待填写
        "confidence": "",  # HIGH/MEDIUM/LOW
        "evidence_level": "",  # A/B/C/D
        "gemini_result": {
            "classification": "",
            "data": "",
            "pmid": "",
            "summary": ""
        },
        "codex_result": {
            "classification": "",
            "fe_distance": "",
            "binding_mode": "",
            "summary": ""
        },
        "consensus": False,  # 待填写
        "needs_human_review": False,  # 待填写
        "reasoning": ""  # 待填写：详细推理过程
    }

def generate_gemini_prompt(enzyme: str, ligand: str, pdb: str) -> str:
    """生成Gemini搜索提示词"""
    return f"""Literature search for: {enzyme} + {ligand}

Search for:
1. Is it SUBSTRATE (metabolized), INHIBITOR (inhibits enzyme), or PRODUCT?
2. Quantitative data: Ki, IC50, Km values (with units)
3. PubMed ID (PMID) if available

Return format:
- Classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
- Data: Ki=X nM, IC50=Y uM, etc.
- PMID: XXXXXXXX
- If no literature found: explicitly state "No literature found"

Brief answer only, no long explanations."""

def generate_codex_prompt(enzyme: str, ligand: str, pdb: str) -> str:
    """生成Codex分析提示词"""
    return f"""Analyze PDB structure: {pdb}
Enzyme: {enzyme}
Ligand: {ligand}

Analyze:
1. Fe-ligand distance (if P450 heme is present)
2. Binding mode (active site vs allosteric)
3. Does structure support SUBSTRATE/INHIBITOR/PRODUCT classification?

Return:
- Fe-N distance: X.XX Å (if nitrogen coordinates to Fe)
- Binding mode: [description]
- Structural classification: SUBSTRATE/INHIBITOR/PRODUCT/EXCLUDE
- Reasoning: [brief]

Brief answer only."""

def print_processing_summary():
    """打印处理进度汇总"""
    if not OUTPUT_JSONL.exists():
        print("尚无处理结果")
        return

    results = []
    with open(OUTPUT_JSONL, 'r', encoding='utf-8') as f:
        for line in f:
            results.append(json.loads(line))

    total = len(results)
    changed = sum(1 for r in results if r.get("classification_changed", False))
    needs_review = sum(1 for r in results if r.get("needs_human_review", False))

    class_counts = {}
    for r in results:
        fc = r.get("final_class", "UNKNOWN")
        class_counts[fc] = class_counts.get(fc, 0) + 1

    print(f"\n=== 处理进度汇总 ===")
    print(f"已处理记录数：{total}/57")
    print(f"分类变更数：{changed}")
    print(f"需人工复审：{needs_review}")
    print(f"\n分类分布：")
    for cls, count in sorted(class_counts.items()):
        print(f"  {cls}: {count}")

if __name__ == "__main__":
    # 确保输出目录存在
    OUTPUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    # 加载记录
    records = load_records()
    print(f"加载了 {len(records)} 条记录")

    # 打印前3条记录示例
    print("\n前3条记录示例：")
    for i, rec in enumerate(records[:3], 1):
        print(f"\n记录 {i}: {rec['record_id']}")
        print(f"  酶: {rec['enzyme_name']}")
        print(f"  配体: {rec['ligand_name']}")
        print(f"  PDB: {rec['pdb_id']}")
        print(f"  任务02分类: {rec['verified_class']}")

    # 打印当前处理进度
    print_processing_summary()

    print("\n脚本加载完成。可以使用以下函数：")
    print("- load_records() - 加载所有记录")
    print("- create_result_template(record) - 创建结果模板")
    print("- generate_gemini_prompt(enzyme, ligand, pdb) - 生成Gemini提示词")
    print("- generate_codex_prompt(enzyme, ligand, pdb) - 生成Codex提示词")
    print("- save_result(result) - 保存结果到JSONL")
    print("- print_processing_summary() - 打印进度汇总")
