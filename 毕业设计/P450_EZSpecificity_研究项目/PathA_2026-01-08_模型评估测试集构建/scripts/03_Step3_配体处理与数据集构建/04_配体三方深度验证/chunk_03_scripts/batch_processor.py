#!/usr/bin/env python3
"""
Batch processor for Chunk 03 ligand verification
This script processes remaining records with three-way collaboration
"""

import csv
import json
import sys
from pathlib import Path

# Input/Output paths
INPUT_CSV = Path("c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_03.csv")
OUTPUT_JSONL = Path("c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_03_results/verified_results.jsonl")

def count_processed():
    """Count how many records have been processed"""
    if not OUTPUT_JSONL.exists():
        return 0
    with open(OUTPUT_JSONL, 'r', encoding='utf-8') as f:
        return sum(1 for line in f if line.strip())

def get_remaining_records():
    """Get list of remaining unprocessed records"""
    processed_ids = set()
    if OUTPUT_JSONL.exists():
        with open(OUTPUT_JSONL, 'r', encoding='utf-8') as f:
            for line in f:
                if line.strip():
                    data = json.loads(line)
                    processed_ids.add(data['record_id'])

    remaining = []
    with open(INPUT_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['record_id'] not in processed_ids:
                remaining.append(row)

    return remaining

if __name__ == "__main__":
    processed_count = count_processed()
    remaining_records = get_remaining_records()

    print(f"Processed: {processed_count}/57")
    print(f"Remaining: {len(remaining_records)}/57")
    print(f"\nNext records to process:")
    for i, rec in enumerate(remaining_records[:10], 1):
        print(f"  {i}. {rec['record_id']}: {rec['enzyme_name']} + {rec['ligand_name']}")

    if len(remaining_records) > 10:
        print(f"  ... and {len(remaining_records) - 10} more")
