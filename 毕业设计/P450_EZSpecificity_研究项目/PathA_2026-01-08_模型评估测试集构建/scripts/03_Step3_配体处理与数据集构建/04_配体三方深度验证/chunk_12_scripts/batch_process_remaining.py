#!/usr/bin/env python3
"""
批处理脚本：处理Chunk 12剩余记录
使用三方协作（Gemini+Codex+Claude）进行深度验证
"""

import json
import csv
from pathlib import Path

# 路径配置
INPUT_CSV = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_12.csv"
OUTPUT_JSONL = r"c:/Users/Administrator/Desktop/EZSpecificity_Project/毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_12_results/verified_results.jsonl"

# 读取已处理记录
processed_ids = set()
if Path(OUTPUT_JSONL).exists():
    with open(OUTPUT_JSONL, 'r', encoding='utf-8') as f:
        for line in f:
            record = json.loads(line)
            processed_ids.add(record['record_id'])

print(f"已处理记录数：{len(processed_ids)}")

# 读取所有记录
all_records = []
with open(INPUT_CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        all_records.append(row)

print(f"总记录数：{len(all_records)}")

# 待处理记录
remaining_records = [r for r in all_records if r['record_id'] not in processed_ids]
print(f"剩余待处理记录数：{len(remaining_records)}")

# 输出待处理列表
print("\n待处理记录列表：")
for i, r in enumerate(remaining_records, 1):
    print(f"{i}. {r['record_id']}: {r['enzyme_name']} + {r['ligand_name']} (PDB {r['pdb_id']})")
