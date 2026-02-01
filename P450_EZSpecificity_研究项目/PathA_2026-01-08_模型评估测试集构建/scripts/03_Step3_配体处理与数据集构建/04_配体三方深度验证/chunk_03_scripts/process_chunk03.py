"""
Chunk 03 配体三方深度验证辅助脚本
用于读取CSV并追加JSONL结果
"""
import csv
import json
import os
from typing import Dict, List

INPUT_CSV = "c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks/chunk_03.csv"
OUTPUT_JSONL = "c:/Users/Administrator/Desktop/EZSpecificity_Project/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_03_results/verified_results.jsonl"

def load_csv_records() -> List[Dict]:
    """读取CSV文件返回记录列表"""
    records = []
    with open(INPUT_CSV, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records

def append_jsonl_result(result: Dict):
    """追加一个结果到JSONL文件"""
    with open(OUTPUT_JSONL, 'a', encoding='utf-8') as f:
        f.write(json.dumps(result, ensure_ascii=False) + '\n')

def count_processed_records() -> int:
    """统计已处理的记录数"""
    if not os.path.exists(OUTPUT_JSONL):
        return 0
    with open(OUTPUT_JSONL, 'r', encoding='utf-8') as f:
        return sum(1 for _ in f)

if __name__ == "__main__":
    records = load_csv_records()
    processed = count_processed_records()
    print(f"Total records: {len(records)}")
    print(f"Processed: {processed}")
    print(f"Remaining: {len(records) - processed}")
