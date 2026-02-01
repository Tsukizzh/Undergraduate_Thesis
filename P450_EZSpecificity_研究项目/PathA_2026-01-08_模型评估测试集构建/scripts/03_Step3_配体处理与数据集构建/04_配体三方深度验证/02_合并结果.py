#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
结果合并脚本 - 04_配体三方深度验证
合并12个chunk的验证结果为统一数据集

输入：12个chunk的verified_results.jsonl
输出：合并后的JSONL/CSV文件 + 统计报告
"""

import json
import pandas as pd
from pathlib import Path
from collections import Counter

PROJECT_ROOT = Path(r"c:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建")
INPUT_DIR = PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results"
OUTPUT_DIR = PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/03_合并结果"

def merge_results():
    """合并12个chunk的验证结果"""
    print("=" * 60)
    print("开始合并验证结果")
    print("=" * 60)

    all_records = []
    chunk_stats = []

    for i in range(1, 13):
        chunk_id = f"{i:02d}"
        jsonl_file = INPUT_DIR / f"chunk_{chunk_id}_results/verified_results.jsonl"

        if not jsonl_file.exists():
            print(f"[WARN] Chunk {chunk_id}: 文件不存在 - {jsonl_file}")
            continue

        chunk_records = []
        with open(jsonl_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                try:
                    record = json.loads(line)
                    chunk_records.append(record)
                except json.JSONDecodeError as e:
                    print(f"[WARN] Chunk {chunk_id} 第{line_num}行: JSON解析错误 - {e}")

        # 统计该chunk
        needs_review = sum(1 for r in chunk_records if r.get('needs_human_review', False))
        consensus_count = sum(1 for r in chunk_records if r.get('consensus', False))
        changed_count = sum(1 for r in chunk_records if r.get('classification_changed', False))

        chunk_stats.append({
            "chunk_id": chunk_id,
            "total_records": len(chunk_records),
            "needs_review": needs_review,
            "consensus_count": consensus_count,
            "classification_changed": changed_count
        })

        all_records.extend(chunk_records)
        print(f"[OK] Chunk {chunk_id}: {len(chunk_records)}条, {changed_count}条发现错误, {needs_review}条需复审")

    if not all_records:
        print("\n[ERROR] 没有找到任何记录！")
        return

    # 总体统计
    print("\n" + "-" * 40)
    print("总体统计")
    print("-" * 40)

    total = len(all_records)
    class_dist = Counter(r.get('final_class', 'UNKNOWN') for r in all_records)
    confidence_dist = Counter(r.get('confidence', 'UNKNOWN') for r in all_records)
    evidence_dist = Counter(r.get('evidence_level', 'UNKNOWN') for r in all_records)
    needs_review = sum(1 for r in all_records if r.get('needs_human_review', False))
    consensus = sum(1 for r in all_records if r.get('consensus', False))

    print(f"总记录数: {total}")
    print(f"需人工审核: {needs_review} ({needs_review/total*100:.1f}%)")
    print(f"三方一致: {consensus} ({consensus/total*100:.1f}%)")
    print(f"\n分类分布:")
    for cls, cnt in sorted(class_dist.items()):
        print(f"  {cls}: {cnt} ({cnt/total*100:.1f}%)")
    print(f"\n置信度分布:")
    for conf, cnt in sorted(confidence_dist.items()):
        print(f"  {conf}: {cnt} ({cnt/total*100:.1f}%)")
    print(f"\n证据等级分布:")
    for ev, cnt in sorted(evidence_dist.items()):
        print(f"  {ev}: {cnt} ({cnt/total*100:.1f}%)")

    # 保存合并结果
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # JSONL
    jsonl_output = OUTPUT_DIR / "三方深度验证结果_merged.jsonl"
    with open(jsonl_output, 'w', encoding='utf-8') as f:
        for record in all_records:
            f.write(json.dumps(record, ensure_ascii=False) + '\n')
    print(f"\n[OK] JSONL已保存: {jsonl_output.name}")

    # CSV (展平嵌套字段)
    flat_records = []
    for r in all_records:
        flat_r = {k: v for k, v in r.items() if not isinstance(v, dict)}
        # 展平gemini_result
        if 'gemini_result' in r and isinstance(r['gemini_result'], dict):
            for k, v in r['gemini_result'].items():
                flat_r[f'gemini_{k}'] = v
        # 展平codex_result
        if 'codex_result' in r and isinstance(r['codex_result'], dict):
            for k, v in r['codex_result'].items():
                flat_r[f'codex_{k}'] = v
        flat_records.append(flat_r)

    df = pd.DataFrame(flat_records)
    csv_output = OUTPUT_DIR / "三方深度验证结果_merged.csv"
    df.to_csv(csv_output, index=False, encoding='utf-8-sig')
    print(f"[OK] CSV已保存: {csv_output.name}")

    # 统计JSON
    stats = {
        "total_records": total,
        "needs_human_review": needs_review,
        "consensus_count": consensus,
        "class_distribution": dict(class_dist),
        "confidence_distribution": dict(confidence_dist),
        "evidence_distribution": dict(evidence_dist),
        "chunk_stats": chunk_stats
    }
    stats_output = OUTPUT_DIR / "三方深度验证统计.json"
    with open(stats_output, 'w', encoding='utf-8') as f:
        json.dump(stats, f, indent=2, ensure_ascii=False)
    print(f"[OK] 统计已保存: {stats_output.name}")

    print("\n" + "=" * 60)
    print("合并完成！")
    print("=" * 60)

if __name__ == "__main__":
    merge_results()
