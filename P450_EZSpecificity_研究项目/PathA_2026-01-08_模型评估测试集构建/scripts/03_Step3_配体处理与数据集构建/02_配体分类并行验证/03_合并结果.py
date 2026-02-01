"""
结果合并脚本 (v2.0 - 添加容错机制)
============

修复内容：
1. 处理JSON解析错误（空行、坏行）
2. 处理缺失字段
3. 验证完整性（682条全处理）
4. 记录问题行
"""

import json
import csv
from pathlib import Path
from collections import Counter
from datetime import datetime


def merge_results(num_chunks: int = 6, expected_total: int = 682):
    """合并所有chunk的结果"""

    BASE_DIR = Path(__file__).parent.parent.parent.parent
    DATA_DIR = BASE_DIR / "data" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证"
    RESULTS_DIR = DATA_DIR / "02_输出数据_results"
    OUTPUT_DIR = DATA_DIR / "03_合并结果"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("合并验证结果 v2.0")
    print("=" * 60)

    all_results = []
    chunk_stats = []
    bad_lines = []  # 记录解析失败的行
    missing_chunks = []  # 记录缺失的chunk

    # 读取所有chunk的结果
    for i in range(1, num_chunks + 1):
        chunk_id = f"{i:02d}"
        result_file = RESULTS_DIR / f"chunk_{chunk_id}_results" / "verified_results.jsonl"

        if not result_file.exists():
            print(f"\n[WARNING] Chunk {chunk_id} 结果文件不存在: {result_file}")
            missing_chunks.append(chunk_id)
            continue

        chunk_results = []
        line_num = 0
        with open(result_file, 'r', encoding='utf-8') as f:
            for line in f:
                line_num += 1
                line = line.strip()

                # 跳过空行
                if not line:
                    continue

                # 尝试解析JSON
                try:
                    result = json.loads(line)
                except json.JSONDecodeError as e:
                    bad_lines.append({
                        'chunk_id': chunk_id,
                        'line_num': line_num,
                        'error': str(e),
                        'content': line[:100] + '...' if len(line) > 100 else line
                    })
                    continue

                # 确保chunk_id存在（如果SOP没写入的话）
                result['chunk_id'] = result.get('chunk_id', chunk_id)

                # 验证必需字段（使用默认值容错）
                result.setdefault('record_id', f"UNKNOWN_{chunk_id}_{line_num}")
                result.setdefault('verified_class', 'UNKNOWN')
                result.setdefault('confidence', 'UNKNOWN')
                result.setdefault('evidence_level', 'D')
                result.setdefault('trap_detected', 'none')
                result.setdefault('needs_human_review', False)
                result.setdefault('review_reason', '')
                result.setdefault('reasoning', '')

                chunk_results.append(result)
                all_results.append(result)

        print(f"\nChunk {chunk_id}: {len(chunk_results)} 条记录")

        # 统计该chunk
        stats = {
            'chunk_id': chunk_id,
            'total': len(chunk_results),
            'classifications': Counter(r['verified_class'] for r in chunk_results),
            'confidence': Counter(r['confidence'] for r in chunk_results),
            'evidence_level': Counter(r.get('evidence_level', 'D') for r in chunk_results),
            'traps_detected': Counter(r.get('trap_detected', 'none') for r in chunk_results),
            'needs_review': sum(1 for r in chunk_results if r.get('needs_human_review', False))
        }
        chunk_stats.append(stats)

        print(f"  分类: {dict(stats['classifications'])}")
        print(f"  需人工审核: {stats['needs_review']} 条")

    # 报告解析错误
    if bad_lines:
        print(f"\n[WARNING] 解析失败的行: {len(bad_lines)} 条")
        for bl in bad_lines[:5]:  # 最多显示5条
            print(f"  Chunk {bl['chunk_id']}, Line {bl['line_num']}: {bl['error']}")
        if len(bad_lines) > 5:
            print(f"  ... 还有 {len(bad_lines) - 5} 条错误")

    # 整体统计
    print("\n" + "=" * 60)
    print("整体统计")
    print("=" * 60)

    total = len(all_results)
    print(f"\n总记录数: {total}")

    # 检查完整性
    if total == 0:
        print("[ERROR] 没有找到任何结果！")
        return

    if total < expected_total:
        print(f"[WARNING] 预期 {expected_total} 条，实际 {total} 条，缺少 {expected_total - total} 条")
    elif total > expected_total:
        print(f"[WARNING] 预期 {expected_total} 条，实际 {total} 条，多出 {total - expected_total} 条（可能有重复）")
    else:
        print(f"✓ 记录数完整：{total}/{expected_total}")

    if missing_chunks:
        print(f"[WARNING] 缺失的Chunk: {', '.join(missing_chunks)}")

    # 分类统计
    class_counts = Counter(r['verified_class'] for r in all_results)
    print(f"\n分类分布:")
    for cls, count in class_counts.most_common():
        pct = count / total * 100
        print(f"  {cls}: {count} ({pct:.1f}%)")

    # 置信度统计
    conf_counts = Counter(r['confidence'] for r in all_results)
    print(f"\n置信度分布:")
    for conf, count in conf_counts.most_common():
        pct = count / total * 100
        print(f"  {conf}: {count} ({pct:.1f}%)")

    # 证据等级统计
    ev_counts = Counter(r.get('evidence_level', 'D') for r in all_results)
    print(f"\n证据等级分布:")
    for ev, count in ev_counts.most_common():
        pct = count / total * 100
        print(f"  {ev}: {count} ({pct:.1f}%)")

    # 陷阱检测统计
    trap_counts = Counter(r.get('trap_detected', 'none') for r in all_results)
    print(f"\n陷阱检测:")
    for trap, count in trap_counts.items():
        if trap != 'none' and count > 0:
            print(f"  {trap}: {count}")

    # 需要人工审核的记录
    needs_review = [r for r in all_results if r.get('needs_human_review', False)]
    print(f"\n需人工审核: {len(needs_review)} 条 ({len(needs_review)/total*100:.1f}%)")

    # 保存合并结果
    merged_file = OUTPUT_DIR / "配体分类验证结果_并行处理_merged.jsonl"
    with open(merged_file, 'w', encoding='utf-8') as f:
        for result in all_results:
            f.write(json.dumps(result, ensure_ascii=False) + '\n')
    print(f"\n合并结果已保存: {merged_file}")

    # 转换为CSV格式
    csv_file = OUTPUT_DIR / "配体分类验证结果_并行处理_merged.csv"
    if all_results:
        flat_results = []
        for r in all_results:
            flat = {
                'record_id': r.get('record_id', ''),
                'chunk_id': r.get('chunk_id', ''),
                'uniprot_id': r.get('uniprot_id', ''),
                'enzyme_name': r.get('enzyme_name', ''),
                'ligand_name': r.get('ligand_name', ''),
                'pdb_id': r.get('pdb_id', ''),
                'original_class': r.get('original_class', ''),
                'verified_class': r.get('verified_class', ''),
                'confidence': r.get('confidence', ''),
                'evidence_level': r.get('evidence_level', ''),
                'is_correct': r.get('is_correct', ''),
                'is_mutant': r.get('is_mutant', False),
                'trap_detected': r.get('trap_detected', 'none'),
                'needs_human_review': r.get('needs_human_review', False),
                'review_reason': r.get('review_reason', ''),
                'reasoning': r.get('reasoning', '')[:200]
            }
            flat_results.append(flat)

        with open(csv_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=flat_results[0].keys())
            writer.writeheader()
            writer.writerows(flat_results)
        print(f"CSV格式已保存: {csv_file}")

    # 保存统计报告
    stats_file = OUTPUT_DIR / "配体分类验证统计_并行处理.json"
    with open(stats_file, 'w', encoding='utf-8') as f:
        json.dump({
            'merged_at': datetime.now().isoformat(),
            'expected_total': expected_total,
            'actual_total': total,
            'num_chunks': num_chunks,
            'missing_chunks': missing_chunks,
            'bad_lines_count': len(bad_lines),
            'classifications': dict(class_counts),
            'confidence': dict(conf_counts),
            'evidence_level': dict(ev_counts),
            'traps_detected': dict(trap_counts),
            'needs_review': len(needs_review),
            'chunk_stats': chunk_stats
        }, f, indent=2, ensure_ascii=False)
    print(f"统计报告已保存: {stats_file}")

    # 保存错误日志（如果有的话）
    if bad_lines:
        error_file = OUTPUT_DIR / "配体分类验证_解析错误.json"
        with open(error_file, 'w', encoding='utf-8') as f:
            json.dump(bad_lines, f, indent=2, ensure_ascii=False)
        print(f"错误日志已保存: {error_file}")

    # 检查记录ID唯一性
    record_ids = [r['record_id'] for r in all_results]
    unique_ids = set(record_ids)
    print(f"\n记录ID唯一性检查:")
    if len(unique_ids) == total:
        print(f"  ✓ 所有 {total} 条记录ID唯一")
    else:
        duplicates = total - len(unique_ids)
        print(f"  ✗ 发现重复！唯一ID数: {len(unique_ids)}, 总记录数: {total}")
        print(f"  重复记录数: {duplicates}")

        # 找出重复的ID
        id_counts = Counter(record_ids)
        dup_ids = [id for id, count in id_counts.items() if count > 1]
        print(f"  重复的ID: {dup_ids[:10]}{'...' if len(dup_ids) > 10 else ''}")

    # 最终状态
    print("\n" + "=" * 60)
    if total == expected_total and not missing_chunks and not bad_lines:
        print("✓ 合并完成！数据完整无错误")
    else:
        print("⚠ 合并完成，但存在以下问题：")
        if total != expected_total:
            print(f"  - 记录数不匹配: {total}/{expected_total}")
        if missing_chunks:
            print(f"  - 缺失Chunk: {missing_chunks}")
        if bad_lines:
            print(f"  - 解析错误: {len(bad_lines)} 行")
    print("=" * 60)


if __name__ == "__main__":
    import sys
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    merge_results(6, 682)
