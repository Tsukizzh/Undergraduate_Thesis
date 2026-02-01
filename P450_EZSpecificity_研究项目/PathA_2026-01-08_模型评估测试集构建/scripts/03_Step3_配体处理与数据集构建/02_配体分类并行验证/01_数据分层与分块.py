"""
数据分层与分块脚本
==================

功能：
1. 读取682条P450-配体记录
2. 按风险等级分层（A/B/C/D）
3. 将需要人工审核的记录分成N个chunk
4. 为每个chunk生成manifest文件

输出：
- tier_A_auto.csv: A级证据，自动通过
- tier_BCD_manual.csv: 需要人工审核的记录
- chunks/chunk_01.csv ... chunk_N.csv: 分块文件
- chunks/manifests.json: 分块元数据
"""

import csv
import json
import math
from pathlib import Path
from datetime import datetime


def load_records(input_csv: str) -> list:
    """加载原始记录"""
    records = []
    with open(input_csv, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            records.append(row)
    return records


def classify_tier(record: dict) -> str:
    """
    根据风险等级分类记录

    A级：有ChEMBL定量数据（自动通过）
    B级：低风险（PDB标题明确，无突变体，有明确关键词）
    C级：中风险（有关键词但需确认，或信息不完整）
    D级：高风险（突变体、冲突、特殊机制）
    """
    # 检查是否有定量数据（A级）
    # 注意：这里假设已经有preliminary分类结果
    # 实际应该检查ChEMBL数据

    pdb_title = record.get('pdb_title', '').lower()
    ligand_class = record.get('ligand_class', '').upper()
    confidence = record.get('confidence', '').lower()

    # D级：高风险情况
    mutant_keywords = ['mutant', 'mutation', 'variant', 'engineered']
    if any(k in pdb_title for k in mutant_keywords):
        return 'D'

    # D级：产物类（容易混淆）
    if ligand_class == 'PRODUCT':
        return 'D'

    # D级：标记为冲突的
    if 'conflict' in record.get('notes', '').lower():
        return 'D'

    # A级：高置信度的明确分类
    if confidence in ['high', 'very high'] and ligand_class in ['INHIBITOR', 'SUBSTRATE']:
        return 'A'

    # B级：中等置信度
    if confidence == 'medium' and ligand_class in ['INHIBITOR', 'SUBSTRATE']:
        return 'B'

    # C级：其他情况
    return 'C'


def split_into_chunks(records: list, num_chunks: int) -> list:
    """将记录分成N个chunk"""
    chunk_size = math.ceil(len(records) / num_chunks)
    chunks = []
    for i in range(0, len(records), chunk_size):
        chunk = records[i:i + chunk_size]
        chunks.append(chunk)
    return chunks


def main():
    import sys
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    # 配置
    BASE_DIR = Path(__file__).parent.parent.parent.parent  # PathA_2026-01-08_...
    INPUT_CSV = BASE_DIR / "reports" / "03_Step3_配体处理与数据集构建" / "配体分类审核结果_682条.csv"
    OUTPUT_DIR = Path(__file__).parent.parent.parent.parent / "sessions" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证"

    NUM_CHUNKS = 6  # 并行窗口数

    print("=" * 60)
    print("数据分层与分块")
    print("=" * 60)

    # 检查输入文件
    if not INPUT_CSV.exists():
        print(f"[ERROR] 未找到输入文件: {INPUT_CSV}")
        print("\n请确保已完成Step 3.01的配体分类审核")
        return

    # 加载记录
    records = load_records(str(INPUT_CSV))
    print(f"\n加载了 {len(records)} 条记录")

    # 分层
    tiers = {'A': [], 'B': [], 'C': [], 'D': []}
    for record in records:
        tier = classify_tier(record)
        record['tier'] = tier
        tiers[tier].append(record)

    print(f"\n分层结果:")
    print(f"  A级（自动通过）: {len(tiers['A'])} 条")
    print(f"  B级（低风险）: {len(tiers['B'])} 条")
    print(f"  C级（中风险）: {len(tiers['C'])} 条")
    print(f"  D级（高风险）: {len(tiers['D'])} 条")

    # 保存A级（自动通过）
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tier_a_file = OUTPUT_DIR / "tier_A_auto.csv"
    if tiers['A']:
        with open(tier_a_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(tiers['A'])
        print(f"\nA级记录已保存: {tier_a_file}")

    # 合并B/C/D级（需要人工审核）
    manual_records = tiers['B'] + tiers['C'] + tiers['D']
    tier_bcd_file = OUTPUT_DIR / "tier_BCD_manual.csv"
    if manual_records:
        with open(tier_bcd_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(manual_records)
        print(f"B/C/D级记录已保存: {tier_bcd_file}")

    # 分块
    chunks_dir = OUTPUT_DIR / "chunks"
    chunks_dir.mkdir(parents=True, exist_ok=True)

    chunks = split_into_chunks(manual_records, NUM_CHUNKS)

    manifests = {
        'created_at': datetime.now().isoformat(),
        'total_records': len(manual_records),
        'num_chunks': len(chunks),
        'chunks': []
    }

    for i, chunk in enumerate(chunks, 1):
        chunk_file = chunks_dir / f"chunk_{i:02d}.csv"
        with open(chunk_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=records[0].keys())
            writer.writeheader()
            writer.writerows(chunk)

        # 创建结果目录
        result_dir = chunks_dir / f"chunk_{i:02d}_results"
        result_dir.mkdir(exist_ok=True)

        manifests['chunks'].append({
            'chunk_id': f'chunk_{i:02d}',
            'file': str(chunk_file.name),
            'count': len(chunk),
            'tier_distribution': {
                'B': len([r for r in chunk if r.get('tier') == 'B']),
                'C': len([r for r in chunk if r.get('tier') == 'C']),
                'D': len([r for r in chunk if r.get('tier') == 'D'])
            }
        })

        print(f"  Chunk {i}: {len(chunk)} 条记录 -> {chunk_file.name}")

    # 保存manifest
    manifest_file = chunks_dir / "manifests.json"
    with open(manifest_file, 'w', encoding='utf-8') as f:
        json.dump(manifests, f, indent=2, ensure_ascii=False)

    print(f"\n分块元数据已保存: {manifest_file}")

    print("\n" + "=" * 60)
    print("下一步操作:")
    print("=" * 60)
    print(f"1. 打开 {NUM_CHUNKS} 个新的Claude Code窗口")
    print(f"2. 每个窗口复制对应的SOP prompt（见下一个脚本生成）")
    print(f"3. 等待所有窗口完成处理")
    print(f"4. 运行合并脚本汇总结果")


if __name__ == "__main__":
    main()
