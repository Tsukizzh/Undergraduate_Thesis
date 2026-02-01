#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据分割脚本 - 04_配体三方深度验证
将682条记录分割为12个chunk，每个约57条记录

输入：02_配体分类并行验证的merged结果
输出：12个chunk CSV文件
"""

import pandas as pd
import json
from pathlib import Path

# 路径配置
PROJECT_ROOT = Path(r"c:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建")
INPUT_FILE = PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/02_配体分类并行验证/03_合并结果/配体分类验证结果_并行处理_merged.csv"
OUTPUT_DIR = PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/01_输入数据_chunks"

def split_data():
    """分割数据为12个chunk"""
    print("=" * 60)
    print("开始数据分割")
    print("=" * 60)

    # 读取数据
    print(f"\n1. 读取输入文件：{INPUT_FILE}")
    df = pd.read_csv(INPUT_FILE)
    total_records = len(df)
    print(f"   总记录数：{total_records}")

    # 计算分割
    num_chunks = 12
    chunk_size = total_records // num_chunks
    remainder = total_records % num_chunks

    print(f"\n2. 分割配置：")
    print(f"   窗口数：{num_chunks}")
    print(f"   基础大小：{chunk_size}条/窗口")
    print(f"   剩余记录：{remainder}条")
    print(f"   分配策略：前{remainder}个窗口各{chunk_size+1}条，其余各{chunk_size}条")

    # 分割数据
    manifest = {
        "total_records": total_records,
        "num_chunks": num_chunks,
        "chunks": []
    }

    start_idx = 0
    for i in range(num_chunks):
        # 前remainder个chunk多分配1条
        current_chunk_size = chunk_size + (1 if i < remainder else 0)
        end_idx = start_idx + current_chunk_size

        # 提取chunk数据
        chunk_df = df.iloc[start_idx:end_idx].copy()

        # 更新chunk_id列
        chunk_id = f"{i+1:02d}"
        if 'chunk_id' in chunk_df.columns:
            chunk_df['chunk_id'] = chunk_id
        else:
            chunk_df.insert(0, 'chunk_id', chunk_id)

        # 保存CSV
        output_file = OUTPUT_DIR / f"chunk_{chunk_id}.csv"
        chunk_df.to_csv(output_file, index=False, encoding='utf-8')

        # 更新manifest
        manifest["chunks"].append({
            "chunk_id": chunk_id,
            "file": f"chunk_{chunk_id}.csv",
            "record_range": f"REC_{df.iloc[start_idx]['record_id']}_to_REC_{df.iloc[end_idx-1]['record_id']}",
            "start_index": int(start_idx),
            "end_index": int(end_idx),
            "num_records": int(current_chunk_size)
        })

        print(f"   [OK] Chunk {chunk_id}: 记录 {start_idx+1}-{end_idx} ({current_chunk_size}条) -> {output_file.name}")

        start_idx = end_idx

    # 保存manifest
    manifest_file = OUTPUT_DIR / "manifest.json"
    with open(manifest_file, 'w', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)
    print(f"\n3. Manifest文件已保存：{manifest_file.name}")

    print("\n" + "=" * 60)
    print("数据分割完成！")
    print("=" * 60)
    print(f"\n输出目录：{OUTPUT_DIR}")
    print(f"生成文件：")
    print(f"  - 12个chunk CSV文件 (chunk_01.csv ~ chunk_12.csv)")
    print(f"  - 1个manifest.json文件")

    return manifest

if __name__ == "__main__":
    manifest = split_data()
