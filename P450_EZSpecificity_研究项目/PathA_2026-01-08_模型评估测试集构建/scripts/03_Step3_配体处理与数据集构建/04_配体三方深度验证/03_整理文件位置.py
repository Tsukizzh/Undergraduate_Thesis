#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
整理任务04文件位置
将错误位置的文件移动到正确位置
"""

import os
import shutil
from pathlib import Path

# 项目根目录
PROJECT_ROOT = Path(r"c:\Users\Administrator\Desktop\EZSpecificity_Project\P450_EZSpecificity_研究项目\PathA_2026-01-08_模型评估测试集构建")

# 需要移动的文件列表 (源路径, 目标路径)
files_to_move = [
    # Chunk 04
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/session_log.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_04_session/session_log.md"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/file_index.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_04_session/file_index.md"
    ),
    # Chunk 06
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_06_results/session_log.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_06_session/session_log.md"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_06_results/file_index.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_06_session/file_index.md"
    ),
    # Chunk 08
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_08_results/session_log.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_08_session/session_log.md"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_08_results/file_index.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_08_session/file_index.md"
    ),
    # Chunk 11
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_11_results/session_log.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_11_session/session_log.md"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_11_results/file_index.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_11_session/file_index.md"
    ),
    # Chunk 05
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/03_操作日志/chunk_05_session_log.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_05_session/session_log.md"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/03_操作日志/chunk_05_file_index.md",
        PROJECT_ROOT / "sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_05_session/file_index.md"
    ),
]

# 需要移动的脚本文件 (chunk_04的修复脚本)
scripts_to_move = [
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/fix_jsonl.py",
        PROJECT_ROOT / "scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_04_scripts/fix_jsonl.py"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/fix_jsonl_v2.py",
        PROJECT_ROOT / "scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_04_scripts/fix_jsonl_v2.py"
    ),
    (
        PROJECT_ROOT / "data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_04_results/rebuild_jsonl.py",
        PROJECT_ROOT / "scripts/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_04_scripts/rebuild_jsonl.py"
    ),
]

def move_file(src, dst):
    """移动文件，如果源文件存在"""
    if src.exists():
        # 确保目标目录存在
        dst.parent.mkdir(parents=True, exist_ok=True)

        # 如果目标文件已存在，先删除
        if dst.exists():
            print(f"  [WARN] 目标文件已存在，覆盖: {dst.name}")
            dst.unlink()

        # 移动文件
        shutil.move(str(src), str(dst))
        print(f"  [OK] 移动: {src.name} -> {dst.parent.name}/{dst.name}")
        return True
    else:
        print(f"  [WARN] 源文件不存在: {src}")
        return False

def main():
    print("=" * 80)
    print("开始整理任务04文件位置")
    print("=" * 80)

    # 移动session日志和文件索引
    print("\n[步骤1] 移动session日志和文件索引到sessions目录")
    success_count = 0
    for src, dst in files_to_move:
        if move_file(src, dst):
            success_count += 1

    print(f"\n[步骤1完成] 成功移动 {success_count}/{len(files_to_move)} 个文件")

    # 移动脚本文件
    print("\n[步骤2] 移动脚本文件到scripts目录")
    script_success_count = 0
    for src, dst in scripts_to_move:
        if move_file(src, dst):
            script_success_count += 1

    print(f"\n[步骤2完成] 成功移动 {script_success_count}/{len(scripts_to_move)} 个脚本")

    # 检查是否所有chunk都有session文件
    print("\n[步骤3] 检查所有chunk的session文件完整性")
    missing_chunks = []
    for i in range(1, 13):
        chunk_id = f"{i:02d}"
        session_log = PROJECT_ROOT / f"sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_{chunk_id}_session/session_log.md"
        file_index = PROJECT_ROOT / f"sessions/03_Step3_配体处理与数据集构建/04_配体三方深度验证/chunk_{chunk_id}_session/file_index.md"

        if not session_log.exists() or not file_index.exists():
            missing_chunks.append(chunk_id)
            print(f"  [MISS] Chunk {chunk_id}: 缺失session文件")
        else:
            print(f"  [OK] Chunk {chunk_id}: session文件完整")

    if missing_chunks:
        print(f"\n[WARN] 以下chunk缺失session文件: {', '.join(missing_chunks)}")
    else:
        print(f"\n[OK] 所有12个chunk的session文件都已就位")

    # 检查所有chunk的JSONL文件
    print("\n[步骤4] 检查所有chunk的JSONL文件")
    missing_jsonl = []
    for i in range(1, 13):
        chunk_id = f"{i:02d}"
        jsonl_file = PROJECT_ROOT / f"data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/02_输出数据_results/chunk_{chunk_id}_results/verified_results.jsonl"

        if not jsonl_file.exists():
            missing_jsonl.append(chunk_id)
            print(f"  [MISS] Chunk {chunk_id}: 缺失verified_results.jsonl")
        else:
            print(f"  [OK] Chunk {chunk_id}: verified_results.jsonl 存在")

    if missing_jsonl:
        print(f"\n[WARN] 以下chunk缺失JSONL文件: {', '.join(missing_jsonl)}")
    else:
        print(f"\n[OK] 所有12个chunk的JSONL文件都已生成")

    print("\n" + "=" * 80)
    print("文件整理完成！")
    print("=" * 80)

    # 总结
    print("\n【整理总结】")
    print(f"- 移动session文件: {success_count}/{len(files_to_move)}")
    print(f"- 移动脚本文件: {script_success_count}/{len(scripts_to_move)}")
    print(f"- 缺失session文件的chunk: {len(missing_chunks)}/12")
    print(f"- 缺失JSONL文件的chunk: {len(missing_jsonl)}/12")

    if missing_chunks or missing_jsonl:
        print("\n[WARN] 存在缺失文件，请检查对应chunk的处理情况")
    else:
        print("\n[OK] 所有文件已就位，可以运行合并脚本了")

if __name__ == "__main__":
    main()
