"""
Step 5: ESM酶嵌入生成脚本

功能：
    为Enzymes.csv中的酶序列生成ESM-2蛋白质嵌入，存储到LMDB数据库

输入：
    - data/04_Step4_格式修正后数据/Enzymes.csv
    - 必需列：Protein sequence, Enzyme_Index

输出：
    - data/05_Step5_ESM酶嵌入/enzyme_features.lmdb

技术规格：
    - 模型：esm2_t33_650M_UR50D (650M参数, 1280维输出)
    - 提取层：第33层 (最后一层)
    - 存储格式：per-residue embedding (L, 1280)，非Mean Pooled
    - LMDB key：str(Enzyme_Index).encode() (来自CSV的Enzyme_Index列)
    - LMDB value：pickle dict {'embedding': (L,1280), 'active_site': int, 'sequence': (L,)}

三方审核：Claude Code + Codex + Gemini (2026-01-30)
"""

import os
import sys
import pickle
import numpy as np
import pandas as pd
import torch
import lmdb
import esm
from tqdm import tqdm
from pathlib import Path

# ============================================================
# 配置
# ============================================================
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent.parent  # PathA_2026-01-08_模型评估测试集构建

INPUT_CSV = PROJECT_DIR / "data" / "04_Step4_格式修正后数据" / "Enzymes.csv"
OUTPUT_DIR = PROJECT_DIR / "data" / "05_Step5_ESM酶嵌入"
OUTPUT_LMDB = OUTPUT_DIR / "enzyme_features.lmdb"

# ESM配置
MAX_SEQUENCE_LENGTH = 1000  # 超过此长度的序列将被跳过
DEVICE = "cuda:0"  # 硬编码为cuda:0，与原始代码一致

# LMDB配置
LMDB_MAP_SIZE = 4 * 1024 * 1024 * 1024  # 4GB，足够292条酶

# 氨基酸字母表 - 与原始代码 src/Datasets/const.py 完全一致
# 注意：允许 X, Z, U, B 等非标准氨基酸
LETTER_TO_NUM = {
    'C': 4, 'D': 3, 'S': 15, 'Q': 5, 'K': 11, 'I': 9,
    'P': 14, 'T': 16, 'F': 13, 'A': 0, 'G': 7, 'H': 8,
    'E': 6, 'L': 10, 'R': 1, 'W': 17, 'V': 19,
    'N': 2, 'Y': 18, 'M': 12, 'X': 20, 'Z': 21, 'U': 22, 'B': 23
}

# ============================================================
# 辅助函数
# ============================================================

def convert_protein_sequence_to_number(sequence: str) -> np.ndarray:
    """
    将氨基酸序列转换为数字编码
    与原始代码 src/Datasets/const.py:letter_to_num 完全一致

    Args:
        sequence: 氨基酸序列字符串

    Returns:
        numpy数组，每个氨基酸对应一个整数

    Raises:
        ValueError: 如果序列包含未知字符
    """
    result = []
    for aa in sequence.upper():
        if aa not in LETTER_TO_NUM:
            raise ValueError(f"Unknown amino acid: {aa}")
        result.append(LETTER_TO_NUM[aa])

    return np.array(result, dtype=np.int64)


def validate_sequence(sequence: str, max_length: int = MAX_SEQUENCE_LENGTH) -> tuple[bool, str]:
    """
    验证蛋白质序列是否符合要求
    与原始代码保持一致，允许 X, Z, U, B 等字符

    Args:
        sequence: 氨基酸序列
        max_length: 最大允许长度

    Returns:
        (is_valid, reason)
    """
    if not isinstance(sequence, str):
        return False, "Sequence is not a string"

    sequence = sequence.strip()

    if len(sequence) == 0:
        return False, "Empty sequence"

    if len(sequence) > max_length:
        return False, f"Sequence too long ({len(sequence)} > {max_length})"

    # 检查是否所有字符都在允许的字母表中
    invalid_chars = set(sequence.upper()) - set(LETTER_TO_NUM.keys())
    if invalid_chars:
        return False, f"Unknown amino acids: {invalid_chars}"

    return True, "OK"


def generate_esm_embeddings(
    df: pd.DataFrame,
    output_lmdb_path: Path,
    device: str = DEVICE,
    resume: bool = True
) -> dict:
    """
    为DataFrame中的酶序列生成ESM嵌入

    Args:
        df: 包含 'Protein sequence' 列的DataFrame
        output_lmdb_path: 输出LMDB路径
        device: 计算设备
        resume: 是否跳过已处理的序列

    Returns:
        统计信息字典
    """
    print("=" * 60)
    print("Step 5: ESM酶嵌入生成")
    print("=" * 60)

    # 统计
    stats = {
        "total": len(df),
        "processed": 0,
        "skipped_length": 0,
        "skipped_unknown_aa": 0,  # 未知氨基酸（不在LETTER_TO_NUM中）
        "skipped_exists": 0,
        "failed": 0
    }

    # 加载ESM模型
    print(f"\n[1/4] 加载ESM-2模型 (esm2_t33_650M_UR50D)...")
    print(f"      设备: {device}")

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(torch.device(device))
    model.eval()  # 重要：设置为评估模式
    batch_converter = alphabet.get_batch_converter()

    print(f"      模型加载完成")

    # 创建输出目录
    output_lmdb_path.parent.mkdir(parents=True, exist_ok=True)

    # 打开LMDB（使用subdir=False，与原始代码一致）
    print(f"\n[2/4] 初始化LMDB数据库...")
    print(f"      路径: {output_lmdb_path}")

    env = lmdb.open(
        str(output_lmdb_path),
        map_size=LMDB_MAP_SIZE,
        subdir=False,  # 关键：与原始代码一致
        create=True,
        readonly=False
    )

    # 获取已存在的keys（用于resume）
    existing_keys = set()
    if resume:
        with env.begin(write=False) as txn:
            cursor = txn.cursor()
            for key, _ in cursor:
                existing_keys.add(key.decode())
        print(f"      已存在 {len(existing_keys)} 条记录")

    # 处理每条酶序列
    print(f"\n[3/4] 生成ESM嵌入...")
    print(f"      总数: {len(df)}")

    # 检查是否有Enzyme_Index列
    use_enzyme_index_col = "Enzyme_Index" in df.columns
    if use_enzyme_index_col:
        print(f"      使用 'Enzyme_Index' 列作为LMDB key")
    else:
        print(f"      使用DataFrame行索引作为LMDB key")

    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing"):
        # 使用Enzyme_Index列（如果存在），否则使用行索引
        if use_enzyme_index_col:
            key = str(int(row["Enzyme_Index"]))
        else:
            key = str(idx)

        # 检查是否已存在
        if key in existing_keys:
            stats["skipped_exists"] += 1
            continue

        sequence = row["Protein sequence"]

        # 验证序列
        is_valid, reason = validate_sequence(sequence)
        if not is_valid:
            if "too long" in reason:
                stats["skipped_length"] += 1
            else:
                stats["skipped_unknown_aa"] += 1
            print(f"      跳过 key={key}: {reason}")
            continue

        sequence = sequence.strip().upper()

        try:
            # 准备batch数据
            batch_data = [(key, sequence)]
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            batch_tokens = batch_tokens.to(torch.device(device))

            # 生成嵌入
            with torch.no_grad():
                results = model(
                    batch_tokens,
                    repr_layers=[33],
                    return_contacts=False  # 节省显存
                )

            # 提取per-residue嵌入
            # token_representations shape: [batch, seq_len+2, 1280]
            # +2是因为有BOS和EOS token，需要去掉
            token_representations = results["representations"][33]

            # 提取实际序列部分（去掉BOS和EOS）
            # [1:1+len(seq)] 提取从位置1开始的len(seq)个token
            embedding = token_representations[0, 1:1+len(sequence), :].cpu().numpy()

            # 转换序列为数字编码
            sequence_encoded = convert_protein_sequence_to_number(sequence)

            # 构建LMDB value
            # 与原始代码 src/Datasets/create_features.py:335-340 保持一致
            value = {
                "embedding": embedding,  # shape: (L, 1280)
                "active_site": 1,        # 默认值，P450测试集不使用此字段
                "sequence": sequence_encoded  # shape: (L,)
            }

            # 写入LMDB
            with env.begin(write=True) as txn:
                txn.put(key.encode(), pickle.dumps(value))

            stats["processed"] += 1

        except Exception as e:
            stats["failed"] += 1
            print(f"      失败 key={key}: {e}")

    env.close()

    # 输出统计
    print(f"\n[4/4] 完成")
    print(f"      处理成功: {stats['processed']}")
    print(f"      跳过(已存在): {stats['skipped_exists']}")
    print(f"      跳过(过长): {stats['skipped_length']}")
    print(f"      跳过(未知AA): {stats['skipped_unknown_aa']}")
    print(f"      失败: {stats['failed']}")

    return stats


def verify_lmdb(lmdb_path: Path, expected_count: int) -> bool:
    """
    验证LMDB输出

    Args:
        lmdb_path: LMDB文件路径
        expected_count: 预期的记录数

    Returns:
        验证是否通过
    """
    print("\n" + "=" * 60)
    print("验证LMDB输出")
    print("=" * 60)

    if not lmdb_path.exists():
        print(f"错误: LMDB文件不存在: {lmdb_path}")
        return False

    env = lmdb.open(str(lmdb_path), subdir=False, readonly=True)

    with env.begin(write=False) as txn:
        stat = txn.stat()
        actual_count = stat["entries"]
        print(f"记录数: {actual_count} (预期: {expected_count})")

        # 随机检查几条记录
        print("\n抽样检查:")
        cursor = txn.cursor()
        checked = 0
        errors = []

        for key, value in cursor:
            if checked >= 5:
                break

            key_str = key.decode()
            data = pickle.loads(value)

            # 验证结构
            if "embedding" not in data:
                errors.append(f"Key {key_str}: missing 'embedding'")
            elif len(data["embedding"].shape) != 2 or data["embedding"].shape[1] != 1280:
                errors.append(f"Key {key_str}: wrong embedding shape {data['embedding'].shape}")
            elif np.isnan(data["embedding"]).any():
                errors.append(f"Key {key_str}: embedding contains NaN")
            else:
                print(f"  Key={key_str}: embedding shape={data['embedding'].shape}, "
                      f"sequence len={len(data.get('sequence', []))}")

            checked += 1

    env.close()

    if errors:
        print("\n发现错误:")
        for e in errors:
            print(f"  - {e}")
        return False

    print("\n验证通过 [OK]")
    return True


# ============================================================
# 主函数
# ============================================================

def main():
    print("\n" + "=" * 60)
    print("Step 5: ESM酶嵌入生成")
    print("=" * 60)
    print(f"输入: {INPUT_CSV}")
    print(f"输出: {OUTPUT_LMDB}")

    # 检查输入文件
    if not INPUT_CSV.exists():
        print(f"错误: 输入文件不存在: {INPUT_CSV}")
        sys.exit(1)

    # 检查GPU
    if not torch.cuda.is_available():
        print("错误: CUDA不可用，此脚本需要GPU")
        sys.exit(1)

    print(f"\nGPU: {torch.cuda.get_device_name(0)}")
    print(f"CUDA版本: {torch.version.cuda}")

    # 读取输入数据
    print(f"\n读取Enzymes.csv...")
    df = pd.read_csv(INPUT_CSV)
    print(f"  总行数: {len(df)}")
    print(f"  列: {list(df.columns)}")

    # 验证必需列
    if "Protein sequence" not in df.columns:
        print(f"错误: 缺少 'Protein sequence' 列")
        sys.exit(1)

    # 生成嵌入
    stats = generate_esm_embeddings(df, OUTPUT_LMDB)

    # 验证输出
    expected_count = stats["processed"] + stats["skipped_exists"]
    verify_lmdb(OUTPUT_LMDB, expected_count)

    print("\n" + "=" * 60)
    print("Step 5 完成")
    print("=" * 60)


if __name__ == "__main__":
    main()
