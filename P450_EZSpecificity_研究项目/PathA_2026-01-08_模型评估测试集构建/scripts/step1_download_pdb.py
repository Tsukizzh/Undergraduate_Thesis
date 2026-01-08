"""
Step 1: 下载161个PDB文件
- 从ML训练数据集_186个.csv读取PDB ID
- 排除25个与ESIBank重叠的UniProt ID
- 下载剩余PDB文件到data/00_raw/pdb_files/
"""

import pandas as pd
import requests
import os
import json
import time
from pathlib import Path

# 需要排除的25个UniProt ID
EXCLUDE_UNIPROT = {
    'A2TEF2', 'C4B644', 'E3VWI3', 'O24782', 'P00189', 'P00191',
    'P05093', 'P08684', 'P11511', 'P19099', 'P22680', 'P23295',
    'P33261', 'Q06069', 'Q09128', 'Q16850', 'Q6TBX7', 'Q6VVX0',
    'Q6WG30', 'Q83WG3', 'Q8VQF6', 'Q9UNU6', 'Q9Y6A2', 'Q9ZAU3', 'S4UX02'
}

def parse_uniprot_ids(uniprot_str):
    """解析uniprot_ids列（JSON格式）"""
    try:
        ids = json.loads(uniprot_str.replace("'", '"'))
        return set(ids)
    except:
        return set()

def download_pdb(pdb_id, output_dir):
    """下载单个PDB文件"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.exists(output_path):
        return True, "already exists"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(response.text)
            return True, "downloaded"
        else:
            return False, f"HTTP {response.status_code}"
    except Exception as e:
        return False, str(e)

def main():
    # 路径设置
    base_dir = Path(__file__).parent.parent
    csv_path = base_dir / "source_data" / "01_核心数据" / "ML训练数据集_186个.csv"
    output_dir = base_dir / "data" / "00_raw" / "pdb_files"
    log_dir = base_dir / "logs"

    # 确保输出目录存在
    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # 读取CSV
    print(f"读取CSV: {csv_path}")
    df = pd.read_csv(csv_path)
    print(f"总共 {len(df)} 条记录")

    # 过滤掉需要排除的UniProt ID
    def should_exclude(row):
        uniprot_ids = parse_uniprot_ids(row['uniprot_ids'])
        return bool(uniprot_ids & EXCLUDE_UNIPROT)

    df['exclude'] = df.apply(should_exclude, axis=1)
    excluded_df = df[df['exclude']]
    filtered_df = df[~df['exclude']]

    print(f"排除 {len(excluded_df)} 条记录（与ESIBank重叠）")
    print(f"剩余 {len(filtered_df)} 条记录需要下载")

    # 保存排除列表
    excluded_path = log_dir / "excluded_pdbs.csv"
    excluded_df[['pdb_id', 'uniprot_ids']].to_csv(excluded_path, index=False)
    print(f"排除列表已保存: {excluded_path}")

    # 下载PDB文件
    pdb_ids = filtered_df['pdb_id'].tolist()
    success_count = 0
    fail_count = 0
    failed_pdbs = []

    print(f"\n开始下载 {len(pdb_ids)} 个PDB文件...")
    for i, pdb_id in enumerate(pdb_ids, 1):
        success, msg = download_pdb(pdb_id, output_dir)
        if success:
            success_count += 1
            status = "OK"
        else:
            fail_count += 1
            failed_pdbs.append((pdb_id, msg))
            status = "FAIL"

        print(f"[{i}/{len(pdb_ids)}] {pdb_id}: {status} ({msg})")

        # 避免请求过快
        if msg == "downloaded":
            time.sleep(0.3)

    # 保存失败列表
    if failed_pdbs:
        fail_path = log_dir / "errors" / "failed_pdbs.txt"
        fail_path.parent.mkdir(parents=True, exist_ok=True)
        with open(fail_path, 'w') as f:
            for pdb_id, msg in failed_pdbs:
                f.write(f"{pdb_id}\t{msg}\n")
        print(f"\n失败列表已保存: {fail_path}")

    # 总结
    print(f"\n===== 下载完成 =====")
    print(f"成功: {success_count}")
    print(f"失败: {fail_count}")
    print(f"输出目录: {output_dir}")

if __name__ == "__main__":
    main()
