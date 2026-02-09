"""
Step 1 v2: Download PDB files for 627 unique PDBs
- Reuse existing 158 PDBs from archived v1 data
- Download missing PDBs from RCSB
"""

import pandas as pd
import os
import shutil
import requests
from pathlib import Path

# Paths
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # PathA_2026-01-08_模型评估测试集构建
SOURCE_CSV = BASE_DIR / "source_data" / "01_核心数据" / "独立测试集_740条.csv"
OLD_PDB_DIR = BASE_DIR / "_archived_v1_按UniProt去重_Step1至3" / "data_v1" / "01_Step1_PDB原始文件" / "PDB结构文件_158个"
NEW_PDB_DIR = BASE_DIR / "data" / "01_Step1_PDB文件"
LOG_DIR = BASE_DIR / "sessions" / "2026-01-21_Step1_v2_下载PDB文件"

def main():
    # Read the 740 records CSV
    df = pd.read_csv(SOURCE_CSV)
    required_pdb_ids = set(df['pdb_id'].str.upper().unique())
    print(f"Required unique PDB IDs: {len(required_pdb_ids)}")

    # Get existing PDB files
    existing_files = list(OLD_PDB_DIR.glob("*.*"))
    existing_pdb_ids = set()
    for f in existing_files:
        pdb_id = f.stem.upper()
        existing_pdb_ids.add(pdb_id)
    print(f"Existing PDB files in archive: {len(existing_pdb_ids)}")

    # Find overlap and missing
    overlap = required_pdb_ids & existing_pdb_ids
    missing = required_pdb_ids - existing_pdb_ids
    print(f"Overlap (can reuse): {len(overlap)}")
    print(f"Missing (need download): {len(missing)}")

    # Copy existing files
    copied_count = 0
    for f in existing_files:
        pdb_id = f.stem.upper()
        if pdb_id in required_pdb_ids:
            dest = NEW_PDB_DIR / f.name
            shutil.copy2(f, dest)
            copied_count += 1
    print(f"Copied {copied_count} existing PDB files")

    # Download missing files
    downloaded = []
    failed = []

    for i, pdb_id in enumerate(sorted(missing)):
        pdb_id_lower = pdb_id.lower()

        # Try PDB format first
        url_pdb = f"https://files.rcsb.org/download/{pdb_id_lower}.pdb"
        url_cif = f"https://files.rcsb.org/download/{pdb_id_lower}.cif"

        success = False
        for url, ext in [(url_pdb, ".pdb"), (url_cif, ".cif")]:
            try:
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    dest = NEW_PDB_DIR / f"{pdb_id}{ext}"
                    with open(dest, 'wb') as f:
                        f.write(response.content)
                    downloaded.append(pdb_id)
                    success = True
                    break
            except Exception as e:
                continue

        if not success:
            failed.append(pdb_id)

        if (i + 1) % 50 == 0:
            print(f"Progress: {i + 1}/{len(missing)} (downloaded: {len(downloaded)}, failed: {len(failed)})")

    print(f"\nDownload complete!")
    print(f"Downloaded: {len(downloaded)}")
    print(f"Failed: {len(failed)}")

    # Save log
    log_content = f"""# Step 1 v2 执行日志

## 执行时间
2026-01-21

## 统计
- 需要的唯一PDB数: {len(required_pdb_ids)}
- 已有的PDB数（从v1归档复制）: {len(overlap)}
- 需要下载的PDB数: {len(missing)}
- 成功下载: {len(downloaded)}
- 下载失败: {len(failed)}
- 最终PDB文件数: {copied_count + len(downloaded)}

## 下载失败的PDB
{chr(10).join(failed) if failed else "无"}
"""

    with open(LOG_DIR / "session_log.md", 'w', encoding='utf-8') as f:
        f.write(log_content)

    # Save PDB list
    all_pdbs = sorted(required_pdb_ids)
    with open(LOG_DIR / "required_pdb_list.txt", 'w') as f:
        f.write('\n'.join(all_pdbs))

    if failed:
        with open(LOG_DIR / "failed_pdbs.txt", 'w') as f:
            f.write('\n'.join(failed))

if __name__ == "__main__":
    main()
