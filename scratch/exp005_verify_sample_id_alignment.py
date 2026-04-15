"""
Step 1: 验证 pt_cache 里的 sample_id 就是 split CSV 的原始行号

目的：
  之后要建 dock sidecar，必须知道能不能通过 sample_id 直接索引 split CSV。
  如果 sample_id 是 CSV 行号，那 sidecar 就是 dock_indices[k] = csv.iloc[index.sample_ids[k]]["Dock Index"]
  如果不是，需要其他 join 策略。

方法：
  对 3 个 split (train/val/test) 各抽 20 个样本做穿透验证：
    base_index = pt_cache_allfix_unified/random/{split}/index.pt
    split_csv  = data/splits/random/{training/val/testing}_datas_0_pt.csv
    对每个抽样 k:
      sample_id = base_index["sample_ids"][k]
      csv_row   = csv.iloc[sample_id]
      验证:
        csv_row["Enzyme Index"] == base_index["enzyme_ids"][k]
        csv_row["Substrate Index"] == base_index["substrate_ids"][k]

  额外检查：
    max(sample_ids) < len(csv)  (超出就说明 sample_id 不是 CSV 行号)
    csv row 行数 和 base_index 长度的关系（base_index 行数 <= csv 行数，差值 = orphan filter 过滤数）

通过此脚本全部在服务器 CPU 运行，本地只读脚本内容传过去执行。
"""
from __future__ import annotations

RECON_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv
import random
from pathlib import Path
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
CACHE_DIR = BASE / "data/pt_cache_allfix_unified/random"
SPLIT_DIR = BASE / "data/splits/random"

# split_name -> csv_filename mapping
SPLIT_CSV = {
    "train": "training_datas_0_pt.csv",
    "val":   "val_datas_0_pt.csv",
    "test":  "testing_datas_0_pt.csv",
}

N_SPOT = 20  # spot-check samples per split
random.seed(42)


def load_csv_rows(path):
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def verify_split(split_name, csv_name):
    print(f"\n=== split: {split_name} ===")
    index_path = CACHE_DIR / split_name / "index.pt"
    csv_path = SPLIT_DIR / csv_name

    # 加载
    idx = torch.load(index_path, map_location="cpu", weights_only=False)
    rows = load_csv_rows(csv_path)
    n_idx = len(idx["sample_ids"])
    n_csv = len(rows)

    print(f"  index.pt rows: {n_idx}")
    print(f"  csv rows:      {n_csv}")
    print(f"  diff (orphan): {n_csv - n_idx}")

    sids = idx["sample_ids"].tolist()
    eids = idx["enzyme_ids"].tolist()
    subids = idx["substrate_ids"].tolist()

    print(f"  sample_id range: min={min(sids)}, max={max(sids)}")
    print(f"  sample_id monotonic increasing? {sids == sorted(sids)}")

    # 基础合理性：max(sample_id) < len(csv)
    if max(sids) >= n_csv:
        print(f"  [ERROR] max(sample_id)={max(sids)} >= csv rows {n_csv}")
        return False

    # 20 个 spot check
    pick = random.sample(range(n_idx), min(N_SPOT, n_idx))
    print(f"  spot-checking {len(pick)} random positions:")

    n_match = 0
    failures = []
    for k in pick:
        sid = sids[k]
        row = rows[sid]
        csv_enz = int(row["Enzyme Index"])
        csv_sub = int(row["Substrate Index"])
        idx_enz = eids[k]
        idx_sub = subids[k]
        ok_enz = csv_enz == idx_enz
        ok_sub = csv_sub == idx_sub
        if ok_enz and ok_sub:
            n_match += 1
        else:
            failures.append({
                "k": k, "sample_id": sid,
                "csv": (csv_enz, csv_sub),
                "idx": (idx_enz, idx_sub),
            })

    print(f"  matched: {n_match}/{len(pick)}")
    if failures:
        print(f"  failures:")
        for f in failures[:5]:
            print(f"    {f}")
        return False

    return True


def main():
    all_ok = True
    for split_name, csv_name in SPLIT_CSV.items():
        ok = verify_split(split_name, csv_name)
        all_ok = all_ok and ok

    print()
    print("=" * 60)
    if all_ok:
        print("OVERALL: PASS  — sample_id IS the CSV row index")
        print("        dock sidecar builder can use:")
        print("        dock_indices[k] = csv.iloc[index.sample_ids[k]]['Dock Index']")
    else:
        print("OVERALL: FAIL — do NOT proceed with the simple sidecar plan")
    print("=" * 60)
    import sys
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
'''


def main():
    """local driver: ssh script to server via stdin, utf-8 encoded"""
    import subprocess, sys
    print("[local driver] sending recon script to server via ssh stdin...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=RECON_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=180,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
