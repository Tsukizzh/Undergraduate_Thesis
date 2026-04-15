"""
Step 2: generate dock_sidecar.pt for 3 splits

Design (codex-approved Option B):
  store sidecars directly in the final overlay dir:
    pt_cache_dualgraph_allfix_unified/random/{train,val,test}/dock_sidecar.pt

Each sidecar file is a dict:
  {
      'sample_ids': int32 [N],   # copy from base index.pt (for integrity check)
      'dock_indices': int32 [N], # csv.iloc[sample_ids[k]]['Dock Index']
      'source_csv':  str,        # which CSV was used
      'source_index': str,       # which base index.pt was used
  }

Verification built-in: before saving, re-check sample_ids match base, and
spot-check 20 rows that (enzyme, substrate) from CSV matches base index.
If any mismatch, abort without writing.

Runs entirely in-server via ssh stdin. No local data pull.
"""
from __future__ import annotations

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, random
from pathlib import Path
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
BASE_CACHE = BASE / "data/pt_cache_allfix_unified/random"
DST_CACHE  = BASE / "data/pt_cache_dualgraph_allfix_unified/random"
SPLIT_CSV_DIR = BASE / "data/splits/random"
POCKET_DIR = BASE / "data/structure/str_tmp_data/pocket"

SPLIT_CSV = {
    "train": "training_datas_0_pt.csv",
    "val":   "val_datas_0_pt.csv",
    "test":  "testing_datas_0_pt.csv",
}

N_SPOT = 20
random.seed(123)

ok_overall = True
for split, csv_name in SPLIT_CSV.items():
    print(f"\n=== {split} ===")
    base_index_path = BASE_CACHE / split / "index.pt"
    csv_path = SPLIT_CSV_DIR / csv_name
    out_path = DST_CACHE / split / "dock_sidecar.pt"

    base_idx = torch.load(base_index_path, map_location="cpu", weights_only=False)
    with open(csv_path, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))

    sids = base_idx["sample_ids"]
    N = len(sids)
    print(f"  base N: {N}")
    print(f"  csv N:  {len(rows)}")

    # Build dock_indices aligned with base row order
    dock_indices = torch.zeros(N, dtype=torch.int32)
    enz_from_csv = torch.zeros(N, dtype=torch.int32)
    sub_from_csv = torch.zeros(N, dtype=torch.int32)
    for k in range(N):
        sid = int(sids[k])
        row = rows[sid]
        dock_indices[k] = int(row["Dock Index"])
        enz_from_csv[k] = int(row["Enzyme Index"])
        sub_from_csv[k] = int(row["Substrate Index"])

    # Full integrity check: every (enz, sub) from CSV must match base index
    base_eids = base_idx["enzyme_ids"].to(torch.int32)
    base_subids = base_idx["substrate_ids"].to(torch.int32)
    enz_mismatch = int((enz_from_csv != base_eids).sum().item())
    sub_mismatch = int((sub_from_csv != base_subids).sum().item())
    print(f"  enzyme_id full-match:    {N - enz_mismatch}/{N}")
    print(f"  substrate_id full-match: {N - sub_mismatch}/{N}")

    if enz_mismatch > 0 or sub_mismatch > 0:
        print(f"  [ERROR] join is BROKEN; aborting save for this split")
        ok_overall = False
        continue

    # dock_index sanity: every dock PDB file must exist
    unique_docks = set(dock_indices.tolist())
    print(f"  unique dock_indices in split: {len(unique_docks)}")
    missing_pdb = [d for d in unique_docks if not (POCKET_DIR / f"{d}.pdb").exists()]
    print(f"  missing pocket PDB files: {len(missing_pdb)}")
    if missing_pdb:
        print(f"    first 5: {missing_pdb[:5]}")
        print(f"  [ERROR] some dock PDB files do not exist")
        ok_overall = False
        continue

    # 20 random spot checks (full chain: k -> sid -> (csv row) -> (enz,sub,dock,pdb))
    pick = random.sample(range(N), N_SPOT)
    for k in pick:
        sid = int(sids[k])
        did = int(dock_indices[k])
        assert (POCKET_DIR / f"{did}.pdb").exists(), f"PDB for dock {did} missing"
        assert int(base_eids[k]) == int(enz_from_csv[k])
        assert int(base_subids[k]) == int(sub_from_csv[k])
    print(f"  spot check 20/20 OK")

    # Save
    sidecar = {
        "sample_ids":   sids.clone().to(torch.int32),
        "dock_indices": dock_indices,
        "source_csv":   str(csv_path),
        "source_index": str(base_index_path),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(sidecar, out_path)
    print(f"  wrote: {out_path}")
    print(f"  file size: {out_path.stat().st_size} bytes")

print()
print("=" * 60)
if ok_overall:
    print("STEP 2 DONE. All 3 dock sidecars written.")
else:
    print("STEP 2 FAILED. Do NOT proceed.")
print("=" * 60)
import sys
sys.exit(0 if ok_overall else 1)
'''


def main():
    import subprocess, sys
    print("[local driver] running dock sidecar build on server...")
    result = subprocess.run(
        ["ssh", "-p", "35822", "root@connect.bjb2.seetacloud.com",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=300,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
