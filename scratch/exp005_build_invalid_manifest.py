"""
Build gvp_invalid_docks.pt manifest from failures.log produced by the GVP cache builder.

Per-dock metadata:
  dock_index -> {
    "reason": str,         # e.g., "empty_pocket" or "too_few_residues=1"
    "enzyme_id": int,
    "split": "train" | "val" | "test",
    "uniprot": str,
  }

Cleanest audit trail per codex recommendation. No placeholder .pt files.
"""

REMOTE_SCRIPT = r'''# -*- coding: utf-8 -*-
import csv, json
from pathlib import Path
from collections import defaultdict
import torch

BASE = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
ENZ_CSV = BASE / "data/Enzymes.csv"
OVERLAY = BASE / "data/pt_cache_dualgraph_allfix_unified"
GVP_CACHE = OVERLAY / "gvp_cache"
FAILURES_LOG = GVP_CACHE / "failures.log"
INVALID_MANIFEST = GVP_CACHE / "gvp_invalid_docks.pt"
SIDECAR_PATHS = {
    "train": OVERLAY / "random/train/dock_sidecar.pt",
    "val":   OVERLAY / "random/val/dock_sidecar.pt",
    "test":  OVERLAY / "random/test/dock_sidecar.pt",
}

# Load Enzymes.csv for uniprot lookup
with open(ENZ_CSV, encoding="utf-8-sig") as f:
    enz_rows = list(csv.DictReader(f))

# Build dock -> split from sidecars
dock_to_split = {}
dock_to_enz = {}
for split, sp in SIDECAR_PATHS.items():
    sc = torch.load(str(sp), weights_only=False)
    bi = torch.load(str(sp).replace("dock_sidecar.pt", "index.pt"), weights_only=False)
    dids = sc["dock_indices"].tolist()
    eids = bi["enzyme_ids"].tolist()
    for e, d in zip(eids, dids):
        di = int(d)
        if di not in dock_to_split:
            dock_to_split[di] = split
            dock_to_enz[di] = int(e)

print(f"sidecar covers {len(dock_to_split)} unique docks across 3 splits")

# Read failures.log
failures = []
with open(FAILURES_LOG) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        failures.append(json.loads(line))
print(f"failures.log has {len(failures)} entries")

# Build manifest
manifest = {}
by_reason = {}
for fr in failures:
    dock = int(fr["dock"])
    enz_id = int(fr["enz"])
    reason = fr["status"]
    uniprot = enz_rows[enz_id]["uniprots"].strip() if enz_id < len(enz_rows) else ""
    split = dock_to_split.get(dock, "unknown")
    manifest[dock] = {
        "reason": reason,
        "enzyme_id": enz_id,
        "split": split,
        "uniprot": uniprot,
    }
    by_reason.setdefault(reason, []).append(dock)

# Verify: every failed dock is in a training split
unknown = [d for d, m in manifest.items() if m["split"] == "unknown"]
if unknown:
    print(f"[WARN] {len(unknown)} failed docks not in any training split: {unknown[:5]}...")

# Save
torch.save(manifest, str(INVALID_MANIFEST))
print(f"\nsaved: {INVALID_MANIFEST}")
print(f"total failed docks: {len(manifest)}")
print(f"by reason:")
for r, ds in sorted(by_reason.items(), key=lambda x: -len(x[1])):
    print(f"  {r}: {len(ds)}")
print(f"\nby split:")
by_split = defaultdict(int)
for m in manifest.values():
    by_split[m["split"]] += 1
for s, n in sorted(by_split.items()):
    print(f"  {s}: {n}")

# Sample 5 entries
print(f"\nfirst 5 entries:")
for d in sorted(manifest.keys())[:5]:
    print(f"  dock={d}: {manifest[d]}")
'''


def main():
    import subprocess, sys
    print("[local driver] building invalid manifest on server...")
    result = subprocess.run(
        ["ssh", "autodl-4x5090-bj",
         "export PATH=/root/miniconda3/bin:$PATH && python"],
        input=REMOTE_SCRIPT.encode("utf-8"),
        capture_output=True,
        timeout=120,
    )
    sys.stdout.write(result.stdout.decode("utf-8", errors="replace"))
    if result.stderr:
        sys.stderr.write("\n[STDERR]\n")
        sys.stderr.write(result.stderr.decode("utf-8", errors="replace"))
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
