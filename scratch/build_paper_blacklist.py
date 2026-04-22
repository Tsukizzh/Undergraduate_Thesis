"""
Build blacklist of enzyme_global_ids in our P450 dataset that overlap with
the 389 ESIBank P450 training set.

Inputs:
  - scratch/P450_Enzymes_server.csv  (1622 rows, columns: Protein sequence, uniprots)
  - 389个P450的PDB映射_完整版.csv     (389 UniProt IDs from ESIBank training)

Outputs:
  - scratch/paper_blacklist.json  (list of enzyme_global_ids to drop)
  - scratch/paper_blacklist_report.json  (detailed stats)
"""
from __future__ import annotations
import csv, json, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ENZ_CSV = ROOT / "scratch/P450_Enzymes_server.csv"
ESIBANK_CSV = ROOT / "毕业设计/P450_EZSpecificity_研究项目/PathA_2026-01-08_模型评估测试集构建/source_data/04_PDB映射/389个P450的PDB映射_完整版.csv"

OUT_BLACKLIST = ROOT / "scratch/paper_blacklist.json"
OUT_REPORT = ROOT / "scratch/paper_blacklist_report.json"


def main():
    # Load 389 ESIBank UniProt set
    with open(ESIBANK_CSV, encoding="utf-8-sig") as f:
        esibank = {row["uniprot_id"].strip() for row in csv.DictReader(f) if row["uniprot_id"].strip()}
    assert len(esibank) == 389, f"Expected 389 ESIBank P450, got {len(esibank)}"

    # Load our P450 Enzymes.csv (row_idx == enzyme_global_id per allfix convention)
    with open(ENZ_CSV, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1622, f"Expected 1622 enzymes, got {len(rows)}"

    # Classify rows
    blacklist: list[int] = []
    matched_pairs: list[dict] = []
    fake_id_rows: list[int] = []  # ENZ_G* synthetic IDs (PlantP450DB etc.)
    unmatched_real: list[int] = []

    for idx, row in enumerate(rows):
        uniprot = row["uniprots"].strip()
        seq = row["Protein sequence"].strip()
        if uniprot.startswith("ENZ_"):
            fake_id_rows.append(idx)
            continue
        if uniprot in esibank:
            blacklist.append(idx)
            matched_pairs.append({
                "enzyme_global_id": idx,
                "uniprot": uniprot,
                "seq_len": len(seq),
                "seq_sha256_8": hashlib.sha256(seq.encode()).hexdigest()[:8],
            })
        else:
            unmatched_real.append(idx)

    # Collision check: any ESIBank UniProt ID NOT in our dataset?
    our_uniprots = {r["uniprots"].strip() for r in rows if not r["uniprots"].startswith("ENZ_")}
    esibank_not_in_ours = sorted(esibank - our_uniprots)

    # Save outputs
    with open(OUT_BLACKLIST, "w", encoding="utf-8") as f:
        json.dump({
            "blacklisted_enzyme_global_ids": sorted(blacklist),
            "count": len(blacklist),
        }, f, indent=2)

    report = {
        "source_files": {
            "our_enzymes_csv": str(ENZ_CSV.name),
            "esibank_p450_csv": str(ESIBANK_CSV.name),
        },
        "esibank_p450_count": len(esibank),
        "our_enzymes_count": len(rows),
        "our_real_uniprot_count": len(rows) - len(fake_id_rows),
        "our_fake_id_count": len(fake_id_rows),
        "blacklist_count": len(blacklist),
        "blacklist_fraction_of_ours": round(len(blacklist) / len(rows), 4),
        "blacklist_fraction_of_esibank": round(len(blacklist) / len(esibank), 4),
        "esibank_not_in_our_dataset_count": len(esibank_not_in_ours),
        "esibank_not_in_our_dataset_ids": esibank_not_in_ours[:20],  # first 20 for visibility
        "matched_pairs_first5": matched_pairs[:5],
        "fake_id_rows_first10": fake_id_rows[:10],
        "unmatched_real_count": len(unmatched_real),
    }
    with open(OUT_REPORT, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)

    print("=" * 70)
    print("Paper Blacklist Build — Summary")
    print("=" * 70)
    print(f"ESIBank P450 total:           {len(esibank)}")
    print(f"Our P450 enzymes total:       {len(rows)}")
    print(f"  with real UniProt ID:       {len(rows) - len(fake_id_rows)}")
    print(f"  with synthetic ENZ_* ID:    {len(fake_id_rows)}")
    print(f"")
    print(f"Blacklist size (overlap):     {len(blacklist)}")
    print(f"  {len(blacklist)}/{len(esibank)} = {100*len(blacklist)/len(esibank):.1f}% of ESIBank P450")
    print(f"  {len(blacklist)}/{len(rows)} = {100*len(blacklist)/len(rows):.1f}% of our enzymes")
    print(f"")
    print(f"ESIBank P450 NOT in our data: {len(esibank_not_in_ours)}")
    print(f"  (these paper-training enzymes are not in our dataset, so cannot cause leakage)")
    print(f"")
    print(f"Outputs:")
    print(f"  {OUT_BLACKLIST}")
    print(f"  {OUT_REPORT}")


if __name__ == "__main__":
    main()
