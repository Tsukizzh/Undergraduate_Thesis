#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from Bio.SVDSuperimposer import SVDSuperimposer

sys.path.insert(0, "/root/autodl-tmp/EZSpecificity/PathD/P450/scripts")

from q01_exp002_make_fe_overlay_samples_20260527 import read_official_pdb_ca, sequence_align_pairs
from q01_exp003_make_alphafill_heme_overlay_samples_20260528 import read_alphafill_ca_and_heme


BASE = Path("/root/autodl-tmp/EZSpecificity/PathD/P450/data/q01_fe_embedding_patch/exp003_full_p450_fe_heme_overlay")
OFFICIAL_PDB = BASE / "official_esibank_p450_pdb_20260528_package/pdb"
OFFICIAL_MANIFEST = BASE / "official_esibank_p450_pdb_20260528_package/exp003_p450_official_esibank_pdb_manifest_20260528.csv"
ALPHAFILL_SUMMARY = BASE / "alphafill_p450_216_20260528/manifests/alphafill_hem_candidate_summary_20260528.csv"
ALPHAFILL_CIF = BASE / "alphafill_p450_216_20260528/cif"
OUT_JSON = BASE / "alphafill_p450_216_20260528/manifests/p0dki7_homolog_template_search_20260528.json"


def best_for_pdb(pdb_path: Path, alpha_rows: list[dict]) -> dict | None:
    official = read_official_pdb_ca(pdb_path)
    best = None
    for official_chain, official_records in official.items():
        if len(official_records) < 80:
            continue
        for row in alpha_rows:
            cif = ALPHAFILL_CIF / f"{row['uniprot']}.cif"
            try:
                alphafill_ca, heme_atoms = read_alphafill_ca_and_heme(cif, row["best_asym_id"])
            except Exception:
                continue
            if not heme_atoms:
                continue
            for af_chain, af_records in alphafill_ca.items():
                if len(af_records) < 80:
                    continue
                pairs, score = sequence_align_pairs(official_records, af_records)
                if len(pairs) < 80:
                    continue
                matches = sum(
                    1
                    for i, j in pairs
                    if official_records[i][0] == af_records[j][0] and official_records[i][0] != "X"
                )
                identity = matches / len(pairs)
                coverage = len(pairs) / max(1, len(official_records))
                official_xyz = np.array([official_records[i][1] for i, _ in pairs], dtype=np.float64)
                af_xyz = np.array([af_records[j][1] for _, j in pairs], dtype=np.float64)
                sup = SVDSuperimposer()
                sup.set(official_xyz, af_xyz)
                sup.run()
                rmsd = float(sup.get_rms())
                candidate = {
                    "template_uniprot": row["uniprot"],
                    "template_pdb": row["best_pdb_id"],
                    "heme_asym_id": row["best_asym_id"],
                    "official_chain": official_chain,
                    "alphafill_chain": af_chain,
                    "pairs": len(pairs),
                    "identity": round(identity, 6),
                    "coverage": round(coverage, 6),
                    "rmsd": round(rmsd, 6),
                    "score": round(float(score), 3),
                    "atom_count": row["best_atom_count"],
                }
                key = (identity, coverage, len(pairs), -rmsd)
                if best is None or key > best[0]:
                    best = (key, candidate)
    return best[1] if best else None


def main() -> int:
    official_rows = list(csv.DictReader(OFFICIAL_MANIFEST.open()))
    p0_rows = [
        row
        for row in official_rows
        if row["p450_uniprots"] == "P0DKI7"
        and row["is_trainable_in_exp003_manifest"].lower() == "true"
        and (row["found_in_drive_metadata"].lower() == "true" or row["found_in_drivefs_metadata"].lower() == "true")
    ]
    alpha_rows = [
        row
        for row in csv.DictReader(ALPHAFILL_SUMMARY.open())
        if int(row.get("hem_candidates") or 0) > 0
    ]
    results = []
    for row in p0_rows:
        dock_index = row["structure_index"]
        best = best_for_pdb(OFFICIAL_PDB / f"{dock_index}.pdb", alpha_rows)
        results.append({"dock_index": dock_index, **(best or {})})

    OUT_JSON.write_text(json.dumps(results, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"p0dki7_rows": len(p0_rows), "results": results}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
