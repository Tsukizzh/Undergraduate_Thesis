#!/usr/bin/env python3
"""Audit downloaded RCSB mmCIF files for EXP002 Fe/HEM preparation.

This script is read-only for raw mmCIF files. It writes versioned manifest files
that say which candidate structures truly contain HEM/Fe atoms and whether the
signal is on the mapped author chain.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


HEME_COMPONENTS = {
    "HEM",
    "HEC",
    "HEA",
    "HEB",
    "HEO",
    "HAS",
    "HDD",
    "DHE",
}
IRON_COMPONENTS = {"FE", "FES", "SF4", "F3S", "FEO"}
HEME_OR_FE_COMPONENTS = HEME_COMPONENTS | IRON_COMPONENTS


def split_chain_ids(value) -> set[str]:
    if value is None:
        return set()
    if isinstance(value, float) and math.isnan(value):
        return set()
    out = set()
    for part in str(value).replace(",", ";").split(";"):
        part = part.strip()
        if part and part.lower() != "nan":
            out.add(part)
    return out


def as_list(value) -> list[str]:
    if isinstance(value, list):
        return [str(x) for x in value]
    return [str(value)]


def get_loop(mm: dict, key: str, n: int, default: str = "") -> list[str]:
    if key not in mm:
        return [default] * n
    vals = as_list(mm[key])
    if len(vals) == n:
        return vals
    if len(vals) == 1:
        return vals * n
    raise ValueError(f"Loop length mismatch for {key}: got {len(vals)}, expected {n}")


def safe_float(value: str):
    try:
        return float(value)
    except Exception:
        return None


def audit_one(cif_path: Path, target_auth_chains: set[str]) -> dict:
    row = {
        "entry_id": cif_path.stem.upper(),
        "target_auth_asym_ids": ";".join(sorted(target_auth_chains)),
        "file_path": str(cif_path),
        "exists": cif_path.exists(),
        "parse_ok": False,
        "n_atoms": 0,
        "n_target_atoms": 0,
        "n_protein_atoms_total": 0,
        "n_protein_atoms_target_chain": 0,
        "heme_comp_ids_present": "",
        "iron_comp_ids_present": "",
        "n_heme_atoms_total": 0,
        "n_heme_atoms_target_chain": 0,
        "n_fe_atoms_total": 0,
        "n_fe_atoms_target_chain": 0,
        "heme_auth_asym_ids": "",
        "fe_auth_asym_ids": "",
        "target_has_heme": False,
        "target_has_fe": False,
        "target_has_heme_fe": False,
        "any_has_heme": False,
        "any_has_fe": False,
        "first_fe_x": "",
        "first_fe_y": "",
        "first_fe_z": "",
        "error": "",
    }
    if not cif_path.exists():
        row["error"] = "missing_file"
        return row

    try:
        mm = MMCIF2Dict(str(cif_path))
        group = as_list(mm["_atom_site.group_PDB"])
        n = len(group)
        comp = get_loop(mm, "_atom_site.label_comp_id", n)
        auth_chain = get_loop(mm, "_atom_site.auth_asym_id", n)
        type_symbol = get_loop(mm, "_atom_site.type_symbol", n)
        x_vals = get_loop(mm, "_atom_site.Cartn_x", n)
        y_vals = get_loop(mm, "_atom_site.Cartn_y", n)
        z_vals = get_loop(mm, "_atom_site.Cartn_z", n)
    except Exception as exc:  # noqa: BLE001
        row["error"] = f"parse_error: {exc}"
        return row

    heme_comp_seen: set[str] = set()
    iron_comp_seen: set[str] = set()
    heme_chains: set[str] = set()
    fe_chains: set[str] = set()
    first_fe = None

    for i in range(n):
        g = group[i].strip().upper()
        c = comp[i].strip().upper()
        ch = auth_chain[i].strip()
        elem = type_symbol[i].strip().upper()
        in_target = (not target_auth_chains) or ch in target_auth_chains

        row["n_atoms"] += 1
        if in_target:
            row["n_target_atoms"] += 1
        if g == "ATOM":
            row["n_protein_atoms_total"] += 1
            if in_target:
                row["n_protein_atoms_target_chain"] += 1

        is_heme = c in HEME_COMPONENTS
        is_fe_atom = elem == "FE" or c in IRON_COMPONENTS

        if is_heme:
            heme_comp_seen.add(c)
            heme_chains.add(ch)
            row["n_heme_atoms_total"] += 1
            if in_target:
                row["n_heme_atoms_target_chain"] += 1
        if is_fe_atom:
            iron_comp_seen.add(c)
            fe_chains.add(ch)
            row["n_fe_atoms_total"] += 1
            if in_target:
                row["n_fe_atoms_target_chain"] += 1
            if first_fe is None:
                first_fe = (safe_float(x_vals[i]), safe_float(y_vals[i]), safe_float(z_vals[i]))

    row["parse_ok"] = True
    row["heme_comp_ids_present"] = ";".join(sorted(heme_comp_seen))
    row["iron_comp_ids_present"] = ";".join(sorted(iron_comp_seen))
    row["heme_auth_asym_ids"] = ";".join(sorted(ch for ch in heme_chains if ch))
    row["fe_auth_asym_ids"] = ";".join(sorted(ch for ch in fe_chains if ch))
    row["target_has_heme"] = row["n_heme_atoms_target_chain"] > 0
    row["target_has_fe"] = row["n_fe_atoms_target_chain"] > 0
    row["target_has_heme_fe"] = row["target_has_heme"] and row["target_has_fe"]
    row["any_has_heme"] = row["n_heme_atoms_total"] > 0
    row["any_has_fe"] = row["n_fe_atoms_total"] > 0
    if first_fe is not None:
        row["first_fe_x"], row["first_fe_y"], row["first_fe_z"] = first_fe
    return row


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-mmcif-dir", required=True)
    parser.add_argument("--candidate-csv", required=True)
    parser.add_argument("--best1-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tag", default="v1_20260527")
    args = parser.parse_args()

    raw_dir = Path(args.raw_mmcif_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidate_csv)
    best1 = pd.read_csv(args.best1_csv)

    entry_to_chains: dict[str, set[str]] = defaultdict(set)
    for _, r in candidates.iterrows():
        entry = str(r.get("entry_id", "")).upper().strip()
        if not entry or entry == "NAN":
            continue
        entry_to_chains[entry].update(split_chain_ids(r.get("auth_asym_ids")))

    all_cifs = sorted(raw_dir.glob("*.cif"))
    rows = []
    t0 = time.time()
    for i, cif_path in enumerate(all_cifs, start=1):
        entry = cif_path.stem.upper()
        rows.append(audit_one(cif_path, entry_to_chains.get(entry, set())))
        if i % 100 == 0:
            print(f"[audit] {i}/{len(all_cifs)}")

    audit_csv = out_dir / f"exp002_mmcif_heme_fe_atom_audit_{args.tag}.csv"
    audit_json = out_dir / f"exp002_mmcif_heme_fe_atom_audit_{args.tag}.json"
    write_csv(audit_csv, rows)

    audit_df = pd.DataFrame(rows)
    best = best1.copy()
    audit_small = audit_df.rename(columns={"entry_id": "best_entry_id"})
    best["best_entry_id"] = best["best_entry_id"].astype(str).str.upper()
    enriched = best.merge(audit_small, on="best_entry_id", how="left", suffixes=("", "_audit"))
    enriched["exp002_atom_audit_usable_target_heme_fe"] = (
        enriched["exists"].fillna(False).astype(bool)
        & enriched["parse_ok"].fillna(False).astype(bool)
        & enriched["target_has_heme_fe"].fillna(False).astype(bool)
    )
    enriched_csv = out_dir / f"exp002_best1_with_mmcif_atom_audit_{args.tag}.csv"
    enriched.to_csv(enriched_csv, index=False)

    summary = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "raw_mmcif_dir": str(raw_dir),
        "candidate_csv": str(args.candidate_csv),
        "best1_csv": str(args.best1_csv),
        "n_cif_files": len(all_cifs),
        "n_audit_rows": len(rows),
        "parse_ok": int(audit_df["parse_ok"].sum()) if len(audit_df) else 0,
        "missing_files": int((~audit_df["exists"]).sum()) if len(audit_df) else 0,
        "any_has_heme": int(audit_df["any_has_heme"].sum()) if len(audit_df) else 0,
        "any_has_fe": int(audit_df["any_has_fe"].sum()) if len(audit_df) else 0,
        "any_has_heme_fe": int((audit_df["any_has_heme"] & audit_df["any_has_fe"]).sum()) if len(audit_df) else 0,
        "target_has_heme_fe": int(audit_df["target_has_heme_fe"].sum()) if len(audit_df) else 0,
        "best1_rows": int(len(enriched)),
        "best1_usable_target_heme_fe_rows": int(enriched["exp002_atom_audit_usable_target_heme_fe"].sum()),
        "best1_usable_unique_uniprot": int(enriched.loc[enriched["exp002_atom_audit_usable_target_heme_fe"], "uniprot"].nunique()),
        "best1_usable_positive_samples": int(enriched.loc[enriched["exp002_atom_audit_usable_target_heme_fe"], "n_positive"].fillna(0).sum()),
        "elapsed_sec": round(time.time() - t0, 2),
        "outputs": {
            "audit_csv": str(audit_csv),
            "audit_json": str(audit_json),
            "best1_enriched_csv": str(enriched_csv),
        },
    }
    audit_json.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
