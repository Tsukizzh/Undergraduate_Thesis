#!/usr/bin/env python3
"""
Build dock_index_mapping.csv for Step 8 structure feature extraction.

This script creates a mapping from Dock_Index to (pdb_id, ligand_ccd, record_id)
which is needed for correctly extracting individual ligands from multi-ligand PDB files.

Usage:
    python build_dock_mapping.py --base-dir <PathA_root>
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

# Reuse helpers from build_datasets.py
UTF8_SIG = "utf-8-sig"


class CanonicalRecord(NamedTuple):
    record_id: str
    pdb_id: str
    ligand_ccd: str
    uniprot_id: str
    final_class: str


def _norm(s: str) -> str:
    return s.strip().upper()


def _read_csv_dicts(path: Path) -> List[Dict[str, str]]:
    with open(path, "r", encoding=UTF8_SIG, newline="") as f:
        return list(csv.DictReader(f))


def _load_chunks_as_canonical(chunks_dir: Path) -> Dict[str, Dict[str, str]]:
    """Load all chunk_*.csv files and build record_id -> {pdb_id, ligand_ccd, uniprot_id}."""
    out: Dict[str, Dict[str, str]] = {}
    for chunk_file in sorted(chunks_dir.glob("chunk_*.csv")):
        rows = _read_csv_dicts(chunk_file)
        for r in rows:
            rid = r.get("record_id", "").strip()
            if not rid:
                continue
            out[rid] = {
                "pdb_id": _norm(r.get("pdb_id", "")),
                "ligand_ccd": _norm(r.get("ligand_ccd", "")),
                "uniprot_id": r.get("uniprot_id", "").strip(),
            }
    return out


def _load_merged_final_class(merged_csv: Path) -> Dict[str, str]:
    """Load merged.csv and build record_id -> final_class map."""
    rows = _read_csv_dicts(merged_csv)
    out: Dict[str, str] = {}
    for r in rows:
        rid = r.get("record_id", "").strip()
        final_class = _norm(r.get("final_class", ""))
        if rid and final_class:
            out[rid] = final_class
    return out


def _load_smiles_map(ccd_smiles_csv: Path) -> Dict[str, str]:
    """Load CCD -> SMILES mapping."""
    rows = _read_csv_dicts(ccd_smiles_csv)
    out: Dict[str, str] = {}
    for r in rows:
        # Try multiple possible column names
        ccd = _norm(r.get("ccd_code", "") or r.get("ligand_ccd", "") or r.get("ccd", ""))
        smiles = r.get("smiles", "").strip()
        if ccd and smiles:
            out[ccd] = smiles
    return out


def _load_enzyme_mapping(mapping_csv: Path) -> Dict[Tuple[str, str], int]:
    """Load (pdb_id, ligand_ccd) -> Enzyme_Index mapping."""
    rows = _read_csv_dicts(mapping_csv)
    out: Dict[Tuple[str, str], int] = {}
    for r in rows:
        pdb_id = _norm(r.get("pdb_id", ""))
        ligand_ccd = _norm(r.get("ligand_ccd", ""))
        try:
            enzyme_idx = int(r.get("Enzyme_Index", ""))
        except ValueError:
            continue
        if pdb_id and ligand_ccd:
            out[(pdb_id, ligand_ccd)] = enzyme_idx
    return out


def _load_enzyme_index_set(enzymes_csv: Path) -> Set[int]:
    """Load valid Enzyme_Index values from Enzymes.csv."""
    rows = _read_csv_dicts(enzymes_csv)
    return {i for i, _ in enumerate(rows)}


def _default_paths(base_dir: Path) -> Dict[str, Path]:
    return {
        "chunks_dir": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证" / "01_输入数据_chunks",
        "merged_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "04_配体三方深度验证" / "03_合并结果" / "三方深度验证结果_merged.csv",
        "ccd_smiles_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "05_提取SMILES" / "ccd_smiles.csv",
        "enzymes_csv": base_dir / "data" / "02_Step2_酶序列" / "Enzymes.csv",
        "record_enzyme_mapping_csv": base_dir / "data" / "02_Step2_酶序列" / "record_enzyme_mapping.csv",
        "output_dir": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "06_数据集切分",
    }


# Scheme definitions (must match build_datasets.py)
B3_INCLUDED = {"SUBSTRATE", "INHIBITOR", "PRODUCT"}
B3_POSITIVE = {"SUBSTRATE"}
B3_NEGATIVE = {"INHIBITOR", "PRODUCT"}


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Build dock_index_mapping.csv for Step 8.")
    p.add_argument("--base-dir", default=".", help="Base directory (PathA root)")
    p.add_argument("--output", default=None, help="Output CSV path (default: output_dir/dock_index_mapping.csv)")
    args = p.parse_args(argv)

    base_dir = Path(args.base_dir).resolve()
    defaults = _default_paths(base_dir)

    chunks_dir = defaults["chunks_dir"]
    merged_csv = defaults["merged_csv"]
    ccd_smiles_csv = defaults["ccd_smiles_csv"]
    enzymes_csv = defaults["enzymes_csv"]
    mapping_csv = defaults["record_enzyme_mapping_csv"]
    output_dir = defaults["output_dir"]

    output_path = Path(args.output) if args.output else output_dir / "dock_index_mapping.csv"

    print(f"[INFO] chunks_dir={chunks_dir}")
    print(f"[INFO] merged_csv={merged_csv}")
    print(f"[INFO] output_path={output_path}")

    # Load data
    chunk_map = _load_chunks_as_canonical(chunks_dir)
    class_map = _load_merged_final_class(merged_csv)
    smiles_by_ccd = _load_smiles_map(ccd_smiles_csv)
    enzyme_by_key = _load_enzyme_mapping(mapping_csv)
    valid_enzyme_indices = _load_enzyme_index_set(enzymes_csv)

    print(f"[INFO] Loaded {len(chunk_map)} records from chunks")
    print(f"[INFO] Loaded {len(class_map)} final_class from merged")
    print(f"[INFO] Loaded {len(smiles_by_ccd)} SMILES mappings")
    print(f"[INFO] Loaded {len(enzyme_by_key)} enzyme mappings")

    # Build canonical records (same logic as build_datasets.py)
    common_ids = set(chunk_map.keys()) & set(class_map.keys())
    canonical_records: List[CanonicalRecord] = []
    for rid in sorted(common_ids):
        cm = chunk_map[rid]
        canonical_records.append(CanonicalRecord(
            record_id=rid,
            pdb_id=cm["pdb_id"],
            ligand_ccd=cm["ligand_ccd"],
            uniprot_id=cm["uniprot_id"],
            final_class=class_map[rid],
        ))

    print(f"[INFO] {len(canonical_records)} canonical records")

    # Filter by B3 scheme (SUBSTRATE + INHIBITOR + PRODUCT)
    included_records = [r for r in canonical_records if r.final_class in B3_INCLUDED]
    print(f"[INFO] {len(included_records)} records in B3 scheme")

    # Build Dock_Index mapping (same order as build_datasets.py)
    # Track unique SMILES for Substrate_Index
    seen_smiles: Dict[str, int] = {}
    substrates: List[str] = []

    mapping_rows: List[Tuple[int, str, str, str, int, int, int]] = []
    # (Dock_Index, record_id, pdb_id, ligand_ccd, Enzyme_Index, Substrate_Index, Label)

    for rec in included_records:
        smi = smiles_by_ccd.get(rec.ligand_ccd, "")
        if not smi:
            continue

        if smi not in seen_smiles:
            seen_smiles[smi] = len(substrates)
            substrates.append(smi)

        substrate_idx = seen_smiles[smi]
        enzyme_key = (rec.pdb_id, rec.ligand_ccd)
        enzyme_idx = enzyme_by_key.get(enzyme_key)

        if enzyme_idx is None:
            continue

        if enzyme_idx not in valid_enzyme_indices:
            continue

        # Determine label
        if rec.final_class in B3_POSITIVE:
            label = 1
        elif rec.final_class in B3_NEGATIVE:
            label = 0
        else:
            continue

        dock_idx = len(mapping_rows)
        mapping_rows.append((
            dock_idx,
            rec.record_id,
            rec.pdb_id,
            rec.ligand_ccd,
            enzyme_idx,
            substrate_idx,
            label,
        ))

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Dock_Index", "record_id", "pdb_id", "ligand_ccd",
            "Enzyme_Index", "Substrate_Index", "Label"
        ])
        for row in mapping_rows:
            writer.writerow(row)

    print(f"[INFO] Wrote {output_path} ({len(mapping_rows)} rows)")

    # Validate against data.csv
    data_csv = output_dir / "B3_完整数据集_272pos267neg" / "data.csv"
    if data_csv.exists():
        with open(data_csv, "r", encoding=UTF8_SIG) as f:
            data_rows = list(csv.DictReader(f))
        if len(data_rows) == len(mapping_rows):
            print(f"[OK] Row count matches data.csv ({len(data_rows)})")
            # Verify indices match
            mismatches = 0
            for i, (data_row, map_row) in enumerate(zip(data_rows, mapping_rows)):
                if int(data_row["Dock Index"]) != map_row[0]:
                    mismatches += 1
                if int(data_row["Enzyme Index"]) != map_row[4]:
                    mismatches += 1
                if int(data_row["Substrate Index"]) != map_row[5]:
                    mismatches += 1
                if int(data_row["Label"]) != map_row[6]:
                    mismatches += 1
            if mismatches == 0:
                print("[OK] All indices match data.csv")
            else:
                print(f"[ERROR] {mismatches} index mismatches with data.csv!")
                return 1
        else:
            print(f"[ERROR] Row count mismatch: mapping={len(mapping_rows)}, data.csv={len(data_rows)}")
            return 1

    # Report multi-ligand PDBs
    pdb_counts: Dict[str, int] = {}
    for row in mapping_rows:
        pdb_id = row[2]
        pdb_counts[pdb_id] = pdb_counts.get(pdb_id, 0) + 1

    multi_ligand = [(pdb, cnt) for pdb, cnt in pdb_counts.items() if cnt > 1]
    if multi_ligand:
        print(f"\n[INFO] Multi-ligand PDBs ({len(multi_ligand)}):")
        for pdb, cnt in sorted(multi_ligand, key=lambda x: -x[1])[:10]:
            print(f"  {pdb}: {cnt} ligands")

    return 0


if __name__ == "__main__":
    sys.exit(main())
