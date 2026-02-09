# -*- coding: utf-8 -*-
"""
Build B1/B2/B3 datasets from P450 enzyme-ligand classification results.
VERSION 2: Fixed record_id alignment bug - use chunk files as canonical source.

Data Flow (CORRECT):
  chunks/*.csv (record_id → pdb_id, ligand_ccd, uniprot_id)
       ↓ JOIN on record_id
  merged.csv (record_id → final_class)
       ↓
  ccd_smiles.csv (ligand_ccd → SMILES)
       ↓
  record_enzyme_mapping.csv ((pdb_id, ligand_ccd) → Enzyme_Index)
       ↓
  Output: Substrates.csv, B1/B2/B3 data.csv

Schemes:
  - B1: SUBSTRATE only (Label=1)
  - B2: SUBSTRATE (Label=1) + PRODUCT (Label=0)
  - B3: SUBSTRATE (Label=1) + PRODUCT (Label=0) + INHIBITOR (Label=0)
  - B4: PRODUCT only (Label=0)
  - B5: INHIBITOR only (Label=0)
  - B6: SUBSTRATE (Label=1) + INHIBITOR (Label=0)
  - EXCLUDE is never included in any scheme.

Note: PRODUCT is labeled as negative (0) because products are the RESULT of catalysis,
      not substrates that can be catalyzed. This aligns with EZSpecificity's training data
      where only substrates are positive samples.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from glob import glob
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

UTF8_SIG = "utf-8-sig"


def _norm_pdb_id(v: str) -> str:
    return (v or "").strip().lower()


def _norm_ccd(v: str) -> str:
    return (v or "").strip().upper()


def _norm_class(v: str) -> str:
    return (v or "").strip().upper()


def _norm_uniprot(v: str) -> str:
    return (v or "").strip().upper()


def _read_csv_dicts(path: Path, *, encoding: str = UTF8_SIG) -> List[Dict[str, str]]:
    with path.open("r", encoding=encoding, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return [dict(row) for row in reader]


def _write_csv(path: Path, header: Sequence[str], rows, *, encoding: str = UTF8_SIG) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding=encoding, newline="") as f:
        w = csv.writer(f)
        w.writerow(list(header))
        for r in rows:
            w.writerow(list(r))


@dataclass(frozen=True)
class Scheme:
    name: str
    included: Set[str]
    positive: Set[str]
    negative: Set[str]


SCHEMES: List[Scheme] = [
    Scheme(name="B1", included={"SUBSTRATE"}, positive={"SUBSTRATE"}, negative=set()),
    Scheme(name="B2", included={"SUBSTRATE", "PRODUCT"}, positive={"SUBSTRATE"}, negative={"PRODUCT"}),
    Scheme(name="B3", included={"SUBSTRATE", "PRODUCT", "INHIBITOR"}, positive={"SUBSTRATE"}, negative={"PRODUCT", "INHIBITOR"}),
    Scheme(name="B4", included={"PRODUCT"}, positive=set(), negative={"PRODUCT"}),
    Scheme(name="B5", included={"INHIBITOR"}, positive=set(), negative={"INHIBITOR"}),
    Scheme(name="B6", included={"SUBSTRATE", "INHIBITOR"}, positive={"SUBSTRATE"}, negative={"INHIBITOR"}),
]

# Map scheme names to descriptive folder names
# Format: BX_内容描述_正样本数pos负样本数neg
SCHEME_FOLDER_NAMES = {
    "B1": "B1_仅底物_271pos",
    "B2": "B2_底物正产物负_271pos23neg",
    "B3": "B3_完整数据集_271pos268neg",
    "B4": "B4_仅产物_23neg",
    "B5": "B5_仅抑制剂_245neg",
    "B6": "B6_底物正抑制剂负_271pos245neg",
}


@dataclass
class CanonicalRecord:
    """A record with all necessary fields from chunk + merged."""
    record_id: str
    pdb_id: str
    ligand_ccd: str
    uniprot_id: str
    final_class: str


def _load_chunks_as_canonical(chunks_dir: Path) -> Dict[str, Dict[str, str]]:
    """Load all chunk_*.csv and build record_id -> {pdb_id, ligand_ccd, uniprot_id} map."""
    chunk_files = sorted(glob(str(chunks_dir / "chunk_*.csv")))
    if not chunk_files:
        raise FileNotFoundError(f"No chunk_*.csv files found in {chunks_dir}")

    record_map: Dict[str, Dict[str, str]] = {}
    for cf in chunk_files:
        rows = _read_csv_dicts(Path(cf))
        for r in rows:
            rid = r.get("record_id", "").strip()
            if not rid:
                continue
            if rid in record_map:
                raise ValueError(f"Duplicate record_id in chunks: {rid}")
            record_map[rid] = {
                "pdb_id": _norm_pdb_id(r.get("pdb_id", "")),
                "ligand_ccd": _norm_ccd(r.get("ligand_ccd", "")),
                "uniprot_id": _norm_uniprot(r.get("uniprot_id", "")),
            }
    return record_map


def _load_merged_final_class(merged_csv: Path) -> Dict[str, str]:
    """Load merged.csv and build record_id -> final_class map."""
    rows = _read_csv_dicts(merged_csv)
    out: Dict[str, str] = {}
    for r in rows:
        rid = r.get("record_id", "").strip()
        final_class = _norm_class(r.get("final_class", ""))
        if rid and final_class:
            out[rid] = final_class
    return out


def _build_canonical_records(
    chunks_dir: Path,
    merged_csv: Path,
) -> Tuple[List[CanonicalRecord], List[str]]:
    """Build list of CanonicalRecord by joining chunks + merged on record_id."""
    chunk_map = _load_chunks_as_canonical(chunks_dir)
    class_map = _load_merged_final_class(merged_csv)

    warnings: List[str] = []

    # Validate record_id sets match - STRICT: must be identical
    chunk_ids = set(chunk_map.keys())
    merged_ids = set(class_map.keys())

    missing_in_merged = chunk_ids - merged_ids
    missing_in_chunks = merged_ids - chunk_ids

    if missing_in_merged or missing_in_chunks:
        msg = f"record_id mismatch: {len(missing_in_merged)} in chunks only, {len(missing_in_chunks)} in merged only"
        warnings.append(msg)
        # For this pipeline, mismatch is a data integrity issue - raise error
        raise ValueError(f"CRITICAL: {msg}. Cannot proceed with inconsistent data.")

    # Build canonical records (only for ids present in both)
    common_ids = chunk_ids & merged_ids
    records: List[CanonicalRecord] = []
    for rid in sorted(common_ids):
        cm = chunk_map[rid]
        records.append(CanonicalRecord(
            record_id=rid,
            pdb_id=cm["pdb_id"],
            ligand_ccd=cm["ligand_ccd"],
            uniprot_id=cm["uniprot_id"],
            final_class=class_map[rid],
        ))

    return records, warnings


def _load_smiles_map(ccd_smiles_csv: Path) -> Dict[str, str]:
    rows = _read_csv_dicts(ccd_smiles_csv)
    out: Dict[str, str] = {}
    for r in rows:
        ccd = _norm_ccd(r.get("ccd_code", ""))
        smi = (r.get("smiles", "") or "").strip()
        if ccd and smi:
            out.setdefault(ccd, smi)
    return out


def _load_enzyme_index_set(enzymes_csv: Path) -> Set[int]:
    rows = _read_csv_dicts(enzymes_csv)
    out: Set[int] = set()
    for r in rows:
        s = (r.get("Enzyme_Index", "") or "").strip()
        if s:
            try:
                out.add(int(s))
            except ValueError:
                pass
    return out


def _load_enzyme_mapping(mapping_csv: Path) -> Dict[Tuple[str, str], int]:
    """Load record_enzyme_mapping.csv, keyed by (pdb_id, ligand_ccd).

    Note: For robustness, we could use (pdb_id, ligand_ccd, uniprot_id) as key,
    but current mapping only has (pdb_id, ligand_ccd). We'll detect conflicts.
    """
    rows = _read_csv_dicts(mapping_csv)
    out: Dict[Tuple[str, str], int] = {}
    conflicts: List[Tuple[Tuple[str, str], int, int]] = []

    for r in rows:
        pdb_id = _norm_pdb_id(r.get("pdb_id", ""))
        ligand_ccd = _norm_ccd(r.get("ligand_ccd", ""))
        enzyme_s = (r.get("Enzyme_Index", "") or "").strip()
        if not pdb_id or not ligand_ccd or not enzyme_s:
            continue
        try:
            enzyme_idx = int(enzyme_s)
        except ValueError:
            continue

        k = (pdb_id, ligand_ccd)
        if k in out and out[k] != enzyme_idx:
            conflicts.append((k, out[k], enzyme_idx))
            continue
        out[k] = enzyme_idx

    if conflicts:
        sys.stderr.write(f"[WARN] Conflicting Enzyme_Index for {len(conflicts)} (pdb_id, ligand_ccd) keys.\n")
        for k, a, b in conflicts[:10]:
            sys.stderr.write(f"  [WARN] {k} -> {a} vs {b} (keeping {a})\n")
        if len(conflicts) > 10:
            sys.stderr.write(f"  [WARN] ... and {len(conflicts) - 10} more\n")

    return out


def _build_substrates_and_data(
    canonical_records: List[CanonicalRecord],
    smiles_by_ccd: Dict[str, str],
    enzyme_by_key: Dict[Tuple[str, str], int],
    valid_enzyme_indices: Set[int],
    scheme: Scheme,
    *,
    strict: bool,
) -> Tuple[List[str], List[Tuple[int, int, int, int]], List[str]]:
    """Build Substrates list and data rows for a given scheme.

    Returns: (substrates_smiles, data_rows, warnings)
    """
    warnings: List[str] = []

    # Filter records by scheme
    included_records = [r for r in canonical_records if r.final_class in scheme.included]

    # Build unique substrates by SMILES
    seen_smiles: Dict[str, int] = {}
    substrates: List[str] = []
    missing_smiles: List[str] = []

    for rec in included_records:
        smi = smiles_by_ccd.get(rec.ligand_ccd, "")
        if not smi:
            missing_smiles.append(f"{rec.record_id}:{rec.ligand_ccd}")
            continue
        if smi not in seen_smiles:
            seen_smiles[smi] = len(substrates)
            substrates.append(smi)

    if missing_smiles:
        msg = f"{scheme.name}: {len(missing_smiles)} records missing SMILES"
        warnings.append(msg)
        if strict:
            raise ValueError(msg + f". Examples: {missing_smiles[:5]}")

    # Build data rows
    data_rows: List[Tuple[int, int, int, int]] = []
    missing_enzyme: List[str] = []
    invalid_enzyme: List[str] = []

    for rec in included_records:
        smi = smiles_by_ccd.get(rec.ligand_ccd, "")
        if not smi:
            continue  # Already warned above

        substrate_idx = seen_smiles[smi]
        enzyme_key = (rec.pdb_id, rec.ligand_ccd)
        enzyme_idx = enzyme_by_key.get(enzyme_key)

        if enzyme_idx is None:
            missing_enzyme.append(f"{rec.record_id}:{enzyme_key}")
            continue

        if enzyme_idx not in valid_enzyme_indices:
            invalid_enzyme.append(f"{rec.record_id}:{enzyme_idx}")
            continue  # Skip invalid enzyme indices

        # Determine label
        if rec.final_class in scheme.positive:
            label = 1
        elif rec.final_class in scheme.negative:
            label = 0
        else:
            raise ValueError(f"Class {rec.final_class} is included but not labeled in {scheme.name}")

        dock_idx = len(data_rows)
        data_rows.append((dock_idx, enzyme_idx, substrate_idx, label))

    if missing_enzyme:
        msg = f"{scheme.name}: {len(missing_enzyme)} records missing Enzyme_Index"
        warnings.append(msg)
        if strict:
            raise ValueError(msg + f". Examples: {missing_enzyme[:5]}")

    if invalid_enzyme:
        warnings.append(f"{scheme.name}: {len(invalid_enzyme)} records have Enzyme_Index not in Enzymes.csv")

    return substrates, data_rows, warnings


def _default_paths(base_dir: Path) -> Dict[str, Path]:
    return {
        "chunks_dir": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "02_配体分类并行验证" / "01_输入数据_chunks",
        "merged_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "04_配体三方深度验证" / "03_合并结果" / "三方深度验证结果_merged.csv",
        "ccd_smiles_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "05_提取SMILES" / "ccd_smiles.csv",
        "enzymes_csv": base_dir / "data" / "02_Step2_酶序列" / "Enzymes.csv",
        "record_enzyme_mapping_csv": base_dir / "data" / "02_Step2_酶序列" / "record_enzyme_mapping.csv",
        "output_dir": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "06_数据集切分",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Build B1/B2/B3 datasets (v2: fixed record_id alignment).")
    p.add_argument("--base-dir", default=".", help="Base directory")
    p.add_argument("--chunks-dir", default=None, help="Directory containing chunk_*.csv files")
    p.add_argument("--merged-csv", default=None, help="Path to 三方深度验证结果_merged.csv")
    p.add_argument("--ccd-smiles-csv", default=None, help="Path to ccd_smiles.csv")
    p.add_argument("--enzymes-csv", default=None, help="Path to Enzymes.csv")
    p.add_argument("--record-enzyme-mapping-csv", default=None, help="Path to record_enzyme_mapping.csv")
    p.add_argument("--output-dir", default=None, help="Output directory")
    p.add_argument("--strict", action="store_true", help="Fail on any missing data")
    args = p.parse_args(argv)

    base_dir = Path(args.base_dir).resolve()
    defaults = _default_paths(base_dir)

    chunks_dir = Path(args.chunks_dir) if args.chunks_dir else defaults["chunks_dir"]
    merged_csv = Path(args.merged_csv) if args.merged_csv else defaults["merged_csv"]
    ccd_smiles_csv = Path(args.ccd_smiles_csv) if args.ccd_smiles_csv else defaults["ccd_smiles_csv"]
    enzymes_csv = Path(args.enzymes_csv) if args.enzymes_csv else defaults["enzymes_csv"]
    mapping_csv = Path(args.record_enzyme_mapping_csv) if args.record_enzyme_mapping_csv else defaults["record_enzyme_mapping_csv"]
    output_dir = Path(args.output_dir) if args.output_dir else defaults["output_dir"]

    print(f"[INFO] chunks_dir={chunks_dir}")
    print(f"[INFO] merged_csv={merged_csv}")
    print(f"[INFO] ccd_smiles_csv={ccd_smiles_csv}")
    print(f"[INFO] enzymes_csv={enzymes_csv}")
    print(f"[INFO] mapping_csv={mapping_csv}")
    print(f"[INFO] output_dir={output_dir}")
    print(f"[INFO] strict={args.strict}")

    # Step 1: Build canonical records from chunks + merged
    print("\n[INFO] Loading canonical records from chunks + merged...")
    canonical_records, load_warnings = _build_canonical_records(chunks_dir, merged_csv)
    print(f"[INFO] Canonical records: {len(canonical_records)}")
    for w in load_warnings:
        sys.stderr.write(f"[WARN] {w}\n")

    # Step 2: Load auxiliary data
    smiles_by_ccd = _load_smiles_map(ccd_smiles_csv)
    print(f"[INFO] SMILES map: {len(smiles_by_ccd)} CCD codes")

    enzyme_by_key = _load_enzyme_mapping(mapping_csv)
    print(f"[INFO] Enzyme mapping: {len(enzyme_by_key)} (pdb_id, ligand_ccd) keys")

    valid_enzyme_indices = _load_enzyme_index_set(enzymes_csv)
    print(f"[INFO] Valid Enzyme_Index: {len(valid_enzyme_indices)}")

    # Step 3: Build B3 first to get unified Substrates.csv
    print("\n[INFO] Building B3 (for unified Substrates.csv)...")
    b3_scheme = next(s for s in SCHEMES if s.name == "B3")
    b3_substrates, b3_rows, b3_warnings = _build_substrates_and_data(
        canonical_records, smiles_by_ccd, enzyme_by_key, valid_enzyme_indices,
        b3_scheme, strict=args.strict
    )
    for w in b3_warnings:
        sys.stderr.write(f"[WARN] {w}\n")

    # Write Substrates.csv (shared across all schemes)
    substrates_path = output_dir / "Substrates.csv"
    _write_csv(
        substrates_path,
        header=["Substrate_Index", "Substrate_SMILES"],
        rows=((i, smi) for i, smi in enumerate(b3_substrates)),
    )
    print(f"[INFO] Wrote {substrates_path} (n={len(b3_substrates)})")

    # Build SMILES -> Substrate_Index lookup for other schemes
    smiles_to_idx = {smi: i for i, smi in enumerate(b3_substrates)}

    # Write B3/data.csv
    b3_folder = SCHEME_FOLDER_NAMES.get("B3", "B3")
    b3_path = output_dir / b3_folder / "data.csv"
    _write_csv(
        b3_path,
        header=["Dock Index", "Enzyme Index", "Substrate Index", "Label"],
        rows=b3_rows,
    )
    pos = sum(1 for *_, y in b3_rows if y == 1)
    neg = sum(1 for *_, y in b3_rows if y == 0)
    print(f"[INFO] Wrote {b3_path} (n={len(b3_rows)}, pos={pos}, neg={neg})")

    # Step 4: Build B1 and B2
    for scheme in SCHEMES:
        if scheme.name == "B3":
            continue  # Already done

        print(f"\n[INFO] Building {scheme.name}...")
        included_records = [r for r in canonical_records if r.final_class in scheme.included]

        data_rows: List[Tuple[int, int, int, int]] = []
        missing = 0

        for rec in included_records:
            smi = smiles_by_ccd.get(rec.ligand_ccd, "")
            if not smi:
                missing += 1
                continue

            substrate_idx = smiles_to_idx.get(smi)
            if substrate_idx is None:
                # This shouldn't happen if B3 is superset
                sys.stderr.write(f"[WARN] {scheme.name}: SMILES not in Substrates.csv: {rec.record_id}\n")
                continue

            enzyme_key = (rec.pdb_id, rec.ligand_ccd)
            enzyme_idx = enzyme_by_key.get(enzyme_key)
            if enzyme_idx is None:
                continue

            if rec.final_class in scheme.positive:
                label = 1
            elif rec.final_class in scheme.negative:
                label = 0
            else:
                continue

            dock_idx = len(data_rows)
            data_rows.append((dock_idx, enzyme_idx, substrate_idx, label))

        out_folder = SCHEME_FOLDER_NAMES.get(scheme.name, scheme.name)
        out_path = output_dir / out_folder / "data.csv"
        _write_csv(
            out_path,
            header=["Dock Index", "Enzyme Index", "Substrate Index", "Label"],
            rows=data_rows,
        )
        pos = sum(1 for *_, y in data_rows if y == 1)
        neg = sum(1 for *_, y in data_rows if y == 0)
        print(f"[INFO] Wrote {out_path} (n={len(data_rows)}, pos={pos}, neg={neg})")

    # Final validation summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    print(f"Canonical records: {len(canonical_records)}")
    print(f"Substrates.csv: {len(b3_substrates)} unique SMILES")

    # Count by class
    class_counts = {}
    for r in canonical_records:
        class_counts[r.final_class] = class_counts.get(r.final_class, 0) + 1
    for cls, cnt in sorted(class_counts.items()):
        print(f"  {cls}: {cnt}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
