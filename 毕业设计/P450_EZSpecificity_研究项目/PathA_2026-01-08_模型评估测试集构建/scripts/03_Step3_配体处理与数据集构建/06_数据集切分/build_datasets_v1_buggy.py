# -*- coding: utf-8 -*-
"""
Build B1/B2/B3 datasets from P450 enzyme-ligand classification results.

Inputs (defaults are project-relative):
  1) source_data/01_核心数据/修复后最终版/独立测试集_682条.csv
  2) data/03_Step3_配体处理与数据集构建/04_配体三方深度验证/03_合并结果/三方深度验证结果_merged.csv
  3) data/03_Step3_配体处理与数据集构建/05_提取SMILES/ccd_smiles.csv
  4) data/02_Step2_酶序列/Enzymes.csv
  5) data/02_Step2_酶序列/record_enzyme_mapping.csv

Outputs (to data/03_Step3_配体处理与数据集构建/06_数据集切分/ by default):
  - Substrates.csv: Substrate_Index, Substrate_SMILES (unique by SMILES; built from B3-eligible records)
  - B1/data.csv, B2/data.csv, B3/data.csv: Dock Index, Enzyme Index, Substrate Index, Label

Schemes:
  - B1 includes: SUBSTRATE                 ; Label=1: SUBSTRATE                   ; Label=0: (none)
  - B2 includes: SUBSTRATE, PRODUCT        ; Label=1: SUBSTRATE, PRODUCT          ; Label=0: (none)
  - B3 includes: SUBSTRATE, PRODUCT, INHIBITOR ; Label=1: SUBSTRATE, PRODUCT      ; Label=0: INHIBITOR

Notes:
  - EXCLUDE is never included in any scheme.
  - record_id REC_0001 corresponds to row 0 in 独立测试集_682条.csv.
  - Use utf-8-sig for I/O (BOM-safe for Chinese filenames).
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


UTF8_SIG = "utf-8-sig"


def _norm_pdb_id(v: str) -> str:
    return (v or "").strip().lower()


def _norm_ccd(v: str) -> str:
    return (v or "").strip().upper()


def _norm_class(v: str) -> str:
    return (v or "").strip().upper()


def _read_csv_dicts(path: Path, *, encoding: str = UTF8_SIG) -> List[Dict[str, str]]:
    with path.open("r", encoding=encoding, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return [dict(row) for row in reader]


def _write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]], *, encoding: str = UTF8_SIG) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding=encoding, newline="") as f:
        w = csv.writer(f)
        w.writerow(list(header))
        for r in rows:
            w.writerow(list(r))


def _parse_record_id_to_index(record_id: str) -> int:
    """REC_0001 -> 0, REC_0002 -> 1, ..."""
    s = (record_id or "").strip()
    if not s.startswith("REC_"):
        raise ValueError(f"Unexpected record_id format: {record_id!r}")
    num_part = s.split("_", 1)[1]
    try:
        n = int(num_part)
    except ValueError as e:
        raise ValueError(f"Unexpected record_id format: {record_id!r}") from e
    if n <= 0:
        raise ValueError(f"record_id must be 1-indexed: {record_id!r}")
    return n - 1


@dataclass(frozen=True)
class Scheme:
    name: str
    included: Set[str]
    positive: Set[str]
    negative: Set[str]


SCHEMES: List[Scheme] = [
    Scheme(name="B1", included={"SUBSTRATE"}, positive={"SUBSTRATE"}, negative=set()),
    Scheme(name="B2", included={"SUBSTRATE", "PRODUCT"}, positive={"SUBSTRATE", "PRODUCT"}, negative=set()),
    Scheme(name="B3", included={"SUBSTRATE", "PRODUCT", "INHIBITOR"}, positive={"SUBSTRATE", "PRODUCT"}, negative={"INHIBITOR"}),
]


def _default_paths(base_dir: Path) -> Dict[str, Path]:
    return {
        "independent_csv": base_dir / "source_data" / "01_核心数据" / "修复后最终版" / "独立测试集_682条.csv",
        "merged_results_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "04_配体三方深度验证" / "03_合并结果" / "三方深度验证结果_merged.csv",
        "ccd_smiles_csv": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "05_提取SMILES" / "ccd_smiles.csv",
        "enzymes_csv": base_dir / "data" / "02_Step2_酶序列" / "Enzymes.csv",
        "record_enzyme_mapping_csv": base_dir / "data" / "02_Step2_酶序列" / "record_enzyme_mapping.csv",
        "output_dir": base_dir / "data" / "03_Step3_配体处理与数据集构建" / "06_数据集切分",
    }


def _load_smiles_map(ccd_smiles_csv: Path) -> Dict[str, str]:
    rows = _read_csv_dicts(ccd_smiles_csv)
    out: Dict[str, str] = {}
    for r in rows:
        ccd = _norm_ccd(r.get("ccd_code", ""))
        smi = (r.get("smiles", "") or "").strip()
        if not ccd or not smi:
            continue
        out.setdefault(ccd, smi)
    return out


def _load_enzyme_index_set(enzymes_csv: Path) -> Set[int]:
    rows = _read_csv_dicts(enzymes_csv)
    out: Set[int] = set()
    for r in rows:
        s = (r.get("Enzyme_Index", "") or "").strip()
        if not s:
            continue
        try:
            out.add(int(s))
        except ValueError:
            continue
    return out


def _load_record_to_enzyme_index_map(mapping_csv: Path) -> Dict[Tuple[str, str], int]:
    """Keyed by (pdb_id, ligand_ccd)"""
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
        sys.stderr.write("[WARN] Conflicting Enzyme_Index mappings found for some (pdb_id, ligand_ccd) keys.\n")
        for (k, a, b) in conflicts[:20]:
            sys.stderr.write(f"  [WARN] {k} -> {a} vs {b} (keeping {a})\n")
        if len(conflicts) > 20:
            sys.stderr.write(f"  [WARN] ... and {len(conflicts) - 20} more conflicts\n")
    return out


def _build_substrates_from_b3(
    independent_rows: Sequence[Dict[str, str]],
    merged_rows: Sequence[Dict[str, str]],
    smiles_by_ccd: Dict[str, str],
    *,
    strict: bool,
) -> Tuple[List[str], Dict[str, int]]:
    """Returns: (substrates_smiles, substrate_index_by_smiles)"""
    b3_included = next(s for s in SCHEMES if s.name == "B3").included
    seen: Dict[str, int] = {}
    substrates: List[str] = []
    missing_smiles_ccd: List[Tuple[str, str]] = []

    for mr in merged_rows:
        final_class = _norm_class(mr.get("final_class", ""))
        if final_class == "EXCLUDE" or final_class not in b3_included:
            continue
        record_id = (mr.get("record_id", "") or "").strip()
        idx = _parse_record_id_to_index(record_id) if record_id else None
        if idx is None or idx < 0 or idx >= len(independent_rows):
            raise IndexError(f"record_id {record_id!r} points outside 独立测试集 rows (n={len(independent_rows)})")
        ligand_ccd = _norm_ccd(independent_rows[idx].get("ligand_ccd", ""))
        if not ligand_ccd:
            missing_smiles_ccd.append((record_id, ligand_ccd))
            continue
        smi = smiles_by_ccd.get(ligand_ccd, "")
        if not smi:
            missing_smiles_ccd.append((record_id, ligand_ccd))
            continue
        if smi in seen:
            continue
        seen[smi] = len(substrates)
        substrates.append(smi)

    if missing_smiles_ccd:
        msg = f"Missing SMILES for {len(missing_smiles_ccd)} B3-eligible records. Example: {missing_smiles_ccd[:5]!r}\n"
        if strict:
            raise ValueError(msg.strip())
        sys.stderr.write("[WARN] " + msg)

    return substrates, seen


def _validate_record_order(merged_rows: Sequence[Dict[str, str]], n_independent: int) -> None:
    if len(merged_rows) != n_independent:
        raise ValueError(f"Row count mismatch: merged_results={len(merged_rows)} vs 独立测试集={n_independent}")
    bad = 0
    for i, mr in enumerate(merged_rows):
        rid = (mr.get("record_id", "") or "").strip()
        expected = f"REC_{i+1:04d}"
        if rid != expected:
            bad += 1
            if bad <= 10:
                sys.stderr.write(f"[WARN] record_id mismatch at merged row {i}: got {rid!r}, expected {expected!r}\n")
    if bad:
        sys.stderr.write(f"[WARN] record_id mismatches: {bad} (will join using record_id -> row index mapping)\n")


def _build_scheme_rows(
    scheme: Scheme,
    independent_rows: Sequence[Dict[str, str]],
    merged_rows: Sequence[Dict[str, str]],
    smiles_by_ccd: Dict[str, str],
    enzyme_index_by_key: Dict[Tuple[str, str], int],
    substrate_index_by_smiles: Dict[str, int],
    valid_enzyme_indices: Optional[Set[int]],
    *,
    strict: bool,
) -> List[Tuple[int, int, int, int]]:
    """Returns rows: (Dock Index, Enzyme Index, Substrate Index, Label)"""
    rows_out: List[Tuple[int, int, int, int]] = []
    missing_smiles: List[Tuple[str, str]] = []
    missing_enzyme: List[Tuple[str, str, str]] = []
    invalid_enzyme: List[Tuple[str, int]] = []

    for mr in merged_rows:
        record_id = (mr.get("record_id", "") or "").strip()
        final_class = _norm_class(mr.get("final_class", ""))
        if final_class == "EXCLUDE" or final_class not in scheme.included:
            continue

        idx = _parse_record_id_to_index(record_id)
        if idx < 0 or idx >= len(independent_rows):
            raise IndexError(f"record_id {record_id!r} points outside 独立测试集 rows (n={len(independent_rows)})")
        ir = independent_rows[idx]

        pdb_id = _norm_pdb_id(ir.get("pdb_id", ""))
        ligand_ccd = _norm_ccd(ir.get("ligand_ccd", ""))

        smi = smiles_by_ccd.get(ligand_ccd, "")
        if not smi:
            missing_smiles.append((record_id, ligand_ccd))
            continue
        substrate_idx = substrate_index_by_smiles.get(smi)
        if substrate_idx is None:
            raise KeyError(f"SMILES not found in Substrates.csv index (record_id={record_id!r}, ccd={ligand_ccd!r})")

        enzyme_key = (pdb_id, ligand_ccd)
        enzyme_idx = enzyme_index_by_key.get(enzyme_key)
        if enzyme_idx is None:
            missing_enzyme.append((record_id, pdb_id, ligand_ccd))
            continue
        if valid_enzyme_indices is not None and enzyme_idx not in valid_enzyme_indices:
            invalid_enzyme.append((record_id, enzyme_idx))

        if final_class in scheme.positive:
            label = 1
        elif final_class in scheme.negative:
            label = 0
        else:
            raise ValueError(f"Class {final_class!r} is included but not labeled in scheme {scheme.name}")

        dock_idx = len(rows_out)
        rows_out.append((dock_idx, enzyme_idx, substrate_idx, label))

    if missing_smiles:
        msg = f"{scheme.name}: missing SMILES for {len(missing_smiles)} included records. Example: {missing_smiles[:5]!r}\n"
        if strict:
            raise ValueError(msg.strip())
        sys.stderr.write("[WARN] " + msg)
    if missing_enzyme:
        msg = f"{scheme.name}: missing Enzyme_Index for {len(missing_enzyme)} included records. Example: {missing_enzyme[:5]!r}\n"
        if strict:
            raise ValueError(msg.strip())
        sys.stderr.write("[WARN] " + msg)
    if invalid_enzyme:
        sys.stderr.write(f"[WARN] {scheme.name}: {len(invalid_enzyme)} rows have Enzyme_Index not in Enzymes.csv. Example: {invalid_enzyme[:5]!r}\n")

    return rows_out


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Build B1/B2/B3 datasets (Substrates + scheme data.csv files).")
    p.add_argument("--base-dir", default=".", help="Base directory for resolving default input/output paths.")
    p.add_argument("--independent-csv", default=None, help="Path to 独立测试集_682条.csv")
    p.add_argument("--merged-results-csv", default=None, help="Path to 三方深度验证结果_merged.csv")
    p.add_argument("--ccd-smiles-csv", default=None, help="Path to ccd_smiles.csv")
    p.add_argument("--enzymes-csv", default=None, help="Path to Enzymes.csv")
    p.add_argument("--record-enzyme-mapping-csv", default=None, help="Path to record_enzyme_mapping.csv")
    p.add_argument("--output-dir", default=None, help="Output directory for Substrates.csv and per-scheme folders.")
    p.add_argument("--strict", action="store_true", help="Fail if any included record is missing SMILES or Enzyme_Index.")
    args = p.parse_args(argv)

    base_dir = Path(args.base_dir).resolve()
    defaults = _default_paths(base_dir)

    independent_csv = Path(args.independent_csv) if args.independent_csv else defaults["independent_csv"]
    merged_results_csv = Path(args.merged_results_csv) if args.merged_results_csv else defaults["merged_results_csv"]
    ccd_smiles_csv = Path(args.ccd_smiles_csv) if args.ccd_smiles_csv else defaults["ccd_smiles_csv"]
    enzymes_csv = Path(args.enzymes_csv) if args.enzymes_csv else defaults["enzymes_csv"]
    record_enzyme_mapping_csv = Path(args.record_enzyme_mapping_csv) if args.record_enzyme_mapping_csv else defaults["record_enzyme_mapping_csv"]
    output_dir = Path(args.output_dir) if args.output_dir else defaults["output_dir"]

    strict = bool(args.strict)

    print(f"[INFO] base_dir={base_dir}")
    print(f"[INFO] independent_csv={independent_csv}")
    print(f"[INFO] merged_results_csv={merged_results_csv}")
    print(f"[INFO] ccd_smiles_csv={ccd_smiles_csv}")
    print(f"[INFO] enzymes_csv={enzymes_csv}")
    print(f"[INFO] record_enzyme_mapping_csv={record_enzyme_mapping_csv}")
    print(f"[INFO] output_dir={output_dir}")
    print(f"[INFO] strict={strict}")

    independent_rows = _read_csv_dicts(independent_csv)
    merged_rows = _read_csv_dicts(merged_results_csv)

    _validate_record_order(merged_rows, n_independent=len(independent_rows))

    smiles_by_ccd = _load_smiles_map(ccd_smiles_csv)
    enzyme_index_by_key = _load_record_to_enzyme_index_map(record_enzyme_mapping_csv)
    valid_enzyme_indices = _load_enzyme_index_set(enzymes_csv)

    substrates_smiles, substrate_index_by_smiles = _build_substrates_from_b3(
        independent_rows, merged_rows, smiles_by_ccd, strict=strict
    )

    # Write Substrates.csv
    substrates_path = output_dir / "Substrates.csv"
    _write_csv(
        substrates_path,
        header=["Substrate_Index", "Substrate_SMILES"],
        rows=((i, smi) for i, smi in enumerate(substrates_smiles)),
        encoding=UTF8_SIG,
    )
    print(f"[INFO] Wrote {substrates_path} (n={len(substrates_smiles)})")

    # Write per-scheme data.csv files
    for scheme in SCHEMES:
        scheme_rows = _build_scheme_rows(
            scheme=scheme,
            independent_rows=independent_rows,
            merged_rows=merged_rows,
            smiles_by_ccd=smiles_by_ccd,
            enzyme_index_by_key=enzyme_index_by_key,
            substrate_index_by_smiles=substrate_index_by_smiles,
            valid_enzyme_indices=valid_enzyme_indices,
            strict=strict,
        )
        out_path = output_dir / scheme.name / "data.csv"
        _write_csv(
            out_path,
            header=["Dock Index", "Enzyme Index", "Substrate Index", "Label"],
            rows=scheme_rows,
            encoding=UTF8_SIG,
        )
        pos = sum(1 for *_, y in scheme_rows if y == 1)
        neg = sum(1 for *_, y in scheme_rows if y == 0)
        print(f"[INFO] Wrote {out_path} (n={len(scheme_rows)}, pos={pos}, neg={neg})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
