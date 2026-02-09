#!/usr/bin/env python3
"""
Extract canonical SMILES for unique PDB CCD (Chemical Component Dictionary) codes.

Workflow (per CCD code):
  1) Primary: query RCSB CCD API: https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code}
  2) Fallback: if API fails or yields invalid SMILES, try extracting from local PDB files
  3) Validate + standardize with RDKit: MolFromSmiles → SanitizeMol → MolToSmiles (canonical)

Output CSV columns:
  - ccd_code, smiles, source (API/PDB_file/FAILED), validation_status (VALID/INVALID), error_message

Progress: print every 50 codes, save every 100 codes, retry/backoff for network errors
"""

import argparse
import csv
import json
import random
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import requests
except Exception:
    requests = None

try:
    from rdkit import Chem
except Exception as e:
    Chem = None
    _RDKIT_IMPORT_ERROR = e
else:
    _RDKIT_IMPORT_ERROR = None

RCSB_CCD_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}"


@dataclass
class SmilesResult:
    ccd_code: str
    smiles: str
    source: str  # API / PDB_file / FAILED
    validation_status: str  # VALID / INVALID
    error_message: str


def _jittered_backoff(base_seconds: float, attempt: int, max_seconds: float) -> float:
    """Exponential backoff with jitter (attempt is 1-based)"""
    delay = base_seconds * (2 ** (attempt - 1))
    delay = min(delay, max_seconds)
    delay += random.uniform(0.0, min(0.25, delay * 0.1))
    return delay


def _safe_str(x) -> str:
    return "" if x is None else str(x)


def _read_unique_ccd_codes_and_pdb_map(
    input_csv: Path, ccd_column: str, pdb_column: str, encoding: str = "utf-8-sig"
) -> Tuple[List[str], Dict[str, List[str]]]:
    """Returns: (sorted unique CCD codes, mapping: CCD -> sorted unique PDB IDs)"""
    unique_ccd = set()
    ccd_to_pdbs: Dict[str, set] = {}

    with input_csv.open("r", encoding=encoding, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {input_csv}")
        if ccd_column not in reader.fieldnames:
            raise ValueError(f"Missing column '{ccd_column}' in {input_csv}")
        if pdb_column not in reader.fieldnames:
            raise ValueError(f"Missing column '{pdb_column}' in {input_csv}")

        for row in reader:
            ccd = _safe_str(row.get(ccd_column)).strip().upper()
            pdb = _safe_str(row.get(pdb_column)).strip().upper()
            if not ccd:
                continue
            unique_ccd.add(ccd)
            if pdb:
                ccd_to_pdbs.setdefault(ccd, set()).add(pdb)

    codes = sorted(unique_ccd)
    pdb_map: Dict[str, List[str]] = {k: sorted(v) for k, v in ccd_to_pdbs.items()}
    return codes, pdb_map


def _load_existing_results(output_csv: Path) -> Dict[str, SmilesResult]:
    """Resume support: load existing results from output CSV if present"""
    if not output_csv.exists():
        return {}

    results: Dict[str, SmilesResult] = {}
    try:
        with output_csv.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                ccd = _safe_str(row.get("ccd_code")).strip().upper()
                if not ccd:
                    continue
                results[ccd] = SmilesResult(
                    ccd_code=ccd,
                    smiles=_safe_str(row.get("smiles")).strip(),
                    source=_safe_str(row.get("source")).strip() or "FAILED",
                    validation_status=_safe_str(row.get("validation_status")).strip() or "INVALID",
                    error_message=_safe_str(row.get("error_message")).strip(),
                )
    except Exception:
        return {}
    return results


def _write_results_csv_atomic(output_csv: Path, rows: Sequence[SmilesResult]) -> None:
    """Atomic write with temp file"""
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_csv.with_suffix(output_csv.suffix + ".tmp")

    fieldnames = ["ccd_code", "smiles", "source", "validation_status", "error_message"]
    with tmp_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(
                {
                    "ccd_code": r.ccd_code,
                    "smiles": r.smiles,
                    "source": r.source,
                    "validation_status": r.validation_status,
                    "error_message": r.error_message,
                }
            )

    tmp_path.replace(output_csv)


def _write_error_report_csv(output_csv: Path, rows: Sequence[SmilesResult]) -> Optional[Path]:
    """Generate error report CSV for failed/invalid records"""
    error_rows = [r for r in rows if r.error_message or r.validation_status != "VALID" or r.source == "FAILED"]
    if not error_rows:
        return None
    err_path = output_csv.with_name(output_csv.stem + "_errors.csv")
    _write_results_csv_atomic(err_path, error_rows)
    return err_path


def _rdkit_require() -> None:
    if Chem is None:
        raise RuntimeError(f"RDKit is required but could not be imported. Import error: {_RDKIT_IMPORT_ERROR!r}")


def _validate_and_canonicalize_smiles(smiles: str) -> Tuple[Optional[str], Optional[str]]:
    """Returns (canonical_smiles, error_message). If invalid -> (None, error_message)"""
    _rdkit_require()

    s = (smiles or "").strip()
    if not s:
        return None, "empty_smiles"

    mol = Chem.MolFromSmiles(s, sanitize=False)
    if mol is None:
        return None, "MolFromSmiles_returned_None"

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None, f"SanitizeMol_failed: {type(e).__name__}: {e}"

    try:
        can = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception as e:
        return None, f"MolToSmiles_failed: {type(e).__name__}: {e}"

    # Double-check canonical SMILES parseability
    mol2 = Chem.MolFromSmiles(can, sanitize=False)
    if mol2 is None:
        return None, "canonical_smiles_not_parseable"
    try:
        Chem.SanitizeMol(mol2)
    except Exception as e:
        return None, f"canonical_smiles_sanitize_failed: {type(e).__name__}: {e}"

    return can, None


def _http_get_json(
    url: str,
    timeout_seconds: float,
    max_retries: int,
    backoff_base_seconds: float,
    backoff_max_seconds: float,
    min_delay_seconds: float,
) -> Tuple[Optional[dict], Optional[str]]:
    """Returns (json_dict, error_message). Handles 429 rate limiting and retries"""
    time.sleep(min_delay_seconds)  # Pacing delay

    last_err = None
    for attempt in range(1, max_retries + 1):
        try:
            if requests is not None:
                resp = requests.get(url, timeout=timeout_seconds, headers={"Accept": "application/json"})
                status = resp.status_code
                if status == 200:
                    return resp.json(), None
                if status == 404:
                    return None, "http_404_not_found"
                if status == 429:
                    retry_after = resp.headers.get("Retry-After")
                    if retry_after and str(retry_after).strip().isdigit():
                        sleep_s = float(str(retry_after).strip())
                    else:
                        sleep_s = _jittered_backoff(backoff_base_seconds, attempt, backoff_max_seconds)
                    last_err = f"http_429_rate_limited (sleep {sleep_s:.2f}s)"
                    time.sleep(sleep_s)
                    continue
                if 500 <= status <= 599:
                    sleep_s = _jittered_backoff(backoff_base_seconds, attempt, backoff_max_seconds)
                    last_err = f"http_{status}_server_error (sleep {sleep_s:.2f}s)"
                    time.sleep(sleep_s)
                    continue
                return None, f"http_{status}"
            else:
                # urllib fallback
                import urllib.request
                import urllib.error

                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=timeout_seconds) as resp:
                    status = getattr(resp, "status", 200)
                    body = resp.read().decode("utf-8", errors="replace")
                    if status != 200:
                        return None, f"http_{status}"
                    return json.loads(body), None
        except Exception as e:
            sleep_s = _jittered_backoff(backoff_base_seconds, attempt, backoff_max_seconds)
            last_err = f"network_error: {type(e).__name__}: {e} (sleep {sleep_s:.2f}s)"
            time.sleep(sleep_s)
            continue

    return None, last_err or "unknown_network_error"


def _extract_smiles_candidates_from_ccd_json(data: dict) -> List[str]:
    """Extract likely SMILES fields from RCSB chemcomp JSON"""
    candidates: List[str] = []

    # rcsb_chem_comp_descriptor: smiles_stereo, smiles
    desc = data.get("rcsb_chem_comp_descriptor") or {}
    if isinstance(desc, dict):
        for k in ("smilesstereo", "smiles_stereo", "smiles", "SMILES_STEREO", "SMILES"):
            v = desc.get(k)
            if isinstance(v, str) and v.strip():
                candidates.append(v.strip())

    # pdbx_chem_comp_descriptor: list of {type, descriptor, ...}
    lst = data.get("pdbx_chem_comp_descriptor")
    if isinstance(lst, list):
        preferred_types = ["SMILES_STEREO", "SMILES_CANONICAL", "SMILES", "OpenEye OEToolkits SMILES"]
        typed: List[Tuple[int, str]] = []
        for item in lst:
            if not isinstance(item, dict):
                continue
            t = _safe_str(item.get("type")).strip()
            smi = _safe_str(item.get("descriptor")).strip()
            if not smi:
                continue
            rank = preferred_types.index(t) if t in preferred_types else len(preferred_types)
            typed.append((rank, smi))
        for _, smi in sorted(typed, key=lambda x: x[0]):
            candidates.append(smi)

    # De-dup while preserving order
    seen = set()
    out: List[str] = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            out.append(c)
    return out


def _try_api_for_ccd(
    ccd_code: str,
    timeout_seconds: float,
    max_retries: int,
    backoff_base_seconds: float,
    backoff_max_seconds: float,
    min_delay_seconds: float,
) -> Tuple[Optional[str], Optional[str]]:
    """Returns (canonical_smiles, error_message)"""
    url = RCSB_CCD_URL.format(ccd=ccd_code)
    data, err = _http_get_json(
        url=url,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        backoff_base_seconds=backoff_base_seconds,
        backoff_max_seconds=backoff_max_seconds,
        min_delay_seconds=min_delay_seconds,
    )
    if data is None:
        return None, f"api_failed: {err}"

    candidates = _extract_smiles_candidates_from_ccd_json(data)
    if not candidates:
        return None, "api_no_smiles_fields"

    # Validate candidates in priority order
    last_err = None
    for smi in candidates:
        can, v_err = _validate_and_canonicalize_smiles(smi)
        if can is not None:
            return can, None
        last_err = v_err

    return None, f"api_smiles_invalid: {last_err or 'unknown'}"


def _parse_pdb_residue_key(line: str) -> Tuple[str, str, str]:
    """Extract chainID, resSeq, iCode from PDB ATOM/HETATM line"""
    chain = (line[21:22] if len(line) >= 22 else "").strip()
    resseq = (line[22:26] if len(line) >= 26 else "").strip()
    icode = (line[26:27] if len(line) >= 27 else "").strip()
    return chain, resseq, icode


def _parse_pdb_atom_serial(line: str) -> Optional[int]:
    """Extract atom serial number from PDB ATOM/HETATM line"""
    try:
        return int(line[6:11].strip())
    except Exception:
        return None


_CONECT_INT_RE = re.compile(r"\d+")


def _filter_conect_line_to_serials(line: str, allowed_serials: set) -> Optional[str]:
    """Keep CONECT records only between atoms within the selected residue"""
    nums = [int(x) for x in _CONECT_INT_RE.findall(line[6:])]
    if not nums:
        return None
    if nums[0] not in allowed_serials:
        return None
    kept = [nums[0]] + [n for n in nums[1:] if n in allowed_serials]
    if len(kept) < 2:
        return None
    return "CONECT" + "".join(f"{n:5d}" for n in kept)


def _extract_smiles_from_pdb_file(pdb_path: Path, ccd_code: str) -> Tuple[Optional[str], Optional[str]]:
    """Try to build ligand molecule from PDB file by selecting largest residue instance with resname == ccd_code"""
    _rdkit_require()

    try:
        text = pdb_path.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return None, f"pdb_read_failed: {type(e).__name__}: {e}"

    residues: Dict[Tuple[str, str, str], List[str]] = {}
    residue_serials: Dict[Tuple[str, str, str], set] = {}
    conect_lines: List[str] = []

    for raw in text.splitlines():
        rec = raw[:6].strip()
        if rec == "HETATM":
            resname = (raw[17:20] if len(raw) >= 20 else "").strip().upper()
            if resname != ccd_code:
                continue
            key = _parse_pdb_residue_key(raw)
            residues.setdefault(key, []).append(raw.rstrip("\n"))
            serial = _parse_pdb_atom_serial(raw)
            if serial is not None:
                residue_serials.setdefault(key, set()).add(serial)
        elif rec == "CONECT":
            conect_lines.append(raw.rstrip("\n"))

    if not residues:
        return None, "pdb_no_matching_HETATM_residue"

    # Choose largest instance (most atoms)
    best_key = max(residues.keys(), key=lambda k: len(residues[k]))
    atom_lines = residues[best_key]
    allowed_serials = residue_serials.get(best_key, set())

    filtered_conect: List[str] = []
    if allowed_serials and conect_lines:
        for cl in conect_lines:
            kept = _filter_conect_line_to_serials(cl, allowed_serials)
            if kept is not None:
                filtered_conect.append(kept)

    pdb_block = "\n".join(atom_lines + filtered_conect + ["END"]) + "\n"

    try:
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    except Exception as e:
        return None, f"MolFromPDBBlock_failed: {type(e).__name__}: {e}"
    if mol is None:
        return None, "MolFromPDBBlock_returned_None"

    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return None, f"SanitizeMol_failed_on_PDB_mol: {type(e).__name__}: {e}"

    try:
        smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception as e:
        return None, f"MolToSmiles_failed_on_PDB_mol: {type(e).__name__}: {e}"

    can, v_err = _validate_and_canonicalize_smiles(smi)
    if can is None:
        return None, f"pdb_smiles_invalid_after_canonicalize: {v_err}"
    return can, None


def _try_pdb_fallback_for_ccd(
    ccd_code: str, pdb_ids: Sequence[str], pdb_dir: Path
) -> Tuple[Optional[str], Optional[str]]:
    """Iterate through associated PDB IDs and attempt to extract ligand SMILES from local PDB files"""
    if not pdb_ids:
        return None, "no_pdb_ids_for_ccd"
    if not pdb_dir.exists():
        return None, f"pdb_dir_missing: {pdb_dir}"

    last_err = None
    for pdb_id in pdb_ids:
        pdb_path = pdb_dir / f"{pdb_id.upper()}.pdb"
        if not pdb_path.exists():
            last_err = f"missing_pdb_file_for_entry: {pdb_id}"
            continue
        smi, err = _extract_smiles_from_pdb_file(pdb_path, ccd_code)
        if smi is not None:
            return smi, None
        last_err = f"{pdb_id}: {err}"

    return None, f"pdb_fallback_failed: {last_err or 'unknown'}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description="Extract canonical SMILES from PDB CCD codes with API+PDB fallback.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("source_data/01_核心数据/修复后最终版/独立测试集_682条.csv"),
        help="Input CSV path (must contain ligand_ccd and pdb_id columns).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/03_Step3_配体处理与数据集构建/05_提取SMILES/ccd_smiles.csv"),
        help="Output CSV path.",
    )
    parser.add_argument("--ccd-column", type=str, default="ligand_ccd", help="Column name for CCD codes.")
    parser.add_argument("--pdb-column", type=str, default="pdb_id", help="Column name for PDB IDs.")
    parser.add_argument(
        "--pdb-dir",
        type=Path,
        default=Path("data/01_Step1_PDB文件"),
        help="Directory containing downloaded PDB structure files (e.g., 1ABC.pdb).",
    )
    parser.add_argument("--no-api", action="store_true", help="Skip API and only try PDB-file fallback.")
    parser.add_argument("--resume", action="store_true", default=True, help="Resume from existing output CSV if present.")
    parser.add_argument("--timeout", type=float, default=30.0, help="HTTP timeout seconds.")
    parser.add_argument("--max-retries", type=int, default=6, help="Max retries for network/API errors.")
    parser.add_argument("--backoff-base", type=float, default=1.0, help="Exponential backoff base seconds.")
    parser.add_argument("--backoff-max", type=float, default=30.0, help="Exponential backoff max seconds.")
    parser.add_argument("--api-min-delay", type=float, default=0.15, help="Minimum delay between API calls (seconds).")
    args = parser.parse_args(argv)

    if Chem is None:
        print(
            "ERROR: RDKit import failed; required for SMILES validation/standardization.\n"
            f"Import error: {_RDKIT_IMPORT_ERROR!r}",
            file=sys.stderr,
        )
        return 2

    if not args.input.exists():
        print(f"ERROR: input CSV not found: {args.input}", file=sys.stderr)
        return 2

    codes, ccd_to_pdbs = _read_unique_ccd_codes_and_pdb_map(
        input_csv=args.input, ccd_column=args.ccd_column, pdb_column=args.pdb_column
    )
    total = len(codes)
    print(f"Unique CCD codes: {total}")

    existing: Dict[str, SmilesResult] = _load_existing_results(args.output) if args.resume else {}
    if existing:
        print(f"Resume enabled: loaded {len(existing)} existing results from {args.output}")

    results: Dict[str, SmilesResult] = dict(existing)
    processed_since_save = 0

    for idx, ccd in enumerate(codes, start=1):
        if ccd in results:
            continue

        source = "FAILED"
        validation_status = "INVALID"
        smiles_out = ""
        err_msg = ""

        try:
            if not args.no_api:
                smi, api_err = _try_api_for_ccd(
                    ccd_code=ccd,
                    timeout_seconds=args.timeout,
                    max_retries=args.max_retries,
                    backoff_base_seconds=args.backoff_base,
                    backoff_max_seconds=args.backoff_max,
                    min_delay_seconds=args.api_min_delay,
                )
                if smi is not None:
                    source = "API"
                    validation_status = "VALID"
                    smiles_out = smi
                    err_msg = ""
                else:
                    err_msg = api_err or "api_failed_unknown"

            # Fallback if API failed or SMILES invalid
            if source != "API":
                pdb_ids = ccd_to_pdbs.get(ccd, [])
                smi, pdb_err = _try_pdb_fallback_for_ccd(ccd_code=ccd, pdb_ids=pdb_ids, pdb_dir=args.pdb_dir)
                if smi is not None:
                    source = "PDB_file"
                    validation_status = "VALID"
                    smiles_out = smi
                    err_msg = ""
                else:
                    if err_msg:
                        err_msg = f"{err_msg}; {pdb_err}"
                    else:
                        err_msg = pdb_err or "pdb_fallback_failed_unknown"

        except Exception as e:
            source = "FAILED"
            validation_status = "INVALID"
            smiles_out = ""
            err_msg = f"unhandled_exception: {type(e).__name__}: {e}"

        results[ccd] = SmilesResult(
            ccd_code=ccd,
            smiles=smiles_out,
            source=source,
            validation_status=validation_status,
            error_message=err_msg,
        )
        processed_since_save += 1

        if idx % 50 == 0 or idx == total:
            n_valid = sum(1 for r in results.values() if r.validation_status == "VALID")
            n_failed = sum(1 for r in results.values() if r.source == "FAILED")
            print(f"Progress: {idx}/{total} | VALID={n_valid} | FAILED={n_failed}")

        if processed_since_save >= 100:
            try:
                rows = [results[k] for k in sorted(results.keys())]
                _write_results_csv_atomic(args.output, rows)
                print(f"Saved intermediate results: {args.output} (rows={len(rows)})")
            except Exception as e:
                print(f"WARNING: failed to save intermediate CSV: {type(e).__name__}: {e}", file=sys.stderr)
            processed_since_save = 0

    # Final write
    rows = [results[k] for k in sorted(results.keys())]
    _write_results_csv_atomic(args.output, rows)
    err_path = _write_error_report_csv(args.output, rows)

    n_valid = sum(1 for r in rows if r.validation_status == "VALID")
    n_invalid = len(rows) - n_valid
    n_api = sum(1 for r in rows if r.source == "API")
    n_pdb = sum(1 for r in rows if r.source == "PDB_file")
    n_failed = sum(1 for r in rows if r.source == "FAILED")

    print("\nDone.")
    print(f"Output CSV: {args.output}")
    if err_path is not None:
        print(f"Error report CSV: {err_path}")
    print(f"Stats: total={len(rows)} valid={n_valid} invalid={n_invalid} api={n_api} pdb_fallback={n_pdb} failed={n_failed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
