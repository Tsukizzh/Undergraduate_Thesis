"""
Step 2: Extract enzyme amino-acid sequences from PDB/mmCIF structure files.

Input:  data/00_raw/pdb_files/*.pdb + *.cif
Output: data/01_processed/Enzymes.csv (single column: "Protein sequence")

Constraints:
  - max sequence length: 1000 aa (reject if longer)
  - only standard 20 amino acids (reject if contains non-standard)
"""

from __future__ import annotations
import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

STANDARD20 = set("ACDEFGHIKLMNPQRSTVWY")

# Modified residues mapping to standard codes
CUSTOM_3TO1 = {
    "MSE": "M",  # selenomethionine
    "SEP": "S",  # phosphoserine
    "TPO": "T",  # phosphothreonine
    "PTR": "Y",  # phosphotyrosine
    "CSO": "C",  # cysteine sulfinic acid
    "CSD": "C",
    "CME": "C",
    "HYP": "P",  # hydroxyproline
    "SEC": "U",  # selenocysteine (will be rejected)
    "PYL": "O",  # pyrrolysine (will be rejected)
}


@dataclass(frozen=True)
class FileResult:
    file_name: str
    ok: bool
    selected_chain: str
    sequence: str
    length: int
    reason: str
    invalid_chars: str
    parser: str


def iter_structure_files(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    files = [p for p in input_dir.iterdir()
             if p.is_file() and p.suffix.lower() in {".pdb", ".cif", ".mmcif"}]
    files.sort(key=lambda p: p.name.lower())
    return files


def parse_structure(file_path: Path):
    suffix = file_path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
        parser_name = "mmcif"
    else:
        parser = PDBParser(QUIET=True)
        parser_name = "pdb"
    structure = parser.get_structure(file_path.stem, str(file_path))
    return structure, parser_name


def chain_to_sequence(chain) -> str:
    letters = []
    for res in chain.get_residues():
        if not is_aa(res, standard=False):
            continue
        letter = seq1(res.get_resname(), custom_map=CUSTOM_3TO1, undef_code="X")
        letters.append(letter)
    return "".join(letters).upper()


def validate_sequence(seq: str, max_len: int) -> tuple[bool, str, str]:
    if not seq:
        return False, "empty_sequence", ""
    if len(seq) > max_len:
        return False, f"too_long>{max_len}", ""
    invalid = sorted({c for c in seq if c not in STANDARD20})
    if invalid:
        return False, f"non_standard:{''.join(invalid)}", "".join(invalid)
    return True, "ok", ""


def extract_best_sequence(file_path: Path, max_len: int) -> FileResult:
    try:
        structure, parser_name = parse_structure(file_path)
    except Exception as e:
        return FileResult(file_path.name, False, "", "", 0,
                         f"parse_error:{type(e).__name__}", "", "unknown")

    try:
        model = next(structure.get_models())
    except StopIteration:
        return FileResult(file_path.name, False, "", "", 0,
                         "no_models", "", parser_name)

    candidates = []
    for chain in model.get_chains():
        seq = chain_to_sequence(chain)
        if not seq:
            continue
        ok, reason, invalid_chars = validate_sequence(seq, max_len)
        candidates.append((ok, len(seq), str(chain.id), seq, reason, invalid_chars))

    if not candidates:
        return FileResult(file_path.name, False, "", "", 0,
                         "no_protein_residues", "", parser_name)

    # Sort: valid first, then longest, then chain id
    candidates.sort(key=lambda t: (not t[0], -t[1], t[2]))
    ok, length, chain_id, seq, reason, invalid_chars = candidates[0]

    return FileResult(file_path.name, ok, chain_id, seq, length,
                     reason, invalid_chars, parser_name)


def write_enzymes_csv(output_csv: Path, sequences: Iterable[str]) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Protein sequence"])
        for seq in sequences:
            writer.writerow([seq])


def write_report_csv(report_csv: Path, results: Iterable[FileResult]) -> None:
    report_csv.parent.mkdir(parents=True, exist_ok=True)
    with report_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["file_name", "ok", "selected_chain", "length",
                        "reason", "invalid_chars", "parser"])
        for r in results:
            writer.writerow([r.file_name, int(r.ok), r.selected_chain,
                           r.length, r.reason, r.invalid_chars, r.parser])


def main() -> int:
    parser = argparse.ArgumentParser(description="Step 2: extract enzyme sequences")
    parser.add_argument("--input-dir", type=Path, default=None)
    parser.add_argument("--output-csv", type=Path, default=None)
    parser.add_argument("--report-csv", type=Path, default=None)
    parser.add_argument("--max-len", type=int, default=1000)
    parser.add_argument("--expect-files", type=int, default=None)
    parser.add_argument("--allow-skips", action="store_true")
    args = parser.parse_args()

    base_dir = Path(__file__).resolve().parent.parent
    input_dir = args.input_dir or (base_dir / "data" / "00_raw" / "pdb_files")
    output_csv = args.output_csv or (base_dir / "data" / "01_processed" / "Enzymes.csv")
    report_csv = args.report_csv or (base_dir / "reports" / "tables" / "step2_sequence_extraction_report.csv")

    files = iter_structure_files(input_dir)
    if args.expect_files is not None and len(files) != args.expect_files:
        print(f"ERROR: expected {args.expect_files} files, found {len(files)}", file=sys.stderr)
        return 2

    results = []
    sequences = []
    for fp in files:
        r = extract_best_sequence(fp, max_len=args.max_len)
        results.append(r)
        if r.ok:
            sequences.append(r.sequence)

    write_report_csv(report_csv, results)

    failed = [r for r in results if not r.ok]
    if failed and not args.allow_skips:
        print(f"ERROR: {len(failed)}/{len(results)} files failed; report: {report_csv}", file=sys.stderr)
        for r in failed[:20]:
            print(f"  - {r.file_name}: {r.reason}", file=sys.stderr)
        return 2

    write_enzymes_csv(output_csv, sequences)
    print(f"Input files: {len(files)} ({input_dir})")
    print(f"Valid sequences: {len(sequences)} ({output_csv})")
    print(f"Report: {report_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
