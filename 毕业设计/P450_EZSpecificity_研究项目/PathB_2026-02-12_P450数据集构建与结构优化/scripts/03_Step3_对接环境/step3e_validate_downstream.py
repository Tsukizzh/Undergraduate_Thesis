"""Step3E: Validate pilot output compatibility with Step 8-10 pipeline."""
from __future__ import annotations

import argparse
import csv
import logging
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

LOG = logging.getLogger(__name__)


def _norm_key(k: str) -> str:
    return k.strip().lower().replace(" ", "_")


def run_cmd(cmd: list[str], timeout: int = 600) -> tuple[int, str, str, float]:
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return proc.returncode, proc.stdout or "", proc.stderr or "", time.perf_counter() - t0


def build_mapping(pair_manifest: Path, output_csv: Path) -> int:
    """Create Dock_Index → Substrate_Index mapping for step8_align_ligand.py."""
    count = 0
    with pair_manifest.open("r", encoding="utf-8-sig", newline="") as fin, \
         output_csv.open("w", encoding="utf-8", newline="") as fout:
        w = csv.writer(fout)
        w.writerow(["Dock_Index", "Substrate_Index"])
        for row in csv.DictReader(fin):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            w.writerow([nr["dock_index"], nr["substrate_index"]])
            count += 1
    return count


def count_success(csv_path: Path) -> tuple[int, int]:
    if not csv_path.exists():
        return 0, 0
    total = ok = 0
    with csv_path.open("r", encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            nr = {_norm_key(k): v.strip() for k, v in row.items()}
            total += 1
            if nr.get("success", "").lower() in {"true", "1", "yes"}:
                ok += 1
    return total, ok


def tail(text: str, n: int = 30) -> str:
    lines = [l for l in text.splitlines() if l.strip()]
    return "\n".join(lines[-n:]) if len(lines) > n else "\n".join(lines)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot_dir", type=Path, required=True, help="Pilot 50 directory")
    p.add_argument("--step2_scripts_dir", type=Path, required=True, help="Step 2 scripts directory")
    p.add_argument("--substrates_csv", type=Path, required=True, help="B6 Substrates.csv")
    p.add_argument("--results_dir", type=Path, required=True, help="Report output directory")
    p.add_argument("--python_exe", type=Path, default=Path(sys.executable))
    return p.parse_args()


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    raw_ligand_dir = args.pilot_dir / "raw_ligand"
    pocket_dir = args.pilot_dir / "pocket"
    align_script = args.step2_scripts_dir / "step8_align_ligand.py"
    lmdb_script = args.step2_scripts_dir / "step8_generate_structure_lmdb.py"

    pair_manifest = args.pilot_dir / "pilot_50_pairs.csv"
    if not pair_manifest.exists():
        pair_manifest = args.pilot_dir / "pilot_5_pairs.csv"

    for label, path in [("raw_ligand_dir", raw_ligand_dir), ("pocket_dir", pocket_dir),
                         ("align_script", align_script), ("lmdb_script", lmdb_script),
                         ("pair_manifest", pair_manifest), ("substrates_csv", args.substrates_csv)]:
        if not path.exists():
            LOG.error("%s not found: %s", label, path)
            return 2

    args.results_dir.mkdir(parents=True, exist_ok=True)

    # Build mapping CSV
    mapping_csv = args.pilot_dir / "mapping_for_align.csv"
    n_mapped = build_mapping(pair_manifest, mapping_csv)
    LOG.info("Mapping CSV: %d rows", n_mapped)

    aligned_dir = args.pilot_dir / "ligand_aligned"
    align_summary = args.pilot_dir / "alignment_summary.csv"
    lmdb_dir = args.pilot_dir / "structure_features_validation"

    # Run step8_align_ligand.py
    cmd1 = [str(args.python_exe), str(align_script),
            "--raw_ligand_dir", str(raw_ligand_dir),
            "--aligned_ligand_dir", str(aligned_dir),
            "--mapping_csv", str(mapping_csv),
            "--substrates_csv", str(args.substrates_csv),
            "--summary_csv", str(align_summary)]
    rc1, out1, err1, dt1 = run_cmd(cmd1)
    LOG.info("step8_align_ligand.py: rc=%d, %.1fs", rc1, dt1)

    # Run step8_generate_structure_lmdb.py
    cmd2 = [str(args.python_exe), str(lmdb_script),
            "--pocket_dir", str(pocket_dir),
            "--ligand_dir", str(aligned_dir),
            "--alignment_summary", str(align_summary),
            "--output_dir", str(lmdb_dir),
            "--experiment_name", "STEP3_VALIDATION"]
    rc2, out2, err2, dt2 = run_cmd(cmd2)
    LOG.info("step8_generate_structure_lmdb.py: rc=%d, %.1fs", rc2, dt2)

    # Count success rates
    align_total, align_ok = count_success(align_summary)
    lmdb_summary = lmdb_dir / "structure_build_summary.csv"
    lmdb_total, lmdb_ok = count_success(lmdb_summary)

    # Write report
    report = args.results_dir / "downstream_validation.md"
    lines = [
        "# Downstream Validation Report", "",
        f"- Pilot dir: `{args.pilot_dir}`",
        f"- Pairs mapped: {n_mapped}", "",
        "## Step 8.2: Ligand Alignment",
        f"- Return code: {rc1}",
        f"- Runtime: {dt1:.1f}s",
        f"- Success: {align_ok}/{align_total}" if align_total else "- No alignment summary",
    ]
    if out1.strip():
        lines.extend(["", "```", tail(out1), "```"])

    lines.extend(["", "## Step 8.3: Structure LMDB",
        f"- Return code: {rc2}",
        f"- Runtime: {dt2:.1f}s",
        f"- Success: {lmdb_ok}/{lmdb_total}" if lmdb_total else "- No LMDB summary",
    ])
    if out2.strip():
        lines.extend(["", "```", tail(out2), "```"])

    overall = rc1 == 0 and rc2 == 0
    lines.extend(["", f"## Overall: {'PASS' if overall else 'FAIL'}"])
    report.write_text("\n".join(lines), encoding="utf-8")

    print(f"\n=== Downstream Validation ===")
    print(f"Alignment: {align_ok}/{align_total}, LMDB: {lmdb_ok}/{lmdb_total}")
    print(f"Overall: {'PASS' if overall else 'FAIL'}")
    print(f"Report: {report}")

    return 0 if overall else 1


if __name__ == "__main__":
    sys.exit(main())
