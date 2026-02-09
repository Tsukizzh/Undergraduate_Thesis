"""
Step 4: CSV Format Correction for EZSpecificity Model

This script fixes the CSV format issues in B1-B6 datasets:
1. Removes UTF-8 BOM (Byte Order Mark)
2. Ensures all index columns are integer type
3. Validates column names match EZSpecificity requirements

Required columns: Enzyme Index, Substrate Index, Label, Dock Index

Usage:
    python step4_csv_format_fix.py [--dry-run]
"""

from pathlib import Path
import pandas as pd
import sys
import shutil
from datetime import datetime


def check_bom(file_path: Path) -> bool:
    """Check if file has UTF-8 BOM."""
    with open(file_path, 'rb') as f:
        return f.read(3) == b'\xef\xbb\xbf'


def fix_csv_format(input_path: Path, output_path: Path) -> dict:
    """
    Fix CSV format: remove BOM and ensure integer types.

    Returns:
        dict with validation results
    """
    result = {
        'input': str(input_path),
        'output': str(output_path),
        'had_bom': check_bom(input_path),
        'rows': 0,
        'columns': [],
        'dtypes_before': {},
        'dtypes_after': {},
        'issues': []
    }

    # Read CSV (pandas handles BOM automatically)
    df = pd.read_csv(input_path)
    result['rows'] = len(df)
    result['columns'] = df.columns.tolist()
    result['dtypes_before'] = {col: str(df[col].dtype) for col in df.columns}

    # Validate required columns
    required_cols = {'Enzyme Index', 'Substrate Index', 'Label', 'Dock Index'}
    missing = required_cols - set(df.columns)
    if missing:
        result['issues'].append(f"Missing columns: {missing}")

    # Ensure integer types for index columns
    int_cols = ['Enzyme Index', 'Substrate Index', 'Label', 'Dock Index']
    for col in int_cols:
        if col in df.columns:
            if df[col].dtype != 'int64':
                df[col] = df[col].astype('int64')

    result['dtypes_after'] = {col: str(df[col].dtype) for col in df.columns}

    # Write without BOM (standard UTF-8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False, encoding='utf-8')

    # Verify no BOM in output
    if check_bom(output_path):
        result['issues'].append("Output still has BOM!")

    return result


def main():
    dry_run = '--dry-run' in sys.argv

    # Base paths
    base_dir = Path(__file__).resolve().parent.parent.parent
    source_dir = base_dir / 'data' / '03_Step3_配体处理与数据集构建' / '06_数据集切分'
    output_dir = base_dir / 'data' / '04_Step4_格式修正后数据'

    # Dataset definitions
    datasets = [
        ('B1_仅底物_272pos', 'data.csv'),
        ('B2_底物正产物负_272pos23neg', 'data.csv'),
        ('B3_完整数据集_272pos267neg', 'data.csv'),
        ('B4_仅产物_23neg', 'data.csv'),
        ('B5_仅抑制剂_244neg', 'data.csv'),
        ('B6_底物正抑制剂负_272pos244neg', 'data.csv'),
    ]

    # Also process Substrates.csv
    substrates_file = ('', 'Substrates.csv')

    print(f"Step 4: CSV Format Correction")
    print(f"{'=' * 50}")
    print(f"Source: {source_dir}")
    print(f"Output: {output_dir}")
    print(f"Dry run: {dry_run}")
    print()

    results = []

    # Process each dataset
    for folder, filename in datasets:
        input_path = source_dir / folder / filename
        output_folder = output_dir / folder
        output_path = output_folder / filename

        if not input_path.exists():
            print(f"[SKIP] {folder}/{filename} - not found")
            continue

        print(f"[PROCESS] {folder}/{filename}")

        if dry_run:
            has_bom = check_bom(input_path)
            df = pd.read_csv(input_path)
            print(f"  - Rows: {len(df)}, BOM: {has_bom}")
            print(f"  - Columns: {df.columns.tolist()}")
        else:
            result = fix_csv_format(input_path, output_path)
            results.append(result)
            print(f"  - Rows: {result['rows']}, Had BOM: {result['had_bom']}")
            if result['issues']:
                print(f"  - Issues: {result['issues']}")

    # Process Substrates.csv
    substrates_input = source_dir / 'Substrates.csv'
    if substrates_input.exists():
        print(f"\n[PROCESS] Substrates.csv")
        substrates_output = output_dir / 'Substrates.csv'

        if dry_run:
            has_bom = check_bom(substrates_input)
            df = pd.read_csv(substrates_input)
            print(f"  - Rows: {len(df)}, BOM: {has_bom}")
        else:
            # Copy Substrates.csv to output (just remove BOM if present)
            df = pd.read_csv(substrates_input)
            output_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(substrates_output, index=False, encoding='utf-8')
            print(f"  - Rows: {len(df)}, Copied to output")

    # Copy Enzymes.csv from Step 2
    enzymes_source = base_dir / 'data' / '02_Step2_酶序列' / 'Enzymes.csv'
    if enzymes_source.exists():
        print(f"\n[COPY] Enzymes.csv from Step 2")
        enzymes_output = output_dir / 'Enzymes.csv'

        if not dry_run:
            df = pd.read_csv(enzymes_source)
            output_dir.mkdir(parents=True, exist_ok=True)
            df.to_csv(enzymes_output, index=False, encoding='utf-8')
            print(f"  - Rows: {len(df)}, Copied to output")
        else:
            df = pd.read_csv(enzymes_source)
            print(f"  - Rows: {len(df)}, Would copy to output")

    print(f"\n{'=' * 50}")
    print(f"Completed! Output directory: {output_dir}")

    if not dry_run:
        # Write summary report
        report_path = output_dir / 'format_fix_report.txt'
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(f"Step 4 CSV Format Fix Report\n")
            f.write(f"Generated: {datetime.now().isoformat()}\n")
            f.write(f"{'=' * 50}\n\n")
            for r in results:
                f.write(f"Dataset: {Path(r['input']).parent.name}\n")
                f.write(f"  Rows: {r['rows']}\n")
                f.write(f"  Had BOM: {r['had_bom']}\n")
                f.write(f"  Columns: {r['columns']}\n")
                f.write(f"  Issues: {r['issues'] or 'None'}\n\n")
        print(f"Report saved: {report_path}")


if __name__ == '__main__':
    main()
