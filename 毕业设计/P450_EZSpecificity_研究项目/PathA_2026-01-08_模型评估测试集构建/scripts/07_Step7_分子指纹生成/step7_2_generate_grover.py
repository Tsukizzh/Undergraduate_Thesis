"""
Step 7.2: GROVER Fingerprint Generation
Orchestrates the 3-step GROVER pipeline:
  1) save_features.py  → atom/bond descriptors (.npz)
  2) build_vocab.py    → molecular vocabulary (.pkl)
  3) main.py fingerprint → neural fingerprint (LMDB)

Output: grover_fingerprint.lmdb with keys "0".."435"
  - 'embedding': (n_atoms, 2400) atom-level
  - 'total_embedding': (4885,) molecule-level

NOTE: LMDB on Windows cannot handle non-ASCII paths.
      We use a temp ASCII path then move the result.
"""
import argparse
import os
import pickle
import shutil
import subprocess
import sys
import time

import lmdb
import numpy as np
import pandas as pd

def _find_project_root(start):
    """Search upward for project root (contains src/ and saved_model/)."""
    p = os.path.abspath(start)
    while True:
        if os.path.isdir(os.path.join(p, "src")) and os.path.isdir(os.path.join(p, "saved_model")):
            return p
        parent = os.path.dirname(p)
        if parent == p:
            break
        p = parent
    raise RuntimeError("Cannot find project root (expected src/ + saved_model/)")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PATHA_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
PROJECT_ROOT = _find_project_root(SCRIPT_DIR)

GROVER_DIR = os.path.join(PROJECT_ROOT, "src", "other_softwares", "grover_software")
PYTHON_EXE = sys.executable

DEFAULT_INPUT = os.path.join(PATHA_ROOT, "data", "04_Step4_格式修正后数据", "Substrates.csv")
DEFAULT_OUTPUT_DIR = os.path.join(PATHA_ROOT, "data", "07_Step7_分子指纹")
DEFAULT_CHECKPOINT = os.path.join(PROJECT_ROOT, "data", "pretrain_model", "grover_large.pt")


def run_command(cmd, env, cwd, step_name):
    print(f"\n{'='*60}")
    print(f"  Running: {step_name}")
    print(f"  CMD: {' '.join(cmd)}")
    print(f"  CWD: {cwd}")
    print(f"{'='*60}")

    result = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True)

    if result.stdout:
        print(result.stdout)
    if result.stderr:
        for line in result.stderr.strip().split("\n"):
            if line.strip():
                print(f"  [stderr] {line}")

    if result.returncode != 0:
        print(f"\nERROR: {step_name} failed with return code {result.returncode}")
        sys.exit(1)

    print(f"  [{step_name}] Done.\n")
    return result


def build_vocab_low_memory(grover_dir, data_path, vocab_save_folder, dataset_name, num_workers=1):
    """
    Custom vocab builder with low memory footprint.
    The original build_vocab.py uses Pool(100) which crashes on Windows.
    We call MolVocab directly with num_workers=1.
    """
    print(f"\n{'='*60}")
    print(f"  Running: Step 3/4: build_vocab (low-memory mode)")
    print(f"  num_workers={num_workers}")
    print(f"{'='*60}")

    # Add GROVER to path for imports
    if grover_dir not in sys.path:
        sys.path.insert(0, grover_dir)

    from grover.data.torchvocab import MolVocab

    os.makedirs(vocab_save_folder, exist_ok=True)

    for vocab_type in ['atom', 'bond']:
        vocab_file = f"{vocab_type}_vocab.pkl"
        if dataset_name:
            vocab_file = f"{dataset_name}_{vocab_file}"
        vocab_save_path = os.path.join(vocab_save_folder, vocab_file)

        print(f"  Building {vocab_type} vocab...")
        vocab = MolVocab(
            file_path=data_path,
            max_size=None,
            min_freq=1,
            num_workers=num_workers,  # Low memory!
            vocab_type=vocab_type
        )
        print(f"    {vocab_type} vocab size: {len(vocab)}")
        vocab.save_vocab(vocab_save_path)
        print(f"    Saved: {vocab_save_path}")

    print(f"  [build_vocab] Done.\n")


def create_grover_csv(input_csv, output_csv):
    """Create single-column CSV for GROVER (SMILES in first column)."""
    df = pd.read_csv(input_csv)
    if "Substrate_Index" in df.columns:
        df = df.sort_values("Substrate_Index").reset_index(drop=True)
    grover_df = pd.DataFrame({"Substrate_SMILES": df["Substrate_SMILES"]})
    grover_df.to_csv(output_csv, index=False)
    print(f"Created {output_csv} with {len(grover_df)} SMILES")
    return len(grover_df)


def validate_lmdb(lmdb_path, expected_count):
    """Validate GROVER LMDB output."""
    db = lmdb.open(lmdb_path, map_size=10*(1024**3), create=False,
                   subdir=False, readonly=True, lock=False)

    with db.begin() as txn:
        keys = list(txn.cursor().iternext(values=False))

    actual_count = len(keys)
    print(f"\nLMDB Validation: {lmdb_path}")
    print(f"  Records: {actual_count} (expected {expected_count})")

    expected_keys = {str(i).encode() for i in range(expected_count)}
    actual_keys = set(keys)
    missing = expected_keys - actual_keys
    extra = actual_keys - expected_keys

    if missing:
        print(f"  WARNING: Missing keys: {sorted([k.decode() for k in missing])[:10]}...")
    if extra:
        print(f"  WARNING: Extra keys: {sorted([k.decode() for k in extra])[:10]}...")

    # Spot-check first 3 records
    print(f"\n  Spot-check:")
    for i in range(min(3, actual_count)):
        with db.begin() as txn:
            value = txn.get(str(i).encode())
            if value is None:
                print(f"    [{i}] MISSING!")
                continue
            data = pickle.loads(value)
            emb_shape = data.get("embedding", np.array([])).shape
            total_shape = data.get("total_embedding", np.array([])).shape
            print(f"    [{i}] embedding={emb_shape}, total_embedding={total_shape}")

    db.close()

    if actual_count != expected_count or missing:
        print(f"\n  VALIDATION FAILED!")
        return False

    print(f"\n  VALIDATION PASSED!")
    return True


def main():
    parser = argparse.ArgumentParser(description="Step 7.2: GROVER Fingerprint Generation")
    parser.add_argument("--input-csv", default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--checkpoint", default=DEFAULT_CHECKPOINT)
    parser.add_argument("--grover-dir", default=GROVER_DIR)
    parser.add_argument("--no-cuda", action="store_true", default=False)
    args = parser.parse_args()

    output_dir = args.output_dir
    grover_dir = args.grover_dir
    checkpoint = args.checkpoint

    os.makedirs(output_dir, exist_ok=True)

    # LMDB on Windows cannot handle non-ASCII paths.
    # Use a simple temp directory in C:\grover_temp for GROVER processing.
    temp_dir = r"C:\grover_temp"
    os.makedirs(temp_dir, exist_ok=True)
    print(f"\nUsing temp directory (ASCII-safe): {temp_dir}")

    # Paths - use temp dir for GROVER intermediate files
    grover_csv = os.path.join(temp_dir, "grover_substrates.csv")
    features_npz = os.path.join(temp_dir, "grover_substrates.npz")
    vocab_dir = os.path.join(temp_dir, "grover_vocab")
    dummy_npz = os.path.join(temp_dir, "grover_fingerprint.npz")
    temp_lmdb = os.path.join(temp_dir, "grover_fingerprint.lmdb")

    # Final output path
    final_lmdb = os.path.join(output_dir, "grover_fingerprint.lmdb")

    # Verify prerequisites
    if not os.path.isfile(checkpoint):
        sys.exit(f"ERROR: Checkpoint not found: {checkpoint}")
    if not os.path.isdir(grover_dir):
        sys.exit(f"ERROR: GROVER directory not found: {grover_dir}")

    # Environment: add GROVER to PYTHONPATH
    env = os.environ.copy()
    env["PYTHONPATH"] = grover_dir + os.pathsep + env.get("PYTHONPATH", "")

    # ---- Step 1: Create GROVER input CSV ----
    print("\n" + "="*60)
    print("  Step 1/4: Create GROVER input CSV")
    print("="*60)
    n_substrates = create_grover_csv(args.input_csv, grover_csv)

    # ---- Step 2: Generate features (save_features.py) ----
    save_features_script = os.path.join(grover_dir, "scripts", "save_features.py")
    run_command([
        PYTHON_EXE, save_features_script,
        "--data_path", grover_csv,
        "--save_path", features_npz,
        "--features_generator", "fgtasklabel",
        "--restart",
        "--sequential",
    ], env=env, cwd=grover_dir, step_name="Step 2/4: save_features.py")

    if not os.path.isfile(features_npz):
        sys.exit(f"ERROR: Features file not created: {features_npz}")

    # ---- Step 3: Build vocabulary (custom low-memory version) ----
    # Original build_vocab.py uses Pool(100) which crashes Windows with OOM
    build_vocab_low_memory(grover_dir, grover_csv, vocab_dir, "test", num_workers=1)

    # ---- Step 4: Generate fingerprints (main.py fingerprint) ----
    main_script = os.path.join(grover_dir, "main.py")
    fp_cmd = [
        PYTHON_EXE, main_script, "fingerprint",
        "--data_path", grover_csv,
        "--features_path", features_npz,
        "--checkpoint_path", checkpoint,
        "--fingerprint_source", "both",
        "--output_path", dummy_npz,
        "--save_lmdb_path", temp_lmdb,
        "--batch_size", "32",
    ]
    if args.no_cuda:
        fp_cmd.append("--no_cuda")

    run_command(fp_cmd, env=env, cwd=grover_dir, step_name="Step 4/4: main.py fingerprint")

    # ---- Move LMDB to final destination ----
    if not os.path.isfile(temp_lmdb):
        sys.exit(f"ERROR: LMDB not created: {temp_lmdb}")

    print(f"\nMoving LMDB from temp to final location...")
    if os.path.exists(final_lmdb):
        os.remove(final_lmdb)
    shutil.copy2(temp_lmdb, final_lmdb)
    print(f"  Copied: {temp_lmdb} -> {final_lmdb}")

    # Also copy intermediate files for reference
    for fname in ["grover_substrates.csv", "grover_substrates.npz"]:
        src = os.path.join(temp_dir, fname)
        dst = os.path.join(output_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, dst)
    vocab_dst = os.path.join(output_dir, "grover_vocab")
    if os.path.exists(vocab_dir):
        if os.path.exists(vocab_dst):
            shutil.rmtree(vocab_dst)
        shutil.copytree(vocab_dir, vocab_dst)

    # ---- Validation ----
    ok = validate_lmdb(final_lmdb, n_substrates)

    # ---- Cleanup temp dir ----
    print(f"\nCleaning up temp directory: {temp_dir}")
    shutil.rmtree(temp_dir, ignore_errors=True)

    # ---- Summary ----
    print(f"\n{'='*60}")
    print(f"GROVER Fingerprint Generation {'Complete' if ok else 'FAILED'}")
    print(f"{'='*60}")
    print(f"Input:       {args.input_csv}")
    print(f"Output LMDB: {final_lmdb}")
    print(f"Substrates:  {n_substrates}")
    print(f"Checkpoint:  {checkpoint}")
    print(f"CUDA:        {'No' if args.no_cuda else 'Auto'}")
    print(f"Status:      {'PASS' if ok else 'FAIL'}")
    print(f"{'='*60}")

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\nTotal elapsed: {time.time() - t0:.1f}s")
