"""
Wait for enzyme_features.lmdb copy to complete, then launch training.
Monitors copy progress and auto-launches main_training.py when ready.
"""
import os
import sys
import time
import subprocess

DATA_DIR = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化\data\09_Step9_AllSplit训练\brenda"
ENZYME_PATH = os.path.join(DATA_DIR, "enzyme_features.lmdb")
EXPECTED_SIZE = 53_916_635_136  # bytes (from Google Drive source)
TRAIN_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "main_training.py")

REQUIRED_FILES = {
    "enzyme_features.lmdb": EXPECTED_SIZE,
    "grover_fingerprint.lmdb": 9_602_629_632,
    "reaction_features.lmdb": 206_630_912,
    "morgan_fingerprint.npy": 324_690_048,
}

def check_all_files():
    """Verify all required data files are present and complete."""
    print("\n[Pre-flight] Checking all required data files...")
    all_ok = True
    for name, expected in REQUIRED_FILES.items():
        path = os.path.join(DATA_DIR, name)
        if not os.path.exists(path):
            print(f"  MISSING: {name}")
            all_ok = False
            continue
        actual = os.path.getsize(path)
        if actual < expected * 0.99:  # allow 1% tolerance
            print(f"  INCOMPLETE: {name} ({actual/1024**3:.1f}/{expected/1024**3:.1f}GB)")
            all_ok = False
        else:
            print(f"  OK: {name} ({actual/1024**3:.2f}GB)")

    # Check structure subdirectory
    str_dir = os.path.join(DATA_DIR, "structure")
    for name in ["structure_features.lmdb", "sequence_features.lmdb"]:
        path = os.path.join(str_dir, name)
        if not os.path.exists(path):
            print(f"  MISSING: structure/{name}")
            all_ok = False
        else:
            print(f"  OK: structure/{name} ({os.path.getsize(path)/1024**3:.2f}GB)")

    return all_ok


def wait_for_enzyme():
    print(f"\nWaiting for enzyme_features.lmdb to complete ({EXPECTED_SIZE/1024**3:.1f}GB)...")
    prev_size = 0
    prev_time = time.time()
    while True:
        if not os.path.exists(ENZYME_PATH):
            print(f"  File not found, waiting...")
            time.sleep(30)
            continue

        cur_size = os.path.getsize(ENZYME_PATH)
        cur_time = time.time()
        pct = cur_size / EXPECTED_SIZE * 100
        dt = cur_time - prev_time
        rate = (cur_size - prev_size) / dt / 1024**2 if prev_size > 0 and dt > 0 else 0
        eta = (EXPECTED_SIZE - cur_size) / (rate * 1024**2) if rate > 0 else 9999

        print(f"  {cur_size/1024**3:.1f}/{EXPECTED_SIZE/1024**3:.1f}GB ({pct:.1f}%) "
              f"rate={rate:.1f}MB/s ETA={eta/60:.1f}min")

        if cur_size >= EXPECTED_SIZE:
            print("  Enzyme file complete!")
            return True

        prev_size = cur_size
        prev_time = cur_time
        time.sleep(15)


def main():
    print("=" * 60)
    print("Step 9: Wait for data + Auto-launch training")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    wait_for_enzyme()
    time.sleep(5)  # flush delay

    if not check_all_files():
        print("\n[ERROR] Some files are missing or incomplete. Aborting.")
        sys.exit(1)

    print()
    print("=" * 60)
    print("Launching main_training.py...")
    print(f"Time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    result = subprocess.run(
        [sys.executable, TRAIN_SCRIPT],
        cwd=r"D:\EZSpecificity_Project\src"
    )

    print(f"\nTraining exited with code: {result.returncode}")
    print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == '__main__':
    main()
