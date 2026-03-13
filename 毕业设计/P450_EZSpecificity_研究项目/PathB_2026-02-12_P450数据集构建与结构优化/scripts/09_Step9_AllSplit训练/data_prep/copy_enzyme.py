"""
Robust chunked copy of enzyme_features.lmdb from Google Drive.
Reads sequentially from byte 0 (no seek), handles network interruptions.
"""
import os
import sys
import time

SRC = r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\brenda\enzyme_features.lmdb"
DST = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化\data\09_Step9_AllSplit训练\brenda\enzyme_features.lmdb"
DST_TMP = DST + ".tmp"

CHUNK_SIZE = 16 * 1024 * 1024  # 16MB chunks (smaller = more resilient)
MAX_RETRIES = 20
RETRY_DELAY = 15


def copy_sequential():
    src_size = os.path.getsize(SRC)
    print(f"Source: {src_size/1024**3:.2f}GB ({src_size:,} bytes)")

    # Remove partial temp if it exists (no seek resume to avoid Google Drive issues)
    if os.path.exists(DST_TMP):
        existing = os.path.getsize(DST_TMP)
        print(f"Removing partial temp file ({existing/1024**3:.1f}GB)")
        os.remove(DST_TMP)

    retries = 0
    while retries < MAX_RETRIES:
        try:
            print(f"\nAttempt {retries+1}/{MAX_RETRIES}: copying from byte 0...")
            sys.stdout.flush()

            with open(SRC, 'rb') as fsrc:
                with open(DST_TMP, 'wb') as fdst:
                    copied = 0
                    t0 = time.time()
                    last_report = t0

                    while True:
                        try:
                            chunk = fsrc.read(CHUNK_SIZE)
                        except (OSError, IOError) as e:
                            # Flush what we have, then raise to trigger retry
                            fdst.flush()
                            raise

                        if not chunk:
                            break
                        fdst.write(chunk)
                        copied += len(chunk)

                        now = time.time()
                        if now - last_report >= 10:
                            elapsed = now - t0
                            rate = copied / elapsed / 1024**2
                            pct = copied / src_size * 100
                            eta = (src_size - copied) / (rate * 1024**2) if rate > 0 else 9999
                            print(f"  {copied/1024**3:.1f}/{src_size/1024**3:.1f}GB "
                                  f"({pct:.1f}%) rate={rate:.1f}MB/s ETA={eta/60:.1f}min")
                            sys.stdout.flush()
                            last_report = now

            # Verify size
            tmp_size = os.path.getsize(DST_TMP)
            if tmp_size != src_size:
                print(f"\n[WARN] Size mismatch: {tmp_size} != {src_size}")
                retries += 1
                continue

            # Success - rename to final
            if os.path.exists(DST):
                os.remove(DST)
            os.rename(DST_TMP, DST)

            elapsed = time.time() - t0
            print(f"\nCopy complete! {src_size/1024**3:.2f}GB in {elapsed:.0f}s "
                  f"({src_size/elapsed/1024**2:.1f}MB/s)")
            return True

        except (OSError, IOError) as e:
            retries += 1
            print(f"\n[ERROR] Copy failed at attempt {retries}/{MAX_RETRIES}: {e}")
            print(f"  Retrying from beginning in {RETRY_DELAY}s...")
            sys.stdout.flush()
            # Clean up partial temp
            if os.path.exists(DST_TMP):
                try:
                    os.remove(DST_TMP)
                except OSError:
                    pass
            time.sleep(RETRY_DELAY)

    print(f"[FATAL] Copy failed after {MAX_RETRIES} retries.")
    return False


if __name__ == '__main__':
    print("=" * 60)
    print("Robust copy: enzyme_features.lmdb (sequential, no seek)")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    success = copy_sequential()

    if success:
        print(f"\nDone at {time.strftime('%Y-%m-%d %H:%M:%S')}")
    else:
        print(f"\nFailed at {time.strftime('%Y-%m-%d %H:%M:%S')}")
        sys.exit(1)
