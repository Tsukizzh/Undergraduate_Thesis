"""Phase 1: Build grover_fingerprint_fixed.lmdb via pure key rewriting.

Bug: grover_fingerprint.lmdb was written with sequential counter keys after
*[H] (Substrate Index 8) was manually removed from grover_substrates.csv.
As a result:
    old key k, k < 8   -> Substrate Index k   (correct)
    old key k, k >= 8  -> Substrate Index k+1 (shifted by +1)
    Substrate Index 8  -> no LMDB entry anywhere
    Substrate Index 2124 -> no LMDB entry (shifted off the end)

Fix (pure rekey, no re-running GROVER):
    new key k = k            if k < 8
    new key k = k + 1        if k >= 8
New LMDB keyset = {0..7, 9..2124}, with b"8" absent.

Codex-reviewed (session 019d88a7). Key review points:
    - bytes(v) copies payload out of mmap (sufficient)
    - txn.put(overwrite=False) returns False on duplicate (not raise) - must assert
    - map_size=2 GB is ample for ~452 MB source
    - one write txn is fine at 2124 entries
    - close() auto-flushes; sync() optional
"""
import os
import shutil
import sys
import lmdb

OLD = "/root/autodl-tmp/EZSpecificity/PathC/P450/data/features/grover_fingerprint.lmdb"
NEW = "/root/autodl-tmp/EZSpecificity/PathC/P450/data/features/grover_fingerprint_fixed.lmdb"
MAP_SIZE = 2 * 1024 ** 3  # 2 GB

assert os.path.exists(OLD), f"missing source: {OLD}"

# Delete existing NEW to avoid stale-data contamination
if os.path.exists(NEW):
    print(f"Removing existing {NEW}")
    os.remove(NEW)
    lock_file = NEW + "-lock"
    if os.path.exists(lock_file):
        os.remove(lock_file)

# Also handle subdir-style stale paths just in case
if os.path.isdir(NEW):
    shutil.rmtree(NEW)


def remap(old_int: int) -> int:
    return old_int if old_int < 8 else old_int + 1


# --- Read old LMDB (keys + payloads into memory) ---
env_old = lmdb.open(
    OLD, subdir=False, readonly=True, lock=False, readahead=False, max_readers=4
)
items = []
with env_old.begin(write=False) as txn:
    old_key_ints = []
    for k, v in txn.cursor():
        old_int = int(k.decode())
        old_key_ints.append(old_int)
        new_int = remap(old_int)
        items.append((str(new_int).encode(), bytes(v)))  # copy out of mmap

print(f"Loaded {len(items)} entries from old LMDB")
print(f"  old key range: {min(old_key_ints)}..{max(old_key_ints)}")

assert len(items) == 2124, f"expected 2124 old keys, got {len(items)}"
assert min(old_key_ints) == 0 and max(old_key_ints) == 2123, (
    f"expected old keys 0..2123, got {min(old_key_ints)}..{max(old_key_ints)}"
)

# Verify the remap produces no collisions
new_keys_seen = set()
for nk, _ in items:
    assert nk not in new_keys_seen, f"duplicate new key: {nk!r}"
    new_keys_seen.add(nk)
print(f"  remap uniqueness check: OK ({len(new_keys_seen)} unique new keys)")

# Sort by integer key for tidy write order
items.sort(key=lambda p: int(p[0].decode()))
print(f"  new key range: {items[0][0].decode()}..{items[-1][0].decode()}")

# --- Write new LMDB ---
env_new = lmdb.open(
    NEW,
    subdir=False,
    readonly=False,
    create=True,
    map_size=MAP_SIZE,
    max_dbs=0,
    lock=True,
)
with env_new.begin(write=True) as txn:
    for nk, nv in items:
        ok = txn.put(nk, nv, overwrite=False)
        assert ok, f"txn.put returned False for key {nk!r} (duplicate?)"

env_old.close()
env_new.close()

print("\n--- Reopening NEW lmdb read-only for sanity check ---")
env_check = lmdb.open(
    NEW, subdir=False, readonly=True, lock=False, readahead=False, max_readers=4
)
env_old_check = lmdb.open(
    OLD, subdir=False, readonly=True, lock=False, readahead=False, max_readers=4
)

with env_check.begin() as tn, env_old_check.begin() as to:
    # Count
    n = tn.stat()["entries"]
    print(f"entries: {n}")
    assert n == 2124, f"expected 2124 entries, got {n}"

    # b"8" absent
    assert tn.get(b"8") is None, "new key '8' should NOT exist (*[H] gap)"
    print("key '8': ABSENT (correct)")

    # Boundary byte-level comparisons
    checks = [
        ("7", "7", True),      # unchanged
        ("8", "9", True),      # shift start: old 8 -> new 9
        ("100", "101", True),  # mid: old 100 -> new 101
        ("2123", "2124", True),  # end: old 2123 -> new 2124
    ]
    for old_k, new_k, expect_match in checks:
        ov = to.get(old_k.encode())
        nv = tn.get(new_k.encode())
        assert ov is not None, f"old['{old_k}'] missing"
        assert nv is not None, f"new['{new_k}'] missing"
        match = ov == nv
        status = "OK" if match == expect_match else "FAIL"
        print(f"  old[{old_k!r}] == new[{new_k!r}]: {match}  [{status}]")
        assert match == expect_match, f"mismatch at {old_k}->{new_k}"

    # Key set completeness: {0..7, 9..2124}
    actual_keys = set()
    for k, _ in tn.cursor():
        actual_keys.add(int(k.decode()))
    expected_keys = set(range(0, 8)) | set(range(9, 2125))
    assert actual_keys == expected_keys, (
        f"key set mismatch: missing={expected_keys-actual_keys} "
        f"extra={actual_keys-expected_keys}"
    )
    print(f"key set: exactly {{0..7, 9..2124}} (2124 keys)")

    # Additional: full-sweep byte equality for the shift region
    #   old k   ->  new (k+1)  for all k in [8..2123]
    print("\nFull-sweep byte check: old[k] == new[k+1] for k in 8..2123 ...")
    mismatches = 0
    for k in range(8, 2124):
        ov = to.get(str(k).encode())
        nv = tn.get(str(k + 1).encode())
        if ov != nv:
            mismatches += 1
    print(f"  mismatches: {mismatches} / 2116")
    assert mismatches == 0

    # And for k in [0..7], old[k] == new[k]
    print("Full-sweep byte check: old[k] == new[k] for k in 0..7 ...")
    mismatches = 0
    for k in range(0, 8):
        ov = to.get(str(k).encode())
        nv = tn.get(str(k).encode())
        if ov != nv:
            mismatches += 1
    print(f"  mismatches: {mismatches} / 8")
    assert mismatches == 0

env_check.close()
env_old_check.close()

print("\n=== Phase 1 DONE ===")
print(f"  source:  {OLD}")
print(f"  output:  {NEW}")
sz = os.path.getsize(NEW)
print(f"  size:    {sz / 1024**2:.1f} MB")
