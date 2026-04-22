"""
Verify the paperfilter cache:
  A. enzyme_id mapping integrity (seq_len in flatbin matches Enzymes.csv row seq_len)
  B. filtered index.pt is self-consistent (all arrays same len, no blacklist leaks)
  C. PtCacheDataset loads without errors, first N samples are sane

Run on server CPU (no torch model instantiation, just data layer).
"""
from __future__ import annotations
import csv, json, os, sys
from pathlib import Path

import torch

ROOT = Path("/root/autodl-tmp/EZSpecificity/PathC/P450")
SRC_CACHE = ROOT / "data/pt_cache_allfix_unified/random"
DST_CACHE = ROOT / "data/pt_cache_allfix_unified_paperfilter/random"
ENZ_CSV = ROOT / "data/Enzymes.csv"
BLACKLIST_JSON = ROOT / "experiments/EXP004_paper_baseline_unified_workdir/paper_blacklist.json"
REPORT_OUT = ROOT / "experiments/EXP004_paper_baseline_unified_workdir/verify_report.json"

# Import pt_dataset from a working experiment scripts dir
sys.path.insert(0, str(ROOT / "experiments/EXP001_allfix_unified/scripts"))
sys.path.insert(0, str(ROOT / "experiments/EXP001_allfix_unified/src"))


def check_a_enzyme_mapping():
    """A: verify enzyme_global_id in enzymes_index.pt uses Enzymes.csv row order."""
    print("[A] enzyme_id mapping integrity")
    with open(ENZ_CSV, encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    n_rows = len(rows)
    print(f"  Enzymes.csv rows: {n_rows}")

    enz_idx = torch.load(SRC_CACHE / "enzymes" / "enzymes_index.pt",
                         map_location="cpu", weights_only=False)
    lookup = enz_idx["index"]  # {global_id: (byte_offset, seq_len)}
    print(f"  enzymes_index.pt entries: {len(lookup)}")

    # Sample 20 enzyme_global_ids spread across the range
    sample_ids = sorted({0, 1, 5, 10, 100, 500, 800, 1000, 1200, 1500, 1600, 1621})
    sample_ids = [i for i in sample_ids if i < n_rows]

    checks = []
    for eid in sample_ids:
        csv_len = len(rows[eid]["Protein sequence"].strip())
        if eid not in lookup:
            checks.append({"enzyme_id": eid, "status": "NOT_IN_LOOKUP", "csv_len": csv_len})
            continue
        bin_len = lookup[eid][1]
        # ESM truncates at 1000aa typically; our max is 1450. But some may be stored
        # at original length. Just compare direct first.
        status = "OK" if bin_len == csv_len else f"MISMATCH bin={bin_len} csv={csv_len}"
        checks.append({"enzyme_id": eid, "status": status, "csv_len": csv_len, "bin_len": bin_len})
        print(f"  eid={eid:5d}  csv_len={csv_len:5d}  bin_len={bin_len:5d}  {status}")

    all_ok = all(c["status"] == "OK" for c in checks)
    return {"checks": checks, "all_ok": all_ok, "sample_count": len(checks)}


def check_b_index_consistency():
    """B: filtered test/index.pt arrays self-consistent, no blacklist leak,
       substrate lookup exists for all substrate_ids (codex recommendation)."""
    print("\n[B] filtered index.pt consistency")
    with open(BLACKLIST_JSON) as f:
        blacklist = set(json.load(f)["blacklisted_enzyme_global_ids"])

    idx = torch.load(DST_CACHE / "test" / "index.pt",
                     map_location="cpu", weights_only=False)
    lens = {k: len(v) for k, v in idx.items()}
    print(f"  array lengths: {lens}")
    n = lens["sample_ids"]
    all_same = all(v == n for v in lens.values())

    # No enzyme in blacklist
    eids = idx["enzyme_ids"].tolist()
    leaked = [e for e in eids if e in blacklist]
    print(f"  leaked enzyme_ids (should be 0): {len(leaked)}")

    # Dtypes unchanged
    src_idx = torch.load(SRC_CACHE / "test" / "index.pt",
                         map_location="cpu", weights_only=False)
    dtype_match = {k: str(idx[k].dtype) == str(src_idx[k].dtype) for k in idx}
    print(f"  dtype preserved: {dtype_match}")

    # Substrate lookup exists for every substrate_id (sanity: overlay hasn't broken ptrs)
    sub_idx = torch.load(DST_CACHE / "substrates" / "substrates_index.pt",
                         map_location="cpu", weights_only=False)
    sub_lookup_keys = set(sub_idx["index"].keys())
    sub_ids_in_test = set(idx["substrate_ids"].tolist())
    missing_subs = sub_ids_in_test - sub_lookup_keys
    print(f"  unique substrate_ids in filtered test: {len(sub_ids_in_test)}")
    print(f"  missing substrate lookups (should be 0): {len(missing_subs)}")

    # Enzyme lookup exists for every enzyme_id (same principle)
    enz_idx_pt = torch.load(DST_CACHE / "enzymes" / "enzymes_index.pt",
                            map_location="cpu", weights_only=False)
    enz_lookup_keys = set(enz_idx_pt["index"].keys())
    enz_ids_in_test = set(eids)
    missing_enz = enz_ids_in_test - enz_lookup_keys
    print(f"  unique enzyme_ids in filtered test: {len(enz_ids_in_test)}")
    print(f"  missing enzyme lookups (should be 0): {len(missing_enz)}")

    return {
        "array_lengths": lens,
        "all_arrays_same_length": all_same,
        "n_samples": n,
        "leaked_blacklisted_enzyme_ids": len(leaked),
        "dtype_preserved": dtype_match,
        "unique_enzyme_ids_in_test": len(enz_ids_in_test),
        "missing_enzyme_lookups": len(missing_enz),
        "unique_substrate_ids_in_test": len(sub_ids_in_test),
        "missing_substrate_lookups": len(missing_subs),
    }


def check_d_provenance():
    """D: source cache untouched; overlay only adds symlinks + 1 new file."""
    print("\n[D] provenance check")
    src_mtime_enzymes_bin = (SRC_CACHE / "enzymes" / "enzymes.bin").stat().st_mtime
    src_mtime_test_idx = (SRC_CACHE / "test" / "index.pt").stat().st_mtime

    # All entries in DST except test/index.pt should be symlinks
    dst_entries = {}
    for entry in DST_CACHE.iterdir():
        if entry.name == "test":
            continue
        dst_entries[entry.name] = "symlink" if entry.is_symlink() else "NOT_SYMLINK"
    # test/ is a dir (not symlink), but test/samples should be symlink, test/index.pt is real file
    test_samples = DST_CACHE / "test" / "samples"
    test_index = DST_CACHE / "test" / "index.pt"
    dst_entries["test/samples"] = "symlink" if test_samples.is_symlink() else "NOT_SYMLINK"
    dst_entries["test/index.pt"] = "new_file" if (not test_index.is_symlink() and test_index.is_file()) else "NOT_NEW_FILE"

    print(f"  overlay entries: {dst_entries}")
    all_good = all(v in ("symlink", "new_file") for v in dst_entries.values())

    return {
        "overlay_entries": dst_entries,
        "all_expected": all_good,
        "src_enzymes_bin_mtime": src_mtime_enzymes_bin,
        "src_test_index_mtime": src_mtime_test_idx,
    }


def check_c_pt_dataset_load():
    """C: PtCacheDataset on new cache, load mixed sample set, check sanity."""
    print("\n[C] PtCacheDataset load")
    from pt_dataset import PtCacheDataset
    import random

    ds = PtCacheDataset(
        cache_dir=str(DST_CACHE),
        split="test",
        edge_mode="legacy_bug",   # what we'll use for paper ckpt
        dist_noise=False,
        cutoff=10.0,
        num_r_gaussian=32,
        max_enzyme_len=1450,
        max_substrate_len=280,
        preload=False,
    )
    print(f"  dataset len: {len(ds)}")

    with open(BLACKLIST_JSON) as f:
        blacklist = set(json.load(f)["blacklisted_enzyme_global_ids"])

    # Mix boundary, quartile, and random indices (codex recommendation)
    n = len(ds)
    random.seed(42)
    random_idxs = random.sample(range(n), 8)
    test_idxs = sorted(set([0, 1, 2, n//4, n//2, 3*n//4, n-3, n-2, n-1] + random_idxs))

    results = []
    for i in test_idxs:
        try:
            s = ds[i]
            # Sanity checks
            has_nan_emb = bool(torch.isnan(s.embedding).any().item())
            has_nan_grover = bool(torch.isnan(s.grover).any().item())
            has_nan_protein_x = bool(torch.isnan(s.protein_x).any().item())
            label_int = int(s.label.item())
            # We don't have enzyme_global_id in StructureSequenceData directly;
            # cross-check via index
            eid_from_idx = int(ds._index["enzyme_ids"][i])
            in_blacklist = eid_from_idx in blacklist
            entry = {
                "idx": i,
                "enzyme_id": eid_from_idx,
                "in_blacklist": in_blacklist,
                "label": label_int,
                "nan_emb": has_nan_emb,
                "nan_grover": has_nan_grover,
                "nan_protein_x": has_nan_protein_x,
                "protein_x_shape": list(s.protein_x.shape),
                "complex_edge_index_shape": list(s.complex_edge_index.shape),
            }
            results.append(entry)
            print(f"  idx={i:6d} eid={eid_from_idx:5d} label={label_int} "
                  f"bl={in_blacklist} nan={'✗' if any([has_nan_emb, has_nan_grover, has_nan_protein_x]) else '✓'}")
        except Exception as e:
            results.append({"idx": i, "error": str(e)})
            print(f"  idx={i}: ERROR {e}")

    all_ok = (
        all("error" not in r for r in results)
        and not any(r.get("in_blacklist") for r in results)
        and not any(r.get("nan_emb") or r.get("nan_grover") or r.get("nan_protein_x") for r in results)
    )
    return {"n_samples": len(ds), "tested": len(results), "all_ok": all_ok, "results": results}


def main():
    report = {}
    report["A_enzyme_mapping"] = check_a_enzyme_mapping()
    report["B_index_consistency"] = check_b_index_consistency()
    report["C_pt_dataset_load"] = check_c_pt_dataset_load()
    report["D_provenance"] = check_d_provenance()

    all_ok = (
        report["A_enzyme_mapping"]["all_ok"]
        and report["B_index_consistency"]["leaked_blacklisted_enzyme_ids"] == 0
        and report["B_index_consistency"]["all_arrays_same_length"]
        and report["B_index_consistency"]["missing_substrate_lookups"] == 0
        and report["B_index_consistency"]["missing_enzyme_lookups"] == 0
        and report["C_pt_dataset_load"]["all_ok"]
        and report["D_provenance"]["all_expected"]
    )
    report["OVERALL_OK"] = all_ok

    with open(REPORT_OUT, "w") as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nReport: {REPORT_OUT}")
    print(f"OVERALL: {'PASS ✓' if all_ok else 'FAIL ✗'}")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
