#!/usr/bin/env python3
"""Prepare EXP003 full-P450 Fe/HEM overlay manifests.

This script is intentionally manifest-only. It does not modify PT caches,
training scripts, checkpoints, or source data. The goal is to make the EXP003
data boundary explicit before cache construction starts.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import time
from collections import Counter, defaultdict
from pathlib import Path


def now_tag() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def as_bool(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text in {"true", "1", "yes", "y"}


def as_int(value: object) -> int:
    text = str(value or "").strip()
    if not text:
        return 0
    return int(float(text))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def sample_path(cache_dir: Path, split: str, sample_id: int) -> Path:
    return cache_dir / split / "samples" / f"{sample_id // 1000:03d}" / f"sample_{sample_id:06d}.pt"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p450-csv", required=True)
    parser.add_argument("--structure-records-csv", required=True)
    parser.add_argument("--best1-audit-csv", required=True)
    parser.add_argument("--google-manifest-csv", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--normalized-splits", default=None)
    parser.add_argument("--base-cache", default=None)
    parser.add_argument("--tag", default=None)
    args = parser.parse_args()

    tag = args.tag or now_tag()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    p450_path = Path(args.p450_csv)
    structure_path = Path(args.structure_records_csv)
    audit_path = Path(args.best1_audit_csv)
    google_path = Path(args.google_manifest_csv)

    p450_rows = read_csv(p450_path)
    structure_rows = read_csv(structure_path)
    audit_rows = read_csv(audit_path)
    google_rows = read_csv(google_path)

    p450_by_uniprot = {}
    duplicate_p450_ids: list[str] = []
    for row in p450_rows:
        uniprot = str(row.get("uniprot_id") or row.get("uniprot") or "").strip()
        if not uniprot:
            continue
        if uniprot in p450_by_uniprot:
            duplicate_p450_ids.append(uniprot)
        p450_by_uniprot[uniprot] = row
    p450_ids = set(p450_by_uniprot)

    google_by_structure = {
        str(as_int(row.get("structure_index"))): row
        for row in google_rows
        if str(row.get("structure_index") or "").strip()
    }

    structure_overlap = []
    structure_by_index = {}
    structure_stats = defaultdict(lambda: {"rows": 0, "samples": 0, "positive": 0})
    for row in structure_rows:
        uniprot = str(row.get("uniprot") or row.get("uniprot_id") or "").strip()
        if uniprot not in p450_ids:
            continue
        out = dict(row)
        structure_overlap.append(out)
        structure_by_index[str(as_int(row.get("structure_index")))] = out
        stat = structure_stats[uniprot]
        stat["rows"] += 1
        stat["samples"] += as_int(row.get("n_samples"))
        stat["positive"] += as_int(row.get("n_positive"))

    target_candidates = []
    all_p450_best1 = []
    target_stats = defaultdict(lambda: {"rows": 0, "samples": 0, "positive": 0})
    for row in audit_rows:
        uniprot = str(row.get("uniprot") or row.get("uniprot_id") or "").strip()
        if uniprot not in p450_ids:
            continue
        structure_key = str(as_int(row.get("structure_index")))
        google = google_by_structure.get(structure_key, {})
        merged = dict(row)
        merged["found_in_drive_metadata"] = google.get("found_in_drive_metadata", "")
        merged["drive_path"] = google.get("drive_path", "")
        merged["pdb_name"] = google.get("pdb_name", "")
        all_p450_best1.append(merged)

        usable = as_bool(row.get("exp002_atom_audit_usable_target_heme_fe"))
        found_official_pdb = as_bool(google.get("found_in_drive_metadata"))
        if usable and found_official_pdb:
            target_candidates.append(merged)
            stat = target_stats[uniprot]
            stat["rows"] += 1
            stat["samples"] += as_int(row.get("n_samples"))
            stat["positive"] += as_int(row.get("n_positive"))

    uniprot_summary = []
    missing_structure = []
    missing_target = []
    for uniprot in sorted(p450_ids):
        p450_row = p450_by_uniprot[uniprot]
        s = structure_stats[uniprot]
        t = target_stats[uniprot]
        record = {
            "uniprot_id": uniprot,
            "has_heme_in_uniprot": p450_row.get("has_heme", ""),
            "has_esibank_structure_records": s["rows"] > 0,
            "n_esibank_structure_records": s["rows"],
            "n_esibank_samples": s["samples"],
            "n_esibank_positive_samples": s["positive"],
            "has_fe_overlay_target": t["rows"] > 0,
            "n_fe_overlay_target_rows": t["rows"],
            "n_fe_overlay_target_samples": t["samples"],
            "n_fe_overlay_target_positive_samples": t["positive"],
            "protein_families": p450_row.get("protein_families", ""),
            "cc_cofactor": p450_row.get("cc_cofactor", ""),
            "xref_interpro": p450_row.get("xref_interpro", ""),
        }
        uniprot_summary.append(record)
        if s["rows"] == 0:
            missing_structure.append(p450_row)
        if t["rows"] == 0:
            missing_target.append(p450_row)

    split_manifest = []
    trainable_split_manifest = []
    split_summary_rows = []
    split_counts = Counter()
    split_positive = Counter()
    split_cache_missing = Counter()
    split_unique_uniprots = defaultdict(set)
    split_unique_substrates = defaultdict(set)
    if args.normalized_splits:
        split_root = Path(args.normalized_splits)
        base_cache = Path(args.base_cache) if args.base_cache else None
        cache_sample_ids: dict[str, set[int]] = {}
        if base_cache:
            try:
                import torch
            except Exception as exc:  # pragma: no cover - environment dependent
                raise SystemExit(f"--base-cache requires torch import: {exc!r}") from exc
            for split in ["train", "val", "test"]:
                index_path = base_cache / split / "index.pt"
                if index_path.exists():
                    index = torch.load(index_path, map_location="cpu", weights_only=False)
                    cache_sample_ids[split] = set(int(x) for x in index["sample_ids"].tolist())
                else:
                    cache_sample_ids[split] = set()

        for split in ["train", "val", "test"]:
            split_csv = split_root / "brenda" / f"{split}.csv"
            if not split_csv.exists():
                continue
            with split_csv.open("r", encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                for sample_id, row in enumerate(reader):
                    dock_index = str(as_int(row.get("Dock Index")))
                    structure = structure_by_index.get(dock_index)
                    if not structure:
                        continue
                    uniprot = str(structure.get("uniprot") or "").strip()
                    label = as_int(row.get("Label"))
                    cache_present = ""
                    cache_file_exists = ""
                    cache_path = ""
                    if base_cache:
                        cache_present = sample_id in cache_sample_ids.get(split, set())
                        path = sample_path(base_cache, split, sample_id)
                        cache_file_exists = path.exists()
                        cache_path = str(path)
                    record = {
                        "split": split,
                        "sample_id": sample_id,
                        "dock_index": dock_index,
                        "enzyme_index": row.get("Enzyme Index", ""),
                        "substrate_index": row.get("Substrate Index", ""),
                        "label": label,
                        "uniprot": uniprot,
                        "structure_record_substrate": structure.get("substrate", ""),
                        "structure_record_n_samples": structure.get("n_samples", ""),
                        "structure_record_n_positive": structure.get("n_positive", ""),
                        "cache_sample_id_present": cache_present,
                        "cache_file_exists": cache_file_exists,
                        "cache_path": cache_path,
                    }
                    split_manifest.append(record)
                    if not base_cache or (cache_present and cache_file_exists):
                        trainable_split_manifest.append(record)
                    split_counts[split] += 1
                    split_positive[split] += label
                    split_unique_uniprots[split].add(uniprot)
                    split_unique_substrates[split].add(str(row.get("Substrate Index", "")))
                    if base_cache and not (cache_present and cache_file_exists):
                        split_cache_missing[split] += 1

        for split in ["train", "val", "test"]:
            total = split_counts[split]
            positive = split_positive[split]
            split_summary_rows.append(
                {
                    "split": split,
                    "samples": total,
                    "positive": positive,
                    "negative": total - positive,
                    "positive_rate": round(positive / total, 6) if total else "",
                    "unique_uniprot": len(split_unique_uniprots[split]),
                    "unique_substrate_index": len(split_unique_substrates[split]),
                    "cache_missing": split_cache_missing[split],
                }
            )

    source_copy = out_dir / f"p450_389_uniprot_source_20260102_copy_{tag}.csv"
    if not source_copy.exists():
        shutil.copy2(p450_path, source_copy)

    def prefixed(name: str) -> Path:
        return out_dir / f"exp003_{name}_{tag}.csv"

    structure_fields = list(structure_overlap[0].keys()) if structure_overlap else []
    best1_fields = list(all_p450_best1[0].keys()) if all_p450_best1 else []
    target_fields = list(target_candidates[0].keys()) if target_candidates else best1_fields
    p450_fields = list(p450_rows[0].keys()) if p450_rows else []
    uniprot_fields = list(uniprot_summary[0].keys()) if uniprot_summary else []

    structure_out = prefixed("p450389_overlap_structure_records")
    all_best1_out = prefixed("p450389_overlap_best1_audited")
    target_out = prefixed("fe_overlay_target_candidates")
    uniprot_out = prefixed("p450389_uniprot_summary")
    missing_structure_out = prefixed("missing_from_esibank_structure_records")
    missing_target_out = prefixed("missing_from_fe_overlay_targets")
    split_manifest_out = prefixed("p450389_split_sample_manifest")
    trainable_split_manifest_out = prefixed("p450389_trainable_split_sample_manifest")
    split_summary_out = prefixed("p450389_split_summary")

    write_csv(structure_out, structure_overlap, structure_fields)
    write_csv(all_best1_out, all_p450_best1, best1_fields)
    write_csv(target_out, target_candidates, target_fields)
    write_csv(uniprot_out, uniprot_summary, uniprot_fields)
    write_csv(missing_structure_out, missing_structure, p450_fields)
    write_csv(missing_target_out, missing_target, p450_fields)
    split_manifest_fields = list(split_manifest[0].keys()) if split_manifest else [
        "split", "sample_id", "dock_index", "enzyme_index", "substrate_index",
        "label", "uniprot", "structure_record_substrate",
        "structure_record_n_samples", "structure_record_n_positive",
        "cache_sample_id_present", "cache_file_exists", "cache_path",
    ]
    split_summary_fields = list(split_summary_rows[0].keys()) if split_summary_rows else [
        "split", "samples", "positive", "negative", "positive_rate",
        "unique_uniprot", "unique_substrate_index", "cache_missing",
    ]
    write_csv(split_manifest_out, split_manifest, split_manifest_fields)
    write_csv(trainable_split_manifest_out, trainable_split_manifest, split_manifest_fields)
    write_csv(split_summary_out, split_summary_rows, split_summary_fields)

    trainable_counts = Counter(row["split"] for row in trainable_split_manifest)

    summary = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "tag": tag,
        "inputs": {
            "p450_csv": str(p450_path),
            "structure_records_csv": str(structure_path),
            "best1_audit_csv": str(audit_path),
            "google_manifest_csv": str(google_path),
            "normalized_splits": args.normalized_splits,
            "base_cache": args.base_cache,
        },
        "outputs": {
            "p450_source_copy": str(source_copy),
            "p450389_overlap_structure_records": str(structure_out),
            "p450389_overlap_best1_audited": str(all_best1_out),
            "fe_overlay_target_candidates": str(target_out),
            "p450389_uniprot_summary": str(uniprot_out),
            "missing_from_esibank_structure_records": str(missing_structure_out),
            "missing_from_fe_overlay_targets": str(missing_target_out),
            "p450389_split_sample_manifest": str(split_manifest_out),
            "p450389_trainable_split_sample_manifest": str(trainable_split_manifest_out),
            "p450389_split_summary": str(split_summary_out),
        },
        "counts": {
            "p450_source_rows": len(p450_rows),
            "p450_unique_uniprot": len(p450_ids),
            "p450_duplicate_uniprot_count": len(duplicate_p450_ids),
            "esibank_structure_records_total": len(structure_rows),
            "p450_overlap_structure_records": len(structure_overlap),
            "p450_with_esibank_structure_records": sum(1 for s in structure_stats.values() if s["rows"] > 0),
            "p450_esibank_samples": sum(s["samples"] for s in structure_stats.values()),
            "p450_esibank_positive_samples": sum(s["positive"] for s in structure_stats.values()),
            "p450_overlap_best1_audited_rows": len(all_p450_best1),
            "fe_overlay_target_candidate_rows": len(target_candidates),
            "p450_with_fe_overlay_targets": sum(1 for t in target_stats.values() if t["rows"] > 0),
            "fe_overlay_target_samples": sum(t["samples"] for t in target_stats.values()),
            "fe_overlay_target_positive_samples": sum(t["positive"] for t in target_stats.values()),
            "p450_missing_from_esibank_structure_records": len(missing_structure),
            "p450_missing_from_fe_overlay_targets": len(missing_target),
            "split_manifest_rows": len(split_manifest),
            "split_manifest_train": split_counts["train"],
            "split_manifest_val": split_counts["val"],
            "split_manifest_test": split_counts["test"],
            "split_manifest_cache_missing": sum(split_cache_missing.values()),
            "trainable_split_manifest_rows": len(trainable_split_manifest),
            "trainable_split_manifest_train": trainable_counts["train"],
            "trainable_split_manifest_val": trainable_counts["val"],
            "trainable_split_manifest_test": trainable_counts["test"],
        },
        "duplicate_p450_ids": sorted(set(duplicate_p450_ids)),
        "interpretation": (
            "target candidates are P450 rows that have ESIBank structure records, "
            "a usable audited RCSB HEM/Fe target, and an official ESIBank PDB file."
        ),
    }
    summary_out = out_dir / f"exp003_full_p450_manifest_summary_{tag}.json"
    write_json(summary_out, summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
