#!/usr/bin/env python3
"""Create EXP009 strict NN60 split candidates for PathD Q2.

EXP009 does not train. It generates several strict NN60 candidate splits with
different target ratios and seed choices, then audits each candidate.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import shutil
import sys
from pathlib import Path
from typing import Any

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from q02_create_strict_nn_split_20260523 import (  # noqa: E402
    SPLITS,
    component_stats,
    duplicate_sequence_audit,
    export_split_files,
    markdown_table,
    now,
    read_edges,
    validate_cross_split,
    write_fasta,
    write_json,
)


RATIO_PLANS = {
    "ratio211": {
        "label": "2:1:1",
        "fractions": {"train": 0.50, "val": 0.25, "test": 0.25},
        "keep_top": 3,
    },
    "ratio71515": {
        "label": "7:1.5:1.5",
        "fractions": {"train": 0.70, "val": 0.15, "test": 0.15},
        "keep_top": 2,
    },
    "ratio811": {
        "label": "8:1:1",
        "fractions": {"train": 0.80, "val": 0.10, "test": 0.10},
        "keep_top": 2,
    },
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--root",
        default="/root/autodl-tmp/EZSpecificity/PathD/P450",
        help="PathD/P450 root on the server.",
    )
    p.add_argument("--exp-name", default="exp009_strict_nn60_candidates")
    p.add_argument("--threshold", type=int, default=60)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--seeds", type=int, default=1200, help="Number of greedy random seeds per ratio.")
    p.add_argument("--force", action="store_true", help="Replace this script's generated EXP009 data dir.")
    p.add_argument(
        "--existing-mmseqs",
        default=None,
        help="Existing all_vs_all_nn60.m8. Defaults to EXP005 strict NN60 file.",
    )
    return p.parse_args()


def ensure_new_dir(path: Path, force: bool, root: Path) -> None:
    path = path.resolve()
    root = root.resolve()
    if path != root and root not in path.parents:
        raise RuntimeError(f"Refusing to write outside root: {path}")
    if path.exists():
        if not force:
            raise RuntimeError(f"Output exists; use --force only for this generated dir: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True)


def balanced_targets(samples: pd.DataFrame, enzymes: pd.DataFrame, fractions: dict[str, float]) -> dict[str, dict[str, int]]:
    total_samples = len(samples)
    total_pos = int((samples["label"] == 1).sum())
    total_enzymes = int(enzymes["enzyme_index"].nunique())
    targets: dict[str, dict[str, int]] = {}
    splits = list(SPLITS)
    for split in splits[:-1]:
        frac = fractions[split]
        targets[split] = {
            "samples": int(round(total_samples * frac)),
            "pos": int(round(total_pos * frac)),
            "enzymes": int(round(total_enzymes * frac)),
        }
    last = splits[-1]
    targets[last] = {
        "samples": total_samples - sum(targets[s]["samples"] for s in splits[:-1]),
        "pos": total_pos - sum(targets[s]["pos"] for s in splits[:-1]),
        "enzymes": total_enzymes - sum(targets[s]["enzymes"] for s in splits[:-1]),
    }
    return targets


def score_stats(stats: dict[str, dict[str, int]], targets: dict[str, dict[str, int]]) -> float:
    score = 0.0
    for split in SPLITS:
        target_samples = max(1, targets[split]["samples"])
        target_pos = max(1, targets[split]["pos"])
        target_enzymes = max(1, targets[split]["enzymes"])
        score += abs(stats[split]["n_samples"] - targets[split]["samples"]) / target_samples
        score += 1.6 * abs(stats[split]["n_positive"] - targets[split]["pos"]) / target_pos
        score += 0.25 * abs(stats[split]["n_enzymes"] - targets[split]["enzymes"]) / target_enzymes
    return float(score)


def assign_for_seed(components: pd.DataFrame, targets: dict[str, dict[str, int]], seed: int) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = random.Random(seed)
    records = [dict(x) for x in components.to_dict("records")]
    rng.shuffle(records)
    records.sort(key=lambda r: (r["n_samples"], r["n_positive"], r["n_enzymes"]), reverse=True)
    stats = {
        split: {"n_components": 0, "n_enzymes": 0, "n_samples": 0, "n_positive": 0, "n_negative": 0}
        for split in SPLITS
    }
    assigned = []
    for comp in records:
        best_split = None
        best_score = math.inf
        for split in SPLITS:
            projected = {name: dict(values) for name, values in stats.items()}
            projected[split]["n_components"] += 1
            projected[split]["n_enzymes"] += int(comp["n_enzymes"])
            projected[split]["n_samples"] += int(comp["n_samples"])
            projected[split]["n_positive"] += int(comp["n_positive"])
            projected[split]["n_negative"] += int(comp["n_negative"])
            candidate_score = score_stats(projected, targets)
            if candidate_score < best_score:
                best_score = candidate_score
                best_split = split
        assert best_split is not None
        for key in ("n_components", "n_enzymes", "n_samples", "n_positive", "n_negative"):
            if key == "n_components":
                stats[best_split][key] += 1
            else:
                stats[best_split][key] += int(comp[key])
        row = dict(comp)
        row["assigned_split"] = best_split
        assigned.append(row)
    score = score_stats(stats, targets)
    return pd.DataFrame(assigned), {"seed": seed, "score": score, "stats": stats}


def assignment_signature(assignment: pd.DataFrame) -> tuple[tuple[str, str], ...]:
    return tuple(sorted((str(r.component_id), str(r.assigned_split)) for r in assignment.itertuples(index=False)))


def split_enzyme_sets(samples: pd.DataFrame) -> dict[str, set[int]]:
    return {split: set(map(int, samples.loc[samples["assigned_split"] == split, "enzyme_index"].unique())) for split in SPLITS}


def jaccard(a: set[int], b: set[int]) -> float:
    if not a and not b:
        return 1.0
    return round(len(a & b) / max(1, len(a | b)), 6)


def summarize_candidate(samples: pd.DataFrame, assignment: pd.DataFrame, targets: dict[str, dict[str, int]]) -> dict[str, dict[str, int | float]]:
    out: dict[str, dict[str, int | float]] = {}
    for split in SPLITS:
        rows = samples[samples["assigned_split"] == split]
        comps = assignment[assignment["assigned_split"] == split]
        n = int(len(rows))
        pos = int((rows["label"] == 1).sum())
        out[split] = {
            "n_components": int(len(comps)),
            "n_enzymes": int(rows["enzyme_index"].nunique()),
            "n_samples": n,
            "n_positive": pos,
            "n_negative": int(n - pos),
            "positive_rate": round(pos / n, 6) if n else 0.0,
            "target_samples": targets[split]["samples"],
            "target_positive": targets[split]["pos"],
            "target_enzymes": targets[split]["enzymes"],
            "sample_delta_vs_target": int(n - targets[split]["samples"]),
            "positive_delta_vs_target": int(pos - targets[split]["pos"]),
            "enzyme_delta_vs_target": int(rows["enzyme_index"].nunique() - targets[split]["enzymes"]),
        }
    return out


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    data_dir = root / "data/q02_sequence_similarity_split" / args.exp_name
    ensure_new_dir(data_dir, args.force, root)
    (data_dir / "splits").mkdir()
    (data_dir / "reports").mkdir()
    (data_dir / "audits").mkdir()
    (data_dir / "inputs").mkdir()

    actual_dir = root / "data/actual_used_baseline/tables"
    enzymes = pd.read_csv(actual_dir / "Enzymes_actual_used.csv")
    source_samples_path = (
        root
        / "data/q02_sequence_similarity_split/exp002_actual_used_cache_valid/splits/id60/samples_with_id60_cluster_split.csv"
    )
    source_samples = pd.read_csv(source_samples_path)
    if int(enzymes["enzyme_index"].nunique()) != 1479 or len(source_samples) != 44090:
        raise RuntimeError("EXP009 must use actual-used 1479 enzymes / 44090 samples.")

    all_fasta = data_dir / "inputs/enzymes_actual_used.fasta"
    write_fasta(all_fasta, enzymes)

    mmseqs_path = (
        Path(args.existing_mmseqs).resolve()
        if args.existing_mmseqs
        else root / "data/q02_sequence_similarity_split/exp005_strict_nn60/mmseqs/all_vs_all_nn60.m8"
    )
    if not mmseqs_path.is_file():
        raise FileNotFoundError(f"Missing existing NN60 all-vs-all file: {mmseqs_path}")

    edges, hits = read_edges(mmseqs_path, args.threshold)
    pd.DataFrame(edges, columns=["enzyme_a", "enzyme_b"]).to_csv(data_dir / "audits/conflict_edges_nn60.csv", index=False)
    components, enzyme_components = component_stats(enzymes, source_samples, edges, "nn60")
    components.to_csv(data_dir / "audits/conflict_components_nn60.csv", index=False)

    exp005_samples_path = (
        root
        / "data/q02_sequence_similarity_split/exp005_strict_nn60/splits/strict_nn60_candidate_001/samples_with_strict_nn60_split.csv"
    )
    exp005_sets: dict[str, set[int]] = {}
    if exp005_samples_path.is_file():
        exp005_samples = pd.read_csv(exp005_samples_path)
        exp005_sets = split_enzyme_sets(exp005_samples)

    candidate_records: list[dict[str, Any]] = []
    selected: list[dict[str, Any]] = []
    for ratio_name, plan in RATIO_PLANS.items():
        targets = balanced_targets(source_samples, enzymes, plan["fractions"])
        seen: set[tuple[tuple[str, str], ...]] = set()
        ranked: list[dict[str, Any]] = []
        for seed in range(args.seeds):
            assignment, meta = assign_for_seed(components, targets, seed)
            sig = assignment_signature(assignment)
            if sig in seen:
                continue
            seen.add(sig)
            ranked.append({"ratio_name": ratio_name, "plan": plan, "targets": targets, "assignment": assignment, "meta": meta})
        ranked.sort(key=lambda x: float(x["meta"]["score"]))
        keep_top = int(plan["keep_top"])
        for rank, item in enumerate(ranked[:keep_top], start=1):
            candidate_name = f"strict_nn60_{ratio_name}_rank{rank:02d}_seed{int(item['meta']['seed']):04d}"
            item["candidate_name"] = candidate_name
            selected.append(item)

    for item in selected:
        candidate_name = str(item["candidate_name"])
        ratio_name = str(item["ratio_name"])
        split_dir = data_dir / "splits" / candidate_name
        split_dir.mkdir(parents=True)
        samples, enzyme_assignment = export_split_files(
            split_dir,
            args.threshold,
            "nn60",
            enzymes,
            source_samples,
            enzyme_components,
            item["assignment"],
        )
        duplicate_rows = duplicate_sequence_audit(enzymes, enzyme_assignment, "strict_nn60_split")
        duplicate_cross = sum(1 for row in duplicate_rows if int(row["n_splits"]) > 1)
        write_json(split_dir / "duplicate_sequence_split_audit.json", duplicate_rows)
        validation = validate_cross_split(split_dir, args.threshold, "nn60", args.threads)
        summary = summarize_candidate(samples, item["assignment"], item["targets"])
        sample_total = sum(int(summary[s]["n_samples"]) for s in SPLITS)
        enzyme_total = int(samples["enzyme_index"].nunique())
        all_passes = (
            sample_total == 44090
            and enzyme_total == 1479
            and duplicate_cross == 0
            and validation["test_vs_train_passes"]
            and validation["val_vs_train_passes"]
        )
        sets = split_enzyme_sets(samples)
        exp005_jaccard = {
            split: jaccard(sets[split], exp005_sets.get(split, set())) if exp005_sets else None for split in SPLITS
        }
        record = {
            "created_at": now(),
            "exp_name": args.exp_name,
            "candidate": candidate_name,
            "ratio_name": ratio_name,
            "ratio_label": item["plan"]["label"],
            "target_fractions": item["plan"]["fractions"],
            "seed": int(item["meta"]["seed"]),
            "score": float(item["meta"]["score"]),
            "n_samples": sample_total,
            "n_enzymes": enzyme_total,
            "n_components": int(len(components)),
            "n_conflict_edges_undirected": int(len(edges)),
            "mmseqs_hits_after_filter": int(len(hits)),
            "split_summary": summary,
            "duplicate_sequence_groups_cross_split": duplicate_cross,
            "validation": validation,
            "exp005_split_enzyme_jaccard": exp005_jaccard,
            "all_passes": bool(all_passes),
            "files": {
                "split_dir": str(split_dir),
                "samples": str(split_dir / "samples_with_strict_nn60_split.csv"),
                "train": str(split_dir / "train_samples_strict_nn60.csv"),
                "val": str(split_dir / "val_samples_strict_nn60.csv"),
                "test": str(split_dir / "test_samples_strict_nn60.csv"),
            },
        }
        write_json(split_dir / "split_summary.json", record)
        write_json(split_dir / "audits/strict_nn60_validation.json", record)
        (split_dir / "README.md").write_text(
            f"# {candidate_name}\n\n"
            f"- ratio: {item['plan']['label']}\n"
            f"- seed: {int(item['meta']['seed'])}\n"
            f"- score: {float(item['meta']['score']):.6f}\n"
            f"- all_passes: {bool(all_passes)}\n\n"
            + markdown_table(
                [{"split": split, **summary[split]} for split in SPLITS],
                [
                    "split",
                    "n_components",
                    "n_enzymes",
                    "n_samples",
                    "n_positive",
                    "n_negative",
                    "positive_rate",
                    "sample_delta_vs_target",
                    "positive_delta_vs_target",
                    "enzyme_delta_vs_target",
                ],
            )
            + "\n",
            encoding="utf-8",
        )
        flat = {
            "candidate": candidate_name,
            "ratio_name": ratio_name,
            "ratio_label": item["plan"]["label"],
            "seed": int(item["meta"]["seed"]),
            "score": float(item["meta"]["score"]),
            "all_passes": bool(all_passes),
            "test_vs_train_hits": int(validation["test_vs_train_hits"]),
            "val_vs_train_hits": int(validation["val_vs_train_hits"]),
            "duplicate_cross": int(duplicate_cross),
            "test_exp005_jaccard": exp005_jaccard["test"],
            "val_exp005_jaccard": exp005_jaccard["val"],
        }
        for split in SPLITS:
            for key, value in summary[split].items():
                flat[f"{split}_{key}"] = value
        candidate_records.append(flat)

    ranking = pd.DataFrame(candidate_records).sort_values(["ratio_name", "score"])
    ranking.to_csv(data_dir / "reports/exp009_candidate_ranking.csv", index=False)
    write_json(data_dir / "reports/exp009_candidate_ranking.json", ranking.to_dict("records"))
    write_json(
        data_dir / "manifest_exp009_candidates.json",
        {
            "created_at": now(),
            "exp_name": args.exp_name,
            "purpose": "Strict NN60 candidate split audit before choosing one EXP010 training split.",
            "threshold": args.threshold,
            "coverage": 0.8,
            "source_samples": str(source_samples_path),
            "existing_mmseqs": str(mmseqs_path),
            "n_enzymes": 1479,
            "n_samples": 44090,
            "ratio_plans": RATIO_PLANS,
            "candidate_count": int(len(ranking)),
        },
    )

    report_rows = ranking.to_dict("records")
    report = "# EXP009 strict NN60 candidate ranking\n\n"
    report += "EXP009 only creates and audits candidate splits. It does not build pt caches or train models.\n\n"
    report += markdown_table(
        report_rows,
        [
            "candidate",
            "ratio_label",
            "seed",
            "score",
            "all_passes",
            "test_n_samples",
            "test_n_positive",
            "test_positive_rate",
            "test_exp005_jaccard",
            "val_n_samples",
            "val_n_positive",
            "train_n_samples",
            "train_n_positive",
        ],
    )
    report += "\n"
    (data_dir / "reports/exp009_candidate_ranking.md").write_text(report, encoding="utf-8")
    print(json.dumps({"data_dir": str(data_dir), "candidates": candidate_records}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
