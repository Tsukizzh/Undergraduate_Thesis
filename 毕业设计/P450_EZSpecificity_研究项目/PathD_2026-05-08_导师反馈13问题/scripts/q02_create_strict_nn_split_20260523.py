#!/usr/bin/env python3
"""Create strict nearest-neighbor sequence splits for PathD Q2.

The script follows the EXP005 strict NN60 idea and only changes the sequence
identity threshold. It keeps the actual-used data universe and the original
feature cache unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import re
import shutil
import subprocess
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import pandas as pd


SPLITS = ("train", "val", "test")
SEQ_ID_RE = re.compile(r"enzyme_index=(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build strict NN threshold split from actual-used PathD Q2 data."
    )
    parser.add_argument(
        "--root",
        default="/root/autodl-tmp/EZSpecificity/PathD/P450",
        help="PathD/P450 root on the server.",
    )
    parser.add_argument(
        "--threshold",
        type=int,
        required=True,
        help="Sequence identity threshold percentage, e.g. 40, 60, or 80.",
    )
    parser.add_argument(
        "--exp-name",
        required=True,
        help="Data experiment directory name, e.g. exp006_strict_nn40.",
    )
    parser.add_argument(
        "--candidate",
        default=None,
        help="Split candidate name. Defaults to strict_nnXX_candidate_001.",
    )
    parser.add_argument(
        "--source-samples",
        default=None,
        help="Source id60 sample CSV used for EXP004 target counts.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="MMseqs2 threads.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=800,
        help="Random seeds for greedy component assignment.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow replacing this script's own generated output directory.",
    )
    return parser.parse_args()


def now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


def write_json(path: Path, obj: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding="utf-8")


def markdown_table(rows: list[dict[str, object]], columns: list[str]) -> str:
    if not rows:
        return "_empty_"
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(row.get(col, "")) for col in columns) + " |")
    return "\n".join(lines)


def ensure_inside(root: Path, path: Path) -> Path:
    resolved = path.resolve()
    root = root.resolve()
    if resolved != root and root not in resolved.parents:
        raise RuntimeError(f"Refusing to write outside PathD/P450: {resolved}")
    return resolved


def prepare_output_dir(root: Path, out_dir: Path, force: bool) -> None:
    out_dir = ensure_inside(root, out_dir)
    if out_dir.exists():
        if not force:
            raise RuntimeError(f"Output exists; use --force only for this generated dir: {out_dir}")
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)


def fasta_header(row: pd.Series) -> str:
    uniprots = str(row.get("uniprots", "")).replace(" ", "_")
    return f"enzyme_index={int(row['enzyme_index'])}|uniprots={uniprots}"


def write_fasta(path: Path, enzymes: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for _, row in enzymes.sort_values("enzyme_index").iterrows():
            handle.write(f">{fasta_header(row)}\n")
            seq = str(row["Protein_sequence"]).strip()
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def run_mmseqs(
    fasta_query: Path,
    fasta_target: Path,
    out_path: Path,
    log_path: Path,
    tmp_dir: Path,
    threshold: int,
    threads: int,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs",
        "easy-search",
        str(fasta_query),
        str(fasta_target),
        str(out_path),
        str(tmp_dir),
        "--min-seq-id",
        f"{threshold / 100:.3f}",
        "-c",
        "0.8",
        "--cov-mode",
        "0",
        "--format-output",
        "query,target,pident,alnlen,qlen,tlen,qcov,tcov,evalue,bits",
        "--threads",
        str(threads),
    ]
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"MMseqs2 failed with exit code {proc.returncode}: {log_path}")


def enzyme_id_from_header(value: str) -> int:
    match = SEQ_ID_RE.search(str(value))
    if not match:
        raise ValueError(f"Cannot parse enzyme index from MMseqs2 id: {value}")
    return int(match.group(1))


class DSU:
    def __init__(self, items: list[int]) -> None:
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in items}

    def find(self, item: int) -> int:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, a: int, b: int) -> None:
        ra = self.find(a)
        rb = self.find(b)
        if ra == rb:
            return
        if self.size[ra] < self.size[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        self.size[ra] += self.size[rb]


def read_edges(m8_path: Path, threshold: int) -> tuple[list[tuple[int, int]], pd.DataFrame]:
    names = ["query", "target", "pident", "alnlen", "qlen", "tlen", "qcov", "tcov", "evalue", "bits"]
    if not m8_path.exists() or m8_path.stat().st_size == 0:
        empty = pd.DataFrame(columns=names + ["query_enzyme_index", "target_enzyme_index"])
        return [], empty
    hits = pd.read_csv(m8_path, sep="\t", names=names)
    hits["query_enzyme_index"] = hits["query"].map(enzyme_id_from_header)
    hits["target_enzyme_index"] = hits["target"].map(enzyme_id_from_header)
    hits = hits[(hits["pident"] >= threshold) & (hits["qcov"] >= 0.8) & (hits["tcov"] >= 0.8)].copy()
    edges_set: set[tuple[int, int]] = set()
    for row in hits.itertuples(index=False):
        a = int(row.query_enzyme_index)
        b = int(row.target_enzyme_index)
        if a == b:
            continue
        edges_set.add((a, b) if a < b else (b, a))
    return sorted(edges_set), hits


def component_stats(
    enzymes: pd.DataFrame,
    source_samples: pd.DataFrame,
    edges: list[tuple[int, int]],
    tag: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    enzyme_ids = sorted(int(x) for x in enzymes["enzyme_index"].unique())
    dsu = DSU(enzyme_ids)
    for a, b in edges:
        dsu.union(a, b)

    groups: dict[int, list[int]] = defaultdict(list)
    for enzyme_id in enzyme_ids:
        groups[dsu.find(enzyme_id)].append(enzyme_id)

    ordered_groups = sorted(groups.values(), key=lambda xs: (-len(xs), min(xs)))
    component_by_enzyme: dict[int, str] = {}
    component_rows = []
    sample_counts = source_samples.groupby("enzyme_index").agg(
        n_samples=("label", "size"),
        n_positive=("label", lambda s: int((s == 1).sum())),
        n_negative=("label", lambda s: int((s == 0).sum())),
    )
    for idx, enzyme_list in enumerate(ordered_groups, start=1):
        comp_id = f"{tag}_comp_{idx:04d}"
        for enzyme_id in enzyme_list:
            component_by_enzyme[enzyme_id] = comp_id
        sub = sample_counts.loc[enzyme_list]
        component_rows.append(
            {
                "component_id": comp_id,
                "n_enzymes": len(enzyme_list),
                "n_samples": int(sub["n_samples"].sum()),
                "n_positive": int(sub["n_positive"].sum()),
                "n_negative": int(sub["n_negative"].sum()),
                "positive_rate": round(float(sub["n_positive"].sum() / max(1, sub["n_samples"].sum())), 6),
                "enzyme_indices": ";".join(str(x) for x in enzyme_list),
            }
        )
    enzyme_components = pd.DataFrame(
        [{"enzyme_index": enzyme_id, "component_id": comp_id} for enzyme_id, comp_id in component_by_enzyme.items()]
    )
    components = pd.DataFrame(component_rows).sort_values(
        ["n_samples", "n_positive", "n_enzymes"], ascending=False
    )
    return components, enzyme_components


def target_counts_from_source(source_samples: pd.DataFrame) -> dict[str, dict[str, int]]:
    if "assigned_split" not in source_samples.columns:
        raise RuntimeError("Source samples must contain EXP004 assigned_split column.")
    targets: dict[str, dict[str, int]] = {}
    for split in SPLITS:
        rows = source_samples[source_samples["assigned_split"] == split]
        targets[split] = {
            "samples": int(len(rows)),
            "pos": int((rows["label"] == 1).sum()),
            "enzymes": int(rows["enzyme_index"].nunique()),
        }
    return targets


def assignment_score(stats: dict[str, dict[str, int]], targets: dict[str, dict[str, int]]) -> float:
    score = 0.0
    for split in SPLITS:
        score += abs(stats[split]["n_samples"] - targets[split]["samples"]) / max(1, targets[split]["samples"])
        score += 1.5 * abs(stats[split]["n_positive"] - targets[split]["pos"]) / max(1, targets[split]["pos"])
        score += 0.3 * abs(stats[split]["n_enzymes"] - targets[split]["enzymes"]) / max(1, targets[split]["enzymes"])
    return float(score)


def choose_assignment(
    components: pd.DataFrame,
    targets: dict[str, dict[str, int]],
    seeds: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    records = components.to_dict("records")
    best_assignment: pd.DataFrame | None = None
    best_meta: dict[str, object] | None = None
    for seed in range(seeds):
        rng = random.Random(seed)
        shuffled = [dict(x) for x in records]
        rng.shuffle(shuffled)
        shuffled.sort(key=lambda r: (r["n_samples"], r["n_positive"], r["n_enzymes"]), reverse=True)
        stats = {
            split: {"n_components": 0, "n_enzymes": 0, "n_samples": 0, "n_positive": 0, "n_negative": 0}
            for split in SPLITS
        }
        assigned = []
        for comp in shuffled:
            best_split = None
            best_score = math.inf
            for split in SPLITS:
                projected = {name: dict(values) for name, values in stats.items()}
                projected[split]["n_components"] += 1
                projected[split]["n_enzymes"] += int(comp["n_enzymes"])
                projected[split]["n_samples"] += int(comp["n_samples"])
                projected[split]["n_positive"] += int(comp["n_positive"])
                projected[split]["n_negative"] += int(comp["n_negative"])
                score = assignment_score(projected, targets)
                if score < best_score:
                    best_score = score
                    best_split = split
            assert best_split is not None
            stats[best_split]["n_components"] += 1
            stats[best_split]["n_enzymes"] += int(comp["n_enzymes"])
            stats[best_split]["n_samples"] += int(comp["n_samples"])
            stats[best_split]["n_positive"] += int(comp["n_positive"])
            stats[best_split]["n_negative"] += int(comp["n_negative"])
            row = dict(comp)
            row["assigned_split"] = best_split
            assigned.append(row)
        score = assignment_score(stats, targets)
        if best_meta is None or score < float(best_meta["score"]):
            best_assignment = pd.DataFrame(assigned)
            best_meta = {"seed": seed, "score": score, "stats": stats}
    assert best_assignment is not None and best_meta is not None
    return best_assignment, best_meta


def export_split_files(
    split_dir: Path,
    threshold: int,
    tag: str,
    enzymes: pd.DataFrame,
    source_samples: pd.DataFrame,
    enzyme_components: pd.DataFrame,
    component_assignment: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    comp_split = component_assignment.set_index("component_id")["assigned_split"].to_dict()
    enzyme_meta = enzymes.merge(enzyme_components, on="enzyme_index", how="left")
    enzyme_meta[f"strict_{tag}_split"] = enzyme_meta["component_id"].map(comp_split)
    if enzyme_meta[f"strict_{tag}_split"].isna().any():
        raise RuntimeError("Some enzymes did not receive a split.")

    per_enzyme = source_samples.groupby("enzyme_index").agg(
        n_samples=("label", "size"),
        n_positive=("label", lambda s: int((s == 1).sum())),
        n_negative=("label", lambda s: int((s == 0).sum())),
        old_exp004_split=("assigned_split", lambda s: ";".join(sorted(set(map(str, s))))),
        old_id60_cluster_ids=("cluster_id", lambda s: ";".join(sorted(set(map(str, s))))),
    )
    enzyme_assignment = (
        enzyme_meta[
            ["enzyme_index", "uniprots", f"strict_{tag}_split", "component_id", "sequence_length"]
        ]
        .merge(per_enzyme, left_on="enzyme_index", right_index=True, how="left")
        .sort_values("enzyme_index")
    )
    enzyme_assignment.to_csv(split_dir / "enzyme_split_assignment.csv", index=False)

    samples = source_samples.merge(
        enzyme_assignment[["enzyme_index", f"strict_{tag}_split", "component_id"]],
        on="enzyme_index",
        how="left",
    ).copy()
    samples["assigned_split"] = samples[f"strict_{tag}_split"]
    samples[f"{tag}_component_id"] = samples["component_id"]
    samples = samples.drop(columns=["component_id"])
    samples.to_csv(split_dir / f"samples_with_strict_{tag}_split.csv", index=False)

    component_out = component_assignment.rename(columns={"assigned_split": f"strict_{tag}_split"}).copy()
    component_out.to_csv(split_dir / "component_split_assignment.csv", index=False)

    for split in SPLITS:
        split_samples = samples[samples["assigned_split"] == split].copy()
        split_samples.to_csv(split_dir / f"{split}_samples_strict_{tag}.csv", index=False)
        split_enzymes = enzymes[enzymes["enzyme_index"].isin(split_samples["enzyme_index"].unique())]
        write_fasta(split_dir / f"{split}_enzymes_strict_{tag}.fasta", split_enzymes)

    summary_rows = []
    for split in SPLITS:
        rows = samples[samples["assigned_split"] == split]
        comps = component_out[component_out[f"strict_{tag}_split"] == split]
        summary_rows.append(
            {
                "split": split,
                "n_components": int(len(comps)),
                "n_enzymes": int(rows["enzyme_index"].nunique()),
                "n_samples": int(len(rows)),
                "n_positive": int((rows["label"] == 1).sum()),
                "n_negative": int((rows["label"] == 0).sum()),
                "positive_rate": round(float((rows["label"] == 1).sum() / max(1, len(rows))), 6),
            }
        )
    pd.DataFrame(summary_rows).to_csv(split_dir / "split_summary.csv", index=False)
    return samples, enzyme_assignment


def duplicate_sequence_audit(enzymes: pd.DataFrame, enzyme_assignment: pd.DataFrame, split_col: str) -> list[dict[str, object]]:
    rows = []
    split_by_enzyme = enzyme_assignment.set_index("enzyme_index")[split_col].to_dict()
    duplicated = enzymes.groupby("Protein_sequence")["enzyme_index"].apply(list)
    group_id = 0
    for indices in duplicated:
        if len(indices) <= 1:
            continue
        group_id += 1
        splits = sorted({str(split_by_enzyme[int(idx)]) for idx in indices})
        rows.append(
            {
                "duplicate_group_id": group_id,
                "enzyme_indices": ";".join(str(int(x)) for x in sorted(indices)),
                "splits": ";".join(splits),
                "n_splits": len(splits),
                "passes": len(splits) <= 1,
            }
        )
    return rows


def validate_cross_split(
    split_dir: Path,
    threshold: int,
    tag: str,
    threads: int,
) -> dict[str, object]:
    mmseqs_dir = split_dir / "mmseqs"
    audits_dir = split_dir / "audits"
    audits_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {}
    for query_split in ("test", "val"):
        out_path = mmseqs_dir / f"{query_split}_vs_train_{tag}.m8"
        run_mmseqs(
            split_dir / f"{query_split}_enzymes_strict_{tag}.fasta",
            split_dir / f"train_enzymes_strict_{tag}.fasta",
            out_path,
            mmseqs_dir / f"{query_split}_vs_train_{tag}.log",
            mmseqs_dir / f"tmp_{query_split}_vs_train_{tag}_{int(time.time())}",
            threshold,
            threads,
        )
        _, hits = read_edges(out_path, threshold)
        violations = hits[hits["query_enzyme_index"] != hits["target_enzyme_index"]].copy()
        nearest_rows = []
        if not hits.empty:
            nonself = hits[hits["query_enzyme_index"] != hits["target_enzyme_index"]].copy()
            if not nonself.empty:
                idx = nonself.groupby("query_enzyme_index")["pident"].idxmax()
                nearest_rows = nonself.loc[idx].sort_values("pident", ascending=False).to_dict("records")
        pd.DataFrame(nearest_rows).to_csv(audits_dir / f"{query_split}_vs_train_nearest.csv", index=False)
        violations.to_csv(audits_dir / f"{query_split}_vs_train_violations.csv", index=False)
        results[f"{query_split}_vs_train_hits"] = int(len(violations))
        results[f"{query_split}_vs_train_passes"] = int(len(violations)) == 0
    return results


def main() -> None:
    args = parse_args()
    root = Path(args.root).resolve()
    tag = f"nn{args.threshold}"
    strict_tag = f"strict_{tag}"
    candidate = args.candidate or f"{strict_tag}_candidate_001"
    data_dir = root / "data" / "q02_sequence_similarity_split" / args.exp_name
    split_dir = data_dir / "splits" / candidate
    prepare_output_dir(root, data_dir, args.force)
    (data_dir / "logs").mkdir(parents=True, exist_ok=True)
    (data_dir / "mmseqs").mkdir(parents=True, exist_ok=True)
    (data_dir / "audits").mkdir(parents=True, exist_ok=True)

    actual_dir = root / "data" / "actual_used_baseline" / "tables"
    enzymes = pd.read_csv(actual_dir / "Enzymes_actual_used.csv")
    actual_samples = pd.read_csv(actual_dir / "samples_actual_used.csv")
    source_samples_path = (
        Path(args.source_samples).resolve()
        if args.source_samples
        else root
        / "data"
        / "q02_sequence_similarity_split"
        / "exp002_actual_used_cache_valid"
        / "splits"
        / "id60"
        / "samples_with_id60_cluster_split.csv"
    )
    source_samples = pd.read_csv(source_samples_path)
    if len(source_samples) != len(actual_samples):
        raise RuntimeError(f"Source samples size mismatch: {len(source_samples)} != {len(actual_samples)}")
    if source_samples["enzyme_index"].nunique() != 1479:
        raise RuntimeError("Source samples do not cover 1479 actual-used enzymes.")

    inputs_dir = data_dir / "inputs"
    all_fasta = inputs_dir / "enzymes_actual_used.fasta"
    write_fasta(all_fasta, enzymes)
    for split in SPLITS:
        rows = source_samples[source_samples["assigned_split"] == split]
        write_fasta(inputs_dir / f"{split}_enzymes_exp004_id60.fasta", enzymes[enzymes["enzyme_index"].isin(rows["enzyme_index"].unique())])
    write_json(
        data_dir / "manifest_inputs.json",
        {
            "created_at": now(),
            "experiment": args.exp_name,
            "threshold": args.threshold,
            "tag": tag,
            "actual_used_enzymes": int(enzymes["enzyme_index"].nunique()),
            "actual_used_samples": int(len(actual_samples)),
            "source_samples": str(source_samples_path),
            "source_role": "EXP004 id60 target-count reference, same as EXP005 strict NN60 flow.",
        },
    )

    all_vs_all = data_dir / "mmseqs" / f"all_vs_all_{tag}.m8"
    run_mmseqs(
        all_fasta,
        all_fasta,
        all_vs_all,
        data_dir / "mmseqs" / f"all_vs_all_{tag}.log",
        data_dir / "mmseqs" / f"tmp_all_vs_all_{tag}_{int(time.time())}",
        args.threshold,
        args.threads,
    )
    edges, hits = read_edges(all_vs_all, args.threshold)
    pd.DataFrame(edges, columns=["enzyme_a", "enzyme_b"]).to_csv(data_dir / "audits" / f"conflict_edges_{tag}.csv", index=False)

    components, enzyme_components = component_stats(enzymes, source_samples, edges, tag)
    components.to_csv(data_dir / "audits" / f"conflict_components_{tag}.csv", index=False)
    targets = target_counts_from_source(source_samples)
    assignment, assignment_meta = choose_assignment(components, targets, args.seeds)

    split_dir.mkdir(parents=True, exist_ok=True)
    samples, enzyme_assignment = export_split_files(
        split_dir, args.threshold, tag, enzymes, source_samples, enzyme_components, assignment
    )

    duplicate_rows = duplicate_sequence_audit(enzymes, enzyme_assignment, f"strict_{tag}_split")
    write_json(split_dir / "duplicate_sequence_split_audit.json", duplicate_rows)
    duplicate_cross = sum(1 for row in duplicate_rows if int(row["n_splits"]) > 1)

    cross_edges = {}
    for left, right in (("train", "val"), ("train", "test"), ("val", "test")):
        left_comps = set(assignment.loc[assignment["assigned_split"] == left, "component_id"])
        right_comps = set(assignment.loc[assignment["assigned_split"] == right, "component_id"])
        cross_edges[f"{left}_{right}"] = int(len(left_comps & right_comps))

    validation = validate_cross_split(split_dir, args.threshold, tag, args.threads)

    split_summary = {}
    for split in SPLITS:
        rows = samples[samples["assigned_split"] == split]
        comps = assignment[assignment["assigned_split"] == split]
        split_summary[split] = {
            "n_components": int(len(comps)),
            "n_enzymes": int(rows["enzyme_index"].nunique()),
            "n_samples": int(len(rows)),
            "n_positive": int((rows["label"] == 1).sum()),
            "n_negative": int((rows["label"] == 0).sum()),
            "positive_rate": round(float((rows["label"] == 1).sum() / max(1, len(rows))), 6),
            "target_samples": targets[split]["samples"],
            "target_positive": targets[split]["pos"],
            "target_enzymes": targets[split]["enzymes"],
            "sample_delta_vs_target": int(len(rows) - targets[split]["samples"]),
            "positive_delta_vs_target": int((rows["label"] == 1).sum() - targets[split]["pos"]),
            "enzyme_delta_vs_target": int(rows["enzyme_index"].nunique() - targets[split]["enzymes"]),
        }

    all_passes = (
        int(len(samples)) == 44090
        and int(samples["enzyme_index"].nunique()) == 1479
        and duplicate_cross == 0
        and all(bool(validation[f"{split}_vs_train_passes"]) for split in ("test", "val"))
    )
    summary = {
        "created_at": now(),
        "candidate": candidate,
        "source_samples": str(source_samples_path),
        "strict_rule": (
            f"No MMseqs2 pident>={args.threshold} hit passing -c 0.8 --cov-mode 0 may cross "
            f"train/val/test because assignment is by connected component of the NN{args.threshold} conflict graph."
        ),
        "target_counts": targets,
        "score": assignment_meta["score"],
        "n_enzymes": int(enzymes["enzyme_index"].nunique()),
        "n_samples": int(len(samples)),
        "n_components": int(len(components)),
        "n_conflict_edges_undirected": int(len(edges)),
        "mmseqs_hits_after_filter": int(len(hits)),
        "largest_components": assignment.head(10).to_dict("records"),
        "split_summary": split_summary,
        "cross_split_conflict_components": cross_edges,
        "duplicate_sequence_groups": len(duplicate_rows),
        "duplicate_sequence_groups_cross_split": duplicate_cross,
        f"passes_graph_{tag}_no_cross_split": all(v == 0 for v in cross_edges.values()),
        "validation": validation,
        "all_passes": all_passes,
        "files": {
            "component_split_assignment": str(split_dir / "component_split_assignment.csv"),
            "enzyme_split_assignment": str(split_dir / "enzyme_split_assignment.csv"),
            f"samples_with_strict_{tag}_split": str(split_dir / f"samples_with_strict_{tag}_split.csv"),
            "split_summary_csv": str(split_dir / "split_summary.csv"),
        },
    }
    write_json(split_dir / "split_summary.json", summary)
    write_json(split_dir / "audits" / f"strict_{tag}_validation.json", summary)

    split_rows = [
        {"split": split, **values}
        for split, values in split_summary.items()
    ]
    report = f"""# {args.exp_name} {strict_tag} split validation

创建时间：{summary['created_at']}

## 范围

- actual-used enzyme：{summary['n_enzymes']}
- actual-used sample：{summary['n_samples']}
- 序列相似度阈值：{args.threshold}%
- 覆盖度要求：80%，`--cov-mode 0`

## split 分布

{markdown_table(split_rows, ['split', 'n_components', 'n_enzymes', 'n_samples', 'n_positive', 'n_negative', 'positive_rate', 'sample_delta_vs_target', 'positive_delta_vs_target', 'enzyme_delta_vs_target'])}

## 审计结果

- conflict graph component 数：{summary['n_components']}
- 无向 conflict edge 数：{summary['n_conflict_edges_undirected']}
- test 到 train 的 >= {args.threshold}% 命中数：{validation['test_vs_train_hits']}
- val 到 train 的 >= {args.threshold}% 命中数：{validation['val_vs_train_hits']}
- 精确重复序列跨 split 组数：{duplicate_cross}
- 总体通过：{all_passes}

## 文件

- split summary：`{split_dir / 'split_summary.json'}`
- samples：`{split_dir / f'samples_with_strict_{tag}_split.csv'}`
- train/val/test CSV：`{split_dir}`
"""
    (split_dir / "README.md").write_text(report, encoding="utf-8")
    (split_dir / "audits" / f"strict_{tag}_validation.md").write_text(report, encoding="utf-8")

    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
