#!/usr/bin/env python3
"""Run row-level Q10 inference with the Q2 random best checkpoint."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import torch
import yaml
from easydict import EasyDict
from torch_geometric.loader import DataLoader


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    tmp_path.replace(path)


def load_config(path: Path) -> EasyDict:
    with path.open("r", encoding="utf-8") as f:
        return EasyDict(yaml.safe_load(f))


def choose_device(requested: str) -> torch.device:
    if requested == "cpu":
        return torch.device("cpu")
    if requested.startswith("cuda"):
        return torch.device(requested)
    return torch.device("cuda:0" if torch.cuda.is_available() else "cpu")


def load_manifest(cache_dir: Path, split: str) -> Dict[int, Dict[str, str]]:
    manifest_path = cache_dir / "manifests" / f"q10_model_cache_manifest_{split}_msa_m1r1_v1.csv"
    rows = read_csv(manifest_path)
    return {int(row["sample_id"]): row for row in rows}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--exp008-scripts", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--split", default="test")
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda:0 ...")
    parser.add_argument("--edge-mode", default="fixed", choices=["fixed", "legacy_bug"])
    parser.add_argument("--use-gdtable", action="store_true")
    parser.add_argument("--dense-dtype", default="fp16", choices=["fp16", "fp32"])
    args = parser.parse_args()

    exp008_scripts = Path(args.exp008_scripts).resolve()
    exp008_src = exp008_scripts.parent / "src"
    if exp008_src.exists():
        sys.path.insert(0, str(exp008_src))
    sys.path.insert(0, str(exp008_scripts))

    from Models.ss import SS
    from pt_dataset import PtCacheDataset
    from pt_dataset_gdtable import GraphOnlyPtCacheDataset
    from main_training_gdtable import GdTableSS

    cache_dir = Path(args.cache_dir).resolve()
    checkpoint = Path(args.checkpoint).resolve()
    config = load_config(Path(args.config))
    device = choose_device(args.device)
    use_gdtable = bool(args.use_gdtable and device.type == "cuda")

    dataset_cls = GraphOnlyPtCacheDataset if use_gdtable else PtCacheDataset
    extra = {"graph_fp16": False} if use_gdtable else {}
    dataset = dataset_cls(
        cache_dir=cache_dir,
        split=args.split,
        edge_mode=args.edge_mode,
        dist_noise=False,
        cutoff=config.transform.cutoff,
        num_r_gaussian=config.transform.num_r_gaussian,
        max_enzyme_len=config.data.max_enzyme_length,
        max_substrate_len=config.data.max_substrate_length,
        preload=False,
        **extra,
    )
    loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers,
        follow_batch=["ligand_index"],
    )

    ckpt = torch.load(checkpoint, map_location="cpu", weights_only=False)
    if use_gdtable:
        model = GdTableSS(config, cache_dir=str(cache_dir), edge_mode=args.edge_mode, dense_dtype=args.dense_dtype)
        model.gdtable.activate(args.split, device, rank=0)
    else:
        model = SS(config)
    load_result = model.load_state_dict(ckpt.get("state_dict", ckpt), strict=False)
    model.to(device).eval()

    sample_ids = [int(x) for x in dataset.valid_idx]
    manifest_by_sample = load_manifest(cache_dir, args.split)
    rows: List[Dict[str, object]] = []
    offset = 0
    with torch.inference_mode():
        for batch in loader:
            n_graphs = int(batch.num_graphs)
            batch_sample_ids = sample_ids[offset: offset + n_graphs]
            offset += n_graphs
            if use_gdtable:
                model.gdtable.add_row_indices(batch)
            batch = batch.to(device, non_blocking=True)
            if use_gdtable:
                model.gdtable.attach(batch)
            output = model(batch)
            logits = output[0] if isinstance(output, tuple) else output
            logits = logits.squeeze(-1).detach().cpu()
            scores = torch.sigmoid(logits)
            for sid, logit, score in zip(batch_sample_ids, logits.tolist(), scores.tolist()):
                meta = manifest_by_sample[int(sid)]
                rows.append({
                    "sample_id": int(sid),
                    "candidate_id": meta["candidate_id"],
                    "input_list": meta["input_list"],
                    "enzyme_id": meta["enzyme_id"],
                    "substrate_id": meta["substrate_id"],
                    "model_logit": float(logit),
                    "model_score_sigmoid": float(score),
                    "docking_score": meta.get("docking_score", ""),
                    "box_source": meta.get("box_source", ""),
                    "heme_fe_quality": meta.get("heme_fe_quality", ""),
                    "nearest_cys_sg_fe_distance": meta.get("nearest_cys_sg_fe_distance", ""),
                    "template_pdb_id": meta.get("template_pdb_id", ""),
                    "fit_rmsd": meta.get("fit_rmsd", ""),
                    "pocket_method": meta.get("pocket_method", ""),
                    "n_ligand_atoms": meta.get("n_ligand_atoms", ""),
                    "n_protein_atoms_pocket": meta.get("n_protein_atoms_pocket", ""),
                    "pose_sdf_path": meta.get("pose_sdf_path", ""),
                    "receptor_pdb_path": meta.get("receptor_pdb_path", ""),
                })

    rows.sort(key=lambda r: float(r["model_score_sigmoid"]), reverse=True)
    for i, row in enumerate(rows, start=1):
        row["rank_overall"] = i

    by_list: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_list[str(row["input_list"])].append(row)
    for group_rows in by_list.values():
        group_rows.sort(key=lambda r: float(r["model_score_sigmoid"]), reverse=True)
        for i, row in enumerate(group_rows, start=1):
            row["rank_within_input_list"] = i

    out_dir = Path(args.output_dir)
    fieldnames = [
        "rank_overall", "rank_within_input_list", "candidate_id", "input_list",
        "sample_id", "enzyme_id", "substrate_id", "model_logit", "model_score_sigmoid",
        "docking_score", "box_source", "heme_fe_quality", "nearest_cys_sg_fe_distance",
        "template_pdb_id", "fit_rmsd", "pocket_method", "n_ligand_atoms",
        "n_protein_atoms_pocket", "pose_sdf_path", "receptor_pdb_path",
    ]
    all_path = out_dir / f"q10_model_scores_all_{args.split}_msa_m1r1_v1.csv"
    write_csv(all_path, rows, fieldnames)
    split_paths = {}
    for input_list, group_rows in sorted(by_list.items()):
        safe = input_list.replace("/", "_")
        p = out_dir / f"q10_model_scores_{safe}_{args.split}_msa_m1r1_v1.csv"
        write_csv(p, group_rows, fieldnames)
        split_paths[input_list] = str(p)

    score_values = [float(r["model_score_sigmoid"]) for r in rows]
    summary = {
        "cache_dir": str(cache_dir),
        "checkpoint": str(checkpoint),
        "config": str(Path(args.config).resolve()),
        "device": str(device),
        "use_gdtable": use_gdtable,
        "n_samples": len(rows),
        "by_input_list": {k: len(v) for k, v in sorted(by_list.items())},
        "score_min": min(score_values) if score_values else math.nan,
        "score_max": max(score_values) if score_values else math.nan,
        "score_mean": sum(score_values) / len(score_values) if score_values else math.nan,
        "load_state_dict_missing_keys": list(load_result.missing_keys),
        "load_state_dict_unexpected_keys": list(load_result.unexpected_keys),
        "load_state_dict_missing_count": len(load_result.missing_keys),
        "load_state_dict_unexpected_count": len(load_result.unexpected_keys),
        "all_scores_csv": str(all_path),
        "split_score_csvs": split_paths,
        "note": "Q10 has no labels here; model_score_sigmoid is a ranking score, not validated wet-lab probability.",
    }
    summary_path = out_dir / f"q10_model_scores_summary_{args.split}_msa_m1r1_v1.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
