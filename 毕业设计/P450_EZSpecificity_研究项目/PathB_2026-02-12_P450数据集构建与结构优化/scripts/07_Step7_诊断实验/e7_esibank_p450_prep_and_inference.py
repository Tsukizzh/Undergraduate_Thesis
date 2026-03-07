"""
E7 Step 0+1: ESIBank P450 Internal Benchmark — Data Preparation & Inference

Filters ESIBank brenda/random_split/testing_datas_0.csv for verified P450 enzymes
(367 domain-filtered UniProt IDs, NOT EC number filtering), then runs model inference
using the same checkpoint as E6 (random_split fold 0).

Output: esibank_p450_predictions.csv + data prep summary
"""
from __future__ import annotations

import sys
import os
import warnings
import json

import numpy as np
import pandas as pd
import torch
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score

warnings.filterwarnings("ignore")

# --- Paths ---
PROJECT_ROOT = Path(r"D:\EZSpecificity_Project")
SRC_DIR = PROJECT_ROOT / "src"
PATHB = PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目" / "PathB_2026-02-12_P450数据集构建与结构优化"
OUTPUT_DIR = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E7_ESIBank_P450_内部基准"
CHECKPOINT_DIR = PROJECT_ROOT / "saved_model" / "model" / "run_0"

# ESIBank brenda paths (Google Drive shortcut)
G_BRENDA = Path(r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\brenda")
ENZYMES_CSV = G_BRENDA / "enzymes.csv"
TESTING_CSV = G_BRENDA / "random_split" / "testing_datas_0.csv"

# P450 verified UniProt IDs (domain-filtered, 367 enzymes with training data)
P450_UNIPROT_CSV = (PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目"
                    / "PathA_2026-01-08_模型评估测试集构建" / "source_data"
                    / "02_底物数据" / "P450酶底物反应详表_汇总版.csv")

# Column mapping: brenda CSV -> model expected
COL_MAP = {
    "reaction": "Substrate Index",
    "enzyme": "Enzyme Index",
    "label": "Label",
    "structure_index": "Dock Index",
}


def step0_data_preparation():
    """Match P450 UniProt IDs to ESIBank enzyme indices, filter test fold."""
    print("=" * 70)
    print("  Step 0: Data Preparation")
    print("=" * 70)

    # 1. Load verified P450 UniProt IDs
    p450_df = pd.read_csv(P450_UNIPROT_CSV, usecols=["uniprot_id"])
    p450_uniprots = set(p450_df["uniprot_id"].str.strip().str.upper().unique())
    print(f"  Verified P450 UniProt IDs: {len(p450_uniprots)}")

    # 2. Load ESIBank enzymes.csv and match
    enzymes = pd.read_csv(ENZYMES_CSV)
    enzymes["enzyme_index"] = enzymes.index  # row number = enzyme index
    enzymes["uniprots"] = enzymes["uniprots"].str.strip().str.upper()

    matched = enzymes[enzymes["uniprots"].isin(p450_uniprots)]
    p450_enzyme_indices = set(matched["enzyme_index"].values)

    unmatched = p450_uniprots - set(matched["uniprots"].values)
    multi_match = matched["uniprots"].value_counts()
    multi_match = multi_match[multi_match > 1]

    print(f"  ESIBank total enzymes: {len(enzymes)}")
    print(f"  Matched P450 enzymes: {len(matched)} (from {len(p450_uniprots)} UniProt IDs)")
    print(f"  Unmatched UniProt IDs: {len(unmatched)}")
    if len(multi_match) > 0:
        print(f"  WARNING: {len(multi_match)} UniProt IDs matched multiple times")
    if len(unmatched) > 0 and len(unmatched) <= 10:
        print(f"  Unmatched IDs: {unmatched}")

    # 3. Filter testing_datas_0.csv for P450 rows
    test_df = pd.read_csv(TESTING_CSV)
    print(f"\n  ESIBank test fold 0 total rows: {len(test_df)}")
    print(f"  Columns: {list(test_df.columns)}")

    p450_test = test_df[test_df["enzyme"].isin(p450_enzyme_indices)].copy()
    print(f"  P450 rows in test fold: {len(p450_test)}")
    print(f"  P450 unique enzymes: {p450_test['enzyme'].nunique()}")
    print(f"  P450 unique substrates: {p450_test['reaction'].nunique()}")
    print(f"  P450 positive: {(p450_test['label'] == 1).sum()}")
    print(f"  P450 negative: {(p450_test['label'] == 0).sum()}")

    # 4. Check structure availability
    n_has_structure = (p450_test["structure_index"] >= 0).sum()
    n_no_structure = (p450_test["structure_index"] < 0).sum()
    print(f"  With structure (structure_index >= 0): {n_has_structure}")
    print(f"  Without structure (structure_index = -1): {n_no_structure}")

    # 5. Rename columns for model compatibility
    p450_renamed = p450_test.rename(columns=COL_MAP)
    for col in ["Substrate Index", "Enzyme Index", "Dock Index", "Label"]:
        if col in p450_renamed.columns:
            p450_renamed[col] = p450_renamed[col].astype(int)

    # 6. Save filtered CSV
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    prep_csv = OUTPUT_DIR / "esibank_p450_testing_renamed.csv"
    p450_renamed.to_csv(prep_csv, index=False)
    print(f"\n  Saved filtered CSV: {prep_csv}")

    # 7. Save matched enzyme mapping for later analysis
    matched_info = matched[["enzyme_index", "uniprots"]].copy()
    matched_info.to_csv(OUTPUT_DIR / "esibank_p450_enzyme_mapping.csv", index=False)

    # 8. Save prep summary
    prep_summary = {
        "n_verified_p450_uniprots": len(p450_uniprots),
        "n_matched_enzymes": len(matched),
        "n_unmatched_uniprots": len(unmatched),
        "n_test_fold_total": len(test_df),
        "n_p450_test_rows": len(p450_test),
        "n_p450_enzymes": int(p450_test["enzyme"].nunique()),
        "n_p450_substrates": int(p450_test["reaction"].nunique()),
        "n_positive": int((p450_test["label"] == 1).sum()),
        "n_negative": int((p450_test["label"] == 0).sum()),
        "n_with_structure": int(n_has_structure),
        "n_without_structure": int(n_no_structure),
        "pos_neg_ratio": f"1:{(p450_test['label'] == 0).sum() / max((p450_test['label'] == 1).sum(), 1):.1f}",
    }
    with open(OUTPUT_DIR / "step0_prep_summary.json", "w", encoding="utf-8") as f:
        json.dump(prep_summary, f, indent=2, ensure_ascii=False)

    return prep_csv, prep_summary


def _sigmoid_np(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out.astype(np.float32)


def build_config_for_brenda_p450(config_yaml, tmp_csv_path):
    """Build config pointing to brenda feature files."""
    sys.path.insert(0, str(SRC_DIR))
    from utils import load_config

    config = load_config(str(config_yaml))

    config.num_cpus = 0
    config.num_gpus = 1

    config.data.tag = "e7_esibank_p450_benchmark"
    config.data.representer = "structure_sequence"

    config.data.train_data_path = str(tmp_csv_path)
    config.data.val_data_path = str(tmp_csv_path)
    config.data.test_data_path = str(tmp_csv_path)

    # Point to brenda feature files
    config.data.enzyme_lmdb_path = str(G_BRENDA / "enzyme_features.lmdb")
    config.data.reaction_lmdb_path = str(G_BRENDA / "reaction_features.lmdb")
    config.data.grover_path = str(G_BRENDA / "grover_fingerprint.lmdb")
    config.data.morgan_path = str(G_BRENDA / "morgan_fingerprint.npy")
    config.data.structure_processed_path = str(G_BRENDA / "structure" / "structure_features.lmdb")
    config.data.sequence_processed_path = "nonexistent_seq_features.lmdb"
    config.data.high_quality_id_path = "nonexistent_hq.txt"

    config.data.full_data = False
    config.data.batch_size = 16
    config.data.sample_weight = [1.0, 1.0]
    config.data.fake_sequence_ratio = 0
    config.data.max_substrate_length = 280
    config.data.max_enzyme_length = 1450
    config.data.features = ["morgan", 1024, "grover_mean", 4885]
    config.data.atom_features = ["grover", 2400]

    return config


@torch.no_grad()
def step1_inference(prep_csv):
    """Run model inference on filtered ESIBank P450 test data."""
    print("\n" + "=" * 70)
    print("  Step 1: Model Inference")
    print("=" * 70)

    sys.path.insert(0, str(SRC_DIR))

    # Suppress RDKit warnings
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    from Datasets.brenda import Singledataset
    from Models.ss import SS

    config_yaml = CHECKPOINT_DIR / "complete-full-random-all-0-complex.yml"
    checkpoint = CHECKPOINT_DIR / "models" / "best-checkpoint.ckpt"
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    print(f"  Device: {device}")
    print(f"  Config: {config_yaml}")
    print(f"  Checkpoint: {checkpoint}")

    config = build_config_for_brenda_p450(config_yaml, prep_csv)

    dm = Singledataset(config)
    test_loader = dm.test_dataloader()

    model = SS.load_from_checkpoint(str(checkpoint), config=config, map_location="cpu")
    model.eval()
    model.to(device)

    logits_chunks = []
    n_processed = 0
    n_total = len(dm.test_prediction_df)
    print(f"  Valid samples after feature filtering: {n_total}")

    for batch_idx, batch in enumerate(test_loader):
        batch = batch.to(device)
        logits, _tags = model(batch)
        logits_np = logits.squeeze(-1).detach().float().cpu().numpy().ravel()
        logits_chunks.append(logits_np)
        n_processed += len(logits_np)
        if (batch_idx + 1) % 20 == 0:
            print(f"    Batch {batch_idx+1}: {n_processed}/{n_total}")

    logits_all = np.concatenate(logits_chunks).astype(np.float32) if logits_chunks else np.zeros(0, dtype=np.float32)
    pred_df = dm.test_prediction_df.copy()
    assert len(pred_df) == len(logits_all), f"Mismatch: df={len(pred_df)}, logits={len(logits_all)}"

    pred_df["logit"] = logits_all
    pred_df["prob"] = _sigmoid_np(logits_all)

    # Save predictions
    pred_csv = OUTPUT_DIR / "esibank_p450_predictions.csv"
    pred_df.to_csv(pred_csv, index=False)

    # Quick metrics
    if len(pred_df) > 0 and pred_df["Label"].nunique() == 2:
        auc_roc = roc_auc_score(pred_df["Label"], pred_df["logit"])
        auc_pr = average_precision_score(pred_df["Label"], pred_df["logit"])
    else:
        auc_roc = auc_pr = float("nan")

    n_pos = int(pred_df["Label"].sum())
    n_neg = len(pred_df) - n_pos

    print(f"\n  --- Inference Results ---")
    print(f"  Valid predictions: {len(pred_df)}")
    print(f"  Positive: {n_pos}, Negative: {n_neg}")
    print(f"  AUC-ROC: {auc_roc:.4f}")
    print(f"  AUC-PR: {auc_pr:.4f}")
    print(f"  Logit mean: {pred_df['logit'].mean():.4f}, std: {pred_df['logit'].std():.4f}")
    print(f"  Saved: {pred_csv}")

    # Save inference summary
    inference_summary = {
        "n_input_csv_rows": pd.read_csv(prep_csv).shape[0],
        "n_valid_predictions": len(pred_df),
        "n_dropped_by_feature_filter": pd.read_csv(prep_csv).shape[0] - len(pred_df),
        "n_positive": n_pos,
        "n_negative": n_neg,
        "auc_roc": float(auc_roc),
        "auc_pr": float(auc_pr),
        "logit_mean": float(pred_df["logit"].mean()),
        "logit_std": float(pred_df["logit"].std()),
        "logit_pos_mean": float(pred_df[pred_df["Label"] == 1]["logit"].mean()) if n_pos > 0 else None,
        "logit_neg_mean": float(pred_df[pred_df["Label"] == 0]["logit"].mean()) if n_neg > 0 else None,
        "score_gap": float(pred_df[pred_df["Label"] == 1]["logit"].mean() - pred_df[pred_df["Label"] == 0]["logit"].mean()) if n_pos > 0 and n_neg > 0 else None,
    }
    with open(OUTPUT_DIR / "step1_inference_summary.json", "w", encoding="utf-8") as f:
        json.dump(inference_summary, f, indent=2, ensure_ascii=False)

    return pred_csv, inference_summary


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Step 0: Data Preparation
    prep_csv, prep_summary = step0_data_preparation()

    # Step 1: Inference
    pred_csv, inf_summary = step1_inference(prep_csv)

    print("\n" + "=" * 70)
    print("  DONE: Prep + Inference Complete")
    print("=" * 70)
    print(f"  Prep summary: {OUTPUT_DIR / 'step0_prep_summary.json'}")
    print(f"  Predictions: {pred_csv}")
    print(f"  Inference summary: {OUTPUT_DIR / 'step1_inference_summary.json'}")
    print(f"\n  Next: Run e7_esibank_p450_analysis.py for E1'-E6'")


if __name__ == "__main__":
    main()
