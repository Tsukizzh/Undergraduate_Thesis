"""
E8 Channel Ablation: zero-out input channels at inference to measure contribution.

Channels:
    SEQ  - Enzyme sequence (ESM embeddings -> protein_mlp)
    MOL  - Molecular features (grover atom + grover_mean + morgan)
    STR  - Structure features (EGNN: x_str atom + str_mean global)

Conditions (8):
    baseline, seq_off, mol_off, str_off, seq_only, mol_only, str_only, all_off

Method: zero hidden tensors AFTER encoder projection (not raw inputs).
When one side of cross-attention is completely off, zero the attention outputs
explicitly rather than letting bias-influenced attention run.

Families: P450 (floor, AUC~0.52) + Phosphatase (healthy control, AUC~0.88)
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    roc_auc_score,
)
from rdkit import RDLogger

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(r"D:\EZSpecificity_Project")
SRC_DIR = PROJECT_ROOT / "src"
PATHB = PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目" / "PathB_2026-02-12_P450数据集构建与结构优化"
CHECKPOINT_DIR = PROJECT_ROOT / "saved_model" / "model" / "run_0"
OUTPUT_DIR = PATHB / "results" / "08_Step8_补充诊断实验" / "E8_通道关闭"

# P450
P450_CSV = PATHB / "data" / "00_shared" / "datasets" / "B6_v1" / "data.csv"
P450_SHARED = PATHB / "data" / "00_shared" / "features"
P450_STRUCTURE = PATHB / "data" / "02_Step2_因子实验" / "EXP01_B6_10A_noHeme" / "structure_features.lmdb"
P450_HQ_IDS = PATHB / "data" / "02_Step2_因子实验" / "EXP01_B6_10A_noHeme" / "high_quality_id.txt"

# Phosphatase (control)
PHOSPHATASE_ROOT = Path(r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\small_family\Phosphatase")
PHOSPHATASE_CSV = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E6_expansion_多家族推理" / "Phosphatase_test_renamed.csv"

# ── Ablation conditions ────────────────────────────────────────────────────
ABLATION_CONDITIONS = {
    "baseline":  {"zero_seq": False, "zero_mol": False, "zero_str": False},
    "seq_off":   {"zero_seq": True,  "zero_mol": False, "zero_str": False},
    "mol_off":   {"zero_seq": False, "zero_mol": True,  "zero_str": False},
    "str_off":   {"zero_seq": False, "zero_mol": False, "zero_str": True},
    "seq_only":  {"zero_seq": False, "zero_mol": True,  "zero_str": True},
    "mol_only":  {"zero_seq": True,  "zero_mol": False, "zero_str": True},
    "str_only":  {"zero_seq": True,  "zero_mol": True,  "zero_str": False},
    "all_off":   {"zero_seq": True,  "zero_mol": True,  "zero_str": True},
}
CONDITION_ORDER = list(ABLATION_CONDITIONS.keys())
BOOTSTRAP_N = 2000
BOOTSTRAP_SEED = 42

if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))


# ── Utilities ──────────────────────────────────────────────────────────────
def sigmoid_np(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    exp_x = np.exp(x[~pos])
    out[~pos] = exp_x / (1.0 + exp_x)
    return out.astype(np.float32)


def masked_mean(h, mask):
    """Mean-pool hidden states with mask, clamping denominator."""
    denom = mask.sum(dim=(1, 2)).clamp_min(1.0).unsqueeze(-1)
    return (h * mask).sum(dim=1) / denom


def safe_auc(labels, scores):
    if np.unique(labels).size < 2:
        return float("nan")
    return float(roc_auc_score(labels, scores))


def safe_aupr(labels, scores):
    if np.unique(labels).size < 2:
        return float("nan")
    return float(average_precision_score(labels, scores))


def pearson_r(a, b):
    if a.size == 0 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def bootstrap_delta_ci(labels, baseline_scores, ablated_scores, scorer):
    if labels.size == 0 or np.unique(labels).size < 2:
        return float("nan"), float("nan")
    rng = np.random.RandomState(BOOTSTRAP_SEED)
    deltas = []
    for _ in range(BOOTSTRAP_N):
        idx = rng.choice(labels.size, size=labels.size, replace=True)
        y = labels[idx]
        if np.unique(y).size < 2:
            continue
        deltas.append(scorer(y, ablated_scores[idx]) - scorer(y, baseline_scores[idx]))
    if not deltas:
        return float("nan"), float("nan")
    return float(np.percentile(deltas, 2.5)), float(np.percentile(deltas, 97.5))


# ── Config builders ────────────────────────────────────────────────────────
def _find_checkpoint():
    yaml_path = CHECKPOINT_DIR / "complete-full-random-all-0-complex.yml"
    ckpt_path = CHECKPOINT_DIR / "models" / "best-checkpoint.ckpt"
    assert yaml_path.is_file(), f"Missing config: {yaml_path}"
    assert ckpt_path.is_file(), f"Missing checkpoint: {ckpt_path}"
    return yaml_path, ckpt_path


def _base_config(config_yaml, tag):
    from utils import load_config
    config = load_config(str(config_yaml))
    config.num_cpus = 0
    config.num_gpus = 1
    config.data.tag = tag
    config.data.representer = "structure_sequence"
    config.data.full_data = False
    config.data.batch_size = 16
    config.data.sample_weight = [1.0, 1.0]
    config.data.fake_sequence_ratio = 0
    config.data.max_substrate_length = 280
    config.data.max_enzyme_length = 1450
    config.data.features = ["morgan", 1024, "grover_mean", 4885]
    config.data.atom_features = ["grover", 2400]
    config.data.sequence_processed_path = "nonexistent.lmdb"
    return config


def build_p450_config(config_yaml):
    config = _base_config(config_yaml, "e8_p450_ablation")
    config.data.train_data_path = str(P450_CSV)
    config.data.val_data_path = str(P450_CSV)
    config.data.test_data_path = str(P450_CSV)
    config.data.enzyme_lmdb_path = str(P450_SHARED / "enzyme_features.lmdb")
    config.data.reaction_lmdb_path = str(P450_SHARED / "reaction_features.lmdb")
    config.data.grover_path = str(P450_SHARED / "grover_fingerprint.lmdb")
    config.data.morgan_path = str(P450_SHARED / "morgan_fingerprint.npy")
    config.data.structure_processed_path = str(P450_STRUCTURE)
    config.data.high_quality_id_path = str(P450_HQ_IDS)
    return config


def build_phosphatase_config(config_yaml):
    config = _base_config(config_yaml, "e8_phosphatase_ablation")
    config.data.train_data_path = str(PHOSPHATASE_CSV)
    config.data.val_data_path = str(PHOSPHATASE_CSV)
    config.data.test_data_path = str(PHOSPHATASE_CSV)
    config.data.enzyme_lmdb_path = str(PHOSPHATASE_ROOT / "enzyme_features.lmdb")
    config.data.reaction_lmdb_path = str(PHOSPHATASE_ROOT / "reaction_features.lmdb")
    config.data.grover_path = str(PHOSPHATASE_ROOT / "grover_fingerprint.lmdb")
    config.data.morgan_path = str(PHOSPHATASE_ROOT / "morgan_fingerprint.npy")
    config.data.structure_processed_path = str(PHOSPHATASE_ROOT / "structure" / "af2" / "structure_features.lmdb")
    config.data.high_quality_id_path = "nonexistent_hq.txt"
    return config


# ── Ablated forward pass ───────────────────────────────────────────────────
def ablated_forward(model, G, zero_seq=False, zero_mol=False, zero_str=False):
    """Replicate SS.forward() with selective channel zeroing after projection.

    Zeroing strategy (Codex-reviewed):
    - Zero hidden tensors AFTER encoder projection, not raw inputs
    - When one side of cross-attention is completely off, zero attention
      outputs explicitly (don't let bias-influenced attention leak info)
    - When both reaction-atom branches are off, zero x_reaction after aggregator
    """
    # ── Enzyme sequence channel ──
    x_pro = model.protein_mlp(G.embedding)
    x_pro = x_pro.view(-1, 1450, model.hidden_dim)
    if zero_seq:
        x_pro = torch.zeros_like(x_pro)

    # ── Atom-level features (structure + molecular) ──
    atom_features = []
    str_mean = None

    if model.config.model.use_gnn:
        str_mean, x_str = model._get_graph_output(G, x_pro.shape[0])
        if zero_str:
            str_mean = torch.zeros_like(str_mean)
            x_str = torch.zeros_like(x_str)
        atom_features.append(x_str)

    if "grover" in model.config.data.atom_features:
        grover_dim = int(model.config.data.atom_features[
            model.config.data.atom_features.index("grover") + 1
        ])
        x_grover = G.grover.view(-1, 280, grover_dim)
        x_grover = model.atom_feature_mlp_dict["grover"](x_grover)
        if zero_mol:
            x_grover = torch.zeros_like(x_grover)
        atom_features.append(x_grover)

    x_reaction = model.atom_feature_aggregator(atom_features)

    # If both reaction-atom branches are off, zero x_reaction to avoid bias leak
    reaction_atom_off = zero_mol and zero_str
    if reaction_atom_off:
        x_reaction = torch.zeros_like(x_reaction)

    # ── Cross-attention ──
    if model.use_attention:
        # If either side is completely off, zero attention outputs
        interaction_off = zero_seq or reaction_atom_off
        if interaction_off:
            weighted_sum_pro = torch.zeros_like(x_pro)
            weighted_sum_reaction = torch.zeros_like(x_reaction)
        else:
            weighted_sum_pro, _ = model.enzyme_attention(
                x_pro, x_reaction, x_reaction,
                need_weights=False,
                key_padding_mask=G.reaction_padding_mask,
            )
            weighted_sum_reaction, _ = model.reaction_attention(
                x_reaction, x_pro, x_pro,
                need_weights=False,
                key_padding_mask=G.enzyme_padding_mask,
            )

    # ── Masked mean pooling ──
    reaction_mask = (~G.reaction_padding_mask).unsqueeze(-1).float()
    enzyme_mask = (~G.enzyme_padding_mask).unsqueeze(-1).float()

    x_reaction = masked_mean(x_reaction, reaction_mask)
    x_pro = masked_mean(x_pro, enzyme_mask)

    if model.use_attention:
        weighted_sum_reaction = masked_mean(weighted_sum_reaction, reaction_mask)
        weighted_sum_pro = masked_mean(weighted_sum_pro, enzyme_mask)

    # ── Specificity header input ──
    embeddings = [x_pro, x_reaction]
    if model.use_attention:
        embeddings.extend([weighted_sum_pro, weighted_sum_reaction])
    if model.config.model.use_gnn:
        embeddings.append(str_mean)

    if "grover_mean" in model.config.data.features:
        grover_mean_emb = model.feature_mlp_dict["grover_mean"](G.grover_mean)
        if zero_mol:
            grover_mean_emb = torch.zeros_like(grover_mean_emb)
        embeddings.append(grover_mean_emb)

    if "morgan" in model.config.data.features:
        morgan_emb = model.feature_mlp_dict["morgan"](G.morgan)
        if zero_mol:
            morgan_emb = torch.zeros_like(morgan_emb)
        embeddings.append(morgan_emb)

    specificity_output = model.specificity_header(embeddings)
    return specificity_output, [G.tag, G.str_tag]


# ── Inference runner ───────────────────────────────────────────────────────
@torch.inference_mode()
def run_family(family_name, config, ckpt_path, device):
    """Single-pass optimization: load each batch ONCE, run all 8 conditions on it."""
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    dm = Singledataset(config)
    test_loader = dm.test_dataloader()
    template_df = dm.test_prediction_df.copy()

    model = SS.load_from_checkpoint(str(ckpt_path), config=config, map_location="cpu")
    model.eval()
    model.to(device)

    n_total = len(template_df)
    chunks_by_cond = {name: [] for name in CONDITION_ORDER}

    print(f"  [{family_name}] single-pass: 1 loader pass × {len(CONDITION_ORDER)} conditions per batch")
    n_processed = 0
    for batch in test_loader:
        batch = batch.to(device, non_blocking=(device.type == "cuda"))
        for cond_name, flags in ABLATION_CONDITIONS.items():
            logits, _ = ablated_forward(model, batch, **flags)
            chunks_by_cond[cond_name].append(
                logits.squeeze(-1).detach().float().cpu().numpy().ravel()
            )
        n_processed += len(chunks_by_cond["baseline"][-1])
        if n_processed % (16 * 10) < 16:
            print(f"    {n_processed}/{n_total} samples", flush=True)

    outputs = {}
    for cond_name in CONDITION_ORDER:
        chunks = chunks_by_cond[cond_name]
        logits_all = np.concatenate(chunks).astype(np.float32) if chunks else np.zeros(0, np.float32)
        assert len(logits_all) == n_total, f"Mismatch: {len(logits_all)} vs {n_total}"

        df = template_df.copy()
        df["logit"] = logits_all
        df["prob"] = sigmoid_np(logits_all)
        outputs[cond_name] = df

    for cond_name in CONDITION_ORDER:
        auc = safe_auc(outputs[cond_name]["Label"].values, outputs[cond_name]["logit"].values)
        print(f"  [{family_name}] {cond_name}: AUC={auc:.4f}")

    return outputs


# ── Metrics ────────────────────────────────────────────────────────────────
def compute_metrics(df):
    labels = df["Label"].values
    logits = df["logit"].values
    probs = df["prob"].values

    pos_logits = logits[labels == 1]
    neg_logits = logits[labels == 0]
    pos_mean = float(np.mean(pos_logits)) if pos_logits.size else float("nan")
    neg_mean = float(np.mean(neg_logits)) if neg_logits.size else float("nan")

    return {
        "n": int(labels.size),
        "n_pos": int((labels == 1).sum()),
        "n_neg": int((labels == 0).sum()),
        "auc_roc": safe_auc(labels, logits),
        "auc_pr": safe_aupr(labels, logits),
        "brier": float(brier_score_loss(labels, probs)) if labels.size else float("nan"),
        "logit_mean": float(np.mean(logits)),
        "logit_std": float(np.std(logits)),
        "pos_logit_mean": pos_mean,
        "neg_logit_mean": neg_mean,
        "score_gap": pos_mean - neg_mean if pos_logits.size and neg_logits.size else float("nan"),
    }


def summarize_family(family_name, outputs):
    baseline_df = outputs["baseline"]
    bl_labels = baseline_df["Label"].values
    bl_logits = baseline_df["logit"].values
    bl_metrics = compute_metrics(baseline_df)

    rows = []
    for cond_name in CONDITION_ORDER:
        df = outputs[cond_name]
        m = compute_metrics(df)
        logits = df["logit"].values

        ci_low, ci_high = bootstrap_delta_ci(bl_labels, bl_logits, logits, roc_auc_score)
        ci_pr_low, ci_pr_high = bootstrap_delta_ci(
            bl_labels, bl_logits, logits, average_precision_score
        )

        m.update({
            "family": family_name,
            "condition": cond_name,
            "delta_auc": m["auc_roc"] - bl_metrics["auc_roc"] if not np.isnan(m["auc_roc"]) else float("nan"),
            "delta_aupr": m["auc_pr"] - bl_metrics["auc_pr"] if not np.isnan(m["auc_pr"]) else float("nan"),
            "delta_score_gap": m["score_gap"] - bl_metrics["score_gap"] if not np.isnan(m["score_gap"]) else float("nan"),
            "delta_logit_std": m["logit_std"] - bl_metrics["logit_std"],
            "mean_abs_delta_logit": float(np.mean(np.abs(logits - bl_logits))),
            "pearson_r_vs_baseline": pearson_r(logits, bl_logits),
            "delta_auc_ci_low": ci_low,
            "delta_auc_ci_high": ci_high,
            "delta_aupr_ci_low": ci_pr_low,
            "delta_aupr_ci_high": ci_pr_high,
        })
        rows.append(m)

    return pd.DataFrame(rows)


# ── Visualization ──────────────────────────────────────────────────────────
def plot_results(p450_summary, phos_summary, output_dir):
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))

    for ax_idx, (family_name, summary, color) in enumerate([
        ("P450", p450_summary, "crimson"),
        ("Phosphatase", phos_summary, "steelblue"),
    ]):
        # AUC by condition
        ax = axes[0][ax_idx]
        conds = summary["condition"].tolist()
        aucs = summary["auc_roc"].tolist()
        bars = ax.bar(range(len(conds)), aucs, color=color, alpha=0.7)
        ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("AUC-ROC")
        ax.set_title(f"{family_name}: AUC by Condition")
        ax.set_ylim(0.3, 1.0)
        for bar, auc in zip(bars, aucs):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{auc:.3f}", ha="center", fontsize=7, fontweight="bold")

        # Score gap by condition
        ax = axes[1][ax_idx]
        gaps = summary["score_gap"].tolist()
        bars = ax.bar(range(len(conds)), gaps, color=color, alpha=0.7)
        ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
        ax.set_xticks(range(len(conds)))
        ax.set_xticklabels(conds, rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Score Gap (pos - neg mean logit)")
        ax.set_title(f"{family_name}: Score Gap by Condition")
        for bar, gap in zip(bars, gaps):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.01 if gap >= 0 else bar.get_height() - 0.15,
                    f"{gap:.2f}", ha="center", fontsize=7)

    plt.suptitle("E8 Channel Ablation: AUC-ROC and Score Gap", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "e8_channel_ablation.png", dpi=150)
    plt.close()
    print(f"  Plot saved: {output_dir / 'e8_channel_ablation.png'}")


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config_yaml, ckpt_path = _find_checkpoint()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    print(f"Device: {device}")
    print(f"Config: {config_yaml}")
    print(f"Checkpoint: {ckpt_path}")

    # ── P450 ──
    print(f"\n{'='*60}")
    print("  P450 (our data, floor performance)")
    print(f"{'='*60}")
    p450_config = build_p450_config(config_yaml)
    p450_outputs = run_family("P450", p450_config, ckpt_path, device)

    p450_dir = OUTPUT_DIR / "p450"
    p450_dir.mkdir(parents=True, exist_ok=True)
    p450_summary = summarize_family("P450", p450_outputs)
    p450_summary.to_csv(p450_dir / "summary.csv", index=False)
    for cond_name, df in p450_outputs.items():
        df.to_csv(p450_dir / f"predictions_{cond_name}.csv", index=False)

    # ── Phosphatase (control) ──
    print(f"\n{'='*60}")
    print("  Phosphatase (ESIBank control, healthy performance)")
    print(f"{'='*60}")
    phos_config = build_phosphatase_config(config_yaml)
    phos_outputs = run_family("Phosphatase", phos_config, ckpt_path, device)

    phos_dir = OUTPUT_DIR / "phosphatase"
    phos_dir.mkdir(parents=True, exist_ok=True)
    phos_summary = summarize_family("Phosphatase", phos_outputs)
    phos_summary.to_csv(phos_dir / "summary.csv", index=False)
    for cond_name, df in phos_outputs.items():
        df.to_csv(phos_dir / f"predictions_{cond_name}.csv", index=False)

    # ── Combined summary ──
    combined = pd.concat([p450_summary, phos_summary], ignore_index=True)
    combined.to_csv(OUTPUT_DIR / "summary_all.csv", index=False)

    # ── JSON results ──
    results_json = {}
    for fam, summary in [("P450", p450_summary), ("Phosphatase", phos_summary)]:
        results_json[fam] = {}
        for _, row in summary.iterrows():
            cond = row["condition"]
            results_json[fam][cond] = {
                k: (None if (isinstance(v, float) and np.isnan(v)) else v)
                for k, v in row.to_dict().items()
            }
    with open(OUTPUT_DIR / "e8_results.json", "w", encoding="utf-8") as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)

    # ── Visualization ──
    plot_results(p450_summary, phos_summary, OUTPUT_DIR)

    # ── Print summary table ──
    print(f"\n\n{'='*100}")
    print(f"{'Family':<14} {'Condition':<12} {'AUC':>7} {'ΔAUC':>7} {'AUPR':>7} {'Gap':>7} {'LogitStd':>9} {'r_base':>7}")
    print("-" * 100)
    for _, row in combined.iterrows():
        print(f"{row['family']:<14} {row['condition']:<12} "
              f"{row['auc_roc']:>7.4f} {row['delta_auc']:>+7.4f} "
              f"{row['auc_pr']:>7.4f} {row['score_gap']:>7.2f} "
              f"{row['logit_std']:>9.3f} {row['pearson_r_vs_baseline']:>7.3f}")

    # ── Sanity checks ──
    all_off_p450 = p450_summary[p450_summary["condition"] == "all_off"].iloc[0]
    all_off_phos = phos_summary[phos_summary["condition"] == "all_off"].iloc[0]
    print(f"\nSanity check (all_off should have ~zero logit_std):")
    print(f"  P450 all_off: logit_std={all_off_p450['logit_std']:.6f}")
    print(f"  Phosphatase all_off: logit_std={all_off_phos['logit_std']:.6f}")

    baseline_p450 = p450_summary[p450_summary["condition"] == "baseline"].iloc[0]
    baseline_phos = phos_summary[phos_summary["condition"] == "baseline"].iloc[0]
    print(f"\nBaseline verification:")
    print(f"  P450 baseline AUC: {baseline_p450['auc_roc']:.4f} (expected ~0.517)")
    print(f"  Phosphatase baseline AUC: {baseline_phos['auc_roc']:.4f} (expected ~0.877)")

    print(f"\nAll results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
