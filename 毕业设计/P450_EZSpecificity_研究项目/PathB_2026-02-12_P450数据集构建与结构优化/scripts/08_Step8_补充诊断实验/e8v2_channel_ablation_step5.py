"""
E8v2 Channel Ablation on Step 5 data (random negatives + Vina docked structures).

Redo of E8 on the actual problem data (AUC~0.517) instead of
EXP01 (inhibitor negatives, crystal structures, AUC~0.71).

Channels:
    SEQ  - Enzyme sequence (ESM embeddings -> protein_mlp)
    MOL  - Molecular features (grover atom + grover_mean + morgan)
    STR  - Structure features (EGNN: x_str atom + str_mean global)

Conditions (8):
    baseline, seq_off, mol_off, str_off, seq_only, mol_only, str_only, all_off

Method: zero hidden tensors AFTER encoder projection (not raw inputs).
When one side of cross-attention is completely off, zero the attention outputs
explicitly rather than letting bias-influenced attention run.

Performance: encoder caching — compute expensive encoder outputs (protein_mlp,
EGNN, Grover MLP) ONCE per batch, then apply 8 ablation masks cheaply.

Data:
    P450: Step 5 random negatives + Vina docked, 2766 samples, expected AUC~0.517
    Phosphatase: ESIBank control (same as E8v1), 6291 samples, expected AUC~0.877
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

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
OUTPUT_DIR = PATHB / "results" / "08_Step8_补充诊断实验" / "E8v2_通道关闭_Step5"

# P450 — Step 5 data (random negatives + Vina docked structures)
P450_CSV = PATHB / "results" / "04_Step4_批量对接" / "data.csv"
P450_SHARED = PATHB / "data" / "00_shared" / "features"
P450_STRUCTURE = PATHB / "data" / "05_Step5_重构评估" / "structure_features.lmdb"
P450_HQ_IDS = PATHB / "data" / "05_Step5_重构评估" / "high_quality_id.txt"

# Phosphatase (control — same as E8v1)
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
    config = _base_config(config_yaml, "e8v2_p450_step5")
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
    config = _base_config(config_yaml, "e8v2_phosphatase")
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


# ── Inference runner with encoder caching ─────────────────────────────────
@torch.inference_mode()
def run_family(family_name, config, ckpt_path, device):
    """Optimized: encode once per batch, decode 8 conditions with cached outputs.

    Expensive ops (protein_mlp, EGNN, Grover MLP) run 1× per batch.
    Cheap ops (mask, aggregator, cross-attention, header) run 8× per batch.
    """
    from Datasets.brenda import Singledataset
    from Models.ss import SS

    dm = Singledataset(config)
    test_loader = dm.test_dataloader()
    template_df = dm.test_prediction_df.copy()

    # ── Preflight: detect silent row filtering (Codex HIGH fix) ──
    raw_df = pd.read_csv(config.data.test_data_path)
    raw_n = len(raw_df)
    kept_n = len(template_df)
    if kept_n != raw_n:
        dropped = raw_n - kept_n
        print(f"  [{family_name}] WARNING: dataloader filtered {dropped} rows "
              f"({raw_n} -> {kept_n})", flush=True)
        if family_name == "P450":
            raise RuntimeError(
                f"P450: {dropped}/{raw_n} rows silently filtered. "
                "Check shared LMDB coverage and Step 5 structure/HQ alignment."
            )

    model = SS.load_from_checkpoint(str(ckpt_path), config=config, map_location="cpu")
    model.eval()
    model.to(device)

    n_total = len(template_df)
    batch_size = int(config.data.batch_size)
    n_batches = (n_total + batch_size - 1) // batch_size
    chunks_by_cond = {name: [] for name in CONDITION_ORDER}

    # Pre-check model features
    use_gnn = model.config.model.use_gnn
    use_attention = model.use_attention
    has_grover_atom = "grover" in model.config.data.atom_features
    has_grover_mean = "grover_mean" in model.config.data.features
    has_morgan = "morgan" in model.config.data.features
    grover_dim = 0
    if has_grover_atom:
        grover_dim = int(model.config.data.atom_features[
            model.config.data.atom_features.index("grover") + 1
        ])

    print(f"  [{family_name}] {n_total} samples, ~{n_batches} batches × "
          f"{len(CONDITION_ORDER)} conditions (encoder cached)")

    t0 = time.time()
    t_encode = 0.0
    t_decode = 0.0
    n_processed = 0

    for batch in test_loader:
        batch = batch.to(device, non_blocking=(device.type == "cuda"))

        # ══ ENCODE ONCE (expensive: protein_mlp, EGNN, Grover MLP) ══
        te0 = time.time()

        x_pro_full = model.protein_mlp(batch.embedding)
        x_pro_full = x_pro_full.view(-1, 1450, model.hidden_dim)

        str_mean_full = x_str_full = None
        if use_gnn:
            str_mean_full, x_str_full = model._get_graph_output(batch, x_pro_full.shape[0])

        x_grover_full = None
        if has_grover_atom:
            x_grover_raw = batch.grover.view(-1, 280, grover_dim)
            x_grover_full = model.atom_feature_mlp_dict["grover"](x_grover_raw)

        grover_mean_emb_full = None
        if has_grover_mean:
            grover_mean_emb_full = model.feature_mlp_dict["grover_mean"](batch.grover_mean)

        morgan_emb_full = None
        if has_morgan:
            morgan_emb_full = model.feature_mlp_dict["morgan"](batch.morgan)

        reaction_mask = (~batch.reaction_padding_mask).unsqueeze(-1).float()
        enzyme_mask = (~batch.enzyme_padding_mask).unsqueeze(-1).float()

        t_encode += time.time() - te0

        # ══ DECODE PER CONDITION (cheap: mask + aggregator + attention + header) ══
        td0 = time.time()

        for cond_name, flags in ABLATION_CONDITIONS.items():
            zero_seq = flags["zero_seq"]
            zero_mol = flags["zero_mol"]
            zero_str = flags["zero_str"]

            # Apply channel masks to cached encoder outputs
            x_pro = torch.zeros_like(x_pro_full) if zero_seq else x_pro_full

            atom_features = []
            str_mean = None
            if use_gnn:
                x_str = torch.zeros_like(x_str_full) if zero_str else x_str_full
                str_mean = torch.zeros_like(str_mean_full) if zero_str else str_mean_full
                atom_features.append(x_str)
            if has_grover_atom:
                x_grover = torch.zeros_like(x_grover_full) if zero_mol else x_grover_full
                atom_features.append(x_grover)

            x_reaction = model.atom_feature_aggregator(atom_features)

            reaction_atom_off = zero_mol and zero_str
            if reaction_atom_off:
                x_reaction = torch.zeros_like(x_reaction)

            # Cross-attention
            if use_attention:
                interaction_off = zero_seq or reaction_atom_off
                if interaction_off:
                    weighted_sum_pro = torch.zeros_like(x_pro)
                    weighted_sum_reaction = torch.zeros_like(x_reaction)
                else:
                    weighted_sum_pro, _ = model.enzyme_attention(
                        x_pro, x_reaction, x_reaction,
                        need_weights=False,
                        key_padding_mask=batch.reaction_padding_mask,
                    )
                    weighted_sum_reaction, _ = model.reaction_attention(
                        x_reaction, x_pro, x_pro,
                        need_weights=False,
                        key_padding_mask=batch.enzyme_padding_mask,
                    )

            # Masked mean pooling
            x_reaction_pooled = masked_mean(x_reaction, reaction_mask)
            x_pro_pooled = masked_mean(x_pro, enzyme_mask)

            embeddings = [x_pro_pooled, x_reaction_pooled]
            if use_attention:
                ws_reaction_pooled = masked_mean(weighted_sum_reaction, reaction_mask)
                ws_pro_pooled = masked_mean(weighted_sum_pro, enzyme_mask)
                embeddings.extend([ws_pro_pooled, ws_reaction_pooled])
            if use_gnn:
                embeddings.append(str_mean)
            if has_grover_mean:
                emb = torch.zeros_like(grover_mean_emb_full) if zero_mol else grover_mean_emb_full
                embeddings.append(emb)
            if has_morgan:
                emb = torch.zeros_like(morgan_emb_full) if zero_mol else morgan_emb_full
                embeddings.append(emb)

            logits = model.specificity_header(embeddings)
            chunks_by_cond[cond_name].append(
                logits.squeeze(-1).detach().float().cpu().numpy().ravel()
            )

        t_decode += time.time() - td0

        n_processed += len(chunks_by_cond["baseline"][-1])
        if n_processed % (batch_size * 20) < batch_size:
            elapsed = time.time() - t0
            rate = n_processed / elapsed if elapsed > 0 else 0
            eta = (n_total - n_processed) / rate if rate > 0 else 0
            print(f"    {n_processed}/{n_total} samples "
                  f"({elapsed:.0f}s, ~{eta:.0f}s left, "
                  f"encode={t_encode:.0f}s decode={t_decode:.0f}s)", flush=True)

    elapsed_total = time.time() - t0
    print(f"  [{family_name}] done in {elapsed_total:.1f}s "
          f"({n_total / elapsed_total:.1f} samples/s)")
    print(f"  [{family_name}] timing: encode={t_encode:.1f}s ({t_encode/elapsed_total*100:.0f}%) "
          f"decode={t_decode:.1f}s ({t_decode/elapsed_total*100:.0f}%)")

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
        ("P450 (Step5 Vina)", p450_summary, "crimson"),
        ("Phosphatase", phos_summary, "steelblue"),
    ]):
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

    plt.suptitle("E8v2 Channel Ablation (Step5 Vina Data): AUC-ROC and Score Gap",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "e8v2_channel_ablation.png", dpi=150)
    plt.close()
    print(f"  Plot saved: {output_dir / 'e8v2_channel_ablation.png'}")


# ── Main ───────────────────────────────────────────────────────────────────
def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    config_yaml, ckpt_path = _find_checkpoint()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    print(f"Device: {device}")
    print(f"Config: {config_yaml}")
    print(f"Checkpoint: {ckpt_path}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"\nE8v2: Channel ablation on Step 5 data (random negatives + Vina)")
    print(f"  P450 CSV: {P450_CSV}")
    print(f"  P450 structure: {P450_STRUCTURE}")
    print(f"  P450 HQ IDs: {P450_HQ_IDS}")

    # ── P450 (Step 5 data) ──
    print(f"\n{'='*60}")
    print("  P450 (Step 5: random negatives + Vina docking)")
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
    with open(OUTPUT_DIR / "e8v2_results.json", "w", encoding="utf-8") as f:
        json.dump(results_json, f, indent=2, ensure_ascii=False)

    # ── Visualization ──
    plot_results(p450_summary, phos_summary, OUTPUT_DIR)

    # ── Print summary table ──
    print(f"\n\n{'='*100}")
    print(f"{'Family':<14} {'Condition':<12} {'AUC':>7} {'dAUC':>7} {'AUPR':>7} {'Gap':>7} {'LogitStd':>9} {'r_base':>7}")
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

    print(f"\n{'='*60}")
    print("  E8v1 (EXP01) vs E8v2 (Step5) comparison notes")
    print(f"{'='*60}")
    print("  Key question: does STR dominance hold on Vina-docked data?")
    print("  If str_off IMPROVES AUC -> Vina structures are actively hurting.")

    print(f"\nAll results saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
