"""
E7 Step 2: ESIBank P450 Internal Benchmark — Full Analysis (E1'-E6')

Runs all 6 diagnostic experiments on ESIBank P450 predictions and compares
with our independent P450 (Step 5, AUC=0.517).

Requires: esibank_p450_predictions.csv from e7_esibank_p450_prep_and_inference.py
"""
from __future__ import annotations

import sys
import os
import warnings
import json
import pickle
import lmdb

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score, r2_score
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

# --- Paths ---
PROJECT_ROOT = Path(r"D:\EZSpecificity_Project")
PATHB = PROJECT_ROOT / "毕业设计" / "P450_EZSpecificity_研究项目" / "PathB_2026-02-12_P450数据集构建与结构优化"
E7_DIR = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E7_ESIBank_P450_内部基准"
E6_DIR = PATHB / "results" / "07_Step7_Tier1_诊断实验" / "E6_expansion_多家族推理"

G_BRENDA = Path(r"G:\.shortcut-targets-by-id\173a36NiOLgXcvzvJjRDH29y2xd7Ey3Pr\ESIBank\brenda")

# ESIBank P450 predictions (from Step 1)
ESIBANK_PRED = E7_DIR / "esibank_p450_predictions.csv"
# Our P450 predictions (from Step 5)
OUR_PRED = PATHB / "results" / "05_Step5_重构评估" / "predictions.csv"
# Our P450 ESM embeddings
OUR_ESM_LMDB = PATHB / "data" / "00_shared" / "features" / "enzyme_features.lmdb"
# Our Morgan fingerprints
OUR_MORGAN = PATHB / "data" / "00_shared" / "features" / "morgan_fingerprint.npy"
# ESIBank Morgan
ESIBANK_MORGAN = G_BRENDA / "morgan_fingerprint.npy"
# ESIBank enzyme features
ESIBANK_ESM_LMDB = G_BRENDA / "enzyme_features.lmdb"
# E6 existing results
E6_RESULTS_JSON = E6_DIR / "e6_expansion_results.json"

ALL_RESULTS = {}


# ============================================================
# Utility functions
# ============================================================

def _safe_json(obj):
    """Convert numpy types and NaN/Inf for JSON serialization."""
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        if np.isnan(v) or np.isinf(v):
            return None
        return v
    if isinstance(obj, float):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return obj
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, dict):
        return {k: _safe_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_safe_json(v) for v in obj]
    return obj


def _load_esm_embeddings_from_lmdb(lmdb_path, enzyme_indices):
    """Load and mean-pool ESM-2 embeddings from LMDB."""
    env = lmdb.open(str(lmdb_path), readonly=True, lock=False, subdir=False)
    embeddings = {}
    with env.begin() as txn:
        for idx in enzyme_indices:
            key = str(idx).encode()
            val = txn.get(key)
            if val is not None:
                data = pickle.loads(val)
                emb = data["embedding"] if isinstance(data, dict) else data
                if isinstance(emb, np.ndarray):
                    if emb.ndim == 2:
                        emb = emb.mean(axis=0)
                    embeddings[idx] = emb.astype(np.float32)
    env.close()
    return embeddings


# ============================================================
# E1': Score Distribution Probe
# ============================================================

def run_e1_prime(esb_df, our_df):
    """E1': Compare score distributions between ESIBank P450 and our P450."""
    print("\n" + "=" * 70)
    print("  E1': Score Distribution Probe")
    print("=" * 70)

    results = {}

    # --- Analysis A: Logit distributions ---
    esb_pos = esb_df[esb_df["Label"] == 1]["logit"].values
    esb_neg = esb_df[esb_df["Label"] == 0]["logit"].values
    our_pos = our_df[our_df["Label"] == 1]["logit"].values
    our_neg = our_df[our_df["Label"] == 0]["logit"].values

    esb_gap = esb_pos.mean() - esb_neg.mean()
    our_gap = our_pos.mean() - our_neg.mean()

    results["analysis_A"] = {
        "esibank": {
            "pos_mean": float(esb_pos.mean()), "pos_std": float(esb_pos.std()),
            "neg_mean": float(esb_neg.mean()), "neg_std": float(esb_neg.std()),
            "gap": float(esb_gap), "n_pos": len(esb_pos), "n_neg": len(esb_neg),
        },
        "ours": {
            "pos_mean": float(our_pos.mean()), "pos_std": float(our_pos.std()),
            "neg_mean": float(our_neg.mean()), "neg_std": float(our_neg.std()),
            "gap": float(our_gap), "n_pos": len(our_pos), "n_neg": len(our_neg),
        },
    }

    # KS tests: ESIBank pos vs our pos, ESIBank neg vs our neg
    ks_pos = stats.ks_2samp(esb_pos, our_pos)
    ks_neg = stats.ks_2samp(esb_neg, our_neg)
    results["analysis_A"]["ks_pos_vs_pos"] = {"statistic": float(ks_pos.statistic), "pvalue": float(ks_pos.pvalue)}
    results["analysis_A"]["ks_neg_vs_neg"] = {"statistic": float(ks_neg.statistic), "pvalue": float(ks_neg.pvalue)}

    print(f"  ESIBank: pos_mean={esb_pos.mean():.3f}, neg_mean={esb_neg.mean():.3f}, gap={esb_gap:.3f}")
    print(f"  Ours:    pos_mean={our_pos.mean():.3f}, neg_mean={our_neg.mean():.3f}, gap={our_gap:.3f}")
    print(f"  KS(pos vs pos): stat={ks_pos.statistic:.3f}, p={ks_pos.pvalue:.2e}")
    print(f"  KS(neg vs neg): stat={ks_neg.statistic:.3f}, p={ks_neg.pvalue:.2e}")

    # --- Analysis B: Tanimoto chemical space ---
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs

        esb_morgan_all = np.load(str(ESIBANK_MORGAN))
        our_morgan_all = np.load(str(OUR_MORGAN))

        esb_sub_indices = esb_df["Substrate Index"].unique()
        our_sub_indices = our_df["Substrate Index"].unique()

        esb_fps = esb_morgan_all[esb_sub_indices]
        our_fps = our_morgan_all[our_sub_indices]

        def _tanimoto_mean(fps_a, fps_b=None):
            """Compute mean Tanimoto similarity. If fps_b is None, compute within fps_a."""
            sims = []
            if fps_b is None:
                for i in range(len(fps_a)):
                    for j in range(i + 1, len(fps_a)):
                        inter = np.sum(fps_a[i] * fps_a[j])
                        union = np.sum(fps_a[i]) + np.sum(fps_a[j]) - inter
                        sims.append(inter / union if union > 0 else 0)
            else:
                # Sample if too many pairs
                n_pairs = len(fps_a) * len(fps_b)
                if n_pairs > 500000:
                    idx_a = np.random.choice(len(fps_a), min(500, len(fps_a)), replace=False)
                    idx_b = np.random.choice(len(fps_b), min(500, len(fps_b)), replace=False)
                    fps_a = fps_a[idx_a]
                    fps_b = fps_b[idx_b]
                for i in range(len(fps_a)):
                    for j in range(len(fps_b)):
                        inter = np.sum(fps_a[i] * fps_b[j])
                        union = np.sum(fps_a[i]) + np.sum(fps_b[j]) - inter
                        sims.append(inter / union if union > 0 else 0)
            return float(np.mean(sims)) if sims else 0.0

        esb_intra = _tanimoto_mean(esb_fps) if len(esb_fps) > 1 else 0.0
        our_intra = _tanimoto_mean(our_fps) if len(our_fps) > 1 else 0.0
        cross = _tanimoto_mean(esb_fps, our_fps) if len(esb_fps) > 0 and len(our_fps) > 0 else 0.0

        results["analysis_B"] = {
            "esibank_intra_tanimoto": esb_intra,
            "ours_intra_tanimoto": our_intra,
            "cross_tanimoto": cross,
            "n_esibank_substrates": len(esb_sub_indices),
            "n_our_substrates": len(our_sub_indices),
        }
        print(f"  Tanimoto intra (ESIBank): {esb_intra:.4f} (n={len(esb_sub_indices)})")
        print(f"  Tanimoto intra (ours):    {our_intra:.4f} (n={len(our_sub_indices)})")
        print(f"  Tanimoto cross:           {cross:.4f}")
    except Exception as e:
        print(f"  WARNING: Tanimoto analysis failed: {e}")
        results["analysis_B"] = {"error": str(e)}

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    ax = axes[0]
    ax.hist(esb_pos, bins=50, alpha=0.5, label=f"ESIBank pos (n={len(esb_pos)})", color="blue", density=True)
    ax.hist(esb_neg, bins=50, alpha=0.5, label=f"ESIBank neg (n={len(esb_neg)})", color="red", density=True)
    ax.set_xlabel("Logit")
    ax.set_ylabel("Density")
    ax.set_title(f"ESIBank P450 (gap={esb_gap:.2f})")
    ax.legend(fontsize=8)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

    ax = axes[1]
    ax.hist(our_pos, bins=50, alpha=0.5, label=f"Our pos (n={len(our_pos)})", color="blue", density=True)
    ax.hist(our_neg, bins=50, alpha=0.5, label=f"Our neg (n={len(our_neg)})", color="red", density=True)
    ax.set_xlabel("Logit")
    ax.set_ylabel("Density")
    ax.set_title(f"Our P450 (gap={our_gap:.2f})")
    ax.legend(fontsize=8)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

    plt.suptitle("E1': Score Distribution Comparison", fontweight="bold")
    plt.tight_layout()
    plt.savefig(E7_DIR / "e1_prime_score_distribution.png", dpi=150)
    plt.close()

    ALL_RESULTS["E1_prime"] = results
    print("  Done.")


# ============================================================
# E2': ESM Embedding Similarity
# ============================================================

def run_e2_prime(esb_df, our_df):
    """E2': ESM-2 embedding similarity analysis."""
    print("\n" + "=" * 70)
    print("  E2': ESM Embedding Similarity")
    print("=" * 70)

    results = {}

    esb_enzyme_indices = sorted(esb_df["Enzyme Index"].unique())
    our_enzyme_indices = sorted(our_df["Enzyme Index"].unique())

    print(f"  Loading ESIBank ESM embeddings ({len(esb_enzyme_indices)} enzymes)...")
    esb_embs = _load_esm_embeddings_from_lmdb(ESIBANK_ESM_LMDB, esb_enzyme_indices)
    print(f"  Loaded: {len(esb_embs)} / {len(esb_enzyme_indices)}")

    print(f"  Loading our ESM embeddings ({len(our_enzyme_indices)} enzymes)...")
    our_embs = _load_esm_embeddings_from_lmdb(OUR_ESM_LMDB, our_enzyme_indices)
    print(f"  Loaded: {len(our_embs)} / {len(our_enzyme_indices)}")

    # --- Analysis A: ESIBank P450 intra-family cosine similarity ---
    esb_vecs = np.array([esb_embs[i] for i in sorted(esb_embs.keys())])
    if len(esb_vecs) > 1:
        from sklearn.metrics.pairwise import cosine_similarity
        sim_matrix = cosine_similarity(esb_vecs)
        upper = sim_matrix[np.triu_indices(len(esb_vecs), k=1)]
        results["analysis_A"] = {
            "n_enzymes": len(esb_vecs),
            "cosine_mean": float(upper.mean()),
            "cosine_std": float(upper.std()),
            "above_090": float((upper > 0.90).mean()),
            "above_095": float((upper > 0.95).mean()),
        }
        print(f"  ESIBank intra cosine: mean={upper.mean():.4f}, std={upper.std():.4f}")
        print(f"  >0.90: {(upper > 0.90).mean()*100:.1f}%, >0.95: {(upper > 0.95).mean()*100:.1f}%")

    # --- Analysis B: PCA ---
    if len(esb_vecs) > 2:
        pca = PCA()
        pca.fit(esb_vecs)
        cumvar = np.cumsum(pca.explained_variance_ratio_)
        n90 = int(np.searchsorted(cumvar, 0.90) + 1)
        n95 = int(np.searchsorted(cumvar, 0.95) + 1)
        results["analysis_B_pca"] = {"n_dims_90pct": n90, "n_dims_95pct": n95}
        print(f"  PCA: 90% variance in {n90} dims, 95% in {n95} dims")

    # --- Analysis C: Cross-group similarity ---
    our_vecs = np.array([our_embs[i] for i in sorted(our_embs.keys())])
    if len(esb_vecs) > 0 and len(our_vecs) > 0:
        from sklearn.metrics.pairwise import cosine_similarity as cos_sim
        cross_sim = cos_sim(esb_vecs, our_vecs)
        results["analysis_C_cross"] = {
            "n_esibank": len(esb_vecs),
            "n_ours": len(our_vecs),
            "cross_cosine_mean": float(cross_sim.mean()),
            "cross_cosine_std": float(cross_sim.std()),
        }
        print(f"  Cross-group cosine: mean={cross_sim.mean():.4f}, std={cross_sim.std():.4f}")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    if len(esb_vecs) > 1:
        ax = axes[0]
        ax.hist(upper, bins=80, alpha=0.7, color="steelblue", density=True, label="ESIBank P450 intra")
        if len(our_vecs) > 1:
            our_sim = cosine_similarity(our_vecs)
            our_upper = our_sim[np.triu_indices(len(our_vecs), k=1)]
            ax.hist(our_upper, bins=80, alpha=0.5, color="crimson", density=True, label="Our P450 intra")
        ax.set_xlabel("Cosine Similarity")
        ax.set_ylabel("Density")
        ax.set_title("Intra-Family ESM-2 Cosine Similarity")
        ax.legend(fontsize=8)

    if len(esb_vecs) > 0 and len(our_vecs) > 0:
        ax = axes[1]
        ax.hist(cross_sim.ravel(), bins=80, alpha=0.7, color="green", density=True)
        ax.set_xlabel("Cosine Similarity")
        ax.set_ylabel("Density")
        ax.set_title("Cross-Group (ESIBank vs Our P450)")

    plt.suptitle("E2': ESM-2 Embedding Similarity", fontweight="bold")
    plt.tight_layout()
    plt.savefig(E7_DIR / "e2_prime_esm_similarity.png", dpi=150)
    plt.close()

    ALL_RESULTS["E2_prime"] = results
    print("  Done.")


# ============================================================
# E3': Per-Substrate / Per-Enzyme AUC + Direction Analysis
# ============================================================

def run_e3_prime(esb_df):
    """E3': Per-substrate and per-enzyme AUC analysis."""
    print("\n" + "=" * 70)
    print("  E3': Per-Substrate / Per-Enzyme AUC")
    print("=" * 70)

    results = {}

    # --- Analysis A: Per-substrate AUC ---
    per_sub = []
    for sub_idx, grp in esb_df.groupby("Substrate Index"):
        n_pos = (grp["Label"] == 1).sum()
        n_neg = (grp["Label"] == 0).sum()
        if n_pos >= 1 and n_neg >= 5:
            auc = roc_auc_score(grp["Label"], grp["logit"])
            per_sub.append({"Substrate Index": sub_idx, "auc": auc, "n_pos": n_pos, "n_neg": n_neg})
    per_sub_df = pd.DataFrame(per_sub)

    if len(per_sub_df) > 0:
        results["per_substrate"] = {
            "n_qualifying": len(per_sub_df),
            "median_auc": float(per_sub_df["auc"].median()),
            "mean_auc": float(per_sub_df["auc"].mean()),
            "std_auc": float(per_sub_df["auc"].std()),
            "above_075": float((per_sub_df["auc"] > 0.75).mean()),
            "above_065": float((per_sub_df["auc"] > 0.65).mean()),
            "below_035": float((per_sub_df["auc"] < 0.35).mean()),
        }
        print(f"  Per-substrate ({len(per_sub_df)} qualifying):")
        print(f"    median={per_sub_df['auc'].median():.3f}, mean={per_sub_df['auc'].mean():.3f}, std={per_sub_df['auc'].std():.3f}")
        print(f"    >0.75: {(per_sub_df['auc'] > 0.75).mean()*100:.1f}%, <0.35: {(per_sub_df['auc'] < 0.35).mean()*100:.1f}%")

    per_sub_df.to_csv(E7_DIR / "e3_prime_per_substrate_auc.csv", index=False)

    # --- Analysis B: Per-enzyme AUC ---
    per_enz = []
    for enz_idx, grp in esb_df.groupby("Enzyme Index"):
        n_pos = (grp["Label"] == 1).sum()
        n_neg = (grp["Label"] == 0).sum()
        if n_pos >= 1 and n_neg >= 5:
            auc = roc_auc_score(grp["Label"], grp["logit"])
            per_enz.append({"Enzyme Index": enz_idx, "auc": auc, "n_pos": n_pos, "n_neg": n_neg})
    per_enz_df = pd.DataFrame(per_enz)

    if len(per_enz_df) > 0:
        results["per_enzyme"] = {
            "n_qualifying": len(per_enz_df),
            "median_auc": float(per_enz_df["auc"].median()),
            "mean_auc": float(per_enz_df["auc"].mean()),
        }
        print(f"  Per-enzyme ({len(per_enz_df)} qualifying): median={per_enz_df['auc'].median():.3f}, mean={per_enz_df['auc'].mean():.3f}")

    per_enz_df.to_csv(E7_DIR / "e3_prime_per_enzyme_auc.csv", index=False)

    # --- Analysis C: Direction analysis ---
    dir_a = 0
    dir_b = 0
    for sub_idx, grp in esb_df.groupby("Substrate Index"):
        n_pos = (grp["Label"] == 1).sum()
        n_neg = (grp["Label"] == 0).sum()
        if n_pos >= 1 and n_neg >= 1:
            pos_mean = grp[grp["Label"] == 1]["logit"].mean()
            neg_mean = grp[grp["Label"] == 0]["logit"].mean()
            if pos_mean > neg_mean:
                dir_b += 1
            else:
                dir_a += 1
    total_dir = dir_a + dir_b
    results["direction"] = {
        "dir_a": dir_a, "dir_b": dir_b,
        "dir_a_pct": float(dir_a / total_dir * 100) if total_dir > 0 else 0,
        "dir_b_pct": float(dir_b / total_dir * 100) if total_dir > 0 else 0,
    }
    print(f"  Direction: A={dir_a} ({dir_a/total_dir*100:.1f}%), B={dir_b} ({dir_b/total_dir*100:.1f}%)")

    # --- Plot ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    if len(per_sub_df) > 0:
        axes[0].hist(per_sub_df["auc"], bins=30, alpha=0.7, color="steelblue", edgecolor="white")
        axes[0].axvline(0.5, color="gray", linestyle="--", alpha=0.5)
        axes[0].axvline(per_sub_df["auc"].median(), color="red", linestyle="-", alpha=0.8, label=f"median={per_sub_df['auc'].median():.3f}")
        axes[0].set_xlabel("AUC-ROC")
        axes[0].set_ylabel("Count")
        axes[0].set_title(f"Per-Substrate AUC (n={len(per_sub_df)})")
        axes[0].legend()

    if len(per_enz_df) > 0:
        axes[1].hist(per_enz_df["auc"], bins=30, alpha=0.7, color="darkorange", edgecolor="white")
        axes[1].axvline(0.5, color="gray", linestyle="--", alpha=0.5)
        axes[1].axvline(per_enz_df["auc"].median(), color="red", linestyle="-", alpha=0.8, label=f"median={per_enz_df['auc'].median():.3f}")
        axes[1].set_xlabel("AUC-ROC")
        axes[1].set_ylabel("Count")
        axes[1].set_title(f"Per-Enzyme AUC (n={len(per_enz_df)})")
        axes[1].legend()

    plt.suptitle("E3': Per-Substrate / Per-Enzyme AUC", fontweight="bold")
    plt.tight_layout()
    plt.savefig(E7_DIR / "e3_prime_per_sub_enz_auc.png", dpi=150)
    plt.close()

    ALL_RESULTS["E3_prime"] = results
    print("  Done.")


# ============================================================
# E4': Score Decomposition (Core)
# ============================================================

def run_e4_prime(esb_df):
    """E4': LOO score decomposition — R²(enzyme), R²(substrate), residual AUC."""
    print("\n" + "=" * 70)
    print("  E4': Score Decomposition (Core)")
    print("=" * 70)

    logits = esb_df["logit"].values.astype(np.float64)
    labels = esb_df["Label"].values
    enzymes = esb_df["Enzyme Index"].values
    substrates = esb_df["Substrate Index"].values
    grand_mean = logits.mean()
    n = len(logits)

    # Precompute group sums and counts for LOO
    enz_sum = {}
    enz_cnt = {}
    for i in range(n):
        e = enzymes[i]
        enz_sum[e] = enz_sum.get(e, 0.0) + logits[i]
        enz_cnt[e] = enz_cnt.get(e, 0) + 1

    sub_sum = {}
    sub_cnt = {}
    for i in range(n):
        s = substrates[i]
        sub_sum[s] = sub_sum.get(s, 0.0) + logits[i]
        sub_cnt[s] = sub_cnt.get(s, 0) + 1

    # LOO predictions
    enz_pred_loo = np.empty(n, dtype=np.float64)
    sub_pred_loo = np.empty(n, dtype=np.float64)
    for i in range(n):
        e = enzymes[i]
        s = substrates[i]
        if enz_cnt[e] > 1:
            enz_pred_loo[i] = (enz_sum[e] - logits[i]) / (enz_cnt[e] - 1)
        else:
            enz_pred_loo[i] = grand_mean
        if sub_cnt[s] > 1:
            sub_pred_loo[i] = (sub_sum[s] - logits[i]) / (sub_cnt[s] - 1)
        else:
            sub_pred_loo[i] = grand_mean

    additive_pred_loo = enz_pred_loo + sub_pred_loo - grand_mean
    residual = logits - additive_pred_loo

    # R² scores
    ss_tot = np.sum((logits - grand_mean) ** 2)
    r2_enz = float(1 - np.sum((logits - enz_pred_loo) ** 2) / ss_tot) if ss_tot > 0 else 0.0
    r2_sub = float(1 - np.sum((logits - sub_pred_loo) ** 2) / ss_tot) if ss_tot > 0 else 0.0
    r2_add = float(1 - np.sum((logits - additive_pred_loo) ** 2) / ss_tot) if ss_tot > 0 else 0.0

    # AUCs
    auc_enz = float(roc_auc_score(labels, enz_pred_loo)) if len(set(labels)) == 2 else float("nan")
    auc_sub = float(roc_auc_score(labels, sub_pred_loo)) if len(set(labels)) == 2 else float("nan")
    auc_add = float(roc_auc_score(labels, additive_pred_loo)) if len(set(labels)) == 2 else float("nan")
    auc_res = float(roc_auc_score(labels, residual)) if len(set(labels)) == 2 else float("nan")

    results = {
        "r2_enzyme": r2_enz,
        "r2_substrate": r2_sub,
        "r2_additive": r2_add,
        "auc_enzyme_only": auc_enz,
        "auc_substrate_only": auc_sub,
        "auc_additive": auc_add,
        "auc_residual": auc_res,
        "n_unique_enzymes": len(enz_cnt),
        "n_unique_substrates": len(sub_cnt),
    }

    print(f"  R2(enzyme):    {r2_enz:.4f}  AUC: {auc_enz:.4f}")
    print(f"  R2(substrate): {r2_sub:.4f}  AUC: {auc_sub:.4f}")
    print(f"  R2(additive):  {r2_add:.4f}  AUC: {auc_add:.4f}")
    print(f"  Residual AUC:  {auc_res:.4f}")
    print(f"  Unique enzymes: {len(enz_cnt)}, substrates: {len(sub_cnt)}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(8, 5))
    components = ["Enzyme\nonly", "Substrate\nonly", "Additive", "Residual"]
    r2_vals = [r2_enz, r2_sub, r2_add, None]
    auc_vals = [auc_enz, auc_sub, auc_add, auc_res]

    x = np.arange(len(components))
    bars = ax.bar(x, auc_vals, color=["steelblue", "darkorange", "green", "gray"], alpha=0.7)
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5, label="Random (0.5)")
    ax.set_xticks(x)
    ax.set_xticklabels(components)
    ax.set_ylabel("AUC-ROC")
    ax.set_title("E4': Score Decomposition — ESIBank P450")
    for i, (bar, auc) in enumerate(zip(bars, auc_vals)):
        r2_text = f"R2={r2_vals[i]:.3f}\n" if r2_vals[i] is not None else ""
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{r2_text}AUC={auc:.3f}", ha="center", fontsize=9, fontweight="bold")
    ax.legend()
    ax.set_ylim(0.3, max(auc_vals) + 0.15)

    plt.tight_layout()
    plt.savefig(E7_DIR / "e4_prime_score_decomposition.png", dpi=150)
    plt.close()

    ALL_RESULTS["E4_prime"] = results
    print("  Done.")


# ============================================================
# E5': Fingerprint Classifier Baseline
# ============================================================

def run_e5_prime(esb_df):
    """E5': Morgan fingerprint classifier with GroupKFold."""
    print("\n" + "=" * 70)
    print("  E5': Fingerprint Classifier Baseline")
    print("=" * 70)

    results = {}

    # Load Morgan fingerprints
    morgan_all = np.load(str(ESIBANK_MORGAN))
    sub_indices = esb_df["Substrate Index"].values
    X = morgan_all[sub_indices]
    y = esb_df["Label"].values
    groups = sub_indices

    print(f"  Samples: {len(X)}, Features: {X.shape[1]}")
    print(f"  Pos: {y.sum()}, Neg: {(1-y).sum()}")
    print(f"  Unique substrates (groups): {len(np.unique(groups))}")

    # Dir A ratio
    sub_pos_rates = esb_df.groupby("Substrate Index")["Label"].mean()
    global_pos_rate = y.mean()
    dir_a_approx = ((sub_pos_rates - global_pos_rate).abs() < 0.05).mean()
    results["dir_a_approx_pct"] = float(dir_a_approx * 100)
    print(f"  Dir A approx (pos rate within 5% of global): {dir_a_approx*100:.1f}%")

    classifiers = {
        "LogisticRegression": LogisticRegression(C=1.0, class_weight="balanced", max_iter=500, random_state=42),
        "RandomForest": RandomForestClassifier(n_estimators=100, max_depth=5, class_weight="balanced", random_state=42),
    }

    cv_strategies = {
        "StratifiedKFold": StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
        "StratifiedGroupKFold": StratifiedGroupKFold(n_splits=5),
    }

    for clf_name, clf_template in classifiers.items():
        for cv_name, cv in cv_strategies.items():
            aucs = []
            pr_aucs = []
            splits = cv.split(X, y, groups) if "Group" in cv_name else cv.split(X, y)
            for train_idx, test_idx in splits:
                clf = type(clf_template)(**clf_template.get_params())
                clf.fit(X[train_idx], y[train_idx])
                if hasattr(clf, "predict_proba"):
                    probs = clf.predict_proba(X[test_idx])[:, 1]
                else:
                    probs = clf.decision_function(X[test_idx])
                if len(set(y[test_idx])) == 2:
                    aucs.append(roc_auc_score(y[test_idx], probs))
                    pr_aucs.append(average_precision_score(y[test_idx], probs))

            key = f"{clf_name}_{cv_name}"
            mean_auc = float(np.mean(aucs)) if aucs else float("nan")
            mean_pr = float(np.mean(pr_aucs)) if pr_aucs else float("nan")
            results[key] = {"auc_roc": mean_auc, "pr_auc": mean_pr}
            print(f"  {clf_name:25s} + {cv_name:25s}: AUC={mean_auc:.4f}, PR-AUC={mean_pr:.4f}")

    # Random baseline
    results["random_baseline"] = {"auc_roc": 0.5, "pr_auc": float(y.mean())}

    ALL_RESULTS["E5_prime"] = results
    print("  Done.")


# ============================================================
# E6': Family Positioning + Intra-Family Analysis
# ============================================================

def run_e6_prime(esb_df):
    """E6': Place ESIBank P450 in multi-family comparison."""
    print("\n" + "=" * 70)
    print("  E6': Family Positioning + Intra-Family Analysis")
    print("=" * 70)

    results = {}

    # --- Load existing E6 results ---
    with open(E6_RESULTS_JSON, "r") as f:
        e6_results = json.load(f)

    esb_auc = roc_auc_score(esb_df["Label"], esb_df["logit"])
    esb_aupr = average_precision_score(esb_df["Label"], esb_df["logit"])

    # --- Analysis 1: Family attribute comparison ---
    results["analysis_1_attributes"] = {
        "esibank_p450": {
            "n_enzymes": int(esb_df["Enzyme Index"].nunique()),
            "n_pairs": len(esb_df),
            "n_pos": int(esb_df["Label"].sum()),
            "n_neg": int((esb_df["Label"] == 0).sum()),
            "auc_roc": float(esb_auc),
            "auc_pr": float(esb_aupr),
        },
        "our_p450": {
            "auc_roc": float(e6_results.get("P450", {}).get("auc_roc", 0.517)),
        },
    }
    print(f"  ESIBank P450 AUC: {esb_auc:.4f}")

    # --- Analysis 2: Intra-P450 promiscuity stratification ---
    # Compute per-enzyme positive count
    enz_pos = esb_df[esb_df["Label"] == 1].groupby("Enzyme Index").size()
    enz_all = esb_df.groupby("Enzyme Index").size()

    low_prom_enzymes = set(enz_pos[enz_pos <= 1].index) | (set(enz_all.index) - set(enz_pos.index))
    high_prom_enzymes = set(enz_pos[enz_pos > 1].index)

    low_df = esb_df[esb_df["Enzyme Index"].isin(low_prom_enzymes)]
    high_df = esb_df[esb_df["Enzyme Index"].isin(high_prom_enzymes)]

    low_auc = roc_auc_score(low_df["Label"], low_df["logit"]) if low_df["Label"].nunique() == 2 else float("nan")
    high_auc = roc_auc_score(high_df["Label"], high_df["logit"]) if len(high_df) > 0 and high_df["Label"].nunique() == 2 else float("nan")

    results["analysis_2_stratification"] = {
        "low_promiscuity": {"n_enzymes": len(low_prom_enzymes), "auc": float(low_auc)},
        "high_promiscuity": {"n_enzymes": len(high_prom_enzymes), "auc": float(high_auc)},
    }
    print(f"  Low promiscuity (<=1 substrate): {len(low_prom_enzymes)} enzymes, AUC={low_auc:.4f}")
    print(f"  High promiscuity (>1 substrate): {len(high_prom_enzymes)} enzymes, AUC={high_auc:.4f}")

    # --- Analysis 3: Score gap comparison ---
    esb_pos_mean = esb_df[esb_df["Label"] == 1]["logit"].mean()
    esb_neg_mean = esb_df[esb_df["Label"] == 0]["logit"].mean()
    esb_gap = esb_pos_mean - esb_neg_mean

    results["analysis_3_score_gap"] = {
        "esibank_p450_gap": float(esb_gap),
        "esibank_pos_mean": float(esb_pos_mean),
        "esibank_neg_mean": float(esb_neg_mean),
    }
    print(f"  Score gap: {esb_gap:.3f} (pos={esb_pos_mean:.3f}, neg={esb_neg_mean:.3f})")

    # --- Analysis D': Updated 8-family table ---
    family_table = {}
    for fam_name, fam_data in e6_results.items():
        family_table[fam_name] = {
            "auc_roc": fam_data.get("auc_roc"),
            "promiscuity_mean": fam_data.get("promiscuity_mean"),
            "n_valid": fam_data.get("n_valid_predictions"),
            "distribution": "in-dist" if fam_name != "P450" else "enzyme-no-overlap+pipeline-mismatch",
        }

    # Compute ESIBank P450 promiscuity from brenda data.csv
    try:
        brenda_data = pd.read_csv(G_BRENDA / "data.csv")
        brenda_pos = brenda_data[brenda_data["label"] == 1]
        # Filter to P450 enzymes
        enzyme_mapping = pd.read_csv(E7_DIR / "esibank_p450_enzyme_mapping.csv")
        p450_indices = set(enzyme_mapping["enzyme_index"].values)
        p450_pos = brenda_pos[brenda_pos["enzyme"].isin(p450_indices)]
        esb_promiscuity = p450_pos.groupby("reaction")["enzyme"].nunique().mean()
    except Exception:
        esb_promiscuity = 2.25  # fallback

    family_table["P450_ESIBank"] = {
        "auc_roc": float(esb_auc),
        "promiscuity_mean": float(esb_promiscuity),
        "n_valid": len(esb_df),
        "distribution": "in-dist",
    }

    # Rename original P450 for clarity
    if "P450" in family_table:
        family_table["P450_independent"] = family_table.pop("P450")

    results["analysis_D_family_table"] = family_table
    print(f"  ESIBank P450 promiscuity: {esb_promiscuity:.2f}")

    # --- Analysis E': Level 1 correlation (8 points) ---
    points = []
    for fam_name, fam_info in family_table.items():
        if fam_info["auc_roc"] is not None and fam_info["promiscuity_mean"] is not None:
            points.append({
                "family": fam_name,
                "promiscuity": fam_info["promiscuity_mean"],
                "auc": fam_info["auc_roc"],
                "distribution": fam_info["distribution"],
            })
    points_df = pd.DataFrame(points)

    if len(points_df) >= 4:
        sp = stats.spearmanr(points_df["promiscuity"], points_df["auc"])
        results["analysis_E_level1"] = {
            "n_points": len(points_df),
            "spearman_rho": float(sp.statistic) if hasattr(sp, 'statistic') else float(sp[0]),
            "spearman_p": float(sp.pvalue) if hasattr(sp, 'pvalue') else float(sp[1]),
        }
        rho_val = float(sp.statistic) if hasattr(sp, 'statistic') else float(sp[0])
        p_val = float(sp.pvalue) if hasattr(sp, 'pvalue') else float(sp[1])
        print(f"  Level 1 Spearman (n={len(points_df)}): rho={rho_val:.3f}, p={p_val:.3f}")

        # In-dist only
        in_dist = points_df[points_df["distribution"] == "in-dist"]
        if len(in_dist) >= 4:
            sp_in = stats.spearmanr(in_dist["promiscuity"], in_dist["auc"])
            rho_in = float(sp_in.statistic) if hasattr(sp_in, 'statistic') else float(sp_in[0])
            p_in = float(sp_in.pvalue) if hasattr(sp_in, 'pvalue') else float(sp_in[1])
            results["analysis_E_level1"]["in_dist_only"] = {
                "n_points": len(in_dist),
                "spearman_rho": rho_in,
                "spearman_p": p_in,
            }
            print(f"  Level 1 (in-dist only, n={len(in_dist)}): rho={rho_in:.3f}, p={p_in:.3f}")

    # --- Analysis F': Level 2 intra-family binning ---
    # Compute per-substrate promiscuity from FULL brenda data (not just test predictions)
    try:
        brenda_full = pd.read_csv(G_BRENDA / "data.csv")
        enzyme_mapping_df = pd.read_csv(E7_DIR / "esibank_p450_enzyme_mapping.csv")
        p450_idx_set = set(enzyme_mapping_df["enzyme_index"].values)
        brenda_p450 = brenda_full[brenda_full["enzyme"].isin(p450_idx_set)]
        brenda_p450_pos = brenda_p450[brenda_p450["label"] == 1]
        per_sub_prom = brenda_p450_pos.groupby("reaction")["enzyme"].nunique().reset_index()
        per_sub_prom.columns = ["Substrate Index", "promiscuity"]
        print(f"  Level 2: promiscuity from full data ({len(per_sub_prom)} substrates)")
    except Exception as e:
        print(f"  WARNING: Could not load full data for promiscuity, using predictions: {e}")
        per_sub_prom = esb_df[esb_df["Label"] == 1].groupby("Substrate Index")["Enzyme Index"].nunique().reset_index()
        per_sub_prom.columns = ["Substrate Index", "promiscuity"]

    per_sub_auc = []
    for sub_idx, grp in esb_df.groupby("Substrate Index"):
        n_pos = (grp["Label"] == 1).sum()
        n_neg = (grp["Label"] == 0).sum()
        if n_pos >= 2 and n_neg >= 2:
            auc = roc_auc_score(grp["Label"], grp["logit"])
            per_sub_auc.append({"Substrate Index": sub_idx, "auc": auc, "n_pos": n_pos, "n_neg": n_neg})
    per_sub_auc_df = pd.DataFrame(per_sub_auc)

    if len(per_sub_auc_df) > 0:
        per_sub_auc_df = per_sub_auc_df.merge(per_sub_prom, on="Substrate Index", how="left")
        per_sub_auc_df["promiscuity"] = per_sub_auc_df["promiscuity"].fillna(0)

        # Log-scale bins: 1, 2, 3-4, 5-8, 9-16, ...
        bins_edges = [0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 9999]
        bin_labels = ["1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129-256", "257+"]
        per_sub_auc_df["prom_bin"] = pd.cut(per_sub_auc_df["promiscuity"], bins=bins_edges, labels=bin_labels, right=True)

        level2_rows = []
        for b in bin_labels:
            bdata = per_sub_auc_df[per_sub_auc_df["prom_bin"] == b]
            if len(bdata) == 0:
                continue
            total_pos = bdata["n_pos"].sum()
            total_neg = bdata["n_neg"].sum()
            if total_pos >= 1 and total_neg >= 1:
                # Compute AUC on pooled data
                pooled = esb_df[esb_df["Substrate Index"].isin(bdata["Substrate Index"])]
                if pooled["Label"].nunique() == 2:
                    bin_auc = roc_auc_score(pooled["Label"], pooled["logit"])
                    confidence = "high" if total_pos >= 10 and total_neg >= 10 else "low"
                    level2_rows.append({
                        "bin": b, "n_pos": int(total_pos), "n_neg": int(total_neg),
                        "n_substrates": len(bdata), "auc": float(bin_auc), "confidence": confidence,
                    })

        results["analysis_F_level2"] = level2_rows
        for row in level2_rows:
            flag = " *" if row["confidence"] == "low" else ""
            print(f"    [{row['bin']:>7}] pos={row['n_pos']:>4} neg={row['n_neg']:>5} subs={row['n_substrates']:>3} AUC={row['auc']:.3f}{flag}")

    # --- Plot: 8-family scatter ---
    fig, ax = plt.subplots(figsize=(10, 7))
    for _, row in points_df.iterrows():
        if row["family"] == "P450_independent":
            color, marker, size = "crimson", "D", 150
        elif row["family"] == "P450_ESIBank":
            color, marker, size = "darkorange", "*", 200
        else:
            color, marker, size = "steelblue", "o", 100
        ax.scatter(row["promiscuity"], row["auc"], s=size, c=color, marker=marker, zorder=5, edgecolors="black", linewidths=0.5)
        offset = (8, 8) if row["family"] != "P450_independent" else (8, -12)
        ax.annotate(row["family"], (row["promiscuity"], row["auc"]),
                    textcoords="offset points", xytext=offset, fontsize=8,
                    fontweight="bold" if "P450" in row["family"] else "normal")
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel("Mean Substrate Promiscuity")
    ax.set_ylabel("AUC-ROC")
    ax.set_title("E6': 8-Family AUC vs Promiscuity (ESIBank P450 added)")
    ax.set_xscale("symlog", linthresh=1)

    plt.tight_layout()
    plt.savefig(E7_DIR / "e6_prime_family_scatter.png", dpi=150)
    plt.close()

    ALL_RESULTS["E6_prime"] = results
    print("  Done.")


# ============================================================
# Main
# ============================================================

def main():
    os.makedirs(E7_DIR, exist_ok=True)
    sys.path.insert(0, str(PROJECT_ROOT / "src"))

    # Load predictions
    print("Loading predictions...")
    esb_df = pd.read_csv(ESIBANK_PRED)
    our_df = pd.read_csv(OUR_PRED)
    print(f"  ESIBank P450: {len(esb_df)} rows (pos={esb_df['Label'].sum()}, neg={(esb_df['Label']==0).sum()})")
    print(f"  Our P450: {len(our_df)} rows (pos={our_df['Label'].sum()}, neg={(our_df['Label']==0).sum()})")

    esb_auc = roc_auc_score(esb_df["Label"], esb_df["logit"])
    our_auc = roc_auc_score(our_df["Label"], our_df["logit"])
    print(f"  ESIBank AUC: {esb_auc:.4f}, Our AUC: {our_auc:.4f}")

    # Run all experiments
    run_e1_prime(esb_df, our_df)
    run_e2_prime(esb_df, our_df)
    run_e3_prime(esb_df)
    run_e4_prime(esb_df)
    run_e5_prime(esb_df)
    run_e6_prime(esb_df)

    # --- Comparison summary ---
    print("\n" + "=" * 70)
    print("  COMPARISON SUMMARY: ESIBank P450 vs Our P450")
    print("=" * 70)

    e4_esb = ALL_RESULTS.get("E4_prime", {})
    print(f"  {'Metric':<30} {'ESIBank P450':>15} {'Our P450':>15} {'Delta':>10}")
    print(f"  {'-'*70}")
    comparisons = [
        ("AUC-ROC", f"{esb_auc:.4f}", f"{our_auc:.4f}", f"{esb_auc-our_auc:+.4f}"),
        ("Score gap", f"{ALL_RESULTS.get('E1_prime',{}).get('analysis_A',{}).get('esibank',{}).get('gap','?'):.3f}" if isinstance(ALL_RESULTS.get('E1_prime',{}).get('analysis_A',{}).get('esibank',{}).get('gap'), (int,float)) else "?",
         "0.260", ""),
        ("R2(enzyme)", f"{e4_esb.get('r2_enzyme','?'):.4f}" if isinstance(e4_esb.get('r2_enzyme'), (int,float)) else "?",
         "-0.0600", ""),
        ("R2(substrate)", f"{e4_esb.get('r2_substrate','?'):.4f}" if isinstance(e4_esb.get('r2_substrate'), (int,float)) else "?",
         "0.3700", ""),
        ("Residual AUC", f"{e4_esb.get('auc_residual','?'):.4f}" if isinstance(e4_esb.get('auc_residual'), (int,float)) else "?",
         "0.5090", ""),
    ]
    for name, esb_val, our_val, delta in comparisons:
        print(f"  {name:<30} {esb_val:>15} {our_val:>15} {delta:>10}")

    # Save all results
    with open(E7_DIR / "e7_analysis_results.json", "w", encoding="utf-8") as f:
        json.dump(_safe_json(ALL_RESULTS), f, indent=2, ensure_ascii=False)
    print(f"\n  All results saved to {E7_DIR / 'e7_analysis_results.json'}")


if __name__ == "__main__":
    main()
