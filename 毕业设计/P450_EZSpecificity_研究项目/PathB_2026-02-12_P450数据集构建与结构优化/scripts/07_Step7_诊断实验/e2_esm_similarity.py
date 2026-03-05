"""
E2: ESM Embedding Similarity Analysis
Tests: M1 (sequence channel compression) — How similar are P450 ESM-2 embeddings?
Method: Read ESM-2 embeddings from LMDB, mean-pool to per-enzyme vectors, compute pairwise cosine similarity.
Control: Generate ESM-2 embeddings for toy_example enzymes (nitrilases) as non-P450 baseline.
"""
import lmdb
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy import stats
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
ENZYME_LMDB = os.path.join(PATHB, "data", "00_shared", "features", "enzyme_features.lmdb")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E2_ESM嵌入相似度")


def load_enzyme_embeddings(lmdb_path, indices=None):
    """Load ESM-2 embeddings from LMDB and mean-pool to per-enzyme vectors."""
    env = lmdb.open(lmdb_path, subdir=False, readonly=True, lock=False)
    embeddings = {}
    with env.begin() as txn:
        cursor = txn.cursor()
        for key, value in cursor:
            idx = int(key)
            if indices is not None and idx not in indices:
                continue
            data = pickle.loads(value)
            emb = data["embedding"]  # shape: (seq_len, 1280)
            mean_emb = emb.mean(axis=0)  # mean pooling -> (1280,)
            embeddings[idx] = mean_emb
    env.close()
    return embeddings


def generate_control_embeddings():
    """Generate ESM-2 embeddings for toy_example enzymes (non-P450) as control group."""
    import pandas as pd
    import sys
    sys.path.insert(0, r"D:\EZSpecificity_Project\src")

    toy_enzymes = pd.read_csv(r"D:\EZSpecificity_Project\data\toy_example\Enzymes.csv")
    sequences = toy_enzymes["Protein sequence"].tolist()
    print(f"Generating ESM-2 embeddings for {len(sequences)} control enzymes...")

    try:
        import torch
        import esm
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        batch_converter = alphabet.get_batch_converter()
        model.eval()
        if torch.cuda.is_available():
            model = model.cuda()

        embeddings = {}
        for i, seq in enumerate(sequences):
            if len(seq) > 1000:
                continue
            data = [("protein", seq)]
            _, _, tokens = batch_converter(data)
            if torch.cuda.is_available():
                tokens = tokens.cuda()
            with torch.no_grad():
                results = model(tokens, repr_layers=[33])
            emb = results["representations"][33][0, 1:-1, :].cpu().numpy()  # remove BOS/EOS
            embeddings[i] = emb.mean(axis=0)
            if (i + 1) % 5 == 0:
                print(f"  {i + 1}/{len(sequences)} done")
        print(f"Generated {len(embeddings)} control embeddings")
        return embeddings
    except Exception as e:
        print(f"Failed to generate control embeddings: {e}")
        return None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load P450 embeddings
    print("Loading P450 ESM-2 embeddings...")
    p450_emb = load_enzyme_embeddings(ENZYME_LMDB)
    print(f"Loaded {len(p450_emb)} P450 enzyme embeddings")

    # Convert to matrix
    p450_indices = sorted(p450_emb.keys())
    p450_matrix = np.array([p450_emb[i] for i in p450_indices])

    # L2 normalize
    norms = np.linalg.norm(p450_matrix, axis=1, keepdims=True)
    p450_normed = p450_matrix / norms

    # Pairwise cosine similarity
    cos_sim = cosine_similarity(p450_normed)
    n = len(cos_sim)
    upper_tri = cos_sim[np.triu_indices(n, k=1)]

    print(f"\nP450 pairwise cosine similarity ({n} enzymes, {len(upper_tri)} pairs):")
    print(f"  Mean:   {upper_tri.mean():.4f}")
    print(f"  Std:    {upper_tri.std():.4f}")
    print(f"  Median: {np.median(upper_tri):.4f}")
    print(f"  Min:    {upper_tri.min():.4f}")
    print(f"  Max:    {upper_tri.max():.4f}")
    print(f"  Q25:    {np.percentile(upper_tri, 25):.4f}")
    print(f"  Q75:    {np.percentile(upper_tri, 75):.4f}")
    print(f"  > 0.90: {(upper_tri > 0.9).sum()} pairs ({(upper_tri > 0.9).mean():.1%})")
    print(f"  > 0.95: {(upper_tri > 0.95).sum()} pairs ({(upper_tri > 0.95).mean():.1%})")

    # Effective dimensionality / concentration
    embedding_dim = p450_matrix.shape[1]
    eigenvalues = np.linalg.eigvalsh(np.cov(p450_normed.T))[::-1]
    cum_var = np.cumsum(eigenvalues) / eigenvalues.sum()
    n_90 = np.searchsorted(cum_var, 0.9) + 1
    n_95 = np.searchsorted(cum_var, 0.95) + 1
    print(f"\nDimensionality analysis:")
    print(f"  Embedding dim: {embedding_dim}")
    print(f"  PCA dims for 90% var: {n_90}")
    print(f"  PCA dims for 95% var: {n_95}")

    # Control group: toy_example enzymes (non-P450)
    ctrl_emb = generate_control_embeddings()
    ctrl_stats = None
    if ctrl_emb and len(ctrl_emb) > 5:
        ctrl_matrix = np.array([ctrl_emb[i] for i in sorted(ctrl_emb.keys())])
        ctrl_norms = np.linalg.norm(ctrl_matrix, axis=1, keepdims=True)
        ctrl_normed = ctrl_matrix / ctrl_norms
        ctrl_cos = cosine_similarity(ctrl_normed)
        n_ctrl = len(ctrl_cos)
        ctrl_upper = ctrl_cos[np.triu_indices(n_ctrl, k=1)]

        print(f"\nControl (non-P450) pairwise cosine ({n_ctrl} enzymes, {len(ctrl_upper)} pairs):")
        print(f"  Mean:   {ctrl_upper.mean():.4f}")
        print(f"  Std:    {ctrl_upper.std():.4f}")
        print(f"  Min:    {ctrl_upper.min():.4f}")
        print(f"  Max:    {ctrl_upper.max():.4f}")

        # Cross-group similarity
        cross_cos = cosine_similarity(p450_normed, ctrl_normed)
        print(f"\nCross-group (P450 vs control):")
        print(f"  Mean:   {cross_cos.mean():.4f}")
        print(f"  Std:    {cross_cos.std():.4f}")

        # KS test
        ks_stat, ks_p = stats.ks_2samp(upper_tri, ctrl_upper)
        print(f"\nKS test (P450 vs control): stat={ks_stat:.4f}, p={ks_p:.2e}")

        ctrl_stats = {
            "n_control": n_ctrl,
            "ctrl_cos_mean": float(ctrl_upper.mean()),
            "ctrl_cos_std": float(ctrl_upper.std()),
            "cross_cos_mean": float(cross_cos.mean()),
            "ks_stat": float(ks_stat),
            "ks_p": float(ks_p),
        }

    # M1 verdict
    p450_mean_cos = upper_tri.mean()
    if p450_mean_cos > 0.95:
        m1_verdict = "M1 STRONGLY CONFIRMED: P450 ESM embeddings are extremely compressed (mean cosine > 0.95). Sequence channel is effectively useless."
    elif p450_mean_cos > 0.90:
        m1_verdict = "M1 CONFIRMED: P450 ESM embeddings are highly similar (mean cosine > 0.90). Sequence channel has very limited discriminative power."
    elif p450_mean_cos > 0.80:
        m1_verdict = "M1 PARTIALLY SUPPORTED: P450 ESM embeddings are moderately similar (mean cosine > 0.80). Sequence channel may retain some signal."
    else:
        m1_verdict = "M1 NOT CONFIRMED: P450 ESM embeddings show sufficient diversity (mean cosine < 0.80). Sequence channel should provide discriminative signal."

    if ctrl_stats:
        delta = p450_mean_cos - ctrl_stats["ctrl_cos_mean"]
        m1_verdict += f" Delta vs control = {delta:+.4f}."

    print(f"\nM1 verdict: {m1_verdict}")

    # Visualization
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Histogram of pairwise cosine similarity
    ax = axes[0]
    ax.hist(upper_tri, bins=50, alpha=0.7, color="steelblue", density=True, label=f"P450 (n={n})")
    if ctrl_stats:
        ax.hist(ctrl_upper, bins=50, alpha=0.5, color="coral", density=True, label=f"Control (n={n_ctrl})")
    ax.axvline(x=upper_tri.mean(), color="blue", linestyle="--", alpha=0.8)
    ax.set_xlabel("Cosine Similarity")
    ax.set_ylabel("Density")
    ax.set_title("Pairwise Cosine Similarity Distribution")
    ax.legend()

    # Heatmap (subset for visibility)
    ax = axes[1]
    n_show = min(50, n)
    im = ax.imshow(cos_sim[:n_show, :n_show], cmap="RdYlBu_r", vmin=0.5, vmax=1.0)
    ax.set_title(f"Cosine Similarity Heatmap (first {n_show})")
    ax.set_xlabel("Enzyme Index")
    ax.set_ylabel("Enzyme Index")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # PCA cumulative variance
    ax = axes[2]
    ax.plot(range(1, min(51, len(cum_var) + 1)), cum_var[:50], "b-o", markersize=3)
    ax.axhline(y=0.9, color="red", linestyle="--", label=f"90% at {n_90} dims")
    ax.axhline(y=0.95, color="orange", linestyle="--", label=f"95% at {n_95} dims")
    ax.set_xlabel("Number of PCA Components")
    ax.set_ylabel("Cumulative Variance Explained")
    ax.set_title("Effective Dimensionality")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e2_esm_similarity.png"), dpi=150)
    plt.close()

    # Save results
    results = {
        "experiment": "E2",
        "name": "ESM Embedding Similarity",
        "tests": "M1 - sequence channel compression",
        "n_enzymes": n,
        "n_pairs": len(upper_tri),
        "cos_mean": float(upper_tri.mean()),
        "cos_std": float(upper_tri.std()),
        "cos_median": float(np.median(upper_tri)),
        "cos_min": float(upper_tri.min()),
        "cos_max": float(upper_tri.max()),
        "cos_q25": float(np.percentile(upper_tri, 25)),
        "cos_q75": float(np.percentile(upper_tri, 75)),
        "pct_above_090": float((upper_tri > 0.9).mean()),
        "pct_above_095": float((upper_tri > 0.95).mean()),
        "pca_dims_90pct": int(n_90),
        "pca_dims_95pct": int(n_95),
        "control": ctrl_stats,
        "m1_verdict": m1_verdict,
    }
    with open(os.path.join(OUTPUT_DIR, "e2_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
