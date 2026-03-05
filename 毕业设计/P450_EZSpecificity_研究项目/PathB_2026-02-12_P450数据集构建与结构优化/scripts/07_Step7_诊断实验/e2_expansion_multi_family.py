"""
E2 Expansion: Multi-Family ESM-2 Cosine Similarity Analysis
Computes intra-family and inter-family cosine similarity for 8 enzyme families.
Confirms whether ESM-2 anisotropy is universal or P450-specific.
"""
import lmdb
import pickle
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import gaussian_kde
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
LOCAL_LMDB = r"D:\EZSpecificity_Project\tmp_lmdb"
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E2_expansion_多家族ESM相似度")

# Use local copies for speed (copied from G: drive ESIBank/small_family)
# Halogenase excluded: 929MB LMDB too large for Google Drive streaming
FAMILIES = {
    "P450": os.path.join(PATHB, "data", "00_shared", "features", "enzyme_features.lmdb"),
    "Duf": os.path.join(LOCAL_LMDB, "Duf_enzyme_features.lmdb"),
    "Esterase": os.path.join(LOCAL_LMDB, "Esterase_enzyme_features.lmdb"),
    "Gt_acceptor": os.path.join(LOCAL_LMDB, "Gt_acceptor_enzyme_features.lmdb"),
    "Nitrilase": os.path.join(LOCAL_LMDB, "Nitrilase_enzyme_features.lmdb"),
    "Phosphatase": os.path.join(LOCAL_LMDB, "Phosphatase_enzyme_features.lmdb"),
    "Thiolase": os.path.join(LOCAL_LMDB, "Thiolase_enzyme_features.lmdb"),
}


def load_enzyme_embeddings(lmdb_path):
    """Load ESM-2 embeddings from LMDB, mean-pool to per-enzyme 1280-dim vectors."""
    env = lmdb.open(lmdb_path, subdir=False, readonly=True, lock=False, readahead=False)
    embeddings = {}
    with env.begin() as txn:
        for key, value in txn.cursor():
            idx = int(key)
            data = pickle.loads(value)
            emb = data["embedding"]  # (seq_len, 1280)
            embeddings[idx] = emb.mean(axis=0)  # mean pooling -> (1280,)
    env.close()
    return embeddings


def normalize_matrix(matrix):
    """L2-normalize rows with zero-norm safety."""
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms = np.clip(norms, 1e-12, None)
    return matrix / norms


def compute_intra_stats(cos_values):
    """Compute summary statistics for a cosine similarity distribution."""
    return {
        "mean": float(np.mean(cos_values)),
        "std": float(np.std(cos_values)),
        "median": float(np.median(cos_values)),
        "q25": float(np.percentile(cos_values, 25)),
        "q75": float(np.percentile(cos_values, 75)),
        "min": float(np.min(cos_values)),
        "max": float(np.max(cos_values)),
        "pct_above_090": float(np.mean(cos_values > 0.90)),
        "pct_above_095": float(np.mean(cos_values > 0.95)),
    }


def compute_pca_dims(normed_matrix):
    """Compute number of PCA dims for 90% and 95% variance explained."""
    n = normed_matrix.shape[0]
    if n < 3:
        return {"n90": -1, "n95": -1}
    cov = np.cov(normed_matrix.T)
    eigenvalues = np.linalg.eigvalsh(cov)[::-1]
    cum_var = np.cumsum(eigenvalues) / eigenvalues.sum()
    n90 = int(np.searchsorted(cum_var, 0.9) + 1)
    n95 = int(np.searchsorted(cum_var, 0.95) + 1)
    return {"n90": n90, "n95": n95}


def plot_ridge(family_results, output_path):
    """Ridge plot of intra-family cosine similarity distributions."""
    sorted_families = sorted(family_results.items(), key=lambda x: x[1]["stats"]["mean"], reverse=True)
    n_families = len(sorted_families)

    fig, axes = plt.subplots(n_families, 1, figsize=(9, 1.2 + 0.9 * n_families), sharex=True)
    if n_families == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=-0.3)

    cmap = plt.get_cmap("tab10")

    # Determine x range from data
    all_vals = np.concatenate([data["cosine_values"] for _, data in sorted_families])
    x_min = max(0.75, np.min(all_vals) - 0.02)
    x_max = min(1.0, np.max(all_vals) + 0.01)
    x_eval = np.linspace(x_min, x_max, 500)

    for i, (ax, (name, data)) in enumerate(zip(axes, sorted_families)):
        ax.set_zorder(n_families - i)
        ax.patch.set_alpha(0)

        cos_vals = data["cosine_values"]
        n_enz = data["n_enzymes"]
        is_p450 = name == "P450"

        kde = gaussian_kde(cos_vals)
        y_eval = kde(x_eval)

        color = "crimson" if is_p450 else cmap(i % 10)
        alpha_fill = 0.7 if is_p450 else 0.4
        fontweight = "bold" if is_p450 else "normal"

        ax.plot(x_eval, y_eval, color=color, lw=1.5)
        ax.fill_between(x_eval, 0, y_eval, color=color, alpha=alpha_fill)

        label = f"{name}\n(n={n_enz})"
        ax.set_ylabel(label, rotation=0, ha="right", va="center",
                       fontweight=fontweight, color="crimson" if is_p450 else "black", fontsize=9)
        ax.set_yticks([])

        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
        if i < n_families - 1:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", length=0)

        ax.axvline(0.90, color="gray", linestyle="--", alpha=0.5, zorder=0)
        ax.axvline(0.95, color="gray", linestyle="--", alpha=0.5, zorder=0)

    axes[-1].set_xlabel("ESM-2 Cosine Similarity", fontsize=11)
    axes[-1].set_xlim(x_min, x_max)
    fig.suptitle("Intra-family ESM-2 Cosine Similarity Distribution", fontsize=13, fontweight="bold", y=0.98)

    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def plot_heatmap(sim_matrix, family_names, output_path):
    """8x8 heatmap of inter/intra-family mean cosine similarity."""
    n = len(family_names)
    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.imshow(sim_matrix, cmap="YlOrRd", vmin=0.85, vmax=1.0, aspect="equal")
    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.04)
    cbar.ax.set_ylabel("Mean Cosine Similarity", rotation=-90, va="bottom", fontsize=11, labelpad=15)

    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(family_names, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_yticklabels(family_names)

    for i, name in enumerate(family_names):
        if name == "P450":
            ax.get_xticklabels()[i].set_fontweight("bold")
            ax.get_xticklabels()[i].set_color("crimson")
            ax.get_yticklabels()[i].set_fontweight("bold")
            ax.get_yticklabels()[i].set_color("crimson")

    for i in range(n):
        for j in range(n):
            val = sim_matrix[i, j]
            text_color = "white" if val > 0.94 else "black"
            text_weight = "bold" if i == j else "normal"
            ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                    color=text_color, fontweight=text_weight, fontsize=9)

    ax.tick_params(which="both", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title("Inter- and Intra-family ESM-2 Cosine Similarity", fontsize=13, fontweight="bold", pad=20)

    fig.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- Step 1: Load all families ---
    family_matrices = {}
    for name, path in FAMILIES.items():
        print(f"Loading {name} from {path}...")
        emb = load_enzyme_embeddings(path)
        indices = sorted(emb.keys())
        matrix = np.array([emb[i] for i in indices])
        normed = normalize_matrix(matrix)
        family_matrices[name] = normed
        print(f"  {name}: {len(indices)} enzymes, shape={matrix.shape}")

    # --- Step 2: Intra-family analysis ---
    family_results = {}
    print("\n=== Intra-family Cosine Similarity ===")
    print(f"{'Family':<15} {'N':>5} {'Mean':>7} {'Std':>7} {'>0.90':>7} {'>0.95':>7} {'PCA90':>6}")
    print("-" * 60)

    for name in FAMILIES:
        normed = family_matrices[name]
        n = normed.shape[0]
        cos_sim = cosine_similarity(normed)
        upper = cos_sim[np.triu_indices(n, k=1)]
        stats = compute_intra_stats(upper)
        pca = compute_pca_dims(normed)

        print(f"{name:<15} {n:>5} {stats['mean']:>7.4f} {stats['std']:>7.4f} "
              f"{stats['pct_above_090']:>6.1%} {stats['pct_above_095']:>6.1%} {pca['n90']:>6}")

        family_results[name] = {
            "n_enzymes": n,
            "n_pairs": len(upper),
            "stats": stats,
            "pca": pca,
            "cosine_values": upper,  # kept for plotting
        }

    # --- Step 3: Inter-family analysis (28 pairs) ---
    family_names = sorted(FAMILIES.keys())
    n_fam = len(family_names)
    sim_matrix = np.zeros((n_fam, n_fam))
    cross_details = {}

    print("\n=== Inter-family Cosine Similarity ===")
    for i in range(n_fam):
        for j in range(n_fam):
            ni = family_names[i]
            nj = family_names[j]
            if i == j:
                sim_matrix[i, j] = family_results[ni]["stats"]["mean"]
            else:
                cross = cosine_similarity(family_matrices[ni], family_matrices[nj])
                mean_cross = float(cross.mean())
                sim_matrix[i, j] = mean_cross
                if i < j:
                    cross_details[f"{ni}_vs_{nj}"] = {
                        "mean": mean_cross,
                        "median": float(np.median(cross)),
                        "std": float(np.std(cross)),
                    }
                    print(f"  {ni} vs {nj}: mean={mean_cross:.4f}")

    # --- Step 4: Visualization ---
    print("\nGenerating plots...")
    plot_ridge(family_results, os.path.join(OUTPUT_DIR, "e2_exp_ridge_plot.png"))
    plot_heatmap(sim_matrix, family_names, os.path.join(OUTPUT_DIR, "e2_exp_heatmap.png"))
    print("Plots saved.")

    # --- Step 5: Save results ---
    results = {
        "experiment": "E2_expansion",
        "name": "Multi-Family ESM-2 Cosine Similarity",
        "families": {},
        "cross_family": cross_details,
        "similarity_matrix": sim_matrix.tolist(),
        "family_order": family_names,
    }
    for name in FAMILIES:
        r = family_results[name]
        results["families"][name] = {
            "n_enzymes": r["n_enzymes"],
            "n_pairs": r["n_pairs"],
            **r["stats"],
            "pca_n90": r["pca"]["n90"],
            "pca_n95": r["pca"]["n95"],
        }

    with open(os.path.join(OUTPUT_DIR, "e2_expansion_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
