"""
E1: Label Semantic Probe
Tests: C2 (label semantic mismatch) — does the model learn binding affinity or catalytic activity?
Method: Compare negative sample logit distributions between ABL-04 (Vina+inhibitor neg) and ABL-01 (Vina+random neg).
Both use Vina structures, so logit difference is purely from negative sample identity.
Also: Chemical space analysis of substrates vs inhibitors using Morgan fingerprints.
"""
import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy import stats
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
ABL04_PATH = os.path.join(PATHB, "results", "06_Step6_消融实验", "ABL-04_vina_inhibitor", "predictions.csv")
ABL01_DIR = os.path.join(PATHB, "results", "06_Step6_消融实验", "ABL-01_vina_1to1")
SUBSTRATES_PATH = os.path.join(PATHB, "data", "00_shared", "datasets", "B6_v1", "Substrates.csv")
DATA_PATH = os.path.join(PATHB, "data", "00_shared", "datasets", "B6_v1", "data.csv")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E1_标签语义探测")


def cliffs_delta(x, y):
    """Compute Cliff's delta effect size."""
    n_x, n_y = len(x), len(y)
    more = sum(1 for xi in x for yi in y if xi > yi)
    less = sum(1 for xi in x for yi in y if xi < yi)
    return (more - less) / (n_x * n_y)


def run_logit_comparison():
    """Part 1: Compare inhibitor vs random negative logit distributions."""
    print("=" * 60)
    print("E1 Part 1: Logit Distribution Comparison")
    print("=" * 60)

    # Load ABL-04 (Vina + inhibitor negatives)
    df04 = pd.read_csv(ABL04_PATH)
    neg04 = df04[df04["Label"] == 0]["logit"].values
    print(f"ABL-04 negatives (inhibitors): n={len(neg04)}, mean={neg04.mean():.4f}, std={neg04.std():.4f}")

    # Load all 5 ABL-01 seeds (Vina + random negatives)
    seed_files = sorted([f for f in os.listdir(ABL01_DIR) if f.startswith("predictions_seed")])
    all_neg01 = {}
    for sf in seed_files:
        seed = sf.replace("predictions_seed", "").replace(".csv", "")
        df01 = pd.read_csv(os.path.join(ABL01_DIR, sf))
        all_neg01[seed] = df01[df01["Label"] == 0]["logit"].values
        print(f"ABL-01 seed {seed} negatives (random): n={len(all_neg01[seed])}, "
              f"mean={all_neg01[seed].mean():.4f}, std={all_neg01[seed].std():.4f}")

    # Aggregate random negatives across seeds for main comparison
    neg01_concat = np.concatenate(list(all_neg01.values()))

    # Statistical tests
    ks_stat, ks_p = stats.ks_2samp(neg04, neg01_concat)
    mwu_stat, mwu_p = stats.mannwhitneyu(neg04, neg01_concat, alternative="two-sided")
    cd = cliffs_delta(neg04, neg01_concat)

    # AUC: Can model scores distinguish inhibitor neg from random neg?
    labels = np.concatenate([np.ones(len(neg04)), np.zeros(len(neg01_concat))])
    scores = np.concatenate([neg04, neg01_concat])
    dist_auc = roc_auc_score(labels, scores)

    # Per-seed AUC
    per_seed_aucs = []
    per_seed_ks = []
    for seed, neg_vals in all_neg01.items():
        labels_s = np.concatenate([np.ones(len(neg04)), np.zeros(len(neg_vals))])
        scores_s = np.concatenate([neg04, neg_vals])
        per_seed_aucs.append(roc_auc_score(labels_s, scores_s))
        ks_s, _ = stats.ks_2samp(neg04, neg_vals)
        per_seed_ks.append(ks_s)

    print(f"\n--- Statistical Tests ---")
    print(f"KS test:        stat={ks_stat:.4f}, p={ks_p:.2e}")
    print(f"Mann-Whitney U: stat={mwu_stat:.0f}, p={mwu_p:.2e}")
    print(f"Cliff's delta:  {cd:.4f}")
    print(f"Distinguishing AUC (inhib=1 vs random=0): {dist_auc:.4f}")
    print(f"  Per-seed AUCs: {[f'{a:.4f}' for a in per_seed_aucs]}")
    print(f"  Per-seed KS:   {[f'{k:.4f}' for k in per_seed_ks]}")

    # Interpretation
    if neg04.mean() < neg01_concat.mean():
        direction = "inhibitor negatives get LOWER scores than random negatives"
        c2_verdict = ("Model does NOT treat inhibitors as 'near-positive'. "
                      "Instead, inhibitors are more confidently classified as negative, "
                      "likely via molecular fingerprint shortcut (substrates vs drug-like inhibitors).")
    else:
        direction = "inhibitor negatives get HIGHER scores than random negatives"
        c2_verdict = ("Model treats inhibitors as closer to positives (binding affinity signal). "
                      "EXP01's 0.71 baseline may be inflated by this.")

    print(f"\nDirection: {direction}")
    print(f"C2 verdict: {c2_verdict}")

    # Visualization: Violin/density plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Density plot
    ax = axes[0]
    pos04 = df04[df04["Label"] == 1]["logit"].values
    ax.hist(pos04, bins=40, alpha=0.5, density=True, label=f"Positive (n={len(pos04)})", color="green")
    ax.hist(neg04, bins=40, alpha=0.5, density=True, label=f"Inhibitor neg (n={len(neg04)})", color="red")
    # Use first seed for random negatives in plot
    first_neg01 = list(all_neg01.values())[0]
    ax.hist(first_neg01, bins=40, alpha=0.5, density=True, label=f"Random neg (n={len(first_neg01)})", color="blue")
    ax.set_xlabel("Logit Score")
    ax.set_ylabel("Density")
    ax.set_title("Score Distributions by Sample Type")
    ax.legend()
    ax.axvline(x=0, color="gray", linestyle="--", alpha=0.5)

    # Box plot comparing all negative types
    ax = axes[1]
    data_box = [neg04, neg01_concat]
    bp = ax.boxplot(data_box, labels=["Inhibitor neg", "Random neg"], patch_artist=True)
    bp["boxes"][0].set_facecolor("red")
    bp["boxes"][0].set_alpha(0.5)
    bp["boxes"][1].set_facecolor("blue")
    bp["boxes"][1].set_alpha(0.5)
    ax.set_ylabel("Logit Score")
    ax.set_title(f"Negative Sample Logit Comparison\nKS={ks_stat:.3f}, Cliff's d={cd:.3f}")

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e1_logit_comparison.png"), dpi=150)
    plt.close()
    print(f"\nPlot saved: e1_logit_comparison.png")

    return {
        "neg04_mean": float(neg04.mean()),
        "neg04_std": float(neg04.std()),
        "neg01_mean": float(neg01_concat.mean()),
        "neg01_std": float(neg01_concat.std()),
        "ks_stat": float(ks_stat),
        "ks_p": float(ks_p),
        "mwu_p": float(mwu_p),
        "cliffs_delta": float(cd),
        "distinguishing_auc": float(dist_auc),
        "per_seed_aucs": [float(a) for a in per_seed_aucs],
        "direction": direction,
        "c2_verdict": c2_verdict,
    }


def run_chemical_space_analysis():
    """Part 2: Chemical space comparison of substrates vs inhibitors."""
    print("\n" + "=" * 60)
    print("E1 Part 2: Chemical Space Analysis")
    print("=" * 60)

    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
        from sklearn.decomposition import PCA
    except ImportError:
        print("RDKit or sklearn not available, skipping chemical space analysis")
        return None

    # Load substrate SMILES and classify as substrate vs inhibitor
    subs_df = pd.read_csv(SUBSTRATES_PATH)
    data_df = pd.read_csv(DATA_PATH)

    pos_sub_idx = set(data_df[data_df["Label"] == 1]["Substrate Index"].unique())
    neg_sub_idx = set(data_df[data_df["Label"] == 0]["Substrate Index"].unique())
    only_pos = pos_sub_idx - neg_sub_idx
    only_neg = neg_sub_idx - pos_sub_idx

    print(f"Substrate-only indices: {len(only_pos)}")
    print(f"Inhibitor-only indices: {len(only_neg)}")

    # Generate Morgan fingerprints
    fps_data = []
    for _, row in subs_df.iterrows():
        idx = row["Substrate_Index"]
        smiles = row["Substrate_SMILES"]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
        fp_arr = np.array(fp)
        label = "substrate" if idx in only_pos else ("inhibitor" if idx in only_neg else "both")
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        fps_data.append({"idx": idx, "label": label, "fp": fp_arr, "MW": mw, "LogP": logp})

    fps_df = pd.DataFrame([{k: v for k, v in d.items() if k != "fp"} for d in fps_data])
    fp_matrix = np.array([d["fp"] for d in fps_data])

    print(f"Total molecules with valid fingerprints: {len(fps_data)}")
    print(f"  Substrates: {sum(1 for d in fps_data if d['label'] == 'substrate')}")
    print(f"  Inhibitors: {sum(1 for d in fps_data if d['label'] == 'inhibitor')}")
    print(f"  Both:       {sum(1 for d in fps_data if d['label'] == 'both')}")

    # Property comparison
    subs_mw = [d["MW"] for d in fps_data if d["label"] == "substrate"]
    inhib_mw = [d["MW"] for d in fps_data if d["label"] == "inhibitor"]
    subs_logp = [d["LogP"] for d in fps_data if d["label"] == "substrate"]
    inhib_logp = [d["LogP"] for d in fps_data if d["label"] == "inhibitor"]

    print(f"\nProperty comparison:")
    print(f"  Substrates MW:  mean={np.mean(subs_mw):.1f}, std={np.std(subs_mw):.1f}")
    print(f"  Inhibitors MW:  mean={np.mean(inhib_mw):.1f}, std={np.std(inhib_mw):.1f}")
    print(f"  Substrates LogP: mean={np.mean(subs_logp):.2f}, std={np.std(subs_logp):.2f}")
    print(f"  Inhibitors LogP: mean={np.mean(inhib_logp):.2f}, std={np.std(inhib_logp):.2f}")

    # PCA on Morgan fingerprints
    pca = PCA(n_components=2)
    fp_2d = pca.fit_transform(fp_matrix)

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # PCA scatter
    ax = axes[0]
    colors = {"substrate": "green", "inhibitor": "red", "both": "orange"}
    for label_type in ["substrate", "inhibitor", "both"]:
        mask = [d["label"] == label_type for d in fps_data]
        if sum(mask) == 0:
            continue
        pts = fp_2d[mask]
        ax.scatter(pts[:, 0], pts[:, 1], c=colors[label_type], alpha=0.6,
                   s=30, label=f"{label_type} (n={sum(mask)})", edgecolors="white", linewidths=0.3)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")
    ax.set_title("Chemical Space (Morgan FP PCA)")
    ax.legend()

    # Property distributions
    ax = axes[1]
    ax.scatter(subs_mw, subs_logp, c="green", alpha=0.5, s=30, label=f"Substrates (n={len(subs_mw)})")
    ax.scatter(inhib_mw, inhib_logp, c="red", alpha=0.5, s=30, label=f"Inhibitors (n={len(inhib_mw)})")
    ax.set_xlabel("Molecular Weight")
    ax.set_ylabel("LogP")
    ax.set_title("Physicochemical Properties")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e1_chemical_space.png"), dpi=150)
    plt.close()
    print(f"Plot saved: e1_chemical_space.png")

    # Tanimoto similarity within and between groups
    from rdkit import DataStructs
    subs_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(subs_df.loc[subs_df.Substrate_Index == idx, "Substrate_SMILES"].iloc[0]), 2, nBits=1024)
                for idx in list(only_pos)[:50] if Chem.MolFromSmiles(subs_df.loc[subs_df.Substrate_Index == idx, "Substrate_SMILES"].iloc[0]) is not None]
    inhib_fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(subs_df.loc[subs_df.Substrate_Index == idx, "Substrate_SMILES"].iloc[0]), 2, nBits=1024)
                 for idx in list(only_neg)[:50] if Chem.MolFromSmiles(subs_df.loc[subs_df.Substrate_Index == idx, "Substrate_SMILES"].iloc[0]) is not None]

    within_subs = [DataStructs.TanimotoSimilarity(subs_fps[i], subs_fps[j])
                   for i in range(len(subs_fps)) for j in range(i+1, len(subs_fps))]
    within_inhib = [DataStructs.TanimotoSimilarity(inhib_fps[i], inhib_fps[j])
                    for i in range(len(inhib_fps)) for j in range(i+1, len(inhib_fps))]
    between = [DataStructs.TanimotoSimilarity(subs_fps[i], inhib_fps[j])
               for i in range(len(subs_fps)) for j in range(len(inhib_fps))]

    print(f"\nTanimoto similarity:")
    print(f"  Within substrates:  mean={np.mean(within_subs):.4f}")
    print(f"  Within inhibitors:  mean={np.mean(within_inhib):.4f}")
    print(f"  Between sub-inhib:  mean={np.mean(between):.4f}")

    return {
        "n_substrates": len(only_pos),
        "n_inhibitors": len(only_neg),
        "substrates_mw_mean": float(np.mean(subs_mw)),
        "inhibitors_mw_mean": float(np.mean(inhib_mw)),
        "substrates_logp_mean": float(np.mean(subs_logp)),
        "inhibitors_logp_mean": float(np.mean(inhib_logp)),
        "tanimoto_within_substrates": float(np.mean(within_subs)),
        "tanimoto_within_inhibitors": float(np.mean(within_inhib)),
        "tanimoto_between": float(np.mean(between)),
        "pca_var_explained": [float(v) for v in pca.explained_variance_ratio_[:2]],
    }


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    logit_results = run_logit_comparison()
    chem_results = run_chemical_space_analysis()

    # Combined results
    results = {
        "experiment": "E1",
        "name": "Label Semantic Probe",
        "tests": "C2 — does model learn binding or catalysis?",
        "logit_comparison": logit_results,
        "chemical_space": chem_results,
    }

    with open(os.path.join(OUTPUT_DIR, "e1_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nAll results saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
