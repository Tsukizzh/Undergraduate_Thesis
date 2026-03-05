"""
E5 (New): Substrate Fingerprint Classifier Baseline
Tests: Can substrate Morgan fingerprints alone predict labels under Direction A?
Replaces flawed old E5 (per-substrate positive rate as score was tautological under Direction A).
Method: Train LogisticRegression + RandomForest on 1024-bit Morgan fingerprints, 2 CV modes.
"""
import pandas as pd
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, StratifiedGroupKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, average_precision_score
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PATHB = r"D:\EZSpecificity_Project\毕业设计\P450_EZSpecificity_研究项目\PathB_2026-02-12_P450数据集构建与结构优化"
PREDICTIONS_PATH = os.path.join(PATHB, "results", "05_Step5_重构评估", "predictions.csv")
MORGAN_PATH = os.path.join(PATHB, "data", "00_shared", "features", "morgan_fingerprint.npy")
OUTPUT_DIR = os.path.join(PATHB, "results", "07_Step7_Tier1_诊断实验", "E5_new_指纹分类器基线")


def run_cv(X, y, groups, clf_name, clf, cv, cv_name):
    """Run cross-validation and return per-fold AUC and AP scores."""
    roc_aucs, pr_aucs = [], []

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X, y, groups)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Skip fold if test has no positive or no negative
        if y_test.sum() == 0 or y_test.sum() == len(y_test):
            print(f"  {clf_name}/{cv_name} fold {fold_idx}: skipped (no pos/neg in test)")
            continue

        clf_instance = clf.__class__(**clf.get_params())
        clf_instance.fit(X_train, y_train)

        if hasattr(clf_instance, "predict_proba"):
            y_score = clf_instance.predict_proba(X_test)[:, 1]
        else:
            y_score = clf_instance.decision_function(X_test)

        roc = roc_auc_score(y_test, y_score)
        ap = average_precision_score(y_test, y_score)
        roc_aucs.append(roc)
        pr_aucs.append(ap)
        print(f"  {clf_name}/{cv_name} fold {fold_idx}: ROC-AUC={roc:.4f}, PR-AUC={ap:.4f} "
              f"(test: {len(y_test)} samples, {y_test.sum()} pos)")

    return roc_aucs, pr_aucs


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load data
    df = pd.read_csv(PREDICTIONS_PATH)
    morgan = np.load(MORGAN_PATH)
    print(f"Loaded {len(df)} predictions, morgan shape={morgan.shape}")

    # Validate indices
    sub_indices = df["Substrate Index"].values
    assert sub_indices.min() >= 0, "Negative substrate index"
    assert sub_indices.max() < morgan.shape[0], f"Substrate index {sub_indices.max()} >= morgan rows {morgan.shape[0]}"

    X = morgan[sub_indices].astype(np.float32)
    y = df["Label"].values
    groups = sub_indices

    n_pos = y.sum()
    n_neg = len(y) - n_pos
    prevalence = n_pos / len(y)
    print(f"Samples: {len(y)} ({n_pos} pos, {n_neg} neg, prevalence={prevalence:.4f})")
    print(f"Unique substrates: {len(np.unique(groups))}")
    print(f"Feature dim: {X.shape[1]}")

    # Define classifiers
    classifiers = {
        "LogisticRegression": LogisticRegression(
            C=1.0, penalty="l2", class_weight="balanced", max_iter=1000, random_state=42
        ),
        "RandomForest": RandomForestClassifier(
            n_estimators=100, max_depth=5, class_weight="balanced", random_state=42
        ),
    }

    # Define CV strategies
    cv_strategies = {
        "StratifiedKFold": StratifiedKFold(n_splits=5, shuffle=True, random_state=42),
        "StratifiedGroupKFold": StratifiedGroupKFold(n_splits=5),
    }

    # Run all combinations
    results = {
        "experiment": "E5_new",
        "name": "Substrate Fingerprint Classifier Baseline",
        "n_samples": len(y),
        "n_pos": int(n_pos),
        "n_neg": int(n_neg),
        "prevalence": float(prevalence),
        "n_substrates": int(len(np.unique(groups))),
        "feature_dim": int(X.shape[1]),
        "combinations": {},
    }

    all_results = []
    print("\n=== Cross-Validation Results ===")

    for clf_name, clf in classifiers.items():
        for cv_name, cv in cv_strategies.items():
            print(f"\n{clf_name} + {cv_name}:")
            roc_aucs, pr_aucs = run_cv(X, y, groups, clf_name, clf, cv, cv_name)

            if roc_aucs:
                combo_key = f"{clf_name}_{cv_name}"
                combo_result = {
                    "roc_auc_mean": float(np.mean(roc_aucs)),
                    "roc_auc_std": float(np.std(roc_aucs)),
                    "pr_auc_mean": float(np.mean(pr_aucs)),
                    "pr_auc_std": float(np.std(pr_aucs)),
                    "n_folds": len(roc_aucs),
                    "roc_aucs": [float(x) for x in roc_aucs],
                    "pr_aucs": [float(x) for x in pr_aucs],
                }
                results["combinations"][combo_key] = combo_result
                all_results.append({
                    "clf": clf_name, "cv": cv_name,
                    "roc_mean": combo_result["roc_auc_mean"],
                    "roc_std": combo_result["roc_auc_std"],
                    "pr_mean": combo_result["pr_auc_mean"],
                    "pr_std": combo_result["pr_auc_std"],
                })
                print(f"  => ROC-AUC: {combo_result['roc_auc_mean']:.4f} +/- {combo_result['roc_auc_std']:.4f}")
                print(f"  => PR-AUC:  {combo_result['pr_auc_mean']:.4f} +/- {combo_result['pr_auc_std']:.4f}")

    # Baselines
    results["baselines"] = {
        "random_roc_auc": 0.500,
        "random_pr_auc": float(prevalence),
    }

    # Summary table
    print("\n=== Summary ===")
    print(f"{'Classifier':<22} {'CV Method':<25} {'ROC-AUC':>15} {'PR-AUC':>15}")
    print("-" * 80)
    for r in all_results:
        print(f"{r['clf']:<22} {r['cv']:<25} {r['roc_mean']:.4f}+/-{r['roc_std']:.4f} "
              f"{r['pr_mean']:.4f}+/-{r['pr_std']:.4f}")
    print(f"{'Constant baseline':<22} {'—':<25} {'0.5000':>15} {prevalence:.4f}")

    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Bar chart of ROC-AUC
    ax = axes[0]
    labels = [f"{r['clf'][:6]}\n{r['cv'][:10]}" for r in all_results]
    roc_means = [r["roc_mean"] for r in all_results]
    roc_stds = [r["roc_std"] for r in all_results]
    colors = ["steelblue" if "Logistic" in r["clf"] else "coral" for r in all_results]
    bars = ax.bar(range(len(labels)), roc_means, yerr=roc_stds, capsize=5, color=colors, alpha=0.7)
    ax.axhline(0.5, color="gray", linestyle="--", label="Random (0.5)")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("ROC-AUC")
    ax.set_title("Substrate Fingerprint Classifier: ROC-AUC")
    ax.set_ylim(0.3, 0.7)
    ax.legend()

    # Bar chart of PR-AUC
    ax = axes[1]
    pr_means = [r["pr_mean"] for r in all_results]
    pr_stds = [r["pr_std"] for r in all_results]
    bars = ax.bar(range(len(labels)), pr_means, yerr=pr_stds, capsize=5, color=colors, alpha=0.7)
    ax.axhline(prevalence, color="gray", linestyle="--", label=f"Prevalence ({prevalence:.3f})")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("PR-AUC")
    ax.set_title("Substrate Fingerprint Classifier: PR-AUC")
    ax.set_ylim(0.0, 0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "e5_new_fingerprint_classifier.png"), dpi=150)
    plt.close()

    # Save results
    with open(os.path.join(OUTPUT_DIR, "e5_new_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
    print(f"\nResults saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
