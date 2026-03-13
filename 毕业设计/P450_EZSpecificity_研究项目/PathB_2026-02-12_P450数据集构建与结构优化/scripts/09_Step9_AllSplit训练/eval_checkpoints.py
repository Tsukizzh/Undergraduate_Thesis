"""
Evaluate all checkpoints: compute AUC/AUPR/loss + per-family metrics,
save to CSV, and generate comprehensive training report figures.

Best checkpoint selected by AUC-ROC (highest overall val AUC).

Figures generated:
  Fig 1 (3x2): Train/Val loss, AUC+AUPR curves, LR schedule, Gradient norm,
                Per-family AUC bars, Per-family AUC trend
  Fig 2 (3x2): Score distribution, PR curve, ROC curve, Threshold sweep,
                Confusion matrix (best F1 threshold), Confusion matrix (0.5 threshold)
  Fig 3:       Family performance heatmap (AUC-ROC x AUPR)

Usage:
    cd D:/EZSpecificity_Project/src
    D:/anaconda3/envs/torch/python.exe "../毕业设计/.../scripts/09_Step9_AllSplit训练/eval_checkpoints.py"
"""
import os, sys, glob, time, csv, math
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", "src"))
PATHB_DIR = os.path.normpath(os.path.join(SCRIPT_DIR, "..", ".."))
RESULT_DIR = os.path.join(PATHB_DIR, "results", "09_Step9_AllSplit训练")
CKPT_DIR = os.path.join(RESULT_DIR, "checkpoints")
CONFIG_PATH = os.path.join(SCRIPT_DIR, "train_allsplit_config.yml")
CACHE_DIR = os.path.normpath(os.path.join(SRC_DIR, "..", "data", "step9_structure_cache"))

for d in (SCRIPT_DIR, SRC_DIR):
    if d not in sys.path:
        sys.path.insert(0, d)

import warnings
warnings.filterwarnings("ignore")

import lmdb
_original_lmdb_open = lmdb.open
def _patched_lmdb_open(path, **kwargs):
    if kwargs.get("readonly", False) or not kwargs.get("create", True):
        kwargs["map_size"] = 256 * 1024 ** 3
    return _original_lmdb_open(path, **kwargs)
lmdb.open = _patched_lmdb_open

FAMILY_NAMES = {
    "0": "brenda", "1": "Duf", "2": "Esterase", "3": "Gt_acceptor",
    "4": "Nitrilase", "5": "Phosphatase", "6": "Thiolase",
}


def evaluate_checkpoint(model_cls, config, dm, ckpt_path):
    """Run validation and return (metrics_dict, logits, labels)."""
    import torch
    import pytorch_lightning as pl
    from sklearn import metrics as skmetrics

    basename = os.path.basename(ckpt_path)
    ckpt_data = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    epoch = ckpt_data.get("epoch", -1)
    print(f"\nEvaluating: {basename} (epoch {epoch})")

    model = model_cls(config)
    trainer = pl.Trainer(
        accelerator="gpu", devices=1, precision=16,
        enable_progress_bar=True, logger=False,
    )
    trainer.validate(model, datamodule=dm, ckpt_path=ckpt_path)

    m = trainer.callback_metrics
    result = {"epoch": epoch, "checkpoint": basename}

    def get(key):
        v = m.get(key)
        if v is None:
            return None
        v = v.item()
        return v if not math.isnan(v) else None

    result["val_auc"] = get("auc/val")
    result["val_aupr"] = get("aupr/val")
    result["val_loss"] = get("sp_loss/val")

    for i in range(7):
        result[f"val_auc_{i}"] = get(f"auc/{i}/val")
        result[f"val_aupr_{i}"] = get(f"aupr/{i}/val")

    # Extract raw logits/labels for detailed analysis
    logits = getattr(model, "logits", None)
    labels = getattr(model, "labels", None)

    return result, logits, labels


def _load_training_csv():
    """Load metrics_history.csv from results dir if available (has train_loss, lr, grad_norm)."""
    csv_path = os.path.join(RESULT_DIR, "metrics_history.csv")
    if not os.path.exists(csv_path):
        return {}
    import csv as csv_mod
    rows = {}
    with open(csv_path, encoding="utf-8") as f:
        for row in csv_mod.DictReader(f):
            ep = int(row.get("epoch", -1))
            rows[ep] = row
    return rows


def generate_figures(results, best_logits, best_labels, best_epoch, output_dir=RESULT_DIR):
    """Generate all report figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn import metrics as skmetrics

    results = sorted(results, key=lambda x: x["epoch"])
    epochs = [r["epoch"] for r in results]

    # Try to load training CSV for train_loss / lr / grad_norm
    csv_data = _load_training_csv()

    def csv_col(key):
        vals = []
        for ep in epochs:
            row = csv_data.get(ep)
            if row and row.get(key):
                try:
                    vals.append(float(row[key]))
                except (ValueError, TypeError):
                    vals.append(np.nan)
            else:
                vals.append(np.nan)
        return vals

    train_loss = csv_col("train_loss")
    lr_vals = csv_col("lr")
    grad_norm = csv_col("grad_norm")

    # ================================================================
    # Fig 1: Training Dynamics (3x2)
    # ================================================================
    fig, axes = plt.subplots(3, 2, figsize=(14, 15))

    # 1a: Train/Val Loss
    ax = axes[0, 0]
    val_loss = [r.get("val_loss") or np.nan for r in results]
    if not all(np.isnan(v) for v in train_loss):
        ax.plot(epochs, train_loss, "o-", label="Train", color="#FF5722", markersize=5, linewidth=2)
    ax.plot(epochs, val_loss, "s-", label="Val", color="#2196F3", markersize=5, linewidth=2)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Loss"); ax.set_title("Train / Val Loss")
    ax.legend(); ax.grid(True, alpha=0.3)

    # 1b: AUC + AUPR
    ax = axes[0, 1]
    val_auc = [r.get("val_auc") or np.nan for r in results]
    val_aupr = [r.get("val_aupr") or np.nan for r in results]
    ax.plot(epochs, val_auc, "o-", label="AUC-ROC", color="#2196F3", markersize=5, linewidth=2)
    ax.plot(epochs, val_aupr, "s-", label="AUPR", color="#4CAF50", markersize=5, linewidth=2)
    ax.set_xlabel("Epoch"); ax.set_ylabel("Score"); ax.set_title("Val AUC-ROC & AUPR")
    ax.legend(); ax.grid(True, alpha=0.3)
    bi = int(np.nanargmax(val_auc))
    ax.annotate(f"Best AUC: {val_auc[bi]:.4f} (ep{epochs[bi]})",
                xy=(epochs[bi], val_auc[bi]), xytext=(0, 15), textcoords="offset points",
                arrowprops=dict(arrowstyle="->", color="red"),
                color="red", fontweight="bold", ha="center", fontsize=9)

    # 1c: Learning Rate
    ax = axes[1, 0]
    if not all(np.isnan(v) for v in lr_vals):
        ax.plot(epochs, lr_vals, "D-", color="#FF9800", markersize=5, linewidth=2)
        ax.set_yscale("log")
        ax.set_xlabel("Epoch"); ax.set_ylabel("Learning Rate")
        ax.set_title("Learning Rate Schedule")
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, "No LR data (pre-upgrade epochs)", ha="center", va="center",
                transform=ax.transAxes, fontsize=12, color="gray")
        ax.set_title("Learning Rate Schedule")

    # 1d: Gradient Norm
    ax = axes[1, 1]
    if not all(np.isnan(v) for v in grad_norm):
        ax.plot(epochs, grad_norm, "^-", color="#9C27B0", markersize=5, linewidth=2)
        ax.set_xlabel("Epoch"); ax.set_ylabel("Avg Gradient L2 Norm")
        ax.set_title("Gradient Norm (per epoch avg)")
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, "No gradient norm data (pre-upgrade epochs)", ha="center",
                va="center", transform=ax.transAxes, fontsize=12, color="gray")
        ax.set_title("Gradient Norm")

    # 1e: Per-Family AUC bars at best epoch
    ax = axes[2, 0]
    best_r = [r for r in results if r["epoch"] == best_epoch][0]
    names, aucs_bar = [], []
    fam_colors = ["#E53935", "#1E88E5", "#43A047", "#FB8C00", "#8E24AA", "#00ACC1", "#6D4C41"]
    for i in range(7):
        v = best_r.get(f"val_auc_{i}")
        name = FAMILY_NAMES.get(str(i), str(i))
        if v is not None:
            names.append(name)
            aucs_bar.append(v)
    if names:
        bars = ax.bar(range(len(names)), aucs_bar, color=fam_colors[:len(names)], alpha=0.85)
        ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=30, ha="right")
        for bar, v in zip(bars, aucs_bar):
            ax.text(bar.get_x() + bar.get_width()/2, v + 0.01, f"{v:.3f}",
                    ha="center", va="bottom", fontsize=8)
        overall = best_r.get("val_auc")
        if overall:
            ax.axhline(y=overall, color="black", linestyle="--", alpha=0.5,
                        label=f"Overall={overall:.3f}")
            ax.legend(fontsize=9)
    ax.set_ylabel("AUC-ROC"); ax.set_title(f"Per-Family AUC (Best ep{best_epoch})")
    ax.grid(True, alpha=0.3, axis="y")

    # 1f: Per-Family AUC trend
    ax = axes[2, 1]
    for i in range(7):
        name = FAMILY_NAMES.get(str(i), str(i))
        vals = [r.get(f"val_auc_{i}") or np.nan for r in results]
        if not all(np.isnan(v) for v in vals):
            ax.plot(epochs, vals, "-", label=name, color=fam_colors[i], alpha=0.7, linewidth=1.5)
    ax.set_xlabel("Epoch"); ax.set_ylabel("AUC-ROC"); ax.set_title("Per-Family AUC Trend")
    ax.legend(fontsize=7, ncol=2); ax.grid(True, alpha=0.3)

    plt.suptitle("Step 9: AllSplit Training Report", fontsize=16, fontweight="bold")
    plt.tight_layout()
    p = os.path.join(output_dir, "fig1_training_dynamics.png")
    plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[Fig1] Saved: {p}")

    # ================================================================
    # Fig 2: Detailed Analysis of Best Checkpoint (requires logits)
    # ================================================================
    if best_logits is not None and best_labels is not None:
        logits = np.array(best_logits).ravel()
        labels = np.array(best_labels).ravel()
        probs = 1.0 / (1.0 + np.exp(-logits))

        fig, axes = plt.subplots(3, 2, figsize=(14, 15))

        # 2a: Score distribution
        ax = axes[0, 0]
        pos_scores = probs[labels == 1]
        neg_scores = probs[labels == 0]
        bins = np.linspace(0, 1, 50)
        ax.hist(neg_scores, bins=bins, alpha=0.6, label=f"Neg (n={len(neg_scores)})",
                color="#2196F3", density=True)
        ax.hist(pos_scores, bins=bins, alpha=0.6, label=f"Pos (n={len(pos_scores)})",
                color="#E53935", density=True)
        ax.set_xlabel("Predicted Probability"); ax.set_ylabel("Density")
        ax.set_title("Score Distribution"); ax.legend(); ax.grid(True, alpha=0.3)

        # 2b: PR Curve
        ax = axes[0, 1]
        precision, recall, _ = skmetrics.precision_recall_curve(labels, probs)
        aupr_val = skmetrics.average_precision_score(labels, probs)
        prevalence = labels.mean()
        ax.plot(recall, precision, color="#4CAF50", linewidth=2, label=f"AUPR={aupr_val:.4f}")
        ax.axhline(y=prevalence, color="gray", linestyle="--", alpha=0.7,
                    label=f"Prevalence={prevalence:.3f}")
        ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
        ax.set_title("Precision-Recall Curve"); ax.legend(); ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)

        # 2c: ROC Curve
        ax = axes[1, 0]
        fpr, tpr, _ = skmetrics.roc_curve(labels, probs)
        auroc = skmetrics.auc(fpr, tpr)
        ax.plot(fpr, tpr, color="#2196F3", linewidth=2, label=f"AUC={auroc:.4f}")
        ax.plot([0, 1], [0, 1], "k--", alpha=0.3)
        ax.set_xlabel("FPR"); ax.set_ylabel("TPR")
        ax.set_title("ROC Curve"); ax.legend(); ax.grid(True, alpha=0.3)

        # 2d: Threshold sweep (F1, MCC, Precision, Recall)
        ax = axes[1, 1]
        thresholds = np.linspace(0.01, 0.99, 200)
        f1s, mccs, precs, recs = [], [], [], []
        for t in thresholds:
            preds = (probs >= t).astype(int)
            tp = ((preds == 1) & (labels == 1)).sum()
            fp = ((preds == 1) & (labels == 0)).sum()
            fn = ((preds == 0) & (labels == 1)).sum()
            tn = ((preds == 0) & (labels == 0)).sum()
            p = tp / (tp + fp) if (tp + fp) > 0 else 0
            r = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * p * r / (p + r) if (p + r) > 0 else 0
            denom = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) if (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) > 0 else 1
            mcc = (tp*tn - fp*fn) / denom
            f1s.append(f1); mccs.append(mcc); precs.append(p); recs.append(r)
        ax.plot(thresholds, precs, label="Precision", color="#FF9800", linewidth=1.5)
        ax.plot(thresholds, recs, label="Recall", color="#03A9F4", linewidth=1.5)
        ax.plot(thresholds, f1s, label="F1", color="#4CAF50", linewidth=2)
        ax.plot(thresholds, mccs, label="MCC", color="#9C27B0", linewidth=2)
        best_f1_idx = np.argmax(f1s)
        ax.axvline(x=thresholds[best_f1_idx], color="red", linestyle="--", alpha=0.5,
                    label=f"Best F1={f1s[best_f1_idx]:.3f} @t={thresholds[best_f1_idx]:.2f}")
        ax.set_xlabel("Threshold"); ax.set_ylabel("Score")
        ax.set_title("Threshold Sweep"); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        # 2e: Confusion Matrix (at best F1 threshold)
        ax = axes[2, 0]
        best_t = thresholds[best_f1_idx]
        preds = (probs >= best_t).astype(int)
        cm = skmetrics.confusion_matrix(labels, preds)
        im = ax.imshow(cm, interpolation="nearest", cmap="Blues")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Neg", "Pos"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["Neg", "Pos"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("True")
        ax.set_title(f"Confusion Matrix (threshold={best_t:.2f})")
        for i_cm in range(2):
            for j_cm in range(2):
                color = "white" if cm[i_cm, j_cm] > cm.max() / 2 else "black"
                ax.text(j_cm, i_cm, f"{cm[i_cm, j_cm]}", ha="center", va="center",
                        fontsize=16, fontweight="bold", color=color)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # 2f: Confusion matrix at threshold=0.5 for comparison
        ax = axes[2, 1]
        preds_05 = (probs >= 0.5).astype(int)
        cm05 = skmetrics.confusion_matrix(labels, preds_05)
        im2 = ax.imshow(cm05, interpolation="nearest", cmap="Oranges")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Neg", "Pos"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["Neg", "Pos"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("True")
        ax.set_title("Confusion Matrix (threshold=0.50)")
        for i_cm in range(2):
            for j_cm in range(2):
                color = "white" if cm05[i_cm, j_cm] > cm05.max() / 2 else "black"
                ax.text(j_cm, i_cm, f"{cm05[i_cm, j_cm]}", ha="center", va="center",
                        fontsize=16, fontweight="bold", color=color)
        plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.04)

        plt.suptitle(f"Best Checkpoint Analysis (Epoch {best_epoch})", fontsize=16, fontweight="bold")
        plt.tight_layout()
        p = os.path.join(output_dir, "fig2_best_checkpoint_analysis.png")
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"[Fig2] Saved: {p}")

        # ================================================================
        # Fig 3: Family Performance Heatmap
        # ================================================================
        fig, ax = plt.subplots(figsize=(10, 5))
        fam_data = []
        fam_names_valid = []
        for i in range(7):
            name = FAMILY_NAMES.get(str(i), str(i))
            auc_v = best_r.get(f"val_auc_{i}")
            aupr_v = best_r.get(f"val_aupr_{i}")
            if auc_v is not None and aupr_v is not None:
                fam_names_valid.append(name)
                fam_data.append([auc_v, aupr_v])
        if fam_data:
            data = np.array(fam_data)
            # Add overall
            oa = best_r.get("val_auc")
            op = best_r.get("val_aupr")
            if oa and op:
                fam_names_valid.append("OVERALL")
                data = np.vstack([data, [oa, op]])

            im = ax.imshow(data, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["AUC-ROC", "AUPR"])
            ax.set_yticks(range(len(fam_names_valid)))
            ax.set_yticklabels(fam_names_valid)
            for i in range(data.shape[0]):
                for j in range(data.shape[1]):
                    color = "white" if data[i, j] < 0.5 else "black"
                    ax.text(j, i, f"{data[i,j]:.3f}", ha="center", va="center",
                            fontsize=11, fontweight="bold", color=color)
            plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
            ax.set_title(f"Family Performance Heatmap (Epoch {best_epoch})", fontsize=14)

        plt.tight_layout()
        p = os.path.join(output_dir, "fig3_family_heatmap.png")
        plt.savefig(p, dpi=150, bbox_inches="tight"); plt.close()
        print(f"[Fig3] Saved: {p}")


def main():
    import torch
    import pytorch_lightning as pl
    import yaml
    from easydict import EasyDict
    from Models.ss import SS
    from main_training_cached import CachedTrainingDataModule, resolve_cache_paths

    # Edge-mode fix
    from Datasets.Structure.transforms import EdgeConnection
    _orig_call = EdgeConnection.__call__
    def _fixed_call(self, data):
        result = _orig_call(self, data)
        if hasattr(result, 'edge_attr') and hasattr(result, 'edge_index'):
            if result.edge_attr.shape[0] != result.edge_index.shape[1]:
                n = min(result.edge_attr.shape[0], result.edge_index.shape[1])
                result.edge_attr = result.edge_attr[:n]
                result.edge_index = result.edge_index[:, :n]
        return result
    EdgeConnection.__call__ = _fixed_call

    with open(CONFIG_PATH, encoding="utf-8") as f:
        config = EasyDict(yaml.safe_load(f))
    cache_paths = resolve_cache_paths(config, CACHE_DIR)

    # Locate checkpoints
    ckpt_files = sorted(glob.glob(os.path.join(CKPT_DIR, "allsplit-fold0-*.ckpt")))
    last_ckpt = os.path.join(CKPT_DIR, "last.ckpt")
    if os.path.exists(last_ckpt):
        last_epoch = torch.load(last_ckpt, map_location="cpu", weights_only=False).get("epoch", -1)
        named_epochs = []
        for f in ckpt_files:
            e = torch.load(f, map_location="cpu", weights_only=False).get("epoch", -2)
            named_epochs.append(e)
        if last_epoch not in named_epochs:
            ckpt_files.append(last_ckpt)

    if not ckpt_files:
        print("No checkpoints found!"); return

    dm = CachedTrainingDataModule(config, cache_paths, edge_mode="fixed", num_workers=2,
                                   use_prefetch_wrapper=False)

    # Evaluate all checkpoints
    results = []
    best_logits, best_labels, best_epoch = None, None, -1
    best_auc = -1

    for ckpt_path in ckpt_files:
        r, logits, labels = evaluate_checkpoint(SS, config, dm, ckpt_path)
        results.append(r)
        auc = r.get("val_auc") or 0
        if auc > best_auc:
            best_auc = auc
            best_logits = logits
            best_labels = labels
            best_epoch = r["epoch"]

    # Print summary table
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    header = f"{'Epoch':<6} {'AUC':<8} {'AUPR':<8} {'Loss':<8}"
    for i in range(7):
        header += f" {FAMILY_NAMES[str(i)][:6]:<7}"
    header += "  (family=AUC)"
    print(header)
    print("-"*80)
    for r in sorted(results, key=lambda x: x["epoch"]):
        line = f"{r['epoch']:<6} {r.get('val_auc',0) or 0:.4f}  {r.get('val_aupr',0) or 0:.4f}  {r.get('val_loss',0) or 0:.4f}"
        for i in range(7):
            v = r.get(f"val_auc_{i}")
            line += f"  {v:.3f}" if v else "    N/A"
        print(line)

    # Determine max epoch for output folder name
    max_epoch = max(r["epoch"] for r in results)
    output_dir = os.path.join(RESULT_DIR, f"eval_epoch{max_epoch:02d}")
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput dir: {output_dir}")

    # Save CSV
    csv_path = os.path.join(output_dir, "metrics_history.csv")
    fieldnames = ["epoch", "timestamp", "lr", "grad_norm", "train_loss", "val_loss", "val_auc", "val_aupr"]
    for i in range(7):
        fieldnames += [f"val_auc_{i}", f"val_aupr_{i}"]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for r in sorted(results, key=lambda x: x["epoch"]):
            r["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
            writer.writerow(r)
    print(f"CSV saved: {csv_path}")

    # Generate figures
    generate_figures(results, best_logits, best_labels, best_epoch, output_dir=output_dir)
    print("\nDone!")


if __name__ == "__main__":
    main()
