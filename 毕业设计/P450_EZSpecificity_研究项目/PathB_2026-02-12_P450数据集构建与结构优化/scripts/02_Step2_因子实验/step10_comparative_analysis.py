"""
PathB Step 2 - Comparative analysis across 2x2 factorial experiments.

Compares EXP01-EXP04 predictions and produces Gate A recommendation.

Usage:
    python step10_comparative_analysis.py \
        --experiments_dir <parent of EXP01-EXP04> \
        --shared_datasets <path_to_shared/datasets> \
        --output_dir      <analysis output dir>

Outputs:
    comparative_metrics.csv
    gate_a_recommendation.txt
    figures/comparative_roc.png
    figures/factorial_heatmap.png
"""

from __future__ import annotations

import argparse
import re
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score, average_precision_score, roc_auc_score, roc_curve,
)

BASELINE_AUC_ROC = 0.6636
DELTA_THRESHOLD = 0.005
KEY_COLS = ["Dock Index", "Enzyme Index", "Substrate Index"]


def _resolve_dataset_csv(shared_datasets: Path) -> Path:
    for c in [shared_datasets / "B6_v1" / "data.csv", shared_datasets / "data.csv"]:
        if c.is_file():
            return c
    raise FileNotFoundError(f"data.csv not found under {shared_datasets}")


def _discover_experiments(parent: Path) -> List[Path]:
    dirs = [d for d in parent.iterdir()
            if d.is_dir() and re.match(r"^EXP\d+", d.name, re.IGNORECASE)]
    return sorted(dirs, key=lambda p: int(re.match(r"EXP(\d+)", p.name, re.IGNORECASE).group(1)))


def _parse_factors(name: str) -> Tuple[bool, float]:
    low = name.lower()
    include_heme = "heme" in low and "noheme" not in low
    m = re.search(r"(\d+(?:\.\d+)?)a", low)
    if not m:
        raise ValueError(f"Cannot parse radius from {name}")
    return include_heme, float(m.group(1))


def _sigmoid(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    e = np.exp(x[~pos])
    out[~pos] = e / (1.0 + e)
    return out


def _get_probs(df: pd.DataFrame) -> np.ndarray:
    if "prob" in df.columns:
        return df["prob"].astype(float).values
    if "logit" in df.columns:
        return _sigmoid(df["logit"].astype(float).values)
    if "score" in df.columns:
        s = df["score"].astype(float).values
        return s if (np.nanmin(s) >= 0 and np.nanmax(s) <= 1) else _sigmoid(s)
    raise ValueError("No prob/logit/score column in predictions.csv")


def _evaluate(exp_dir: Path, labels_df: pd.DataFrame):
    name = exp_dir.name
    heme, radius = _parse_factors(name)
    empty = {"experiment_name": name, "include_heme": heme, "pocket_radius": radius,
             "auc_roc": np.nan, "auc_pr": np.nan, "accuracy": np.nan, "n_samples": 0}

    pred_path = exp_dir / "predictions.csv"
    if not pred_path.exists():
        return empty, None, "predictions.csv not found"

    pred_df = pd.read_csv(pred_path)
    if pred_df.empty:
        return empty, None, "predictions.csv is empty"

    if all(c in pred_df.columns for c in KEY_COLS):
        merged = pred_df.drop(columns=["Label"], errors="ignore").merge(labels_df, on=KEY_COLS, how="left")
        if "Label" in pred_df.columns:
            merged["Label"] = merged["Label"].fillna(pred_df["Label"])
        eval_df = merged.dropna(subset=["Label"]).copy()
    elif "Label" in pred_df.columns:
        eval_df = pred_df.copy()
    else:
        return empty, None, "no label info"

    eval_df["Label"] = eval_df["Label"].astype(int)
    y_true = eval_df["Label"].values
    y_prob = _get_probs(eval_df)

    mask = np.isfinite(y_true.astype(float)) & np.isfinite(y_prob)
    y_true, y_prob = y_true[mask], y_prob[mask]
    empty["n_samples"] = int(len(y_true))

    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return empty, None, "insufficient classes"

    empty["auc_roc"] = float(roc_auc_score(y_true, y_prob))
    empty["auc_pr"] = float(average_precision_score(y_true, y_prob))
    empty["accuracy"] = float(accuracy_score(y_true, (y_prob >= 0.5).astype(int)))

    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc = {"name": name, "fpr": fpr, "tpr": tpr, "auc": empty["auc_roc"]}
    return empty, roc, None


def _plot_roc(roc_list, path: Path):
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for i, item in enumerate(roc_list):
        c = colors[i % len(colors)]
        ax.plot(item["fpr"], item["tpr"], lw=2, color=c,
                label=f"{item['name']} (AUC={item['auc']:.4f})")
    ax.plot([0, 1], [0, 1], "--", color="gray", lw=1)
    ax.set(xlim=(0, 1), ylim=(0, 1.02), xlabel="FPR", ylabel="TPR",
           title="PathB Step 2: Comparative ROC")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _plot_heatmap(df: pd.DataFrame, path: Path):
    matrix = np.full((2, 2), np.nan)  # rows: 10A, 6A; cols: noHeme, Heme
    for _, r in df.iterrows():
        if pd.isna(r["auc_roc"]):
            continue
        ri = 0 if abs(float(r["pocket_radius"]) - 10.0) < 1e-6 else 1
        ci = 1 if r["include_heme"] else 0
        matrix[ri, ci] = float(r["auc_roc"])

    fig, ax = plt.subplots(figsize=(6, 5))
    cmap = plt.cm.viridis.copy()
    cmap.set_bad("#e0e0e0")

    valid = matrix[np.isfinite(matrix)]
    vmin = float(valid.min()) if valid.size else 0.0
    vmax = float(valid.max()) if valid.size else 1.0
    if np.isclose(vmin, vmax):
        vmin, vmax = max(0, vmin - 0.01), min(1, vmax + 0.01)

    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["noHeme", "Heme"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["10\u00c5", "6\u00c5"])
    ax.set(xlabel="Heme Factor", ylabel="Pocket Radius",
           title="2\u00d72 Factorial Heatmap (AUC-ROC)")

    for i in range(2):
        for j in range(2):
            v = matrix[i, j]
            txt = "NA" if not np.isfinite(v) else f"{v:.4f}"
            color = "white" if np.isfinite(v) and v > (vmin + vmax) / 2 else "black"
            ax.text(j, i, txt, ha="center", va="center", color=color, fontsize=12, fontweight="bold")

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="AUC-ROC")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300)
    plt.close(fig)


def _gate_a_report(df: pd.DataFrame, notes: Dict[str, str]) -> str:
    lines = [
        "Gate A Recommendation Report",
        "=" * 70,
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"PathA baseline AUC-ROC: {BASELINE_AUC_ROC:.4f}",
        f"Decision threshold: delta > {DELTA_THRESHOLD:.4f}",
        "",
    ]

    valid = df.dropna(subset=["auc_roc"])
    if valid.empty:
        lines.append("Decision: INSUFFICIENT DATA")
        lines.append("Run Step 9 inference for all experiments first.")
    else:
        best = valid.loc[valid["auc_roc"].idxmax()]
        delta = float(best["auc_roc"]) - BASELINE_AUC_ROC
        lines += [
            f"Best experiment: {best['experiment_name']}",
            f"  include_heme={bool(best['include_heme'])}, pocket_radius={float(best['pocket_radius']):.1f}A",
            f"  AUC-ROC={float(best['auc_roc']):.4f}",
            f"  Delta vs baseline: {delta:+.4f} ({delta*100:+.2f}%)",
            "",
        ]
        if delta > DELTA_THRESHOLD:
            lines.append("Gate A Decision: PASS")
            lines.append("Recommendation: adopt this config for downstream PathB steps.")
        else:
            lines.append("Gate A Decision: HOLD")
            lines.append("Improvement does not exceed 0.5%. Keep baseline config.")

    # Main effects
    def _mean(mask):
        v = df.loc[mask, "auc_roc"].dropna()
        return float(v.mean()) if len(v) else np.nan

    heme_on = _mean(df["include_heme"] == True)
    heme_off = _mean(df["include_heme"] == False)
    r10 = _mean(np.isclose(df["pocket_radius"], 10.0))
    r6 = _mean(np.isclose(df["pocket_radius"], 6.0))

    heme_eff = heme_on - heme_off if np.isfinite(heme_on) and np.isfinite(heme_off) else np.nan
    rad_eff = r10 - r6 if np.isfinite(r10) and np.isfinite(r6) else np.nan

    fmt = lambda v: "NA" if pd.isna(v) else f"{v:+.4f}"
    lines += [
        "",
        "Main Effects:",
        f"  Heme effect  = mean(Heme) - mean(noHeme)  = {fmt(heme_eff)}",
        f"  Radius effect = mean(10A) - mean(6A)       = {fmt(rad_eff)}",
        "",
        "Per-experiment:",
    ]
    for _, r in df.iterrows():
        f4 = lambda v: "NA" if pd.isna(v) else f"{v:.4f}"
        lines.append(
            f"  {r['experiment_name']}: heme={bool(r['include_heme'])}, "
            f"radius={float(r['pocket_radius']):.0f}A, "
            f"AUC-ROC={f4(r['auc_roc'])}, AUC-PR={f4(r['auc_pr'])}, "
            f"acc={f4(r['accuracy'])}, n={int(r['n_samples'])}"
        )

    if notes:
        lines += ["", "Notes:"]
        for k, v in notes.items():
            lines.append(f"  {k}: {v}")

    lines.append("\n" + "=" * 70)
    return "\n".join(lines)


def main() -> int:
    p = argparse.ArgumentParser(description="PathB Step 2 comparative analysis.")
    p.add_argument("--experiments_dir", required=True)
    p.add_argument("--shared_datasets", required=True)
    p.add_argument("--output_dir", required=True)
    args = p.parse_args()

    exp_dir = Path(args.experiments_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    labels_df = pd.read_csv(_resolve_dataset_csv(Path(args.shared_datasets).resolve()),
                            usecols=KEY_COLS + ["Label"])

    exp_dirs = _discover_experiments(exp_dir)
    if not exp_dirs:
        raise FileNotFoundError(f"No EXP* dirs under {exp_dir}")

    rows, rocs, notes = [], [], {}
    for d in exp_dirs:
        row, roc, note = _evaluate(d, labels_df)
        rows.append(row)
        if roc:
            rocs.append(roc)
        if note:
            notes[d.name] = note
            print(f"[WARN] {d.name}: {note}")

    metrics_df = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)

    metrics_df.to_csv(out_dir / "comparative_metrics.csv", index=False)
    _plot_roc(rocs, out_dir / "figures" / "comparative_roc.png")
    _plot_heatmap(metrics_df, out_dir / "figures" / "factorial_heatmap.png")
    (out_dir / "gate_a_recommendation.txt").write_text(
        _gate_a_report(metrics_df, notes), encoding="utf-8")

    print(f"[OK] Metrics: {out_dir / 'comparative_metrics.csv'}")
    print(f"[OK] ROC figure: {out_dir / 'figures' / 'comparative_roc.png'}")
    print(f"[OK] Heatmap: {out_dir / 'figures' / 'factorial_heatmap.png'}")
    print(f"[OK] Gate A report: {out_dir / 'gate_a_recommendation.txt'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}")
        import traceback
        traceback.print_exc()
        raise SystemExit(1)
