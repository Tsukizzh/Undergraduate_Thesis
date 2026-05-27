#!/usr/bin/env python3
"""Plot PathD Q2 sequence-similarity summary figures.

Inputs are compact CSV files exported from the server split CSVs. The script
does not read or modify the original experiment data.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "results" / "q02_sequence_similarity_figures" / "figure_data_20260525_v1"
OUT_DIR = ROOT / "results" / "q02_sequence_similarity_figures" / "figures_20260525_v1"


COLORS = {
    "train": "#4C78A8",
    "val": "#F58518",
    "test": "#54A24B",
    "auc": "#4C78A8",
    "aupr": "#E45756",
    "accent": "#7A5195",
    "light": "#DDE7F2",
}


def setup_fonts() -> None:
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT_DIR / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


def annotate_bars(ax: plt.Axes, bars, fmt: str = "{:.0f}") -> None:
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            fmt.format(height),
            ha="center",
            va="bottom",
            fontsize=8,
        )


def plot_threshold_structure(overall: pd.DataFrame) -> None:
    subset = overall[overall["exp"].isin(["EXP006", "EXP005", "EXP007"])].copy()
    subset["threshold"] = pd.Categorical(
        subset["exp"].map({"EXP006": "NN40", "EXP005": "NN60", "EXP007": "NN80"}),
        categories=["NN40", "NN60", "NN80"],
        ordered=True,
    )
    subset = subset.sort_values("threshold")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    bars = axes[0].bar(subset["threshold"], subset["total_groups"], color=["#7A5195", "#4C78A8", "#54A24B"])
    axes[0].set_title("阈值越低，组件合并越明显")
    axes[0].set_ylabel("组件总数")
    annotate_bars(axes[0], bars)

    bars = axes[1].bar(
        subset["threshold"],
        subset["largest_group_enzymes"],
        color=["#7A5195", "#4C78A8", "#54A24B"],
    )
    axes[1].set_title("最大组件包含的酶数")
    axes[1].set_ylabel("最大组件酶数")
    annotate_bars(axes[1], bars)

    fig.suptitle("Q2 序列相似度阈值与组件结构", fontsize=15, fontweight="bold")
    save(fig, "fig01_threshold_component_structure.png")


def plot_split_distribution(split_df: pd.DataFrame) -> None:
    order = ["EXP004", "EXP006", "EXP005", "EXP010", "EXP011", "EXP007"]
    labels = {
        "EXP004": "EXP004\nid60聚类",
        "EXP006": "EXP006\nNN40",
        "EXP005": "EXP005\nNN60 2:1:1",
        "EXP010": "EXP010\nNN60 7:1.5:1.5",
        "EXP011": "EXP011\nNN60 8:1:1",
        "EXP007": "EXP007\nNN80",
    }
    fig, axes = plt.subplots(2, 1, figsize=(13.5, 8), sharex=True)
    x = np.arange(len(order))
    width = 0.24

    for offset, split in zip([-width, 0, width], ["train", "val", "test"]):
        rows = split_df[split_df["split"] == split].set_index("exp").loc[order]
        bars = axes[0].bar(
            x + offset,
            rows["groups"],
            width=width,
            label=split,
            color=COLORS[split],
        )
        annotate_bars(axes[0], bars)
        bars = axes[1].bar(
            x + offset,
            rows["enzymes"],
            width=width,
            label=split,
            color=COLORS[split],
        )
        annotate_bars(axes[1], bars)

    axes[0].set_title("每个实验的组数分配")
    axes[0].set_ylabel("组数")
    axes[1].set_title("每个实验的酶数分配")
    axes[1].set_ylabel("酶数")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([labels[x] for x in order])
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    axes[1].legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    fig.suptitle("Q2 train / val / test 分到了多少组和酶", fontsize=15, fontweight="bold")
    save(fig, "fig02_split_group_enzyme_distribution.png")


def plot_test_concentration(test_df: pd.DataFrame) -> None:
    subset = test_df[test_df["exp"].isin(["EXP005", "EXP010", "EXP011"])].copy()
    order = ["EXP005", "EXP010", "EXP011"]
    subset = subset.set_index("exp").loc[order].reset_index()
    x = np.arange(len(subset))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    series = [
        ("top5_group_sample_share", "Top5组件样本占比", "#4C78A8"),
        ("top5_group_positive_share", "Top5组件正样本占比", "#F58518"),
        ("top_positive_enzyme_share", "最大正样本酶占比", "#E45756"),
    ]
    for i, (col, label, color) in enumerate(series):
        values = subset[col].to_numpy() * 100
        bars = ax.bar(x + (i - 1) * width, values, width=width, label=label, color=color)
        annotate_bars(ax, bars, "{:.1f}%")

    ax.set_xticks(x)
    ax.set_xticklabels(["EXP005\n2:1:1", "EXP010\n7:1.5:1.5", "EXP011\n8:1:1"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("占 test 集比例")
    ax.set_title("strict NN60 下 test 集集中度")
    ax.legend(ncol=1, loc="upper left")
    ax.grid(axis="y", alpha=0.2)
    save(fig, "fig03_test_concentration_strict_nn60.png")


def plot_metric_comparison(metrics: pd.DataFrame) -> None:
    order = ["EXP008", "EXP004", "EXP007", "EXP010", "EXP005", "EXP011", "EXP006"]
    labels = {
        "EXP008": "EXP008\n随机",
        "EXP004": "EXP004\nid60聚类",
        "EXP007": "EXP007\nNN80",
        "EXP010": "EXP010\nNN60 7:1.5:1.5",
        "EXP005": "EXP005\nNN60 2:1:1",
        "EXP011": "EXP011\nNN60 8:1:1",
        "EXP006": "EXP006\nNN40",
    }
    df = metrics.set_index("exp").loc[order].reset_index()
    x = np.arange(len(df))
    width = 0.36

    fig, ax = plt.subplots(figsize=(13, 5.2))
    auc_bars = ax.bar(x - width / 2, df["test_auc"], width=width, label="Test AUC", color=COLORS["auc"])
    aupr_bars = ax.bar(x + width / 2, df["test_aupr"], width=width, label="Test AUPR", color=COLORS["aupr"])
    annotate_bars(ax, auc_bars, "{:.3f}")
    annotate_bars(ax, aupr_bars, "{:.3f}")
    ax.set_xticks(x)
    ax.set_xticklabels([labels[e] for e in df["exp"]])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("分数")
    ax.set_title("Q2 不同划分方式下的测试表现")
    ax.legend()
    ax.grid(axis="y", alpha=0.2)
    save(fig, "fig04_test_auc_aupr_comparison.png")


def main() -> None:
    setup_fonts()
    overall = pd.read_csv(DATA_DIR / "overall_group_structure.csv")
    split_df = pd.read_csv(DATA_DIR / "split_group_distribution.csv")
    test_df = pd.read_csv(DATA_DIR / "test_concentration.csv")
    metrics = pd.read_csv(DATA_DIR / "test_metrics.csv")

    plot_threshold_structure(overall)
    plot_split_distribution(split_df)
    plot_test_concentration(test_df)
    plot_metric_comparison(metrics)
    print(f"Wrote figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
