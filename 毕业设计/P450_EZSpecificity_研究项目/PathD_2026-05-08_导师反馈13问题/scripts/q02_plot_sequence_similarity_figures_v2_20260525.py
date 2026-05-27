#!/usr/bin/env python3
"""Draw cleaner Q2 sequence-similarity figures for presentation use."""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "results" / "q02_sequence_similarity_figures" / "figure_data_20260525_v1"
OUT_DIR = ROOT / "results" / "q02_sequence_similarity_figures" / "figures_20260525_v2"


INK = "#202124"
MUTED = "#667085"
GRID = "#E6E8EB"
BLUE = "#2F5D8C"
BLUE_LIGHT = "#D8E3EF"
RED = "#B94A48"
GRAY = "#A7ADB5"
GRAY_LIGHT = "#F3F5F7"
PURPLE = "#6E5A8A"


def setup() -> None:
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["axes.edgecolor"] = "#D0D5DD"
    plt.rcParams["axes.labelcolor"] = INK
    plt.rcParams["xtick.color"] = INK
    plt.rcParams["ytick.color"] = INK


def save(fig: plt.Figure, name: str) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / name, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def node(ax: plt.Axes, xy: tuple[float, float], label: str, color: str = BLUE) -> None:
    circ = patches.Circle(xy, radius=0.055, facecolor="white", edgecolor=color, linewidth=2)
    ax.add_patch(circ)
    ax.text(*xy, label, ha="center", va="center", fontsize=10, color=INK, fontweight="bold")


def rounded_box(
    ax: plt.Axes,
    xy: tuple[float, float],
    width: float,
    height: float,
    label: str,
    edge: str = BLUE,
    face: str = "white",
    linestyle: str = "-",
) -> None:
    box = patches.FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.018",
        facecolor=face,
        edgecolor=edge,
        linewidth=1.8,
        linestyle=linestyle,
    )
    ax.add_patch(box)
    ax.text(xy[0] + width / 2, xy[1] + height - 0.04, label, ha="center", va="top", fontsize=10, color=INK)


def edge(
    ax: plt.Axes,
    a: tuple[float, float],
    b: tuple[float, float],
    color: str = RED,
    linestyle: str = "-",
    rad: float = 0.0,
) -> None:
    arrow = patches.FancyArrowPatch(
        a,
        b,
        connectionstyle=f"arc3,rad={rad}",
        arrowstyle="-",
        linewidth=2.0,
        linestyle=linestyle,
        color=color,
        alpha=0.9,
    )
    ax.add_patch(arrow)


def plot_concept_network() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    ax = axes[0]
    ax.set_title("EXP004：按 id60 cluster 分", fontsize=16, fontweight="bold", color=INK, pad=10)
    rounded_box(ax, (0.08, 0.45), 0.32, 0.32, "cluster c0461", edge=BLUE, face=GRAY_LIGHT, linestyle="--")
    rounded_box(ax, (0.60, 0.45), 0.32, 0.32, "cluster c0675", edge=BLUE, face=GRAY_LIGHT, linestyle="--")
    node(ax, (0.19, 0.58), "1191")
    node(ax, (0.30, 0.58), "1190")
    node(ax, (0.76, 0.58), "1520")
    edge(ax, (0.24, 0.58), (0.71, 0.58), color=RED, linestyle="--", rad=0.05)
    ax.text(0.50, 0.65, "跨 cluster 仍达到 NN60\n相似度 78.0%，覆盖度 100%", ha="center", va="center", fontsize=10, color=RED)
    ax.text(0.24, 0.39, "EXP004: test", ha="center", va="center", fontsize=10, color=MUTED)
    ax.text(0.76, 0.39, "EXP004: train", ha="center", va="center", fontsize=10, color=MUTED)
    rounded_box(ax, (0.12, 0.10), 0.76, 0.18, "只保证同一 cluster 不拆开；跨 cluster 的近邻可能落在 train/test 两边", edge=RED, face="#FFF5F5")

    ax = axes[1]
    ax.set_title("EXP005：按所有 NN60 连接分", fontsize=16, fontweight="bold", color=INK, pad=10)
    rounded_box(ax, (0.13, 0.38), 0.74, 0.42, "strict NN60 组件", edge=PURPLE, face="#F6F2FA", linestyle="--")
    node(ax, (0.33, 0.58), "1191", color=PURPLE)
    node(ax, (0.45, 0.58), "1190", color=PURPLE)
    node(ax, (0.66, 0.58), "1520", color=PURPLE)
    edge(ax, (0.38, 0.58), (0.61, 0.58), color=RED)
    ax.text(0.50, 0.69, "只要存在 >=60% 且覆盖度 >=80% 的连接\n就合并到同一组件", ha="center", va="center", fontsize=10, color=INK)
    ax.text(0.50, 0.32, "同一组件整体进入一个 split；train-test NN60 连接数 = 0", ha="center", va="center", fontsize=11, color=PURPLE, fontweight="bold")
    rounded_box(ax, (0.24, 0.12), 0.52, 0.12, "结果：这组不能再拆成 train 和 test", edge=PURPLE, face="white")

    fig.suptitle("为什么 id60 cluster 和 strict NN60 组件不同", fontsize=18, fontweight="bold", color=INK, y=0.98)
    save(fig, "fig01_cluster_vs_strict_nn60_network.png")


def plot_component_size_spectrum(component_sizes: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(12.5, 6.2))
    configs = [
        ("EXP006", "NN40 最严格", PURPLE),
        ("EXP005", "NN60 主实验", BLUE),
        ("EXP007", "NN80 较宽松", "#4A6F4A"),
        ("EXP004", "id60 cluster", GRAY),
    ]
    for exp, label, color in configs:
        vals = component_sizes.loc[component_sizes["exp"] == exp, "n_enzymes"].sort_values(ascending=False).to_numpy()
        x = np.arange(1, len(vals) + 1)
        ax.plot(x, vals, color=color, linewidth=2.2, label=f"{label}  最大={int(vals.max())}")
        ax.scatter(x[:25], vals[:25], s=12, color=color, alpha=0.55)
    ax.set_yscale("log")
    ax.set_xlabel("组排名，按包含酶数从大到小")
    ax.set_ylabel("每组包含的酶数，对数刻度")
    ax.set_title("组件大小谱：阈值越低，越容易形成大组件", fontsize=16, fontweight="bold", color=INK)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.legend(frameon=False, loc="upper right")
    ax.text(
        0.02,
        0.05,
        "读取服务器 split CSV 后按组重新统计。NN40 最大组件包含 433 个酶，说明远缘阈值会把大量酶连成大块。",
        transform=ax.transAxes,
        fontsize=10,
        color=MUTED,
    )
    save(fig, "fig02_component_size_spectrum.png")


def plot_split_bubble_matrix(split_df: pd.DataFrame) -> None:
    order = ["EXP004", "EXP006", "EXP005", "EXP010", "EXP011", "EXP007"]
    labels = {
        "EXP004": "EXP004  id60聚类",
        "EXP006": "EXP006  NN40",
        "EXP005": "EXP005  NN60 2:1:1",
        "EXP010": "EXP010  NN60 7:1.5:1.5",
        "EXP011": "EXP011  NN60 8:1:1",
        "EXP007": "EXP007  NN80",
    }
    splits = ["train", "val", "test"]
    fig, ax = plt.subplots(figsize=(13.5, 6.5))
    ax.set_xlim(-0.45, len(splits) - 0.05)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.invert_yaxis()
    ax.set_xticks(range(len(splits)))
    ax.set_xticklabels(["训练集", "验证集", "测试集"], fontsize=12)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([labels[x] for x in order], fontsize=11)
    ax.set_title("每个划分中，train / val / test 分到了多少组件和酶", fontsize=16, fontweight="bold")

    max_groups = split_df["groups"].max()
    for yi, exp in enumerate(order):
        rows = split_df[split_df["exp"] == exp].set_index("split")
        for xi, split in enumerate(splits):
            row = rows.loc[split]
            size = 120 + 950 * (row["groups"] / max_groups)
            face = BLUE if split != "test" else RED
            ax.scatter(xi, yi, s=size, facecolor=face, edgecolor="white", linewidth=1.6, alpha=0.78)
            ax.text(
                xi + 0.12,
                yi,
                f"{int(row['groups'])}组 / {int(row['enzymes'])}酶",
                ha="left",
                va="center",
                fontsize=9.5,
                color=INK,
            )
    for x in range(len(splits)):
        ax.axvline(x, color=GRID, linewidth=0.8, zorder=0)
    for y in range(len(order)):
        ax.axhline(y, color=GRID, linewidth=0.8, zorder=0)
    ax.text(
        0.5,
        -0.12,
        "圆越大表示组件数越多；测试集用红色强调。EXP010/EXP011 的测试集组件明显减少。",
        transform=ax.transAxes,
        ha="center",
        fontsize=10,
        color=MUTED,
    )
    for spine in ax.spines.values():
        spine.set_visible(False)
    save(fig, "fig03_split_component_bubble_matrix.png")


def plot_performance_landscape(metrics: pd.DataFrame, test_conc: pd.DataFrame) -> None:
    order = ["EXP008", "EXP004", "EXP007", "EXP010", "EXP005", "EXP011", "EXP006"]
    labels = {
        "EXP008": "随机",
        "EXP004": "id60\ncluster",
        "EXP007": "NN80",
        "EXP010": "NN60\n7:1.5:1.5",
        "EXP005": "NN60\n2:1:1",
        "EXP011": "NN60\n8:1:1",
        "EXP006": "NN40",
    }
    df = metrics.set_index("exp").loc[order].reset_index()
    comp_map = test_conc.set_index("exp")["test_groups"].to_dict()
    df["test_groups"] = df["exp"].map(comp_map).fillna(0)
    x = np.arange(len(df))

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.plot(x, df["test_auc"], color=BLUE, linewidth=2.4, marker="o", markersize=7, label="Test AUC")
    ax.plot(x, df["test_aupr"], color=RED, linewidth=2.0, marker="o", markersize=6, label="Test AUPR")
    for i, row in df.iterrows():
        if row["test_groups"] > 0:
            size = 60 + row["test_groups"] * 2.2
            ax.scatter(i, row["test_auc"], s=size, facecolor="none", edgecolor=BLUE, linewidth=1.4, alpha=0.75)
            ax.text(i, row["test_auc"] + 0.035, f"{row['test_auc']:.3f}", ha="center", fontsize=9, color=INK)
            ax.text(i, 0.035, f"test {int(row['test_groups'])}组", ha="center", fontsize=8.5, color=MUTED)
        else:
            ax.text(i, row["test_auc"] + 0.035, f"{row['test_auc']:.3f}", ha="center", fontsize=9, color=INK)
            ax.text(i, 0.035, "随机", ha="center", fontsize=8.5, color=MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels([labels[e] for e in df["exp"]], fontsize=11)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("测试集指标")
    ax.set_title("约束越接近远缘泛化，测试表现整体下降", fontsize=16, fontweight="bold")
    ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False, loc="upper right")
    ax.text(
        0.02,
        0.90,
        "空心圆大小表示 test 集组件数；AUPR 受正样本数量和测试集组成影响更大。",
        transform=ax.transAxes,
        fontsize=10,
        color=MUTED,
    )
    save(fig, "fig04_performance_landscape.png")


def main() -> None:
    setup()
    component_sizes = pd.read_csv(DATA_DIR / "component_size_long.csv")
    split_df = pd.read_csv(DATA_DIR / "split_group_distribution.csv")
    test_conc = pd.read_csv(DATA_DIR / "test_concentration.csv")
    metrics = pd.read_csv(DATA_DIR / "test_metrics.csv")
    plot_concept_network()
    plot_component_size_spectrum(component_sizes)
    plot_split_bubble_matrix(split_df)
    plot_performance_landscape(metrics, test_conc)
    print(f"Wrote v2 figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
