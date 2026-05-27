#!/usr/bin/env python3
"""Draw one presentation-style Q2 sequence-similarity overview figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "q02_sequence_similarity_figures" / "figures_20260525_v3"


INK = "#1F2933"
MUTED = "#667085"
LINE = "#D7DCE2"
PANEL = "#FBFCFD"
BLUE = "#2E5E8C"
BLUE_SOFT = "#EAF1F7"
RED = "#B85252"
RED_SOFT = "#FFF2F2"
GREEN = "#4F6F52"
GREEN_SOFT = "#EEF4EE"
PURPLE = "#66517A"
PURPLE_SOFT = "#F4F0F7"
GRAY = "#6B7280"


def setup() -> None:
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Arial Unicode MS",
        "DejaVu Sans",
    ]
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 150


def add_panel(ax: plt.Axes, x: float, y: float, w: float, h: float, title: str, subtitle: str) -> None:
    panel = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=PANEL,
        edgecolor=LINE,
        linewidth=1.4,
    )
    ax.add_patch(panel)
    ax.text(x + 0.03, y + h - 0.045, title, ha="left", va="top", fontsize=14.5, color=INK, fontweight="bold")
    ax.text(x + 0.03, y + h - 0.087, subtitle, ha="left", va="top", fontsize=9.2, color=MUTED)


def add_node(ax: plt.Axes, x: float, y: float, label: str, edge: str = BLUE, fill: str = "white", r: float = 0.021) -> None:
    circ = patches.Circle((x, y), r, facecolor=fill, edgecolor=edge, linewidth=1.8)
    ax.add_patch(circ)
    ax.text(x, y, label, ha="center", va="center", fontsize=9.2, color=INK, fontweight="bold")


def add_split_box(ax: plt.Axes, x: float, y: float, text: str, color: str, w: float = 0.09) -> None:
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        0.042,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor="white",
        edgecolor=color,
        linewidth=1.3,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + 0.021, text, ha="center", va="center", fontsize=9.2, color=color, fontweight="bold")


def connect(ax: plt.Axes, a: tuple[float, float], b: tuple[float, float], color: str, dashed: bool = False, lw: float = 2.0) -> None:
    line = patches.FancyArrowPatch(
        a,
        b,
        arrowstyle="-",
        connectionstyle="arc3,rad=0.05",
        color=color,
        linewidth=lw,
        linestyle="--" if dashed else "-",
        alpha=0.92,
    )
    ax.add_patch(line)


def callout(ax: plt.Axes, x: float, y: float, w: float, h: float, text: str, color: str, fill: str) -> None:
    box = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.01,rounding_size=0.014",
        facecolor=fill,
        edgecolor=color,
        linewidth=1.4,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10.2, color=color, fontweight="bold")


def draw_random_panel(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    add_panel(ax, x, y, w, h, "随机划分", "相似酶可能同时出现在训练集和测试集")
    add_split_box(ax, x + 0.045, y + 0.265, "train", BLUE, w=0.10)
    add_split_box(ax, x + w - 0.145, y + 0.265, "test", RED, w=0.10)
    add_node(ax, x + 0.095, y + 0.185, "E1")
    add_node(ax, x + 0.17, y + 0.285, "E2")
    add_node(ax, x + w - 0.095, y + 0.205, "E3", edge=RED)
    connect(ax, (x + 0.17, y + 0.285), (x + w - 0.095, y + 0.205), RED, dashed=True)
    ax.text(x + w / 2, y + 0.335, "相似酶可能跨集合", ha="center", fontsize=10.5, color=RED, fontweight="bold")
    callout(ax, x + 0.055, y + 0.055, w - 0.11, 0.07, "Test AUC 0.934", BLUE, BLUE_SOFT)


def draw_cluster_panel(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    add_panel(ax, x, y, w, h, "EXP004：id60 聚类簇划分", "同一 cluster 不拆开，但跨 cluster 近邻仍可能存在")
    c1 = patches.FancyBboxPatch((x + 0.045, y + 0.185), 0.13, 0.145, boxstyle="round,pad=0.01,rounding_size=0.012", facecolor=BLUE_SOFT, edgecolor=BLUE, linestyle="--", linewidth=1.5)
    c2 = patches.FancyBboxPatch((x + w - 0.175, y + 0.185), 0.13, 0.145, boxstyle="round,pad=0.01,rounding_size=0.012", facecolor=BLUE_SOFT, edgecolor=BLUE, linestyle="--", linewidth=1.5)
    ax.add_patch(c1)
    ax.add_patch(c2)
    ax.text(x + 0.11, y + 0.335, "cluster 1", ha="center", fontsize=9.0, color=BLUE)
    ax.text(x + w - 0.11, y + 0.335, "cluster 2", ha="center", fontsize=9.0, color=BLUE)
    add_node(ax, x + 0.085, y + 0.245, "E1")
    add_node(ax, x + 0.135, y + 0.245, "E2")
    add_node(ax, x + w - 0.11, y + 0.245, "E3")
    connect(ax, (x + 0.145, y + 0.245), (x + w - 0.11, y + 0.245), RED, dashed=True)
    ax.text(x + w / 2, y + 0.385, "跨 cluster 仍可相似", ha="center", va="center", fontsize=10.5, color=RED, fontweight="bold")
    add_split_box(ax, x + 0.07, y + 0.145, "test", RED, w=0.07)
    add_split_box(ax, x + w - 0.145, y + 0.145, "train", BLUE, w=0.08)
    callout(ax, x + 0.055, y + 0.055, w - 0.11, 0.07, "train-test NN60：19 对", RED, RED_SOFT)


def draw_strict_panel(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    add_panel(ax, x, y, w, h, "EXP005：strict NN60", "把所有 NN60 近邻连接合并成组件，再整体分配")
    comp = patches.FancyBboxPatch((x + 0.06, y + 0.185), w - 0.12, 0.165, boxstyle="round,pad=0.01,rounding_size=0.012", facecolor=PURPLE_SOFT, edgecolor=PURPLE, linestyle="--", linewidth=1.7)
    ax.add_patch(comp)
    ax.text(x + w / 2, y + 0.35, "NN60 组件", ha="center", fontsize=9.5, color=PURPLE, fontweight="bold")
    pts = [(x + 0.12, y + 0.25, "E1"), (x + 0.19, y + 0.25, "E2"), (x + w - 0.13, y + 0.25, "E3")]
    for px, py, lab in pts:
        add_node(ax, px, py, lab, edge=PURPLE)
    connect(ax, (x + 0.14, y + 0.25), (x + 0.17, y + 0.25), RED)
    connect(ax, (x + 0.21, y + 0.25), (x + w - 0.15, y + 0.25), RED)
    ax.text(x + w / 2, y + 0.385, "所有 NN60 连接合并", ha="center", va="center", fontsize=10.5, color=PURPLE, fontweight="bold")
    callout(ax, x + 0.055, y + 0.055, w - 0.11, 0.07, "train-test NN60：0 对", PURPLE, PURPLE_SOFT)


def draw_metric_axis(ax: plt.Axes) -> None:
    x0, x1, y = 0.08, 0.92, 0.12
    ax.plot([x0, x1], [y, y], color=LINE, linewidth=2)
    entries = [
        ("随机", 0.934, 0.08),
        ("id60 cluster", 0.808, 0.25),
        ("NN80", 0.820, 0.40),
        ("NN60", 0.671, 0.58),
        ("NN40", 0.638, 0.78),
        ("更严格远缘", None, 0.92),
    ]
    for label, auc, xpos in entries:
        color = RED if label in {"NN60", "NN40"} else BLUE if auc is not None else MUTED
        ax.scatter(xpos, y, s=160 if auc is not None else 60, facecolor="white", edgecolor=color, linewidth=2.0, zorder=3)
        if auc is not None:
            ax.text(xpos, y + 0.045, f"{auc:.3f}", ha="center", fontsize=11, color=INK, fontweight="bold")
        ax.text(xpos, y - 0.055, label, ha="center", va="top", fontsize=10.5, color=color)
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="->", color=MUTED, linewidth=1.8))
    ax.text(0.50, 0.02, "序列隔离越严格，测试更接近远缘酶泛化；随机划分的高分不应直接代表远缘泛化能力", ha="center", fontsize=11, color=INK)


def main() -> None:
    setup()
    fig = plt.figure(figsize=(14.8, 8.3), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.5, 0.955, "Q2 序列相似度划分：模型真的能泛化到远缘酶吗？", ha="center", va="top", fontsize=22, color=INK, fontweight="bold")
    ax.text(0.5, 0.905, "核心变化：从随机样本划分，转为按酶序列相似关系整体划分", ha="center", va="top", fontsize=13, color=MUTED)

    draw_random_panel(ax, 0.055, 0.38, 0.27, 0.42)
    draw_cluster_panel(ax, 0.365, 0.38, 0.27, 0.42)
    draw_strict_panel(ax, 0.675, 0.38, 0.27, 0.42)
    draw_metric_axis(ax)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / "fig01_q2_sequence_split_story_overview.png"
    fig.savefig(out, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
