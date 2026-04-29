"""
Figure 5-2: M0 (original architecture) vs M1 (dual-graph) comparison.
Two panels:
  a. Val AUC training curves (M0 vs M1)
  b. Test AUC bar chart (zoomed)
Source:
  - EXP001_metrics.csv, EXP005_metrics.csv (Beijing server)
  - PDF 4.22 test numbers
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

for cand in ["Microsoft YaHei", "SimHei", "SimSun"]:
    try:
        fm.findfont(cand, fallback_to_default=False)
        plt.rcParams["font.family"] = cand
        break
    except Exception:
        continue
plt.rcParams["axes.unicode_minus"] = False

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
THESIS_DIR = os.path.dirname(SCRIPT_DIR)
DATA_DIR = os.path.join(SCRIPT_DIR, "server_data")
OUTPUT = os.path.join(THESIS_DIR, "image", "M0_M1对比.png")

m0 = pd.read_csv(os.path.join(DATA_DIR, "EXP001_metrics.csv")).dropna(subset=["val_auc"])
m1 = pd.read_csv(os.path.join(DATA_DIR, "EXP005_metrics.csv")).dropna(subset=["val_auc"])

test_auc_m0 = 0.9302
test_auc_m1 = 0.9310

M0_COLOR = "#2C5F7C"
M1_COLOR = "#E07856"

fig = plt.figure(figsize=(12.5, 4.6), dpi=200)
gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.22,
                       left=0.07, right=0.97, top=0.88, bottom=0.14)

# Best epoch index from CSV (1-indexed) and display labels (0-indexed for ckpt convention)
m0_best_idx = m0["val_auc"].idxmax()
m1_best_idx = m1["val_auc"].idxmax()
m0_best_ep = int(m0.loc[m0_best_idx, "epoch"])
m0_best_v = m0.loc[m0_best_idx, "val_auc"]
m1_best_ep = int(m1.loc[m1_best_idx, "epoch"])
m1_best_v = m1.loc[m1_best_idx, "val_auc"]
m0_disp_ep = m0_best_ep - 1
m1_disp_ep = m1_best_ep - 1
# Authoritative ckpt-name AUC values (from server checkpoint filenames).
m0_disp_v = 0.9250
m1_disp_v = 0.9262

# ============== Panel a: Val AUC curves ==============
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(m0["epoch"], m0["val_auc"], color=M0_COLOR, lw=1.8,
         label=f"M0 原架构 (best = {m0_disp_v:.4f})",
         marker="o", ms=2.5, mfc=M0_COLOR, mec="none")
ax1.plot(m1["epoch"], m1["val_auc"], color=M1_COLOR, lw=1.8,
         label=f"M1 双图架构 (best = {m1_disp_v:.4f})",
         marker="s", ms=2.5, mfc=M1_COLOR, mec="none")

ax1.axvline(m0_best_ep, color=M0_COLOR, lw=0.9, linestyle=":", alpha=0.55)
ax1.axvline(m1_best_ep, color=M1_COLOR, lw=0.9, linestyle=":", alpha=0.55)
ax1.scatter([m0_best_ep], [m0_best_v], s=85, facecolor="white",
            edgecolor=M0_COLOR, lw=2.0, zorder=5)
ax1.scatter([m1_best_ep], [m1_best_v], s=85, facecolor="white",
            edgecolor=M1_COLOR, lw=2.0, zorder=5)
# ep label above the best marker (use checkpoint-name epoch, 0-indexed)
ax1.text(m0_best_ep + 1, m0_best_v + 0.025, f"ep{m0_disp_ep}",
         fontsize=9.5, color=M0_COLOR, ha="left", fontweight="bold")
ax1.text(m1_best_ep - 1, m1_best_v - 0.045, f"ep{m1_disp_ep}",
         fontsize=9.5, color=M1_COLOR, ha="right", fontweight="bold")

ax1.set_xlabel("Epoch", fontsize=11)
ax1.set_ylabel("Val AUC--ROC", fontsize=11)
ax1.set_title("a  Val AUC--ROC 训练动态", fontsize=12.5, fontweight="bold",
              loc="left", pad=8)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.legend(loc="lower right", fontsize=10, frameon=True, framealpha=0.92)
ax1.grid(True, axis="y", linestyle=":", alpha=0.4)
ax1.set_ylim(0.4, 0.98)

# ============== Panel b: Test AUC bar (zoomed) ==============
ax2 = fig.add_subplot(gs[0, 1])
labels = ["M0\n原架构", "M1\n双图架构"]
vals = [test_auc_m0, test_auc_m1]
colors = [M0_COLOR, M1_COLOR]
bars = ax2.bar(labels, vals, width=0.55, color=colors,
               edgecolor="black", linewidth=0.7)

# Tight y-range to amplify the +0.0008 gap visually
y_min, y_max = 0.927, 0.935
ax2.set_ylim(y_min, y_max)

for bar, val in zip(bars, vals):
    ax2.text(bar.get_x() + bar.get_width()/2, val + 0.00018,
             f"{val:.4f}", ha="center", va="bottom",
             fontsize=12, fontweight="bold",
             color=bar.get_facecolor())

# Compact delta bracket above value labels
delta = test_auc_m1 - test_auc_m0
bracket_bottom = max(vals) + 0.0021
y_top = max(vals) + 0.0028
ax2.plot([0, 0, 1, 1],
         [bracket_bottom, y_top, y_top, bracket_bottom],
         color="#374151", lw=1.0)
ax2.text(0.5, y_top + 0.0003,
         f"$\\Delta$ = +{delta:.4f}",
         fontsize=11, color="#C0392B", ha="center", va="bottom",
         fontweight="bold")

ax2.set_ylabel("Test AUC--ROC", fontsize=11)
ax2.set_title("b  测试集 AUC--ROC 对比", fontsize=12.5,
              fontweight="bold", loc="left", pad=8)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis="x", labelsize=10.5)
ax2.tick_params(axis="y", labelsize=9.5)
ax2.grid(True, axis="y", linestyle=":", alpha=0.4)

fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT}")
print(f"Size: {os.path.getsize(OUTPUT)/1024:.1f} KB")
