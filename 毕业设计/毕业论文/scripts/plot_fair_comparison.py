"""
Figure 4-1: Paper baseline vs Ours under leakage-controlled fair comparison.
Two panels: AUC-ROC and AUPR; each shows two pairs of bars
(filtered vs unfiltered test set; paper vs ours).
Source: 组会 4.22 PDF table.
"""
import os
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

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
OUTPUT = os.path.join(THESIS_DIR, "image", "公平对比.png")

# Numbers per the meeting PDF (4.22), four-decimal precision
groups = ["过滤后\n(n=7,963)", "未过滤\n(n=10,999)"]
paper_auc = [0.5586, 0.5860]
ours_auc = [0.9205, 0.9302]
paper_aupr = [0.1004, 0.1124]
ours_aupr = [0.6403, 0.6749]

base_rate = [680 / 7963, 984 / 10999]  # 0.0854, 0.0895; mean 0.0875

PAPER_COLOR = "#9CA3AF"
OURS_COLOR = "#2C5F7C"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 5.4), dpi=200,
                                gridspec_kw={"wspace": 0.22})

x = np.arange(len(groups))
bar_w = 0.36

# ============== Panel a: AUC-ROC ==============
b1 = ax1.bar(x - bar_w/2, paper_auc, bar_w, color=PAPER_COLOR, edgecolor="black",
             linewidth=0.7, label="参考 ckpt")
b2 = ax1.bar(x + bar_w/2, ours_auc, bar_w, color=OURS_COLOR, edgecolor="black",
             linewidth=0.7, label="本研究 M0")
ax1.axhline(0.5, color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
ax1.text(1.45, 0.515, "随机基线 0.5", fontsize=9, color="gray", ha="right")
ax1.set_ylim(0, 1.05)
ax1.set_xticks(x)
ax1.set_xticklabels(groups, fontsize=10.5)
ax1.set_ylabel("AUC--ROC", fontsize=12)
ax1.set_title("a  AUC--ROC", fontsize=13, fontweight="bold", loc="left", pad=8)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis="y", labelsize=10)

for bar, val in zip(b1, paper_auc):
    ax1.text(bar.get_x() + bar.get_width()/2, val + 0.020,
             f"{val:.4f}", ha="center", va="bottom",
             fontsize=11, color="#4B5563")
for bar, val in zip(b2, ours_auc):
    ax1.text(bar.get_x() + bar.get_width()/2, val + 0.020,
             f"{val:.4f}", ha="center", va="bottom",
             fontsize=11.5, color=OURS_COLOR, fontweight="bold")

def _delta_badge(ax, x, y, delta, fontsize=10, color="white", bg="#C0392B"):
    """Compact badge above an 'ours' bar showing the improvement."""
    ax.text(x, y, f"$\\uparrow$ +{delta:.4f}",
            fontsize=fontsize, color=color, ha="center", va="bottom",
            fontweight="bold",
            bbox=dict(facecolor=bg, edgecolor="none",
                      boxstyle="round,pad=0.30"))

# delta values are now inline above each ours bar


# ============== Panel b: AUPR ==============
b3 = ax2.bar(x - bar_w/2, paper_aupr, bar_w, color=PAPER_COLOR, edgecolor="black",
             linewidth=0.7, label="参考 ckpt")
b4 = ax2.bar(x + bar_w/2, ours_aupr, bar_w, color=OURS_COLOR, edgecolor="black",
             linewidth=0.7, label="本研究 M0")
ax2.axhline(np.mean(base_rate), color="gray", linestyle="--", linewidth=1.0, alpha=0.6)
ax2.text(1.45, np.mean(base_rate) + 0.012,
         f"阳性基线 $\\approx$ {np.mean(base_rate):.4f}",
         fontsize=9, color="gray", ha="right")
ax2.set_ylim(0, 0.78)
ax2.set_xticks(x)
ax2.set_xticklabels(groups, fontsize=10.5)
ax2.set_ylabel("AUPR", fontsize=12)
ax2.set_title("b  AUPR", fontsize=13, fontweight="bold", loc="left", pad=8)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(axis="y", labelsize=10)

for bar, val in zip(b3, paper_aupr):
    ax2.text(bar.get_x() + bar.get_width()/2, val + 0.016,
             f"{val:.4f}", ha="center", va="bottom",
             fontsize=11, color="#4B5563")
for bar, val in zip(b4, ours_aupr):
    ax2.text(bar.get_x() + bar.get_width()/2, val + 0.016,
             f"{val:.4f}", ha="center", va="bottom",
             fontsize=11.5, color=OURS_COLOR, fontweight="bold")

# delta values are now inline above each ours bar


# Shared legend at top
handles = [b1[0], b2[0]]
labels_l = ["参考 ckpt", "本研究 M0"]
fig.legend(handles, labels_l, loc="upper center", ncol=2, fontsize=11,
           frameon=False, bbox_to_anchor=(0.5, 1.005))

# Bottom strip
fig.text(0.5, -0.02,
         "测试集说明：过滤后 = 剔除 ESIBank 训练集中出现过的 P450 酶；未过滤 = 全部 P450 测试样本",
         ha="center", fontsize=10, color="#374151")

plt.subplots_adjust(top=0.88)
fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT}")
print(f"Size: {os.path.getsize(OUTPUT)/1024:.1f} KB")
