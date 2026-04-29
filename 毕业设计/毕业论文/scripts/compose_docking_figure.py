"""
Compose Figure 3-3: overview + zoom-in for docking complex.
Reads docking_overview.png and docking_zoom.png, produces 对接复合物.png.
"""
import os
import sys
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import ConnectionPatch
import matplotlib.font_manager as fm

# Chinese font
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
# Inputs may be in script dir (after copy from /c/temp) or co-located
OVERVIEW = os.path.join(SCRIPT_DIR, "docking_overview.png")
ZOOM = os.path.join(SCRIPT_DIR, "docking_zoom.png")
if not os.path.exists(OVERVIEW):
    OVERVIEW = os.path.join(os.getcwd(), "docking_overview.png")
    ZOOM = os.path.join(os.getcwd(), "docking_zoom.png")
OUTPUT = os.path.join(THESIS_DIR, "image", "对接复合物.png")

img_o = Image.open(OVERVIEW)
img_z = Image.open(ZOOM)

fig = plt.figure(figsize=(14, 6.5), dpi=200)
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.04, left=0.02, right=0.98, top=0.93, bottom=0.04)

# Left: overview
ax_o = fig.add_subplot(gs[0, 0])
ax_o.imshow(img_o)
ax_o.set_xticks([])
ax_o.set_yticks([])
for s in ax_o.spines.values():
    s.set_visible(False)
ax_o.text(0.02, 0.98, "a  整体结构", transform=ax_o.transAxes,
          fontsize=15, fontweight="bold", va="top", ha="left",
          bbox=dict(facecolor="white", edgecolor="none", pad=2))

# Box on overview indicating active site (in image pixel coords)
ow, oh = img_o.size
box_cx, box_cy = ow * 0.50, oh * 0.43
box_w, box_h = ow * 0.18, oh * 0.18
rect = patches.Rectangle(
    (box_cx - box_w / 2, box_cy - box_h / 2),
    box_w, box_h,
    linewidth=1.6, edgecolor="black", facecolor="none", linestyle=(0, (4, 3)))
ax_o.add_patch(rect)

# Right: zoom
ax_z = fig.add_subplot(gs[0, 1])
ax_z.imshow(img_z)
ax_z.set_xticks([])
ax_z.set_yticks([])
for s in ax_z.spines.values():
    s.set_visible(True)
    s.set_linewidth(1.6)
    s.set_linestyle((0, (4, 3)))
    s.set_edgecolor("black")
ax_z.text(0.02, 0.98, "b  活性口袋放大", transform=ax_z.transAxes,
          fontsize=15, fontweight="bold", va="top", ha="left",
          bbox=dict(facecolor="white", edgecolor="none", pad=2))

# Dashed connector lines (top-right and bottom-right of box -> top-left and bottom-left of zoom)
con_top = ConnectionPatch(
    xyA=(box_cx + box_w / 2, box_cy - box_h / 2), coordsA=ax_o.transData,
    xyB=(0, 0), coordsB=ax_z.transAxes,
    color="black", linewidth=1.2, linestyle=(0, (4, 3)))
con_bot = ConnectionPatch(
    xyA=(box_cx + box_w / 2, box_cy + box_h / 2), coordsA=ax_o.transData,
    xyB=(0, 1), coordsB=ax_z.transAxes,
    color="black", linewidth=1.2, linestyle=(0, (4, 3)))
fig.add_artist(con_top)
fig.add_artist(con_bot)

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
fig.savefig(OUTPUT, dpi=300, bbox_inches="tight", facecolor="white")
print(f"Saved: {OUTPUT}")
print(f"Size: {os.path.getsize(OUTPUT)/1024:.1f} KB")
