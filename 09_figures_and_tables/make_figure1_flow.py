"""Figure 1 — STARD-AI participant flow diagram.

Draws the screening/exclusion flow of the MGMT multivendor radiomics study
from the counts in supplementary_appendices/Appendix_B_reporting_checklists.md
(section B2) and saves it as the Figure 1 file of the main figure set.

Output: paper1/figures/Figure1_participant_flow_STARD_AI.png (300 dpi)

Style follows the other main figures: matplotlib defaults, PNG at 300 dpi,
drawn at 7.0 in width so that 8 pt type is not scaled down when the figure
is placed on the page at 16.5 cm.
"""

import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   "..", "..", "figures")
OUT = os.path.join(FIG, "Figure1_participant_flow_STARD_AI.png")

plt.rcParams.update({'font.size': 8.0, 'savefig.dpi': 300,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})

EDGE = "#333333"
FILL_MAIN = "#ffffff"
FILL_EXCL = "#f0f0f0"
FILL_FINAL = "#dbe9f6"

W_IN, H_IN = 7.0, 5.2
FS = 8.0

fig = plt.figure(figsize=(W_IN, H_IN))
ax = fig.add_axes((0, 0, 1, 1))       # axes fill the figure: 1 x-unit = W_IN/100 in
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

UY = H_IN / 100.0                     # one y-unit in inches


def box_height(text, fs=FS, pad=0.16):
    """Box height in y-units, auto-sized to the number of text lines."""
    return (len(text.split("\n")) * fs / 72.0 * 1.5 + pad) / UY


def box(x, y, w, text, fill=FILL_MAIN, bold=False, fs=FS):
    ax.add_patch(plt.Rectangle((x - w / 2, y - box_height(text, fs) / 2), w,
                               box_height(text, fs),
                               facecolor=fill, edgecolor=EDGE, linewidth=0.9))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", linespacing=1.5)


def arrow(x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                                 arrowstyle="-|>", mutation_scale=9,
                                 linewidth=0.9, color=EDGE,
                                 shrinkA=0, shrinkB=0))


# ---- texts ----------------------------------------------------------------
t1 = "109 patients screened:\n" \
     "histopathologically confirmed glioblastoma and\n" \
     "archived preoperative brain MRI, July 2019 to May 2026"
t2 = "93 patients\n" \
     "underwent MGMT promoter methylation testing\n" \
     "of the resection or biopsy specimen"
t3 = "90 patients\n" \
     "with a complete preoperative three-sequence MRI\n" \
     "(T1-weighted, contrast-enhanced T1-weighted, FLAIR)"
t4 = "89 patients analysed\n" \
     "MGMT methylated 45 / unmethylated 44\n" \
     "Siemens 57 / GE 32; all 3.0 T; all IDH-wildtype"

e1 = "16 patients excluded:\n" \
     "no MGMT promoter\n" \
     "methylation result available"
e2 = "3 patients excluded:\n" \
     "preoperative three-sequence\n" \
     "MRI unavailable or unusable"
e3 = "1 patient excluded:\n" \
     "single 1.5 T United Imaging\n" \
     "examination; no platform-level\n" \
     "inference possible"

# ---- geometry -------------------------------------------------------------
LX, LW = 31, 58          # main-flow column: centre x, box width
RX, RW = 81, 34          # exclusion column
GAP = 8                  # vertical gap between main boxes (y-units)

H1, H2, H3, H4 = (box_height(t) for t in (t1, t2, t3, t4))
HE1, HE2, HE3 = (box_height(t) for t in (e1, e2, e3))

Y4 = H4 / 2 + 4
Y3 = Y4 + H4 / 2 + GAP + H3 / 2
Y2 = Y3 + H3 / 2 + GAP + H2 / 2
Y1 = Y2 + H2 / 2 + GAP + H1 / 2

# ---- main flow ------------------------------------------------------------
box(LX, Y1, LW, t1, bold=True)
box(LX, Y2, LW, t2)
box(LX, Y3, LW, t3)
box(LX, Y4, LW, t4, fill=FILL_FINAL, bold=True)

# ---- exclusions -----------------------------------------------------------
box(RX, (Y1 + Y2) / 2, RW, e1, fill=FILL_EXCL)
box(RX, (Y2 + Y3) / 2, RW, e2, fill=FILL_EXCL)
box(RX, (Y3 + Y4) / 2, RW, e3, fill=FILL_EXCL)

# ---- arrows ---------------------------------------------------------------
for ya, yb, ha, hb in [(Y1, Y2, H1, H2), (Y2, Y3, H2, H3), (Y3, Y4, H3, H4)]:
    arrow(LX, ya - ha / 2, LX, yb + hb / 2)
    ymid = (ya - ha / 2 + yb + hb / 2) / 2
    arrow(LX + LW / 2, ymid, RX - RW / 2, ymid)

# crop the canvas to the content with a small margin
top = Y1 + H1 / 2 + 2
bot = Y4 - H4 / 2 - 2
ax.set_ylim(bot, top)

os.makedirs(FIG, exist_ok=True)
fig.savefig(OUT)
print("saved", os.path.abspath(OUT))
