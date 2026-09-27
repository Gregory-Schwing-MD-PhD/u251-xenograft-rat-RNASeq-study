# -*- coding: utf-8 -*-
"""Study design in the form of a journal Figure 1a: (a) procedure per arm, (b) where tissue was taken, (c) sample matrix.

Source of truth (Nagaraja lab records, 2023): the RNA-seq sample key (metadata.csv: IL-64, IL-67, IL-68, IL-69 Primary;
IL-65, IL-66, IL-70, IL-71 Recurrent; N1, N2 Control 'Contralateral'; C1, C2 culture) and the methylation-array sample sheet
'1724 (Raj-8).xlsx' (IL66B..IL71B, '69B control N2', C2B). Procedure and timing from Nagaraja et al., J Neurosurg 2026 (LITT
when the tumour reached about 4 mm; recurrence imaged at 2 and 4 weeks). Harvest time per animal is not recorded, so no
per-animal time axis is drawn. Human-read shares: xengsort graft fractions (SLIDES/figures_cns/chart_sorting.csv).
Methylation marks: the EPIC v1 chip 205648300021 as listed in ANALYSIS/methylation/assets/samplesheet.csv (eight arrays,
each joined to its RNA library by the rna_library column). ANALYSIS/SAMPLE_KEY.md has the table and what is not known.

    python SLIDES/12_design_figure.py      -> SLIDES/figures_cns/fig_design.png
"""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

from matplotlib.patches import Ellipse, FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon, Rectangle
from matplotlib.path import Path as MPath

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("cns04", HERE / "04_make_cns_figures.py")
f4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f4)
plt, save, NAVY, GREY = f4.plt, f4.save, f4.NAVY, f4.GREY
ORANGE, INKK, RULE, FILL = "#E08214", "#1A1A1A", "#9A9A9A", "#F3F3F3"

graft = {r["sample"]: float(r["graft"]) for r in csv.DictReader(open(HERE / "figures_cns" / "chart_sorting.csv", encoding="utf-8"))}
assert graft["IL64B"] < 1 and graft["N168B"] < 1 and 4 < graft["N269B"] < 6
assert min(graft[s] for s in ("IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B")) > 25

sheet = list(csv.DictReader(open(HERE.parent / "ANALYSIS" / "methylation" / "assets" / "samplesheet.csv", encoding="utf-8")))
assert len(sheet) == 8 and {r["sentrix_id"] for r in sheet} == {"205648300021"}
array_of = {r["rna_library"]: r["sample"] for r in sheet}          # RNA library -> array name
assert set(array_of) <= set(graft), set(array_of) - set(graft)

W = 14.3                                                            # 13.2 before the two assay columns; a and b keep their size
fig = plt.figure(figsize=(W, 4.7))


def ax_at(x_in, w_in):
    return fig.add_axes([x_in / W, 0.0, w_in / W, 1.0])


def letter(ax, s):
    ax.text(0.0, 1.0, s, transform=ax.transAxes, ha="left", va="top", fontsize=17, fontweight="bold", color=INKK)


def box(ax, x, y, w, h, text, edge=RULE, face=FILL, color=INKK, bold=False, size=11.5):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.018",
                                fc=face, ec=edge, lw=1.0, transform=ax.transAxes))
    ax.text(x, y, text, transform=ax.transAxes, ha="center", va="center", fontsize=size, color=color,
            fontweight="bold" if bold else "normal", linespacing=1.15)


def arrow(ax, x0, y0, x1, y1, color=INKK):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=11,
                                 lw=1.1, color=color, shrinkA=0, shrinkB=0))


def rat(ax, x, y, s=1.0, color="#4D4D4D"):
    """A flat side-view rat, facing right, centred near (x, y) in axes units."""
    t = ax.transAxes
    ax.add_patch(Ellipse((x, y), 0.105 * s, 0.13 * s, fc=color, ec="none", transform=t))                   # body
    ax.add_patch(Ellipse((x + 0.058 * s, y + 0.012 * s), 0.06 * s, 0.075 * s, fc=color, ec="none", transform=t))  # head
    ax.add_patch(Polygon([(x + 0.075 * s, y + 0.03 * s), (x + 0.108 * s, y + 0.004 * s), (x + 0.078 * s, y - 0.012 * s)],
                         closed=True, fc=color, ec="none", transform=t))                                        # snout
    ax.add_patch(Ellipse((x + 0.05 * s, y + 0.058 * s), 0.024 * s, 0.035 * s, fc=color, ec="white", lw=0.8, transform=t))  # ear
    ax.add_patch(Ellipse((x + 0.074 * s, y + 0.02 * s), 0.006 * s, 0.009 * s, fc="white", ec="none", transform=t))  # eye
    verts = [(x - 0.05 * s, y - 0.01 * s), (x - 0.1 * s, y - 0.06 * s), (x - 0.14 * s, y + 0.02 * s), (x - 0.12 * s, y + 0.06 * s)]
    ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]), fc="none", ec=color,
                           lw=1.6, transform=t))                                                               # tail
    for dx in (-0.025, 0.02):                                                                                  # feet
        ax.add_patch(Ellipse((x + dx * s, y - 0.062 * s), 0.02 * s, 0.012 * s, fc=color, ec="none", transform=t))


# ---------------------------------------------------------------- a  procedure
a = ax_at(0.0, 6.6); a.axis("off"); letter(a, "a")
rat(a, 0.075, 0.52, 0.95)
a.text(0.07, 0.38, "athymic\nRNU/RNU rat", transform=a.transAxes, ha="center", va="top", fontsize=10.5, color=GREY)
box(a, 0.31, 0.50, 0.19, 0.2, "U251N cells\nimplanted in\none hemisphere")
box(a, 0.52, 0.50, 0.15, 0.2, "MRI:\ntumour\n≈ 4 mm")
arrow(a, 0.18, 0.50, 0.215, 0.50)
arrow(a, 0.405, 0.50, 0.445, 0.50)
yp, yr = 0.80, 0.20
xs, xe = 0.595, 0.612                                             # square elbow from the MRI box into each arm
a.plot([xs, xe], [0.50, 0.50], color=INKK, lw=1.1, transform=a.transAxes, solid_capstyle="butt")
for yy, col, x_in in ((yp, NAVY, 0.65), (yr, ORANGE, 0.625)):
    a.plot([xe, xe], [0.50, yy], color=col, lw=1.1, transform=a.transAxes, solid_capstyle="butt")
    arrow(a, xe, yy, x_in, yy, col)
# primary arm
box(a, 0.735, yp, 0.17, 0.15, "no treatment", edge=NAVY, face="white", color=NAVY)
arrow(a, 0.82, yp, 0.875, yp, NAVY)
box(a, 0.935, yp, 0.12, 0.15, "harvest", edge=NAVY, face="white", color=NAVY)
a.text(0.625, yp + 0.12, "Primary   4 rats (64, 67, 68, 69)", transform=a.transAxes, ha="left", va="bottom", fontsize=12,
       color=NAVY, fontweight="bold")
# recurrent arm
box(a, 0.685, yr, 0.12, 0.17, "MRI-guided\nLITT", edge=ORANGE, face="white", color=ORANGE, size=10)
arrow(a, 0.745, yr, 0.765, yr, ORANGE)
box(a, 0.82, yr, 0.11, 0.17, "regrowth\non MRI\n(2–4 wk)", edge=ORANGE, face="white", color=ORANGE, size=10.5)
arrow(a, 0.875, yr, 0.893, yr, ORANGE)
box(a, 0.938, yr, 0.09, 0.17, "harvest", edge=ORANGE, face="white", color=ORANGE, size=11)
a.text(0.625, yr - 0.12, "Recurrent   4 rats (65, 66, 70, 71)", transform=a.transAxes, ha="left", va="top", fontsize=12,
       color=ORANGE, fontweight="bold")

# ---------------------------------------------------------------- b  tissue sampled
b = ax_at(6.666, 2.178); b.axis("off"); letter(b, "b")
b.set_xlim(0, 1); b.set_ylim(0, 1)
cx, cy = 0.5, 0.54                                              # coronal section, true proportions (axes are 2.2 x 4.7 in)
b.add_patch(Ellipse((cx, cy), 0.92, 0.25, fc="#FAFAFA", ec=INKK, lw=1.2))   # whole section, ~1.6 : 1 as in a rat
b.plot([cx, cx], [cy - 0.125, cy + 0.125], color=INKK, lw=1.0)
b.add_patch(Ellipse((cx - 0.2, cy + 0.01), 0.17, 0.075, fc="#5A5A5A", ec="none"))
b.add_patch(Rectangle((cx - 0.34, cy - 0.065), 0.29, 0.14, fill=False, ec=INKK, lw=1.1, ls=(0, (3, 2))))
b.add_patch(Rectangle((cx + 0.05, cy - 0.065), 0.29, 0.14, fill=False, ec=GREY, lw=1.1, ls=(0, (3, 2))))
b.text(cx - 0.215, cy + 0.14, "implanted\nside", ha="center", va="bottom", fontsize=11, color=INKK)
b.text(cx + 0.215, cy + 0.14, "opposite\nside", ha="center", va="bottom", fontsize=11, color=GREY)
b.text(cx - 0.215, cy - 0.14, "tumour\n(all rats)", ha="center", va="top", fontsize=10.5, color=INKK)
b.text(cx + 0.215, cy - 0.14, "rats 68, 69\nonly", ha="center", va="top", fontsize=10.5, color=GREY)
b.text(cx, 0.06, "2-mm coronal slices,\nfrozen, then dissected", ha="center", va="bottom", fontsize=10, color=GREY, style="italic")

# ---------------------------------------------------------------- c  sample matrix
CX, CW = 9.042, W - 9.042                                          # panel c in inches; columns placed in inches below
c = fig.add_axes([CX / W, 0.02, CW / W, 0.96]); c.axis("off"); letter(c, "c")
c.set_xlim(0, CW); c.set_ylim(0, 1)
ROWS = [("IL67B", 67, "tumour", "P", "cmp"), ("IL68B", 68, "tumour", "P", "cmp"), ("IL69B", 69, "tumour", "P", "cmp"),
        ("IL66B", 66, "tumour", "R", "cmp"), ("NL70B", 70, "tumour", "R", "cmp"), ("NL71B", 71, "tumour", "R", "cmp"),
        ("IL64B", 64, "no tumour", "P", "ctl"), ("N168B", 68, "opposite", "P", "ctl"), ("N269B", 69, "opposite", "P", "ctl"),
        ("C2B", "–", "culture", "", "ref")]
assert set(array_of) <= {r[0] for r in ROWS}
colx = {"lib": 0.08, "rat": 1.04, "tis": 1.50, "hum": 2.74, "rna": 3.13, "dna": 3.55, "cmp": 4.13, "ctl": 4.57, "ref": 5.0}
L, R_ = 0.05, CW - 0.03
top, step = 0.80, 0.072
heads = [("lib", "sample", "left"), ("rat", "rat", "center"), ("tis", "tissue", "left"), ("hum", "human\nreads", "right"),
         ("rna", "RNA-\nseq", "center"), ("dna", "EPIC\narray", "center"),
         ("cmp", "3 v 3", "center"), ("ctl", "control", "center"), ("ref", "ref.", "center")]
for k, h, ha in heads:
    c.text(colx[k], top + 0.045, h, ha=ha, va="bottom", fontsize=10.5, color=GREY, linespacing=1.0)
for x0, x1, h in ((colx["rna"] - 0.2, colx["dna"] + 0.2, "assay"), (colx["cmp"] - 0.2, colx["ref"] + 0.16, "use")):
    c.text((x0 + x1) / 2, top + 0.13, h, ha="center", va="bottom", fontsize=10.5, color=GREY)
    c.plot([x0, x1], [top + 0.122, top + 0.122], color=RULE, lw=0.7)
c.plot([L, R_], [top + 0.035, top + 0.035], color=INKK, lw=0.9)
S = 0.18                                                            # mark size, inches wide; height in axes units below
for i, (lib, r, tis, arm, use) in enumerate(ROWS):
    y = top - (i + 0.5) * step
    if i in (3, 6, 9):
        c.plot([L, R_], [y + step / 2, y + step / 2], color=RULE, lw=0.6)
    col = NAVY if arm == "P" and use == "cmp" else ORANGE if arm == "R" else INKK
    c.text(colx["lib"], y, lib, ha="left", va="center", fontsize=11.5, color=col, fontweight="bold")
    c.text(colx["rat"], y, str(r), ha="center", va="center", fontsize=11.5, color=INKK)
    c.text(colx["tis"], y, tis, ha="left", va="center", fontsize=11, color=INKK)
    g = graft[lib]
    c.text(colx["hum"], y, f"{g:.0f} %" if g >= 1 else f"{g:.1f} %", ha="right", va="center", fontsize=11, color=INKK)
    for k in ["rna"] + (["dna"] if lib in array_of else []) + [use]:
        c.add_patch(Rectangle((colx[k] - S / 2, y - 0.022), S, 0.044, fc=INKK, ec=INKK, lw=0.8))
c.plot([L, R_], [top - 10 * step, top - 10 * step], color=INKK, lw=0.9)
renamed = [(lib, array_of[lib]) for lib, *_ in ROWS if lib in array_of and array_of[lib] != lib]
note = ("Array names " + ", ".join(arr for _, arr in renamed) + " = libraries " + ", ".join(lib for lib, _ in renamed)
        + ".\n" if renamed else "")
c.text(L, top - 10 * step - 0.03, note + "Not profiled: rat 65 (recurrent), culture C1.", ha="left", va="top",
       fontsize=10, color=GREY, linespacing=1.3)
save(fig, "fig_design.png", pad=0.04)
