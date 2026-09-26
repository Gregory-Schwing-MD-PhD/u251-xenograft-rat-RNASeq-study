# -*- coding: utf-8 -*-
"""Figures for the CNS 2026 Abstract 418 deck on the official CNS speaker template (04_build_cns_template_deck.py).

Same data as 02_make_figures.py (Supplementary_Data.xlsx sheets, figures/facts.json, the digitised PCA), redrawn in
the CNS template palette (navy 191D63, blue 4C70B7, red EF4627, gold E0A800) at the size each chart occupies on the
slide, with 14-20 pt type so it reads from the back of the room. The two horizontal bar charts of the first deck
(leading gene sets, drug ranking) are drawn as vertical columns. New here: the margin diagram and its known / not
known callouts, the theme and convergence charts (facts.json#themes), the three tables as images, the overview
cartoon cropped from ASSETS/Visual_abstract.png, and the ablation video of Nagaraja 2021 (64 x 64 px, 7.9 s, from
'Ablation video_TNN_18Sept26.pptx') re-encoded at 512 px with its first frame as the poster. Writes figures_cns/.

No number is typed here: every value is read from facts.json, the supplementary sheets or the deck's own tables.

    python SLIDES/04_make_cns_figures.py
"""
import importlib.util
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Ellipse, FancyBboxPatch, Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
from PIL import Image  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "figures_cns"
OUT.mkdir(exist_ok=True)
spec = importlib.util.spec_from_file_location("mk02", HERE / "02_make_figures.py")
mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)      # data loaders + label maps, no drawing
F, sheet, spread, LABELS, NICE_SUBTYPE, digitize_pca = mk.F, mk.sheet, mk.spread, mk.LABELS, mk.NICE_SUBTYPE, mk.digitize_pca

NAVY, BLUE, PALE, RED, GOLD, GREY, INK, LIGHT, FAINT, GRID = "#191D63", "#4C70B7", "#C9D4EC", "#EF4627", "#E0A800", "#8A8A8A", "#191D63", "#EEF1F8", "#C9CCCB", "#DDDDDD"
plt.rcParams.update({"font.family": "Arial", "font.size": 15, "axes.edgecolor": "#BFBFBF", "axes.linewidth": 0.9,
                     "xtick.color": GREY, "ytick.color": GREY, "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                     "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "mathtext.default": "regular",
                     "figure.dpi": 220, "savefig.dpi": 220})
MINUS = lambda s: s.replace("-", "−")  # noqa: E731


def save(fig, name, pad=0.05):
    fig.savefig(OUT / name, facecolor="white", bbox_inches="tight", pad_inches=pad)
    plt.close(fig)
    with Image.open(OUT / name) as im:
        print(f"  {name:28s} {im.size}  {im.size[0] / im.size[1]:.2f}")


def bare(ax, grid="both"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=14)


import matplotlib.text as mtext  # noqa: E402

DIRS = [(1, 0, "left", "center"), (-1, 0, "right", "center"), (0, 1, "center", "bottom"), (0, -1, "center", "top"),
        (0.8, 0.8, "left", "bottom"), (-0.8, 0.8, "right", "bottom"), (0.8, -0.8, "left", "top"), (-0.8, -0.8, "right", "top")]


def place_labels(ax, items, obstacles=(), fs=12, radii=(9, 16, 26, 38, 52), leader_from=16, halo=True, dirs=None, skip=False, **kw):
    """Greedy label placement. items: (x, y, text, extra-kwargs) in data coordinates, placed in the order given; an item's
    extra-kwargs may carry its own 'dirs' order. For each, try offsets around its point at increasing radii; keep the
    first whose text box stays inside the axes and overlaps neither a label already placed nor any obstacle point (x, y,
    radius in points). A leader line is drawn from leader_from points on. With skip=True a label with no clean slot is
    left off (and named on stdout) rather than drawn over something. Returns the annotations."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pt = fig.dpi / 72.0
    axbb = ax.get_window_extent(rend)
    obs = [(ax.transData.transform((ox, oy)), orad * pt) for ox, oy, orad in obstacles]
    placed, out = [], []
    for it in items:
        x, y, text = it[:3]
        extra = dict(kw, **(it[3] if len(it) > 3 else {}))
        fsz = extra.pop("fontsize", fs)
        own_dirs = extra.pop("dirs", None)
        best = None
        for rad in radii:
            for dx, dy, ha, va in (own_dirs or dirs or DIRS):
                ap = dict(arrowstyle="-", color=extra.get("color", "#8A8A8A"), lw=0.8, shrinkA=1, shrinkB=4, alpha=0.8) if rad >= leader_from else None
                t = ax.annotate(text, (x, y), xytext=(dx * rad, dy * rad), textcoords="offset points", ha=ha, va=va, fontsize=fsz,
                                arrowprops=ap, zorder=7, **extra)
                if halo:
                    t.set_path_effects([pe.withStroke(linewidth=3.2, foreground="white")])
                t.update_positions(rend)                      # resolve the offset before measuring the text box
                bb = mtext.Text.get_window_extent(t, rend)
                grow = bb.expanded(1.04, 1.12)
                inside = grow.x0 >= axbb.x0 - 2 and grow.x1 <= axbb.x1 + 2 and grow.y0 >= axbb.y0 - 2 and grow.y1 <= axbb.y1 + 2
                clash = any(grow.overlaps(p) for p in placed)
                hit = any((grow.x0 - r_ <= o[0] <= grow.x1 + r_) and (grow.y0 - r_ <= o[1] <= grow.y1 + r_) for o, r_ in obs)
                if inside and not clash and not hit:
                    best = (t, grow)
                    break
                t.remove()
            if best:
                break
        if best is None and skip:
            print(f"    place_labels: no clean slot for {text!r}; left off")
            continue
        if best is None:           # nothing clean: take the nearest right-hand slot and say so
            print(f"    place_labels: no clean slot for {text!r}")
            t = ax.annotate(text, (x, y), xytext=(radii[0], 0), textcoords="offset points", ha="left", va="center", fontsize=fsz, zorder=7, **extra)
            if halo:
                t.set_path_effects([pe.withStroke(linewidth=3.2, foreground="white")])
            t.update_positions(rend)
            best = (t, mtext.Text.get_window_extent(t, rend))
        placed.append(best[1]); out.append(best[0])
    return out


# ====================================================================== the margin diagram
def rings(ax, cx=0.0, cy=0.0, s=1.0, labels=True, fs=13):
    """Three zones; the fiber tip sits just above the centre of the core so no label lies under it."""
    ax.add_patch(Ellipse((cx, cy), 5.4 * s, 4.2 * s, fc=LIGHT, ec=BLUE, lw=2.0, zorder=1))
    ax.add_patch(Ellipse((cx, cy), 3.8 * s, 2.9 * s, fc="#FFF1C2", ec=GOLD, lw=2.4, zorder=2))
    ax.add_patch(Ellipse((cx, cy), 1.9 * s, 1.3 * s, fc="#FBD2C8", ec=RED, lw=2.4, zorder=3))
    ax.plot([cx - 2.95 * s, cx - 0.1 * s], [cy + 1.75 * s, cy + 0.25 * s], color=INK, lw=3.5, zorder=4, solid_capstyle="round")
    if labels:
        ax.text(cx + 0.05 * s, cy - 0.2 * s, "ablation\ncore", ha="center", va="center", fontsize=fs, color=RED, fontweight="bold", zorder=5,
                linespacing=1.05)
        ax.text(cx, cy - 1.05 * s, "sublethal margin", ha="center", va="center", fontsize=fs - 1, color="#9A7300", fontweight="bold", zorder=5)
        ax.text(cx, cy - 1.78 * s, "peritumoral brain", ha="center", va="center", fontsize=fs - 1, color=BLUE, fontweight="bold", zorder=5)
        ax.text(cx - 2.95 * s, cy + 1.9 * s, "laser fiber", ha="left", va="bottom", fontsize=fs, color=INK, zorder=5)


def fig_margin():
    fig = plt.figure(figsize=(6.4, 4.9))
    ax = fig.add_axes([0.0, 0.0, 0.58, 1.0]); ax.set_xlim(-3.0, 2.8); ax.set_ylim(-2.3, 2.4); ax.set_aspect("equal"); ax.set_axis_off()
    rings(ax, fs=13)
    # thermal dose against distance from the fiber, the ablative threshold and the sublethal band
    bx = fig.add_axes([0.66, 0.16, 0.32, 0.72])
    d = np.linspace(0, 1, 200); dose = np.exp(-3.2 * d)
    bx.plot(d, dose, color=INK, lw=3)
    thr = float(np.exp(-3.2 * 0.33)); sub = float(np.exp(-3.2 * 0.68))
    bx.axhline(thr, color=RED, lw=1.6, ls=(0, (4, 3)))
    bx.axhspan(sub, thr, color="#FFF1C2", zorder=0)
    bx.axvspan(0, 0.33, color="#FBD2C8", zorder=0, alpha=0.7)
    bx.text(0.98, thr + 0.03, "ablative\nthreshold", ha="right", va="bottom", fontsize=13, color=RED, linespacing=1.1)
    bx.text(0.98, (thr + sub) / 2, "sublethal", ha="right", va="center", fontsize=13, color="#9A7300", fontweight="bold")
    bx.set_xticks([]); bx.set_yticks([])
    for sd in ("top", "right"):
        bx.spines[sd].set_visible(False)
    bx.set_xlabel("distance from fiber", fontsize=14, color=INK); bx.set_ylabel("thermal dose", fontsize=14, color=INK)
    save(fig, "fig_margin.png", pad=0.02)


def fig_margin_known():
    s = 0.95
    fig = plt.figure(figsize=(11.6, 4.3))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.set_xlim(-7.6, 7.6); ax.set_ylim(-2.6, 2.6); ax.set_aspect("equal"); ax.set_axis_off()
    rings(ax, s=s, labels=False)
    ax.text(0.05 * s, -0.25 * s, "core", ha="center", va="center", fontsize=14, color=RED, fontweight="bold", zorder=5)
    # a few surviving cells in the margin ring
    rng = np.random.default_rng(3)
    for _ in range(70):
        t = rng.uniform(0, 2 * np.pi); r = rng.uniform(0.62, 0.92)
        ax.plot(r * 1.9 * s * np.cos(t), r * 1.45 * s * np.sin(t), "o", ms=5, color="#9A7300", zorder=4, alpha=0.8)
    # described (navy, left) and not yet done (red, right)
    def call(x, y, text, col, px, py, ha):
        ax.annotate(text, (px, py), xytext=(x, y), fontsize=13.5, color=col, ha=ha, va="center", linespacing=1.22,
                    fontweight="bold" if col == RED else "normal",
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.5, shrinkA=0, shrinkB=2), zorder=6)
    ax.text(-7.4, 2.5, "described", fontsize=16, color=NAVY, fontweight="bold", ha="left", va="top")
    call(-7.4, 1.45, "peri-ablation vessels open the\nblood–brain barrier for a time\n(Cleary 2026)", NAVY, -2.3, 0.95, "left")
    call(-7.4, 0.1, "the immune microenvironment of\nthe sublethal zone, profiled in a\nrodent model (Tao 2026)", NAVY, -1.55, -0.55, "left")
    call(-7.4, -1.55, "a first-pass analysis of these same\nlibraries, recurrent vs primary:\ncell cycle, motility, inflammation\n(Nagaraja 2026)", NAVY, -0.75, -0.45, "left")
    ax.text(7.4, 2.5, "not yet done", fontsize=16, color=RED, fontweight="bold", ha="right", va="top")
    call(7.4, 0.55, "the regrowing tumor's own\nprogram, sorted from host\nreads, read gene-set-wide,\nand whether it points\nto a drug", RED, 1.75, 0.55, "right")
    ax.text(7.4, -1.45, "this analysis", fontsize=15, color=RED, ha="right", va="center", style="italic")
    save(fig, "fig_margin_known.png", pad=0.02)


# ====================================================================== PCA
def fig_pca():
    g = digitize_pca()
    fig, ax = plt.subplots(figsize=(5.0, 4.9))
    bare(ax)
    ax.axhline(0, color=GRID, lw=0.9, zorder=1); ax.axvline(0, color=GRID, lw=0.9, zorder=1)
    style = {"Primary": (NAVY, "^", "Primary\n(pre-LITT)"), "Recurrent": (BLUE, "s", "Recurrent\n(post-LITT)")}   # as on the sorting and trajectory charts
    for name, pts in g.items():
        col, mark, _ = style[name]
        ax.scatter(pts[:, 0], pts[:, 1], s=300, marker=mark, color=col, edgecolor="white", linewidth=1.8, zorder=4)
    ax.text(-7.4, 4.4, style["Primary"][2], color=NAVY, fontsize=16, fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.text(7.0, 4.6, style["Recurrent"][2], color=BLUE, fontsize=16, fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.set_xlabel("PC1  (%.1f %% of variance)" % F["pc1"], fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("PC2  (%.1f %%)" % F["pc2"], fontsize=16, color=INK, labelpad=4)
    ax.set_xlim(-11.5, 11.5); ax.set_ylim(-10.5, 10.5)
    ax.set_xticks([-10, -5, 0, 5, 10]); ax.set_yticks([-10, -5, 0, 5, 10])      # v4: even ticks (2.5 steps printed as 8, 5, 2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    save(fig, "chart_pca.png")


# ====================================================================== volcano
def fig_volcano():
    d = sheet("S1_DE_all").dropna(subset=["padj", "log2FoldChange"])
    d["y"] = -np.log10(d.padj.clip(lower=1e-16))
    sig = d[(d.padj < 0.05) & (d.log2FoldChange.abs() > 1)]
    ns = d.drop(sig.index)
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    bare(ax)
    ax.scatter(ns.log2FoldChange, ns.y, s=8, color=FAINT, lw=0, zorder=2)
    ax.scatter(sig.log2FoldChange, sig.y, s=70, color=np.where(sig.log2FoldChange > 0, "#E08214", NAVY), lw=0.8, edgecolor="white", zorder=4)   # v4: up orange, down navy
    for v in (-1, 1):
        ax.axvline(v, color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.axhline(-np.log10(0.05), color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.set_xlim(-11.8, 11.8); ax.set_ylim(-0.6, sig.y.max() * 1.18)
    ax.set_xlabel("log$_2$ fold change   (recurrent vs primary)", fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("−log$_{10}$ adjusted P", fontsize=16, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    top = ax.get_ylim()[1]
    ax.text(-11.2, top * 0.955, "%d down" % F["de_down"], color=NAVY, fontsize=19, fontweight="bold", ha="left", va="top")
    ax.text(11.2, top * 0.955, "%d up" % F["de_up"], color="#E08214", fontsize=19, fontweight="bold", ha="right", va="top")
    ax.text(-11.2, -np.log10(0.05) + top * 0.02, "FDR 0.05", color=GREY, fontsize=13, ha="left", va="bottom")
    # the six largest movers (the table beside the chart) and the network hubs, largest change placed first
    big6 = sig.reindex(sig.log2FoldChange.abs().sort_values(ascending=False).index).head(6)
    hubs = sig[sig.symbol.isin(["BGN", "COL1A1", "IGFBP3", "CALB1", "GPC3"])]
    lab = pd.concat([big6, hubs]).drop_duplicates("symbol")
    lab = lab.reindex(lab.log2FoldChange.abs().sort_values(ascending=False).index)
    obstacles = [(r.log2FoldChange, r.y, 5.5) for _, r in sig.iterrows()]
    place_labels(ax, [(r.log2FoldChange, r.y, str(r.symbol)) for _, r in lab.iterrows()], obstacles=obstacles, fs=12.5,
                 color=INK, style="italic", radii=(11, 17, 26, 38, 52), leader_from=11)
    save(fig, "chart_volcano.png")


# ====================================================================== leading sets, vertical
SHORT = {"KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION": "translation\ninitiation",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION": "translation\nelongation",
         "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "GCN2\namino-acid\nstress",
         "KEGG_RIBOSOME": "ribosome", "REACTOME_SELENOAMINO_ACID_METABOLISM": "seleno-\namino-acid\nmetabolism",
         "REACTOME_CELLULAR_RESPONSE_TO_STARVATION": "starvation\nresponse"}


def hanging(ax, xs, vals, cols, labels, fs_lab=14, q=None, ylim=(-2.45, 0.35), yl="normalized enrichment score", fs_val=14):
    ax.bar(xs, vals, width=0.62, color=cols, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    for x, v in zip(xs, vals):
        ax.text(x, v - 0.05, MINUS(f"{v:.2f}"), ha="center", va="top", fontsize=fs_val, fontweight="bold", color=INK)
    if q is not None:
        for x, v, qq in zip(xs, vals, q):
            ax.text(x, v - 0.30, (f"q {qq:.3f}" if qq < 0.05 else f"q {qq:.2f}"), ha="center", va="top", fontsize=12.5,
                    color=INK if qq < 0.05 else GREY, fontweight="bold" if qq < 0.05 else "normal")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=fs_lab, color=INK, linespacing=1.1)
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=6)
    ax.set_ylim(*ylim); ax.set_ylabel(yl, fontsize=15, color=INK)
    ax.set_yticks([0, -0.5, -1.0, -1.5, -2.0]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.1f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)


def fig_gsea_sets():
    sets = [s for s in F["gsea_down_top"] if s["name"] in LABELS][:6]
    assert len(sets) == 6
    pmax = max(s["p"] for s in sets)
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    hanging(ax, list(range(6)), [s["nes"] for s in sets], [NAVY if s["q"] < 0.05 else BLUE for s in sets],
            [SHORT[s["name"]] for s in sets], q=[s["q"] for s in sets], fs_lab=11.5, ylim=(-2.9, 0.35))
    ax.text(-0.5, -2.86, "nominal p ≤ %.3f for all six  ·  %s gene sets tested" % (pmax, format(F["gsea_n_sets"], ",")),
            fontsize=12.5, color=GREY, style="italic", ha="left", va="bottom")
    save(fig, "chart_gsea_sets.png")
    json.dump({"p_max": pmax, "p": {s["name"]: s["p"] for s in sets}}, open(OUT / "chart_gsea_sets.json", "w"), indent=1)


def theme(name):
    th = F["themes"]
    if name in th:
        return th[name]
    hits = [k for k in th if name.lower().replace(" ", "") in k.lower().replace(" ", "")]
    assert len(hits) == 1, (name, hits, list(th))
    return th[hits[0]]


S3 = None


def gs(name):
    """One gene set from the pipeline's full GSEA table (NES, nominal p, FDR q); F["themes"] names are accepted too."""
    global S3
    if S3 is None:
        S3 = sheet("S3_GSEA_pipeline").set_index("gene_set")
    if name in S3.index:
        r = S3.loc[name]
        return {"nes": float(r.NES), "p": float(r.nominal_p), "q": float(r.FDR_q)}
    return theme(name)


# growth signalling, glycolysis, respiration, the TCA cycle and sterol synthesis, each in the collections that carry it:
# glycolysis falls in all three; respiration does not agree across collections; sterol synthesis rises
THEMES_ENERGY = [("growth signalling", NAVY, [("ribosome\nbiogenesis", "GOBP_RIBOSOME_BIOGENESIS"), ("MYC", "HALLMARK_MYC_TARGETS_V1"),
                                              ("mTORC1", "HALLMARK_MTORC1_SIGNALING")]),
                 ("glycolysis", BLUE, [("Hallmark", "HALLMARK_GLYCOLYSIS"), ("Reactome", "REACTOME_GLYCOLYSIS"),
                                       ("KEGG", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")]),
                 ("respiration", GREY, [("Hallmark\nOXPHOS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"), ("KEGG\nOXPHOS", "KEGG_OXIDATIVE_PHOSPHORYLATION"),
                                        ("Reactome\nETC", "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT")]),
                 ("TCA", GREY, [("TCA\ncycle", "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE")]),
                 ("sterols", RED, [("cholesterol\nsynthesis", "REACTOME_CHOLESTEROL_BIOSYNTHESIS")])]


def grouped_bars(groups, name, figsize, pitch=1.9, gapg=1.2, width=1.0, fs_tick=11.5, ylim=(-3.2, 2.8), head_y=None, lookup=None,
                 up_note=True):
    """Bars that may go either way, grouped, with NES and nominal p on each bar and the group name under the axis."""
    lookup = lookup or theme
    fig, ax = plt.subplots(figsize=figsize)
    xs, vals, cols, labels, ps, heads, rows = [], [], [], [], [], [], []
    x = 0.0
    for head, col, members in groups:
        x0 = x
        for lab, key in members:
            t = lookup(key); xs.append(x); vals.append(t["nes"]); cols.append(col); labels.append(lab); ps.append(t["p"]); x += pitch
            rows.append({"group": head, "label": lab.replace("\n", " "), "set": key, "nes": t["nes"], "p": t["p"], "q": t.get("q")})
        heads.append(((x0 + x - pitch) / 2, head, col)); x += gapg
    ax.bar(xs, vals, width=width, color=cols, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    for xx, v, p in zip(xs, vals, ps):
        up = v > 0
        ax.text(xx, v + (0.05 if up else -0.05), MINUS(f"{v:+.2f}"), ha="center", va="bottom" if up else "top", fontsize=13, fontweight="bold", color=INK)
        ax.text(xx, v + (0.32 if up else -0.32), f"p {p:.3f}" if p < 0.01 else f"p {p:.2f}", ha="center", va="bottom" if up else "top", fontsize=11.5, color=GREY)
    hy = head_y if head_y is not None else ylim[0] + 0.35
    for xc, head, col in heads:
        ax.text(xc, hy, head, ha="center", va="top", fontsize=15, color=col, fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=fs_tick, color=INK, linespacing=1.1)
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=4)
    ax.set_ylim(*ylim); ax.set_ylabel("normalized enrichment score", fontsize=14, color=INK)
    ax.set_yticks([v for v in (2, 1, 0, -1, -2) if ylim[0] < v < ylim[1]]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    if up_note:
        ax.text(-0.7, ylim[1] - 0.1, "up in recurrence", ha="left", va="top", fontsize=13, color=GREY, style="italic")
        ax.text(xs[-1] + 0.6, -0.15, "down in recurrence", ha="right", va="top", fontsize=13, color=GREY, style="italic")
    save(fig, name)
    json.dump(rows, open(OUT / name.replace(".png", ".json"), "w"), indent=1)
    return rows


def fig_themes():
    grouped_bars(THEMES_ENERGY, "chart_themes.png", figsize=(13.4, 4.6), pitch=2.85, gapg=1.1, width=1.35, fs_tick=12,
                 ylim=(-3.2, 2.6), lookup=gs, up_note=False)


# ====================================================================== landscape
def fig_landscape():
    g = sheet("S3_GSEA_pipeline").copy()
    g["y"] = -np.log10(g.FDR_q.clip(lower=1e-3))
    hit = g[g.FDR_q < 0.25].sort_values("NES"); rest = g.drop(hit.index)
    fig, ax = plt.subplots(figsize=(6.9, 4.6))
    bare(ax)
    ax.scatter(rest.NES, rest.y, s=7, color=FAINT, lw=0, zorder=2)
    ax.scatter(hit.NES, hit.y, s=110, color=NAVY, edgecolor="white", lw=1.0, zorder=5)
    ax.axhline(-np.log10(0.25), color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.axvline(0, color=GRID, lw=1.0, zorder=1)
    ax.set_xlim(-2.55, 2.55); ax.set_ylim(-0.06, 2.42); ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    ax.set_xlabel("normalized enrichment score", fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("−log$_{10}$ FDR q", fontsize=16, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.text(0.3, -np.log10(0.25) - 0.03, "q = 0.25", color=GREY, fontsize=13, ha="center", va="top", zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none"))
    short = dict(LABELS, **{"REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "GCN2 amino-acid stress"})
    short = {k: v.replace("\n", " ") for k, v in short.items()}
    # one label column to the right of the hits, in the order of their heights, so no leader crosses a label
    hit = hit.sort_values("y")
    ys = spread(hit.y.tolist(), 0.74, 1.86, 0.215)
    col_x = hit.NES.max() + 0.42
    for (_, row), y in zip(hit.iterrows(), ys):
        ax.annotate(short[row.gene_set], (row.NES, row.y), xytext=(col_x, y), ha="left", va="center", fontsize=13.5, color=INK, zorder=6,
                    arrowprops=dict(arrowstyle="-", color="#AFB6B4", lw=0.8, shrinkA=2, shrinkB=4, relpos=(0.0, 0.5)))
    ax.annotate("nothing on this side\nclears q = %.2f" % F["gsea_up_best_q"], xy=(1.56, 0.06), xytext=(1.72, 1.2), ha="center", va="center",
                fontsize=14, color=GREY, linespacing=1.35, arrowprops=dict(arrowstyle="->", color=GREY, lw=1.1, connectionstyle="arc3,rad=0.22"))
    ax.text(-2.50, 2.40, "down in recurrence", color=NAVY, fontsize=15, fontweight="bold", ha="left", va="top")
    ax.text(2.50, 2.40, "up in recurrence", color=GREY, fontsize=15, ha="right", va="top")
    save(fig, "chart_gsea_landscape.png")


# ====================================================================== subtypes
def fig_subtypes():
    rows = sorted(F["gsva"], key=lambda r: r["change"])
    q = bh([r["p"] for r in rows])                        # facts.json carries a mis-ordered BH; recomputed here from the p values
    for r, qq in zip(rows, q):
        r["q_bh"] = qq
    fig, ax = plt.subplots(figsize=(6.6, 4.7))
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.plot([-0.13, 1.0], [0, 0], color=GRID, lw=1.0, zorder=1)
    ends = []
    for r in rows:
        on = r["p"] < 0.05; mtc = False   # v4: only the nominally significant signature is highlighted
        col = NAVY if on else GOLD if mtc else "#C4CBC9"
        a, b = r["primary"], r["recurrent"]
        ax.plot([0, 1], [a, b], color=col, lw=3.2 if (on or mtc) else 1.5, zorder=4 if (on or mtc) else 2, solid_capstyle="round")
        ax.scatter([0, 1], [a, b], s=100 if (on or mtc) else 40, color=col, edgecolor="white", lw=1.3, zorder=5 if (on or mtc) else 3)
        ends.append((b, r, on, mtc, col))
    ends.sort(key=lambda e: e[0])
    for (b, r, on, mtc, col), y in zip(ends, spread([e[0] for e in ends], -0.62, 0.62, 0.105)):
        lab = NICE_SUBTYPE[r["name"]] + ("  *" if on else "")
        ax.plot([1.02, 1.09, 1.13], [b, y, y], color=col if (on or mtc) else "#C4CBC9", lw=0.9, zorder=2)
        ax.text(1.15, y, lab, fontsize=14 if (on or mtc) else 12.5, color=NAVY if on else "#9A7300" if mtc else GREY,
                fontweight="bold" if (on or mtc) else "normal", va="center")
    ax.set_xlim(-0.13, 2.15); ax.set_ylim(-0.66, 0.66)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Primary\n(pre-LITT)", "Recurrent\n(post-LITT)"], fontsize=15, color=INK, linespacing=1.3)
    ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5]); ax.set_yticklabels([MINUS("-0.5"), MINUS("-0.25"), "0", "0.25", "0.5"], fontsize=13)
    ax.tick_params(length=0)
    ax.set_ylabel("GSVA enrichment score", fontsize=15, color=INK)
    best = min(rows, key=lambda r: r["p"])
    ax.text(-0.13, 0.63, "*  nominal p < 0.05; none clears FDR across the ten (lowest q = %.3f)" % best["q_bh"], fontsize=12, color=GREY,
            style="italic", va="center")
    save(fig, "chart_subtypes.png")
    json.dump({r["name"]: {"p": r["p"], "q_bh": r["q_bh"], "change": r["change"]} for r in rows}, open(OUT / "chart_subtypes.json", "w"), indent=1)


def bh(p):
    """Benjamini-Hochberg in p order, monotone from the top."""
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); prev = 1.0
    for k in range(n - 1, -1, -1):
        prev = min(prev, p[o[k]] * n / (k + 1)); q[o[k]] = prev
    return q.tolist()


# ====================================================================== PPI
def fig_ppi():
    e = sheet("S6_PPI_edges").drop_duplicates(subset=["gene1", "gene2"])
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.gene1, r.gene2, w=float(r.STRING_combined_score))
    deg = dict(G.degree())
    hubs = {n for n, d in deg.items() if d >= 3}                   # a cutoff, so ties cannot fall either way
    json.dump({"degree": deg, "hubs_degree_ge3": sorted(hubs)}, open(OUT / "chart_ppi.json", "w"), indent=1)
    lfc = F["ppi_nodes"]
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    pos = {}
    big = G.subgraph(comps[0])
    p = nx.spring_layout(big, seed=7, k=1.6, iterations=900, weight=None)
    pts = np.array([p[n] for n in big]); pts = pts - pts.mean(0); pts = pts / np.abs(pts).max() * 1.6
    for n, xy in zip(big, pts):
        pos[n] = xy + np.array([0.0, 0.95])
    # the four-node component is a star on RYR2, placed by hand under the large one so no label crosses a node
    small = comps[1]
    centre = max(small, key=lambda n: deg[n])
    others = sorted(n for n in small if n != centre)
    star = {centre: (0.15, -1.85)}
    for n, xy in zip(others, [(-0.8, -2.15), (1.05, -1.45), (1.05, -2.35)]):
        star[n] = xy
    for n, xy in star.items():
        pos[n] = np.array(xy)
    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    ax.set_axis_off(); ax.set_aspect("equal")
    scores = np.array([d["w"] for *_, d in G.edges(data=True)])
    for u, v, d in G.edges(data=True):
        lw = 1.6 + 4.0 * (d["w"] - scores.min()) / (scores.max() - scores.min())
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#B7C3DE", lw=lw, zorder=1, solid_capstyle="round")
    for n in G.nodes():
        hub = n in hubs; col = RED if lfc.get(n, 0) > 0 else NAVY
        ax.scatter(*pos[n], s=470 if hub else 250, color=col if hub else "white", edgecolor=col, linewidth=0 if hub else 2.8, zorder=3)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 1.05, max(xs) + 1.05); ax.set_ylim(min(ys) - 0.62, max(ys) + 0.55)   # the legend is slide text
    # labels last, clear of every node; hubs first so they get the nearest slots
    obstacles = [(pos[n][0], pos[n][1], 12.5 if n in hubs else 9.5) for n in G.nodes()]
    order = sorted(G.nodes(), key=lambda n: (-deg[n], n))
    cen = {}
    for comp in comps:
        c = np.mean([pos[m] for m in comp], axis=0)
        for m in comp:
            cen[m] = c

    def outward(n):                 # try the directions that point away from the node's own cluster first
        v = pos[n] - cen[n]
        v = v / (np.linalg.norm(v) + 1e-9)
        return sorted(DIRS, key=lambda d: -(d[0] * v[0] + d[1] * v[1]) / np.hypot(d[0], d[1]))
    place_labels(ax, [(pos[n][0], pos[n][1], n, {"fontsize": 14 if n in hubs else 13, "fontweight": "bold", "dirs": outward(n),
                                                   "color": RED if lfc.get(n, 0) > 0 else NAVY}) for n in order],
                 obstacles=obstacles, radii=(11, 17, 26, 36), leader_from=26)
    save(fig, "chart_ppi.png")


# ====================================================================== drugs, vertical
def overall_ranks():
    """Rank of every clinically available compound by the BBB-weighted score, as in ESM_1 S13 (checked against it)."""
    x = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx")
    s12 = x.parse("S12_Drug_ranking_full")
    s12 = s12[s12["Reached a clinic"].astype(str).str.lower() == "yes"]
    s12 = s12.sort_values("Integrated score (|NES|^1.5 x ADMET-AI BBB)", ascending=False).drop_duplicates("Drug")
    rank = {str(d).lower(): i + 1 for i, d in enumerate(s12["Drug"])}
    s13 = x.parse("S13_Drug_candidates_top20")
    for _, r in s13.iterrows():
        assert rank[str(r.Drug).lower()] == int(r["Rank (BBB-weighted)"]), (r.Drug, rank[str(r.Drug).lower()])
    return rank, s13


def fig_drugs():
    tier = F["drug_tierA"][:10]
    rank, _ = overall_ranks()
    fig, ax = plt.subplots(figsize=(5.9, 4.6))
    xs = list(range(len(tier)))
    cols = [NAVY if d["drug"] == "ciclopirox" else BLUE for d in tier]
    ax.bar(xs, [d["score"] for d in tier], width=0.66, color=cols, zorder=3)
    for x, d in zip(xs, tier):
        ax.text(x, d["score"] + 0.05, "%.2f" % d["score"], ha="center", va="bottom", fontsize=13.5,
                fontweight="bold" if d["drug"] == "ciclopirox" else "normal", color=INK)
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{d['drug'].replace(' bromide', '')} ({rank[d['drug'].lower()]})" for d in tier], fontsize=13, color=INK, rotation=35, ha="right")
    bare(ax, grid="y")
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, max(d["score"] for d in tier) * 1.22)
    ax.set_ylabel("|NES|$^{1.5}$ × predicted BBB permeability", fontsize=15, color=INK)
    ax.text(len(tier) - 0.5, max(d["score"] for d in tier) * 1.2,
            "top 10 of the %d compounds that clear\nboth barrier models; overall rank in brackets" % len(F["drug_tierA"]),
            ha="right", va="top", fontsize=12.5, color=GREY, style="italic", linespacing=1.25)
    save(fig, "chart_drugs.png")


# ====================================================================== conclusions: the arms move together
ARMS = [("translation", NAVY, [("initiation", "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"), ("elongation", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"),
                               ("ribosome", "KEGG_RIBOSOME")]),
        ("growth", NAVY, [("rRNA\nprocessing", "rRNA processing"), ("ribosome\nbiogenesis", "ribosome biogenesis"), ("mTORC1", "mTORC1 signalling"),
                                ("MYC", "MYC targets")]),
        ("glycolysis", BLUE, [("Hallmark", "HALLMARK_GLYCOLYSIS"), ("Reactome", "REACTOME_GLYCOLYSIS"), ("KEGG", "KEGG_GLYCOLYSIS_GLUCONEOGENESIS")]),
        ("hypoxia", BLUE, [("HIF1\ntargets", "HIF1 targets"), ("hypoxia\nmetagene", "hypoxia metagene"), ("hypoxia", "hypoxia")]),
        ("unchanged", GREY, [("TCA\ncycle", "TCA cycle")])]


def fig_arms():
    """One column per arm = the mean normalized enrichment score of its member sets; each member a dot on the column."""
    top = {s["name"]: s["nes"] for s in F["gsea_down_top"]}
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    xs, means, cols, labels, lows = [], [], [], [], []
    members_txt = []
    for i, (head, col, members) in enumerate(ARMS):
        vals = [top[key] if key in top else gs(key)["nes"] for _, key in members]
        x = i * 1.25
        xs.append(x); means.append(float(np.mean(vals))); cols.append(col); lows.append(min(vals))
        labels.append(head if head != "unchanged" else "TCA cycle\n(control)")
        offs = np.linspace(-0.18, 0.18, len(vals)) if len(vals) > 1 else [0.0]
        for v, o in zip(sorted(vals), offs):
            ax.scatter(x + o, v, s=70, color="white", edgecolor=col, linewidth=2.0, zorder=5)
        members_txt.append(f"{head}: " + ", ".join(f"{lab.replace(chr(10), ' ')} {v:.2f}" for (lab, _), v in zip(members, vals)))
    ax.bar(xs, means, width=0.72, color=cols, zorder=3, alpha=0.9)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    for x, m, lo, (head, _, members) in zip(xs, means, lows, ARMS):
        ax.text(x, lo - 0.09, MINUS(f"{m:.2f}"), ha="center", va="top", fontsize=15, fontweight="bold", color=INK)
        ax.text(x, lo - 0.33, f"mean of {len(members)}" if len(members) > 1 else "1 set", ha="center", va="top", fontsize=12.5, color=GREY)
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=14, color=INK, linespacing=1.1)
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=6)
    ax.set_ylim(-2.85, 0.25); ax.set_ylabel("enrichment score (NES)", fontsize=14, color=INK)
    ax.set_yticks([0, -0.5, -1.0, -1.5, -2.0]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.1f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.text(-0.55, -2.84, "columns, the mean of each arm  ·  dots, its member gene sets", fontsize=12, color=GREY, style="italic",
            ha="left", va="bottom")
    save(fig, "chart_arms.png")
    json.dump(members_txt, open(OUT / "chart_arms_members.json", "w"), indent=1)


# ====================================================================== tables as images
def table_image(name, header, rows, colw, fig_w, fs=16, hl=None, row_h=0.6, align=None):
    W = sum(colw); HEAD = 0.72; Hh = HEAD + row_h * len(rows) + 0.06
    align = align or (["left"] + ["center"] * (len(header) - 1))
    fig = plt.figure(figsize=(fig_w, fig_w * Hh / W))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(-Hh, 0); ax.set_axis_off()
    ax.add_patch(Rectangle((0, -HEAD), W, HEAD, facecolor=NAVY, edgecolor="none"))

    def cell(x, w, y, h, text, al, **kw):
        xx = x + 0.12 if al == "left" else x + w - 0.12 if al == "right" else x + w / 2
        ax.text(xx, y + h / 2, text, ha=al, va="center", **kw)
    x = 0
    for j, (h_, w_) in enumerate(zip(header, colw)):
        cell(x, w_, -HEAD, HEAD, h_, align[j], fontsize=fs, fontweight="bold", color="white"); x += w_
    y = -HEAD
    for i, row in enumerate(rows):
        y -= row_h
        if hl and i in hl:
            ax.add_patch(Rectangle((0, y), W, row_h, facecolor=PALE, edgecolor="none"))
        elif i % 2 == 0:
            ax.add_patch(Rectangle((0, y), W, row_h, facecolor="#F4F6FB", edgecolor="none"))
        x = 0
        for j, (c, w_) in enumerate(zip(row, colw)):
            col = RED if (j == 1 and str(c).startswith("+")) else INK
            cell(x, w_, y, row_h, c, align[j], fontsize=fs, color=col, fontweight="bold" if (hl and i in hl) else "normal",
                 style="italic" if (j == 0 and name != "table_drugs.png") else "normal"); x += w_
    save(fig, name, pad=0.0)


def tables():
    deck = json.load(open(HERE / "deck_2026-09-08_tables.json", encoding="utf-8"))   # the three tables of the 8 Sept deck
    fold, hubs, drugs = deck["fold"], deck["hubs"], deck["drugs"]
    for tb in (fold, hubs):                                   # Arial has no subscript-two glyph: mathtext instead
        tb[0] = [c.replace("log₂FC", "log$_2$FC") for c in tb[0]]
    fold[0][0], fold[0][2] = "gene", "locus type"
    table_image("table_fold.png", fold[0], fold[1:], [3.4, 1.25, 3.2], 4.6, fs=13, row_h=0.56)
    hubs[0][0] = "gene"
    table_image("table_hubs.png", hubs[0], hubs[1:], [1.6, 1.1, 3.6], 5.4, fs=15, row_h=0.46)
    # the six highest-ranked compounds on which both barrier models agree, with their OVERALL rank among the 54 (S13)
    _, s13 = overall_ranks()
    agree = s13[s13["Both BBB models agree"].astype(str).str.lower() == "yes"].head(6)
    rows = [[str(int(r["Rank (BBB-weighted)"])), str(r.Drug), str(int(r["Max clinical phase"])), f"{abs(r.NES):.2f}",
             f"{r['ADMET-AI BBB probability']:.2f}", f"{r['Integrated score (|NES|^1.5 x ADMET-AI BBB)']:.2f}"] for _, r in agree.iterrows()]
    table_image("table_drugs.png", ["rank", "compound", "phase", "|NES|", "BBB", "score"], rows, [0.75, 2.3, 0.9, 0.9, 0.9, 1.0], 5.0,
                fs=15, row_h=0.52, hl={0}, align=["center", "left", "center", "center", "center", "center"])
    json.dump(rows, open(OUT / "table_drugs.json", "w"), indent=1)
    s13j = [{"rank": int(r["Rank (BBB-weighted)"]), "rank_unweighted": int(r["Rank (unweighted)"]), "drug": str(r.Drug).lower(),
             "score": float(r["Integrated score (|NES|^1.5 x ADMET-AI BBB)"]), "nes": float(r.NES),
             "both_agree": str(r["Both BBB models agree"]).lower() == "yes", "egg": str(r["BOILED-Egg yolk"])} for _, r in s13.iterrows()]
    json.dump(s13j, open(OUT / "s13_top20.json", "w"), indent=1)


# ====================================================================== overview cartoon + video
def cartoon():
    im = Image.open(ROOT / "ASSETS" / "Visual_abstract.png").convert("RGB")
    w, h = im.size
    crop = im.crop((int(0.34 * w), int(0.36 * h), int(0.97 * w), int(0.96 * h)))     # the two in-vivo stages
    crop.save(OUT / "overview_cartoon.png")
    print("  overview_cartoon.png", crop.size)


def video():
    import imageio_ffmpeg
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    src = OUT / "ablation_64.mp4"
    with zipfile.ZipFile(ROOT / "Ablation video_TNN_18Sept26.pptx") as z:
        media = [n for n in z.namelist() if n.startswith("ppt/media/") and n.lower().endswith(".mp4")]
        assert len(media) == 1, media
        src.write_bytes(z.read(media[0]))
    up = OUT / "ablation_512.mp4"
    subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-y", "-i", str(src), "-vf", "scale=512:512:flags=lanczos,format=yuv420p",
                    "-c:v", "libx264", "-preset", "slow", "-crf", "18", "-movflags", "+faststart", "-an", str(up)], check=True)
    subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-y", "-i", str(up), "-vf", "select=eq(n\\,0)", "-vframes", "1",
                    str(OUT / "video_poster.png")], check=True)
    print("  ablation_512.mp4", up.stat().st_size // 1024, "kB; poster", Image.open(OUT / "video_poster.png").size)


# ====================================================================== the study in one strip (overview)
def fig_pipeline():
    """Five stages as boxes; the MRI-guided LITT stage carries the first frame of the ablation video."""
    fig = plt.figure(figsize=(11.5, 2.1))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 11.5); ax.set_ylim(0.95, 2.9); ax.set_axis_off()
    stages = [("U251N cells", "human\nglioblastoma line", NAVY), ("orthotopic tumor", "athymic RNU/RNU\nrat brain", NAVY),
              ("MRI-guided LITT", "", RED), ("recurrence", "regrowth after\nablation", GOLD),
              ("bulk RNA-seq", "3 primary + 3 recurrent;\nhuman and shared\nreads analysed", BLUE)]
    w, gap = 1.95, 0.36
    x = 0.12
    for i, (head, sub, col) in enumerate(stages):
        ax.add_patch(FancyBboxPatch((x, 1.05), w, 1.75, boxstyle="round,pad=0,rounding_size=0.12", fc=LIGHT, ec=col, lw=2.5))
        ax.add_patch(FancyBboxPatch((x, 2.25), w, 0.55, boxstyle="round,pad=0,rounding_size=0.12", fc=col, ec=col, lw=0))
        ax.text(x + w / 2, 2.52, head, ha="center", va="center", fontsize=15.5, fontweight="bold", color="white")
        if head == "MRI-guided LITT":
            im = Image.open(OUT / "video_poster.png").convert("L")
            ax.imshow(im, extent=(x + 0.1, x + 0.8, 1.2, 1.9), cmap="gray", zorder=3, aspect="auto")
            ax.text(x + 0.88, 1.55, "ablation\nwatched by\ndiffusion\nMRI", ha="left", va="center", fontsize=12, color=INK, linespacing=1.15)
        else:
            ax.text(x + w / 2, 1.62, sub, ha="center", va="center", fontsize=12.5, color=INK, linespacing=1.2)
        if i < len(stages) - 1:
            ax.annotate("", (x + w + gap - 0.05, 1.92), (x + w + 0.05, 1.92), arrowprops=dict(arrowstyle="-|>", color=GOLD, lw=3, mutation_scale=22))
        x += w + gap
    save(fig, "overview_pipeline.png", pad=0.02)


# ====================================================================== what else moves: hypoxia, iron, stress down; cell cycle up
THEMES2 = [("hypoxia", NAVY, [("HIF1\ntargets", "HIF1 targets"), ("Buffa\nmetagene", "hypoxia metagene"), ("Hallmark\nhypoxia", "hypoxia"),
                              ("angio-\ngenesis", "angiogenesis")]),
           ("iron, stress", BLUE, [("iron\nuptake", "iron uptake and transport"), ("senes-\ncence", "senescence"), ("EMT", "EMT")]),
           ("cell cycle", RED, [("spindle", "mitotic spindle"), ("G2M", "G2M checkpoint"), ("E2F", "E2F targets")])]


def fig_themes2():
    grouped_bars(THEMES2, "chart_themes2.png", figsize=(10.4, 4.8), pitch=2.1, gapg=1.2, width=1.05, fs_tick=11.5, ylim=(-3.2, 2.8))


# ====================================================================== the drug filter as a funnel
def fig_funnel():
    """93 compounds named by the 100 most opposing signatures -> 54 with a clinical phase and a structure the barrier
    model can score -> 15 that clear both barrier models. The prior-art audit covered the top 20 of the 54, not these 15,
    so it is not a step of this funnel."""
    x = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx")
    s12 = x.parse("S12_Drug_ranking_full")
    n_sig, n_comp = len(s12), s12["Drug"].astype(str).str.lower().nunique()      # 'Cobalt' and 'COBALT' are one compound
    assert (n_sig, n_comp) == (100, 92), (n_sig, n_comp)
    steps = [(n_comp, f"compounds named\nby the {n_sig} most\nopposing signatures"),
             (F["drug_clinical_n"], "with a clinical\nphase and a\nscorable structure"),
             (len(F["drug_tierA"]), "clear both\nbarrier\nmodels")]
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    cols = [BLUE, BLUE, NAVY]
    xs = [0, 1.55, 3.1]
    ax.bar(xs, [s[0] for s in steps], width=1.0, color=cols, zorder=3)
    for xx, (n, lab), col in zip(xs, steps, cols):
        ax.text(xx, n + 2, str(n), ha="center", va="bottom", fontsize=24, fontweight="bold", color=col)
        ax.text(xx, -4, lab, ha="center", va="top", fontsize=12.5, color=INK, linespacing=1.2)
    ax.set_xticks([]); ax.set_ylim(0, 110); ax.set_xlim(-0.75, 3.85)
    for sd in ("top", "right", "left"):
        ax.spines[sd].set_visible(False)
    ax.set_yticks([])
    ax.text(3.1, 36, "ciclopirox\nranks first", ha="center", va="bottom", fontsize=13.5, fontweight="bold", color=NAVY, linespacing=1.15)
    fig.subplots_adjust(bottom=0.25)
    save(fig, "chart_funnel.png")
    json.dump({"n_signatures": int(n_sig), "n_compounds": int(n_comp), "n_clinical": int(F["drug_clinical_n"]), "n_both_bbb": len(F["drug_tierA"])},
              open(OUT / "chart_funnel.json", "w"), indent=1)


# ====================================================================== DepMap: does U-251 depend on what ciclopirox hits?
def fig_depmap():
    rows = [r.split(",") for r in (ROOT / "depmap" / "depmap_u251_targets.csv").read_text().strip().splitlines()[1:]]
    want = [("RRM1", "ribonucleotide\nreductase"), ("RRM2", "ribonucleotide\nreductase"), ("DOHH", "eIF5A\nhypusination"),
            ("DHPS", "eIF5A\nhypusination"), ("EIF5A", "eIF5A\nhypusination"), ("VHL", "PHD–HIF\naxis"), ("HIF1A", "PHD–HIF\naxis"),
            ("FTH1", "iron\nhandling"), ("GPX4", "iron\nhandling"), ("TFRC", "iron\nhandling")]
    by = {r[1]: r for r in rows}
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    xs, vals, cols, labels, ce = [], [], [], [], []
    x = 0.0; last = None; heads = []
    for gene, grp in want:
        r = by[gene]; eff = float(r[2]); med = float(r[3])
        if last is not None and grp != last:
            heads.append(((x0 + x - 1.45) / 2, last)); x += 0.9
        if grp != last:
            x0 = x
        xs.append(x); vals.append(eff); labels.append(gene); ce.append(med)
        cols.append(NAVY if eff < -0.5 else "#C4CBC9")
        last = grp; x += 1.45
    heads.append(((x0 + x - 1.45) / 2, last))
    ax.bar(xs, vals, width=0.8, color=cols, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    ax.axhline(-0.5, color=RED, lw=1.6, ls=(0, (4, 3)), zorder=4)
    for xx, v, g, med in zip(xs, vals, labels, ce):         # the median across 1,178 lines as a black tick on each bar
        ax.plot([xx - 0.52, xx + 0.52], [med, med], color=INK, lw=2.6, zorder=5, solid_capstyle="butt")
        lo = min(v, med)
        ax.text(xx, lo - 0.06 if v < 0 else v + 0.05, MINUS(f"{v:.2f}"), ha="center", va="top" if v < 0 else "bottom", fontsize=12.5,
                fontweight="bold" if v < -0.5 else "normal", color=INK, zorder=6,
                bbox=dict(boxstyle="square,pad=0.08", fc="white", ec="none"))
    for xc, head in heads:
        ax.text(xc, 0.62, head, ha="center", va="bottom", fontsize=13, color=NAVY, fontweight="bold", linespacing=1.1)
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=12.5, color=INK, style="italic")
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=4)
    ax.set_ylim(-4.35, 1.35); ax.set_ylabel("gene effect in U-251 MG", fontsize=14, color=INK)
    ax.set_yticks([0, -1, -2, -3]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.text(xs[-1] + 0.5, -4.32, "bar, U-251 MG (navy where it depends, below −0.5)  ·  black tick, median of 1,178 lines\n"
                                  "dashed red, the dependency threshold (−0.5)",
            ha="right", va="bottom", fontsize=11.5, color=GREY, style="italic", linespacing=1.3)
    save(fig, "chart_depmap.png")


if __name__ == "__main__":
    print("figures ->", OUT)
    fig_margin(); fig_margin_known(); fig_pca(); fig_volcano(); fig_gsea_sets(); fig_themes(); fig_landscape(); fig_subtypes()
    fig_ppi(); fig_drugs(); fig_arms(); tables(); cartoon(); video(); fig_pipeline(); fig_themes2(); fig_funnel(); fig_depmap()
    print("done")
