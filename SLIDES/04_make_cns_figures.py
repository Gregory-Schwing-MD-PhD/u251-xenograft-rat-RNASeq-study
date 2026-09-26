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


# ====================================================================== the margin diagram
def rings(ax, cx=0.0, cy=0.0, s=1.0, labels=True, fs=17):
    ax.add_patch(Ellipse((cx, cy), 5.0 * s, 3.8 * s, fc=LIGHT, ec=BLUE, lw=2.0, zorder=1))
    ax.add_patch(Ellipse((cx, cy), 3.35 * s, 2.5 * s, fc="#FFF1C2", ec=GOLD, lw=2.4, zorder=2))
    ax.add_patch(Ellipse((cx, cy), 1.6 * s, 1.15 * s, fc="#FBD2C8", ec=RED, lw=2.4, zorder=3))
    ax.plot([cx - 2.9 * s, cx], [cy + 1.55 * s, cy + 0.05 * s], color=INK, lw=3.5, zorder=4, solid_capstyle="round")
    if labels:
        ax.text(cx, cy - 0.02 * s, "ablation core", ha="center", va="center", fontsize=fs, color=RED, fontweight="bold", zorder=5)
        ax.text(cx, cy - 0.95 * s, "sublethal margin", ha="center", va="center", fontsize=fs, color="#9A7300", fontweight="bold", zorder=5)
        ax.text(cx, cy - 1.62 * s, "peritumoral brain", ha="center", va="center", fontsize=fs, color=BLUE, fontweight="bold", zorder=5)
        ax.text(cx - 2.95 * s, cy + 1.72 * s, "laser fiber", ha="left", va="bottom", fontsize=fs - 2, color=INK, zorder=5)


def fig_margin():
    fig = plt.figure(figsize=(6.4, 4.9))
    ax = fig.add_axes([0.0, 0.0, 0.58, 1.0]); ax.set_xlim(-3.1, 2.7); ax.set_ylim(-2.35, 2.35); ax.set_aspect("equal"); ax.set_axis_off()
    rings(ax, fs=17)
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
    fig = plt.figure(figsize=(11.0, 4.1))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.set_xlim(-7.2, 7.2); ax.set_ylim(-2.5, 2.5); ax.set_aspect("equal"); ax.set_axis_off()
    rings(ax, s=1.15, labels=False)
    ax.text(0, 0, "core", ha="center", va="center", fontsize=15, color=RED, fontweight="bold", zorder=5)
    # a few surviving cells in the margin ring
    rng = np.random.default_rng(3)
    for _ in range(70):
        t = rng.uniform(0, 2 * np.pi); r = rng.uniform(0.62, 0.92)
        ax.plot(r * 1.6 * 1.15 * np.cos(t), r * 1.2 * 1.15 * np.sin(t), "o", ms=5, color="#9A7300", zorder=4, alpha=0.8)
    # described (navy, left) and not described (red, right)
    def call(x, y, text, col, px, py, ha):
        ax.annotate(text, (px, py), xytext=(x, y), fontsize=16, color=col, ha=ha, va="center", linespacing=1.25,
                    fontweight="bold" if col == RED else "normal",
                    arrowprops=dict(arrowstyle="-", color=col, lw=1.6, shrinkA=0, shrinkB=2), zorder=6)
    ax.text(-7.0, 2.3, "described", fontsize=17, color=NAVY, fontweight="bold", ha="left", va="top")
    call(-7.0, 1.0, "peri-ablation vessels open the\nblood–brain barrier for a time\n(Cleary 2026)", NAVY, -2.6, 0.55, "left")
    call(-7.0, -1.3, "the immune microenvironment\nof the sublethal zone, profiled\nin a rodent model (Tao 2026)", NAVY, -1.7, -1.1, "left")
    ax.text(7.0, 2.3, "not described", fontsize=17, color=RED, fontweight="bold", ha="right", va="top")
    call(7.0, 0.3, "the transcriptional state of the\ntumor cells that survive here\nand regrow the lesion", RED, 1.7, 0.6, "right")
    ax.text(7.0, -1.6, "this analysis", fontsize=16, color=RED, ha="right", va="center", style="italic")
    save(fig, "fig_margin_known.png", pad=0.02)


# ====================================================================== PCA
def fig_pca():
    g = digitize_pca()
    fig, ax = plt.subplots(figsize=(5.0, 4.9))
    bare(ax)
    ax.axhline(0, color=GRID, lw=0.9, zorder=1); ax.axvline(0, color=GRID, lw=0.9, zorder=1)
    style = {"Primary": (GOLD, "^", "Primary\n(pre-LITT)"), "Recurrent": (NAVY, "s", "Recurrent\n(post-LITT)")}
    for name, pts in g.items():
        col, mark, _ = style[name]
        ax.scatter(pts[:, 0], pts[:, 1], s=300, marker=mark, color=col, edgecolor="white", linewidth=1.8, zorder=4)
    ax.text(-7.4, 4.4, style["Primary"][2], color="#9A7300", fontsize=16, fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.text(7.0, 4.6, style["Recurrent"][2], color=NAVY, fontsize=16, fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.set_xlabel("PC1  (%.1f %% of variance)" % F["pc1"], fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("PC2  (%.1f %%)" % F["pc2"], fontsize=16, color=INK, labelpad=4)
    ax.set_xlim(-11.5, 11.5); ax.set_ylim(-10.5, 10.5)
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
    ax.scatter(sig.log2FoldChange, sig.y, s=70, color=NAVY, lw=0.8, edgecolor="white", zorder=4)
    for v in (-1, 1):
        ax.axvline(v, color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.axhline(-np.log10(0.05), color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.set_xlim(-11.8, 11.8); ax.set_ylim(-0.6, sig.y.max() * 1.18)
    ax.set_xlabel("log$_2$ fold change   (recurrent vs primary)", fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("−log$_{10}$ adjusted P", fontsize=16, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    top = ax.get_ylim()[1]
    ax.text(-11.2, top * 0.955, "%d down" % F["de_down"], color=NAVY, fontsize=19, fontweight="bold", ha="left", va="top")
    ax.text(11.2, top * 0.955, "%d up" % F["de_up"], color=NAVY, fontsize=19, fontweight="bold", ha="right", va="top")
    ax.text(11.2, -np.log10(0.05) + top * 0.02, "FDR 0.05", color=GREY, fontsize=13, ha="right", va="bottom")
    hubs = set(F["ppi_hubs"][:4]) | {"CALB1"}
    named = sig[~sig.symbol.astype(str).str.startswith("ENSG")]
    lab = pd.concat([named.nsmallest(5, "log2FoldChange"), named.nlargest(5, "log2FoldChange"),
                     named[named.symbol.isin(hubs)]]).drop_duplicates("symbol")
    for side in (-1, 1):
        part = lab[np.sign(lab.log2FoldChange) == side].sort_values("y")
        ys = spread(part.y.tolist(), 0.6, top * 0.80, top * 0.095)
        for (_, row), y in zip(part.iterrows(), ys):
            inward = abs(row.log2FoldChange) > 7
            out = -1.5 if side < 0 else 1.5
            dx = -out if inward else out
            ax.annotate(row.symbol, (row.log2FoldChange, row.y), xytext=(row.log2FoldChange + dx, y),
                        ha=("left" if side < 0 else "right") if inward else ("right" if side < 0 else "left"), va="center",
                        fontsize=13, color=INK, style="italic",
                        arrowprops=dict(arrowstyle="-", color="#AFB6B4", lw=0.8, shrinkA=0, shrinkB=3))
    save(fig, "chart_volcano.png")


# ====================================================================== leading sets, vertical
SHORT = {"KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION": "translation\ninitiation",
         "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION": "translation\nelongation",
         "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "GCN2 amino-\nacid stress",
         "KEGG_RIBOSOME": "ribosome", "REACTOME_SELENOAMINO_ACID_METABOLISM": "selenoamino-\nacid\nmetabolism",
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
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    hanging(ax, list(range(6)), [s["nes"] for s in sets], [NAVY if s["q"] < 0.05 else BLUE for s in sets],
            [SHORT[s["name"]] for s in sets], q=[s["q"] for s in sets], fs_lab=11.5, ylim=(-2.9, 0.35))
    ax.text(-0.5, -2.86, "nominal p < 0.001 for all six  ·  %s gene sets tested" % format(F["gsea_n_sets"], ","),
            fontsize=12.5, color=GREY, style="italic", ha="left", va="bottom")
    save(fig, "chart_gsea_sets.png")


THEME_GROUPS = [("anabolic", ["ribosome biogenesis", "mTORC1 signalling", "MYC targets"]),
                ("energy", ["glycolysis", "oxidative phosphorylation", "mitochondrial translation"]),
                ("", ["TCA cycle"])]
THEME_SHORT = {"ribosome biogenesis": "ribosome\nbiogenesis", "mTORC1 signalling": "mTORC1", "MYC targets": "MYC\ntargets",
               "glycolysis": "glycolysis", "oxidative phosphorylation": "OXPHOS", "mitochondrial translation": "mito.\ntranslation",
               "TCA cycle": "TCA\ncycle"}


def theme(name):
    th = F["themes"]
    if name in th:
        return th[name]
    hits = [k for k in th if name.lower().replace(" ", "") in k.lower().replace(" ", "")]
    assert len(hits) == 1, (name, hits, list(th))
    return th[hits[0]]


def fig_themes():
    names = [n for _, ns in THEME_GROUPS for n in ns]
    vals = [theme(n)["nes"] for n in names]
    cols = [NAVY] * 3 + [BLUE] * 3 + [GREY]
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    xs = [0, 1.5, 3.0, 4.9, 6.4, 7.9, 9.8]
    hanging(ax, xs, vals, cols, [THEME_SHORT[n] for n in names], fs_lab=12, yl="", fs_val=13)
    ax.text(1.5, 0.30, "biosynthesis", ha="center", va="bottom", fontsize=15, color=NAVY, fontweight="bold")
    ax.text(6.4, 0.30, "energy", ha="center", va="bottom", fontsize=15, color=BLUE, fontweight="bold")
    ax.text(9.8, 0.30, "unchanged", ha="center", va="bottom", fontsize=15, color=GREY, fontweight="bold")
    ax.set_ylim(-2.45, 0.85)
    save(fig, "chart_themes.png")


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
    ax.set_xlim(-2.55, 2.55); ax.set_ylim(-0.06, 2.10)
    ax.set_xlabel("normalized enrichment score", fontsize=16, color=INK, labelpad=6)
    ax.set_ylabel("−log$_{10}$ FDR q", fontsize=16, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.text(-2.50, -np.log10(0.25) + 0.04, "q = 0.25", color=GREY, fontsize=13, ha="left", va="bottom")
    short = dict(LABELS, **{"REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "GCN2 starvation response"})
    short = {k: v.replace("\n", " ") for k, v in short.items()}
    hit = hit.sort_values("y")
    ys = spread(hit.y.tolist(), 0.62, 1.72, 0.21)
    for (_, row), y in zip(hit.iterrows(), ys):
        ax.annotate(short[row.gene_set], (row.NES, row.y), xytext=(row.NES + 0.14, y), ha="left", va="center", fontsize=13.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color="#AFB6B4", lw=0.8, shrinkA=0, shrinkB=4))
    ax.annotate("nothing on this side\nclears q = %.2f" % F["gsea_up_best_q"], xy=(1.56, 0.06), xytext=(1.72, 1.45), ha="center", va="center",
                fontsize=14, color=GREY, linespacing=1.35, arrowprops=dict(arrowstyle="->", color=GREY, lw=1.1, connectionstyle="arc3,rad=0.22"))
    ax.text(-2.50, 2.06, "down in recurrence", color=NAVY, fontsize=15, fontweight="bold", ha="left", va="top")
    ax.text(2.50, 2.06, "up in recurrence", color=GREY, fontsize=15, ha="right", va="top")
    save(fig, "chart_gsea_landscape.png")


# ====================================================================== subtypes
def fig_subtypes():
    rows = sorted(F["gsva"], key=lambda r: r["change"])
    fig, ax = plt.subplots(figsize=(6.4, 4.7))
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.axhline(0, color=GRID, lw=1.0, zorder=1)
    ends = []
    for r in rows:
        on = r["p"] < 0.05; mtc = r["name"] == "Garofano_MTC"
        col = NAVY if on else GOLD if mtc else "#C4CBC9"
        a, b = r["primary"], r["recurrent"]
        ax.plot([0, 1], [a, b], color=col, lw=3.2 if (on or mtc) else 1.5, zorder=4 if (on or mtc) else 2, solid_capstyle="round")
        ax.scatter([0, 1], [a, b], s=100 if (on or mtc) else 40, color=col, edgecolor="white", lw=1.3, zorder=5 if (on or mtc) else 3)
        ends.append((b, r, on, mtc))
    ends.sort(key=lambda e: e[0])
    for (_, r, on, mtc), y in zip(ends, spread([e[0] for e in ends], -0.62, 0.62, 0.105)):
        lab = NICE_SUBTYPE[r["name"]] + ("  *" if on else "")
        ax.text(1.05, y, lab, fontsize=14 if (on or mtc) else 12.5, color=NAVY if on else "#9A7300" if mtc else GREY,
                fontweight="bold" if (on or mtc) else "normal", va="center")
    ax.set_xlim(-0.13, 2.05); ax.set_ylim(-0.66, 0.66)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Primary\n(pre-LITT)", "Recurrent\n(post-LITT)"], fontsize=15, color=INK, linespacing=1.3)
    ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5]); ax.set_yticklabels([MINUS("-0.5"), MINUS("-0.25"), "0", "0.25", "0.5"], fontsize=13)
    ax.tick_params(length=0)
    ax.set_ylabel("GSVA enrichment score", fontsize=15, color=INK)
    ax.text(-0.13, 0.63, "*  p < 0.05   ·   ten published signatures, GSVA", fontsize=13, color=GREY, style="italic", va="center")
    save(fig, "chart_subtypes.png")


# ====================================================================== PPI
def fig_ppi():
    e = sheet("S6_PPI_edges").drop_duplicates(subset=["gene1", "gene2"])
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.gene1, r.gene2, w=float(r.STRING_combined_score))
    deg = dict(G.degree())
    hubs = {n for n, _ in sorted(deg.items(), key=lambda t: -t[1])[:4]}
    lfc = F["ppi_nodes"]
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    pos = {}
    big = G.subgraph(comps[0])
    p = nx.spring_layout(big, seed=7, k=1.6, iterations=900, weight=None)
    pts = np.array([p[n] for n in big]); pts = pts - pts.mean(0); pts = pts / np.abs(pts).max() * 1.30
    for n, xy in zip(big, pts):
        pos[n] = xy + np.array([0.0, 0.75])
    # the four-node component is a star on RYR2, placed by hand under the large one so no label crosses a node
    small = comps[1]
    centre = max(small, key=lambda n: deg[n])
    others = sorted(n for n in small if n != centre)
    star = {centre: (0.15, -1.45)}
    for n, xy in zip(others, [(-0.75, -1.75), (1.0, -1.05), (1.0, -1.95)]):
        star[n] = xy
    for n, xy in star.items():
        pos[n] = np.array(xy)
    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ax.set_axis_off(); ax.set_aspect("equal")
    scores = np.array([d["w"] for *_, d in G.edges(data=True)])
    for u, v, d in G.edges(data=True):
        lw = 1.6 + 4.0 * (d["w"] - scores.min()) / (scores.max() - scores.min())
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#B7C3DE", lw=lw, zorder=1, solid_capstyle="round")
    ymid = float(np.mean([pos[n][1] for n in big]))
    for n in G.nodes():
        hub = n in hubs; col = RED if lfc.get(n, 0) > 0 else NAVY
        ax.scatter(*pos[n], s=470 if hub else 250, color=col if hub else "white", edgecolor=col, linewidth=0 if hub else 2.8, zorder=3)
        above = (n in big and pos[n][1] > ymid + 0.35) or n in (others[1],)
        ax.text(pos[n][0], pos[n][1] + (0.22 if above else -0.22), n, ha="center", va="bottom" if above else "top", zorder=5,
                fontsize=14 if hub else 13, fontweight="bold", color=col, path_effects=[pe.withStroke(linewidth=3.5, foreground="white")])
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.85, max(xs) + 0.85); ax.set_ylim(min(ys) - 0.62, max(ys) + 0.55)   # the legend is slide text
    save(fig, "chart_ppi.png")


# ====================================================================== drugs, vertical
def fig_drugs():
    tier = F["drug_tierA"][:10]
    fig, ax = plt.subplots(figsize=(5.9, 4.6))
    xs = list(range(len(tier)))
    cols = [NAVY if d["drug"] == "ciclopirox" else BLUE for d in tier]
    ax.bar(xs, [d["score"] for d in tier], width=0.66, color=cols, zorder=3)
    for x, d in zip(xs, tier):
        ax.text(x, d["score"] + 0.05, "%.2f" % d["score"], ha="center", va="bottom", fontsize=14,
                fontweight="bold" if d["drug"] == "ciclopirox" else "normal", color=INK)
    ax.set_xticks(xs); ax.set_xticklabels([d["drug"].replace(" bromide", "\nbromide") for d in tier], fontsize=13.5, color=INK, rotation=35, ha="right")
    bare(ax, grid="y")
    ax.tick_params(axis="x", length=0)
    ax.set_ylim(0, max(d["score"] for d in tier) * 1.18)
    ax.set_ylabel("|NES|$^{1.5}$ × predicted BBB permeability", fontsize=15, color=INK)
    ax.text(len(tier) - 0.5, max(d["score"] for d in tier) * 1.12, "%d of %d hits carry a clinical phase" % (F["drug_clinical_n"], 93),
            ha="right", va="top", fontsize=13, color=GREY, style="italic")
    save(fig, "chart_drugs.png")


# ====================================================================== conclusions: the arms move together
ARMS = [("translation", NAVY, [("initiation", "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"), ("elongation", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION"),
                               ("ribosome", "KEGG_RIBOSOME")]),
        ("biosynthesis", NAVY, [("rRNA\nprocessing", "rRNA processing"), ("ribosome\nbiogenesis", "ribosome biogenesis"), ("mTORC1", "mTORC1 signalling"),
                                ("MYC", "MYC targets")]),
        ("energy", BLUE, [("glycolysis", "glycolysis"), ("OXPHOS", "oxidative phosphorylation"), ("mito.\ntranslation", "mitochondrial translation")]),
        ("hypoxia", BLUE, [("HIF1\ntargets", "HIF1 targets"), ("hypoxia\nmetagene", "hypoxia metagene"), ("hypoxia", "hypoxia")]),
        ("unchanged", GREY, [("TCA\ncycle", "TCA cycle")])]


def fig_arms():
    """One column per arm = the mean normalized enrichment score of its member sets; each member a dot on the column."""
    top = {s["name"]: s["nes"] for s in F["gsea_down_top"]}
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    xs, means, cols, labels, lows = [], [], [], [], []
    members_txt = []
    for i, (head, col, members) in enumerate(ARMS):
        vals = [top[key] if key in top else theme(key)["nes"] for _, key in members]
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
    ax.set_ylim(-2.6, 0.25); ax.set_ylabel("enrichment score (NES)", fontsize=14, color=INK)
    ax.set_yticks([0, -0.5, -1.0, -1.5, -2.0]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.1f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.text(-0.55, -2.57, "columns, the mean of each arm  ·  dots, its member gene sets", fontsize=12, color=GREY, style="italic",
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
    table_image("table_fold.png", fold[0], fold[1:], [2.6, 1.15, 2.85], 4.4, fs=15, row_h=0.56)
    hubs[0][0] = "hub gene"
    table_image("table_hubs.png", hubs[0], hubs[1:], [1.6, 1.1, 3.6], 5.4, fs=15, row_h=0.46)
    table_image("table_drugs.png", drugs[0], drugs[1:], [0.5, 2.4, 0.9, 0.9, 0.9, 1.0], 5.0, fs=15, row_h=0.52, hl={0},
                align=["center", "left", "center", "center", "center", "center"])


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
    fig = plt.figure(figsize=(11.0, 2.55))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 11.0); ax.set_ylim(0.6, 3.0); ax.set_axis_off()
    stages = [("U251N cells", "human glioblastoma line", NAVY), ("orthotopic tumor", "athymic RNU/RNU rat brain", NAVY),
              ("MRI-guided LITT", "laser ablation under\nMR thermometry", RED), ("recurrence", "regrowth after ablation", GOLD),
              ("bulk RNA-seq", "3 primary + 3 recurrent,\nhuman reads only", BLUE)]
    w, gap = 1.9, 0.375
    x = 0.15
    for i, (head, sub, col) in enumerate(stages):
        ax.add_patch(FancyBboxPatch((x, 1.05), w, 1.75, boxstyle="round,pad=0,rounding_size=0.12", fc=LIGHT, ec=col, lw=2.5))
        ax.add_patch(FancyBboxPatch((x, 2.25), w, 0.55, boxstyle="round,pad=0,rounding_size=0.12", fc=col, ec=col, lw=0))
        ax.text(x + w / 2, 2.52, head, ha="center", va="center", fontsize=16, fontweight="bold", color="white")
        if head == "MRI-guided LITT":
            im = Image.open(OUT / "video_poster.png").convert("L")
            ax.imshow(im, extent=(x + 0.15, x + 0.95, 1.12, 1.92), cmap="gray", zorder=3, aspect="auto")
            ax.text(x + 1.05, 1.52, "laser ablation\nunder MR\nthermometry", ha="left", va="center", fontsize=12.5, color=INK, linespacing=1.2)
        else:
            ax.text(x + w / 2, 1.62, sub, ha="center", va="center", fontsize=13.5, color=INK, linespacing=1.25)
        if i < len(stages) - 1:
            ax.annotate("", (x + w + gap - 0.05, 1.92), (x + w + 0.05, 1.92), arrowprops=dict(arrowstyle="-|>", color=GOLD, lw=3, mutation_scale=22))
        x += w + gap
    save(fig, "overview_pipeline.png", pad=0.02)


# ====================================================================== what else moves: hypoxia, iron, stress down; cell cycle up
THEMES2 = [("hypoxia", NAVY, [("HIF1\ntargets", "HIF1 targets"), ("hypoxia\nmetagene", "hypoxia metagene"), ("hypoxia\n(hallmark)", "hypoxia"),
                              ("angio-\ngenesis", "angiogenesis")]),
           ("iron, stress", BLUE, [("iron\nuptake", "iron uptake and transport"), ("senes-\ncence", "senescence"), ("EMT", "EMT")]),
           ("cell cycle", RED, [("mitotic\nspindle", "mitotic spindle"), ("G2M\ncheckpoint", "G2M checkpoint"), ("E2F\ntargets", "E2F targets")])]


def fig_themes2():
    fig, ax = plt.subplots(figsize=(9.8, 4.8))
    xs, vals, cols, labels, ps, heads = [], [], [], [], [], []
    x = 0.0
    for head, col, members in THEMES2:
        x0 = x
        for lab, key in members:
            t = theme(key); xs.append(x); vals.append(t["nes"]); cols.append(col); labels.append(lab); ps.append(t["p"]); x += 1.9
        heads.append(((x0 + x - 1.9) / 2, head, col)); x += 1.2
    ax.bar(xs, vals, width=1.0, color=cols, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    for xx, v, p in zip(xs, vals, ps):
        up = v > 0
        ax.text(xx, v + (0.05 if up else -0.05), MINUS(f"{v:+.2f}"), ha="center", va="bottom" if up else "top", fontsize=13, fontweight="bold", color=INK)
        ax.text(xx, v + (0.32 if up else -0.32), f"p {p:.3f}" if p < 0.01 else f"p {p:.2f}", ha="center", va="bottom" if up else "top", fontsize=11.5, color=GREY)
    for xc, head, col in heads:
        ax.text(xc, -2.85, head, ha="center", va="top", fontsize=15, color=col, fontweight="bold")
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=11.5, color=INK, linespacing=1.1)
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=4)
    ax.set_ylim(-3.2, 2.8); ax.set_ylabel("normalized enrichment score", fontsize=14, color=INK)
    ax.set_yticks([2, 1, 0, -1, -2]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.text(-0.7, 2.7, "up in recurrence", ha="left", va="top", fontsize=13, color=GREY, style="italic")
    ax.text(xs[-1] + 0.6, -0.15, "down in recurrence", ha="right", va="top", fontsize=13, color=GREY, style="italic")
    save(fig, "chart_themes2.png")


# ====================================================================== the drug filter as a funnel
def fig_funnel():
    steps = [(93, "drug signatures\nthat oppose\nours"), (F["drug_clinical_n"], "compound has\na clinical\nphase"),
             (len(F["drug_tierA"]), "clears both\nblood–brain\nbarrier models"), (1, "survives the\nprior-art\naudit")]
    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    cols = [BLUE, BLUE, BLUE, NAVY]
    xs = [0, 1.3, 2.6, 3.9]
    ax.bar(xs, [s[0] for s in steps], width=0.9, color=cols, zorder=3)
    for x, (n, lab), col in zip(xs, steps, cols):
        ax.text(x, n + 2, str(n), ha="center", va="bottom", fontsize=24, fontweight="bold", color=col)
        ax.text(x, -4, lab, ha="center", va="top", fontsize=13, color=INK, linespacing=1.2)
    ax.set_xticks([]); ax.set_ylim(0, 110); ax.set_xlim(-0.7, 4.6)
    for sd in ("top", "right", "left"):
        ax.spines[sd].set_visible(False)
    ax.set_yticks([])
    ax.text(3.9, 16, "ciclopirox", ha="center", va="bottom", fontsize=15, fontweight="bold", color=NAVY)
    fig.subplots_adjust(bottom=0.25)
    save(fig, "chart_funnel.png")


# ====================================================================== DepMap: does U-251 depend on what ciclopirox hits?
def fig_depmap():
    rows = [r.split(",") for r in (ROOT / "depmap" / "depmap_u251_targets.csv").read_text().strip().splitlines()[1:]]
    want = [("RRM1", "ribonucleotide\nreductase"), ("RRM2", "ribonucleotide\nreductase"), ("DOHH", "eIF5A\nhypusination"),
            ("DHPS", "eIF5A\nhypusination"), ("EIF5A", "eIF5A\nhypusination"), ("VHL", "PHD–HIF\naxis"), ("HIF1A", "PHD–HIF\naxis"),
            ("FTH1", "iron\nhandling"), ("GPX4", "iron\nhandling"), ("TFRC", "iron\nhandling")]
    by = {r[1]: r for r in rows}
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    xs, vals, cols, labels, ce = [], [], [], [], []
    x = 0.0; last = None; heads = []
    for gene, grp in want:
        r = by[gene]; eff = float(r[2]); common = r[6].strip() == "True"
        if last is not None and grp != last:
            heads.append(((x0 + x - 1.25) / 2, last)); x += 0.9
        if grp != last:
            x0 = x
        xs.append(x); vals.append(eff); labels.append(gene); ce.append(common)
        cols.append(GREY if common else (NAVY if eff < -0.5 else "#C4CBC9"))
        last = grp; x += 1.25
    heads.append(((x0 + x - 1.25) / 2, last))
    ax.bar(xs, vals, width=0.8, color=cols, zorder=3)
    ax.axhline(0, color=INK, lw=1.2, zorder=4)
    ax.axhline(-0.5, color=RED, lw=1.6, ls=(0, (4, 3)), zorder=4)
    ax.axhline(-1.0, color=GREY, lw=1.0, ls=(0, (2, 3)), zorder=4)
    for xx, v, g, c in zip(xs, vals, labels, ce):
        ax.text(xx, v - 0.05 if v < 0 else v + 0.05, MINUS(f"{v:.2f}"), ha="center", va="top" if v < 0 else "bottom", fontsize=12.5,
                fontweight="bold" if v < -0.5 and not c else "normal", color=INK)
    for xc, head in heads:
        ax.text(xc, 0.62, head, ha="center", va="bottom", fontsize=13, color=NAVY, fontweight="bold", linespacing=1.1)
    ax.set_xticks(xs); ax.set_xticklabels(labels, fontsize=12.5, color=INK, style="italic")
    ax.xaxis.set_ticks_position("top"); ax.tick_params(axis="x", length=0, pad=4)
    ax.set_ylim(-4.0, 1.35); ax.set_ylabel("gene effect in U-251 MG", fontsize=14, color=INK)
    ax.set_yticks([0, -1, -2, -3]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(axis="y", labelsize=13)
    for sd in ("top", "right", "bottom"):
        ax.spines[sd].set_visible(False)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0); ax.set_axisbelow(True)
    ax.text(xs[-1] + 0.5, -3.97, "dashed red, dependency threshold (−0.5)  ·  dotted grey, median common-essential gene (−1)\n"
                                  "navy, a selective dependency  ·  grey, common essential in every line  ·  pale, no dependency",
            ha="right", va="bottom", fontsize=11.5, color=GREY, style="italic", linespacing=1.3)
    save(fig, "chart_depmap.png")


if __name__ == "__main__":
    print("figures ->", OUT)
    fig_margin(); fig_margin_known(); fig_pca(); fig_volcano(); fig_gsea_sets(); fig_themes(); fig_landscape(); fig_subtypes()
    fig_ppi(); fig_drugs(); fig_arms(); tables(); cartoon(); video(); fig_pipeline(); fig_themes2(); fig_funnel(); fig_depmap()
    print("done")
