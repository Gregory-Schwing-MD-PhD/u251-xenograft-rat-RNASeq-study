# -*- coding: utf-8 -*-
"""Slide figures for the CNS 2026 deck.

The published 9-panel composite was drawn for a 10 x 10 inch print figure; shrunk
onto a slide its type lands around 6-7 pt, which is unreadable from the back of a
session room. So every panel that has source data behind it is REDRAWN here at
the exact size it occupies on the slide (figure inches == slide inches, so a
13 pt label in matplotlib is a 13 pt label on the slide).

Source data:
  Supplementary_Data.xlsx  S1 (all genes) S2 (significant) S3 (GSEA)
                           S4 (subtype scores) S6 (STRING edges) S7 (drug ranking)
  figures/facts.json       scalars, written by 01_extract_facts.py

The one exception is the PCA: the variance-stabilized count matrix is not in this
repo, so the six sample coordinates are recovered from the published panel by
locating the plot gridlines and the marker centroids (digitize_pca below). That
is still derived from data at build time, not typed in.

No numbers are hand-entered in this file.
"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter
from PIL import Image, ImageChops
from scipy import ndimage

Image.MAX_IMAGE_PIXELS = None

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(HERE, "figures")
XLSX = os.path.join(ROOT, "manuscript", "Supplementary_Data.xlsx")
COMPOSITE = os.path.join(ROOT, "publication_figure",
                         "Publication_Figure_9Panel_VOLCANO_COMPLETE.png")
F = json.load(open(os.path.join(FIG, "facts.json")))

# Deck palette. GREEN is the Wayne State mark used throughout the template, so
# data marks stay on it rather than on a generic categorical ramp. The one
# two-category chart (PCA) pairs it with OCHRE: that pair clears CVD separation
# (dE 17.6 protan / 25.4 normal) and 3:1 contrast, and carries distinct marker
# shapes plus direct labels so identity is never colour-alone.
GREEN, OCHRE = "#095945", "#B07A1E"
MIDG, PALE = "#6E9E92", "#E8F0EE"
INK, MUTED, FAINT = "#1A1A1A", "#5A5A5A", "#C9CCCB"
GRID = "#DCDCDC"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 13,
    "axes.edgecolor": "#BFBFBF", "axes.linewidth": 0.9,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
    "figure.dpi": 300, "savefig.dpi": 300,
})
MINUS = lambda s: s.replace("-", "−")


def sheet(name):
    return pd.read_excel(XLSX, name)


def save(fig, name):
    path = os.path.join(FIG, name)
    fig.savefig(path, facecolor="white", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print("  %-28s %s" % (name, Image.open(path).size))


def bare(ax, grid="both"):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    if grid:
        ax.grid(axis=grid, color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, labelsize=12)


def spread(pos, lo, hi, gap):
    """Push a sorted list of label positions apart so none collide."""
    pos = list(pos)
    for _ in range(400):
        moved = False
        for i in range(len(pos) - 1):
            d = pos[i + 1] - pos[i]
            if d < gap:
                shift = (gap - d) / 2
                pos[i] -= shift
                pos[i + 1] += shift
                moved = True
        pos = [min(max(p, lo), hi) for p in pos]
        if not moved:
            break
    return pos


# ===================================================================== PCA ===
def digitize_pca():
    """Recover the six PCA sample coordinates from the published panel A.

    ggplot draws one major gridline per axis tick, so the three evenly spaced
    light-grey full-height/width lines inside the panel are the -5 / 0 / +5 ticks;
    that calibrates pixels to PC units. Markers are the two solid fills used by
    the panel (orange triangles = primary, red squares = recurrent).
    """
    comp = Image.open(COMPOSITE).convert("RGB")
    w, h = comp.width // 3, comp.height // 3
    # panel A is the top-left cell; its left ~70% is the PCA, the rest the scree
    a = np.asarray(comp.crop((0, 0, int(w * 0.70), h))).astype(int)

    grey = (abs(a[..., 0] - a[..., 1]) < 6) & (abs(a[..., 1] - a[..., 2]) < 6)
    dark = grey & (a[..., 0] < 160)
    light = grey & (a[..., 0] > 190) & (a[..., 0] < 240)

    def runs(counts, floor):
        hits = np.where(counts > floor)[0]
        out, cur = [], [hits[0]]
        for i in hits[1:]:
            if i - cur[-1] <= 3:
                cur.append(i)
            else:
                out.append(sum(cur) / len(cur))
                cur = [i]
        out.append(sum(cur) / len(cur))
        return out

    vb = runs(dark.sum(0), dark.shape[0] * 0.75)          # panel border, vertical
    hb = runs(dark.sum(1), dark.shape[1] * 0.30)          # panel border, horizontal
    x_lo, x_hi, y_lo, y_hi = vb[0], vb[-1], hb[0], hb[-1]

    def evenly_spaced_triple(cands):
        """Ticks are equally spaced, so take the first consecutive triple whose
        two gaps match; that drops axis furniture caught by the same threshold."""
        for i in range(len(cands) - 2):
            a_, b_, c_ = cands[i:i + 3]
            if abs((b_ - a_) - (c_ - b_)) < 0.02 * (c_ - a_):
                return [a_, b_, c_]
        raise AssertionError("no evenly spaced gridline triple in %s" % cands)

    gv = evenly_spaced_triple(
        [x for x in runs(light.sum(0), dark.shape[0] * 0.75) if x_lo < x < x_hi])
    gh = evenly_spaced_triple(
        [y for y in runs(light.sum(1), dark.shape[1] * 0.30) if y_lo < y < y_hi])
    px2pc1 = lambda px: -5.0 + (px - gv[0]) * 10.0 / (gv[2] - gv[0])
    px2pc2 = lambda py: 5.0 + (py - gh[0]) * -10.0 / (gh[2] - gh[0])

    groups = {}
    for label, rgb in (("Primary", (255, 127, 14)), ("Recurrent", (214, 39, 40))):
        m = np.all(abs(a - np.array(rgb)) < 12, axis=2)
        lab, n = ndimage.label(m)
        sizes = ndimage.sum(m, lab, range(1, n + 1))
        keep = [i + 1 for i, s in enumerate(sizes) if s > 200]
        pts = []
        for cy, cx in ndimage.center_of_mass(m, lab, keep):
            if y_lo < cy < y_hi and x_lo < cx < x_hi:      # drop the legend key
                pts.append((px2pc1(cx), px2pc2(cy)))
        assert len(pts) == 3, (label, pts)
        groups[label] = np.array(pts)
    return groups


def fig_pca():
    g = digitize_pca()
    fig, ax = plt.subplots(figsize=(4.9, 4.7))
    bare(ax)
    ax.axhline(0, color=GRID, lw=0.9, zorder=1)
    ax.axvline(0, color=GRID, lw=0.9, zorder=1)
    style = {"Primary": (OCHRE, "^", "Primary\n(pre-LITT)"),
             "Recurrent": (GREEN, "s", "Recurrent\n(post-LITT)")}
    for name, pts in g.items():
        col, mark, _ = style[name]
        ax.scatter(pts[:, 0], pts[:, 1], s=250, marker=mark, color=col,
                   edgecolor="white", linewidth=1.8, zorder=4)
    ax.text(-7.4, 4.4, style["Primary"][2], color=OCHRE, fontsize=13.5,
            fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.text(7.0, 4.6, style["Recurrent"][2], color=GREEN, fontsize=13.5,
            fontweight="bold", ha="center", va="center", linespacing=1.25)
    ax.set_xlabel("PC1  (%.1f%% of variance)" % F["pc1"], fontsize=13.5, color=INK,
                  labelpad=6)
    ax.set_ylabel("PC2  (%.1f%%)" % F["pc2"], fontsize=13.5, color=INK, labelpad=4)
    ax.set_xlim(-11.5, 11.5)
    ax.set_ylim(-10.5, 10.5)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(labelsize=12.5)
    save(fig, "chart_pca.png")


# ================================================================= volcano ===
def fig_volcano():
    d = sheet("S1_DE_all").dropna(subset=["padj", "log2FoldChange"])
    d["y"] = -np.log10(d.padj.clip(lower=1e-16))
    sig = d[(d.padj < 0.05) & (d.log2FoldChange.abs() > 1)]
    ns = d.drop(sig.index)

    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    bare(ax)
    ax.scatter(ns.log2FoldChange, ns.y, s=7, color=FAINT, lw=0, zorder=2)
    ax.scatter(sig.log2FoldChange, sig.y, s=60, color=GREEN, lw=0.8,
               edgecolor="white", zorder=4)
    for v in (-1, 1):
        ax.axvline(v, color=MUTED, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.axhline(-np.log10(0.05), color=MUTED, lw=0.9, ls=(0, (4, 4)), zorder=1)

    ax.set_xlim(-11.8, 11.8)
    ax.set_ylim(-0.6, sig.y.max() * 1.18)
    ax.set_xlabel("log$_2$ fold change   (recurrent vs primary)", fontsize=13.5,
                  color=INK, labelpad=6)
    ax.set_ylabel("−log$_{10}$ adjusted P", fontsize=13.5, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(labelsize=12.5)

    top = ax.get_ylim()[1]
    ax.text(-11.2, top * 0.955, "%d down" % F["de_down"], color=GREEN, fontsize=15,
            fontweight="bold", ha="left", va="top")
    ax.text(11.2, top * 0.955, "%d up" % F["de_up"], color=GREEN, fontsize=15,
            fontweight="bold", ha="right", va="top")
    ax.text(11.2, -np.log10(0.05) + top * 0.02, "FDR 0.05", color=MUTED,
            fontsize=11.5, ha="right", va="bottom")

    # name the extremes on each side plus the network hubs
    hubs = set(F["ppi_hubs"][:4]) | {"CALB1"}
    named = sig[~sig.symbol.astype(str).str.startswith("ENSG")]   # skip unnamed loci
    lab = pd.concat([named.nsmallest(5, "log2FoldChange"),
                     named.nlargest(5, "log2FoldChange"),
                     named[named.symbol.isin(hubs)]]).drop_duplicates("symbol")
    for side in (-1, 1):
        part = lab[np.sign(lab.log2FoldChange) == side].sort_values("y")
        ys = spread(part.y.tolist(), 0.6, top * 0.80, top * 0.088)
        for (_, row), y in zip(part.iterrows(), ys):
            # points out at the edges get their label turned inward so it stays
            # inside the axes instead of hanging off the figure
            inward = abs(row.log2FoldChange) > 7
            out = -1.5 if side < 0 else 1.5
            dx = -out if inward else out
            ax.annotate(row.symbol, (row.log2FoldChange, row.y),
                        xytext=(row.log2FoldChange + dx, y),
                        ha=("left" if side < 0 else "right") if inward
                        else ("right" if side < 0 else "left"), va="center",
                        fontsize=12, color=INK, style="italic",
                        arrowprops=dict(arrowstyle="-", color="#AFB6B4", lw=0.8,
                                        shrinkA=0, shrinkB=3))
    save(fig, "chart_volcano.png")


# ================================================================ GSEA bar ===
LABELS = {
    "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION": "Translation initiation",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION": "Translation elongation",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY":
        "EIF2AK4 (GCN2) response to\namino-acid deficiency",
    "KEGG_RIBOSOME": "Ribosome",
    "REACTOME_SELENOAMINO_ACID_METABOLISM": "Selenoamino-acid metabolism",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION": "Cellular response to starvation",
}


def fig_gsea_bar():
    sets = [s for s in F["gsea_down_top"] if s["name"] in LABELS][:6]
    assert len(sets) == 6, [s["name"] for s in sets]
    fig, ax = plt.subplots(figsize=(7.9, 4.0))
    y = list(range(len(sets)))
    ax.barh(y, [s["nes"] for s in sets], height=0.62, zorder=3,
            color=[GREEN if s["q"] < 0.05 else MIDG for s in sets])
    ax.set_xlim(0, -2.95)
    ax.set_ylim(len(sets) - 0.5, -0.95)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS[s["name"]] for s in sets], color=INK, fontsize=13,
                       linespacing=1.15)
    ax.set_xticks([0, -0.5, -1.0, -1.5, -2.0])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.1f}")))
    bare(ax, grid="x")
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(labelsize=12.5)
    ax.set_xlabel("Normalized enrichment score   (negative = down in recurrence)",
                  fontsize=13, color=INK, labelpad=8)
    for i, s in enumerate(sets):
        ax.text(s["nes"] - 0.07, i, MINUS(f"{s['nes']:.2f}"), va="center",
                ha="left", fontsize=13, color=INK, fontweight="bold")
        ax.text(-2.93, i, f"q = {s['q']:.3f}" if s["q"] < 0.05 else f"q = {s['q']:.2f}",
                va="center", ha="right", fontsize=12.5,
                color=INK if s["q"] < 0.05 else MUTED,
                fontweight="bold" if s["q"] < 0.05 else "normal")
    ax.text(-2.93, -0.9,
            "nominal p < 0.001 for all six   ·   %s gene sets tested"
            % format(F["gsea_n_sets"], ","),
            va="center", ha="right", fontsize=12, color=MUTED, style="italic")
    save(fig, "chart_gsea.png")


# =========================================================== GSEA landscape ==
def fig_gsea_landscape():
    g = sheet("S3_GSEA_pipeline").copy()
    g["y"] = -np.log10(g.FDR_q.clip(lower=1e-3))
    hit = g[g.FDR_q < 0.25].sort_values("NES")
    rest = g.drop(hit.index)

    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    bare(ax)
    ax.scatter(rest.NES, rest.y, s=6, color=FAINT, lw=0, zorder=2)
    ax.scatter(hit.NES, hit.y, s=95, color=GREEN, edgecolor="white", lw=1.0,
               zorder=5)
    ax.axhline(-np.log10(0.25), color=MUTED, lw=0.9, ls=(0, (4, 4)), zorder=1)
    ax.axvline(0, color=GRID, lw=1.0, zorder=1)

    ax.set_xlim(-2.55, 2.55)
    ax.set_ylim(-0.06, 2.10)
    ax.set_xlabel("Normalized enrichment score", fontsize=13.5, color=INK,
                  labelpad=6)
    ax.set_ylabel("−log$_{10}$ FDR q", fontsize=13.5, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.tick_params(labelsize=12.5)
    ax.text(-2.50, -np.log10(0.25) + 0.04, "q = 0.25", color=MUTED, fontsize=11.5,
            ha="left", va="bottom")

    short = dict(LABELS, **{
        "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY":
            "EIF2AK4 (GCN2) starvation response"})
    short = {k: v.replace("\n", " ") for k, v in short.items()}
    hit = hit.sort_values("y")
    ys = spread(hit.y.tolist(), 0.62, 1.72, 0.195)
    for (_, row), y in zip(hit.iterrows(), ys):
        ax.annotate(short[row.gene_set], (row.NES, row.y),
                    xytext=(row.NES + 0.14, y), ha="left", va="center",
                    fontsize=12, color=INK,
                    arrowprops=dict(arrowstyle="-", color="#AFB6B4", lw=0.8,
                                    shrinkA=0, shrinkB=4))
    ax.annotate("nothing on this side\nclears q = %.2f" % F["gsea_up_best_q"],
                xy=(1.56, 0.06), xytext=(1.72, 1.45), ha="center", va="center",
                fontsize=12.5, color=MUTED, linespacing=1.35,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.1,
                                connectionstyle="arc3,rad=0.22"))
    ax.text(-2.50, 2.06, "down in recurrence", color=GREEN, fontsize=13,
            fontweight="bold", ha="left", va="top")
    ax.text(2.50, 2.06, "up in recurrence", color=MUTED, fontsize=13,
            ha="right", va="top")
    save(fig, "chart_gsea_landscape.png")


# ================================================================ subtypes ===
NICE_SUBTYPE = {
    "Garofano_MTC": "Garofano mitochondrial", "Garofano_GPM": "Garofano glycolytic",
    "Garofano_NEU": "Garofano neuronal", "Garofano_PPR": "Garofano proliferative",
    "Neftel_AC": "Neftel astrocyte-like", "Neftel_OPC": "Neftel OPC-like",
    "Neftel_NPC1": "Neftel NPC1-like", "Neftel_NPC2": "Neftel NPC2-like",
    "Neftel_MES1": "Neftel MES1-like", "Neftel_MES2": "Neftel MES2-like",
}


def fig_subtypes():
    """GSVA over the published signatures; the earlier version scored 5-6 gene
    marker panels by mean z, which is what produced the spurious shift."""
    rows = sorted(F["gsva"], key=lambda r: r["change"])
    fig, ax = plt.subplots(figsize=(7.3, 4.4))
    for side in ("top", "right", "bottom", "left"):
        ax.spines[side].set_visible(False)
    ax.axhline(0, color=GRID, lw=1.0, zorder=1)

    ends = []
    for r in rows:
        on = r["p"] < 0.05
        a, b = r["primary"], r["recurrent"]
        ax.plot([0, 1], [a, b], color=GREEN if on else "#C4CBC9",
                lw=3.0 if on else 1.5, zorder=4 if on else 2,
                solid_capstyle="round")
        ax.scatter([0, 1], [a, b], s=90 if on else 40,
                   color=GREEN if on else "#C4CBC9", edgecolor="white",
                   lw=1.3 if on else 0.7, zorder=5 if on else 3)
        ends.append((b, r, on))
    ends.sort()
    for (_, r, on), y in zip(ends, spread([e[0] for e in ends], -0.62, 0.62,
                                          0.098)):
        lab = NICE_SUBTYPE[r["name"]] + ("  *" if on else "")
        ax.text(1.05, y, lab, fontsize=12.5 if on else 11.5,
                color=GREEN if on else MUTED,
                fontweight="bold" if on else "normal", va="center")

    ax.set_xlim(-0.13, 2.05)
    ax.set_ylim(-0.66, 0.66)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Primary\n(pre-LITT)", "Recurrent\n(post-LITT)"],
                       fontsize=13, color=INK, linespacing=1.3)
    ax.set_yticks([-0.5, -0.25, 0, 0.25, 0.5])
    ax.set_yticklabels([MINUS("-0.5"), MINUS("-0.25"), "0", "0.25", "0.5"],
                       fontsize=12.5)
    ax.tick_params(length=0)
    ax.set_ylabel("GSVA enrichment score", fontsize=13, color=INK)
    ax.text(-0.13, 0.63, "*  p < 0.05   ·   published signatures, GSVA",
            fontsize=11.5, color=MUTED, style="italic", va="center")
    save(fig, "chart_subtypes.png")


# ===================================================================== PPI ===
def fig_ppi():
    e = sheet("S6_PPI_edges").drop_duplicates(subset=["gene1", "gene2"])
    G = nx.Graph()
    for _, r in e.iterrows():
        G.add_edge(r.gene1, r.gene2, w=float(r.STRING_combined_score))
    deg = dict(G.degree())
    hubs = {n for n, _ in sorted(deg.items(), key=lambda t: -t[1])[:4]}

    # lay each connected component out on its own, then place them side by side,
    # so a single spring layout does not squash the small one into a blob
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    boxes = [(-1.42, 0.0, 1.30), (1.30, 0.0, 0.72)]   # cx, cy, half-width
    pos = {}
    for comp, (cx, cy, half) in zip(comps, boxes):
        sub = G.subgraph(comp)
        p = nx.spring_layout(sub, seed=7, k=1.6, iterations=900, weight=None)
        pts = np.array([p[n] for n in sub])
        pts = pts - pts.mean(0)
        pts = pts / np.abs(pts).max() * half
        for n, xy in zip(sub, pts):
            pos[n] = xy + np.array([cx, cy])

    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    ax.set_axis_off()
    ax.set_aspect("equal")
    scores = np.array([d["w"] for *_, d in G.edges(data=True)])
    for u, v, d in G.edges(data=True):
        lw = 1.6 + 4.0 * (d["w"] - scores.min()) / (scores.max() - scores.min())
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color="#9FB8B1", lw=lw, zorder=1, solid_capstyle="round")
    for n in G.nodes():
        hub = n in hubs
        ax.scatter(*pos[n], s=430 if hub else 230,
                   color=GREEN if hub else "white",
                   edgecolor=GREEN, linewidth=0 if hub else 2.6, zorder=3)
        # label sits below the node with a white halo so edges cannot cut it
        ax.text(pos[n][0], pos[n][1] - 0.20, n, ha="center", va="top", zorder=5,
                fontsize=13 if hub else 12, fontweight="bold",
                color=GREEN if hub else INK,
                path_effects=[pe.withStroke(linewidth=3.5, foreground="white")])
    xs = [p[0] for p in pos.values()]
    ys = [p[1] for p in pos.values()]
    ax.set_xlim(min(xs) - 0.75, max(xs) + 0.75)
    ax.set_ylim(min(ys) - 0.60, max(ys) + 0.95)
    ax.text(0.5, 0.995,
            "%d of the %d differentially expressed genes have a STRING interaction\n"
            "filled node = highest degree   ·   line weight = STRING combined score"
            % (G.number_of_nodes(), F["de_total"]),
            transform=ax.transAxes, ha="center", va="top", fontsize=12.5,
            color=MUTED, linespacing=1.5)
    save(fig, "chart_ppi.png")


# =================================================================== drugs ===
def fig_drugs():
    """Tier A: reached a clinic, and both BBB predictors call it permeant."""
    tier = F["drug_tierA"][:10][::-1]
    names = [d["drug"] for d in tier]
    score = [d["score"] for d in tier]
    bbb = [d["bbb"] for d in tier]

    fig, ax = plt.subplots(figsize=(5.48, 4.3))
    y = list(range(len(tier)))
    cmap = plt.get_cmap("BuGn")
    ax.barh(y, score, height=0.66, zorder=3,
            color=[cmap(0.30 + 0.55 * b) for b in bbb],
            edgecolor=GREEN, linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=12.5, color=INK)
    bare(ax, grid="x")
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, max(score) * 1.18)
    ax.set_xlabel("|NES|$^{1.5}$ × predicted BBB permeability", fontsize=13,
                  color=INK, labelpad=6)
    ax.tick_params(axis="y", length=0, pad=2)
    for yi, (v, d) in enumerate(zip(score, tier)):
        ax.text(v + max(score) * 0.015, yi, "%.2f" % v, va="center", ha="left",
                fontsize=12, color=INK,
                fontweight="bold" if d["drug"] == "ciclopirox" else "normal")
    # the compound-count caption lives in the slide text, not on the plot
    save(fig, "chart_drugs.png")


# ================================================================== driver ===
if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    print("figures ->", FIG)
    fig_pca()
    fig_volcano()
    fig_gsea_bar()
    fig_gsea_landscape()
    fig_subtypes()
    fig_ppi()
    fig_drugs()
