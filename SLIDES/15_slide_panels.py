# -*- coding: utf-8 -*-
"""Single-panel slide images (deck v6) from the v5 figures.

Greg's v5 feedback: a four-subplot figure cannot be read on a slide, so every slide gets ONE panel and that
panel fills the slide. This script does not draw anything and computes no number. It

  (a) re-runs MBR/nature_draft/figures/make_figures_rna.py -- the manuscript's own matplotlib code, with the same
      raised font constants 13_v5_figures.py uses -- at a higher save dpi, and writes, besides the whole figure,
      one PNG per axes cropped to that axes' tight bounding box (legend included, because the legend is a child
      of the axes);

  (b) cuts the R/ggplot rasters, which cannot be re-drawn, into their panel grid at the white gutters, and
      re-attaches -- as pixels lifted from the same PNG, never redrawn -- the strips a facetted panel shares with
      its siblings: the y-axis title and tick labels held in the left margin, a legend strip above the row, the
      single x-axis title centred under the middle panel.

Cut positions are hard-coded from the ink/white structure measured on the v5 PNGs; each entry asserts the source
image size, so a re-made source that changed shape fails loudly instead of cropping the wrong place.

    python SLIDES/15_slide_panels.py

Re-runnable: figures_cns/v6/ is deleted and rebuilt every time.
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
from matplotlib.transforms import Bbox
from PIL import Image

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DRAFT = ROOT / "MBR" / "nature_draft" / "figures"        # untracked on this machine
V5 = HERE / "figures_cns" / "v5"
OUT = HERE / "figures_cns" / "v6"

MIN_LONG_SIDE = 1200                                     # a slide panel below this is not worth projecting
UPSCALE_TO = 1300                                        # LANCZOS target when a native crop falls short

# ------------------------------------------------------------------ (a) the manuscript's matplotlib panels
# type exactly as v5 drew it; only the save dpi rises, so a cropped single axes still carries >1200 px
SLIDE_RC = {
    "font.size": 12,
    "axes.titlesize": 12.5,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "lines.linewidth": 1.6,
    "lines.markersize": 7,
    "savefig.dpi": 500,
}
FIGSIZE_SCALE = 1.25
PAD_IN = 0.08                                            # padding round an axes' tight bbox, inches
WRAP_AT = 30                                             # v5 wrapped at 46 and an axis label then ran off the
#                                                          canvas and was clipped -- see the note in save()

MULTI_AXES = {"fig2b_splits", "fig2d_mes_human_vs_virtual", "fig5a_dig_arms"}
# fig2c is two axes (image + colour bar) but one panel: it is copied through whole.

# ------------------------------------------------------------------ (b) the R/ggplot rasters
# xbands/ybands are explicit (start, end) pixel ranges -- not cut lines, so a band can skip a strip that belongs
# to neither neighbour. Cells are read left-to-right, top-to-bottom and numbered p1, p2, ...
# left_margin  (x0, x1)  strip holding the y-axis title + tick labels that a facet row shares; pasted to the
#                        left of every panel that is not in the first column, cropped to that panel's rows.
# top / bottom  [(y0, y1), ...]  full-width bands (shared legend, single shared x-axis title) that are tight-
#                        cropped and stacked above / below every panel of the figure.
# legend_cell  (row, col)  a grid cell that holds only a legend: not emitted as a panel, stacked under each one.
# skip         [(row, col), ...]  cells that hold no result.
RASTER = {
    "fig1b_rna_vs_dna.png": dict(
        size=(3000, 1260),
        xbands=[(0, 1087), (1087, 2038), (2038, 3000)],   # panels 145-1077, 1097-2028, 2048-2980
        ybands=[(150, 1188)],                             # below the legend, above the shared x-axis title
        left_margin=(0, 145),
        top=[(20, 132)], bottom=[(1190, 1252)],
    ),
    "fig1c_share_labellings.png": dict(
        size=(3000, 1500),
        xbands=[(0, 1650), (1650, 3000)],                 # the two facets abut; the seam is the 3 px gap at 1648
        ybands=[(20, 1430)],
        left_margin=(0, 315),                             # the twenty trio names + the y-axis title
        top=[], bottom=[(1432, 1492)],
    ),
    "fig_a_nullA.png": dict(
        size=(3240, 1740),
        # the left panel's last x tick "0.16" runs to x = 1072, past its plot area (ends 1050); the middle
        # panel's y-axis title starts at 1088. Cutting at 1080 keeps the tick whole and leaves no stray glyph.
        xbands=[(0, 1080), (1080, 2124), (2124, 3240)],
        ybands=[(0, 868), (868, 1740)],
        legend_cell=(1, 2),
    ),
    "fig_b_similarity.png": dict(
        size=(3240, 1140),
        xbands=[(0, 1685), (1685, 3240)],
        ybands=[(0, 1140)],
    ),
    "fig_xval_dosage.png": dict(
        size=(3450, 1150),
        xbands=[(0, 1153), (1153, 2303), (2303, 3450)],
        ybands=[(0, 1150)],
    ),
    "fig_xval_genotype.png": dict(
        size=(2400, 1150),
        xbands=[(0, 1204), (1204, 2400)],
        ybands=[(0, 1150)],
    ),
    "fig_xval_meth.png": dict(
        size=(2400, 1150),
        xbands=[(0, 1201), (1201, 2400)],
        ybands=[(0, 1150)],
    ),
    "fig_curve_mouse_validation.png": dict(
        size=(2400, 1500),
        xbands=[(0, 857), (857, 1618), (1618, 2400)],
        # each row band stops above its x-axis title, which is centred under the MIDDLE column only; the two
        # rows carry the same title, so one copy is re-attached under every panel
        ybands=[(195, 968), (1017, 1440)],
        left_margin=(0, 105),
        top=[(108, 190)],                                 # the Brain / Fat / Spleen key; the title band is dropped
        bottom=[(1440, 1482)],
    ),
    "fig6_rats_needed.png": dict(
        size=(2160, 2370),
        xbands=[(0, 722), (722, 1459), (1459, 2160)],
        ybands=[(0, 848), (848, 1604), (1604, 2370)],
        skip=[(1, 1), (1, 2)],                            # "T4 not computed", "T5 not computed"
    ),
}

# single-panel rasters: copied through untouched
COPY_THROUGH = [
    "fig4a_amplitude_dilution.png", "fig4b_dmp_counts.png", "fig4c_mgmt_stp27.png",
    "fig5c_fcpg_w.png", "fig5d_s1_detection.png",
    "fig_c_states.png", "fig_c1_N2_over_core_bins.png", "fig_c1_events.png",
    "fig_curve_dna_vs_rna.png",
]


# ------------------------------------------------------------------------------------------------ helpers
def load_manuscript_module():
    path = DRAFT / "make_figures_rna.py"
    if not path.exists():
        raise SystemExit(f"missing {path}: the draft figure folder is untracked; restore it before building v6")
    spec = importlib.util.spec_from_file_location("mfr_v6", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["mfr_v6"] = mod
    spec.loader.exec_module(mod)
    return mod


def ink(arr: np.ndarray) -> np.ndarray:
    """True where a pixel carries anything but paper white."""
    return (arr[:, :, :3] < 250).any(axis=2)


def tight(arr: np.ndarray, y0: int, y1: int, x0: int = 0, x1: int | None = None, pad: int = 10):
    """Tight bounding box of the ink inside a band, padded; None when the band is blank."""
    x1 = arr.shape[1] if x1 is None else x1
    m = ink(arr[y0:y1, x0:x1])
    if not m.any():
        return None
    rr = np.where(m.any(axis=1))[0]
    cc = np.where(m.any(axis=0))[0]
    return (max(x0, x0 + int(cc[0]) - pad), max(y0, y0 + int(rr[0]) - pad),
            min(x1, x0 + int(cc[-1]) + 1 + pad), min(y1, y0 + int(rr[-1]) + 1 + pad))


def stack(blocks, gap: int = 12) -> Image.Image:
    """Vertically stack images, each centred, on white."""
    blocks = [b for b in blocks if b is not None]
    w = max(b.width for b in blocks)
    h = sum(b.height for b in blocks) + gap * (len(blocks) - 1)
    out = Image.new("RGB", (w, h), (255, 255, 255))
    y = 0
    for b in blocks:
        out.paste(b, ((w - b.width) // 2, y))
        y += b.height + gap
    return out


def beside(left: Image.Image | None, right: Image.Image) -> Image.Image:
    if left is None:
        return right
    h = max(left.height, right.height)
    out = Image.new("RGB", (left.width + right.width, h), (255, 255, 255))
    out.paste(left, (0, 0))
    out.paste(right, (left.width, 0))
    return out


def write(img: Image.Image, path: Path, log: list):
    if max(img.size) < MIN_LONG_SIDE:
        s = UPSCALE_TO / max(img.size)
        img = img.resize((round(img.width * s), round(img.height * s)), Image.LANCZOS)
        note = f"upscaled x{s:.2f}"
    else:
        note = ""
    img.save(path)
    log.append((path.name, img.size, note))


# ------------------------------------------------------------------------------------ (a) matplotlib panels
def build_matplotlib(log: list, notes: dict):
    m = load_manuscript_module()
    m.plt.rcParams.update(SLIDE_RC)
    orig_subplots = m.plt.subplots

    def subplots(*a, **k):
        if k.get("figsize"):
            k["figsize"] = tuple(v * FIGSIZE_SCALE for v in k["figsize"])
        return orig_subplots(*a, **k)

    m.plt.subplots = subplots

    def wrap_long_labels(fig):
        """At slide type an axis label can outrun its axes; wrap it over two lines. The words are unchanged."""
        for ax in fig.axes:
            for getter, setter in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel)):
                t = getter()
                if len(t) > WRAP_AT and "\n" not in t:
                    words = t.split(" ")
                    half, n, best, cut = len(t) // 2, 0, 0, len(words)
                    for i, w in enumerate(words[:-1]):
                        n += len(w) + 1
                        if abs(n - half) < abs(best - half):
                            best, cut = n, i + 1
                    setter(" ".join(words[:cut]) + "\n" + " ".join(words[cut:]))

    def save(fig, stem):
        # v5 wrapped an axis label only past 46 characters. At 12 pt a 36-character label is wider than its
        # axes, runs off the canvas and is drawn clipped -- fig2b lost "...|LFC|)" and "...RNA share|" in v5.
        # Wrapping earlier, then re-running the layout the two-column figures already use, fixes it.
        wrap_long_labels(fig)
        if stem in MULTI_AXES:
            fig.tight_layout()
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        whole = OUT / f"{stem}.png"
        fig.savefig(whole, bbox_inches="tight")
        log.append((whole.name, Image.open(whole).size, "whole figure" if stem in MULTI_AXES else ""))
        if stem in MULTI_AXES:
            axes = [a for a in fig.axes if a.get_visible()]
            axes.sort(key=lambda a: (-round(a.get_position().y1, 3), round(a.get_position().x0, 3)))
            for i, ax in enumerate(axes, 1):
                bb = ax.get_tightbbox(r).transformed(fig.dpi_scale_trans.inverted())
                # deliberately NOT clipped to the canvas: savefig re-renders the requested region, so a label
                # that overhangs the figure edge comes back instead of being cropped away
                bb = Bbox.from_extents(bb.x0 - PAD_IN, bb.y0 - PAD_IN, bb.x1 + PAD_IN, bb.y1 + PAD_IN)
                q = OUT / f"{stem}_p{i}.png"
                fig.savefig(q, bbox_inches=bb)
                lost = []
                if not ax.get_ylabel():
                    lost.append("no y-axis title of its own")
                if not any(t.get_text() and t.get_visible() for t in ax.get_yticklabels()):
                    lost.append("no y tick labels (shared with the sibling panel)")
                if not ax.get_xlabel():
                    lost.append("no x-axis title of its own")
                notes[q.name] = "; ".join(lost)
                log.append((q.name, Image.open(q).size, notes[q.name]))
        m.plt.close(fig)

    m.save = save
    for fn in (m.fig1a, m.fig1d, m.fig2a, m.fig2b, m.fig2c, m.fig2d, m.fig5a, m.fig5e):
        fn()


# ---------------------------------------------------------------------------------------- (b) raster panels
def build_rasters(log: list, problems: list):
    for name, spec in RASTER.items():
        src = V5 / name
        if not src.exists():
            problems.append(f"{name}: not in figures_cns/v5")
            continue
        im = Image.open(src).convert("RGB")
        if im.size != tuple(spec["size"]):
            problems.append(f"{name}: is {im.size}, the cut lines were measured on {tuple(spec['size'])} -- skipped")
            continue
        arr = np.asarray(im)
        stem = Path(name).stem
        xc, yc = spec["xcuts"], spec["ycuts"]
        lm = spec.get("left_margin")
        skip = set(map(tuple, spec.get("skip", [])))
        legend_cell = spec.get("legend_cell")

        tops = [im.crop(t) for t in (tight(arr, *b) for b in spec.get("top", [])) if t]
        bots = [im.crop(t) for t in (tight(arr, *b) for b in spec.get("bottom", [])) if t]

        if legend_cell:
            r, c = legend_cell
            box = tight(arr, yc[r], yc[r + 1], xc[c], xc[c + 1])
            if box:
                leg = im.crop(box)
                bots = bots + [leg]
                write(leg, OUT / f"{stem}_legend.png", log)
            else:
                problems.append(f"{name}: the cell named as the legend block is blank")

        n = 0
        for r in range(len(yc) - 1):
            for c in range(len(xc) - 1):
                n += 1
                if (r, c) in skip or (legend_cell and (r, c) == tuple(legend_cell)):
                    continue
                cell = im.crop((xc[c], yc[r], xc[c + 1], yc[r + 1]))
                margin = None
                if lm and c > 0:
                    margin = im.crop((lm[0], yc[r], lm[1], yc[r + 1]))
                write(stack(tops + [beside(margin, cell)] + bots), OUT / f"{stem}_p{n}.png", log)


def main() -> int:
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)
    log, problems, notes = [], [], {}

    build_matplotlib(log, notes)
    build_rasters(log, problems)

    for name in COPY_THROUGH:
        s = V5 / name
        if not s.exists():
            problems.append(f"{name}: not in figures_cns/v5")
            continue
        shutil.copyfile(s, OUT / name)
        log.append((name, Image.open(OUT / name).size, "single panel, copied through"))

    print(f"\nwrote {len(log)} PNGs to {OUT}\n")
    for nm, size, note in sorted(log):
        flag = "  <-- UNDER 1200 px" if max(size) < MIN_LONG_SIDE else ""
        print(f"  {nm:44s} {size[0]:5d} x {size[1]:5d}  {note}{flag}")
    if problems:
        print("\nproblems:")
        for p in problems:
            print("   ", p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
