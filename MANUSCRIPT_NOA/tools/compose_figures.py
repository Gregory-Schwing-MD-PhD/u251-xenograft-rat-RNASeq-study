#!/usr/bin/env python
"""Compose the manuscript's six display items from the single-panel figures, with panel letters.

Neuro-Oncology Advances: composite figures grouped under one number, subparts labelled A, B, C in the upper left on
the face of the illustration, no more than 8 panels, 600 dpi for bitmaps, and embedded text large enough to survive
reduction. Panels are never redrawn here: each is a file produced by the analysis that owns it, pasted at its own
resolution onto a white canvas.

    python tools/compose_figures.py            # write figures/fig1..fig6 (+ .json manifest per figure)
    python tools/compose_figures.py --check    # verify each figure is newer than its panels; exit 1 if stale
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
PANELS = [REPO / "SLIDES" / "figures_cns" / "v6", REPO / "SLIDES" / "figures_cns" / "v5",
          REPO / "MBR" / "nature_draft" / "figures"]
OUT = ROOT / "figures"
DPI = 600
GAP = 26
LETTER = 64          # px height of the panel letter band

# figure -> rows of panels. Each panel: (file, caption-letter). A row is laid out left to right and scaled to a
# common height; rows are stacked and scaled to a common width.
FIGURES = {
    "fig1_composition": [
        ["fig1a_read_composition.png", "fig1d_lesion_accounting.png"],
        ["fig1b_rna_vs_dna_p1.png", "fig1c_share_labellings_p1.png"],
    ],
    "fig2_decomposition": [
        ["fig2a_shapley_de_genes.png"],
        ["fig2b_splits_p1.png", "fig2b_splits_p2.png"],
    ],
    "fig3_programmes": [
        ["fig2c_categories_by_view.png", "fig2d_mes_human_vs_virtual_p1.png"],
        ["fig2d_mes_human_vs_virtual_p2.png", "fig5a_dig_arms_p2.png"],
    ],
    "fig4_methylation": [
        ["fig4b_dmp_counts.png"],
        ["fig4c_mgmt_stp27.png", "fig4a_amplitude_dilution.png"],
    ],
    "fig5_contralateral": [
        ["fig_c1_N2_over_core_bins.png"],
        ["fig_a_nullA_p1.png", "fig_b_similarity_p2.png", "fig_c1_events.png"],
    ],
    "fig6_rats_needed": [
        ["fig6_rats_needed_p1.png", "fig6_rats_needed_p2.png", "fig6_rats_needed_p3.png"],
        ["fig6_rats_needed_p4.png", "fig6_rats_needed_p8.png", "fig6_rats_needed_p9.png"],
    ],
}


def find(name: str) -> Path:
    for d in PANELS:
        p = d / name
        if p.exists():
            return p
    raise SystemExit(f"panel not found in {[str(d) for d in PANELS]}: {name}")


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def font(size: int):
    for cand in ("arialbd.ttf", "Arialbd.ttf", "DejaVuSans-Bold.ttf"):
        try:
            return ImageFont.truetype(cand, size)
        except OSError:
            continue
    return ImageFont.load_default()


def compose(name: str, rows: list[list[str]]) -> dict:
    letters = iter("ABCDEFGH")
    loaded = [[(Image.open(find(n)).convert("RGB"), n, next(letters)) for n in row] for row in rows]
    # scale each row to a common height, then all rows to a common width
    row_imgs = []
    for row in loaded:
        h = min(im.height for im, _, _ in row)
        scaled = [(im.resize((max(1, round(im.width * h / im.height)), h), Image.LANCZOS), n, L) for im, n, L in row]
        w = sum(im.width for im, _, _ in scaled) + GAP * (len(scaled) - 1)
        canvas = Image.new("RGB", (w, h + LETTER), "white")
        d = ImageDraw.Draw(canvas)
        x = 0
        for im, _, L in scaled:
            d.text((x + 6, 4), L, fill=(0, 0, 0), font=font(LETTER - 12))
            canvas.paste(im, (x, LETTER))
            x += im.width + GAP
        row_imgs.append(canvas)
    W = max(im.width for im in row_imgs)
    row_imgs = [im if im.width == W else im.resize((W, round(im.height * W / im.width)), Image.LANCZOS)
                for im in row_imgs]
    H = sum(im.height for im in row_imgs) + GAP * (len(row_imgs) - 1)
    fig = Image.new("RGB", (W, H), "white")
    y = 0
    for im in row_imgs:
        fig.paste(im, (0, y))
        y += im.height + GAP
    OUT.mkdir(exist_ok=True)
    p = OUT / f"{name}.png"
    fig.save(p, dpi=(DPI, DPI))
    man = {"figure": name, "px": list(fig.size), "dpi": DPI,
           "panels": [{"letter": L, "file": n, "source": str(find(n).relative_to(REPO)), "sha256": sha256(find(n))}
                      for row in loaded for _, n, L in row]}
    (OUT / f"{name}.json").write_text(json.dumps(man, indent=1), encoding="utf-8", newline="\n")
    return man


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    stale = []
    for name, rows in FIGURES.items():
        if a.check:
            p, m = OUT / f"{name}.png", OUT / f"{name}.json"
            if not p.exists() or not m.exists():
                stale.append((name, "not built"))
                continue
            old = {x["file"]: x["sha256"] for x in json.load(open(m, encoding="utf-8"))["panels"]}
            for row in rows:
                for n in row:
                    if old.get(n) != sha256(find(n)):
                        stale.append((name, f"panel changed: {n}"))
            continue
        man = compose(name, rows)
        print(f"{name}: {len(man['panels'])} panels, {man['px'][0]}x{man['px'][1]} px at {DPI} dpi")
    for name, why in stale:
        print(f"   STALE {name}: {why}")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())
