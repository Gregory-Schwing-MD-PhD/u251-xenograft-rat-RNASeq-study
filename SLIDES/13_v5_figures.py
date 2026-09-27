# -*- coding: utf-8 -*-
"""Slide-size versions of the manuscript's RNA panels for the CNS 2026 deck, version 5.

The manuscript figures are drawn by MBR/nature_draft/figures/make_figures_rna.py at 8 pt from the frozen grid
tables. A slide needs the same panel with type that survives projection ([[chart-preferences]]: chart type 17-30 pt,
one manuscript panel per slide), so this script IMPORTS that module, raises its font constants and re-saves each
panel into SLIDES/figures_cns/v5/. It draws nothing of its own and computes no number: the panels are the
manuscript's, at a different size.

The DNA and contralateral panels were drawn on the grid in R at 2,100-3,240 px with large type already; they are
copied, not redrawn (fetch_grid_outputs.sh / MBR/nature_draft/figures/grid_*).

    python SLIDES/13_v5_figures.py
"""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DRAFT = ROOT / "MBR" / "nature_draft" / "figures"      # untracked; the deck build needs it, as it needs MBR/ESM_1.xlsx
OUT = HERE / "figures_cns" / "v5"

# type as drawn; each panel is placed 1.7-1.9x its drawn width on the slide, so 12 pt reads as about 21 pt
SLIDE_RC = {
    "font.size": 12,
    "axes.titlesize": 12.5,
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11,
    "lines.linewidth": 1.6,
    "lines.markersize": 7,
    "savefig.dpi": 200,
}
# the manuscript panels are 2.1-2.8 in tall; widen a little so the larger type does not collide
FIGSIZE_SCALE = 1.25


def load_manuscript_module():
    path = DRAFT / "make_figures_rna.py"
    if not path.exists():
        raise SystemExit(f"missing {path}: the draft figure folder is untracked; restore it before building v5")
    spec = importlib.util.spec_from_file_location("mfr_v5", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["mfr_v5"] = mod
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    m = load_manuscript_module()
    m.plt.rcParams.update(SLIDE_RC)

    orig_subplots = m.plt.subplots

    def subplots(*a, **k):                       # same layout, a little more canvas for the larger type
        if "figsize" in k and k["figsize"]:
            k["figsize"] = tuple(v * FIGSIZE_SCALE for v in k["figsize"])
        return orig_subplots(*a, **k)

    m.plt.subplots = subplots
    written = []

    def wrap_long_labels(fig):
        """At slide type an axis label can outrun its axes. Wrap it over two lines; the words are unchanged."""
        for ax in fig.axes:
            for getter, setter in ((ax.get_xlabel, ax.set_xlabel), (ax.get_ylabel, ax.set_ylabel)):
                t = getter()
                if len(t) > 46 and "\n" not in t:
                    words = t.split(" ")
                    half = len(t) // 2
                    n, best = 0, 0
                    for i, w in enumerate(words[:-1]):
                        n += len(w) + 1
                        if abs(n - half) < abs(best - half):
                            best, cut = n, i + 1
                    setter(" ".join(words[:cut]) + "\n" + " ".join(words[cut:]))

    def save(fig, stem):                          # PNG only; the deck does not read PDFs
        wrap_long_labels(fig)
        p = OUT / f"{stem}.png"
        fig.savefig(p, bbox_inches="tight")
        m.plt.close(fig)
        written.append(p.name)

    m.save = save
    for fn in (m.fig1a, m.fig1d, m.fig2a, m.fig2b, m.fig2c, m.fig2d, m.fig5a, m.fig5e):
        fn()

    # the R-drawn panels and the power figure: copied at their published size
    copies = {
        "grid_dna/fig1b_rna_vs_dna.png": "fig1b_rna_vs_dna.png",
        "grid_dna/fig1c_share_labellings.png": "fig1c_share_labellings.png",
        "grid_dna/fig4a_amplitude_dilution.png": "fig4a_amplitude_dilution.png",
        "grid_dna/fig4b_dmp_counts.png": "fig4b_dmp_counts.png",
        "grid_dna/fig4c_mgmt_stp27.png": "fig4c_mgmt_stp27.png",
        "grid_dna/fig5c_fcpg_w.png": "fig5c_fcpg_w.png",
        "grid_dna/fig5d_s1_detection.png": "fig5d_s1_detection.png",
        "grid_spread/fig_a_nullA.png": "fig_a_nullA.png",
        "grid_spread/fig_b_similarity.png": "fig_b_similarity.png",
        "grid_spread/fig_c_states.png": "fig_c_states.png",
        "grid_spread_dna/fig_c1_N2_over_core_bins.png": "fig_c1_N2_over_core_bins.png",
        "grid_spread_dna/fig_c1_events.png": "fig_c1_events.png",
        "grid_xval/fig_xval_dosage.png": "fig_xval_dosage.png",
        "grid_xval/fig_xval_genotype.png": "fig_xval_genotype.png",
        "grid_xval/fig_xval_meth.png": "fig_xval_meth.png",
        "grid_curve/fig_a_mouse_validation.png": "fig_curve_mouse_validation.png",
        "grid_curve/fig_c_dna_vs_rna.png": "fig_curve_dna_vs_rna.png",
        "fig6_rats_needed.png": "fig6_rats_needed.png",
    }
    for src, dst in copies.items():
        s = DRAFT / src
        if not s.exists():
            print(f"  MISSING (skipped): {src}")
            continue
        shutil.copyfile(s, OUT / dst)
        written.append(dst)

    print(f"wrote {len(written)} figures to {OUT}")
    for n in sorted(written):
        print("   ", n)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
