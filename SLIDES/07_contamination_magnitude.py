# -*- coding: utf-8 -*-
"""How much of the translational fall could rat contamination or tumour purity explain? Two estimates, for the
translation-initiation, elongation and ribosome sets, written to figures_cns/contamination_magnitude.json.

1. Mixture bound. The human stream of a tumour library is graft reads plus the 'both' reads (main.nf concatenates
   *graft* and *both*). Take, as an upper bound, every 'both' read of a tumour as rat-derived with the expression profile
   of the rat-brain controls' human stream (their human-stream CPM). Then a gene's expected CPM in sample i is
   (1 - m_i) * T_g + m_i * C_g, m_i = both_i / (graft_i + both_i), T_g the gene's tumour CPM, C_g its control CPM. DESeq2's
   median-of-ratios normalisation divides by the typical gene, so the predicted log2 shift of gene g between arms is the
   mean over recurrent minus mean over primary of log2(((1 - m) + m r_g) / ((1 - m) + m r_med)), r = C_g / T_g.
2. Matched background. For every member of a set, the 25 tested genes outside all translation/ribosome sets nearest in
   (log tumour CPM, log control/tumour ratio); the set's mean log2FC minus the matched genes' mean log2FC is what the set
   does beyond its abundance and its rat-brain share. Unadjusted RUVSeq-run DESeq2 (the run the contamination slide
   uses); the manuscript's shrunk fold changes are smaller throughout.

    python SLIDES/07_contamination_magnitude.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from u251_paths import COUNTS, DECONTAM_DIR, GSEA_DIR, METADATA, RUVSEQ_DIR, check, gsea_file  # noqa: E402,F401
check()
CTRL = ["IL64B", "N168B", "N269B"]; PRI = ["IL67B", "IL68B", "IL69B"]; REC = ["IL66B", "NL70B", "NL71B"]
SETS = {"initiation": "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "elongation": "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "ribosome": "KEGG_RIBOSOME"}
ALLTRANS = ["KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "KEGG_RIBOSOME",
            "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "REACTOME_SELENOAMINO_ACID_METABOLISM", "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
            "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION"]

cnt = pd.read_csv(COUNTS, sep="\t").set_index("gene_id")
X = cnt[CTRL + PRI + REC].astype(float)
cpm = X / X.sum(0) * 1e6
T = cpm[PRI + REC].mean(1); C = cpm[CTRL].mean(1)
bias = pd.read_csv(DECONTAM_DIR / "contamination_bias_per_gene.csv").set_index("gene_id")
tested = bias.index
lfc = bias["log2FoldChange"]
base = pd.read_csv(RUVSEQ_DIR / "ruvseq_baseline_de.tsv", sep="\t").set_index("gene_id")
sym = base["symbol"].astype(str)
meta = pd.read_csv(METADATA).set_index("sample")
m = (meta["both_pct"] / (meta["graft_pct"] + meta["both_pct"]))
r = (C / T.replace(0, np.nan)).reindex(tested)
r_med = float(r.median())


def members(key):
    s = set(pd.read_csv(gsea_file(f"{key}.tsv"), sep="\t")["SYMBOL"].astype(str))
    return [g for g in tested if sym.get(g) in s]


trans_all = set().union(*[set(members(k)) for k in ALLTRANS if (gsea_file(f"{k}.tsv")).exists()])
bg = [g for g in tested if g not in trans_all and np.isfinite(r.get(g, np.nan)) and T.get(g, 0) > 0]
feat = lambda ids: np.column_stack([np.log10(T.reindex(ids).to_numpy() + 1), np.log10(r.reindex(ids).to_numpy() + 1e-3)])  # noqa: E731
Fbg = feat(bg); lbg = lfc.reindex(bg).to_numpy()
out = {"m_upper_bound": {s: round(float(m[s]), 4) for s in PRI + REC}, "median_ratio_tested": r_med, "sets": {}}
print("m (upper bound on the rat share of each tumour's human stream):", {s: round(float(m[s]), 3) for s in PRI + REC})
for name, key in SETS.items():
    ids = [g for g in members(key) if np.isfinite(r.get(g, np.nan))]
    rg = r.reindex(ids).to_numpy()
    pred = []
    for g_r in rg:
        f = lambda s: np.log2(((1 - m[s]) + m[s] * g_r) / ((1 - m[s]) + m[s] * r_med))  # noqa: E731
        pred.append(np.mean([f(s) for s in REC]) - np.mean([f(s) for s in PRI]))
    Fs = feat(ids)
    matched = []
    for row in Fs:
        d = np.sqrt(((Fbg - row) ** 2).sum(1))
        nn = np.argsort(d)[:25]
        matched.append(lbg[nn].mean())
    obs = float(lfc.reindex(ids).mean())
    rec = {"n": len(ids), "median_ratio": float(np.median(rg)), "observed_mean_lfc": obs, "mixture_bound_mean_lfc": float(np.mean(pred)),
           "matched_background_mean_lfc": float(np.mean(matched)), "excess_over_matched": obs - float(np.mean(matched))}
    out["sets"][name] = rec
    print(f"{name:11s} n {len(ids):3d}  ratio {rec['median_ratio']:.3f} (median gene {r_med:.3f})  observed {obs:+.3f}  mixture bound {rec['mixture_bound_mean_lfc']:+.4f}  "
          f"matched background {rec['matched_background_mean_lfc']:+.3f}  excess {rec['excess_over_matched']:+.3f}")
# ---------------------------------------------------------------------- leave one tumour out
# Fold changes recomputed from median-of-ratios normalised counts of the six tumours (genes with >= 10 reads in >= 3
# tumours), so a tumour can be dropped; with all six they reproduce the DESeq2 set means to 0.003. One primary, IL68B,
# is the highest of the six on nearly every translation gene; this is the check the adversarial review asked for.
X6 = cnt[PRI + REC].astype(float)
X6 = X6[(X6 >= 10).sum(1) >= 3]
lg6 = np.log(X6.replace(0, np.nan)); geo6 = lg6.mean(1); ok6 = geo6.notna()
sf6 = np.exp(lg6[ok6].sub(geo6[ok6], axis=0).median(0))
N6 = X6 / sf6
lfc6 = lambda pri, rec, ids: np.log2((N6.loc[ids, rec].mean(1) + 0.5) / (N6.loc[ids, pri].mean(1) + 0.5))  # noqa: E731
t6 = [g for g in tested if g in N6.index]
bg6 = [g for g in bg if g in N6.index]
Fbg6 = feat(bg6)
out["leave_one_out"] = {"size_factors": {k: float(v) for k, v in sf6.items()}, "sets": {}}
for name, key in SETS.items():
    ids = [g for g in members(key) if g in N6.index and np.isfinite(r.get(g, np.nan))]
    rec_ = {"all_six": float(lfc6(PRI, REC, ids).mean())}
    for s in PRI + REC:
        p_ = [q for q in PRI if q != s]; r_ = [q for q in REC if q != s]
        rec_["without_" + s] = float(lfc6(p_, r_, ids).mean())
    L6 = np.log2(N6.loc[ids] + 0.5)
    rec_["highest_tumour_counts"] = {k: int(v) for k, v in L6.idxmax(1).value_counts().items()}
    rec_["n"] = len(ids)
    for lab, pri in (("all_six", PRI), ("without_IL68B", [q for q in PRI if q != "IL68B"])):
        lbg = lfc6(pri, REC, bg6).to_numpy()
        rec_["matched_" + lab] = float(np.mean([lbg[np.argsort(np.sqrt(((Fbg6 - row) ** 2).sum(1)))[:25]].mean() for row in feat(ids)]))
    out["leave_one_out"]["sets"][name] = rec_
    print(f"{name:11s} all six {rec_['all_six']:+.3f}  without IL68B {rec_['without_IL68B']:+.3f} (matched {rec_['matched_without_IL68B']:+.3f})  "
          f"without others {min(v for k, v in rec_.items() if k.startswith('without_') and k != 'without_IL68B'):+.3f} to "
          f"{max(v for k, v in rec_.items() if k.startswith('without_') and k != 'without_IL68B'):+.3f}  highest {rec_['highest_tumour_counts']}")
out["leave_one_out"]["all_tested_without_IL68B"] = float(lfc6([q for q in PRI if q != "IL68B"], REC, t6).mean())
json.dump(out, open(HERE / "figures_cns" / "contamination_magnitude.json", "w"), indent=1)
print("wrote contamination_magnitude.json")


# ---------------------------------------------------------------------- the slide chart: magnitude, set by set
import importlib.util  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

spec = importlib.util.spec_from_file_location("mk04", HERE / "04_make_cns_figures.py")
mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
cc = json.load(open(HERE / "figures_cns" / "contamination_check.json", encoding="utf-8"))["sets"]
LOO = out["leave_one_out"]["sets"]
bars = [("observed", mk.NAVY, lambda nm: out["sets"][nm]["observed_mean_lfc"]),
        ("without one primary\ntumor (IL68B)", mk.RED, lambda nm: LOO[nm]["without_IL68B"]),
        ("adjusted for control-\nderived factors (k = 2)", mk.BLUE, lambda nm: cc[SETS[nm]]["k2"]["mean_lfc"]),
        ("genes matched on abundance\nand control/tumor ratio", "#A9B4CC", lambda nm: out["sets"][nm]["matched_background_mean_lfc"]),
        ("the most rat reads\ncould produce", mk.GOLD, lambda nm: out["sets"][nm]["mixture_bound_mean_lfc"])]
names = [("initiation", "translation\ninitiation"), ("elongation", "translation\nelongation"), ("ribosome", "ribosome")]
fig, ax = plt.subplots(figsize=(9.6, 5.0))
wbar, pitch = 0.8, 6.0
for gi, (nm, lab) in enumerate(names):
    for bi, (_, col, f) in enumerate(bars):
        v = f(nm); x = gi * pitch + bi
        ax.bar(x, v, width=wbar, color=col, zorder=3)
        ax.text(x, v - 0.012, mk.MINUS(f"{v:+.2f}"), ha="center", va="top", fontsize=13, fontweight="bold" if bi == 0 else "normal", color=mk.INK)
    ax.text(gi * pitch + 2.0, 0.025, f"{lab}\n{out['sets'][nm]['n']} genes", ha="center", va="bottom", fontsize=13.5,
            fontweight="bold", color=mk.INK, linespacing=1.15)
ax.axhline(0, color=mk.INK, lw=1.2, zorder=4)
ax.set_xticks([]); ax.set_xlim(-0.7, 2 * pitch + 4.7); ax.set_ylim(-0.5, 0.17)
ax.set_yticks([0, -0.1, -0.2, -0.3, -0.4]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: mk.MINUS(f"{v:.1f}")))
ax.set_ylabel("mean log$_2$ fold change of the set", fontsize=14, color=mk.INK)
for sd in ("top", "right", "bottom"):
    ax.spines[sd].set_visible(False)
ax.grid(axis="y", color=mk.GRID, lw=0.8, zorder=0); ax.set_axisbelow(True); ax.tick_params(axis="y", labelsize=13, length=0)
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c, _ in bars]
ax.legend(handles, [b[0] for b in bars], loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=11.5, frameon=False,
          handlelength=1.2, columnspacing=1.2)
mk.save(fig, "chart_contamination_magnitude.png")
