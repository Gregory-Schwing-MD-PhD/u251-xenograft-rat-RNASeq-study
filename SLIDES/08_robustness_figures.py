# -*- coding: utf-8 -*-
"""Robustness and confound charts for the v4 talk, read from the analysis outputs (no number is typed here).

  chart_graft_scatter.png   graft % (xengsort human share) against the translation-initiation and Neftel AC GSVA
                            scores, one dot per tumour (ANALYSIS/graft_relation/results/per_sample.tsv, score_correlations.tsv)
  chart_holdout.png         three findings across the full pipeline re-runs: all six tumours, IL68B held out, IL66B held
                            out (ANALYSIS/holdout_IL68B|IL66B/comparison.json): translation-initiation FDR q, ciclopirox
                            final rank, Neftel AC p
  chart_states.png          CIBERSORT v1.04 Neftel fractions per tumour and the culture (ANALYSIS/cibersort/results/v104)
  robustness.json           every number the slides quote from these analyses

    python SLIDES/08_robustness_figures.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
A = HERE.parent / "ANALYSIS"
OUT = HERE / "figures_cns"
NAVY, BLUE, RED, GOLD, GREY, GRID = "#191D63", "#4C70B7", "#EF4627", "#E0A800", "#8A8A8A", "#DDDDDD"
plt.rcParams.update({"font.family": "Arial", "font.size": 15, "axes.edgecolor": "#BFBFBF", "axes.linewidth": 0.9,
                     "xtick.color": GREY, "ytick.color": GREY, "figure.dpi": 220, "savefig.dpi": 220})
PRI, REC = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]
TI = "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"


def bare(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def save(fig, name):
    fig.savefig(OUT / name, facecolor="white", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print("  ", name)


facts = {}

# ---- graft % against the two scores
ps = pd.read_csv(A / "graft_relation/results/per_sample.tsv", sep="\t").set_index("sample")
cr = pd.read_csv(A / "graft_relation/results/score_correlations.tsv", sep="\t").set_index("score")
de = pd.read_csv(A / "graft_relation/results/de_models.tsv", sep="\t")
tl = pd.read_csv(A / "graft_relation/results/translation_sets_lfc_by_model.tsv", sep="\t").set_index("set")
facts["graft"] = {
    "graft_pct": ps.graft_pct.round(1).to_dict(),
    "r_group_graft": float(np.corrcoef((ps.group == "Recurrent").astype(int), ps.graft_pct)[0, 1]),
    "r_translation_initiation": float(cr.loc[TI, "r_graft"]), "p_translation_initiation": float(cr.loc[TI, "p_graft"]),
    "r_AC": float(cr.loc["Neftel_AC", "r_graft"]), "p_AC": float(cr.loc["Neftel_AC", "p_graft"]),
    "r_MES1": float(cr.loc["Neftel_MES1", "r_graft"]), "p_MES1": float(cr.loc["Neftel_MES1", "p_graft"]),
    "AC_group_p": float(cr.loc["Neftel_AC", "group_p"]), "AC_group_p_adj": float(cr.loc["Neftel_AC", "group_p_adj"]),
    "TI_group_p": float(cr.loc[TI, "group_p"]), "TI_group_p_adj": float(cr.loc[TI, "group_p_adj"]),
    "de_group": int(de.loc[de.model == "~ group", "padj05"].iloc[0]),
    "de_group_adj": int(de.loc[de.model.str.startswith("~ graft_z + group (group"), "padj05"].iloc[0]),
    "TI_lfc_group": float(tl.loc[TI, "mean_lfc_group"]), "TI_lfc_group_adj": float(tl.loc[TI, "mean_lfc_group_adj"]),
    "TI_effect": float(cr.loc[TI, "group_effect"]), "TI_effect_adj": float(cr.loc[TI, "group_effect_adj_graft"]),
    "AC_effect": float(cr.loc["Neftel_AC", "group_effect"]), "AC_effect_adj": float(cr.loc["Neftel_AC", "group_effect_adj_graft"]),
    "pair": {s: {"graft": float(ps.loc[s, "graft_pct"]), "TI": float(ps.loc[s, TI]), "AC": float(ps.loc[s, "Neftel_AC"])} for s in ("IL67B", "NL70B")},
    "assigned_pct_min": float((ps.graft_pct + ps.host_pct).min()), "assigned_pct_max": float((ps.graft_pct + ps.host_pct).max()),
    "size_factor_r": float(np.corrcoef(ps.size_factor, ps.graft_pct)[0, 1]),
    "scores_input": "GSVA on DESeq2 vst counts (graft_relation.R)",
}

# ---- v4 slide 7: the translation panel alone, the graft-matched pair ringed
fig, ax = plt.subplots(figsize=(7.4, 4.6))
for grp, c, mk, lab in ((PRI, NAVY, "o", "primary"), (REC, RED, "s", "recurrent")):
    ax.scatter(ps.loc[grp, "graft_pct"], ps.loc[grp, TI], s=190, color=c, marker=mk, edgecolor="white", lw=1.5, zorder=3, label=lab)
    for s in grp:
        off = {"IL67B": (10, -4), "NL70B": (10, -4), "IL66B": (10, -4), "NL71B": (10, 4)}.get(s, (10, 4))
        ax.annotate(s, (ps.loc[s, "graft_pct"], ps.loc[s, TI]), xytext=off, textcoords="offset points", fontsize=13, color=GREY)
pr = ps.loc[["IL67B", "NL70B"]]
cx, cy = pr.graft_pct.mean(), pr[TI].mean()
ax.add_patch(matplotlib.patches.Ellipse((cx, cy), 6.0, abs(pr[TI].diff().iloc[-1]) + 0.35, fill=False, ec=GOLD, lw=2.2, ls="--", zorder=2))
ax.annotate("same human share", (cx, cy + (abs(pr[TI].diff().iloc[-1]) + 0.35) / 2), xytext=(0, 8), textcoords="offset points",
            ha="center", fontsize=12.5, color="#9A7300")
ax.set_xlabel("human reads in the library (%)", fontsize=15)
ax.set_ylabel("translation initiation (GSVA)", fontsize=15)
ax.tick_params(labelsize=13)
ax.legend(frameon=False, fontsize=13, loc="upper left")
bare(ax)
fig.tight_layout()
save(fig, "chart_graft_translation.png")

# ---- v4 slide 6: translation-initiation FDR q in every leave-one-tumour-out run (41 runs)
loo = pd.read_csv(A / "gsea_leave_one_out/loo_sets.tsv", sep="\t")
ti = loo[loo.set == TI].copy()
conds = ["full", "drop_IL67B", "drop_IL68B", "drop_IL69B", "drop_IL66B", "drop_NL70B", "drop_NL71B"]
labs = ["all six", "−IL67B", "−IL68B", "−IL69B", "−IL66B", "−NL70B", "−NL71B"]
fig, ax = plt.subplots(figsize=(9.6, 4.5))
for i, cnd in enumerate(conds):
    g = ti[ti.condition == cnd].sort_values("seed")
    col = GREY if cnd == "full" else NAVY if cnd in ("drop_IL67B", "drop_IL68B", "drop_IL69B") else RED
    ax.scatter(i + np.linspace(-0.16, 0.16, len(g)), g.fdr_q, s=80, color=col, edgecolor="white", lw=1.0, zorder=3)
    nf = ti[ti.condition == cnd + "_nofilter"]
    if len(nf):
        ax.scatter([i + 0.3] * len(nf), nf.fdr_q, s=70, facecolor="white", edgecolor=col, lw=1.6, marker="D", zorder=3)
ax.axhline(0.05, color=RED, ls="--", lw=1.4)
ax.text(len(conds) - 0.45, 0.056, "q = 0.05", color=RED, fontsize=12, ha="right", va="bottom")
ax.set_yscale("log"); ax.set_ylim(0.0007, 1.3)
ax.set_xticks(range(len(conds))); ax.set_xticklabels(labs, fontsize=13)
ax.set_ylabel("translation initiation, FDR q", fontsize=14)
ax.tick_params(labelsize=12)
bare(ax)
fig.tight_layout()
save(fig, "chart_loo_q.png")
main = ti[~ti.condition.str.endswith("_nofilter")]
facts["loo"] = {"runs": int(len(ti)), "nes_min": float(ti.NES.min()), "nes_max": float(ti.NES.max()),
                "nom_p_max": float(ti.nom_p.max()), "q_min": float(ti.fdr_q.min()), "q_max": float(ti.fdr_q.max()),
                "n_q05": int((ti.fdr_q < 0.05).sum()),
                "q_full_min": float(ti[ti.condition == "full"].fdr_q.min()), "q_full_max": float(ti[ti.condition == "full"].fdr_q.max()),
                "q_by_condition": {c: [round(float(x), 4) for x in main[main.condition == c].sort_values("seed").fdr_q] for c in conds}}
fig, axs = plt.subplots(1, 2, figsize=(10.8, 4.4))
for ax, col, lab, key in ((axs[0], TI, "translation initiation (GSVA)", "r_translation_initiation"),
                          (axs[1], "Neftel_AC", "astrocyte-like, AC (GSVA)", "r_AC")):
    for grp, c, mk in ((PRI, NAVY, "o"), (REC, RED, "s")):
        ax.scatter(ps.loc[grp, "graft_pct"], ps.loc[grp, col], s=150, color=c, marker=mk, edgecolor="white", lw=1.5, zorder=3)
        for s in grp:
            off = {"IL67B": (8, -14), "NL70B": (8, 6), "IL66B": (8, 4), "NL71B": (8, 6)}.get(s, (8, 5))
            ax.annotate(s, (ps.loc[s, "graft_pct"], ps.loc[s, col]), xytext=off, textcoords="offset points", fontsize=11, color=GREY)
    ax.set_xlabel("human reads in the library (%)", fontsize=14)
    ax.set_ylabel(lab, fontsize=14)
    ax.set_title(f"r = {facts['graft'][key]:+.2f}", fontsize=14, color=NAVY, loc="left")
    bare(ax)
axs[0].scatter([], [], s=90, color=NAVY, marker="o", label="primary")
axs[0].scatter([], [], s=90, color=RED, marker="s", label="recurrent")
axs[0].legend(frameon=False, fontsize=12, loc="upper left")
fig.tight_layout()
save(fig, "chart_graft_scatter.png")

# ---- the three findings across the hold-outs
arms = {}
for h in ("IL68B", "IL66B"):
    d = json.load(open(A / f"holdout_{h}/comparison.json", encoding="utf-8"))
    lead = {g["set"]: g for g in d["gsea"]["lead"]}
    arms[h] = {"q_TI": lead[TI]["q_ho"], "cic_rank": d["ranking"]["ciclopirox"]["holdout"]["rank_bbb"],
               "cic_rank_nobbb": d["ranking"]["ciclopirox"]["holdout"]["rank_nobbb"],
               "AC_p": d["subtypes"]["holdout"]["Neftel_AC"]["p"], "AC_change": d["subtypes"]["holdout"]["Neftel_AC"]["change"],
               "de_n": d["de"]["holdout_n"], "sets_q05": None}
    pub = {"q_TI": lead[TI]["q_pub"], "cic_rank": d["ranking"]["ciclopirox"]["published"]["rank_bbb"],
           "cic_rank_nobbb": d["ranking"]["ciclopirox"]["published"]["rank_nobbb"],
           "AC_p": d["subtypes"]["published"]["Neftel_AC"]["p"], "AC_change": d["subtypes"]["published"]["Neftel_AC"]["change"],
           "de_n": d["de"]["published_n"]}
facts["holdout"] = {"all_six": pub, **arms}
labels = ["all six", "without\nIL68B", "without\nIL66B"]
rows = [pub, arms["IL68B"], arms["IL66B"]]
fig, axs = plt.subplots(1, 3, figsize=(12.6, 4.2))
x = np.arange(3)
axs[0].plot(x, [r["q_TI"] for r in rows], color=NAVY, lw=2.5, marker="o", ms=11)
axs[0].axhline(0.05, color=RED, ls="--", lw=1.4)
axs[0].set_yscale("log"); axs[0].set_ylim(0.005, 1)
axs[0].set_title("translation initiation, FDR q (seed 1234)", fontsize=14, color=NAVY, loc="left")
axs[1].plot(x, [r["cic_rank"] for r in rows], color=NAVY, lw=2.5, marker="o", ms=11)
axs[1].set_ylim(5.5, 0.5); axs[1].set_yticks([1, 2, 3, 4, 5])
axs[1].set_title("ciclopirox, rank among candidates", fontsize=14, color=NAVY, loc="left")
axs[2].plot(x, [r["AC_p"] for r in rows], color=NAVY, lw=2.5, marker="o", ms=11)
axs[2].axhline(0.05, color=RED, ls="--", lw=1.4)
axs[2].set_ylim(0, 0.1)
axs[2].set_title("astrocyte-like score falls, p", fontsize=14, color=NAVY, loc="left")
for ax, key, fmt in ((axs[0], "q_TI", "{:.3f}"), (axs[1], "cic_rank", "{:d}"), (axs[2], "AC_p", "{:.3f}")):
    for i, r in enumerate(rows):
        ax.annotate(fmt.format(r[key]), (i, r[key]), xytext=(10, -16) if i == 1 else (0, 11), textcoords="offset points",
                    ha="left" if i == 1 else "center", fontsize=12, color=NAVY)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=12); ax.set_xlim(-0.4, 2.4)
    bare(ax)
fig.tight_layout()
save(fig, "chart_holdout.png")

# ---- deconvolution: Neftel fractions per tumour
v = pd.read_csv(A / "cibersort/results/v104/fractions_v104_neftel4_confident.tsv", sep="\t").set_index("sample")
facts["deconv"] = {"MES_primary": float(v.loc[PRI, "MES"].mean()), "MES_recurrent": float(v.loc[REC, "MES"].mean()),
                   "AC_primary": float(v.loc[PRI, "AC"].mean()), "AC_recurrent": float(v.loc[REC, "AC"].mean()),
                   "fitR_min": float(v.loc[PRI + REC, "Correlation"].min()), "fitR_max": float(v.loc[PRI + REC, "Correlation"].max()),
                   "MES_C2B": float(v.loc["C2B", "MES"])}
order = PRI + REC + ["C2B"]
fig, ax = plt.subplots(figsize=(7.6, 4.2))
cols = {"MES": NAVY, "AC": BLUE, "NPC": GOLD, "OPC": GREY}
xi = np.arange(len(order))
for j, (st, c) in enumerate(cols.items()):
    ax.scatter(xi + (j - 1.5) * 0.12, v.loc[order, st], color=c, s=70, zorder=3, label=f"{st}-like")
ax.axvspan(2.5, 5.5, color="#FDECE8", zorder=0)
ax.set_xticks(xi); ax.set_xticklabels(order, fontsize=11); ax.set_ylim(0, 0.8)
ax.set_ylabel("fraction (CIBERSORT)", fontsize=13)
ax.legend(frameon=False, fontsize=11, ncol=4, loc="upper center")
bare(ax)
save(fig, "chart_states.png")

json.dump(facts, open(OUT / "robustness.json", "w", encoding="utf-8"), indent=1, default=float)
print(json.dumps(facts, indent=1, default=float)[:1500])
