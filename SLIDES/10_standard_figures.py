# -*- coding: utf-8 -*-
"""The v4 talk's result figures in the forms this literature uses, drawn at the size they occupy on the slide.

  fig_gsea_ti.png          GSEA enrichment plot in the Broad layout (running enrichment score, hit barcode, ranked-list
                           colour bar, ranked-list metric), translation initiation. The running score is recomputed from
                           the pipeline's ranked list (classic weighted, p = 1) and checked against the pipeline's own
                           running-score column at every hit before drawing.
  fig_le_heatmap.png       row z-scored rlog expression of that set's leading-edge (core enrichment) genes, six tumours
                           grouped by arm, with an annotation row for the human share of each library
  fig_subtype_heatmap.png  GSVA scores of the ten Neftel and Garofano signatures per tumour (DESeq2 vst counts,
                           ANALYSIS/graft_relation), with the two-group p and a Benjamini-Hochberg q across the ten
  fig_graft_corr.png       each score against the human share of the library: ordinary least squares fit, 95 %
                           confidence band, Pearson r and p (two panels: translation initiation, astrocyte-like)
  fig_loo_heatmap.png      translation-initiation FDR q in all 41 leave-one-tumour-out runs (condition x permutation
                           seed), with the NES range of each condition
  fig_cibersort.png        CIBERSORT v1.04 fractions of the four Neftel states per library, stacked, with the fit r
  fig_pca.png              PCA of the six tumours as DESeq2 plotPCA draws it (rlog, 500 most variable genes)
  fig_sample_dist.png      sample-to-sample Euclidean distances on rlog, hierarchically clustered (DESeq2 vignette QC)
  fig_gsea_dotplot.png     the strongest gene sets in both directions: NES, FDR q, set size
  fig_depmap.png           DepMap 24Q4 Chronos gene effect across every screened line (violin) with U-251 MG marked,
                           for the genes of the ciclopirox axes (needs depmap/CRISPRGeneEffect.csv + Model.csv, figshare
                           27993248; skipped with a message when absent)
  standard.json            every number the slides quote from these figures

    python SLIDES/10_standard_figures.py
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from scipy import stats  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
from u251_paths import A, COUNTS, GSEA_LOO, RLOG, gsea_file  # noqa: E402

OUT = HERE / "figures_cns"
PRI, REC = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]
TUM = PRI + REC
TI = "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"
C_PRI, C_REC, INK = "#191D63", "#E08214", "#000000"          # one arm pair on every slide: primary navy, recurrent orange
plt.rcParams.update({"font.family": "Arial", "font.size": 16, "axes.edgecolor": INK, "axes.linewidth": 1.1,
                     "xtick.color": INK, "ytick.color": INK, "xtick.direction": "out", "ytick.direction": "out",
                     "xtick.major.width": 1.1, "ytick.major.width": 1.1, "xtick.major.size": 5, "ytick.major.size": 5,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
                     "legend.frameon": False, "figure.dpi": 220, "savefig.dpi": 220, "axes.unicode_minus": True})
MINUS = lambda s: s.replace("-", "−")  # noqa: E731
facts = {}


def save(fig, name):
    fig.savefig(OUT / name, facecolor="white", bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    print("  ", name)


def zrow(m):
    return (m - m.mean(axis=1, keepdims=True)) / m.std(axis=1, ddof=1, keepdims=True)


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


meta = pd.read_csv(A / "metadata_full.csv").set_index("sample")
graft = meta.loc[TUM, "graft_pct"].astype(float)
RDBU = plt.get_cmap("RdBu_r")

# ====================================================================== GSEA enrichment plot, Broad layout
ranked = pd.read_csv(gsea_file("ranked_gene_list_Recurrent_U2_versus_Primary_U2.tsv"), sep="\t")[["NAME", "SCORE"]]
rep = pd.concat([pd.read_csv(gsea_file(f"gsea_report_for_{g}_U2.tsv"), sep="\t") for g in ("Primary", "Recurrent")])
row = rep[rep.NAME == TI].iloc[0]
st = pd.read_csv(gsea_file(f"{TI}.tsv"), sep="\t").sort_values("RANK IN GENE LIST")
N = len(ranked); score = ranked.SCORE.to_numpy(float)
hits = st["RANK IN GENE LIST"].to_numpy(int)
assert all(ranked.NAME.iloc[h] == s for h, s in zip(hits, st.SYMBOL)), "set-table ranks are not positions in the ranked list"
isin = np.zeros(N, bool); isin[hits] = True
w = np.where(isin, np.abs(score), 0.0)
es = np.cumsum(w) / w.sum() - np.cumsum(~isin) / (N - isin.sum())
err = float(np.max(np.abs(es[hits] - st["RUNNING ES"].to_numpy(float))))
assert err < 1e-3, f"recomputed running score differs from the pipeline's by {err}"
lead = st[st["CORE ENRICHMENT"].astype(str).str.strip() == "Yes"].SYMBOL.astype(str).tolist()
RP = st.SYMBOL.astype(str).str.match(r"^(RPL|RPS|FAU$|UBA52$)")
assert RP.all(), "the KEGG MEDICUS translation-initiation set is no longer all ribosomal-protein genes"
rp_other = {}
for k in ("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "KEGG_RIBOSOME", "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
          "REACTOME_SELENOAMINO_ACID_METABOLISM", "REACTOME_CELLULAR_RESPONSE_TO_STARVATION"):
    t_ = pd.read_csv(gsea_file(f"{k}.tsv"), sep="\t").SYMBOL.astype(str)
    rp_other[k] = [int(t_.str.match(r"^(RPL|RPS|FAU$|UBA52$)").sum()), int(len(t_))]
loo = pd.read_csv(GSEA_LOO / "loo_sets.tsv", sep="\t")
ti_full = loo[(loo.set == TI) & (loo.condition == "full")]
facts["gsea"] = {"N_ranked": N, "n_set": int(len(st)), "n_leading_edge": len(lead), "NES": float(row.NES), "nom_p": float(row["NOM p-val"]),
                 "q_published": float(row["FDR q-val"]), "q_seed_min": float(ti_full.fdr_q.min()), "q_seed_max": float(ti_full.fdr_q.max()),
                 "es_min": float(es.min()), "first_hit_rank1": int(hits.min()) + 1, "es_recompute_max_err": err,
                 "zero_cross_rank1": int(np.argmax(score < 0)) + 1, "metric": "Diff_of_Classes",
                 "ti_all_ribosomal_protein": bool(RP.all()), "rp_in_other_leading_sets": rp_other}

fig = plt.figure(figsize=(5.9, 4.4))
gs = GridSpec(4, 1, height_ratios=[3.2, 0.55, 0.22, 1.7], hspace=0.0, figure=fig)
a0, a1, a2, a3 = (fig.add_subplot(gs[i]) for i in range(4))
x = np.arange(N)
a0.plot(x, es, color="#1B9E3E", lw=2.4)
a0.axhline(0, color=INK, lw=0.8)
a0.set_xlim(0, N - 1); a0.set_xticks([])
a0.set_ylabel("Enrichment score", fontsize=17)
a0.spines["bottom"].set_visible(False)
a0.set_ylim(min(-1.05, es.min() - 0.08), 0.12)
a0.set_title("Translation initiation: %d ribosomal-protein genes" % len(st), fontsize=16, loc="left", pad=6)
a0.text(0.03, 0.10, f"NES = {MINUS(f'{row.NES:.2f}')}\nNominal P < 0.002\nFDR q = {facts['gsea']['q_seed_min']:.3f}–{facts['gsea']['q_seed_max']:.2f} (5 seeds)",
        transform=a0.transAxes, fontsize=16, va="bottom", ha="left", linespacing=1.3)
a1.vlines(hits, 0, 1, color=INK, lw=0.8)
a1.set_xlim(0, N - 1); a1.set_ylim(0, 1); a1.axis("off")
a2.imshow(np.linspace(1, -1, 512)[None, :], aspect="auto", cmap=RDBU, extent=(0, N - 1, 0, 1))
a2.axis("off")
a3.fill_between(x, score, 0, color="#8C8C8C", lw=0)
a3.axhline(0, color=INK, lw=0.8)
zc = facts["gsea"]["zero_cross_rank1"]
a3.axvline(zc, color=INK, lw=0.8, ls=(0, (3, 2)))
a3.set_xlim(0, N - 1)
lim = np.percentile(np.abs(score), 99.5)
a3.set_ylim(-lim, lim)
a3.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: MINUS(f"{v / 1000:.0f}") + ("k" if v else "")))
a3.set_ylabel("Metric\n(clipped)", fontsize=14)
a3.set_xlabel("Rank in ordered gene list", fontsize=17)
a3.set_xticks([0, 5000, 10000, 15000, N - 1]); a3.set_xticklabels(["1", "5k", "10k", "15k", f"{N / 1000:.1f}k"], fontsize=15)
a3.tick_params(axis="y", labelsize=14)
a3.text(0.01, 0.95, "Higher in recurrent", transform=a3.transAxes, fontsize=15, color=C_REC, va="top")
a3.text(0.99, 0.05, "Higher in primary", transform=a3.transAxes, fontsize=15, color="#2166AC", va="bottom", ha="right")
a0.tick_params(axis="y", labelsize=15)
save(fig, "fig_gsea_ti.png")
facts["gsea"]["metric_clip_abs"] = float(lim)

# ====================================================================== leading-edge heatmap
rl = pd.read_csv(RLOG, sep="\t").set_index("gene_id")[TUM]
cnt = pd.read_csv(COUNTS, sep="\t", usecols=["gene_id", "gene_name"]).set_index("gene_id")
rl["symbol"] = cnt.reindex(rl.index).gene_name.values
rl = rl.dropna(subset=["symbol"])
rl["mean"] = rl[TUM].mean(axis=1)
rl = rl.sort_values("mean", ascending=False).drop_duplicates("symbol").set_index("symbol")
le = [g for g in lead if g in rl.index]
Z = zrow(rl.loc[le, TUM].to_numpy(float))
order = leaves_list(linkage(Z, "average", metric="correlation"))
Z = Z[order]
top_col = int(np.argmax(np.bincount(np.argmax(rl.loc[le, TUM].to_numpy(float), axis=1), minlength=6)))
facts["leading_edge"] = {"n_genes": len(le), "n_missing_in_rlog": len(lead) - len(le),
                         "highest_tumour": {s: int(n) for s, n in zip(TUM, np.bincount(np.argmax(rl.loc[le, TUM].to_numpy(float), axis=1), minlength=6))},
                         "mean_z": {s: float(v) for s, v in zip(TUM, Z.mean(axis=0))}}
fig = plt.figure(figsize=(4.7, 4.5))
gs = GridSpec(3, 2, height_ratios=[0.6, 0.6, 8], width_ratios=[1, 0.06], hspace=0.08, wspace=0.08, figure=fig)
ag, ah, am, cb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])
ag.imshow([[0, 0, 0, 1, 1, 1]], aspect="auto", cmap=LinearSegmentedColormap.from_list("g", [C_PRI, C_REC]), vmin=0, vmax=1)
ag.set_xticks([]); ag.set_yticks([0]); ag.set_yticklabels(["Arm"], fontsize=15)
gv = graft.to_numpy()[None, :]
ah.imshow(gv, aspect="auto", cmap="Greys", vmin=0, vmax=80)
for i, v in enumerate(gv[0]):
    ah.text(i, 0, f"{v:.0f}", ha="center", va="center", fontsize=14, color="white" if v > 45 else INK)
ah.set_xticks([]); ah.set_yticks([0]); ah.set_yticklabels(["Human %"], fontsize=15)
im = am.imshow(Z, aspect="auto", cmap=RDBU, norm=TwoSlopeNorm(0, -2, 2), interpolation="nearest")
am.set_yticks([]); am.set_xticks(range(6)); am.set_xticklabels(TUM, rotation=45, ha="right", fontsize=15)
am.set_ylabel(f"{len(le)} leading-edge genes\n(ribosomal proteins)", fontsize=16)
for a in (ag, ah, am):
    for s in a.spines.values():
        s.set_visible(False)
    a.tick_params(length=0)
am.axvline(2.5, color="white", lw=3); ag.axvline(2.5, color="white", lw=3); ah.axvline(2.5, color="white", lw=3)
ag.text(1, 0, "Primary", ha="center", va="center", color="white", fontsize=15, fontweight="bold")
ag.text(4, 0, "Recurrent", ha="center", va="center", color="white", fontsize=15, fontweight="bold")
c = fig.colorbar(im, cax=cb); c.set_label("Row z-score", fontsize=16); c.ax.tick_params(labelsize=14)
c.set_ticks([-2, -1, 0, 1, 2]); c.set_ticklabels([MINUS("-2"), MINUS("-1"), "0", "1", "2"])
save(fig, "fig_le_heatmap.png")

# the same set on a log-scale ranking (difference of rlog means), since Diff_of_Classes on counts weights abundant genes
dl = (rl[REC].mean(axis=1) - rl[PRI].mean(axis=1)).sort_values(ascending=False)
pos_ = {g_: i for i, g_ in enumerate(dl.index)}
rk = np.array([pos_[g_] for g_ in st.SYMBOL.astype(str) if g_ in pos_])
facts["gsea"]["log_scale"] = {"n_genes_ranked": int(len(dl)), "n_set_found": int(len(rk)),
                              "share_in_last_6_6pct": float((rk >= (1 - 0.066) * len(dl)).mean()), "share_in_lower_half": float((rk >= 0.5 * len(dl)).mean()),
                              "median_percentile_from_top": float(np.median(rk) / len(dl) * 100)}
facts["gsea"]["share_in_last_6_6pct_counts_metric"] = float((hits >= (1 - 0.066) * N).mean())

# ====================================================================== PCA, as DESeq2 plotPCA (rlog, top 500 by variance, centred)
r0 = pd.read_csv(RLOG, sep="\t").set_index("gene_id")[TUM]
top = r0.loc[r0.var(axis=1).sort_values(ascending=False).index[:500]]
Xc = top.T.to_numpy(float); Xc = Xc - Xc.mean(axis=0)
U_, S_, Vt_ = np.linalg.svd(Xc, full_matrices=False)
pcs = U_ * S_; ve = S_ ** 2 / np.sum(S_ ** 2)
if pcs[[TUM.index(s_) for s_ in REC], 0].mean() < pcs[[TUM.index(s_) for s_ in PRI], 0].mean():
    pcs[:, 0] *= -1                                        # recurrent to the right, as on every earlier version
pc = pd.DataFrame(pcs[:, :2], index=TUM, columns=["PC1", "PC2"])
prefix = np.array([s_[:2] == "NL" for s_ in TUM], float); arm = np.array([s_ in REC for s_ in TUM], float)
r2 = lambda a_, b_: float(np.corrcoef(a_, b_)[0, 1] ** 2)  # noqa: E731
facts["pca"] = {"var_pc1": float(ve[0]), "var_pc2": float(ve[1]), "coords": {s_: [float(a_), float(b_)] for s_, (a_, b_) in pc.iterrows()},
                "margin_pc1": float(pc.loc[REC, "PC1"].min() - pc.loc[PRI, "PC1"].max()), "pc1_range": float(pc.PC1.max() - pc.PC1.min()),
                "r2_pc1_arm": r2(pc.PC1, arm), "r2_pc1_prefix": r2(pc.PC1, prefix), "input": "rlog (all.rlog.tsv), 500 most variable genes, centred"}
fig, ax = plt.subplots(figsize=(5.2, 4.4))
for grp, col, mk, lab in ((PRI, C_PRI, "o", "Primary"), (REC, C_REC, "s", "Recurrent")):
    ax.scatter(pc.loc[grp, "PC1"], pc.loc[grp, "PC2"], s=150, color=col, marker=mk, edgecolor="white", lw=1.2, zorder=3, label=lab)
for s_, (a_, b_) in pc.iterrows():
    ax.annotate(s_, (a_, b_), xytext=(0, 11), textcoords="offset points", ha="center", fontsize=14, color="#444444")
ax.axhline(0, color="#BDBDBD", lw=0.8, zorder=1); ax.axvline(0, color="#BDBDBD", lw=0.8, zorder=1)
ax.set_xlabel(f"PC1: {100 * ve[0]:.0f} % variance", fontsize=17); ax.set_ylabel(f"PC2: {100 * ve[1]:.0f} % variance", fontsize=17)
lim_ = float(np.abs(pcs[:, :2]).max()) * 1.25
ax.set_xlim(-lim_, lim_); ax.set_ylim(-lim_, lim_); ax.set_aspect("equal")
ax.tick_params(labelsize=15)
ax.legend(fontsize=14, loc="upper left", handletextpad=0.2, borderaxespad=0.2)
save(fig, "fig_pca.png")

# ====================================================================== sample-to-sample distances (DESeq2 vignette), clustered
from scipy.cluster.hierarchy import dendrogram, fcluster  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402
Xall = r0.T.to_numpy(float)
dv = pdist(Xall, "euclidean"); Dm = squareform(dv)
Lk = linkage(dv, "complete")
lv = leaves_list(Lk); names_ = [TUM[i] for i in lv]
cl2 = fcluster(Lk, 2, "maxclust")
facts["sample_dist"] = {"order": names_, "k2_clusters": {s_: int(c_) for s_, c_ in zip(TUM, cl2)},
                        "k2_recovers_arms": bool(len({c_ for s_, c_ in zip(TUM, cl2) if s_ in PRI}) == 1 and len({c_ for s_, c_ in zip(TUM, cl2) if s_ in REC}) == 1
                                                 and {c_ for s_, c_ in zip(TUM, cl2) if s_ in PRI} != {c_ for s_, c_ in zip(TUM, cl2) if s_ in REC}),
                        "method": "Euclidean distance on rlog, all genes; complete linkage"}
fig = plt.figure(figsize=(6.4, 5.4))
gs = GridSpec(4, 2, height_ratios=[1.3, 0.35, 0.35, 6], width_ratios=[1, 0.05], hspace=0.05, wspace=0.06, figure=fig)
adn, ag, ah, am, cb = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[3, 0]), fig.add_subplot(gs[3, 1]))
dendrogram(Lk, ax=adn, color_threshold=0, above_threshold_color="#444444", no_labels=True)
adn.axis("off")
armv = [[1.0 if s_ in REC else 0.0 for s_ in names_]]
ag.imshow(armv, aspect="auto", cmap=LinearSegmentedColormap.from_list("g", [C_PRI, C_REC]), vmin=0, vmax=1)
ag.set_xticks([]); ag.set_yticks([0]); ag.set_yticklabels(["Arm"], fontsize=14)
gvo = np.array([[graft[s_] for s_ in names_]])
ah.imshow(gvo, aspect="auto", cmap="Greys", vmin=0, vmax=80)
for i, v_ in enumerate(gvo[0]):
    ah.text(i, 0, f"{v_:.0f}", ha="center", va="center", fontsize=13, color="white" if v_ > 45 else INK)
ah.set_xticks([]); ah.set_yticks([0]); ah.set_yticklabels(["Human %"], fontsize=14)
im = am.imshow(Dm[np.ix_(lv, lv)], cmap="Blues_r", interpolation="nearest", aspect="auto")
am.set_xticks(range(6)); am.set_xticklabels(names_, rotation=45, ha="right", fontsize=15)
am.set_yticks(range(6)); am.set_yticklabels(names_, fontsize=15)
for a in (ag, ah, am):
    for s in a.spines.values():
        s.set_visible(False)
    a.tick_params(length=0)
c = fig.colorbar(im, cax=cb); c.set_label("Euclidean distance (rlog)", fontsize=15); c.ax.tick_params(labelsize=13)
from matplotlib.patches import Patch  # noqa: E402
am.legend([Patch(color=C_PRI), Patch(color=C_REC)], ["Primary", "Recurrent"], loc="upper left", bbox_to_anchor=(1.18, 1.32), fontsize=13, handlelength=1)
save(fig, "fig_sample_dist.png")

# ====================================================================== subtype GSVA heatmap
ps = pd.read_csv(A / "graft_relation/results/per_sample.tsv", sep="\t").set_index("sample").loc[TUM]
cr = pd.read_csv(A / "graft_relation/results/score_correlations.tsv", sep="\t").set_index("score")
SIG = [("Neftel_MES1", "Neftel MES1-like"), ("Neftel_MES2", "Neftel MES2-like"), ("Neftel_AC", "Neftel AC-like"),
       ("Neftel_OPC", "Neftel OPC-like"), ("Neftel_NPC1", "Neftel NPC1-like"), ("Neftel_NPC2", "Neftel NPC2-like"),
       ("Garofano_GPM", "Garofano GPM"), ("Garofano_MTC", "Garofano MTC"), ("Garofano_NEU", "Garofano NEU"),
       ("Garofano_PPR", "Garofano PPR")]
keys = [k for k, _ in SIG]
p = np.array([cr.loc[k, "group_p"] for k in keys]); q = bh(p)
padj = np.array([cr.loc[k, "group_p_adj"] for k in keys])
diff = np.array([ps.loc[REC, k].mean() - ps.loc[PRI, k].mean() for k in keys])
facts["subtypes"] = {k: {"primary": float(ps.loc[PRI, k].mean()), "recurrent": float(ps.loc[REC, k].mean()), "diff": float(d),
                         "p": float(pp), "q": float(qq), "p_adj_graft": float(pa), "r_graft": float(cr.loc[k, "r_graft"])}
                     for k, d, pp, qq, pa in zip(keys, diff, p, q, padj)}
facts["subtypes_input"] = "GSVA (Gaussian kernel) on DESeq2 vst counts of the six tumours; p from score ~ arm (two-sample t), q BH over the ten"
M = ps[keys].T.to_numpy(float)
vmax = float(np.ceil(np.abs(M).max() * 10) / 10)
fig = plt.figure(figsize=(7.4, 5.2))
gs = GridSpec(3, 3, height_ratios=[0.7, 0.7, 10], width_ratios=[6, 2.1, 0.18], hspace=0.06, wspace=0.05, figure=fig)
ag, ah, am, at, cb = (fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1]), fig.add_subplot(gs[2, 2]))
ah.imshow(graft.to_numpy()[None, :], aspect="auto", cmap="Greys", vmin=0, vmax=80)
for i, v_ in enumerate(graft.to_numpy()):
    ah.text(i, 0, f"{v_:.0f}", ha="center", va="center", fontsize=14, color="white" if v_ > 45 else INK)
ah.set_xticks([]); ah.set_yticks([0]); ah.set_yticklabels(["Human %"], fontsize=15)
for s in ah.spines.values():
    s.set_visible(False)
ah.tick_params(length=0); ah.axvline(2.5, color="white", lw=3)
ag.imshow([[0, 0, 0, 1, 1, 1]], aspect="auto", cmap=LinearSegmentedColormap.from_list("g", [C_PRI, C_REC]), vmin=0, vmax=1)
ag.text(1, 0, "Primary", ha="center", va="center", color="white", fontsize=15, fontweight="bold")
ag.text(4, 0, "Recurrent", ha="center", va="center", color="white", fontsize=15, fontweight="bold")
ag.set_xticks([]); ag.set_yticks([])
im = am.imshow(M, aspect="auto", cmap=RDBU, norm=TwoSlopeNorm(0, -vmax, vmax), interpolation="nearest")
am.set_xticks(range(6)); am.set_xticklabels(TUM, rotation=45, ha="right", fontsize=15)
am.set_yticks(range(len(SIG))); am.set_yticklabels([n for _, n in SIG], fontsize=16)
for t_, k in zip(am.get_yticklabels(), keys):
    if k == "Neftel_AC":
        t_.set_fontweight("bold")
am.axhline(5.5, color="white", lw=3)
for a in (ag, am):
    for s in a.spines.values():
        s.set_visible(False)
    a.tick_params(length=0)
am.axvline(2.5, color="white", lw=3); ag.axvline(2.5, color="white", lw=3)
at.set_xlim(0, 2); at.set_ylim(len(SIG) - 0.5, -0.5); at.axis("off")
at.text(0.55, -0.62, "P", ha="center", va="bottom", fontsize=16, style="italic"); at.text(1.5, -0.62, "q", ha="center", va="bottom", fontsize=16, style="italic")
at.text(1.0, -1.35, "t-test, 3 v 3", ha="center", va="bottom", fontsize=12, color="#555555")
for i, (pp, qq) in enumerate(zip(p, q)):
    b = "bold" if pp < 0.05 else "normal"
    at.text(0.55, i, f"{pp:.3f}" if pp < 0.01 else f"{pp:.2f}", ha="center", va="center", fontsize=16, fontweight=b)
    at.text(1.5, i, f"{qq:.3f}" if qq < 0.1 else f"{qq:.2f}", ha="center", va="center", fontsize=16, fontweight=b)
c = fig.colorbar(im, cax=cb); c.set_label("GSVA score", fontsize=16); c.ax.tick_params(labelsize=14)
save(fig, "fig_subtype_heatmap.png")

# ====================================================================== score against the human share: OLS + 95 % band
fig, axs = plt.subplots(1, 2, figsize=(9.8, 4.2))
facts["graft_corr"] = {}
OFF = {TI: {"IL66B": (0, -22, "center"), "NL71B": (0, 11, "center"), "IL67B": (-9, 9, "right"), "NL70B": (9, -18, "left"),
            "IL68B": (9, 5, "left"), "IL69B": (0, -22, "center")},
       "Neftel_AC": {"IL66B": (0, -22, "center"), "NL71B": (0, 11, "center"), "IL67B": (-8, 11, "right"), "IL68B": (8, 9, "left"),
                     "NL70B": (10, -5, "left"), "IL69B": (-9, 10, "right")}}
for ax, (k, lab) in zip(axs, ((TI, "Ribosomal-protein set\n(GSVA score)"), ("Neftel_AC", "Neftel AC-like\n(GSVA score)"))):
    xv, yv = graft.to_numpy(), ps[k].to_numpy(float)
    res = stats.linregress(xv, yv)
    xx = np.linspace(xv.min() - 3, xv.max() + 3, 100)
    n = len(xv); tq = stats.t.ppf(0.975, n - 2)
    s_err = np.sqrt(np.sum((yv - (res.intercept + res.slope * xv)) ** 2) / (n - 2))
    band = tq * s_err * np.sqrt(1 / n + (xx - xv.mean()) ** 2 / np.sum((xv - xv.mean()) ** 2))
    yy = res.intercept + res.slope * xx
    ax.fill_between(xx, yy - band, yy + band, color="#BDBDBD", alpha=0.45, lw=0)
    ax.plot(xx, yy, color="#555555", lw=1.6)
    for grp, col, mk in ((PRI, C_PRI, "o"), (REC, C_REC, "s")):
        ax.scatter(graft[grp], ps.loc[grp, k], s=110, color=col, marker=mk, edgecolor="white", lw=1.2, zorder=3)
        for s_ in grp:
            dx, dy, ha = OFF[k].get(s_, (8, 6, "left"))
            ax.annotate(s_, (graft[s_], ps.loc[s_, k]), xytext=(dx, dy), textcoords="offset points", fontsize=14, color="#444444", ha=ha)
    ax.set_xlabel("Human reads in library (%)", fontsize=17)
    ax.set_ylabel(lab, fontsize=17)
    ax.text(0.03, 0.97, f"r = {res.rvalue:.2f}, P = {res.pvalue:.2f}", transform=ax.transAxes, va="top", fontsize=16)
    ax.tick_params(labelsize=15)
    facts["graft_corr"][k] = {"r": float(res.rvalue), "p": float(res.pvalue), "slope": float(res.slope)}
axs[1].scatter([], [], s=90, color=C_PRI, marker="o", label="Primary")
axs[1].scatter([], [], s=90, color=C_REC, marker="s", label="Recurrent")
axs[1].legend(fontsize=15, loc="lower right", handletextpad=0.3)
fig.tight_layout(w_pad=2.5)
save(fig, "fig_graft_corr.png")

# ====================================================================== leave-one-tumour-out FDR q
ti = loo[loo.set == TI]
conds = ["full"] + [f"drop_{s}" for s in TUM]
names = ["All six"] + [f"Without {s}" for s in TUM]
seeds = [1234, 1, 2, 3, 4]
Q = np.array([[float(ti[(ti.condition == c) & (ti.seed == sd)].fdr_q.iloc[0]) for sd in seeds] for c in conds])
nes_rng = [(ti[ti.condition == c].NES.min(), ti[ti.condition == c].NES.max()) for c in conds]
facts["loo"] = {"runs": int(len(ti)), "nes_min": float(ti.NES.min()), "nes_max": float(ti.NES.max()), "nom_p_max": float(ti.nom_p.max()),
                "n_q05": int((ti.fdr_q < 0.05).sum()), "q_min": float(ti.fdr_q.min()), "q_max": float(ti.fdr_q.max()),
                "grid": {n: [float(v) for v in r] for n, r in zip(names, Q)},
                "all_seeds_q05": [n for n, r in zip(names, Q) if (r < 0.05).all()], "no_seed_q05": [n for n, r in zip(names, Q) if (r >= 0.05).all()],
                "nofilter_runs": int(ti.condition.str.endswith("_nofilter").sum())}
fig = plt.figure(figsize=(7.2, 4.9))
gs = GridSpec(1, 2, width_ratios=[5, 1.75], wspace=0.04, figure=fig)
ax, at = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
SIGC = "#2166AC"
ax.imshow((Q < 0.05).astype(float), aspect="auto", cmap=LinearSegmentedColormap.from_list("q", ["#F0F0F0", SIGC]), vmin=0, vmax=1)
for i in range(Q.shape[0]):
    for j in range(Q.shape[1]):
        v = Q[i, j]
        ax.text(j, i, f"{v:.3f}" if v < 0.1 else f"{v:.2f}", ha="center", va="center", fontsize=15,
                color="white" if v < 0.05 else INK, fontweight="bold" if v < 0.05 else "normal")
ax.set_xticks(np.arange(-0.5, len(seeds)), minor=True); ax.set_yticks(np.arange(-0.5, len(names)), minor=True)
ax.grid(which="minor", color="white", lw=1.5); ax.tick_params(which="minor", length=0)
ax.set_xticks(range(len(seeds))); ax.set_xticklabels([f"{s}" for s in seeds], fontsize=15)
ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top")
ax.set_xlabel("Permutation seed (1234 = published)", fontsize=16, labelpad=8)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=16)
ax.axhline(0.5, color="white", lw=3); ax.axhline(3.5, color="white", lw=3)
for s in ax.spines.values():
    s.set_visible(False)
ax.tick_params(length=0)
at.set_ylim(len(names) - 0.5, -0.5); at.set_xlim(0, 1); at.axis("off")
at.text(0.5, -0.62, "NES", ha="center", va="bottom", fontsize=16)
for i, (lo_, hi_) in enumerate(nes_rng):
    at.text(0.5, i, MINUS(f"{lo_:.2f} to {hi_:.2f}"), ha="center", va="center", fontsize=14)
save(fig, "fig_loo_heatmap.png")

# ====================================================================== CIBERSORT v1.04, stacked
v = pd.read_csv(A / "cibersort/results/v104/fractions_v104_neftel4_confident.tsv", sep="\t").set_index("sample")
CTL = ["IL64B", "N168B", "N269B"]
order = TUM + ["C2B"] + CTL
STATES = [("MES", "MES-like", "#C0392B"), ("AC", "AC-like", "#E6A817"), ("NPC", "NPC-like", "#2E86AB"), ("OPC", "OPC-like", "#6BAA75")]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
xs = np.array([0, 1, 2, 3.5, 4.5, 5.5, 7.0, 8.5, 9.5, 10.5])
bottom = np.zeros(len(order))
for k, lab, col in STATES:
    vals = v.loc[order, k].to_numpy(float)
    ax.bar(xs, vals, bottom=bottom, width=0.78, color=col, edgecolor="white", lw=0.8, label=lab)
    bottom += vals
ax.set_xticks(xs); ax.set_xticklabels(order, fontsize=14, rotation=45, ha="right")
ax.set_ylim(0, 1.0); ax.set_ylabel("Estimated fraction", fontsize=16)
for x_, s_ in zip(xs, order):
    ax.text(x_, -0.25, f"{v.loc[s_, 'Correlation']:.2f}", ha="center", va="top", fontsize=14, color="#444444", transform=ax.get_xaxis_transform())
ax.text(-0.9, -0.25, "fit r", ha="right", va="top", fontsize=14, color="#444444", style="italic", transform=ax.get_xaxis_transform())
for x_, s_ in zip(xs, order):
    ax.text(x_, -0.34, f"{v.loc[s_, 'P-value']:.3f}"[1:], ha="center", va="top", fontsize=12, color="#444444", transform=ax.get_xaxis_transform())
ax.text(-0.9, -0.34, "P", ha="right", va="top", fontsize=14, color="#444444", style="italic", transform=ax.get_xaxis_transform())
for xc, t_ in ((1, "Primary"), (4.5, "Recurrent"), (7.0, "Culture"), (9.5, "Rat brain")):
    ax.text(xc, 1.02, t_, ha="center", va="bottom", fontsize=15, transform=ax.get_xaxis_transform(), fontweight="bold",
            color=C_PRI if t_ == "Primary" else C_REC if t_ == "Recurrent" else "#666666")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=15, handlelength=1.2, reverse=True)
ax.tick_params(labelsize=14)
save(fig, "fig_cibersort.png")
facts["cibersort"] = {s: {k: float(v.loc[s, k]) for k in ("MES", "AC", "NPC", "OPC", "Correlation", "P-value") if k in v.columns} for s in order}

# ====================================================================== DepMap, every line with U-251 MG marked
GENES = [("RRM1", "Ribonucleotide\nreductase"), ("RRM2", "Ribonucleotide\nreductase"), ("DOHH", "eIF5A\nhypusination"),
         ("DHPS", "eIF5A\nhypusination"), ("EIF5A", "eIF5A\nhypusination"), ("FTH1", "Iron\nhandling"), ("GPX4", "Iron\nhandling"), ("TFRC", "Iron\nhandling"),
         ("HIF1A", "PHD–HIF"), ("VHL", "PHD–HIF")]
dm_path = ROOT / "depmap" / "CRISPRGeneEffect.csv"
if dm_path.exists():
    hdr = pd.read_csv(dm_path, nrows=0).columns
    want = {g: next(c for c in hdr if c.split(" (")[0] == g) for g, _ in GENES}
    eff = pd.read_csv(dm_path, usecols=[hdr[0]] + list(want.values()), index_col=0)
    eff.columns = [c.split(" (")[0] for c in eff.columns]
    U = "ACH-000232"
    assert U in eff.index
    fig, ax = plt.subplots(figsize=(10.0, 4.3))
    pos, x0, last, heads = [], 0.0, None, []
    for g, grp_ in GENES:
        if last is not None and grp_ != last:
            x0 += 0.7
        pos.append(x0); last = grp_; x0 += 1.45
    data = [eff[g].dropna().to_numpy() for g, _ in GENES]
    vp = ax.violinplot(data, positions=pos, widths=0.8, showextrema=False)
    for b in vp["bodies"]:
        b.set_facecolor("#C9C9C9"); b.set_edgecolor("#7F7F7F"); b.set_alpha(1); b.set_lw(0.8)
    for x_, d in zip(pos, data):
        q1, med, q3 = np.percentile(d, [25, 50, 75])
        ax.plot([x_, x_], [q1, q3], color="#333333", lw=3.5, solid_capstyle="butt")
        ax.scatter([x_], [med], s=22, color="white", zorder=3)
    uv = [float(eff.loc[U, g]) for g, _ in GENES]
    mdl = pd.read_csv(ROOT / "depmap" / "Model.csv", usecols=["ModelID", "OncotreeLineage"]).set_index("ModelID")
    cns = [m_ for m_ in eff.index if m_ in mdl.index and mdl.loc[m_, "OncotreeLineage"] == "CNS/Brain"]
    rng_ = np.random.default_rng(7)
    for j_, (x_, (g, _)) in enumerate(zip(pos, GENES)):
        yv_ = eff.loc[cns, g].dropna().to_numpy()
        ax.scatter(x_ + rng_.uniform(-0.16, 0.16, len(yv_)), yv_, s=7, color="#3A3A3A", alpha=0.55, lw=0, zorder=2,
                   label="CNS/brain lines" if j_ == 0 else None)
    ax.scatter(pos, uv, s=130, marker="D", color="#D62728", edgecolor="white", lw=1.0, zorder=4, label="U-251 MG")
    ax.axhline(0, color="#7F7F7F", lw=0.9)
    ax.axhline(-0.5, color="#7F7F7F", lw=1.1, ls=(0, (4, 3)))
    ax.axhline(-1, color="#7F7F7F", lw=0.9, ls=(0, (1, 2)))
    ax.set_xticks(pos); ax.set_xticklabels([g for g, _ in GENES], fontsize=15, style="italic")
    ax.set_ylabel("Gene effect (Chronos)", fontsize=17)
    ax.tick_params(labelsize=15)
    lo_ = min(float(np.concatenate(data).min()), min(uv)) - 0.1
    ax.set_ylim(lo_, 1.05)
    seen = {}
    for x_, (g, grp_) in zip(pos, GENES):
        seen.setdefault(grp_, []).append(x_)
    for grp_, xl in seen.items():
        ax.text(np.mean(xl), 1.02, grp_, ha="center", va="bottom", fontsize=15, transform=ax.get_xaxis_transform())
    ax.legend(loc="lower right", fontsize=14, handletextpad=0.2, markerscale=1.0)
    save(fig, "fig_depmap.png")
    facts["depmap"] = {"n_lines": int(eff.shape[0]), "release": "DepMap 24Q4 (figshare 27993248)",
                       "n_cns": len(cns),
                       "genes": {g: {"u251": float(eff.loc[U, g]), "median": float(eff[g].median()), "n": int(eff[g].notna().sum()),
                                     "percentile_u251": float((eff[g].dropna() < eff.loc[U, g]).mean() * 100),
                                     "share_below_minus05": float((eff[g].dropna() < -0.5).mean()),
                                     "cns_median": float(eff.loc[cns, g].median()), "cns_n": int(eff.loc[cns, g].notna().sum()),
                                     "cns_more_dependent_than_pct": float((eff.loc[cns, g].dropna() > eff.loc[U, g]).mean() * 100)} for g, _ in GENES}}
else:
    print("  fig_depmap.png skipped: depmap/CRISPRGeneEffect.csv absent")

# ====================================================================== GSEA summary dot plot, both directions
def short(n_):
    col_, rest = n_.split("_", 1)
    t_ = rest.replace("_", " ").lower()
    t_ = t_[0].upper() + t_[1:]
    return (t_ if len(t_) <= 40 else t_[:39] + "…") + f" ({col_.split('_')[0].title() if col_ != 'KEGG' else 'KEGG'})"
repn = rep.copy(); repn["NES"] = repn.NES.astype(float); repn["q"] = repn["FDR q-val"].astype(float)
dn = repn.sort_values("NES").head(8); up = repn.sort_values("NES", ascending=False).head(8)
dd = pd.concat([up.sort_values("NES"), dn.sort_values("NES", ascending=False)]).reset_index(drop=True)
dd = pd.concat([dn.sort_values("NES", ascending=False), up.sort_values("NES")]).reset_index(drop=True)
fig, ax = plt.subplots(figsize=(9.6, 5.4))
from matplotlib.colors import LogNorm  # noqa: E402
sc = ax.scatter(dd.NES, range(len(dd)), s=dd.SIZE.astype(float) * 1.6, c=dd.q.clip(lower=1e-3), cmap="viridis", norm=LogNorm(1e-3, 1), edgecolor="#333333", lw=0.6, zorder=3)
ax.set_yticks(range(len(dd))); ax.set_yticklabels([short(n_) for n_ in dd.NAME], fontsize=12.5)
ax.axvline(0, color="#7F7F7F", lw=0.9); ax.axhline(7.5, color="#BDBDBD", lw=0.9, ls=(0, (3, 3)))
ax.set_xlabel("Normalized enrichment score (recurrent vs primary)", fontsize=15); ax.tick_params(axis="x", labelsize=14)
ax.set_xlim(-2.4, 2.4)
cbx = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.04); cbx.set_label("FDR q", fontsize=14); cbx.ax.tick_params(labelsize=12)
for sz in (50, 100, 200):
    ax.scatter([], [], s=sz * 1.6, color="white", edgecolor="#333333", lw=0.6, label=f"{sz}")
ax.legend(title="Set size", fontsize=12, title_fontsize=12, loc="lower right", labelspacing=1.1, borderpad=0.8)
save(fig, "fig_gsea_dotplot.png")
les = {}
for k in [TI] + list(rp_other):
    t_ = pd.read_csv(gsea_file(f"{k}.tsv"), sep="\t")
    les[k] = set(t_[t_["CORE ENRICHMENT"].astype(str).str.strip() == "Yes"].SYMBOL.astype(str))
facts["gsea"]["leading_edges"] = {"shared_by_all_six": len(set.intersection(*les.values())), "union": len(set.union(*les.values())),
                                  "union_rp": int(sum(1 for g_ in set.union(*les.values()) if g_.startswith(("RPL", "RPS")) or g_ in ("FAU", "UBA52")))}
down6 = repn[repn.NAME.isin([TI] + list(rp_other))]
facts["gsea"]["six_q"] = {n_: float(q_) for n_, q_ in zip(down6.NAME, down6.q)}
facts["gsea"]["best_up_q"] = float(repn[repn.NES > 0].q.min())
facts["gsea"]["n_down_q25"] = int(((repn.NES < 0) & (repn.q < 0.25)).sum()); facts["gsea"]["n_up_q25"] = int(((repn.NES > 0) & (repn.q < 0.25)).sum())

json.dump(facts, open(OUT / "standard.json", "w", encoding="utf-8"), indent=1, default=float)
print(json.dumps({k: facts[k] for k in facts if k in ("gsea", "leading_edge", "loo")}, indent=1, default=float)[:3000])
