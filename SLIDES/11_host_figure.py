# -*- coding: utf-8 -*-
"""Backup B7 (preliminary): the rat host tissue, recurrent vs primary. Volcano in the deck's style plus host.json.

Inputs (ANALYSIS/host/results_de/, copied from /rs/rs_grp_oschome/go2432/u251_host/results_de on 2026-09-26):
  host_therapy.deseq2.results.tsv       nf-core/differentialabundance DESeq2, Classification Recurrent_U2 vs Primary_U2
  host_tumour_vs_control.deseq2.results.tsv   Group Tumour vs Control (rat-brain controls)
  host_therapy.deseq2.sizefactors.tsv   DESeq2 size factors (for the per-sample consistency check)
  host_genes.tsv                        gene names and raw salmon counts (xengsort host bin, nf-core/rnaseq 3.22.2, mRatBN7.2)
  human_side_check.tsv                  the same genes in the HUMAN salmon counts (ANALYSIS/results_human_final)
The run crashed after DESeq2, in the exploratory plots (subset-to-contrast with two contrast variables), so there is no
gene-set test and no plot from the pipeline; the DESeq2 tables themselves were written before the crash.

    python SLIDES/11_host_figure.py
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
spec = importlib.util.spec_from_file_location("cns04", HERE / "04_make_cns_figures.py")
f4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f4)
plt, bare, save, place_labels = f4.plt, f4.bare, f4.save, f4.place_labels
NAVY, GREY, INK, FAINT, MINUS = f4.NAVY, f4.GREY, f4.INK, f4.FAINT, f4.MINUS
ORANGE = "#E08214"
HD = ROOT / "ANALYSIS" / "host" / "results_de"
PRI, REC = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]


def results(name):
    out = []
    for r in csv.DictReader(open(HD / f"{name}.deseq2.results.tsv"), delimiter="\t"):
        if r["padj"] not in ("NA", ""):
            out.append(dict(gene=r["gene_id"], lfc=float(r["log2FoldChange"]), padj=float(r["padj"]), base=float(r["baseMean"])))
    return out


genes = {r["gene_id"]: r for r in csv.DictReader(open(HD / "host_genes.tsv"), delimiter="\t")}
sf = {r["sample"]: float(r["sizeFactor"]) for r in csv.DictReader(open(HD / "host_therapy.deseq2.sizefactors.tsv"), delimiter="\t")}
name = lambda g: genes[g]["gene_name"] or g  # noqa: E731
norm = lambda g, s: float(genes[g][f"count_{s}"]) / sf[s]  # noqa: E731

th, tc = results("host_therapy"), results("host_tumour_vs_control")
sig = [r for r in th if r["padj"] < 0.05]
up, down = [r for r in sig if r["lfc"] > 0], [r for r in sig if r["lfc"] < 0]
sig2 = [r for r in sig if abs(r["lfc"]) >= 1]                     # the tumor volcano's rule: FDR < 0.05 and two-fold
up2, down2 = [r for r in sig2 if r["lfc"] > 0], [r for r in sig2 if r["lfc"] < 0]


def separates(r):
    p, q = [norm(r["gene"], s) for s in PRI], [norm(r["gene"], s) for s in REC]
    return min(q) > max(p) if r["lfc"] > 0 else max(q) < min(p)


reads = {s: sum(float(g[f"count_{s}"]) for g in genes.values() if g[f"count_{s}"]) / 1e6 for s in PRI + REC}
# salmon counts are read pairs assigned to genes, 30-45 % of them rRNA; the slide quotes protein-coding counts
pc = {s: sum(float(g[f"count_{s}"]) for g in genes.values() if g[f"count_{s}"] and g["gene_biotype"] == "protein_coding") / 1e6
      for s in PRI + REC}
raw = lambda g, s: float(genes[g][f"count_{s}"])  # noqa: E731
hum = {r["gene_name"]: r for r in csv.DictReader(open(HD / "human_side_check.tsv"), delimiter="\t")}
THEMES = {"blood": ["Hbb", "Hba-a1", "Alas2"], "fibroblast / matrix": ["Mmp13", "Fap", "Lrrc15"],
          "interferon / lymphocyte": ["Cxcl10", "Batf2", "Il2ra", "Slamf6", "Ighm"]}
by = {name(r["gene"]): r for r in sig}
for t, gs in THEMES.items():
    assert all(g in by for g in gs), (t, [g for g in gs if g not in by])
nl_only = [g for g in THEMES["interferon / lymphocyte"]
           if norm(by[g]["gene"], "IL66B") < 2 * max(norm(by[g]["gene"], s) for s in PRI)]   # IL66B barely above the primaries

facts = dict(
    tested=len(th), sig=len(sig), up=len(up), down=len(down), separating=sum(map(separates, sig)),
    sig2=len(sig2), up2=len(up2), down2=len(down2),
    down_all_three=sum(separates(r) for r in down), pc_m_min=round(min(pc.values()), 1), pc_m_max=round(max(pc.values()), 1),
    tumour_vs_control_sig=sum(r["padj"] < 0.05 for r in tc), tumour_vs_control_tested=len(tc),
    reads_m_min=round(min(reads.values()), 1), reads_m_max=round(max(reads.values()), 1),
    mmp13_fold=round(2 ** -by["Mmp13"]["lfc"], 1), mmp13_padj=by["Mmp13"]["padj"],
    human_MMP13_tumour_max=max(int(hum["MMP13"][s]) for s in PRI + REC),
    rat_Mmp13_primary_min=round(min(raw(by["Mmp13"]["gene"], s) for s in PRI)),          # raw reads, like the human side
    rat_Mmp13_fold_range=[round(min(norm(by["Mmp13"]["gene"], p) / norm(by["Mmp13"]["gene"], r) for p in PRI for r in REC)),
                          round(max(norm(by["Mmp13"]["gene"], p) / norm(by["Mmp13"]["gene"], r) for p in PRI for r in REC))],
    il66b_weak_immune=nl_only,
    themes={t: {g: round(by[g]["lfc"], 2) for g in gs} for t, gs in THEMES.items()},
)
# leave-one-tumour-out on the same host counts (ANALYSIS/holdout_separation/loo_separation.R, job 40468363): DESeq2 fits
# with one tumour removed; padj05 = genes at padj < 0.05. Human-side counterpart: loo_separation_human.tsv.
loo = {r["held_out"]: int(r["padj05"]) for r in csv.DictReader(open(HD / "loo_separation_host.tsv"), delimiter="\t")}
loo_h = {r["held_out"]: int(r["padj05"]) for r in csv.DictReader(
    open(ROOT / "ANALYSIS" / "holdout_separation" / "loo_separation_human.tsv"), delimiter="\t")}
facts["loo_host"] = loo
facts["loo_human"] = loo_h
# The pipeline's own Broad GSEA for host_therapy (it finished before the plotting crash; rescued from the work dir)
BG = ROOT / "ANALYSIS" / "host" / "verify" / "broad_gsea_existing"
gs = {}
for side, f in (("up", "host_therapy.combined_rat.gsea_report_for_Recurrent_U2.tsv"),
                ("down", "host_therapy.combined_rat.gsea_report_for_Primary_U2.tsv")):
    rows = list(csv.DictReader(open(BG / f), delimiter="\t"))
    qs = [float(r["FDR q-val"]) for r in rows if r["FDR q-val"] not in ("", "---")]
    gs[side] = dict(n=len(rows), min_fdr=min(qs), n_fdr25=sum(q < 0.25 for q in qs),
                    top=[(r["NAME"], float(r["NES"]), float(r["NOM p-val"])) for r in rows[:3]])
facts["broad_gsea"] = gs
# Local rerun of the same DESeq2 fits (DESeq2 1.50.2; reproduces the pipeline's 38 genes exactly), ANALYSIS/host/verify/prelim_local
PL = ROOT / "ANALYSIS" / "host" / "verify" / "prelim_local"
splits = list(csv.DictReader(open(PL / "r2_relabel_splits.tsv"), delimiter="\t"))
cnt = sorted(((int(r["padj05"]), r["target"], r["true_split"] == "TRUE") for r in splits), reverse=True)
true_n = [n for n, _, t in cnt if t][0]
others = sorted(n for n, _, t in cnt if not t)
facts["relabel"] = dict(n_splits=len(cnt), true_n=true_n, true_rank=[t for _, _, t in cnt].index(True) + 1,
                        largest_n=cnt[0][0], largest_target=cnt[0][1], median_others=others[len(others) // 2])
lg = list(csv.DictReader(open(PL / "r1_loo_genes.tsv"), delimiter="\t"))
facts["loo_all_six"] = [(r["gene_name"], r["direction"]) for r in lg if r["n_loo_kept"] == "6"]
lc = list(csv.DictReader(open(PL / "r1_loo_counts.tsv"), delimiter="\t"))
kept = {r["held_out"]: int(r["of_pipeline_genes_kept_same_sign"]) for r in lc if r["held_out"] != "none"}
facts["loo_kept_of_38"] = kept
assert facts["relabel"]["true_n"] == 38 and gs["up"]["n_fdr25"] == 0 and gs["down"]["n_fdr25"] == 0
assert loo["none"] == 38 and max(loo, key=loo.get) == "IL66B" and loo_h["none"] != loo["none"]
assert facts["sig"] == 38 and facts["up"] == 15 and facts["down"] == 23 and facts["tested"] == 17345
assert (facts["sig2"], facts["up2"], facts["down2"]) == (33, 14, 19) and facts["down_all_three"] == len(down)
assert all(abs(by[g]["lfc"]) >= 1 for gs in THEMES.values() for g in gs)          # every named gene passes the two-fold rule
assert facts["human_MMP13_tumour_max"] <= 7 < facts["rat_Mmp13_primary_min"]
json.dump(facts, open(f4.OUT / "host.json", "w"), indent=1)
print(json.dumps(facts, indent=1))

# ---- volcano, as chart_volcano.png: FDR < 0.05 coloured, up orange, down navy; themed genes labelled.
# The y axis is capped at CAP so one gene (Mmp13, 1e-26) does not flatten the rest; off-scale points are triangles.
CAP = 7.0
x = np.array([r["lfc"] for r in th]); y = -np.log10(np.clip([r["padj"] for r in th], 1e-30, None))
s_ = np.array([r["padj"] < 0.05 and abs(r["lfc"]) >= 1 for r in th]); off = y > CAP
yc = np.where(off, CAP - 0.25, y)
fig, ax = plt.subplots(figsize=(6.6, 4.6))
bare(ax)
ax.scatter(x[~s_], yc[~s_], s=8, color=FAINT, lw=0, zorder=2)
ok = s_ & ~off
ax.scatter(x[ok], yc[ok], s=70, color=np.where(x[ok] > 0, ORANGE, NAVY), lw=0.8, edgecolor="white", zorder=4)
ax.scatter(x[s_ & off], yc[s_ & off], s=110, marker="^", color=np.where(x[s_ & off] > 0, ORANGE, NAVY), lw=0.8,
           edgecolor="white", zorder=4)
ax.axhline(-np.log10(0.05), color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
for v in (-1, 1):
    ax.axvline(v, color=GREY, lw=0.9, ls=(0, (4, 4)), zorder=1)
ax.set_xlim(-7.5, 7.5); ax.set_ylim(-0.3, CAP + 1.7)
ax.set_yticks(range(0, int(CAP) + 1))
ax.set_xlabel("log$_2$ fold change   (recurrent vs primary)", fontsize=16, color=INK, labelpad=6)
ax.set_ylabel("−log$_{10}$ adjusted P", fontsize=16, color=INK)
ax.xaxis.set_major_formatter(f4.FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
top = ax.get_ylim()[1]
ax.text(-7.1, top * 0.975, "%d down" % len(down2), color=NAVY, fontsize=19, fontweight="bold", ha="left", va="top")
ax.text(7.1, top * 0.975, "%d up" % len(up2), color=ORANGE, fontsize=19, fontweight="bold", ha="right", va="top")
ax.text(-7.1, -np.log10(0.05) + 0.08, "FDR 0.05", color=GREY, fontsize=13, ha="left", va="bottom")
assert f"{by['Mmp13']['padj']:.0e}" == "4e-27" and int(off[s_].sum()) == 1 and name(th[int(np.argmax(y))]["gene"]) == "Mmp13"
lab = [by[g] for gs in THEMES.values() for g in gs if g != "Hba-a1"]   # its only free slot sits beside Emcn; Hbb and Alas2 carry blood
lab.sort(key=lambda r: r["padj"])
ycap = lambda r: min(-np.log10(max(r["padj"], 1e-30)), CAP - 0.25)  # noqa: E731
obstacles = [(r["lfc"], ycap(r), 5.5) for r in sig]
place_labels(ax, [(r["lfc"], ycap(r), name(r["gene"]) + ("\nadj. P 4 × 10$^{-27}$" if name(r["gene"]) == "Mmp13" else "")) for r in lab],
             obstacles=obstacles, fs=12.5, color=INK, style="italic", radii=(11, 17, 26, 38, 52), leader_from=11)
save(fig, "fig_host_volcano.png")
