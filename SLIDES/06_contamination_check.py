# -*- coding: utf-8 -*-
"""Is the translational shutdown rat contamination? A read of the repository's own contamination outputs (Phase 2B,
mirrored from the grid) for the six leading gene sets of the reported GSEA.

Definitions are the repository's, not new ones:
  * "detected in rat brain" = CPM >= 1 in >= 2 of the three rat-brain controls (IL64B, N168B, N269B), the rule of
    ANALYSIS/test_contamination_bias.R and of the absolute mode of filter_contaminated_genes.R; per gene in
    ANALYSIS/results_decontamination/contamination_bias_per_gene.csv (column ctrl_detected), over the 14,849 genes the
    unadjusted DESeq2 of the RUVSeq run tested.
  * contaminated_genes.tsv is NOT used: it is the ratio-mode output (mean control CPM / mean tumour CPM >= 0.5), which
    flags genes the tumour barely expresses, not genes that carry rat signal.
  * fold changes: ANALYSIS/results_ruvseq/ruvseq_{baseline_de,adjusted_de_k1,adjusted_de_k2}.tsv (DESeq2 ~ Classification,
    and with k = 1, 2 RUVs factors estimated from the controls); this unadjusted run is the RUVSeq script's own baseline
    (no ashr shrinkage), not the manuscript's DESeq2 table.
  * set members: the pipeline's per-set GSEA tables (work dir of the reported run), joined to the DESeq2 tables by
    Ensembl id through the baseline table's gene_id <-> symbol pairs.

The direction argument: recurrent tumours carry more rat than primaries, so misassigned rat reads on conserved genes
would push those genes UP in recurrence. Per set and transcriptome-wide: share of genes detected in rat brain; mean
log2FC of detected vs undetected members; the set's mean and share-negative before and after RUVs adjustment.
Also: the 35 differentially expressed genes of the manuscript (ESM_1 S3) against the same rule.
Writes figures_cns/chart_contamination.png and figures_cns/contamination_check.json.

    python SLIDES/06_contamination_check.py
"""
import importlib.util
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "figures_cns"
GRID = Path(r"C:\Users\grego\AppData\Local\Temp\claude\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\scratchpad\u251_grid")
spec = importlib.util.spec_from_file_location("mk04", HERE / "04_make_cns_figures.py")
mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
save = mk.save
NAVY, BLUE, PALE, RED, GOLD, GREY, INK, GRIDC = mk.NAVY, mk.BLUE, mk.PALE, mk.RED, mk.GOLD, mk.GREY, mk.INK, mk.GRID
MINUS = mk.MINUS

SETS = [("KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "translation\ninitiation"), ("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "translation\nelongation"),
        ("REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "GCN2 amino-\nacid stress"), ("KEGG_RIBOSOME", "ribosome"),
        ("REACTOME_SELENOAMINO_ACID_METABOLISM", "selenoamino-\nacid metabolism"), ("REACTOME_CELLULAR_RESPONSE_TO_STARVATION", "starvation\nresponse")]
TAB = {name: pd.read_csv(GRID / "ruvseq" / f, sep="\t").set_index("gene_id")
       for name, f in (("baseline", "ruvseq_baseline_de.tsv"), ("k1", "ruvseq_adjusted_de_k1.tsv"), ("k2", "ruvseq_adjusted_de_k2.tsv"))}
bias = pd.read_csv(GRID / "decontam" / "contamination_bias_per_gene.csv").set_index("gene_id")
det = bias["ctrl_detected"].astype(str).str.upper() == "TRUE"
sym2id = TAB["baseline"].reset_index()[["gene_id", "symbol"]].dropna()
sym2id = sym2id[sym2id.gene_id.isin(bias.index)]          # the tested genes only: one id per symbol there
dup = sym2id.symbol.duplicated(keep=False)
print(f"tested genes {len(bias):,}; detected in rat brain {int(det.sum()):,} ({100 * det.mean():.1f} %); symbols mapping to >1 tested id: {int(dup.sum())}")
sym2id = sym2id.drop_duplicates("symbol").set_index("symbol")["gene_id"]


def lfc(name, ids):
    return TAB[name].reindex(ids)["log2FoldChange"].dropna()


tested = bias.index
allg = {name: float(lfc(name, tested).mean()) for name in TAB}
wide = {"n_tested": int(len(bias)), "n_detected": int(det.sum()), "share_detected": float(det.mean()),
        "mean_lfc_detected": float(bias.loc[det, "log2FoldChange"].mean()), "mean_lfc_clean": float(bias.loc[~det, "log2FoldChange"].mean()),
        "median_lfc_detected": float(bias.loc[det, "log2FoldChange"].median()), "median_lfc_clean": float(bias.loc[~det, "log2FoldChange"].median())}
print("transcriptome:", {k: round(v, 3) if isinstance(v, float) else v for k, v in wide.items()}, "all-gene mean lfc", {k: round(v, 3) for k, v in allg.items()})

res = {}
for key, lab in SETS:
    t = pd.read_csv(GRID / "enplots" / f"therapy_impact.combined_human.{key}.tsv", sep="\t")
    syms = t["SYMBOL"].astype(str).tolist()
    ids = [sym2id[s] for s in syms if s in sym2id.index]
    d_ids = [i for i in ids if det.get(i, False)]; c_ids = [i for i in ids if not det.get(i, False)]
    r = {"n_members": len(syms), "n_tested": len(ids), "n_detected": len(d_ids), "share_detected": len(d_ids) / max(len(ids), 1),
         "mean_lfc_detected": float(lfc("baseline", d_ids).mean()) if d_ids else None,
         "mean_lfc_undetected": float(lfc("baseline", c_ids).mean()) if c_ids else None}
    for name in TAB:
        v = lfc(name, ids)
        r[name] = {"mean_lfc": float(v.mean()), "median_lfc": float(v.median()), "share_negative": float((v < 0).mean()), "n": int(len(v))}
    res[key] = r
    print(f"{key[:46]:46s} tested {len(ids):3d}/{len(syms):3d}  detected {len(d_ids):3d} ({100 * r['share_detected']:.0f} %)  "
          f"lfc det {r['mean_lfc_detected'] if d_ids else float('nan'):+.2f} undet {r['mean_lfc_undetected'] if c_ids else float('nan'):+.2f}  "
          f"base {r['baseline']['mean_lfc']:+.2f} k1 {r['k1']['mean_lfc']:+.2f} k2 {r['k2']['mean_lfc']:+.2f}  neg base/k1/k2 "
          f"{r['baseline']['share_negative']:.2f}/{r['k1']['share_negative']:.2f}/{r['k2']['share_negative']:.2f}")

# the manuscript's 35 differentially expressed genes against the same rule
esm = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx").parse("S3_DE_significant")
idcol = [c for c in esm.columns if "gene" in c.lower() and "id" in c.lower()][0]
de_ids = esm[idcol].astype(str).tolist()
de_in = [i for i in de_ids if i in bias.index]
de_det = [i for i in de_in if det[i]]
de35 = {"n": len(de_ids), "n_in_tested": len(de_in), "n_detected": len(de_det)}
print("manuscript DE genes:", de35)

# the direction on the log scale in the MANUSCRIPT's DESeq2 table (ashr-shrunk, Online Resource 1 S2), per set:
# the GSEA ranking metric (Diff_of_Classes) is a linear-scale difference that puts abundant genes at the ends
s2 = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx").parse("S2_DE_all_genes")
symcol = [c for c in s2.columns if c.lower() in ("symbol", "gene symbol", "gene_name", "gene name")][0]
lfccol = [c for c in s2.columns if "log2" in c.lower()][0]
s2 = s2.dropna(subset=[lfccol]).drop_duplicates(symcol).set_index(symcol)
for key, _ in SETS:
    t = pd.read_csv(GRID / "enplots" / f"therapy_impact.combined_human.{key}.tsv", sep="\t")
    v = s2.reindex(t["SYMBOL"].astype(str))[lfccol].dropna()
    res[key]["manuscript_s2"] = {"n": int(len(v)), "share_negative": float((v < 0).mean()), "mean_lfc": float(v.mean())}
    print(f"  manuscript S2 {key[:40]:40s} n {len(v):3d}  share negative {100 * (v < 0).mean():5.1f} %  mean {v.mean():+.3f}")
conc = pd.read_csv(GRID / "ruvseq" / "ruvseq_concordance_summary.csv").to_dict(orient="records")
graft = pd.read_csv(GRID / "metadata_full.csv").to_dict(orient="records")

# ---- chart: per set, detected vs undetected members (left pair), then the adjusted means
fig, ax = plt.subplots(figsize=(10.0, 5.6))
xs = []
for i, (key, lab) in enumerate(SETS):
    r = res[key]; x0 = i * 4.2
    t = pd.read_csv(GRID / "enplots" / f"therapy_impact.combined_human.{key}.tsv", sep="\t")
    ids = [sym2id[s] for s in t["SYMBOL"].astype(str) if s in sym2id.index]
    for j, (name, col, lab2) in enumerate((("baseline", NAVY, "unadjusted"), ("k1", BLUE, "adjusted, k = 1"), ("k2", "#9CB3DC", "adjusted, k = 2"))):
        x = x0 + j
        v = lfc(name, ids)
        ax.bar(x, v.mean(), width=0.8, color=col, zorder=3, label=lab2 if i == 0 else None)
        g = allg[name]
        ax.plot([x - 0.45, x + 0.45], [g, g], color=RED, lw=2.4, zorder=6, solid_capstyle="round", label="all tested genes, mean" if (i == 0 and j == 0) else None)
        rng = np.random.default_rng(i * 3 + j)
        dmask = np.array([det.get(k_, False) for k_ in v.index])
        ax.scatter(x + rng.uniform(-0.28, 0.28, len(v)), v.to_numpy(), s=9, lw=0, zorder=4, alpha=0.45,
                   c=np.where(dmask, GOLD, INK), label=None)
        ax.text(x, min(v.mean(), float(np.percentile(v, 4))) - 0.05, MINUS(f"{v.mean():+.2f}"), ha="center", va="top", fontsize=10.5, fontweight="bold",
                color=INK, bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.85), zorder=7)
    ax.text(x0 + 1.0, 1.02, lab, ha="center", va="bottom", fontsize=12, color=INK, linespacing=1.1)
    ax.text(x0 + 1.0, -1.18, f"{r['n_detected']} of {r['n_tested']}\ndetected in\nrat brain", ha="center", va="top", fontsize=11, color=GREY, linespacing=1.1)
    xs.append(x0 + 1.0)
ax.scatter([], [], s=30, c=GOLD, label="gene detected in rat brain")
ax.axhline(0, color=INK, lw=1.2, zorder=5)
ax.set_xticks([]); ax.set_ylim(-1.75, 1.35); ax.set_ylabel("log$_2$ fold change, recurrent vs primary", fontsize=13.5, color=INK)
ax.set_yticks([-1.0, -0.5, 0, 0.5, 1.0]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.1f}")))
for sd in ("top", "right", "bottom"):
    ax.spines[sd].set_visible(False)
ax.grid(axis="y", color=GRIDC, lw=0.8, zorder=0); ax.set_axisbelow(True); ax.tick_params(axis="y", labelsize=12.5)
ax.legend(loc="lower left", fontsize=11.5, frameon=False, ncol=5, bbox_to_anchor=(-0.02, -0.17), columnspacing=1.2, handletextpad=0.5)
fig.subplots_adjust(bottom=0.14)
save(fig, "chart_contamination.png")
json.dump({"definition": "detected in rat brain = CPM >= 1 in >= 2 of 3 controls (contamination_bias_per_gene.csv ctrl_detected)",
           "sets": res, "transcriptome": wide, "all_genes_mean_lfc": allg, "manuscript_de35": de35, "concordance": conc, "graft_pct": graft},
          open(OUT / "contamination_check.json", "w"), indent=1)
print("wrote contamination_check.json")
