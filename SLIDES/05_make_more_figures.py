# -*- coding: utf-8 -*-
"""More results for the CNS 2026 Abstract 418 deck (Greg 2026-09-26: 19 slides is not enough), drawn only from material
that survives the September re-analysis (Molecular Biology Reports manuscript, MBR/ESM_1.xlsx, the pipeline outputs on
read from the repository result folders named in SLIDES/u251_paths.py). The May nine-panel figure is NOT used: its subtype and pathway panels are the
superseded marker-panel scoring.

  chart_sorting.png        the ten libraries and what fraction of their reads is human graft, rat host or ambiguous
                           (README sample table, generated from the xengsort logs)
  chart_trajectory_pca.png exploratory PCA of log2 CPM over the 500 most variable genes with the culture sample beside
                           the six tumours (salmon.merged.gene_counts.tsv); the culture sample is n = 1 and outside the
                           differential model, said on the slide
  chart_threshold_sweep.png how many genes clear FDR at each fold-change threshold (ESM_1 S4)
  chart_drug_scatter.png   every clinically available candidate: signature reversal against predicted barrier
                           permeability (ESM_1 S12), the MBR Figure 2A at slide size
  table_prior_art.png      the prior-art audit of the ranked candidates (ESM_1 S15)
  enplots_composite.png    the pipeline's own enrichment running-sum plots for the six leading sets

    python SLIDES/05_make_more_figures.py
"""
import importlib.util
import json
import sys
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "figures_cns"
sys.path.insert(0, str(HERE))
from u251_paths import COUNTS, DECONTAM_DIR, GSEA_DIR, METADATA, RUVSEQ_DIR, check, gsea_file  # noqa: E402,F401
check()
spec = importlib.util.spec_from_file_location("mk04", HERE / "04_make_cns_figures.py")
mk = importlib.util.module_from_spec(spec); spec.loader.exec_module(mk)
save, table_image, bare = mk.save, mk.table_image, mk.bare
NAVY, BLUE, PALE, RED, GOLD, GREY, INK, LIGHT, FAINT, GRIDC = mk.NAVY, mk.BLUE, mk.PALE, mk.RED, mk.GOLD, mk.GREY, mk.INK, mk.LIGHT, mk.FAINT, mk.GRID
MINUS = mk.MINUS
ESM = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx")


# ====================================================================== 1. the ten libraries
def sample_table():
    """The README's generated sample table (xengsort classification percentages), parsed, never retyped."""
    t = (ROOT / "README.md").read_text(encoding="utf-8")
    block = t[t.index("<!-- BEGIN sample-table"):t.index("<!-- END sample-table -->")]
    rows = []
    for line in block.splitlines():
        if line.startswith("|") and not line.startswith("| Sample") and not line.startswith("| ---"):
            c = [x.strip() for x in line.strip("|").split("|")]
            rows.append({"sample": c[0], "cohort": c[1], "env": c[2], "graft": float(c[3]), "host": float(c[4]), "role": c[5]})
    return pd.DataFrame(rows)


def fig_sorting():
    d = sample_table()
    order = {"In vitro Culture": 0, "Primary (Pre-LITT)": 1, "Recurrent (Post-LITT)": 2, "Control (failed graft)": 3, "Control (procedural)": 3}
    d["g"] = d.cohort.map(order); d = d.sort_values(["g", "sample"]).reset_index(drop=True)
    heads = {0: ("culture", GOLD), 1: ("primary", NAVY), 2: ("recurrent", BLUE), 3: ("rat brain controls", GREY)}
    fig, ax = plt.subplots(figsize=(10.6, 4.6))
    xs, x, last = [], 0.0, None
    for _, r in d.iterrows():
        if last is not None and r.g != last:
            x += 0.75
        xs.append(x); last = r.g; x += 1.12
    d["x"] = xs
    both = 100 - d.graft - d.host        # reads shared by both genomes (analysed with the human reads) or unresolved
    ax.bar(d.x, d.graft, width=0.78, color=[heads[g][1] for g in d.g], zorder=3, label="human")
    ax.bar(d.x, d.host, bottom=d.graft, width=0.78, color="#D9D9D9", zorder=3, label="rat")
    ax.bar(d.x, both, bottom=d.graft + d.host, width=0.78, color="#F2F2F2", edgecolor="#D9D9D9", lw=0.6, zorder=3, label="shared by both, or unresolved")
    for _, r in d.iterrows():
        if r.graft > 20:
            ax.text(r.x, r.graft - 2, f"{r.graft:.0f} %", ha="center", va="top", fontsize=13, fontweight="bold", color="white", zorder=5)
        else:
            ax.text(r.x, r.graft + 1.5, f"{r.graft:.1f} %", ha="center", va="bottom", fontsize=12, color=INK, zorder=5)
    for g, (name, col) in heads.items():
        sub = d[d.g == g]
        ax.text(sub.x.mean(), 104, name, ha="center", va="bottom", fontsize=14.5, fontweight="bold", color=col)
    ax.set_xticks(d.x); ax.set_xticklabels(d["sample"], fontsize=10.5, color=INK)
    ax.set_ylim(0, 116); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel("share of reads, %", fontsize=15, color=INK)
    bare(ax, grid="y"); ax.tick_params(axis="x", length=0)
    from matplotlib.patches import Patch
    ax.legend([Patch(color=NAVY), Patch(color="#D9D9D9"), Patch(facecolor="#F2F2F2", edgecolor="#D9D9D9")],
              ["human (coloured by group)", "rat", "shared by both, or unresolved"],
              loc="upper left", bbox_to_anchor=(0.0, -0.12), ncol=3, fontsize=12.5, frameon=False)
    fig.subplots_adjust(bottom=0.2)
    save(fig, "chart_sorting.png")
    meta = pd.read_csv(METADATA)[["sample", "both_pct"]]          # the shared ('both') share, in-vivo libraries only
    d.merge(meta, on="sample", how="left").to_csv(OUT / "chart_sorting.csv", index=False)


# ====================================================================== 2. trajectory PCA with the culture sample
def fig_trajectory():
    cnt = pd.read_csv(COUNTS, sep="\t")
    samples = ["C2B", "IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B"]
    X = cnt[samples].to_numpy(dtype=float)
    cpm = X / X.sum(0, keepdims=True) * 1e6
    keep = cpm.mean(1) >= 1
    L = np.log2(cpm[keep] + 1)
    var = L.var(1); top = np.argsort(var)[-500:]
    Z = L[top] - L[top].mean(1, keepdims=True)
    U, S, Vt = np.linalg.svd(Z.T, full_matrices=False)
    pcs = U * S; ev = S ** 2 / (S ** 2).sum() * 100
    fig, ax = plt.subplots(figsize=(6.0, 4.9))
    bare(ax)
    ax.axhline(0, color=GRIDC, lw=0.9); ax.axvline(0, color=GRIDC, lw=0.9)
    style = {"C2B": (GOLD, "D", "culture"), "IL67B": (NAVY, "^", "primary"), "IL68B": (NAVY, "^", "primary"), "IL69B": (NAVY, "^", "primary"),
             "IL66B": (BLUE, "s", "recurrent"), "NL70B": (BLUE, "s", "recurrent"), "NL71B": (BLUE, "s", "recurrent")}
    seen = set()
    for i, s in enumerate(samples):
        col, mk_, lab = style[s]
        ax.scatter(pcs[i, 0], pcs[i, 1], s=260, marker=mk_, color=col, edgecolor="white", linewidth=1.5, zorder=4, label=lab if lab not in seen else None)
        seen.add(lab)
    ax.set_xlabel(f"PC1  ({ev[0]:.1f} % of variance)", fontsize=15, color=INK)
    ax.set_ylabel(f"PC2  ({ev[1]:.1f} %)", fontsize=15, color=INK)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}"))); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: MINUS(f"{v:.0f}")))
    ax.legend(loc="best", fontsize=13, frameon=False)
    # arrows: culture -> primary centroid -> recurrent centroid
    c0 = pcs[0, :2]; c1 = pcs[1:4, :2].mean(0); c2 = pcs[4:7, :2].mean(0)
    for a, b in ((c0, c1), (c1, c2)):
        ax.annotate("", b, a, arrowprops=dict(arrowstyle="-|>", color=GREY, lw=2, mutation_scale=18, shrinkA=14, shrinkB=14), zorder=2)
    save(fig, "chart_trajectory_pca.png")
    pd.DataFrame({"sample": samples, "PC1": pcs[:, 0], "PC2": pcs[:, 1]}).assign(var1=ev[0], var2=ev[1]).to_csv(OUT / "chart_trajectory_pca.csv", index=False)
    Zs = Z.T                                                   # samples x 500 genes, centred
    cen = lambda idx: Zs[idx].mean(0)                          # noqa: E731
    d_dish = float(np.linalg.norm(Zs[0] - cen([1, 2, 3, 4, 5, 6])))
    d_arms = float(np.linalg.norm(cen([1, 2, 3]) - cen([4, 5, 6])))
    json.dump({"n_genes_cpm_ge1": int(keep.sum()), "pc1": float(ev[0]), "pc2": float(ev[1]), "dist_culture_to_tumour_centroid": d_dish,
               "dist_primary_to_recurrent_centroid": d_arms, "ratio": d_dish / d_arms, "pc1_culture": float(pcs[0, 0]),
               "pc1_primary_mean": float(pcs[1:4, 0].mean()), "pc1_recurrent_mean": float(pcs[4:7, 0].mean())},
              open(OUT / "chart_trajectory_pca.json", "w"), indent=1)
    print(f"  distances in the 500-gene space: culture to tumours {d_dish:.1f}, primary to recurrent {d_arms:.1f}, ratio {d_dish / d_arms:.1f}")
    print("  trajectory PCA: genes", int(keep.sum()), "top-500 var; PC1", round(ev[0], 1), "PC2", round(ev[1], 1))


# ====================================================================== 3. the threshold sweep
def fig_threshold():
    s4 = ESM.parse("S4_DE_threshold_sweep")
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    x = np.arange(len(s4))
    ax.bar(x - 0.18, s4["Up in recurrence"], width=0.34, color=RED, zorder=3, label="up in recurrence")
    ax.bar(x + 0.18, s4["Down in recurrence"], width=0.34, color=NAVY, zorder=3, label="down")
    for xi, (u, dn) in enumerate(zip(s4["Up in recurrence"], s4["Down in recurrence"])):
        ax.text(xi - 0.18, u + 1, str(u), ha="center", va="bottom", fontsize=11, color=INK)
        ax.text(xi + 0.18, dn + 1, str(dn), ha="center", va="bottom", fontsize=11, color=INK)
    ax.set_xticks(x); ax.set_xticklabels([f"> {v:g}×" for v in (1.5, 2, 2.83, 4)], fontsize=11.5, color=INK)
    ax.set_xlabel("fold-change threshold (FDR < 0.05 throughout)", fontsize=11.5, color=INK)
    ax.set_ylabel("genes", fontsize=11.5, color=INK); ax.set_ylim(0, 70)
    bare(ax, grid="y"); ax.tick_params(axis="x", length=0, labelsize=11)
    ax.legend(loc="upper right", fontsize=10.5, frameon=False)
    save(fig, "chart_threshold_sweep.png")
    json.dump({"fold": [1.5, 2, 2.83, 4], "up": [int(v) for v in s4["Up in recurrence"]], "down": [int(v) for v in s4["Down in recurrence"]]},
              open(OUT / "chart_threshold_sweep.json", "w"), indent=1)


# ====================================================================== 4. the drug scatter
def fig_drug_scatter():
    d = ESM.parse("S12_Drug_ranking_full")
    d = d[d["Reached a clinic"].astype(str).str.lower() == "yes"].copy()
    d["absnes"] = d["NES"].abs()
    d = d.sort_values("Integrated score (|NES|^1.5 x ADMET-AI BBB)", ascending=False).drop_duplicates("Drug")
    agree = d["Both BBB models agree"].astype(str).str.lower() == "yes"
    sc = d["Integrated score (|NES|^1.5 x ADMET-AI BBB)"]; lo, hi = sc.min(), sc.max()
    sz = 40 + 380 * (sc - lo) / (hi - lo)
    fig, ax = plt.subplots(figsize=(7.2, 4.9))
    bare(ax)
    top = d["Drug"].iloc[0]
    for m, col, lab in ((agree & (d.Drug != top), BLUE, "both models predict crossing"), (~agree, "white", "at least one predicts no crossing")):
        ax.scatter(d.absnes[m], d["ADMET-AI BBB probability"][m], s=sz[m], color=col, edgecolor=BLUE, linewidth=1.6, zorder=3, label=lab, alpha=0.95)
    t = d[d.Drug == top]
    ax.scatter(t.absnes, t["ADMET-AI BBB probability"], s=sz[t.index], color=NAVY, edgecolor="white", linewidth=1.6, zorder=5, label=top)
    ax.axhline(0.5, color=GREY, lw=1.2, ls=(0, (4, 3)), zorder=1)
    ax.text(d.absnes.max() + 0.11, 0.515, "threshold 0.5", fontsize=11.5, color=GREY, ha="right", va="bottom")
    ax.set_xlabel("|NES|  (how strongly the drug signature opposes recurrence)", fontsize=13.5, color=INK)
    ax.set_ylabel("predicted barrier permeability", fontsize=13.5, color=INK)
    ax.set_ylim(-0.03, 1.1); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlim(d.absnes.min() - 0.08, d.absnes.max() + 0.12)
    # labels placed clear of every marker and of each other (radius from the marker area, in points)
    obstacles = [(float(a), float(b), float(np.sqrt(s_ / np.pi)) + 1.5) for a, b, s_ in zip(d.absnes, d["ADMET-AI BBB probability"], sz)]
    items = [(float(t.absnes.iloc[0]), float(t["ADMET-AI BBB probability"].iloc[0]), top, {"fontsize": 15, "fontweight": "bold", "color": NAVY})]
    for name in ("pentetrazol", "metformin", "progesterone"):      # the crowded top row is named in the slide text instead
        r = d[d.Drug.str.lower() == name]
        if len(r):
            items.append((float(r.absnes.iloc[0]), float(r["ADMET-AI BBB probability"].iloc[0]), name, {"fontsize": 11.5, "color": INK}))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3, fontsize=12, frameon=False, scatterpoints=1, markerscale=0.5, columnspacing=1.6)
    fig.subplots_adjust(bottom=0.24)                    # final axes geometry first, then place the labels in it
    mk.place_labels(ax, items, obstacles=obstacles, radii=(8, 14, 22, 32, 44), leader_from=14, skip=True)
    save(fig, "chart_drug_scatter.png")
    allr = ESM.parse("S12_Drug_ranking_full").assign(absnes=lambda t: t["NES"].abs()).sort_values("absnes", ascending=False).drop_duplicates("Drug")
    strongest = [{"drug": str(r.Drug), "absnes": round(float(r.absnes), 3), "clinic": str(r["Reached a clinic"]).lower() == "yes",
                  "phase": None if pd.isna(r["Max clinical phase"]) else float(r["Max clinical phase"])} for _, r in allr.head(6).iterrows()]
    (OUT / "chart_drug_scatter.json").write_text(json.dumps({"n_with_phase": int(len(d)), "top": top, "strongest_overall": strongest}, indent=1),
                                                 encoding="utf-8")


# ====================================================================== 5. the prior-art audit
SHORT_EV = {"ganciclovir": "Phase 3, only as an HSV-tk prodrug or for CMV",          # S15 text, shortened to fit the column
            "amiodarone": "In vitro + murine xenograft; worse survival in patients"}


def table_prior_art():
    """The audited compounds are S13's top 20 (the twenty highest-ranked of the 54 clinically available compounds); the
    audit itself is S15. Joined on the compound name, because S15's own 'In the ranked list' flag misses thioridazine
    (rank 12 in S12 and S13)."""
    s13 = ESM.parse("S13_Drug_candidates_top20")
    s15 = ESM.parse("S15_Prior_art")
    s15 = s15.assign(key=s15.Compound.astype(str).str.lower()).drop_duplicates("key").set_index("key")
    tier_word = {"A": "A: in vivo or clinical glioma evidence", "B": "B: in vitro or contested", "C": "C: none, failed, or not a therapy"}
    rows, tiers = [], {}
    for _, r in s13.sort_values("Rank (BBB-weighted)").iterrows():
        k = str(r.Drug).lower()
        assert k in s15.index, f"{k} has no prior-art row"
        a = s15.loc[k]
        ev = a["Highest level of glioma evidence"]
        ev = "none found" if pd.isna(ev) else re.sub(r"\s+", " ", str(ev))
        ev = SHORT_EV.get(k, ev)
        tier = str(a["Evidence tier"])
        tiers.setdefault(tier, []).append((int(r["Rank (BBB-weighted)"]), k))
        rows.append([k, f"{int(r['Rank (BBB-weighted)'])}", tier, (ev[:58] + "…") if len(ev) > 60 else ev])
    assert len(rows) == 20 and [int(x[1]) for x in rows] == list(range(1, 21)), [x[1] for x in rows]
    table_image("table_prior_art.png", ["compound", "rank", "tier", "highest glioma evidence"], rows, [2.2, 0.8, 0.7, 6.3], 10.0, fs=13.5, row_h=0.42,
                align=["left", "center", "center", "left"], hl={0})
    (OUT / "table_prior_art_tiers.txt").write_text("\n".join(tier_word.values()), encoding="utf-8")
    json.dump(tiers, open(OUT / "table_prior_art.json", "w"), indent=1)
    print("  prior-art tiers:", {t: [n for _, n in v] for t, v in sorted(tiers.items())})


# ====================================================================== 6. the enrichment running-sum plots
def enplot_composite():
    names = [("KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "translation initiation  (q = 0.022)"), ("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "translation elongation"),
             ("REACTOME_RESPONSE_OF_EIF2AK4_GCN2", "GCN2 amino-acid stress response"), ("KEGG_RIBOSOME", "ribosome"),
             ("REACTOME_SELENOAMINO_ACID_METABOLISM", "selenoamino-acid metabolism"), ("REACTOME_CELLULAR_RESPONSE_TO_STARVATION", "response to starvation")]
    tiles = []
    for key, lab in names:
        p = sorted(GSEA_DIR.glob(f"*enplot_{key}*.png"))
        if not p:
            print("  missing enplot", key); continue
        im = Image.open(p[0]).convert("RGB")
        tiles.append((im, lab))
    w, h = tiles[0][0].size
    F = ImageFont.truetype(r"C:\Windows\Fonts\arialbd.ttf", int(h * 0.075))
    pad = int(h * 0.14)
    grid = Image.new("RGB", (3 * w + 2 * 20, 2 * (h + pad) + 20), "white")
    for i, (im, lab) in enumerate(tiles):
        r, c = divmod(i, 3)
        x, y = c * (w + 20), r * (h + pad + 20)
        d = ImageDraw.Draw(grid)
        tw = d.textlength(lab, font=F)
        d.text((x + (w - tw) / 2, y + int(pad * 0.15)), lab, fill=(0x19, 0x1D, 0x63), font=F)
        grid.paste(im, (x, y + pad))
    grid.save(OUT / "enplots_composite.png"); print("  enplots_composite.png", grid.size)


# ====================================================================== 7. the running-sum plots, redrawn at slide size
def fig_running_sum():
    """The GSEA running enrichment score for the six leading sets, from the pipeline's per-set tables (rank of every member
    gene in the ordered list and the running ES at that rank); the curve between hits is the linear miss penalty."""
    rep = pd.read_csv(gsea_file("gsea_report_for_Primary_U2.tsv"), sep="\t")
    ranked = pd.read_csv(gsea_file("ranked_gene_list_Recurrent_U2_versus_Primary_U2.tsv"), sep="\t")
    pos = {g: i for i, g in enumerate(ranked["NAME"].astype(str))}          # 0-based position in the ranked list
    N_RANKED = int(len(ranked))
    names = [("KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "translation initiation"), ("REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION", "translation elongation"),
             ("REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "GCN2 amino-acid stress"), ("KEGG_RIBOSOME", "ribosome"),
             ("REACTOME_SELENOAMINO_ACID_METABOLISM", "selenoamino-acid metabolism"), ("REACTOME_CELLULAR_RESPONSE_TO_STARVATION", "response to starvation")]
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 5.6), sharex=True, sharey=True)
    n_ranked = None
    for ax, (key, lab) in zip(axes.ravel(), names):
        t = pd.read_csv(gsea_file(f"{key}.tsv"), sep="\t").sort_values("RANK IN GENE LIST")
        r = t["RANK IN GENE LIST"].to_numpy(); es = t["RUNNING ES"].to_numpy()
        n_ranked = N_RANKED
        base0 = all(pos.get(s) == rk for s, rk in zip(t["SYMBOL"].astype(str), t["RANK IN GENE LIST"]) if s in pos)
        assert base0, "set-table ranks are not 0-based positions in the ranked list"
        xs = np.concatenate([[0], r, [n_ranked]]); ys = np.concatenate([[0], es, [0]])
        row = rep[rep.NAME == key].iloc[0]
        ax.fill_between(xs, ys, 0, color=PALE, alpha=0.7, zorder=1)
        ax.plot(xs, ys, color=NAVY, lw=2.2, zorder=3)
        ax.vlines(r, -1.28, -1.08, color=INK, lw=0.7, zorder=2)
        ax.axhline(0, color=GRIDC, lw=1)
        ax.set_title(lab, fontsize=13.5, color=INK, fontweight="bold", pad=6)
        ax.text(0.03, 0.30, f"NES {MINUS(f'{row.NES:.2f}')}   {int(row.SIZE)} genes\nq = {row['FDR q-val']:.3f}" if row["FDR q-val"] < 0.05 else
                f"NES {MINUS(f'{row.NES:.2f}')}   {int(row.SIZE)} genes\nq = {row['FDR q-val']:.2f}", transform=ax.transAxes, fontsize=12, color=INK,
                ha="left", va="bottom", linespacing=1.2, fontweight="bold" if row["FDR q-val"] < 0.05 else "normal",
                bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.9), zorder=6)
        for sd in ("top", "right"):
            ax.spines[sd].set_visible(False)
        ax.tick_params(labelsize=10.5)
    for ax in axes[1]:
        ax.set_xlabel("rank in the list, up in recurrence → down", fontsize=11.5, color=INK)
        ax.set_xticks([0, 5000, 10000, 15000, n_ranked - 1]); ax.set_xticklabels(["1", "5,000", "10,000", "15,000", f"{n_ranked:,}"])
    for ax in axes[:, 0]:
        ax.set_ylabel("running enrichment score", fontsize=11.5, color=INK)
    axes[0, 0].set_ylim(-1.3, 0.15)
    fig.tight_layout()
    save(fig, "chart_running_sum.png")
    out = {}
    for key, lab in names:
        t = pd.read_csv(gsea_file(f"{key}.tsv"), sep="\t")
        row = rep[rep.NAME == key].iloc[0]
        first0 = int(t["RANK IN GENE LIST"].min())                    # 0-based, checked in the drawing loop
        rk = t["RANK IN GENE LIST"].to_numpy()
        out[key] = {"label": lab, "n": int(len(t)), "first_rank0": first0, "n_ranked": n_ranked, "n_from_first_to_end": n_ranked - first0,
                    "share_in_bottom_tenth": float((rk >= 0.9 * n_ranked).mean()), "share_in_bottom_quarter": float((rk >= 0.75 * n_ranked).mean()),
                    "share_in_bottom_half": float((rk >= 0.5 * n_ranked).mean()),
                    "bottom_share": (n_ranked - first0) / n_ranked, "es_min": float(t["RUNNING ES"].min()), "nes": float(row.NES), "q": float(row["FDR q-val"]),
                    "metric": "Diff_of_Classes (difference of class means, linear scale); permute gene_set; nperm 1000"}
    json.dump(out, open(OUT / "chart_running_sum.json", "w"), indent=1)


if __name__ == "__main__":
    print("figures ->", OUT)
    fig_sorting(); fig_trajectory(); fig_threshold(); fig_drug_scatter(); table_prior_art(); enplot_composite(); fig_running_sum()
    print("done")
