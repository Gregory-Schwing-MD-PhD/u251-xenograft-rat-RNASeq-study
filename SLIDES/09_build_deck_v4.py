# -*- coding: utf-8 -*-
"""The CNS 2026 Abstract 418 talk, version 4 (2026-09-26): 14 slides and 7 backups, one story, figures in the forms this
literature uses (Broad GSEA enrichment plot, leading-edge and GSVA heatmaps, PCA, volcano, correlation with a fitted
line and 95 % band, CIBERSORT stacked fractions, DepMap distributions), 16-18 pt slide text with the detail in the notes.

The story, as the robustness work left it:
  1. LITT kills the core; recurrence grows from the sublethal margin.
  2. In a human-in-rat model the tumour can be read on its own, but recurrent libraries are less human (a confound).
  3. Few single genes move; whole programs do. Translation initiation is the most depleted program in recurrence.
  4. Its direction survives leaving any tumour out; its FDR does not, and one primary (IL68B) carries its size.
  5. Scores track tumour content, but two tumours with the same human share still differ, and contamination is too small.
  6. Of ten published subtypes only the astrocyte-like score separates the arms, and it survives adjustment.
  7. Ciclopirox opposes the recurrent state; U-251 depends on the hypusination axis it inhibits; LITT opens the window.

Figures: 04_make_cns_figures.py, 05_make_more_figures.py, 07_contamination_magnitude.py, 08_robustness_figures.py,
10_standard_figures.py. Numbers: figures/facts.json, figures_cns/*.json, the supplement (MBR/ESM_1.xlsx, S12). Helpers,
title block and reference list are imported from 04_build_cns_template_deck.py.

    python SLIDES/09_build_deck_v4.py [--out PATH]   (default SLIDES/CNS2026_Schwing_Abstract418_CNStemplate_v4.pptx)
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.oxml.ns import qn
from pptx.util import Emu, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import MSO_ANCHOR

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location("deck3", HERE / "04_build_cns_template_deck.py")
d3 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d3)
from cns_template import (BLUE, BOTTOM, CAP, GAP, GREY, L, NAVY, R, TEMPLATE, TOP, WHITE,  # noqa: E402
                          add_text, content_slide, fit, navy_bar, notes, picture, refs)

FIG = HERE / "figures_cns"
J = lambda name: json.load(open(FIG / name, encoding="utf-8"))  # noqa: E731
m, pv, fig_and_text, fig_top, head_box, flat = d3.m, d3.pv, d3.fig_and_text, d3.fig_top, d3.head_box, d3.flat
F = d3.F
d3.REFS.update({
    "subramanian": "Subramanian A, et al. Proc Natl Acad Sci USA 2005;102:15545–15550.",
    "hanzelmann": "Hänzelmann S, et al. BMC Bioinformatics 2013;14:7.",
    "newman": "Newman AM, et al. Nat Methods 2015;12:453–457.",
    "love": "Love MI, et al. Genome Biol 2014;15:550.",
    "varn": "Varn FS, et al. Cell 2022;185:2184–2199.",
})
C = d3.Refs()
BLACK = RGBColor(0, 0, 0)
DESIGN = d3.DESIGN
TXT = 17          # body text on the slides


def table(s, left, top, widths, rows, size=14, bold_rows=(), row_h=330000):
    """A plain black table (no theme fill): header row bold, a thin grid."""
    shp = s.shapes.add_table(len(rows), len(rows[0]), Emu(left), Emu(top), Emu(sum(widths)), Emu(row_h * len(rows)))
    tbl = shp.table
    tblPr = tbl._tbl.tblPr
    tblPr.set("firstRow", "0"); tblPr.set("bandRow", "0")
    sid = tblPr.find(qn("a:tableStyleId"))
    if sid is None:
        sid = tblPr.makeelement(qn("a:tableStyleId"), {}); tblPr.append(sid)
    sid.text = "{5940675A-B579-460E-94D1-54222C63F5DA}"          # No Style, Table Grid
    for j, w in enumerate(widths):
        tbl.columns[j].width = Emu(w)
    for i, row in enumerate(rows):
        tbl.rows[i].height = Emu(row_h)
        for j, val in enumerate(row):
            cell = tbl.cell(i, j)
            cell.margin_left = cell.margin_right = Emu(60000); cell.margin_top = cell.margin_bottom = Emu(20000)
            tf = cell.text_frame; tf.word_wrap = True; cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = tf.paragraphs[0]; p.text = ""
            r = p.add_run(); r.text = str(val)
            r.font.size = Pt(size); r.font.color.rgb = BLACK
            r.font.bold = (i == 0) or (i in bold_rows)
    return shp


def build(out: Path):
    prs = Presentation(str(TEMPLATE))
    sld = prs.slides._sldIdLst
    RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    for e in list(sld):
        sld.remove(e); prs.part.drop_rel(e.get(RID))

    # ---------------------------------------------------------------- numbers, read once
    srt = list(csv.DictReader(open(FIG / "chart_sorting.csv", encoding="utf-8")))
    grp = lambda coh, col="graft": [float(r[col]) for r in srt if r["cohort"].startswith(coh)]  # noqa: E731
    pri, rec = grp("Primary"), grp("Recurrent")
    assigned = [float(r["graft"]) + float(r["host"]) for r in srt if r["cohort"].startswith(("Primary", "Recurrent"))]
    st = J("standard.json"); rb = J("robustness.json"); mag = J("contamination_magnitude.json"); mag3 = J("contamination_magnitude_ctrl3.json");
    fun = J("chart_funnel.json"); s13 = J("s13_top20.json"); sets_p = J("chart_gsea_sets.json")
    g, lo, gc, sub, cs, dm = st["gsea"], st["loo"], st["graft_corr"], st["subtypes"], st["cibersort"], st["depmap"]
    gr, ho = rb["graft"], rb["holdout"]
    le = st["leading_edge"]
    pvn = F["permanova"]
    TIK = "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"
    assert le["highest_tumour"]["IL68B"] == le["n_genes"]
    assert lo["runs"] == 41 and lo["nom_p_max"] < 0.001
    assert lo["all_seeds_q05"] == ["Without IL68B", "Without NL70B"] and set(lo["no_seed_q05"]) == {"Without IL67B", "Without IL69B", "Without IL66B"}
    ac = sub["Neftel_AC"]; mes1 = sub["Neftel_MES1"]
    assert ac["p"] < 0.05 and all(v["p"] >= 0.05 for k, v in sub.items() if k != "Neftel_AC")
    pair = gr["pair"]
    ini = mag["sets"]["initiation"]; ini3 = mag3["sets"]["initiation"]
    assert mag["controls"] == ["IL64B", "N168B"] and "N269B" in mag3["controls"]   # tumour-free controls (SAMPLE_KEY.md)
    fitr = [cs[k]["Correlation"] for k in ("IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B")]
    cpx = s13[0]; nxt = [x for x in s13[1:] if x["both_agree"]][0]
    assert cpx["drug"] == "ciclopirox" and cpx["rank"] == 1 and cpx["rank_unweighted"] == 1
    assert ho["IL68B"]["cic_rank"] == 1 and ho["IL66B"]["cic_rank"] == 3
    dg = dm["genes"]
    esm = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx").parse("S12_Drug_ranking_full")
    esm = esm[esm["Reached a clinic"].astype(str).str.lower() == "yes"].sort_values(
        "Integrated score (|NES|^1.5 x ADMET-AI BBB)", ascending=False).drop_duplicates("Drug").reset_index(drop=True)
    assert len(esm) == fun["n_clinical"] and esm.Drug.iloc[0] == "ciclopirox"

    # ---- 1 title (the v3 title block)
    s = prs.slides.add_slide(prs.slide_layouts[0])
    ph = {p.placeholder_format.idx: p for p in s.placeholders}
    t = ph[0]
    t.left, t.top, t.width, t.height = Emu(5185558), Emu(950000), Emu(6600000), Emu(2000000)
    t.text_frame.word_wrap = True; t.text_frame.text = d3.TITLE
    for p in t.text_frame.paragraphs:
        for r in p.runs:
            r.font.size = Pt(24); r.font.bold = True; r.font.color.rgb = WHITE
    sb = ph[1]
    sb.left, sb.top, sb.width, sb.height = Emu(5185558), Emu(3050000), Emu(6600000), Emu(1350000)
    sb.text_frame.word_wrap = True; sb.text_frame.text = d3.AUTHORS[0]
    for a in d3.AUTHORS[1:] + d3.AFFIL:
        p = sb.text_frame.add_paragraph(); p.text = a
    for i, p in enumerate(sb.text_frame.paragraphs):
        p.space_after = Pt(1 if i == 0 else 5 if i == 1 else 1); p.level = 0; flat(p)
        for r in p.runs:
            r.font.size = Pt(15 if i < 2 else 11); r.font.color.rgb = WHITE
    body = ph[10]
    body.left, body.top, body.width, body.height = Emu(5185558), Emu(4550000), Emu(6600000), Emu(1450000)
    body.text_frame.word_wrap = True; body.text_frame.text = "CNS 2026 Annual Meeting, Washington, DC"
    for line in d3.SESSION:
        p = body.text_frame.add_paragraph(); p.text = line
    for p in body.text_frame.paragraphs:
        p.space_after = Pt(3); p.level = 0; flat(p)
        for r in p.runs:
            r.font.size = Pt(12); r.font.color.rgb = WHITE
    notes(s, [
        "0:00–0:10  Title. Name, institutions, a rat model of MRI-guided LITT. Do not read the title.",
        "Disclosures, read aloud (CNS displays the list; no disclosure slide per the ACCME speaker guidelines): I.Y.L. has consulting "
        "agreements with Medtronic, Inc. (Minneapolis, MN) and Monteris Medical, Inc. (Plymouth, MN). All other authors declare no "
        "conflicts of interest. Supported by a Henry Ford Health Physician Scientist Award A20050 (I.Y.L.). No commercial support was "
        "received for this analysis. Drug candidates discussed here are computational predictions. None is approved for glioblastoma "
        "and none has been tested in this model.",
        "Timing: 14 slides in six minutes, about 25 s each (slide 3, the design, adds 20 s: everything after it runs 10 s later, taken back on slide 5); seven backups follow the acknowledgments (B7, the rat host, is preliminary). If running long, shorten slide 5 "
        "(PCA and volcano) to one sentence."])
    C.line()

    # ---- numbers added after the adversarial check (ANALYSIS files read directly)
    A_ = ROOT / "ANALYSIS"
    psd = pd.read_csv(A_ / "graft_relation/results/per_sample.tsv", sep="\t").set_index("sample")
    PRI_, REC_ = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]
    ac_sep = bool(psd.loc[REC_, "Neftel_AC"].max() < psd.loc[PRI_, "Neftel_AC"].min())
    assert ac_sep
    both_ = [float(r["both_pct"]) for r in srt if r["cohort"].startswith(("Primary", "Recurrent"))]
    unres = [100 - float(r["graft"]) - float(r["host"]) - float(r["both_pct"]) for r in srt if r["cohort"].startswith(("Primary", "Recurrent"))]
    loo_t = pd.read_csv(A_ / "gsea_leave_one_out/loo_sets.tsv", sep="\t")
    loo_t = loo_t[loo_t.set == TIK]
    filt = loo_t[~loo_t.condition.str.endswith("_nofilter")]; nof = loo_t[loo_t.condition.str.endswith("_nofilter")]
    assert len(filt) == 35 and len(nof) == 6
    q66_nof = float(nof[nof.condition == "drop_IL66B_nofilter"].fdr_q.iloc[0])
    pca = st["pca"]; sdist = st["sample_dist"]
    s12 = pd.ExcelFile(ROOT / "MBR" / "ESM_1.xlsx").parse("S12_Drug_ranking_full")
    n_fdr05 = int((s12["FDR q"].astype(float) < 0.05).sum())
    dfo = s12[s12.Drug.astype(str).str.lower() == "deferoxamine"]
    dfo_nes = float(dfo.NES.astype(float).min()); dfo_unw = abs(dfo_nes) ** 1.5
    cpx_unw = float(esm.loc[esm.Drug == "ciclopirox", "Score without BBB weighting (|NES|^1.5)"].iloc[0])
    assert dfo_unw > cpx_unw and str(dfo["Reached a clinic"].iloc[0]).lower() == "no"
    from scipy import stats as _st
    runs = {}
    for f_, cols_ in (("neftel4_confident", ["MES", "AC", "OPC", "NPC"]), ("neftel4_allcells", ["MES", "AC", "OPC", "NPC"]),
                      ("neftel4_confident_plus_nonmalignant", ["MES", "AC", "OPC", "NPC"]),
                      ("varn2022_glass_tumor3_noRMRP", ["stemcell_tumor", "prolif_stemcell_tumor", "differentiated_tumor"]),
                      ("varn2022_glass_sc12_noRMRP", ["stemcell_tumor", "prolif_stemcell_tumor", "differentiated_tumor"])):
        fr = pd.read_csv(A_ / f"cibersort/results/v104/fractions_v104_{f_}.tsv", sep="\t").set_index("sample")
        runs[f_] = {c_: float(_st.ttest_ind(fr.loc[PRI_, c_], fr.loc[REC_, c_]).pvalue) for c_ in cols_ if fr.loc[PRI_ + REC_, c_].std() > 0}
    p_state_min = min(min(v_.values()) for v_ in runs.values())
    glass_p = [runs[k_]["prolif_stemcell_tumor"] for k_ in ("varn2022_glass_tumor3_noRMRP", "varn2022_glass_sc12_noRMRP")]
    assert p_state_min >= 0.05
    for k_ in ("varn2022_glass_tumor3_noRMRP", "varn2022_glass_sc12_noRMRP"):
        fr = pd.read_csv(A_ / f"cibersort/results/v104/fractions_v104_{k_}.tsv", sep="	").set_index("sample")
        assert fr.loc[REC_, "prolif_stemcell_tumor"].min() > fr.loc[PRI_, "prolif_stemcell_tumor"].max(), k_

    # ---- 2 the problem
    s = content_slide(prs, "LITT kills the core; glioblastoma returns from the margin")
    img = FIG / "fig_margin.png"
    w, h = fit(img, 6200000, BOTTOM - TOP - 1000000)
    picture(s, img, L, TOP, w, h)
    add_text(s, L + w + GAP, TOP + 150000, R - L - w - GAP, h, [
        f"LITT is increasingly used for deep and recurrent glioblastoma.{C('chen')}",
        "Past the ablative threshold lies a sublethal margin.",
        "Recurrence is thought to grow from the cells that survive there."], size=18, gap=14)
    navy_bar(s, L, TOP + h + 150000, R - L, BOTTOM - (TOP + h + 150000),
             "Does the tumor that regrows after LITT carry its own program, and does that program point to a drug?", size=19)
    refs(s, C.line("Diagram: thermal dose falls with distance from the fiber; schematic, not to scale."))
    notes(s, ["0:10–0:35  The one idea: LITT kills the core, but the cells that come back sat in the sublethal margin. A first-pass analysis of these "
              "libraries (Nagaraja 2026, next slide) found cell-cycle, motility and inflammatory genes up in recurrence. What is new here: the reads "
              "are sorted from the host, read gene-set-wide, stress-tested, and taken to a drug. End on the question."])

    # ---- 3 the design (SLIDES/12_design_figure.py; lab sample key and array sheet, ANALYSIS/SAMPLE_KEY.md)
    s = content_slide(prs, "The design: eight rats, ten RNA libraries")
    img = FIG / "fig_design.png"
    w, h = fit(img, R - L, 3700000)
    picture(s, img, L + (R - L - w) // 2, TOP, w, h)
    add_text(s, L, TOP + h + 100000, R - L, BOTTOM - (TOP + h + 100000), [
        "Three primaries (rats 67, 68, 69) are compared with three recurrences (66, 70, 71). Rat 64's sample held almost no "
        "tumor and rat 65 has no library.",
        "The two brain controls are the opposite hemispheres of rats 68 and 69, not separate animals."], size=16, gap=6)
    refs(s, C.line("Nagaraja lab sample key (2023) and methylation-array sample sheet. The sequencing files call rats 70 and 71 "
                   "NL70B and NL71B; the lab records call them IL-70 and IL-71."))
    notes(s, ["0:35–0:55  The design, because it is easy to misread. Eight rats, all implanted. Four primaries, never ablated; four "
              "ablated with LITT and allowed to regrow. We compare three against three: rat 64's primary sample had almost no tumor in it, "
              "and there is no library for rat 65. The two brain controls come from the other side of the brain of two primary rats, "
              "68 and 69, so they are not independent animals. The culture is a reference only.",
              "Why 3 v 3 when Nagaraja 2026 reports 4 v 4: rat 64 (no tumor in the sample) and rat 65 (no library). "
              "Harvest time after LITT (2 or 4 weeks) for each recurrence is not in the records we have."])

    # ---- 4 the model and the read sorting
    s = content_slide(prs, "A human tumor in a rat brain, so the tumor can be read on its own")
    vs = 3000000
    s.shapes.add_movie(str(FIG / "ablation_512.mp4"), Emu(L), Emu(TOP), Emu(vs), Emu(vs),
                       poster_frame_image=str(FIG / "video_poster.png"), mime_type="video/mp4")
    add_text(s, L, TOP + vs + 60000, vs, 600000, [f"MRI-guided LITT in this model, diffusion-weighted MRI (from Nagaraja 2021{C('nag21')}); click to play."],
             size=CAP, color=GREY, gap=0)
    img = FIG / "chart_sorting.png"
    x2 = L + vs + GAP
    w, h = fit(img, R - x2, 3000000)
    picture(s, img, x2, TOP, w, h)
    add_text(s, x2, TOP + h + 150000, R - x2, BOTTOM - (TOP + h + 150000), [
        f"Reads sorted by species{C('xengsort')}: {min(assigned):.0f}–{max(assigned):.0f} % go to one species, {min(both_):.0f}–{max(both_):.0f} % are shared by both "
        f"and kept with the human reads, under {max(unres):.1f} % are unresolved.",
        [(f"Recurrent tumors are less human on average ({sum(rec) / 3:.0f} % against {sum(pri) / 3:.0f} %; ranges overlap): a confound, tested on slide 8.", BLUE)]],
        size=16, gap=10)
    C.only("nag26")
    refs(s, C.line(DESIGN))
    notes(s, ["0:55–1:10  The model the Henry Ford group built (Nagaraja 2026): U251N in the athymic rat, MRI-guided LITT (design on the previous slide). "
              "The graft is human and the host is rat, so xengsort sorts every read by species before anything is counted; human reads plus the "
              "shared reads enter the analysis.",
              "If asked why three per arm when Nagaraja 2026 describes four and four: the libraries come from animals numbered 64 to 71; IL64B's "
              "sample holds almost no tumor (it is used only as a control) and there is no library for animal 65 (a recurrence in the lab's key).",
              "Numbers to reconcile before the talk: the J Neurosurg paper reports four per arm for RNA-seq; this analysis has three per arm plus "
              "IL64B, recorded by the lab as a primary but with almost no tumor in its sample (0.3 % human reads), used as a control. The lab's "
              "sample key lists IL-65 as the fourth recurrence; no library for it exists in the data.",
              "The two other controls are the contralateral hemispheres of tumor-bearing rats 68 and 69 (lab key: 'Contralateral'; array "
              "sheet: '69B control N2'), not separate animals. N269B, from rat 69, the most tumor-rich animal, carries human Y-linked reads, "
              "so it holds some tumor cells."])

    # ---- 5 PCA + volcano
    s = content_slide(prs, "The arms separate modestly, and few single genes move")
    a, b = FIG / "fig_pca.png", FIG / "chart_volcano.png"
    H = 3700000
    wa, _ = fit(a, 10 ** 9, H); wb, _ = fit(b, 10 ** 9, H)
    x0 = L + (R - L - wa - wb - GAP) // 2
    picture(s, a, x0, TOP, wa, H); picture(s, b, x0 + wa + GAP, TOP, wb, H)
    assert pca["margin_pc1"] > 0 and not sdist["k2_recovers_arms"]
    add_text(s, L, TOP + H + 150000, R - L, BOTTOM - (TOP + H + 150000), [
        f"The recurrent tumors lie to the right on PC1, IL66B only just; clustering does not separate the arms (backup B5). "
        f"DESeq2{C('love')}: {F['de_total']} genes at FDR < 0.05 and two-fold. So we read whole programs."], size=16, gap=0)
    refs(s, C.line(f"PCA of the 500 most variable rlog genes (DESeq2 plotPCA); PERMANOVA on all genes R² {pvn['r2']:.2f}, p = {pvn['p']:.2f}, the smallest 3 v 3 allows; "
                   "Wald test with ashr shrinkage. " + DESIGN))
    notes(s, ["Printed abstract: PC1 41.8 % / PC2 26.2 % came from the May figure script (vst); this plot uses rlog (35 % / 25 %). PERMANOVA R2 0.278 and 35 genes (23 up, 12 down) are unchanged; the abstract omits p = 0.10. Its hubs (BGN, COL1A1, IGFBP3, CALB1): only 11 of the 35 genes connect, and CALB1 has one interaction (SLIDES/ABSTRACT_VS_TALK.md).",
              "1:00–1:25  Global structure in one breath. PC1 carries " + f"{100 * pca['var_pc1']:.0f} % of the variance in the 500 genes; the recurrent tumors sit to the right, "
              f"the closest (IL66B) {pca['margin_pc1']:.1f} units past the nearest primary on a {pca['pc1_range']:.0f}-unit axis. PC1 also splits the NL animals from the IL "
              f"animals (R² {pca['r2_pc1_prefix']:.2f} against {pca['r2_pc1_arm']:.2f} for arm): find out what the prefix denotes before the talk. With three per arm, ten "
              "splits exist, so PERMANOVA cannot go below 0.10. Only 35 genes pass, mostly readthrough and non-coding loci, most of them zero in one arm."])

    # ---- 6 GSEA + leading edge
    s = content_slide(prs, "After LITT, the regrown tumor turns its ribosomal-protein genes down")
    a, b = FIG / "fig_gsea_ti.png", FIG / "fig_le_heatmap.png"
    H = 3900000
    wa, _ = fit(a, 10 ** 9, H); wb, _ = fit(b, 10 ** 9, H)
    x0 = L + (R - L - wa - wb - GAP) // 2
    picture(s, a, x0, TOP, wa, H); picture(s, b, x0 + wa + GAP, TOP, wb, H)
    q6 = sorted(g["six_q"].values())
    assert g["n_down_q25"] == 6 and g["n_up_q25"] == 0 and q6[1] >= 0.05
    add_text(s, L, TOP + H + 120000, R - L, BOTTOM - (TOP + H + 120000), [
        f"Of {F['gsea_n_sets']:,} gene sets{C('subramanian')}, six reach q < 0.25, all down in recurrence and all built on the same {g['leading_edges']['shared_by_all_six']} "
        f"ribosomal-protein genes; only this set is below 0.05 (q = {q6[0]:.3f}). No up-regulated set gets below q = {g['best_up_q']:.2f}.",
        [(f"The heatmap shows the catch: one primary, IL68B, is the highest of the six on all {le['n_genes']} leading-edge genes.", BLUE)]],
        size=15, gap=4)
    refs(s, C.line("Broad GSEA, gene-set permutation; metric Diff_of_Classes on normalised counts (lower panel clipped at the 99.5th percentile); "
                   "leading edge = core-enrichment genes, rlog row z-scores; n = 3 per arm."))
    assert g["ti_all_ribosomal_protein"] and all(v[0] >= 83 for v in g["rp_in_other_leading_sets"].values())
    ls_ = g["log_scale"]
    notes(s, ["Printed abstract: elongation NES -3.17, FDR 1.6e-26 came from a first-pass gene-permutation test (clusterProfiler/fgsea, 25 May report). The pipeline's Broad GSEA (gene-set permutation) gives the same direction, NES near -2, and one set at q = 0.022. Its 'integrated stress' sets are the same ribosomal-protein genes (85 of 99 in GCN2). If asked: 'the abstract's FDR came from a first-pass test that treats genes as independent; the permutation test gives the same direction'.",
              "1:25–2:05  The core result. Say 'ribosomal-protein genes', not 'translation initiation': the KEGG MEDICUS set of that name is 80 of 80 "
              "cytosolic ribosomal-protein genes (RPL, RPS, FAU, UBA52) with no initiation factor; the other five sets are the same block counted again "
              f"(their leading edges share {g['leading_edges']['shared_by_all_six']} genes, union {g['leading_edges']['union']});. "
              "Whether translation itself falls is the puromycin / polysome experiment on slide 13.",
              "Left, the standard GSEA plot: genes ranked from higher in recurrent (left) to higher in primary (right); the 80 genes are the ticks, "
              f"NES {m(g['NES'])}. The ranking is a difference of means on counts, which weights abundant genes; on a log-scale ranking all 80 are "
              f"still in the lower half and {100 * ls_['share_in_last_6_6pct']:.0f} % in the last 6.6 %, so the direction does not depend on the metric.",
              "Right, the leading edge per tumor. Say it plainly: IL68B is high on every gene, the other two primaries sit near the middle, and the "
              "three recurrences sit below them on average. The direction is shared; the size leans on one tumor. Next slide asks whether it survives without it."])

    # ---- 7 leave one out
    s = content_slide(prs, "Leave any tumor out: the direction holds, the FDR does not")
    nq = int((filt.fdr_q < 0.05).sum())
    fig_and_text(s, FIG / "fig_loo_heatmap.png", [
        f"The published GSEA re-run {len(filt)} times: each tumor left out, five permutation seeds.",
        f"Every run: NES {m(float(filt.NES.max()))} to {m(float(filt.NES.min()))}, nominal p < 0.002.",
        f"FDR q < 0.05 in {nq} of {len(filt)}; with all six it moves with the seed alone ({lo['grid']['All six'][0]:.3f}–{max(lo['grid']['All six']):.2f}).",
        [("Dropping IL68B keeps q < 0.05 at every seed: it carries the size of the fall, not its direction.", BLUE)]],
        img_w=6500000, size=16, gap=12, text_top=60000)
    refs(s, C.line("ANALYSIS/gsea_leave_one_out: the published command, gene filter re-applied after each drop; seed 1234 reproduces the published "
                   "report exactly; blue cells q < 0.05; NES range over the five seeds."))
    notes(s, ["Printed abstract: its FDR is not a Broad GSEA q; here q runs 0.022-0.31 with all six tumors, depending on the seed alone.",
              "2:05–2:35  Robustness, stated honestly. The direction (NES) never moves. FDR is fragile: with gene-set permutation and six tumors, q "
              "depends on the random seed. With the filter re-applied, q never clears 0.05 without IL67B, IL69B or IL66B, and always does without "
              f"IL68B or NL70B. Six further runs skip the re-filtering (seed 1234): q < 0.05 in {int((nof.fdr_q < 0.05).sum())} of 6, including without IL66B "
              f"(q = {q66_nof:.3f}). Without IL68B the set is no longer the top-ranked one; a stromal set leads."])

    # ---- 8 purity
    s = content_slide(prs, "Recurrence after LITT, or less tumor in the sample?")
    fig_top(s, FIG / "fig_graft_corr.png", [
        f"Both scores rise with the human share (r = {gc[TIK]['r']:.2f} and {gc['Neftel_AC']['r']:.2f}; P = {gc[TIK]['p']:.2f}, {gc['Neftel_AC']['p']:.2f}), and arm and "
        f"share are correlated (r = {m(gr['r_group_graft'])}): six tumors cannot fully separate the two.",
        [(f"IL67B and NL70B have the same share ({pair['IL67B']['graft']:.1f} and {pair['NL70B']['graft']:.1f} %) yet differ on the astrocyte-like score "
          f"({m(pair['IL67B']['AC'], True)} against {m(pair['NL70B']['AC'], True)}). Rat reads could shift the ribosomal-protein genes by about 0.01 log2; they fall "
          f"{abs(ini['observed_mean_lfc']):.2f}.", BLUE)]],
        reserve=1050000, size=15, gap=4)
    refs(s, C.line("GSVA on DESeq2 variance-stabilised counts; OLS fit with 95 % confidence band; human share from xengsort; contamination bound: every "
                   "shared read treated as rat. " + DESIGN))
    notes(s, ["2:35–3:10  The confound, head on. Adjusting for human share leaves the ribosomal-protein score's arm difference at p = "
              f"{gr['TI_group_p_adj']:.2f} (from {gr['TI_group_p']:.2f}); the gene-level estimate does not shrink ({m(gr['TI_lfc_group'])} to {m(gr['TI_lfc_group_adj'])} log2), "
              "though at n = 6 the adjusted estimate is unstable (it leans on IL69B; without IL68B it reverses).",
              f"The matched pair: IL67B (primary) and NL70B (recurrent) are both {pair['IL67B']['graft']:.0f} % human. On the astrocyte-like score they "
              f"differ by {abs(pair['NL70B']['AC'] - pair['IL67B']['AC']):.2f}, far beyond the spread within either arm. On the ribosomal-protein score the gap "
              f"({m(pair['NL70B']['TI'] - pair['IL67B']['TI'])}) is within the spread among the primaries, so the pair says nothing there. One pair, not a test.",
              f"Contamination: counting every shared read as rat, with the two tumor-free controls, gives {ini['mixture_bound_mean_lfc']:+.3f}; with N269B "
              f"(the contralateral hemisphere of rat 69, which holds tumor cells) {ini3['mixture_bound_mean_lfc']:+.3f}, so under 0.01 log2 either way. A second estimate from the controls' shared-read "
              "ratio puts the worst case near 0.05 (ANALYSIS/human_cohorts/PLAN.md), an eighth of the fall. Backup B1."])

    # ---- 9 subtypes
    s = content_slide(prs, "Of ten published subtype signatures, only the astrocyte-like score separates the arms")
    fig_and_text(s, FIG / "fig_subtype_heatmap.png", [
        f"Neftel{C('neftel')} and Garofano{C('garofano')} signatures, scored per tumor by GSVA{C('hanzelmann')}.",
        f"AC-like: {m(ac['primary'], True)} → {m(ac['recurrent'], True)}, every recurrent below every primary (t-test P = {ac['p']:.3f}, q = {ac['q']:.3f}; "
        f"exact floor 0.10); direction kept after adjusting for human share.",
        f"MES1 is lower in recurrence (P = {mes1['p']:.2f}) but tracks human share as closely (r = {mes1['r_graft']:.2f}, P = {gr['p_MES1']:.2f}).",
        f"CIBERSORT{C('newman')} with the Neftel reference fits these cells poorly (r {min(fitr):.2f}–{max(fitr):.2f}; backup B3)."],
        img_w=6600000, size=15, gap=10, text_top=40000)
    refs(s, C.line("GSVA on variance-stabilised counts; two-sample t per signature, Benjamini–Hochberg across the ten. " + DESIGN))
    notes(s, ["Printed abstract: 'toward Garofano mitochondrial, away from AC/NPC' came from May mean z-scores of hand-picked 5-6-gene lists (mitochondrial = CS, ACO2, IDH2, IDH3A, OGDH, SDHA), not GSVA. With the published signatures the mitochondrial score falls (not significant); AC-like agrees; NPC does not move. If asked: 'the mitochondrial shift did not survive the published signatures; the astrocyte-like fall did'.",
              "3:10–3:40  One state transition: the astrocyte-like program score falls, and every recurrent tumor is below every primary. Call it the "
              "strongest group difference and exploratory: 3 v 3 cannot give an exact permutation p below 0.10.",
              f"Adjusted for human share, P = {ac['p_adj_graft']:.3f} on these counts; across other normalised inputs 0.024 to 0.12, and the adjustment leans on "
              "IL69B (graft_relation/SUMMARY.md). The manuscript's run scored TPM (p = 0.007, q = 0.075); 7SK/7SL/Y RNA and rRNA take 40 to 65 % of TPM and "
              "the share varies by library, which variance-stabilised counts remove.",
              f"CIBERSORT puts AC slightly higher in recurrence ({rb['deconv']['AC_primary']:.2f} to {rb['deconv']['AC_recurrent']:.2f}), the opposite direction; at this fit neither is a composition.",
              f"Whole pipeline without IL68B or IL66B (TPM run, from {ho['all_six']['AC_p']:.3f} with all six): AC p = {ho['IL68B']['AC_p']:.3f} and {ho['IL66B']['AC_p']:.3f} (backup B2)."])

    # ---- 10 the drug (exact-name ChEMBL matching, ANALYSIS/drug_rematch; supersedes the first-search-hit ranking of S12)
    RM = ROOT / "ANALYSIS" / "drug_rematch"
    rmj = json.load(open(RM / "rematch.json", encoding="utf-8"))
    ex, oldm = rmj["published"]["exact matching"], rmj["published"]["old matching"]
    cw, dw = ex["watch"]["ciclopirox"], ex["watch"]["deferoxamine"]
    c68, c66 = rmj["IL68B"]["exact matching"]["watch"]["ciclopirox"], rmj["IL66B"]["exact matching"]["watch"]["ciclopirox"]
    assert cw["rank_bbb"] == 1 and c68["rank_bbb"] == 1 and c66["rank_bbb"] == 3 and dw["rank_nobbb"] == 1
    rk = pd.read_csv(RM / "results" / "published_drug_ranking_final.csv").sort_values("rank_bbb").reset_index(drop=True)
    assert rk.Drug.iloc[0] == "ciclopirox" and len(rk) == ex["n_clinical"]
    nxt2 = rk[rk.both_agree.astype(bool) & (rk.Drug != "ciclopirox")].iloc[0]
    sep20 = set(esm.head(20).Drug.str.lower()); new20 = list(rk.head(20).Drug)
    entered20 = [d for d in new20 if d.lower() not in sep20]
    s = content_slide(prs, "A candidate drug against the recurrent state: ciclopirox")
    rows = [["Rank", "Compound", "NES", "BBB probability", "Both BBB models"]]
    for _, r in rk.head(5).iterrows():
        rows.append([str(int(r.rank_bbb)), str(r.Drug), m(float(r.NES)), f"{float(r.BBB_Martins):.2f}", "yes" if bool(r.both_agree) else "no"])
    widths = [680000, 1900000, 950000, 1400000, 1450000]
    table(s, L, TOP + 50000, widths, rows, size=17, bold_rows=(1,), row_h=470000)
    tx = L + sum(widths) + GAP
    add_text(s, tx, TOP + 50000, R - tx, BOTTOM - TOP, [
        f"DSigDB drug signatures{C('yoo')} against the recurrence ranking, used as a ranking: the 100 most opposing name {fun['n_compounds']} compounds; "
        f"{ex['n_clinical']} have a clinical-phase ChEMBL record; {ex['n_both_agree']} pass both barrier models{C('swanson', 'daina')}.",
        f"Ciclopirox ranks first with the barrier weight (second without it); first without IL68B, third without IL66B. Another group's reversal screen also nominated it.{C('sun')}",
        f"Approved as a topical antifungal; an oral form completed phase 1 in hematologic cancer.{C('minden')}"],
        size=15, gap=12)
    add_text(s, L, TOP + 50000 + 470000 * len(rows) + 150000, sum(widths), 1200000, [
        f"Without the barrier weight, deferoxamine, another approved iron chelator, leads (NES {m(dw['NES'])}); the barrier models disagree on it "
        f"(ADMET-AI {dw['BBB_Martins']:.2f}, borderline; BOILED-Egg, no). Ranks 2–5 share one ribosomal-protein gene set: one signal, not four candidates.",
        [("Computational predictions: nothing here has been dosed in this model.", BLUE)]], size=15, gap=8)
    refs(s, C.line(f"Rank among the {ex['n_clinical']} compounds with a clinical-phase ChEMBL record (exact-name matching, audited) by |NES|^1.5 × ADMET-AI BBB "
                   "probability; BBB models ADMET-AI and BOILED-Egg; supersedes Online Resource 1, S12."))
    notes(s, ["3:40–4:15  The drug. Each compound's DSigDB gene set is tested for sitting among the genes that fall in recurrence, then weighted by "
              f"predicted barrier penetration. Ciclopirox scores {cw['score_bbb']:.2f} against {float(nxt2.score_bbb):.2f} for {nxt2.Drug}, the next compound both "
              f"barrier models pass. Without IL66B it is third, {c66['rank_nobbb']}th without the barrier weight: say so if asked.",
              "What changed since the manuscript's S12: compound names are now matched to ChEMBL by exact name or synonym and mapped to the parent "
              "molecule, with every changed or unresolved name audited (ANALYSIS/drug_rematch). The old pipeline took the first free-text search "
              f"hit, which dropped approved drugs (deferoxamine, vandetanib, mefloquine, bromocriptine) and used salt records whose counterions made "
              f"escitalopram, paroxetine, metoprolol, trimipramine and propantheline look barrier-impermeable. Clinical compounds {oldm['n_clinical']} -> "
              f"{ex['n_clinical']}; both barrier models {oldm['n_both_agree']} -> {ex['n_both_agree']} (same software, same day). Ciclopirox's rank with the "
              "barrier weight did not move in any of the three runs. The manuscript's S12 and prior-art table need the same update.",
              f"Deferoxamine: first without the barrier weight (score {dw['score_nobbb']:.2f} against ciclopirox's {cw['score_nobbb']:.2f}), {dw['rank_bbb']}th with it; "
              f"ADMET-AI {dw['BBB_Martins']:.2f} (just over the 0.5 cut) and BOILED-Egg 'out', so the models disagree. Two approved iron chelators at the top of the unweighted list is itself worth saying.",
              f"New entries to the top twenty since the prior-art audit (backup B4): {', '.join(entered20) if entered20 else 'none'}.",
              "Caveats if asked: scored per tumor by GSVA, the ciclopirox gene set does not separate the arms (P = 0.89), so the rank rests on the "
              "whole-list ordering; its DSigDB sets hold at most one ribosomal-protein gene, so the rank is not the ribosomal block again. Brain "
              "pharmacokinetics of ciclopirox have not been measured; the barrier call is a model prediction.",
              "What carries ciclopirox's rank: its induced genes are HIF/hypoxia targets (21 of 36 leading-edge genes are hallmark hypoxia; "
              "without them NES -1.57, p 0.005). Ranks 2-5 (pentetrazol, nilutamide, primidone, pyrantel) share one leading edge of "
              "ribosomal-protein genes plus PGK1 and ENO1 (pairwise Jaccard 0.55-0.86): the ribosomal block counted four times "
              "(ANALYSIS/gsc_drugs/REVIEW_2026-09-26.md).",
              "Stem-cell programs: restricted to six published glioma stem-cell gene sets, ciclopirox's induced genes fall in recurrence only "
              "within the radial-glia-like program (Bhaduri 2020; first of 2,212 drug sets, NES -2.08), and those genes are hypoxia targets, "
              "not stem-cell genes. No drug opposed recurrence inside the stem-specific sets. Consistent with this slide, not independent of it.",
              "The printed abstract: its 3.40 for ciclopirox is this enrichment with a rule-of-thumb barrier score capped at 1.0; with ADMET-AI "
              "it is 3.27 and still first. DMOG and LY-294002, named in the abstract, have no clinical-phase record and drop out. The drug NES "
              "come from gene-permutation GSEA, so they are quoted as a ranking, not with FDRs (SLIDES/ABSTRACT_VS_TALK.md)."])

    # ---- 11 why ciclopirox
    s = content_slide(prs, "Why ciclopirox: U-251 is past the threshold on the axis it inhibits")
    dd_ = dg["DOHH"]
    fig_top(s, FIG / "fig_depmap.png", [
        f"Ciclopirox chelates iron and inhibits DOHH, which activates eIF5A for translation elongation. In DepMap{C('depmap')}, U-251 MG is past the "
        f"dependency threshold for DOHH ({m(dd_['u251'])}), more dependent than {dd_['cns_more_dependent_than_pct']:.0f} % of {dd_['cns_n']} CNS lines (median {m(dd_['cns_median'])}).",
        [(f"In mice, LITT opens the barrier around the ablation for one to three weeks{C('cleary')}: a possible delivery window.", BLUE)]],
        reserve=1000000, size=15, gap=4)
    refs(s, C.line("Violins: all screened lines (n per gene 1,043–1,178); dots, CNS/brain lines; bar, interquartile range; white dot, median; dashed, −0.5; "
                   "dotted, −1 (median common-essential gene)."))
    notes(s, ["4:15–4:45  The mechanism argument, kept modest. DHPS and EIF5A are needed by most lines "
              f"({100 * dg['DHPS']['share_below_minus05']:.0f} % and {100 * dg['EIF5A']['share_below_minus05']:.0f} % below −0.5), so they are not selective; DOHH is the step most lines can do "
              f"without ({100 * dg['DOHH']['share_below_minus05']:.0f} % below −0.5), and U-251 needs it more than the typical glioma line, by a modest margin. U-251 is broadly "
              "dependent, so a percentile like this is ordinary for it. Ribonucleotide reductase "
              f"(RRM1 {m(dg['RRM1']['u251'])}, RRM2 {m(dg['RRM2']['u251'])}) is needed by nearly every line. HIF1A is not a dependency ({m(dg['HIF1A']['u251'], True)}). A public screen of "
              "the line in culture, not our data.",
              "Barrier window, from Cleary 2026 in SB28 mouse glioma: tight junctions open for 7 days up to 100 µm from the ablation; transcytosis "
              "peaks at day 14 and is back to baseline by day 21. No human post-LITT time course exists."])

    # ---- 12 what this means
    s = content_slide(prs, "What this means")
    add_text(s, L, TOP + 100000, R - L, 3000000, [
        "After LITT, the regrown tumor carries a program of its own: lower ribosomal-protein genes and a lower astrocyte-like score.",
        "Its direction holds with any tumor left out; its size leans on one tumor, its FDR on the seed, and tumor content is not yet separated from it.",
        "It ranks first a topically approved iron chelator predicted to cross the barrier (third without one recurrent tumor), on an axis this line depends on."],
        size=19, gap=18)
    navy_bar(s, L, TOP + 3350000, R - L, BOTTOM - (TOP + 3350000),
             "A hypothesis with a mechanism and a window, from six tumors: worth the experiments on the next slide.", size=19)
    refs(s, C.line("Evidence status: one cell line (U251N) in one rat xenograft model, n = 3 per arm; drug candidates are computational predictions; GEO GSE338105."))
    notes(s, ["Say, if time allows: 'Three things differ from the printed abstract: lower ribosomal-protein genes rather than translational shutdown; the mitochondrial shift did not survive the published signatures; and of its three drugs only ciclopirox is clinically available.'",
              "4:45–5:15  Three sentences, then the bar. Do not add numbers here."])

    # ---- 13 limitations -> experiments
    s = content_slide(prs, "Limitations, and the experiments that answer them")
    half = (R - L - GAP) // 2
    left = ["Limitations",
            "Three tumors per arm; the size of the ribosomal-protein fall leans on one.",
            "Lower ribosomal-protein mRNA is not yet lower translation.",
            "Arm and human share are partly collinear; no sham or non-thermal regrowth arm.",
            "One cell line in an athymic rat: no T cells, no patient tissue."]
    right = ["Next",
             "Burden-matched primaries, sham-fiber and non-thermal cytoreduction arms, more animals.",
             "Puromycin incorporation and polysome profiling: is translation suppressed?",
             f"Paired IDH-wild-type primary and recurrent glioblastoma (GLASS{C('varn')}): does it recur in patients?",
             "Ciclopirox dosed into the post-LITT barrier window, survival as the endpoint."]
    for x, block in ((L, left), (L + half + GAP, right)):
        head_box(s, x, TOP, half, 460000, block[0], [], size=22)
        add_text(s, x, TOP + 620000, half, BOTTOM - TOP - 620000, block[1:], size=17, gap=12)
    refs(s, C.line(DESIGN))
    notes(s, ["Printed abstract: its 'translational shutdown' is what the puromycin / polysome experiment tests; it has not been shown.",
              "5:15–5:45  Say the limitations as the reason for each experiment. GLASS recurrences follow standard therapy, not LITT, so it tests "
              "whether the signature accompanies recurrence in general. Bulk RNA: the margin is inferred, not isolated. Methylation on the same six "
              "tumors (EPIC, 866,238 probes) found no probe past FDR; with purity and arm confounded it cannot separate the two. If asked: host (rat) "
              "reads are being analysed separately for the microenvironment; not in this talk."])

    # ---- 14 acknowledgments
    s = content_slide(prs, "Acknowledgments and data availability")
    add_text(s, L, TOP, 6300000, BOTTOM - TOP, [
        "Co-authors Tavarekere N. Nagaraja, PhD, Indrani Datta, DHI, and Ian Y. Lee, MD, Hermelin Brain Tumor Center, Henry Ford Health: the model, the ablations, the tissue and the study design.",
        "Supported by a Henry Ford Health Physician Scientist Award A20050 (I.Y.L.).",
        "Henry Ford Health IACUC protocol #1509; ARRIVE guidelines.",
        "Sequencing data: GEO GSE338105."], size=17, gap=12)
    navy_bar(s, R - 3600000, TOP + 300000, 3600000, 1300000, ["Questions", "go2432@wayne.edu"], size=20)
    notes(s, ["5:45–6:00  Thanks. Leave the GEO accession on screen during questions. Indrani Datta's credential is DHI (confirmed by Greg, 2026-09-26). The JNS 2026 author list prints PhD; do not copy it."])
    assert C.line() == ""

    # ---- backups
    s = content_slide(prs, "Backup B1. Is the fall rat contamination?")
    fig_and_text(s, FIG / "chart_contamination_magnitude.png", [
        f"Counting every shared read as rat would move the ribosomal-protein sets by under 0.01 log2 ({ini['mixture_bound_mean_lfc']:+.3f}; the two tumor-free controls).",
        f"They fall {m(ini['observed_mean_lfc'])}; background genes matched on abundance and control/tumor ratio move {ini['matched_background_mean_lfc']:+.2f}.",
        f"Without IL68B the fall is {m(mag['leave_one_out']['sets']['initiation']['without_IL68B'])}."],
        img_w=6600000, size=16, gap=12, text_top=60000)
    refs(s, C.line(f"Controls: IL64B (no tumor in its sample) and N168B (opposite hemisphere of rat 68); N269B holds tumor cells and is left out. RUVSeq{C('risso')} "
                   "factors (k = 2, from the three control libraries; rerun on the two tumor-free ones pending) also separate the arms."))
    notes(s, ["If the printed abstract is the docx version: its claim that contamination 'would bias ribosomal genes upward' does not match this bound. The size argument holds; the direction depends on N269B, the contralateral hemisphere of rat 69 that holds tumor cells.",
              f"Sensitivity, N269B among the controls: bound {ini3['mixture_bound_mean_lfc']:+.3f}, matched background {ini3['matched_background_mean_lfc']:+.2f}, "
              f"excess {ini3['excess_over_matched']:+.2f} (tumor-free controls: {ini['excess_over_matched']:+.2f}); contamination_magnitude_ctrl3.json.",
              "Backup. Mean log2 fold change (recurrent vs primary) of each set under each analysis; the bound is a model ceiling, not an observed change."])

    s = content_slide(prs, "Backup B2. The whole pipeline re-run without one tumor")
    rows = [["", "DE genes", "Ribosomal-protein set q", "Ciclopirox rank (weighted / not)", "AC-like change", "AC-like p"]]
    for lab, k in (("All six", "all_six"), ("Without IL68B", "IL68B"), ("Without IL66B", "IL66B")):
        r = ho[k]
        cr_ = {"all_six": cw, "IL68B": c68, "IL66B": c66}[k]            # ciclopirox ranks from the exact-name rerun
        rows.append([lab, str(r["de_n"]), f"{r['q_TI']:.3f}", f"{cr_['rank_bbb']} / {cr_['rank_nobbb']}", m(r["AC_change"]), f"{r['AC_p']:.3f}"])
    table(s, L, TOP + 100000, [1900000, 1200000, 2000000, 2350000, 1400000, 1200000], rows, size=16, row_h=520000)
    add_text(s, L, TOP + 100000 + 520000 * 4 + 250000, R - L, 1200000, [
        "Every step re-run from counts: differential expression, GSEA (seed 1234), subtype scoring, drug ranking. IL68B was chosen because it carries the "
        "ribosomal-protein fall; IL66B because it is the library-QC outlier. Both choices were made after seeing the data."], size=15, gap=0)
    refs(s, C.line("ANALYSIS/holdout_IL68B and holdout_IL66B, comparison.json; ciclopirox ranks from the exact-name ChEMBL rerun (ANALYSIS/drug_rematch). Subtype p here is the published TPM-based GSVA run."))
    notes(s, ["Backup. Post hoc holdouts: say so."])

    ctl_r = [cs[k]["Correlation"] for k in ("IL64B", "N168B", "N269B")]
    pv_ = [cs[k]["P-value"] for k in cs]
    s = content_slide(prs, "Backup B3. CIBERSORT deconvolution of the Neftel states")
    fig_and_text(s, FIG / "fig_cibersort.png", [
        f"CIBERSORT v1.04{C('newman')}, signature built from Neftel{C('neftel')} single-cell profiles. P < 0.05 in every library, but fit r "
        f"{min(fitr):.2f}–{max(fitr):.2f} in the tumors and the culture: significant against random mixtures, poor as a description.",
        f"The rat-brain controls, almost no human reads, get fractions too (fit r {min(ctl_r):.2f}–{max(ctl_r):.2f}).",
        f"AC-like {rb['deconv']['AC_primary']:.2f} → {rb['deconv']['AC_recurrent']:.2f}, opposite to the GSVA score; MES-like {rb['deconv']['MES_primary']:.2f} → {rb['deconv']['MES_recurrent']:.2f}.",
        f"Across the five runs with tumor-state columns (three Neftel builds, two GLASS), no state differs at P < 0.05; GLASS proliferating stem-like is "
        f"higher in all three recurrences (P {min(glass_p):.2f}–{max(glass_p):.2f})."],
        img_w=6400000, size=14, gap=8, text_top=20000)
    refs(s, C.line("Relative mode, 1,000 permutations, no quantile normalisation; Ivy GAP, glioma atlas, LM22 and BayesPrism references have no tumor-state columns; ANALYSIS/cibersort."))
    assert max(pv_) < 0.05
    notes(s, ["Backup. CIBERSORT estimates the share of each state in a bulk sample; it does not sort reads."])

    s = content_slide(prs, "Backup B4. What is already known about the top twenty")
    img = FIG / "table_prior_art.png"
    w, h = fit(img, R - L, BOTTOM - TOP)
    picture(s, img, L + (R - L - w) // 2, TOP, w, h)
    refs(s, C.line(f"Tier A: in vivo or clinical glioma evidence (ciclopirox{C('su')}); B: in vitro or contested; C: none, failed or not a therapy. PubMed per compound."))
    notes(s, ["Backup. Prior-art audit of the top twenty of the 54 in the September ranking (S15); the exact-name rerun changes the top twenty "
              "(see slide 10 notes), so new entrants have no audit yet. Known inconsistency to fix in the manuscript: "
              "paroxetine and amiodarone have in vivo evidence but sit in B; diazepam's evidence is contested but sits in C."])

    s = content_slide(prs, "Backup B5. Sample-to-sample distances")
    fig_and_text(s, FIG / "fig_sample_dist.png", [
        "Euclidean distance on rlog expression, all genes, hierarchical clustering (the DESeq2 quality-control view).",
        "The two NL recurrences pair; IL68B stands apart from everything; IL66B joins the other two primaries.",
        "Unsupervised clustering does not recover the arms: whatever separates them is a few coordinated programs, not the whole transcriptome."],
        img_w=6300000, size=16, gap=12, text_top=60000)
    refs(s, C.line(DESIGN))
    notes(s, ["Backup. Show if asked whether the arms separate globally."])

    s = content_slide(prs, "Backup B6. The strongest gene sets in each direction")
    fig_top(s, FIG / "fig_gsea_dotplot.png", [
        f"Down in recurrence: the ribosomal-protein block (six sets at q < 0.25). Up: mitotic and cell-cycle sets lead, none below q = {g['best_up_q']:.2f}."],
        reserve=650000, size=15, gap=4)
    refs(s, C.line(f"Broad GSEA{C('subramanian')}, {F['gsea_n_sets']:,} sets; point size, set size; colour, FDR q."))
    notes(s, ["Backup. The standard GSEA summary: eight strongest sets each way."])
    C.line()

    hs = J("host.json")                                     # SLIDES/11_host_figure.py, from ANALYSIS/host/results_de + verify/
    assert hs["sig2"] == hs["up2"] + hs["down2"] and hs["down_all_three"] == hs["down"]
    assert hs["human_MMP13_tumour_max"] < hs["rat_Mmp13_primary_min"]          # raw reads on both sides
    weak = hs["il66b_weak_immune"]
    rl, bg, lh = hs["relabel"], hs["broad_gsea"], hs["loo_host"]
    robust = hs["loo_all_six"]
    rob_dn = [g for g, d in robust if d == "down"]; rob_up = [g for g, d in robust if d == "up"]
    kept = [v for k, v in hs["loo_kept_of_38"].items() if k != "IL66B"]
    n_sets = bg["up"]["n"] + bg["down"]["n"]
    min_q = min(bg["up"]["min_fdr"], bg["down"]["min_fdr"])
    big = rl["largest_target"].replace("IL69B", "69").replace("NL70B", "70").replace("NL71B", "71").replace(",", ", ")
    s = content_slide(prs, "Backup B7. The rat host tissue in recurrence (preliminary)")
    fig_and_text(s, FIG / "fig_host_volcano.png", [
        f"Rat reads only, 3 v 3, the tumor's rule (FDR < 0.05 and two-fold): {hs['sig2']} genes, {hs['down2']} down and {hs['up2']} up.",
        f"Holds when any tumor is left out: {', '.join(rob_dn)} down (Mmp13 about {hs['mmp13_fold']:.0f}-fold), {', '.join(rob_up)} up. "
        f"Not human reads: human MMP13 at most {hs['human_MMP13_tumour_max']} reads per tumor.",
        f"The rest is unstable: leaving out one tumor keeps {min(kept)}–{max(kept)} of the {hs['sig']}; no gene set passes FDR 0.25 "
        f"({n_sets:,} sets, lowest {min_q:.2f}).",
        f"The true grouping ranks {rl['true_rank']} of the {rl['n_splits']} possible 3 v 3 splits; the largest (rats {big} against the rest, "
        f"{rl['largest_n']:,} genes) follows immune infiltration, not treatment. The interferon/lymphocyte rise comes from rats 70 and 71."],
        img_w=6400000, size=14, gap=8, text_top=20000)
    refs(s, C.line(f"Preliminary: different rats per arm; LITT, time and the ablation wound cannot be separated. Broad GSEA from the pipeline; "
                   f"leave-one-out and relabelling from a rerun that reproduces the pipeline's {hs['sig']} genes. DESeq2{C('love')}."))
    notes(s, ["Backup, preliminary; show only if asked about the host. What holds: " + ", ".join(rob_dn) + " lower in every recurrence and "
              "in every leave-one-out fit (Mmp13 is also not human-read leakage: human MMP13 at most "
              f"{hs['human_MMP13_tumour_max']} reads per tumor against at least {hs['rat_Mmp13_primary_min']} rat reads per primary). "
              f"What does not: the {hs['sig']}-gene list (FDR < 0.05) shrinks to {min(kept)}–{max(kept)} when a primary or rat 70 or 71 is left out; "
              f"leaving IL66B out raises it to {lh['IL66B']} (the tumor side does the same, {hs['loo_human']['none']} to {hs['loo_human']['IL66B']}), "
              "so IL66B is the recurrence most like a primary in both species. "
              f"The true grouping gives {rl['true_n']} genes; the other nine splits give a median of {rl['median_others']}; the largest, "
              f"rats {big} against the rest ({rl['largest_n']:,} genes), is an immune-infiltration (Ptprc, Cd68, Aif1) against "
              "white-matter/neuron (Mbp, Plp1, Pcp4) axis that cuts across treatment. " + ", ".join(weak) + " barely move in IL66B.",
              f"Broad GSEA (the pipeline's run, gene-set permutation): {n_sets:,} sets, none at FDR < 0.25. By nominal p, up: antigen "
              "presentation, allograft rejection, interferon gamma; down: KEGG ribosome and translation initiation (NES about -2.0, "
              "nominal p about 0.04). The rat host's ribosomal-protein genes point the same way as the tumor's; their share tracks the "
              "human read share (Spearman 0.83 over six tumors), which fits read leakage but does not prove it.",
              "Host cell-type deconvolution (CIBERSORT; LM22, Zhang 2014, Bowman 2016 myeloid, and a combined signature): no cell type "
              "differs between primaries and recurrences. Tumor-associated bone-marrow macrophages are high in rats 70 and 71 (0.43, 0.50) "
              "and not in 66 (0.07): the same rats-70-71 pattern, found after looking, exploratory only.",
              "The 'NL' in the file names of rats 70 and 71 is not a group: the lab's sample key calls them IL-70 and IL-71 and its "
              "methylation-array sheet IL70B and IL71B; all ten RNA libraries ran on one lane. Why 70 and 71 differ from 66 (for example "
              "time after LITT) is not recorded. The host is the athymic RNU/RNU rat: do not describe a T-cell response. "
              "Rat reads are the xengsort host bin; nf-core/rnaseq 3.22.2 on mRatBN7.2 (Ensembl 110); nf-core/differentialabundance DESeq2. "
              "Tumor vs control brain (contralateral hemispheres of rats 68 and 69, and IL64B): "
              f"{hs['tumour_vs_control_sig']:,} genes differ."])
    C.line()


    for sl in prs.slides:                                   # a number never parts from its '%'
        for sh in sl.shapes:
            if sh.has_text_frame:
                for p in sh.text_frame.paragraphs:
                    for r in p.runs:
                        if " %" in r.text:
                            r.text = r.text.replace(" %", " %")
    prs.save(str(out))
    print("wrote", out, "slides:", len(prs.slides))
    for k, n in sorted(C.num.items(), key=lambda t: t[1]):
        print(f"  {n:2d}. {d3.REFS[k]}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "CNS2026_Schwing_Abstract418_CNStemplate_v4.pptx"))
    a = ap.parse_args()
    if not TEMPLATE.exists():
        raise SystemExit(f"CNS speaker template missing: {TEMPLATE}")
    build(Path(a.out))
