# -*- coding: utf-8 -*-
"""The CNS 2026 Abstract 418 talk (Characterizing the Transcriptomic Recurrence Signature of Glioblastoma Following
Laser Interstitial Thermal Therapy; Sunrise Science Session, Tumor 2, Monday 2 November 2026, 7:00 AM) on the
official 2026 CNS speaker template, in the style of the SchwingNet deck: one large figure per slide, the result and
its explanation in three to five lines of 15-17 pt, a citation line with the study design on every data slide, no
logos, no company names (ACCME guidelines for CNS speakers; the disclosures are read aloud from the title slide's
notes, CNS displays the list).

Version 2 (Greg 2026-09-26: "a six-minute talk, I need more results and better explanations"): every data slide now
states what the figure shows, the number, and what it means; three results slides are added from material that was
in the notes or the repository but never on a slide (what else moves: hypoxia, iron and stress down, cell cycle up;
the drug filter as a funnel; the DepMap dependency check behind ciclopirox); the overview carries a drawn pipeline
of the study instead of the cartoon.

Content is the 8 September 2026 deck (SLIDES/CNS2026_Schwing_Abstract418.pptx, built by 03_build_deck.py and sent to
the co-authors): its text, tables and timed speaker notes are read from deck_2026-09-08_text.json; every number on
a slide comes from figures/facts.json or the depmap tables. Figures by 04_make_cns_figures.py (figures_cns/).

    python SLIDES/04_make_cns_figures.py
    python SLIDES/04_build_cns_template_deck.py [--out PATH]
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from pptx import Presentation
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.util import Emu, Pt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CNS = Path.home() / "OneDrive" / "Desktop" / "CTSpinoPelvic1K-1" / "paper" / "cns2026"
sys.path.insert(0, str(CNS))
from make_cns_deck_v2 import fit  # noqa: E402
from make_cns_deck_v4 import (BLUE, BOTTOM, CAP, GAP, GREY, L, NAVY, R, REF, REFS_TOP, RED, TEMPLATE, TOP, WHITE,  # noqa: E402
                              add_text, arrow, content_slide, navy_bar, picture, refs)
from make_cns_deck_v6 import notes  # noqa: E402

FIG = HERE / "figures_cns"
DECK = {r["n"]: r for r in json.load(open(HERE / "deck_2026-09-08_text.json", encoding="utf-8"))}
F = json.load(open(HERE / "figures" / "facts.json", encoding="utf-8"))
TH = F["themes"]
TITLE = "Characterizing the Transcriptomic Recurrence Signature of Glioblastoma Following Laser Interstitial Thermal Therapy"
AUTHORS = "Gregory J. Schwing, MD, PhD¹  ·  Tavarekere N. Nagaraja, PhD²  ·  Indrani Datta, DHI²  ·  Ian Y. Lee, MD²"
AFFIL = ("¹ Department of General Surgery, Detroit Medical Center and Wayne State University, Detroit, MI   "
         "² Department of Neurosurgery, Hermelin Brain Tumor Center, Henry Ford Health, Detroit, MI")
SESSION = "Abstract 418  ·  Sunrise Science Session, Tumor 2  ·  Monday 2 November 2026, 7:00 AM"
DESIGN = ("Orthotopic U251N glioblastoma xenograft in the athymic rat, MRI-guided LITT; bulk RNA-seq of primary and recurrent "
          "tumors, n = 3 per arm; GEO GSE338105.")
REF_MINDEN = "14. Minden MD, et al. Am J Hematol 2014;89:363–368."
REF_DEPMAP = "15. DepMap 24Q4, Broad Institute (depmap.org), Chronos gene effect, cell line ACH-000232."
REF_DSIGDB = "16. Yoo M, et al. Bioinformatics 2015;31:3069–3071."


def old_text(n, name):
    for t in DECK[n]["texts"]:
        if t["name"] == name:
            return t["text"]
    raise KeyError((n, name))


def old_notes(*ns):
    return [DECK[n]["notes"].strip() for n in ns if DECK[n]["notes"].strip()]


def paras(text):
    return [p.strip() for p in text.split("\n") if p.strip()]


def nes(name):
    return TH[name]["nes"]


def m(v, signed=False):
    """A number with a typographic minus, so no line breaks after a hyphen."""
    s = f"{v:+.2f}" if signed else f"{v:.2f}"
    return s.replace("-", "−")


def pval(name):
    p = TH[name]["p"]
    return f"p = {p:.3f}" if p < 0.01 else f"p = {p:.2f}"


def head_box(slide, left, top, width, height, head, sub, size=20):
    box = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Emu(left), Emu(top), Emu(width), Emu(height))
    box.fill.solid(); box.fill.fore_color.rgb = NAVY; box.line.color.rgb = NAVY; box.shadow.inherit = False
    tf = box.text_frame; tf.word_wrap = True; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = head; r.font.size = Pt(size); r.font.bold = True; r.font.color.rgb = WHITE
    if sub:
        add_text(slide, left, top + height + 80000, width, 800000, sub, size=15, color=NAVY, gap=0, align=PP_ALIGN.CENTER)


def fig_and_text(s, img, lines, img_w, size=16, gap=10, caption=None, text_top=0):
    """A figure at the left up to img_w wide and the column height, the explanation at the right."""
    w, h = fit(img, img_w, BOTTOM - TOP - (450000 if caption else 0))
    picture(s, img, L, TOP, w, h)
    if caption:
        add_text(s, L, TOP + h + 60000, w, 420000, [caption], size=CAP, color=GREY, gap=0)
    x2 = L + w + GAP
    add_text(s, x2, TOP + text_top, R - x2, BOTTOM - TOP - text_top, lines, size=size, gap=gap)
    return w, h


def build(out: Path):
    prs = Presentation(str(TEMPLATE))
    sld = prs.slides._sldIdLst
    RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
    for e in list(sld):
        sld.remove(e)
        prs.part.drop_rel(e.get(RID))

    # ---- 1 title
    s = prs.slides.add_slide(prs.slide_layouts[0])
    ph = {p.placeholder_format.idx: p for p in s.placeholders}
    t = ph[0]
    t.left, t.top, t.width, t.height = Emu(5185558), Emu(950000), Emu(6600000), Emu(2000000)
    t.text_frame.word_wrap = True
    t.text_frame.text = TITLE
    for p in t.text_frame.paragraphs:
        for r in p.runs:
            r.font.size = Pt(24); r.font.bold = True; r.font.color.rgb = WHITE
    sub = ph[1]
    sub.left, sub.top, sub.width, sub.height = Emu(5185558), Emu(3100000), Emu(6600000), Emu(1250000)
    sub.text_frame.word_wrap = True
    sub.text_frame.text = AUTHORS
    p = sub.text_frame.add_paragraph(); p.text = AFFIL
    for i, p in enumerate(sub.text_frame.paragraphs):
        p.space_after = Pt(4)
        for r in p.runs:
            r.font.size = Pt(15 if i == 0 else 11); r.font.color.rgb = WHITE
    body = ph[10]
    body.left, body.top, body.width, body.height = Emu(5185558), Emu(4500000), Emu(6600000), Emu(1500000)
    body.text_frame.word_wrap = True
    body.text_frame.text = "CNS 2026 Annual Meeting, Washington, DC"
    p = body.text_frame.add_paragraph(); p.text = SESSION
    for p in body.text_frame.paragraphs:
        p.space_after = Pt(5)
        for r in p.runs:
            r.font.size = Pt(12); r.font.color.rgb = WHITE
    notes(s, old_notes(1) + ["Disclosures, read aloud (CNS displays the list; no disclosure slide per the ACCME speaker guidelines): "
                             + old_text(2, "TextBox 24").replace("\n", " ")] + old_notes(2))

    # ---- 2 where this is going: the five beats and the study in one strip
    s = content_slide(prs, "Where this is going")
    items = [("The problem: LITT kills the core; the cells that come back sat in the sublethal margin", 520000),
             ("The model: MRI-guided LITT in an orthotopic rat xenograft, primary and recurrent tumors, three per arm", 520000),
             ("The finding: the recurrent transcriptome switches protein synthesis off, and hypoxia signalling with it", 520000),
             ("The convergence: subtype scores, network hubs and a public dependency screen point the same way", 520000),
             ("The opening: which existing, brain-penetrant drug opposes that state, and why it might work", 520000)]
    y = TOP
    for i, (txt, dy) in enumerate(items, 1):
        box = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Emu(L), Emu(y), Emu(400000), Emu(400000))
        box.fill.solid(); box.fill.fore_color.rgb = NAVY; box.line.color.rgb = NAVY; box.shadow.inherit = False
        tf = box.text_frame; tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
        r = p.add_run(); r.text = str(i); r.font.size = Pt(16); r.font.bold = True; r.font.color.rgb = WHITE
        add_text(s, L + 540000, y - 30000, R - L - 540000, dy, [txt], size=15.5, gap=0)
        y += dy
    img = FIG / "overview_pipeline.png"
    w, h = fit(img, R - L, BOTTOM - y - 120000)
    picture(s, img, L + (R - L - w) // 2, y + 120000, w, h)
    add_text(s, L, y + 120000 + h + 40000, R - L, 300000,
             ["Then: differential expression, gene-set enrichment, published subtype signatures, a STRING network, and the drug signatures that oppose the recurrence signature."],
             size=CAP, color=GREY, gap=0)
    notes(s, old_notes(3) + ["The strip is the study in one line: the human line, the orthotopic tumor, the ablation (the frame is from the video on the model slide), the recurrence, and the sequencing of three primary and three recurrent tumors."])

    # ---- 3 LITT ablates the core; recurrence starts at the margin
    s = content_slide(prs, "LITT ablates the core; recurrence starts at the margin")
    fig_and_text(s, FIG / "fig_margin.png", paras(old_text(4, "TextBox 24")) + [
        "So the question is not whether ablation works at the core. It is what the cells that survived the margin look like when they regrow."],
        img_w=6400000, size=16, gap=12, text_top=100000)
    refs(s, old_text(4, "TextBox 31") + "   Diagram: thermal dose falls with distance from the fiber; schematic, not to scale.")
    notes(s, old_notes(4))

    # ---- 4 what that margin is known to do, and what it is not
    s = content_slide(prs, "What that margin is known to do, and what it is not")
    img = FIG / "fig_margin_known.png"
    w, h = fit(img, R - L, BOTTOM - TOP - 900000)
    picture(s, img, L + (R - L - w) // 2, TOP, w, h)
    navy_bar(s, L, TOP + h + 120000, R - L, 760000,
             "Does post-LITT recurrence carry its own transcriptional program, and does that program point to a drug?", size=18)
    refs(s, old_text(5, "TextBox 29"))
    notes(s, old_notes(5))

    # ---- 5 the model, with the ablation video
    s = content_slide(prs, "The model: MRI-guided LITT in a rat xenograft")
    gapc = 360000
    colw = (R - L - 2 * gapc) // 3
    xs = [L, L + colw + gapc, L + 2 * (colw + gapc)]
    heads = ("U251N xenograft", "MRI-guided LITT", "Tumors re-harvested")
    subs = (["orthotopic human glioblastoma line in the athymic RNU/RNU rat"],
            ["thermal ablation of the primary tumor under real-time MR thermometry"],
            ["primary before ablation, recurrent after regrowth, three animals per arm"])
    for x, head, sub in zip(xs, heads, subs):
        head_box(s, x, TOP, colw, 480000, head, sub, size=20)
    for x in xs[:2]:
        arrow(s, x + colw + (gapc - 320000) // 2, TOP + 80000)
    vy = TOP + 1250000
    vs = BOTTOM - vy
    s.shapes.add_movie(str(FIG / "ablation_512.mp4"), Emu(L), Emu(vy), Emu(vs), Emu(vs), poster_frame_image=str(FIG / "video_poster.png"), mime_type="video/mp4")
    x2 = L + vs + GAP
    add_text(s, x2, vy, R - x2, 2600000, paras(old_text(6, "TextBox 33")) + [
        "The model was established and characterized by the Henry Ford group: the ablation is set by MR thermometry, and diffusion imaging shows the lesion as it forms; the video is one such ablation.",
        "Primary tumors were harvested before ablation and recurrent tumors after regrowth, three animals in each arm; every number that follows rests on those six transcriptomes."],
        size=15, gap=8)
    add_text(s, x2, vy + 2700000, R - x2, 700000,
             ["Video: the ablation under diffusion-weighted MRI in this model; click to play (7.9 s). Nagaraja 2021, Acta Neurochirurgica."],
             size=CAP, color=GREY, gap=0)
    refs(s, old_text(6, "TextBox 34"))
    notes(s, old_notes(6) + ["The video is the intra-MRI LITT with diffusion imaging from the 2021 technical note (Ablation video_TNN_18Sept26.pptx); the source clip is 64 x 64 px, re-encoded at 512 px for the slide."])

    # ---- 6 the analysis, explained
    s = content_slide(prs, "The analysis, and how to read the numbers that follow")
    heads = ("Sort", "Test", "Interpret")
    subs = (["human graft reads sorted from rat host reads (xengsort, k = 25)"],
            ["DESeq2 Wald test with ashr shrinkage; FDR < 0.05, |log₂ fold change| > 1"],
            ["GSEA enrichment, GSVA subtypes, STRING network, DSigDB drug signatures"])
    for x, head, sub in zip(xs, heads, subs):
        head_box(s, x, TOP, colw, 480000, head, sub, size=20)
    for x in xs[:2]:
        arrow(s, x + colw + (gapc - 320000) // 2, TOP + 80000)
    add_text(s, L, TOP + 1450000, R - L, BOTTOM - TOP - 1450000, [
        "Enrichment asks whether a whole gene set drifts toward the down- or up-regulated end of the ranked transcriptome. The normalized enrichment score (NES) is that drift; negative means down in recurrence.",
        "Significance is permuted on gene-set membership, not on samples: with three per arm there are only ten ways to split six samples, so a sample-permuted p can never fall below 0.10. Gene-set permutation is the reproducible choice at this size, and nominal p is shown beside FDR q throughout.⁶",
        "Subtypes are scored per tumor by GSVA over the full published Neftel and Garofano signatures. Drugs are found by asking which compound's own expression signature runs opposite to ours, then weighting by predicted blood–brain barrier permeability.",
        "Raw sequencing data are public: GEO accession GSE338105."], size=14.5, gap=8)
    refs(s, old_text(7, "TextBox 33") + "   " + DESIGN)
    notes(s, old_notes(7))

    # ---- 7 PCA
    s = content_slide(prs, "Primary and recurrent tumors separate globally")
    pv = F["permanova"]
    fig_and_text(s, FIG / "chart_pca.png", [
        f"Each point is one tumor's whole transcriptome. The arms sit on opposite sides of the first axis, which carries {F['pc1']} % of all variance; the second carries {F['pc2']} %.",
        f"PERMANOVA puts {100 * pv['r2']:.1f} % of transcriptome-wide variance between the arms (R² = {pv['r2']}, F = {pv['F']}).",
        f"p = {pv['p']:.2f} is the floor of this test, not a weak result: with three per arm there are only ten ways to split six samples, and the observed split is the most extreme of them.",
        "A large effect that this sample cannot make formally significant, said out loud rather than left in the supplement."],
        img_w=5200000, size=15, gap=10, text_top=80000)
    refs(s, DESIGN + " Principal components of the variance-stabilised counts; PERMANOVA on the transcriptome.")
    notes(s, old_notes(8))

    # ---- 8 only 35 genes move
    s = content_slide(prs, "Only 35 genes move, and the biggest movers are artefacts of low counts")
    img = FIG / "chart_volcano.png"
    w, h = fit(img, 6100000, BOTTOM - TOP)
    picture(s, img, L, TOP, w, h)
    x2 = L + w + GAP
    tb = FIG / "table_fold.png"
    tw_, th_ = fit(tb, R - x2, 2100000)
    picture(s, tb, x2, TOP, tw_, th_)
    add_text(s, x2, TOP + th_ + 120000, R - x2, BOTTOM - (TOP + th_ + 120000), [
        f"{F['de_total']} genes pass FDR < 0.05 and a two-fold change: {F['de_up']} up, {F['de_down']} down; 27 are protein coding.",
        "Every one of the six largest fold changes is a readthrough transcript, a U6 snRNA or an unannotated locus on 14–111 mean counts against a median of 232: a handful of reads makes a large ratio.",
        "The protein-coding changes are near two-fold and cluster in matrix and calcium genes, taken up on the network slide.",
        "A list this short does not carry the biology; the rest of the talk asks what the whole transcriptome is doing."],
        size=13, gap=5)
    refs(s, DESIGN + " Wald test with ashr shrinkage; FDR < 0.05 and |log2 fold change| > 1.")
    notes(s, old_notes(9))

    # ---- 9 translation is switched off
    s = content_slide(prs, "One dominant theme: protein synthesis is switched off")
    a, b = FIG / "chart_gsea_sets.png", FIG / "chart_themes.png"
    hh = BOTTOM - TOP - 1250000
    wa, ha = fit(a, R - L, hh); wb, hb = fit(b, R - L, hh)
    scale = min(1.0, (R - L - 250000) / (wa + wb))
    wa, ha, wb, hb = int(wa * scale), int(ha * scale), int(wb * scale), int(hb * scale)
    picture(s, a, L, TOP, wa, ha)
    picture(s, b, R - wb, TOP, wb, hb)
    top = F["gsea_down_top"]
    add_text(s, L, TOP + max(ha, hb) + 80000, R - L, 1200000, [
        f"Left: of {F['gsea_n_sets']:,} gene sets, the six most depleted in recurrence are all protein synthesis or amino-acid stress. Translation initiation (NES {m(top[0]['nes'])}) is the one set in the screen that clears FDR (q = {top[0]['q']:.3f}); the next five are nominally p < 0.001 and point the same way.",
        f"Right: the same direction runs through the whole anabolic programme (ribosome biogenesis {m(nes('ribosome biogenesis'))}, mTORC1 {m(nes('mTORC1 signalling'))}, MYC targets {m(nes('MYC targets'))}) and through energy (glycolysis {m(nes('glycolysis'))}, OXPHOS {m(nes('oxidative phosphorylation'))}, mitochondrial translation {m(nes('mitochondrial translation'))}), while the TCA cycle does not move ({m(nes('TCA cycle'))}).",
        "A cell that has stopped building, not one that has stopped burning."],
        size=13.5, gap=5)
    refs(s, DESIGN + " Broad GSEA, gene-set permutation, 1,000 permutations; nominal p beside FDR q because n = 3 per arm.")
    notes(s, old_notes(10) + ["Right-hand chart, normalized enrichment scores from facts.json#themes: " + " / ".join(paras(old_text(10, "TextBox 25")))])

    # ---- 10 one-sided
    s = content_slide(prs, "The effect is strikingly one-sided")
    fig_and_text(s, FIG / "chart_gsea_landscape.png", [
        f"Every one of the {F['gsea_n_sets']:,} gene sets is a dot: its enrichment score against how far it clears FDR.",
        "Six reach q < 0.25. All six are down-regulated, and all six are translation or stress sets.",
        f"Nothing on the up-regulated side survives correction (best q = {F['gsea_up_best_q']}) or forms an interpretable group.",
        "The recurrent state is defined by what it switches off, not by a new program it switches on."],
        img_w=6500000, size=15, gap=10, text_top=80000)
    refs(s, DESIGN + " Broad GSEA, gene-set permutation, 1,000 permutations.")
    notes(s, old_notes(11))

    # ---- 11 NEW what else moves
    s = content_slide(prs, "What else moves: hypoxia and iron down, cell cycle up")
    fig_and_text(s, FIG / "chart_themes2.png", [
        f"Hypoxic signalling collapses: HIF1 targets and the hypoxia metagene both {m(nes('HIF1 targets'))} ({pval('HIF1 targets')}, {pval('hypoxia metagene')}), angiogenesis {m(nes('angiogenesis'))}. Ablation removes the hypoxic core; what regrows is better perfused.",
        f"Iron uptake and transport falls ({m(nes('iron uptake and transport'))}, {pval('iron uptake and transport')}), as does the senescence programme ({m(nes('senescence'))}, {pval('senescence')}).",
        f"Cell-cycle sets rise nominally: mitotic spindle {m(nes('mitotic spindle'), True)} ({pval('mitotic spindle')}), G2M checkpoint {m(nes('G2M checkpoint'), True)}, E2F targets {m(nes('E2F targets'), True)}. A regrowing tumor dividing again, read as a direction, not a finding: none of these clears FDR (q ≥ 0.6).",
        "Each of these is what a thermally stressed, reperfused margin would be expected to do, which is why the drug output of the same pipeline is worth acting on."],
        img_w=6300000, size=13.5, gap=8, text_top=60000)
    refs(s, DESIGN + " Broad GSEA on the Hallmark, Reactome, Buffa, Semenza and Fridman sets named; nominal p, gene-set permutation; no set on this slide clears FDR (q ≥ 0.6).")
    notes(s, ["New slide, about 20 s. Hypoxia is the cleanest story: HIF1 targets and the Buffa metagene both at -1.79 with nominal p near 0.01, angiogenesis down with them; the ablated core was the hypoxic core. Iron uptake down is worth one sentence because the top drug is an iron chelator, but do not build on it. The cell-cycle sets are up only nominally and the one-sided slide already said nothing up survives correction: say 'dividing again, as a direction'. Cut this slide if running long."])

    # ---- 12 subtypes
    s = content_slide(prs, "Scored against the published classifiers: the astrocyte-like state is lost")
    gb = F["gsva_best"]
    fig_and_text(s, FIG / "chart_subtypes.png", [
        "Each tumor is scored by GSVA for how much of each published signature it expresses: ten Neftel and Garofano signatures, 39 to 50 genes each, the published definitions rather than a hand-picked panel.",
        f"One transition is significant: the astrocyte-like state is lost, {gb['primary']:+.2f} → {gb['recurrent']:+.2f} (p = {gb['p']:.3f}).",
        "The mitochondrial signature moves the same way (−0.58): 34 of its 35 measurable genes are lower in recurrence, the direction OXPHOS and mitochondrial translation took on the enrichment slide. Two independent computations agreeing is the argument.",
        "Mesenchymal and neural-progenitor scores drift but do not reach significance."],
        img_w=6100000, size=13.5, gap=8, text_top=60000)
    refs(s, old_text(12, "TextBox 26") + "   " + DESIGN)
    notes(s, old_notes(12))

    # ---- 13 hubs
    s = content_slide(prs, "The hubs: matrix genes down, calcium genes up")
    img = FIG / "chart_ppi.png"
    tb = FIG / "table_hubs.png"
    tw_, th_ = fit(tb, 4600000, BOTTOM - TOP)
    picture(s, tb, R - tw_, TOP, tw_, th_)
    cw = R - tw_ - GAP - L
    w, h = fit(img, cw, BOTTOM - TOP - 1750000)
    picture(s, img, L + (cw - w) // 2, TOP, w, h)
    add_text(s, L, TOP + h + 60000, cw, 420000,
             ["11 of the 35 genes have a STRING interaction. Navy down, red up; filled, highest degree; line weight, STRING score."],
             size=CAP, color=GREY, gap=0)
    add_text(s, L, TOP + h + 540000, cw, 1250000, [
        "Down: the structural matrix, every collagen, biglycan, prolargin and hemicentin. Up: neuronal calcium handling and a presynaptic scaffold, calbindin, ryanodine receptor 2 and piccolo.",
        "Whether that echoes the neuron–glioma synapse described in high-grade glioma⁹˒¹⁰, or only reflects which cells survived, six animals cannot tell us."],
        size=13.5, gap=5)
    refs(s, old_text(13, "TextBox 27") + "   STRING interactions among the 35 differentially expressed genes.")
    notes(s, old_notes(13))

    # ---- 14 NEW the drug filter
    s = content_slide(prs, "From 93 opposing drug signatures to one candidate")
    fig_and_text(s, FIG / "chart_funnel.png", [
        "DSigDB holds the expression signatures of drugs.¹⁶ The pipeline asks which drug signature runs opposite to the recurrence signature (GSEA on the drug's genes), then weights by predicted brain penetration (ADMET-AI¹¹), with a second barrier model as a check (BOILED-Egg¹²).",
        f"93 signatures opposed ours. Only {F['drug_clinical_n']} belong to a compound with any clinical phase; the rest include phosphine, tributyltin and a flame retardant.",
        f"{len(F['drug_tierA'])} clear both barrier models. Prior art was then read for each: seven have no glioma literature, pentetrazol causes seizures, and amiodarone and diazepam carry published evidence of harm.¹³",
        "The filter is honest about its limits: it keeps things that are not systemic therapies (oxygen, ozone, magnesium) and drops approved drugs with no recorded phase, deferoxamine among them."],
        img_w=5200000, size=13.5, gap=8, text_top=40000)
    refs(s, old_text(14, "TextBox 27") + "   " + REF_DSIGDB + "   Score: |NES| to the power 1.5, weighted by predicted blood–brain barrier permeability.")
    notes(s, ["New slide, about 20 s. Explain the score in one sentence: how strongly the drug's own signature reverses the tumor signature, weighted by predicted barrier penetration. Then the funnel: 93 to 54 to 15 to one. Own the filter's limits before anyone asks."] + old_notes(14)[:1])

    # ---- 15 the candidate
    s = content_slide(prs, "One candidate survives every filter")
    img = FIG / "chart_drugs.png"
    tb = FIG / "table_drugs.png"
    tw_, th_ = fit(tb, 4900000, 2500000)
    picture(s, tb, R - tw_, TOP, tw_, th_)
    w, h = fit(img, R - tw_ - GAP - L, BOTTOM - TOP)
    picture(s, img, L, TOP, w, h)
    dt = F["drug_top"]; d2 = F["drug_tierA"][1]
    add_text(s, R - tw_, TOP + th_ + 120000, tw_, BOTTOM - (TOP + th_ + 120000), [
        f"Ciclopirox scores {dt['score']:.2f} against {d2['score']:.2f} for the runner-up, and ranks first with or without the permeability weight; both barrier models agree.",
        "It is approved (a topical antifungal), and an oral form has completed phase 1 in cancer patients.¹⁴",
        [("Computational predictions: no agent named here has been tested in this model or is approved for this use.", BLUE)]],
        size=13, gap=6)
    refs(s, old_text(14, "TextBox 27") + "   " + REF_MINDEN)
    notes(s, ["Prior art, said aloud rather than shown: " + paras(old_text(14, "TextBox 26"))[2]] + old_notes(14))

    # ---- 16 NEW why ciclopirox could work
    s = content_slide(prs, "Why ciclopirox: the tumor depends on its target")
    dm = {r.split(",")[1]: float(r.split(",")[2]) for r in (ROOT / "depmap" / "depmap_u251_targets.csv").read_text().strip().splitlines()[1:]}
    fig_and_text(s, FIG / "chart_depmap.png", [
        "Ciclopirox chelates iron and inhibits deoxyhypusine hydroxylase (DOHH), the enzyme that activates eIF5A, a factor the ribosome needs to keep elongating.",
        f"Does this tumor depend on that axis? In DepMap's genome-wide CRISPR screen¹⁵, U-251 MG loses fitness when DOHH ({m(dm['DOHH'])}), DHPS ({m(dm['DHPS'])}) or EIF5A ({m(dm['EIF5A'])}) is knocked out, past the −0.5 dependency threshold, and none is a common-essential gene: a selective dependency, and one independent of our signature.",
        "Ribonucleotide reductase, the other ciclopirox target, is essential in every line; the HIF axis and iron handling are not dependencies here.",
        "Delivery: LITT opens the blood–brain barrier at the margin for weeks after ablation², which is the window.",
        [("Public screen data and a computational prediction; nothing has been dosed in this model.", BLUE)]],
        img_w=6300000, size=13, gap=6, text_top=40000)
    refs(s, "2. Cleary RT, et al. Neuro Oncol 2026;28:1649–1661.   " + REF_MINDEN + "   " + REF_DEPMAP)
    notes(s, ["New slide, about 25 s. This is the answer to 'why should an iron chelator work on a cell that has stopped building': ciclopirox also inhibits DOHH, and DepMap says U-251 depends on DOHH, DHPS and EIF5A selectively (Chronos gene effect below -0.5, not common essential). Say plainly that this is a public screen, not our data, and that nothing has been dosed. The barrier window is the delivery argument."])

    # ---- 17 conclusions
    s = content_slide(prs, "Conclusions")
    img = FIG / "chart_arms.png"
    w, h = fit(img, 5500000, BOTTOM - TOP - 1000000)
    picture(s, img, L, TOP, w, h)
    x2 = L + w + GAP
    add_text(s, x2, TOP, R - x2, h, paras(old_text(15, "TextBox 24")), size=13, gap=5)
    navy_bar(s, L, TOP + h + 120000, R - L, BOTTOM - (TOP + h + 120000), old_text(15, "Rounded Rectangle 25"), size=15)
    add_text(s, L, REFS_TOP, R - L, 400000,
             ["Evidence status: investigational; a single cell line in a single animal model, n = 3 per arm; drug candidates are computational predictions. " + DESIGN],
             size=REF, color=GREY, gap=0)
    notes(s, old_notes(15))

    # ---- 18 limitations and next steps
    s = content_slide(prs, "Limitations and next steps")
    half = (R - L - GAP) // 2
    for x, name in ((L, "TextBox 24"), (L + half + GAP, "TextBox 25")):
        block = paras(old_text(16, name))
        head_box(s, x, TOP, half, 460000, block[0], [], size=22)
        add_text(s, x, TOP + 620000, half, BOTTOM - TOP - 620000, block[1:], size=15, gap=10)
    refs(s, DESIGN)
    notes(s, old_notes(16))

    # ---- 19 acknowledgments and data availability
    s = content_slide(prs, "Acknowledgments and data availability")
    add_text(s, L, TOP, 6300000, BOTTOM - TOP, paras(old_text(17, "TextBox 24")), size=15, gap=10)
    q = paras(old_text(17, "Rounded Rectangle 25"))
    navy_bar(s, R - 3600000, TOP + 300000, 3600000, 1300000, q, size=20)
    add_text(s, R - 3600000, TOP + 1800000, 3600000, 1200000,
             ["Sequencing data: GEO GSE338105", "No commercial support. Animal work under Henry Ford Health IACUC protocol 1509, ARRIVE guidelines."],
             size=13, color=GREY, gap=6, align=PP_ALIGN.CENTER)
    notes(s, old_notes(17))

    prs.save(str(out))
    print("wrote", out, "slides:", len(prs.slides))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(Path.home() / "OneDrive/Desktop/CTSpinoPelvic1K-1/CNS2026_Schwing_Abstract418.pptx"))
    a = ap.parse_args()
    out = Path(a.out)
    bak = HERE / "CNS2026_Schwing_Abstract418_2026-09-07_lab_template.pptx"
    if out.exists() and not bak.exists():
        shutil.copy2(out, bak)
        print("backed up the previous file to", bak)
    build(out)
    keep = HERE / "CNS2026_Schwing_Abstract418_CNStemplate.pptx"
    shutil.copy2(out, keep)
    print("copy:", keep)
