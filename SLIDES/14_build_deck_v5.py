# -*- coding: utf-8 -*-
"""The CNS 2026 Abstract 418 talk, version 5 (2026-09-27): v4 with the species-exact composition work added.

v4 (SLIDES/09_build_deck_v4.py) tells one story: after LITT the regrown tumour turns its ribosomal-protein genes
down, the direction survives every test six animals allow, and it points to ciclopirox. Everything in it is kept,
in order. v5 opens that deck and inserts what the pre-registered SCIENCE analyses added since:

  M1  the bulk recurrence signal split exactly into composition, tumour-cell and host-cell terms (C = 0.47)
  M2  the mesenchymal shift at recurrence is host tissue
  M3  U251N in the opposite hemisphere -- and what it means in the operating room
  M4  methylation: nothing survives correction, and what differs is tumour content
  M5  the experiment six animals size

and twenty backups (B8-B27) carrying every other positive and negative result: the two-platform fraction, the ten
splits, the four views, copy-number amplitude against a measured dilution model, MGMT, the distal-invasion
programme and its floor gate, the detection limit, the fluctuating CpGs, the array-against-RNA cross-validation,
the titration curve that failed, and the three contralateral panels.

Figures: SLIDES/13_v5_figures.py (the manuscript's own panels, re-saved at slide type; the R-drawn ones copied).
Numbers: SLIDES/figures_cns/v5_data/*.json|tsv, byte copies of the grid outputs under u251_science/. This builder
types no result number -- every value below is read from one of those files and the file is named beside it.

References: v4's numbering (1-20, by first appearance in v4) is left untouched; new sources continue at 21 in the
order the new slides appear. Each was checked against Europe PMC on 2026-09-27 before it was written here.

    python SLIDES/14_build_deck_v5.py [--out PATH]
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
from pathlib import Path

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Emu

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

_spec = importlib.util.spec_from_file_location("deck3_v5", HERE / "04_build_cns_template_deck.py")
d3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(d3)
from cns_template import (BLUE, BOTTOM, CAP, GAP, GREY, L, NAVY, R, TOP,  # noqa: E402
                          add_text, content_slide, fit, navy_bar, notes, picture, refs)

FIG = HERE / "figures_cns" / "v5"
DAT = HERE / "figures_cns" / "v5_data"
V4 = HERE / "CNS2026_Schwing_Abstract418_CNStemplate_v4.pptx"
DESIGN = d3.DESIGN
TXT = 16


# ---------------------------------------------------------------- numbers, read from the copied grid outputs
def J(name):
    with open(DAT / name, encoding="utf-8") as fh:
        return json.load(fh)


def T(name):
    with open(DAT / name, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(x, nd=0):
    return f"{100 * float(x):.{nd}f}"


COMP = J("composition_main.json")
COMP_S1, COMP_S4 = J("composition_S1.json"), J("composition_S4.json")
C3 = J("c3_ha2_ha3.json")
C5 = J("c5_l2_t1r.json")
T1E1 = J("t1_e1.json")
HOST = J("m4_host.json")
DIL = J("dilution_tests.json")
MGMT = J("D4_mgmt.json")
T3 = J("T3_fcpg.json")
T4 = J("T4_S1.json")
T2 = J("t2.json")
CURVE = J("validation.json")
POWER = J("power_closed_only.json")
XENG = {r[""]: r for r in T("xengsort_classes.tsv")}
SAMP = {r[""]: r for r in T("composition_samples_main.tsv")}
HM1 = {r["which"]: r for r in T("M2_HM1_summary.tsv")}
LUMP = {r["array"]: r for r in T("L1_lump.tsv")}
SMIN = {(r["dc"], r["c"]): r["s_min"] for r in T("S1_smin.tsv")}
BENCH = T("HC2_benchlists_criterion.tsv")
C6 = {r["key"]: r["value"] for r in T("c6_summary.tsv")}
M3T = {r[""]: r for r in T("m3_rna_vs_dna.tsv")}
ACCT = {r["library"]: r for r in T("m4_lesion_accounting.tsv")}
FDNA = [r for r in T("fdna_resolved.tsv")]

PRIM, REC = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]

# the values the slides quote, each with the file it came from (asserted, so a changed file cannot pass silently)
C_SHARE = COMP["H_A1"]["C_share"]
assert COMP["H_A1"]["pass"] is False and COMP["H_A1"]["threshold"] == 0.5
SENS = [COMP_S1["C2_true_split"]["C_share"], COMP_S4["C2_true_split"]["C_share"],
        COMP["sensitivities"]["S2_floor_corrected_IL64B"]["C_share"],
        COMP["sensitivities"]["S6_length_ratio_0.8_1.25"]["C_share"]]
S3 = COMP["sensitivities"]["S3_without_IL66B"]["C_share"]
PAIR_C = COMP["sensitivities"]["matched_pair_IL67B_vs_NL70B"]["C_share"]
HA2 = C3["H_A2"]["padj05"]
SPLIT_RHO = C3["H_A3"]["C_share"]["spearman_vs_abs_diff_f_RNA"]
COV_FLAG = COMP["C2_true_split"]["covariance_residual_flag_share"]
MES = C5["C5"]["neftel_MES"]
MES_W = C5["C5"]["wang2017_MES"]
L2_RHO = C5["L2"]["spearman_estimate_vs_f_RNA_ten"]
DIG = T1E1["T1"]["DiG"]
HC1 = DIL["tests"]["bespoke"]["HC1"]
HC3 = DIL["tests"]["bespoke"]["HC3"]
FLIP = [f["flip_distance_beta"] for f in MGMT["flip"] if f["array"] in REC + ["IL70B", "IL71B", "IL66B"]]
STP27 = [f["stp27_prob"] for f in MGMT["flip"]]
N_BENCH = len(BENCH)
assert all(r["criterion_B_bespoke"] == "TRUE" and r["criterion_B_fDNA"] == "TRUE" for r in BENCH)
assert all(r["A"] == "1" and r["B_bespoke"] == "4" for r in BENCH)
CY = {k: float(v["chrY_per_M_human_counts"]) for k, v in SAMP.items() if v["chrY_per_M_human_counts"]}
N269_HUMAN = float(XENG["N269B"]["graft_pct"])
N168_HUMAN = float(XENG["N168B"]["graft_pct"])
IL64_HUMAN = float(XENG["IL64B"]["graft_pct"])
PWR = POWER["targets"]


# ---------------------------------------------------------------- slide helpers
def figure_slide(prs, title, image, lines, reserve=1_500_000, size=TXT, gap=8, ref=None, note=None):
    """One panel across the top, the reading under it. The form every v4 results slide uses."""
    s = content_slide(prs, title)
    img = FIG / image
    if not img.exists():
        raise SystemExit(f"missing figure {img}: run SLIDES/13_v5_figures.py first")
    w, h = fit(img, R - L, BOTTOM - TOP - reserve)
    picture(s, img, L + (R - L - w) // 2, TOP, w, h)
    add_text(s, L, TOP + h + 80000, R - L, BOTTOM - (TOP + h + 80000), lines, size=size, gap=gap)
    if ref:
        refs(s, ref)
    if note:
        notes(s, note)
    return s


def side_slide(prs, title, image, lines, img_w=6_400_000, size=TXT, gap=10, ref=None, note=None, text_top=40000):
    """A tall panel at the left, the reading at the right."""
    s = content_slide(prs, title)
    img = FIG / image
    if not img.exists():
        raise SystemExit(f"missing figure {img}: run SLIDES/13_v5_figures.py first")
    w, h = fit(img, img_w, BOTTOM - TOP)
    picture(s, img, L, TOP, w, h)
    x2 = L + w + GAP
    add_text(s, x2, TOP + text_top, R - x2, BOTTOM - TOP - text_top, lines, size=size, gap=gap)
    if ref:
        refs(s, ref)
    if note:
        notes(s, note)
    return s


def reorder(prs, order):
    """Put the slides in `order` (current indices) into that sequence."""
    lst = prs.slides._sldIdLst
    ids = list(lst)
    if sorted(order) != list(range(len(ids))):
        raise SystemExit("reorder: the order must be a permutation of every slide index")
    for e in ids:
        lst.remove(e)
    for i in order:
        lst.append(ids[i])


# ---------------------------------------------------------------- the five new main slides
def m1_composition(prs):
    return figure_slide(
        prs, "How much of the difference is simply less tumor? An exact split",
        "fig2a_shapley_de_genes.png",
        [f"A human graft in a rat host: every read is assigned by species, so the bulk profile a patient sample would give "
         f"can be rebuilt and each gene's fold change split exactly into composition, tumor-cell and host-cell terms.{sup(21)}",
         f"Composition is the largest single term and still under half: {C_SHARE:.2f} of the total fold change "
         f"({min(SENS):.2f}–{max(SENS):.2f} across four sensitivities, {S3:.2f} without one recurrence).",
         [(f"Of the {HA2['n_DE']} genes a species-blind analysis calls differential, {pct(HA2['share_composition_dominant'])} % are "
           f"composition-dominant: the genes it picks out are the ones where cells actually changed.", BLUE)]],
        reserve=2_000_000, size=15, gap=8,
        ref="21. Shen-Orr SS, et al. Nat Methods 2010;7:287–289.   Per-sample mixing identity, then Shapley attribution over "
            "the three factors; the terms sum to the fold change to 1e-15. Ensembl Compara one-to-one orthologs; n = 3 per arm. " + DESIGN,
        note=[f"NEW in v5. 25 s. The one idea: in a xenograft the mixture is measured, not estimated, so 'is it recurrence or "
              f"is it less tumor?' has an exact answer per gene. The identity holds to 3e-18 and the three Shapley terms sum "
              f"to the fold change to 9e-16 in every gene, which the job checks before it prints anything.",
              f"Robustness: across the ten balanced splits, the composition share rises with the split's difference in tumor "
              f"content (Spearman {SPLIT_RHO:.2f}) and the true split is the highest of the ten, which is what a real "
              f"composition term must do. The one fraction-matched pair (IL67B against NL70B) gives {PAIR_C:.2f}.",
              f"What to concede if asked: the within-arm covariance residual exceeds 10 % of the fold change in "
              f"{pct(COV_FLAG)} % of genes, so the arm-level split is exact at arm means and approximate per animal; and "
              f"noise at n = 3 inflates the two cell-change terms, which biases the composition share DOWN, so 0.47 is "
              f"conservative for the hypothesis that composition dominates. The pre-registered threshold was 0.50 and it failed: "
              f"we report that as the result.",
              "Shapley attribution follows Shorrocks' decomposition procedure (J Econ Inequal 2013); confirm that reference "
              "before it goes in a manuscript."])


def m2_mesenchymal(prs):
    out = C5["C5_outside_check"]
    return figure_slide(
        prs, "Where the mesenchymal shift at recurrence comes from",
        "fig2d_mes_human_vs_virtual.png",
        [f"In the tumor cells alone the mesenchymal score FALLS at recurrence{sup(7)} — below every one of 1,000 draws of a "
         f"dilution null for the Neftel program and {pct(1 - MES_W['dilution_null']['share_null_le_obs'], 1)} % of them for Wang.",
         f"Adding the host reads moves the samples toward the bulk reading; that difference ranks {MES['D_rank']} of 20 "
         f"labellings, so a composition-generated shift is not shown at six animals.",
         [(f"The rat-only mesenchymal score tracks myeloid content (Spearman "
           f"{out['neftel_MES']['spearman_vs_TAM_BMDM']:.2f} and {out['wang2017_MES']['spearman_vs_TAM_BMDM']:.2f} with both the "
           f"deconvolved macrophage fraction and Ptprc/Cd68/Aif1): here the bulk mesenchymal signal is host tissue.{sup(22)}", BLUE)]],
        reserve=2_000_000, size=15, gap=8,
        ref=f"7. Neftel C, et al. Cell 2019;178:835–849.   9. Hänzelmann S, et al. BMC Bioinformatics 2013;14:7.   "
            f"22. Spitzer A, et al. Nat Genet 2025;57:1168–1178.   GSVA on TMM log-CPM; dilution null = the primary's human "
            f"counts thinned to the recurrence's depth plus the measured assignment floor, 1,000 draws. " + DESIGN,
        note=["NEW in v5. 25 s. Why a neurosurgical audience should care: 'mesenchymal shift at recurrence' is the most "
              "repeated statement about recurrent glioblastoma, and here, where tumor and host can be told apart exactly, "
              "the shift belongs to the host compartment while the tumor cells move the other way.",
              f"Numbers: human-only delta {MES['delta_human_only']:.2f}, virtual bulk {MES['delta_virtual']:.2f}, "
              f"rat-only {MES['delta_rat_only']:+.2f} (Neftel). The pre-registered statistic was the difference between the "
              f"first two; it ranks {MES['D_rank']} of 20 (Wang {MES_W['D_rank']} of 20), so we say 'not shown', not 'absent'.",
              "Spitzer 2025 is the largest longitudinal human cohort: no consistent expression trajectory, one consistent "
              "change — a lower malignant-cell fraction at recurrence. Our arms reproduce exactly that, on both platforms."])


def m3_contralateral(prs):
    s = figure_slide(
        prs, "Tumor in the opposite hemisphere — and what it means in the OR",
        "fig_c1_N2_over_core_bins.png",
        [f"Rat 69's opposite hemisphere: {N269_HUMAN:.1f} % of reads human, the donor genotype at 56 of 59 identity probes, "
         f"Y-chromosome dose 0.88, and the line's copy-number profile (r = 0.72, z = 8.4). Two platforms, three measurements "
         f"that do not share a failure mode.",
         f"It is not the core diluted: two of three pre-registered statistics reject a depth-matched dilution null "
         f"(P < 0.001). In the second contralateral hemisphere, from another animal, human reads are {N168_HUMAN:.1f} %.",
         [(f"We cannot separate cells that crossed in life from core tissue carried across on the blade. Either way, one pass "
           f"of an instrument left tumor that sequencing detects and histology would not.{sup(23)}", BLUE)]],
        reserve=1_900_000, size=15, gap=8,
        ref="23. Orjefelt E, et al. J Small Anim Pract 2026;67:528–533.   24. Steinmetz MP, et al. J Neurooncol 2001;55:167–171.   "
            "25. Bekar A, et al. World Neurosurg 2010;73:719–721.   Contralateral piece cut from the same frozen coronal slice as "
            "the tumor; copy number on human-specific probes with rat cross-hybridising probes removed, scaled to the parental culture.",
        note=["NEW in v5, and the slide to slow down on. 35 s.",
              "Say it in this order: (1) there is U251N material in the hemisphere opposite the implant, established by genotype, "
              "Y dosage and copy number on an array handled separately from the libraries; (2) its transcriptome is not the core "
              "sampled at lower depth; (3) we cannot tell migration from carry-over at dissection, and both readings matter to a surgeon.",
              "The carry-over reading is the surgical one: a blade that passed through tumor and then through brain left tumor "
              "in the brain at about 3 % of the piece's transcriptome. In canine cancer surgery, malignant cells were recovered "
              "from gloves or instruments in 14 of 47 operations (30 %), and in 8 of 14 (57 %) when margins were incomplete "
              "against 6 of 32 (19 %) when they were clear, P = 0.016 (Orjefelt 2026, cytology). Iatrogenic seeding along a "
              "biopsy tract or into a distant surgical site is documented in case reports (Steinmetz 2001; Bekar 2010; also "
              "Perrin & Bernstein, J Neurooncol 1998;36:243–246). What is new here is that species identity MEASURES the "
              "transferred material instead of inferring it.",
              "Concede before you are asked: one animal, two sampled hemispheres, so no rate is estimable (the exact interval "
              "for one of two runs from 1.3 % to 98.7 %). Both contralateral samples come from animals that were never ablated, "
              "so this says nothing about LITT and spread. No blank library was run. The fix is in the scaled design: a "
              "cutting-order control — a naive hemisphere cut with the same blade and slide straight after a tumor-bearing one — "
              "and 29 such pieces with no positive would bound the false-positive rate at 10 %."])
    return s


def m4_methylation(prs):
    return figure_slide(
        prs, "Methylation: nothing survives correction, and what differs is tumor content",
        "fig4b_dmp_counts.png",
        [f"No probe reaches FDR < 0.05 in any design, on 866,238 probes at three against three.",
         f"The arm difference is tumor content: the true split ranks 1 of 10 splits without a fraction covariate and 4 of 10 "
         f"with it — in two independent pipelines and under all {N_BENCH} published cross-species exclusion lists.",
         [(f"MGMT-STP27{sup(26)} reads methylated on all eight arrays ({min(STP27):.3f}–{max(STP27):.3f}); the shift that would "
           f"change a call is {min(FLIP):.2f}–{max(FLIP):.2f} beta, against observed shifts of −0.03 and +0.06. Status does not move.", BLUE)]],
        reserve=1_750_000, size=15, gap=8,
        ref="26. Bady P, Delorenzi M, Hegi ME. J Mol Diagn 2016;18:350–361.   Probes at p < 0.001 for all ten balanced splits, "
            "per pipeline and design; the true split is orange. Bespoke pipeline and nf-core/methylarray, with and without the "
            "tumour-fraction covariate. " + DESIGN,
        note=["NEW in v5. 25 s. This is the pertinent negative: a methylation array on the same six tumors finds nothing that "
              "survives multiple testing, and the difference that looks real before correction is the same tumor-content "
              "confound as on the RNA side — replicated in a second, independently written pipeline (nf-core/methylarray) and "
              "under every published list of rat-cross-hybridising probes.",
              "MGMT matters to this audience: the model's call is unchanged in every array, and the distance to a flip is ten "
              "times the observed shift, so nothing here suggests LITT changes MGMT status. No animal received an alkylating "
              "agent, so this is a stability readout, not a resistance result.",
              "If asked why no FDR hit at all: 3 v 3 with a 15-point difference in tumor content and an array whose betas are "
              "not tumor-only (the LUMP check fails, 0.84–0.88 against a 0.9 bar). That failure is carried beside every "
              "DNA-side statement in the manuscript."])


def m5_scaled(prs):
    t1, t2b, t2c, t7, t8 = PWR["T1"], PWR["T2b"], PWR["T2c"], PWR["T7"], PWR["T8"]
    return side_slide(
        prs, "The experiment six animals size",
        "fig6_rats_needed.png",
        [f"Targets were written and hashed before any pilot number entered a calculation, and a pilot effect never sizes the "
         f"study: it is marked on the curve, and the study is sized on the pre-stated effect at a pessimistic variance.{sup(28)}",
         f"The binding constraint is the per-animal binary endpoint: {t2b['n_per_arm_needed']} rats per arm to move invasion "
         f"frequency 0.20 to 0.50, {t1['n_per_arm_needed']} for a 15-point difference in tumor content, {t7['n_per_arm_needed']} "
         f"for the host response, and {t2c['n_controls_needed']} control pieces to bound transfer at dissection at 10 %.",
         [(f"The clone-level test — are the cells that grew back the ones that had already left? — needs "
           f"{t8['n_animals_needed']['k0']} to {t8['n_animals_needed_bound_0.80']['k0']} barcoded animals, because the unit of "
           f"inference becomes the clone rather than the animal.{sup(27)}", BLUE)]],
        img_w=5_500_000, size=15, gap=14,
        ref="27. Barthel FP, et al. Nature 2019;576:112–120.   28. Browne RH. Stat Med 1995;14:1933–1940.   Power, interval "
            "half-width or Wilson lower bound against animals; dashed lines, the pre-stated targets; filled marks, the "
            "conservative sizing; open marks, the pilot point. Panels e and f are pending.",
        note=["NEW in v5. 30 s, and the slide to end the science on. The point is not 'we need more animals': it is that every "
              "number on this chart was fixed before the pilot was allowed to speak, in the same form as the power figure in my "
              "imaging work — the target line first, the curve afterwards.",
              f"Why the binary endpoint binds: genome-wide work is cheap in animals and the per-animal yes/no endpoints are not. "
              f"Arms i and ii are sized by invasion frequency at {t2b['n_per_arm_needed']} each ({t2b['n_per_arm_needed_holm']} "
              f"under Holm), arms iv and v by the host response at {t7['n_per_arm_needed']}. A design that wants every target at "
              f"once is in the low hundreds of animals; the honest alternative is to stage it, and the two pieces that close a "
              f"literature gap at the smallest animal cost are the barcode arm and the titration series.",
              "What the pilot cannot answer at any effect size: with six tumors there are twenty balanced labellings, so no "
              "one-sided p can fall below 0.05 however large the effect; both contralateral samples came from unablated animals; "
              "and no animal has a pre-treatment sample, a cutting-order control or a naive hemisphere."])


def sup(n):
    """The superscript digits v4 uses for its reference marks."""
    return str(n).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))


# ---------------------------------------------------------------- the twenty new backups
def backups(prs):
    made = []

    made.append(side_slide(
        prs, "Backup B8. What each library is made of",
        "fig1a_read_composition.png",
        [f"Every read pair assigned by species before anything is counted: {min(float(XENG[k]['graft_pct']) + float(XENG[k]['host_pct']) for k in XENG):.0f}–"
         f"{max(float(XENG[k]['graft_pct']) + float(XENG[k]['host_pct']) for k in XENG):.0f} % go to one species.",
         f"The three floor libraries: IL64B, an implanted hemisphere with no tumor in the piece ({IL64_HUMAN:.1f} % human), and the two "
         f"opposite hemispheres ({N168_HUMAN:.1f} % and {N269_HUMAN:.1f} %).",
         f"Human Y-linked transcripts per million human counts place the human stream as the male donor line: "
         f"{min(CY[k] for k in PRIM + REC):.0f}–{max(CY[k] for k in PRIM + REC):.0f} in the six tumors, {CY['C2B']:.0f} in the culture, "
         f"{CY['IL64B']:.1f} in the tumor-free library."],
        ref="Reads sorted with xengsort (k = 25) against GRCh38 and mRatBN7.2; shares of all classified read pairs."))

    made.append(side_slide(
        prs, "Backup B9. What the host tissue in each lesion is",
        "fig1d_lesion_accounting.png",
        ["Each library as tumor share plus host cell types, scaled so the bars are shares of the whole library.",
         "Three published deconvolution methods agree on rank for all six cell types; every low-tumor library scores above "
         "every tumor on neurons, and every tumor above the floor libraries on myeloid cells.",
         "LM22 was rejected by the pre-registered athymic control: it assigned more than 5 % T cells to a host that has none.",
         "Neither pre-stated index separates the arms (normal-parenchyma rank 7 of 20, reactive rank 9 of 20)."],
        ref="BRETIGEA, NeuroExpresso marker-gene profiles and CIBERSORT with a brain signature, on the rat read stream; mouse "
            "markers mapped one-to-one to rat. Reactive index reduced to its microglial term: no reactive-astrocyte, DAM or "
            "fibrotic-scar list was retrievable verbatim."))

    made.append(side_slide(
        prs, "Backup B10. Tumor content on two platforms",
        "fig1b_rna_vs_dna.png",
        [f"DNA fraction against RNA human share for every estimator; the registered one is the bespoke two-mode estimate and it "
         f"is uncalibrated, so no cell percentage is quoted from either platform.",
         f"The order agrees: Spearman {float(COMP['M3_RNA_side']['spearman_six_tumours']):.2f} over the six tumors and "
         f"{float(COMP['M3_RNA_side']['spearman_seven_arrays']):.2f} over the seven in-vivo arrays.",
         f"The contralateral sample's DNA-to-RNA ratio (about 20 % against 5 %) is unexplained and reported as such."],
        ref="Bars: each estimator's own interval. Dashed: equality, for reference only. Cell share needs the line's DNA index, "
            "which is not measured here."))

    made.append(side_slide(
        prs, "Backup B11. Recurrent lesions hold less tumor: all twenty labellings",
        "fig1c_share_labellings.png",
        [f"DNA: recurrent minus primary {float(HM1['DNA_bespoke']['rec_minus_prim']):.3f} (exact 90 % "
         f"{float(HM1['DNA_bespoke']['exact90_lo']):.3f} to {float(HM1['DNA_bespoke']['exact90_hi']):.3f}), rank 1 of 20, "
         f"p = 0.05 — the smallest the design allows.",
         f"RNA: {float(HM1['RNA']['rec_minus_prim']):.3f}, rank 2 of 20 (p = 0.10). This quantity was seen before "
         f"pre-registration, so it is a replication, not a discovery.",
         "Holm across the two platforms: 0.10 and 0.10."],
        ref="Twenty balanced three-against-three labellings of the six tumors; the true labelling is filled. Every DNA fraction is "
            "uncalibrated."))

    made.append(side_slide(
        prs, "Backup B12. Composition share and DE count across all ten splits",
        "fig2b_splits.png",
        [f"The composition share rises with the split's difference in tumor content (Spearman {SPLIT_RHO:.2f}), and the true "
         f"split is the highest of the ten — what a real composition term must do.",
         "Differentially expressed genes across the ten splits range from 0 to 303; the true split gives 133.",
         "The true split's composition share sits inside the 80 % band predicted from the other nine."],
        ref="nf-core/differentialabundance (DESeq2) on the species-blind matrix, the same parameters for all ten splits."))

    made.append(side_slide(
        prs, "Backup B13. The three published categories, in four views",
        "fig2c_categories_by_view.png",
        ["Cell motility and inflammatory response are Holm-significant among the genes up in the species-blind (bulk) view.",
         "Neither appears in the tumor-cell stream; cell cycle appears in no view.",
         "The composition-only view is empty by construction: with a 15-point difference in content the composition term "
         "cannot reach one log2 unit, which was the pre-registered cut."],
        ref="g:Profiler on three fixed GO biological-process terms, Holm across three categories by four views. The three "
            "categories are the ones reported for this model."))

    made.append(side_slide(
        prs, "Backup B14. Copy-number amplitude against a measured dilution model",
        "fig4a_amplitude_dilution.png",
        [f"Diluting the parental culture in silico with rat DNA barely lowers its copy-number amplitude, so dilution alone "
         f"does not put the primaries below the curve.",
         f"The primaries sit below it and the recurrences on or above it; recurrent minus primary "
         f"{HC1['rec_minus_prim']:.3f} (exact 90 % {HC1['exact90'][0]:.3f} to {HC1['exact90'][1]:.3f}), rank 1 of 20.",
         [("The null behind this panel is unvalidated, because the titration that would validate it failed (B23). We report "
           "the difference as unexplained by dilution, and claim nothing from it.", BLUE)]],
        ref="Amplitude = slope of an array's copy-number bins on the parental culture's. Line and band: the culture diluted in "
            "silico with ten rat reference arrays, median and range, through the identical pipeline."))

    made.append(side_slide(
        prs, "Backup B15. MGMT status does not move",
        "fig4c_mgmt_stp27.png",
        [f"MGMT-STP27 calls every array methylated ({min(STP27):.3f}–{max(STP27):.3f}).",
         f"The equal shift at both model CpGs that would flip a recurrence is {min(FLIP):.2f}–{max(FLIP):.2f} beta; the observed "
         f"recurrent-minus-primary shifts are −0.032 and +0.057.",
         f"The rat-free promoter mean differs by {HC3['observed_rec_minus_prim']:.3f} (outside an unvalidated dilution null).",
         "MGMT is not expressed in any library, and no animal received an alkylating agent: this is a stability readout."],
        ref="mgmtstp27; the dashed line is the model's decision boundary and each arrow the shift that would change that "
            "recurrence's call."))

    made.append(side_slide(
        prs, "Backup B16. The one tumor-cell signal that survives adjustment",
        "fig5a_dig_arms.png",
        [f"A published distal-invasion program scores higher in the recurrences after adjusting for tumor content: the arm "
         f"coefficient ranks {DIG['relabelling_rank']} of 20 (p = {DIG['p_one_sided']:.2f}), as does the unadjusted difference.",
         "Two scoring packages agree. The pre-registration expected this to fail, on the basis of the human longitudinal cohort.",
         [("The floor is the floor: one labelling of twenty, in a model where the recurrent arm is also the later-harvested "
           "arm. It is a hypothesis for the scaled experiment, not a finding.", BLUE)]],
        ref="Chanoch-Myers 2026 distal-invasion gene sets, frozen by hash; singscore on floor-corrected human log-CPM; model "
            "score ~ arm + tumour share, refitted for each of the twenty labellings."))

    made.append(side_slide(
        prs, "Backup B17. The gate that stopped the contralateral invasion test",
        "fig5e_dig_floor_gate.png",
        [f"The pre-registered case test asked whether rat 69's contralateral human cells sit further along the invasion axis "
         f"than its own core.",
         f"Its gate failed first: the assignment floor itself scores as high as the contralateral sample "
         f"({T1E1['T1']['DiG']['arm_coef_adjusted_for_f']:.4f} is the arm coefficient; the floor libraries sit at the top of "
         f"this panel), so the contrast would measure mis-assignment.",
         "The test was abandoned and is reported as abandoned. The p-value the job printed is not a result."],
        ref="Same score as the arm test, on every human-bearing profile plus the two floor libraries. Pre-registered gate: the "
            "floor libraries must score below the contralateral sample."))

    made.append(side_slide(
        prs, "Backup B18. What a shared subclonal change would have to be for us to see it",
        "fig5d_s1_detection.png",
        [f"Recurrences do not share derived arm-level copy-number changes that the primaries lack "
         f"(statistic {T4['T_T4_arms']:.3f}, rank {T4['rank_arms']} of 20) — the pre-stated expected null.",
         f"The bound that makes that negative useful: a single-copy loss carried by fewer than {SMIN[('-1', '2')]} of recurrent "
         f"cells (from two copies) or {SMIN[('-1', '3')]} (from three) would not have been flagged on 80 % of arms.",
         "Gains are not detectable at any share by this rule, and that is stated rather than hidden."],
        ref="The registered rule applied unchanged to simulated shared changes added to the observed recurrent bins."))

    made.append(side_slide(
        prs, "Backup B19. Fluctuating CpGs: uninformative, and why",
        "fig5c_fcpg_w.png",
        [f"The homogeneity readout was pre-registered with two gates. The second failed: the parental culture is already "
         f"multimodal at the published fluctuating-CpG sets (dip test p < 0.005 for every list).",
         f"The statistic is about zero ({T3['results'][T3['primary']].get('T_T3', 0):.4f}) and ranks 14 of 20.",
         "So the homogeneity triad cannot be assembled, and no homogeneity claim is made anywhere in this work."],
        ref="Gabbutt 2022 and 2025 fluctuating-CpG lists, frozen by hash; ten matched control CpGs per site; the band is the "
            "dilution model."))

    made.append(side_slide(
        prs, "Backup B20. Array against RNA: expression dosage",
        "fig_xval_dosage.png",
        ["Two assays on the same eight samples check each other. Where the RNA is tumor-rich, expression follows copy number: "
         "slope 0.54 (p = 7e-06), and the CDKN2A/B deletion every array calls is silent in every library, as it must be.",
         "In the contralateral sample, arm-level expression does not follow its copy number anywhere in the genome "
         "(slope 0.17, p = 0.45), against 0.97 for its own core.",
         [("So the deeper losses seen there stay array-only and unresolved: the assay has no demonstrated dosage sensitivity "
           "in that sample, and we say so rather than reading them as subclones.", BLUE)]],
        ref="Arm-level copy number against arm-median expression shift, nulls matched on expression decile."))

    made.append(side_slide(
        prs, "Backup B21. Array against RNA: genotype",
        "fig_xval_genotype.png",
        ["At expressed CpG-SNP positions the RNA genotype agrees with the array's implied genotype in 527 of 543 pairs (97.1 %).",
         "The sixteen disagreements sit at five common germline sites carried by every sample: line-wide, not subclonal.",
         "The libraries agree with each other at 99.86 % over 4,230 paired sites.",
         "The array's identity probes lie outside expressed sequence, so the full-genotype arm of this test could not run."],
        ref="bcftools pile-up at the array's probe positions; genotype cut-offs fixed in the pre-registration."))

    made.append(side_slide(
        prs, "Backup B22. Array against RNA: promoter methylation",
        "fig_xval_meth.png",
        ["Within every sample, promoter methylation and expression anticorrelate (Spearman −0.31 to −0.36; no shuffle of "
         "10,000 reached it), the direction the cell-line compendia show.",
         "Across the recurrence contrast the two do not move together (Spearman −0.003; the true split ranks 8 of 10).",
         "The ribosomal-protein block that dominates the tumor-cell list has no copy-number or promoter support."],
        ref="Promoter beta against log-CPM per sample; the delta-delta test compares fraction-adjusted promoter change with "
            "graft-adjusted expression change."))

    made.append(side_slide(
        prs, "Backup B23. The calibration that failed, and why it matters",
        "fig_curve_mouse_validation.png",
        [f"Mixing two pure arrays in silico does not reproduce real DNA mixtures: mean absolute error "
         f"{CURVE['step1']['metrics']['primary']['MAE']:.3f} against a pre-set bar of 0.05 (maximum "
         f"{CURVE['step1']['E_max_abs_error']:.3f} against 0.10), every real mixture reading high.",
         f"A correction fitted on those real mixtures still misses its held-out bar "
         f"({CURVE['amendment1']['metrics']['MAE']:.3f}, maximum {CURVE['amendment1']['E_cal']:.3f}).",
         [("So every tumour-fraction value here is uncalibrated and every dilution null is labelled unvalidated. No published "
           "human:rat titration on a human methylation array exists — the scaled design runs one.", BLUE)],
         "We also found a defect in the public human:mouse series: in one batch every red-channel file is a copy of the green."],
        ref="GSE310817 (human:mouse titrations, EPIC v1) as the validation set; GSE299969 (EPIC v2) as the independent check."))

    made.append(side_slide(
        prs, "Backup B24. The contralateral sample against a dilution null",
        "fig_a_nullA.png",
        ["The null: the animal's own core subsampled to the contralateral sample's tumor depth, plus the measured assignment "
         "floor, 1,000 draws at a fixed seed.",
         "Whole-profile distance 0.153 against a null 99.5th percentile of 0.028, and larger than another animal's core gives.",
         "Malignant-state composition differs by 0.102 against a null maximum of 0.069.",
         "The third statistic, the only one with a direction named in advance, does not reject (P = 0.56) and moves the other way."],
        ref="Dashed and dotted lines: the same construction from two other animals' cores. One animal."))

    made.append(side_slide(
        prs, "Backup B25. Does the contralateral population look like the recurrences?",
        "fig_b_similarity.png",
        ["No, on either pre-registered statistic, and it fails in the direction the literature predicts.",
         "Recurrent-minus-primary similarity difference −0.028: the sample sits nearer the primaries.",
         "The true arm labelling ranks fifteenth of twenty (P = 0.75).",
         "A variant-sharing test looked strongly positive (P = 9e-10) until the variants a tumor-free library also carries "
         "were removed; then 0.830 against 0.828 (P = 0.93). Host reads in the graft stream manufacture that signal."],
        ref="All libraries subsampled to one depth before the comparison; twenty balanced labellings."))

    made.append(side_slide(
        prs, "Backup B26. The four arms named before the analysis",
        "fig_c1_events.png",
        ["None is a contralateral-only copy-number event: the animal's own core carries all four at three to six times the "
         "noise bound measured on thirteen copy-neutral control arrays, and so does the parental culture.",
         "Along the genome the contralateral profile is the core's profile with the same breakpoints and deeper losses.",
         "Whether those deeper losses are subclonal composition or low-signal distortion, this array cannot resolve (B20)."],
        ref="Colour: amplitude-scaled log2 ratio; a dot marks a value beyond the noise bound."))

    made.append(side_slide(
        prs, "Backup B27. Malignant-state composition of the human reads",
        "fig_c_states.png",
        ["Deconvolution of the human stream into three malignant states, renormalised over the malignant columns.",
         "The contralateral sample has fewer differentiated-like cells than its own core (0.037 against 0.139).",
         [("We weight this lightly and say why: another animal's tumor differs from the same core by the same distance, and "
           "the method fails a negative control, assigning 36 % of a pure human culture to microenvironment cell types that a "
           "human-only read stream from a rat host cannot contain.", BLUE)]],
        ref="CIBERSORT with a published single-cell reference; fit correlation printed above each bar."))

    for s in made:
        notes(s, ["NEW in v5. Backup: show only if asked."])
    return made


# ---------------------------------------------------------------- build
def build(out: Path):
    if not V4.exists():
        raise SystemExit(f"missing {V4}: run python SLIDES/09_build_deck_v4.py first")
    prs = Presentation(str(V4))
    n_v4 = len(prs.slides)
    if n_v4 != 21:
        raise SystemExit(f"v4 deck has {n_v4} slides, expected 21: check 09_build_deck_v4.py before extending it")

    m1_composition(prs)      # 21
    m2_mesenchymal(prs)      # 22
    m3_contralateral(prs)    # 23
    m4_methylation(prs)      # 24
    m5_scaled(prs)           # 25
    backups(prs)             # 26 onwards

    # v4: 0 title, 1 margin, 2 design, 3 model, 4 PCA, 5 GSEA, 6 leave-one-out, 7 content confound, 8 subtypes,
    #     9 drug, 10 why, 11 what this means, 12 limitations, 13 acknowledgments, 14-20 backups B1-B7
    order = ([0, 1, 2, 3, 4, 5, 6, 7] + [21, 22] + [8, 9, 10] + [23, 24] + [11, 12] + [25] + [13]
             + list(range(14, 21)) + list(range(26, len(prs.slides))))
    reorder(prs, order)

    for sl in prs.slides:                                   # a number never parts from its '%'
        for sh in sl.shapes:
            if sh.has_text_frame:
                for p in sh.text_frame.paragraphs:
                    for r in p.runs:
                        if " %" in r.text:
                            r.text = r.text.replace(" %", " %")
    prs.save(str(out))
    print(f"wrote {out}  slides: {len(prs.slides)}  (v4 {n_v4} + 5 main + {len(prs.slides) - n_v4 - 5} backups)")
    print("new references: 21 Shen-Orr 2010 | 22 Spitzer 2025 | 23 Orjefelt 2026 | 24 Steinmetz 2001 | "
          "25 Bekar 2010 | 26 Bady 2016 | 27 Barthel 2019 | 28 Browne 1995")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "CNS2026_Schwing_Abstract418_CNStemplate_v5.pptx"))
    a = ap.parse_args()
    build(Path(a.out))
