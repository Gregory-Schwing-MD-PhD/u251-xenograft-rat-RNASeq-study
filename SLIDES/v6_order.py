# -*- coding: utf-8 -*-
"""The running order of deck v6, and the slides v6 adds to v4.

Greg's review of v5 (2026-09-27): keep the contralateral slide in the talk, promote more of the backup material into
it in a logical order, remember the talk is six minutes, do not say the instrument reading is what was measured, and
stop putting multi-panel figures and paragraphs on slides.

So: every slide below carries ONE panel, the text is at most three short lines, the deck marks a 6-minute CORE, and
everything else is the expanded talk with a cut order (SLIDES/STORY_v6.md). Every instrument sentence is conditional.

v4's fourteen main slides and seven backups are kept exactly as they are and re-ordered around the new ones; v4's two
side-by-side-figure slides (PCA + volcano, GSEA + leading edge) are left alone because both figures there are already
large and legible - the defect Greg reported was in v5's additions.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from cns_template import BLUE

HERE = Path(__file__).resolve().parent
DAT = HERE / "figures_cns" / "v5_data"


def J(n):
    with open(DAT / n, encoding="utf-8") as fh:
        return json.load(fh)


def T(n):
    with open(DAT / n, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pct(x, nd=0):
    return f"{100 * float(x):.{nd}f}"


def sup(n):
    return str(n).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))


# ---------------------------------------------------------------- the numbers, read from the grid outputs
COMP, COMP_S1, COMP_S4 = J("composition_main.json"), J("composition_S1.json"), J("composition_S4.json")
C3, C5, T1E1 = J("c3_ha2_ha3.json"), J("c5_l2_t1r.json"), J("t1_e1.json")
HOST, DIL, MGMT = J("m4_host.json"), J("dilution_tests.json"), J("D4_mgmt.json")
T3, T4, CURVE, POWER = J("T3_fcpg.json"), J("T4_S1.json"), J("validation.json"), J("power_closed_only.json")
XENG = {r[""]: r for r in T("xengsort_classes.tsv")}
SAMP = {r[""]: r for r in T("composition_samples_main.tsv")}
HM1 = {r["which"]: r for r in T("M2_HM1_summary.tsv")}
C6 = {r["key"]: r["value"] for r in T("c6_summary.tsv")}
LUMP = {r["array"]: float(r["LUMP"]) for r in T("L1_lump.tsv")}
SMIN = {(r["dc"], r["c"]): r["s_min"] for r in T("S1_smin.tsv")}
BENCH = T("HC2_benchlists_criterion.tsv")
ACCT = {r["library"]: r for r in T("m4_lesion_accounting.tsv")}

C_SHARE = COMP["H_A1"]["C_share"]
SENS = [COMP_S1["C2_true_split"]["C_share"], COMP_S4["C2_true_split"]["C_share"],
        COMP["sensitivities"]["S2_floor_corrected_IL64B"]["C_share"],
        COMP["sensitivities"]["S6_length_ratio_0.8_1.25"]["C_share"]]
S3 = COMP["sensitivities"]["S3_without_IL66B"]["C_share"]
HA2 = C3["H_A2"]["padj05"]
SPLIT_RHO = C3["H_A3"]["C_share"]["spearman_vs_abs_diff_f_RNA"]
MES, MES_W = C5["C5"]["neftel_MES"], C5["C5"]["wang2017_MES"]
OUT = C5["C5_outside_check"]
DIG = T1E1["T1"]["DiG"]
HC1 = DIL["tests"]["bespoke"]["HC1"]
FLIP = [f["flip_distance_beta"] for f in MGMT["flip"] if f["array"] in ("IL66B", "IL70B", "IL71B")]
STP27 = [f["stp27_prob"] for f in MGMT["flip"]]
PWR = POWER["targets"]
PRIM, REC = ["IL67B", "IL68B", "IL69B"], ["IL66B", "NL70B", "NL71B"]
CY = {k: float(v["chrY_per_M_human_counts"]) for k, v in SAMP.items() if v["chrY_per_M_human_counts"]}
N269, N168, IL64 = (float(XENG[k]["graft_pct"]) for k in ("N269B", "N168B", "IL64B"))
assert COMP["H_A1"]["pass"] is False and abs(C_SHARE - 0.4657703477667222) < 1e-12
assert all(r["criterion_B_bespoke"] == "TRUE" for r in BENCH) and len(BENCH) == 14
assert HA2["n_DE"] == 133

DESIGN = ("Orthotopic U251N glioblastoma xenograft in the athymic rat, MRI-guided LITT; bulk RNA-seq of primary and "
          "recurrent tumors, n = 3 per arm; GEO GSE338105.")
R21 = "21. Shen-Orr SS, et al. Nat Methods 2010;7:287-289."
R22 = "22. Spitzer A, et al. Nat Genet 2025;57:1168-1178."
R23 = "23. Orjefelt E, et al. J Small Anim Pract 2026;67:528-533."
R24 = "24. Steinmetz MP, et al. J Neurooncol 2001;55:167-171."
R25 = "25. Bekar A, et al. World Neurosurg 2010;73:719-721."
R26 = "26. Bady P, Delorenzi M, Hegi ME. J Mol Diagn 2016;18:350-361."
R27 = "27. Barthel FP, et al. Nature 2019;576:112-120."
R28 = "28. Browne RH. Stat Med 1995;14:1933-1940."

NEW = []          # filled by build_slides, in creation order, for arrange()


def build_slides(prs, slide):
    """Every slide v6 adds. v4's slides are already in the deck; these are appended and then re-ordered."""
    def add(key, *a, **kw):
        NEW.append(key)
        return slide(*a, **kw)

    # ---------------------------------------------------------------- promoted into the talk
    add("libraries", "panel", prs, "What each lesion is actually made of", "fig1d_lesion_accounting.png", [
        "Each library as tumor share plus the host cell types the rat reads name, scaled to the whole library.",
        "Three published methods agree on rank for all six cell types; the athymic control rejects any method that "
        "finds T cells."], core=False, seconds=15,
        ref="BRETIGEA, NeuroExpresso marker-gene profiles and CIBERSORT with a brain signature, on the rat read "
            "stream. " + DESIGN,
        note=["PROMOTED from backup in v6. 15 s. The point: 'less tumor' is not an abstraction - here is what fills "
              "the rest of the lesion, and in the two NL recurrences it is mostly myeloid cells."])

    add("less_tumour", "panel", prs, "Recurrent lesions hold less tumor - on both platforms",
        "fig1c_share_labellings_p1.png", [
            f"DNA: recurrent minus primary {float(HM1['DNA_bespoke']['rec_minus_prim']):.2f} "
            f"(90 % {float(HM1['DNA_bespoke']['exact90_lo']):.2f} to {float(HM1['DNA_bespoke']['exact90_hi']):.2f}), "
            f"the true labelling first of twenty - p = 0.05, the floor of a three-against-three design.",
            f"RNA agrees and ranks second of twenty. Every DNA fraction here is uncalibrated."],
        core=True, seconds=20,
        ref="Twenty balanced three-against-three labellings of the six tumors; the true labelling is filled. " + DESIGN,
        note=["PROMOTED in v6. 20 s. Say: whatever else differs, the recurrent lesions have less tumor in them, and "
              "that is the confound the next three slides deal with. The twenty-labelling floor means p = 0.05 is "
              "the smallest number this design can produce - quote it as a floor, never as significance."])

    add("exact_split", "panel", prs, "How much of the difference is simply less tumor?",
        "fig2a_shapley_de_genes.png", [
            f"A human graft in a rat host: every read is assigned by species, so each gene's fold change splits "
            f"exactly into composition, tumor-cell and host-cell terms.{sup(21)}",
            f"Composition is the largest single term and still under half: {C_SHARE:.2f} "
            f"({min(SENS):.2f}-{max(SENS):.2f} across four sensitivities).",
            [(f"Of the {HA2['n_DE']} genes a species-blind analysis calls differential, "
              f"{pct(HA2['share_composition_dominant'])} % are composition-dominant.", BLUE)]],
        core=True, seconds=28,
        ref=R21 + "   Per-sample mixing identity, then Shapley attribution; the three terms sum to the fold change to "
                  "1e-15. Ensembl Compara one-to-one orthologs. " + DESIGN,
        note=["NEW in v5, kept in v6 as the centre of the talk. 28 s.",
              "The idea in one breath: in a xenograft the mixture is measured, not estimated, so 'recurrence or less "
              "tumor?' has an exact answer per gene. The pre-registered threshold was 0.50 and it FAILED: composition "
              "is the biggest single term and still under half, and the genes a bulk analysis picks out are "
              "disproportionately the ones where cells actually changed.",
              f"If asked how it behaves: across the ten balanced splits the composition share rises with the split's "
              f"difference in tumor content (Spearman {SPLIT_RHO:.2f}) and the true split is highest of the ten. Noise "
              f"at n = 3 inflates the two cell-change terms, so {C_SHARE:.2f} is conservative."])

    add("splits", "panel", prs, "The composition term behaves as a real term should", "fig2b_splits_p1.png", [
        f"Across all ten balanced splits of the six tumors, the composition share rises with the split's difference "
        f"in tumor content (Spearman {SPLIT_RHO:.2f}).",
        "The true split is the highest of the ten, and sits inside the band predicted from the other nine."],
        core=False, seconds=15,
        ref="Ten balanced three-against-three splits; the true split is orange. " + DESIGN,
        note=["PROMOTED in v6. 15 s. This is the check that the composition term is not an artefact of the method: "
              "it tracks the thing it claims to measure."])

    add("categories", "panel", prs, "Motility and inflammation live in the bulk view only",
        "fig2c_categories_by_view.png", [
            "Cell motility and inflammatory response are enriched among the genes up in the species-blind (bulk) "
            "reading - the categories reported for this model.",
            "Neither appears in the tumor cells' own reading. Cell cycle appears in no view."],
        core=False, seconds=18,
        ref="g:Profiler on three fixed GO terms, Holm across three categories by four views. " + DESIGN,
        note=["PROMOTED in v6. 18 s. A neurosurgeon's version: the programs people report for recurrent tumors are, "
              "here, in the mixture rather than in the tumor cells."])

    add("mes", "panel", prs, "Where the mesenchymal shift at recurrence comes from",
        "fig2d_mes_human_vs_virtual_p1.png", [
            f"In the tumor cells alone the mesenchymal score FALLS at recurrence{sup(7)} - below every one of 1,000 "
            f"draws of a dilution null.",
            f"The rat-only score tracks myeloid content (Spearman {OUT['neftel_MES']['spearman_vs_TAM_BMDM']:.2f}), so "
            f"the bulk mesenchymal signal here is host tissue.{sup(22)}"],
        core=False, seconds=22,
        ref="7. Neftel C, et al. Cell 2019;178:835-849.   " + R22 + "   GSVA on TMM log-CPM; dot, human reads; arrow "
            "head, the same sample read as bulk. " + DESIGN,
        note=["NEW in v5, kept in v6. 22 s. 'Mesenchymal shift at recurrence' is the most repeated statement about "
              "recurrent glioblastoma; where tumor and host can be told apart exactly, the shift is in the host "
              f"compartment and the tumor cells move the other way. The composition-generated shift itself ranks "
              f"{MES['D_rank']} of 20, so we say 'not shown', not 'absent'."])

    add("dig", "panel", prs, "One tumor-cell program survives the adjustment", "fig5a_dig_arms_p2.png", [
        f"A published distal-invasion program scores higher in the recurrences after adjusting for tumor content: "
        f"first of twenty labellings (p = {DIG['p_one_sided']:.2f}).",
        [("One labelling of twenty, in a model where the recurrent arm is also the later-harvested arm. A hypothesis "
          "for the next study, not a finding.", BLUE)]],
        core=False, seconds=18,
        ref="Chanoch-Myers 2026 distal-invasion gene sets, frozen by hash; score ~ arm + tumor share, refitted for "
            "each of the twenty labellings. " + DESIGN,
        note=["PROMOTED in v6. 18 s. The pre-registration expected this to fail and it did not. Say the ceiling out "
              "loud: twenty labellings floor the p at 0.05."])

    add("meth", "panel", prs, "Methylation: nothing survives correction", "fig4b_dmp_counts.png", [
        "No probe reaches FDR < 0.05 in any design, on 866,238 probes at three against three.",
        f"What differs is tumor content: the true split ranks 1 of 10 without a fraction covariate and 4 of 10 with "
        f"it - in two independent pipelines and under all {len(BENCH)} published cross-species exclusion lists."],
        core=True, seconds=22,
        ref="Probes at p < 0.001 for all ten balanced splits, per pipeline and design; the true split is orange. " + DESIGN,
        note=["NEW in v5, kept in v6. 22 s. The pertinent negative: an array on the same six tumors finds nothing that "
              "survives multiple testing, and the difference that looks real before correction is the same tumor-content "
              "confound as on the RNA side - replicated in a second, independently written pipeline and under every "
              "published list of rat-cross-hybridising probes."])

    add("mgmt", "panel", prs, "MGMT status does not move", "fig4c_mgmt_stp27.png", [
        f"MGMT-STP27{sup(26)} calls every array methylated ({min(STP27):.3f}-{max(STP27):.3f}).",
        f"The shift that would change a call is {min(FLIP):.2f}-{max(FLIP):.2f} beta; the observed recurrent-minus-"
        f"primary shifts are -0.03 and +0.06.",
        "No animal received an alkylating agent: this is a stability readout, not a resistance result."],
        core=False, seconds=15,
        ref=R26 + "   Dashed line, the model's decision boundary; each arrow is the shift that would change that "
                  "recurrence's call. " + DESIGN,
        note=["PROMOTED in v6, and the slide this audience will ask about. 15 s."])

    add("amplitude", "panel", prs, "One DNA signal dilution does not explain", "fig4a_amplitude_dilution.png", [
        f"Diluting the parental culture in silico with rat DNA barely lowers its copy-number amplitude, so dilution "
        f"alone does not put the primaries below the curve.",
        f"Recurrent minus primary {HC1['rec_minus_prim']:.3f} (90 % {HC1['exact90'][0]:.3f} to "
        f"{HC1['exact90'][1]:.3f}), first of twenty labellings.",
        [("The null behind the band is unvalidated - the titration that would validate it failed - so this is "
          "reported as unexplained by dilution and nothing is claimed from it.", BLUE)]],
        core=False, seconds=18,
        ref="Amplitude = slope of an array's copy-number bins on the parental culture's; line and band, the culture "
            "diluted in silico with ten rat reference arrays. " + DESIGN,
        note=["PROMOTED in v6. 18 s. Say the caveat in the same breath as the result; it is why this is on a slide and "
              "not in the abstract."])

    # ---------------------------------------------------------------- the contralateral observation
    add("contra_identity", "panel", prs, "U251N in the opposite hemisphere", "fig_c1_N2_over_core_bins.png", [
        f"Rat 69's contralateral hemisphere: {N269:.1f} % of reads human, the donor genotype at 56 of 59 identity "
        f"probes, Y-chromosome dose 0.88, and the line's copy-number profile (r = 0.72).",
        f"Three measurements, two platforms, no shared failure mode. The other animal's contralateral hemisphere "
        f"reads {N168:.1f} %."],
        core=True, seconds=25,
        ref="Copy number on human-specific probes with rat cross-hybridising probes removed, scaled to the parental "
            "culture; 31-bin running medians. One animal.",
        note=["KEPT IN THE TALK at Greg's instruction. 25 s.",
              "Establish the identity first and slowly: this is not a statistical whisper, it is the cell line's own "
              "genotype, sex chromosome and copy-number profile in tissue from the other side of the brain, on an "
              "array handled separately from the libraries.",
              "Do not say anything yet about how it got there - that is the next slide but one."])

    add("contra_null", "panel", prs, "It is not simply the core, sampled thinner", "fig_a_nullA_p1.png", [
        "The null: that animal's own tumor subsampled to the contralateral sample's depth, plus the measured "
        "assignment floor, 1,000 draws.",
        "Whole-profile distance 0.153 against a null 99.5th percentile of 0.028, and larger than another animal's "
        "tumor gives against the same reference."],
        core=False, seconds=18,
        ref="Dashed and dotted lines: the same construction from two other animals' tumors. One animal.",
        note=["PROMOTED in v6. 18 s. Two of three pre-registered statistics reject; the third, the only one with a "
              "direction named in advance, does not (P = 0.56) and moves the other way. Say that if asked."])

    add("two_readings", "statement", prs, "Two readings, and what a surgeon can take from it", lines=[
        "Either U251N cells reached the opposite hemisphere in life, or core tissue was carried across on the blade "
        "when the frozen slice was cut. This design cannot separate them, and I am not going to rank them.",
        "What does not depend on which: tumor material is measurable at a site remote from the implant, by two "
        "platforms, at a level routine histology of that piece would miss.",
        [(f"For scale, someone else's measurement: malignant cells were recovered from gloves or instruments in 14 of "
          f"47 canine cancer operations (30 %), and in 57 % when margins were involved.{sup(23)}", BLUE)]],
        core=True, seconds=25,
        bar="If it was carried, then one pass of an instrument moved tumor centimetres - which is a question this "
            "method can answer, and the next study answers it with a naive hemisphere cut straight after a "
            "tumor-bearing one.",
        ref=R23 + "   " + R24 + "   " + R25 + "   Those are other groups' measurements of instrument contamination, "
            "not ours. Both contralateral samples come from animals that were never ablated.",
        note=["NEW in v6, and the slide to be careful on. 25 s.",
              "Greg, 2026-09-27: do NOT say definitively that instrument contamination is what we measured. It is not "
              "established. Say the two readings, say we cannot separate them, and put the instrument sentence in the "
              "conditional - the bar on this slide is phrased that way on purpose.",
              "What is unconditionally true: the material is there, three independent measurements on two platforms "
              "say it is U251N, and its transcriptome is not the core diluted.",
              "What is not: which reading is right; any rate (one animal of two sampled, exact interval 1.3 % to "
              "98.7 %); anything about LITT (both contralateral samples come from animals that were never ablated); "
              "anything about human resections.",
              "The record does not say whether a fresh blade was used between pieces. If asked, say that plainly.",
              "The fix is in the next study: a naive hemisphere cut with the same blade and slide immediately after a "
              "tumor-bearing one, read out by species assignment. 29 such pieces with no positive bound the "
              "false-positive rate at 10 %."])

    add("scaled", "panel", prs, "What the next experiment costs, decided in advance",
        "fig6_rats_needed_p3.png", [
            f"Targets were fixed and hashed before any pilot number entered a calculation.{sup(28)}",
            f"The binding constraint is the per-animal endpoint: {PWR['T2b']['n_per_arm_needed']} rats per arm to move "
            f"invasion frequency 0.20 to 0.50, {PWR['T1']['n_per_arm_needed']} for a 15-point difference in tumor "
            f"content, {PWR['T2c']['n_controls_needed']} control pieces to bound transfer at dissection at 10 %.",
            [(f"The clone-level test - are the cells that grew back the ones that had already left? - needs "
              f"{PWR['T8']['n_animals_needed']['k0']} to {PWR['T8']['n_animals_needed_bound_0.80']['k0']} barcoded "
              f"animals, because the unit becomes the clone.{sup(27)}", BLUE)]],
        core=False, seconds=25,
        ref=R27 + "   " + R28 + "   Power against animals; dashed line, the pre-stated target; filled mark, the "
            "conservative sizing; open marks, the two other contrasts. Panels for the genome-wide targets are pending.",
        note=["NEW in v5, one panel in v6. 25 s. The point is not 'we need more animals': every number was fixed "
              "before the pilot was allowed to speak, in the same form as the power figure in my imaging work."])

    # ---------------------------------------------------------------- new backups, one panel each
    backups = [
        ("b_reads", "Backup. What each library is made of", "fig1a_read_composition.png", [
            "Every read pair assigned by species before anything is counted.",
            f"The three floor libraries: IL64B, an implanted hemisphere with no tumor in the piece ({IL64:.1f} % "
            f"human), and the two contralateral hemispheres ({N168:.1f} % and {N269:.1f} %).",
            f"Human Y-linked transcripts per million human counts: {min(CY[k] for k in PRIM + REC):.0f}-"
            f"{max(CY[k] for k in PRIM + REC):.0f} in the six tumors, {CY['IL64B']:.1f} in the tumor-free library."],
         "xengsort (k = 25) against GRCh38 and mRatBN7.2; shares of all classified read pairs."),
        ("b_fraction", "Backup. Tumor content, three estimators", "fig1b_rna_vs_dna_p1.png", [
            "DNA fraction against RNA human share; the registered estimator is uncalibrated, so no cell percentage is "
            "quoted from either platform.",
            f"The order agrees: Spearman {float(COMP['M3_RNA_side']['spearman_six_tumours']):.2f} over the six tumors."],
         "Bars: the estimator's own interval. Dashed: equality, for reference only."),
        ("b_decount", "Backup. Differentially expressed genes across the ten splits", "fig2b_splits_p2.png", [
            "The true split gives 133 genes; the ten splits range from 0 to 303.",
            "DE count rises with the split's difference in tumor content, as the composition share does."],
         "nf-core/differentialabundance (DESeq2) on the species-blind matrix, identical parameters for all ten splits."),
        ("b_mes_wang", "Backup. The second mesenchymal signature", "fig2d_mes_human_vs_virtual_p2.png", [
            f"Wang 2017 mesenchymal: the human-only fall is below "
            f"{pct(1 - MES_W['dilution_null']['share_null_le_obs'], 1)} % of 1,000 dilution-null draws.",
            f"The composition-generated difference ranks {MES_W['D_rank']} of 20."],
         "GSVA on TMM log-CPM; dot, human reads; arrow head, the same sample read as bulk."),
        ("b_dig_scores", "Backup. The invasion score against tumor content", "fig5a_dig_arms_p1.png", [
            "Each tumor's distal-invasion score against its tumor share; the arm difference is what the model fits.",
            "Two scoring packages agree on the rank."],
         "Chanoch-Myers 2026 gene sets, frozen by hash; singscore on floor-corrected human log-CPM."),
        ("b_floor_gate", "Backup. The gate that stopped the contralateral invasion test", "fig5e_dig_floor_gate.png", [
            "The pre-registered case test asked whether the contralateral human cells sit further along the invasion "
            "axis than their own core.",
            "Its gate failed first: the assignment floor scores as high as the contralateral sample, so the contrast "
            "would measure mis-assignment. The test was abandoned and is reported as abandoned."],
         "Same score as the arm test, on every human-bearing profile plus the two floor libraries."),
        ("b_detection", "Backup. The smallest shared change this design could see", "fig5d_s1_detection.png", [
            f"Recurrences do not share derived arm-level copy-number changes the primaries lack (statistic "
            f"{T4['T_T4_arms']:.3f}, rank {T4['rank_arms']} of 20) - the pre-stated expected null.",
            f"A single-copy loss carried by fewer than {SMIN[('-1', '2')]} of recurrent cells would not have been "
            f"flagged on 80 % of arms."],
         "The registered rule applied unchanged to simulated shared changes added to the observed recurrent bins."),
        ("b_fcpg", "Backup. Fluctuating CpGs: uninformative, and why", "fig5c_fcpg_w.png", [
            "The homogeneity readout was pre-registered with two gates. The second failed: the parental culture is "
            "already multimodal at the published fluctuating-CpG sets.",
            "So no homogeneity claim is made anywhere in this work."],
         "Gabbutt 2022 and 2025 lists, frozen by hash; ten matched control CpGs per site; the band is the dilution model."),
        ("b_xval_dosage", "Backup. Array against RNA: expression follows copy number", "fig_xval_dosage_p1.png", [
            "Where the RNA is tumor-rich, expression follows copy number: slope 0.54 (p = 7e-06).",
            "The CDKN2A/B deletion every array calls is silent in every library, as it must be."],
         "Arm-level copy number against arm-median expression shift; nulls matched on expression decile."),
        ("b_xval_n2", "Backup. Why the contralateral copy-number losses stay unresolved", "fig_xval_dosage_p3.png", [
            "In the contralateral sample, arm-level expression does not follow its copy number anywhere in the genome "
            "(slope 0.17, p = 0.45), against 0.97 for its own core.",
            "So its deeper losses stay array-only: the assay has no demonstrated dosage sensitivity there."],
         "Arm copy number against arm-median expression shift, per sample."),
        ("b_xval_geno", "Backup. Array against RNA: genotype", "fig_xval_genotype_p1.png", [
            "At expressed CpG-SNP positions the RNA genotype agrees with the array in 527 of 543 pairs (97.1 %).",
            "The sixteen disagreements sit at five common germline sites carried by every sample: line-wide, not "
            "subclonal."],
         "bcftools pile-up at the array's probe positions; genotype cut-offs fixed in the pre-registration."),
        ("b_xval_meth", "Backup. Array against RNA: promoter methylation", "fig_xval_meth_p1.png", [
            "Within every sample, promoter methylation and expression anticorrelate (Spearman -0.31 to -0.36).",
            "Across the recurrence contrast they do not move together (-0.003), and the ribosomal-protein fall has no "
            "promoter or copy-number support."],
         "Promoter beta against log-CPM per sample; 10,000 gene-pairing shuffles as the null."),
        ("b_curve", "Backup. The calibration that failed", "fig_curve_mouse_validation_p2.png", [
            f"Mixing two pure arrays in silico does not reproduce real DNA mixtures: mean absolute error "
            f"{CURVE['step1']['metrics']['primary']['MAE']:.3f} against a pre-set bar of 0.05, every real mixture "
            f"reading high.",
            [("So every tumor-fraction value here is uncalibrated and every dilution null is labelled unvalidated. No "
              "published human:rat titration on a human methylation array exists.", BLUE)]],
         "GSE310817 human:mouse titrations (EPIC v1) as the validation set; points, real mixtures; line, in silico."),
        ("b_similarity", "Backup. Does the contralateral population look like the recurrences?",
         "fig_b_similarity_p2.png", [
            "No: the true arm labelling ranks fifteenth of twenty.",
            "A variant-sharing test looked strongly positive until the variants a tumor-free library also carries were "
            "removed; then 0.830 against 0.828. Host reads in the graft stream manufacture that signal."],
         "All libraries subsampled to one depth before the comparison; twenty balanced labellings."),
        ("b_arms", "Backup. The four arms named before the analysis", "fig_c1_events.png", [
            "None is a contralateral-only copy-number event: the animal's own core carries all four, and so does the "
            "parental culture.",
            "Along the genome the contralateral profile is the core's profile with the same breakpoints."],
         "Colour: amplitude-scaled log2 ratio; a dot marks a value beyond the noise bound measured on thirteen normals."),
        ("b_states", "Backup. Malignant-state composition of the human reads", "fig_c_states.png", [
            "The contralateral sample has fewer differentiated-like cells than its own core (0.037 against 0.139).",
            [("Weighted lightly: another animal's tumor differs by the same distance, and the method assigns 36 % of a "
              "pure human culture to cell types a human-only read stream cannot contain.", BLUE)]],
         "CIBERSORT with a published single-cell reference; fit correlation printed above each bar."),
        ("b_t1_power", "Backup. Sizing the tumour-content endpoint", "fig6_rats_needed_p1.png", [
            f"{PWR['T1']['n_per_arm_needed']} rats per arm for a 15-point difference in tumor content at 80 % power, "
            f"{PWR['T1']['n_per_arm_needed_holm']} under Holm; the pilot's own spread sets the variance.",
            "The dotted curve is the optimistic reading at the pilot's variance; the study is sized on the pessimistic one."],
         R28 + "   Power against rats per arm at effects of 5 to 25 points."),
        ("b_t8_power", "Backup. Sizing the clone-level test", "fig6_rats_needed_p9.png", [
            f"Wilson lower bound on the share of barcoded animals rejecting in the pre-stated direction: "
            f"{PWR['T8']['n_animals_needed']['k0']} animals with no failure reach 0.70, "
            f"{PWR['T8']['n_animals_needed_bound_0.80']['k0']} reach 0.80.",
            "The unit of inference is the clone, which is how single digits suffice."],
         "The SchwingNet Figure 5a construction: a Wilson lower bound against the number of units."),
    ]
    for key, title, image, lines, ref in backups:
        add(key, "panel", prs, title, image, lines, core=False, ref=ref + " " + DESIGN,
            note=["NEW in v6. Backup: show only if asked."])


def arrange(prs):
    """v4's slides and the new ones into the v6 running order."""
    n_v4 = 21
    idx = {k: n_v4 + i for i, k in enumerate(NEW)}
    main = [0, 1, 2, 3,                                   # title, margin, design, model
            idx["libraries"], idx["less_tumour"],
            4, 5, 6, 7,                                   # PCA+volcano, GSEA, leave-one-out, the confound
            idx["exact_split"], idx["splits"], idx["categories"], idx["mes"],
            8,                                            # subtypes
            idx["dig"],
            idx["meth"], idx["mgmt"], idx["amplitude"],
            idx["contra_identity"], idx["contra_null"], idx["two_readings"],
            9, 10,                                        # ciclopirox, why ciclopirox
            11,                                           # what this means
            idx["scaled"],
            12, 13]                                       # limitations, acknowledgments
    v4_backups = list(range(14, 21))
    new_backups = [idx[k] for k in NEW if k.startswith("b_")]
    order = main + v4_backups + new_backups
    missing = sorted(set(range(len(prs.slides))) - set(order))
    if missing:
        raise SystemExit(f"arrange: slides not placed: {missing}")

    lst = prs.slides._sldIdLst
    ids = list(lst)
    for e in ids:
        lst.remove(e)
    for i in order:
        lst.append(ids[i])
    print(f"arranged: {len(main)} main ({sum(1 for _ in main)} slides), {len(v4_backups) + len(new_backups)} backups")
