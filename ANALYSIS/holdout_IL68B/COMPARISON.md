# IL68B held out: the full pipeline against the published run

control = the published six-tumour inputs re-run today through the same R and Python stages; holdout = IL68B removed from the sample sheet and every stage re-run from DESeq2.


## 1. Differential expression (DESeq2, FDR < 0.05, |log2FC| > 1)

| | published | held out |
|---|---|---|
| genes | 35 (23 up) | 43 (33 up) |
| shared | 21 (Jaccard 0.37) | |
| Spearman of log2FC, all genes | 0.837 | |


Lost without IL68B (14): ALG10B, CALB1, GDAP1, HTR2A, KCCAT333, MOB3B, MRAP2, NAALADL2, NRG1, PCLO, PIK3AP1, PRSS12, RASL11A, SORL1

Gained (22): AC073316.1, ADGRG6, AP1M2, ARHGEF5, CH507-513H4.6, CYB5R2, ENSG00000287515, ERMN, FAM71F1, GLDC, IL7R, MIR568, MYEOV, NR1H4, RNU6-11P, RNU6-5P, RNVU1-18, RNVU1-7, SPAG17, TSPY8, XACT, ZNF717


## 2. Gene-set enrichment (the pipeline's Broad GSEA, seed 1234)

| set | NES published | q published | NES held out | p held out | q held out |
|---|---|---|---|---|---|
| KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION | -1.99 | 0.022 | -1.93 | 0.000 | 0.027 |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | -1.96 | 0.158 | -1.96 | 0.000 | 0.009 |
| REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY | -1.95 | 0.146 | -1.91 | 0.000 | 0.034 |
| KEGG_RIBOSOME | -1.94 | 0.204 | -1.95 | 0.000 | 0.010 |
| REACTOME_SELENOAMINO_ACID_METABOLISM | -1.94 | 0.169 | -1.93 | 0.000 | 0.021 |
| REACTOME_CELLULAR_RESPONSE_TO_STARVATION | -1.93 | 0.247 | -1.90 | 0.000 | 0.065 |

Sets at q < 0.05 / q < 0.25 (down): published 1 / 6; held out 12 / 29. Up side at q < 0.25: 0 and 0.


## 3. DSigDB drug screen: the 100 most opposing drug gene sets (R stage)

Overlap of the 100 compounds (Jaccard): control vs published 1.00; held out vs published 0.41; held out vs control 0.41.

| arm | ciclopirox NES | FDR | rank (R integrated score) |
|---|---|---|---|
| published | -2.261 | 1.7e-05 | 1 |
| control | -2.261 | 1.7e-05 | 1 |
| holdout | -2.146 | 0.00061 | 1 |

Entered with IL68B held out (41): 2e,4e,6e,8e, 4-oxoretinol, alsterpaullone, amprolium, benzyl butyl phthalate, bisindolylmaleimide ix, chloramphenicol, chlorophenothane, curcumin ii, cytarabine, d-penicillamine, docetaxel, endosulfan, ethanol, ethylene dimethacrylate, fenvalerate, folic acid, gabexate, gemfibrozil, geranylgeranyl pyrophosphate, gly-his-lys, healon, hypochlorous acid, ibmx, iron, letrozole, maleic hydrazide, olanzapine, pamidronate, ribavirin, rosavin, s-1,2-dichlorovinyl-n-acetylcysteine, sodium dichromate, sodium hydroxide, tamoxifen, tetrachloroethylene, trolox, vanadium, vanadium pentoxide, vinblastine

Left (37): 1-naphthyl isothiocyanate, 3-methoxycatechol, 6-thioguanine, 7-diethylamino-4-methylcoumarin, amphotericin b, astemizole, atovaquone, axitinib, bortezomib, carfilzomib, chloranil, cp-690334-01, cumene hydroperoxide, cyperquat chloride, dacarbazine, econazole, escitalopram oxalate, ethaverine, ferric ammonium, gnf-pf-3832, gossypol, lysergide, malathion, mefloquine, metformin, oxygen, paroxetine, propantheline bromide, salbutamol, taurine, terfenadine, thioridazine, tolbutamide, toluidine blue o, tributyltin, trimipramine, wortmannin


## 4. Final ranking (clinical phase, ADMET-AI, BOILED-Egg)

| arm | clinically available | clear both barrier models | ciclopirox rank (weighted / unweighted) | ciclopirox score |
|---|---|---|---|---|
| published | 54 | 15 | 1 / 1 | 3.14 |
| control | 54 | 13 | 1 / 1 | 3.27 |
| holdout | 48 | 9 | 1 / 1 | 3.02 |

Top 20 overlap (Jaccard): control vs published 0.82; held out vs control 0.43.

| rank | published | control | held out |
|---|---|---|---|
| 1 | ciclopirox | ciclopirox | ciclopirox |
| 2 | pentetrazol | pentetrazol | ozone |
| 3 | primidone | nilutamide | magnesium |
| 4 | nilutamide | primidone | d-penicillamine |
| 5 | pyrantel | pyrantel | pentetrazol |
| 6 | ifosfamide | ifosfamide | diazepam |
| 7 | magnesium | magnesium | pyrantel |
| 8 | diazepam | progesterone | l-citrulline |
| 9 | dexverapamil | diazepam | paricalcitol |
| 10 | malathion | trimipramine | hypochlorous acid |
| 11 | progesterone | thioridazine | primidone |
| 12 | thioridazine | malathion | progesterone |
| 13 | oxygen | oxygen | nilutamide |
| 14 | astemizole | nitrofural | ethanol |
| 15 | nitrofural | l-citrulline | amprolium |
| 16 | ganciclovir | ozone | ifosfamide |
| 17 | paroxetine | dacarbazine | nitrofural |
| 18 | ozone | ganciclovir | zinc sulfide |
| 19 | trimipramine | astemizole | chlorophenothane |
| 20 | amiodarone | dexverapamil | vincristine |

## 5. Subtype signatures (GSVA; change recurrent minus primary, Welch p)

| signature | published | control | held out |
|---|---|---|---|
| Garofano_MTC | -0.58 (p 0.142) | -0.58 (p 0.142) | -0.43 (p 0.262) |
| Neftel_AC | -0.53 (p 0.007) | -0.53 (p 0.007) | -0.57 (p 0.044) |
| Neftel_MES1 | -0.48 (p 0.089) | -0.48 (p 0.089) | -0.38 (p 0.204) |
| Neftel_MES2 | -0.36 (p 0.345) | -0.36 (p 0.345) | -0.30 (p 0.524) |
| Neftel_NPC2 | -0.15 (p 0.465) | -0.15 (p 0.465) | -0.06 (p 0.830) |
| Neftel_NPC1 | +0.04 (p 0.876) | +0.04 (p 0.876) | +0.11 (p 0.719) |
| Neftel_OPC | +0.05 (p 0.511) | +0.05 (p 0.511) | +0.05 (p 0.700) |
| Garofano_GPM | +0.11 (p 0.698) | +0.11 (p 0.698) | -0.04 (p 0.844) |
| Garofano_PPR | +0.31 (p 0.419) | +0.31 (p 0.419) | +0.40 (p 0.365) |
| Garofano_NEU | +0.34 (p 0.267) | +0.34 (p 0.267) | +0.29 (p 0.450) |

## 6. STRING network among the DE genes

| arm | nodes | edges | hubs (degree ≥ 3) |
|---|---|---|---|
| published | 11 | 30 | BGN, COL17A1, COL1A1, GPC3, HMCN1, IGFBP3, LRP1, PCLO, PRELP, RYR2 |
| control | 11 | 30 | BGN, COL17A1, COL1A1, GPC3, HMCN1, IGFBP3, LRP1, PCLO, PRELP, RYR2 |
| holdout | 11 | 26 | BGN, COL17A1, COL1A1, GPC3, IGFBP3, LRP1, PRELP |

Not rerun: the prior-art audit (a manual literature review of the top 20; compounds entering the held-out top 20 would need it) and DepMap (a property of the cell line, not of the samples).

