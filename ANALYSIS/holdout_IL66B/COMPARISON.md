# IL66B held out: the full pipeline against the published run

control = the published six-tumour inputs re-run today through the same R and Python stages; holdout = IL66B removed from the sample sheet and every stage re-run from DESeq2.


## 1. Differential expression (DESeq2, FDR < 0.05, |log2FC| > 1)

| | published | held out |
|---|---|---|
| genes | 35 (23 up) | 140 (83 up) |
| shared | 30 (Jaccard 0.21) | |
| Spearman of log2FC, all genes | 0.911 | |


Lost without IL66B (5): ENSG00000289697, GPC3, IGFBP3, NAALADL2, PRSS12

Gained (110): AC084219.4, ACSBG1, ADGRG6, ALDH1A3, ANK3, ANO4, APCDD1, AQP1, ARAP2, ARHGAP19, B4GALNT4, C15orf59, C3, CHAT, CRB2, CSF2RA, CYB5R2, DAAM2, DIRC3-AS1, EDNRB, EGFL8, EGR2, ENO1P4, ENSG00000285920, ENSG00000288709, ENSG00000290842, ERICH3, F3, F8A3, FAM71F1, FAM84A, FLG, FOXD3-AS1, GALNT16, GLDC, HTR1D, IL4R, KLF5, KNL1, LINC01137, LINC01468, LONRF2, LRRC8C, LTF, MAP3K7CL, MFAP5, MIR3648-2, MUC16, MYB, NBEAL2, NPSR1-AS1, NRCAM, NTRK2, OXTR, PAK3, PARP8, PIK3CD, PKDCC, PLEKHH1, PPM1E


## 2. Gene-set enrichment (the pipeline's Broad GSEA, seed 1234)

| set | NES published | q published | NES held out | p held out | q held out |
|---|---|---|---|---|---|
| KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION | -1.99 | 0.022 | -1.95 | 0.000 | 0.309 |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION | -1.96 | 0.158 | -1.87 | 0.005 | 1.000 |
| REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY | -1.95 | 0.146 | -1.91 | 0.000 | 0.785 |
| KEGG_RIBOSOME | -1.94 | 0.204 | -1.93 | 0.000 | 0.586 |
| REACTOME_SELENOAMINO_ACID_METABOLISM | -1.94 | 0.169 | -1.96 | 0.007 | 0.722 |
| REACTOME_CELLULAR_RESPONSE_TO_STARVATION | -1.93 | 0.247 | -1.89 | 0.007 | 1.000 |

Sets at q < 0.05 / q < 0.25 (down): published 1 / 6; held out 0 / 0. Up side at q < 0.25: 0 and 0.


## 3. DSigDB drug screen: the 100 most opposing drug gene sets (R stage)

Overlap of the 100 compounds (Jaccard): control vs published 1.00; held out vs published 0.41; held out vs control 0.41.

| arm | ciclopirox NES | FDR | rank (R integrated score) |
|---|---|---|---|
| published | -2.261 | 1.7e-05 | 1 |
| control | -2.261 | 1.7e-05 | 1 |
| holdout | -1.638 | 0.04 | 2 |

Entered with IL66B held out (38): 1-methylphenanthrene, 2,2',3,4,4',5,5'-heptachlorobiphenyl, 2e,4e,6e,8e, 5-aminosalicylic acid, alitretinoin, anthracene, benzyl butyl phthalate, caffeic acid, calcium, camptothecin, chlorhexidine diacetate, chlorpyrifos, creatine, cytarabine, endosulfan, fenvalerate, flecainide, fluoranthene, gabexate, gadolinium, gentamicin, ionomycin, leflunomide, lynestrenol, menadione bisulfite, methylergometrine, pentamidine, perhexiline, pyrogallol, pyrvinium, raloxifene, rofecoxib, s-1,2-dichlorovinyl-n-acetylcysteine, sodium dichromate, stannic fluoride, tetracycline, vinblastine, zeranol

Left (39): 0297417-0002b, 1,6,7,8,9,11a,12,13,14,14a-decahydro-1,13-dihydroxy-6-methyl-4h-cyclopent [f]oxacyclotridecin-4-one, 7-aminocephalosporanic acid, 7-diethylamino-4-methylcoumarin, amiodarone, astemizole, axitinib, cefotiam, cumene hydroperoxide, cyperquat chloride, dacarbazine, econazole, escitalopram oxalate, estropipate, ethaverine, ferric ammonium, ferrous, gnf-pf-3832, heptachlor, ifosfamide, lysergide, mefloquine, nicotinic acid, nilutamide, nitrofural, osajin, oxygen, paroxetine, pentetrazol, pomiferin, primidone, propantheline bromide, salbutamol, succinylsulfathiazole, terfenadine, tolbutamide, trimipramine, vandetanib, wortmannin


## 4. Final ranking (clinical phase, ADMET-AI, BOILED-Egg)

| arm | clinically available | clear both barrier models | ciclopirox rank (weighted / unweighted) | ciclopirox score |
|---|---|---|---|---|
| published | 54 | 15 | 1 / 1 | 3.14 |
| control | 54 | 13 | 1 / 1 | 3.27 |
| holdout | 49 | 11 | 3 / 14 | 2.02 |

Top 20 overlap (Jaccard): control vs published 0.82; held out vs control 0.33.

| rank | published | control | held out |
|---|---|---|---|
| 1 | ciclopirox | ciclopirox | flecainide |
| 2 | pentetrazol | pentetrazol | magnesium |
| 3 | primidone | nilutamide | ciclopirox |
| 4 | nilutamide | primidone | progesterone |
| 5 | pyrantel | pyrantel | thioridazine |
| 6 | ifosfamide | ifosfamide | diazepam |
| 7 | magnesium | magnesium | perhexiline |
| 8 | diazepam | progesterone | malathion |
| 9 | dexverapamil | diazepam | paricalcitol |
| 10 | malathion | trimipramine | s-1,2-dichlorovinyl-n-acetylcysteine |
| 11 | progesterone | thioridazine | l-citrulline |
| 12 | thioridazine | malathion | ozone |
| 13 | oxygen | oxygen | vincristine |
| 14 | astemizole | nitrofural | leflunomide |
| 15 | nitrofural | l-citrulline | pyrantel |
| 16 | ganciclovir | ozone | rofecoxib |
| 17 | paroxetine | dacarbazine | methylergometrine |
| 18 | ozone | ganciclovir | dexverapamil |
| 19 | trimipramine | astemizole | vinblastine |
| 20 | amiodarone | dexverapamil | alitretinoin |

## 5. Subtype signatures (GSVA; change recurrent minus primary, Welch p)

| signature | published | control | held out |
|---|---|---|---|
| Garofano_MTC | -0.58 (p 0.142) | -0.58 (p 0.142) | -0.73 (p 0.127) |
| Neftel_AC | -0.53 (p 0.007) | -0.53 (p 0.007) | -0.50 (p 0.066) |
| Neftel_MES1 | -0.48 (p 0.089) | -0.48 (p 0.089) | -0.32 (p 0.329) |
| Neftel_MES2 | -0.36 (p 0.345) | -0.36 (p 0.345) | -0.12 (p 0.788) |
| Neftel_NPC2 | -0.15 (p 0.465) | -0.15 (p 0.465) | -0.34 (p 0.087) |
| Neftel_NPC1 | +0.04 (p 0.876) | +0.04 (p 0.876) | -0.12 (p 0.630) |
| Neftel_OPC | +0.05 (p 0.511) | +0.05 (p 0.511) | +0.07 (p 0.694) |
| Garofano_GPM | +0.11 (p 0.698) | +0.11 (p 0.698) | +0.14 (p 0.707) |
| Garofano_PPR | +0.31 (p 0.419) | +0.31 (p 0.419) | +0.58 (p 0.113) |
| Garofano_NEU | +0.34 (p 0.267) | +0.34 (p 0.267) | +0.29 (p 0.489) |

## 6. STRING network among the DE genes

| arm | nodes | edges | hubs (degree ≥ 3) |
|---|---|---|---|
| published | 11 | 30 | BGN, COL17A1, COL1A1, GPC3, HMCN1, IGFBP3, LRP1, PCLO, PRELP, RYR2 |
| control | 11 | 30 | BGN, COL17A1, COL1A1, GPC3, HMCN1, IGFBP3, LRP1, PCLO, PRELP, RYR2 |
| holdout | 49 | 126 | ADGRG6, ANK3, BGN, C3, CALB1, CHAT, COL17A1, COL1A1, CSF2RA, EGR2, F3, FLG, HMCN1, HTR2A, IL4R, LTF, MFAP5, MUC16, NRCAM, NRG1, NTRK2, PCLO, PIK3AP1, PIK3CD, PRELP, ROBO2, RYR2, S100A1, S100B, SCN3A, SEMA3F, SERPINA3, SERPINA5, TACR1 |

Not rerun: the prior-art audit (a manual literature review of the top 20; compounds entering the held-out top 20 would need it) and DepMap (a property of the cell line, not of the samples).

