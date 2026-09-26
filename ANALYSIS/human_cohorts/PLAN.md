# Human replication plan for the U251 LITT xenograft signals

2026-09-26. Inputs: `cohort_inventory.tsv` (27 cohorts), `raw/MANIFEST.tsv` (100 files, 2.95 GB, md5/size-checked), `../power/`, `../graft_relation/`. Two adversarial passes re-checked every local number; where they disagreed with the first draft, the corrected value is used here. Nothing is committed. `raw/` is excluded through `.git/info/exclude`.

## 0. What is being replicated (corrected signals)

| Signal | What the data show | Caveat that changes the plan |
|---|---|---|
| "Translation" fall | Six MSigDB sets, but they are one ~80-gene cytosolic ribosomal-protein (RP) block. KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION is 81/81 RPL/RPS/FAU/UBA52, with no initiation factors. GSEA NES −1.99, q 0.022 (seed 1234); q runs from 0.022 to 0.31 over five seeds. DESeq2 mean log2FC is −0.40. | **IL68B drives the effect size, not the NES.** Without IL68B the NES is −1.92 to −2.00 and q falls to 0.001–0.032, while the DESeq2 set mean shrinks from −0.40 to −0.10 (`graft_relation/SUMMARY.md` C7). Significance breaks when IL67B, IL69B or IL66B is dropped (q 0.08–1.0; `gsea_leave_one_out/loo_sets.tsv`). The RP share of raw counts depends on library composition. |
| Neftel AC-like fall | GSVA 0.261 ± 0.120 → −0.269 ± 0.139 (Welch p 0.0079; `subtypes/subtype_rerun_scores_per_sample.csv`). | All p values are parametric. With 3 vs 3 the exact permutation floor is **p = 0.10**. CIBERSORT AC goes the **other** way (0.221 → 0.271). |
| MES fall | GSVA MES1 −0.48 (p 0.087). CIBERSORT MES 0.51 → 0.39 (p 0.16). | **Drop CIBERSORT MES as an endpoint.** The pure-rat control N168B (graft 0.55 %) returns MES 0.451 at CIBERSORT p 0.045, and IL64B (graft 0.33 %) returns NPC 0.769. Fit R is 0.10–0.16 for the controls against 0.19–0.26 for the tumours (`cibersort/results/v104/fractions_v104_neftel4_confident.tsv`). |
| Ciclopirox | Rank 1 in the published run and with IL68B held out (`holdout_IL68B/COMPARISON.md`). | Drops to rank 3 without IL66B. It has no per-sample signal: DSigDB ciclopirox GSVA group effect +0.03, p 0.89. |

## 1. Your question: does the xengsort graft % relate to the results?

**Yes, through the groups. Rat reads leaking into the human stream do not produce the results.** Recurrences carry less human tissue: graft % is 29.44 / 42.67 / 33.05 against 42.52 / 45.38 / 64.40 (Welch p 0.14, d −1.62; `metadata_full.csv`). Six animals cannot tell "recurrence" apart from "less tumour in the piece sequenced".

Current `graft_relation/results/score_correlations.tsv` (13:19 run, GSVA on DESeq2 vst):

| Signal | vs graft % | Group effect → graft-adjusted | Reading |
|---|---|---|---|
| RP block (6 sets) | r 0.45–0.53 (p ≥ 0.27); Spearman 0.77–0.89 | GSVA −0.76 → −0.83 (p 0.06 → 0.18). DESeq2 set mean −0.40 → −0.57 (`translation_sets_lfc_by_model.tsv`) | Adjustment makes the fall **larger**; only precision is lost. The earlier figures (r 0.76–0.84; −0.85 → −0.49) came from un-normalised log2 TPM, where 7SK/7SL/Y RNA and rRNA take 40–65 % of TPM. The script now uses vst. |
| Neftel AC | r 0.70 (p 0.12) | −0.573 (p 0.0021) → −0.555 (p 0.024) | Survives adjustment. It also holds in the one graft-matched pair, NL70B (42.67 %) vs IL67B (42.52 %): −0.61 (SUMMARY C5). |
| Neftel MES1 | r 0.77 (p 0.076); r with both-share m −0.82 (p 0.048); within-group r 0.55 | −0.478 (p 0.087) → −0.264 (p 0.44) | Tracks sample composition as much as group. Do not present it as a group result. |
| DE genes (FDR < 0.05) | graft term alone: 9 genes | 102 → 42 | The drop is lost precision (df 4 → 3). All 102 genes keep their sign (SUMMARY C4). |
| Ciclopirox | r 0.22 | +0.03 (p 0.89) → +0.19 (p 0.59) | No per-sample signal either way. |

The size of any rat leak:
- **Estimated rat share of the human stream:** 1.9 % (IL69B) to 8.0 % (IL66B). This assumes the pure-rat controls' both/host ratio (5.3 %) carries over to tumour libraries. With N269B included the ratio is 6.5 % and the range becomes 2.3–9.8 %.
- **Where the leak lands:** 70–71.5 % on RN7SL1. RP genes are 0.06–0.24 % of the pure-rat human stream, against about 3–5.5 % in the tumours.
- **Worst-case effect on RP share:** 0.048 log2 out of an observed −0.40, so at most about 12 %. The mixture bound is 0.009 log2, and 0.0008 once N269B (4.9 % graft, not rat-only) is excluded.
- **RUVSeq k = 1 is not decontamination.** W_1 separates the groups completely (r −0.87 with recurrence, p 0.024) but has r 0.33 with graft. Its −0.40 → −0.11 therefore removes the group difference, not contamination.

Corrections to earlier statements:
- **IL66B** is not the most host-rich library: its host % is 59.17 against NL71B's 59.23. It is extreme on read composition instead: both-share 0.246, unique mapping 66.4 % vs 73.6–77.2 %, MapQ0 16 %, and the smallest human library (20.3 M).
- **IL68B** has the second-highest graft % (45.38), not a middling one. It is a composition outlier (7SK/7SL 10.1 %).
- RP share vs graft is monotone in rank but not linear (Spearman 0.89, exact p 0.033; Pearson 0.30). With short ncRNA, rRNA and mitochondrial reads taken out of the denominator, rho falls to 0.71 (p 0.11).

Human data face the same issue. Tumour purity falls at human recurrence: Hoogstrate 2023 (https://doi.org/10.1016/j.ccell.2023.02.019); Spitzer 2025 (https://doi.org/10.1038/s41588-025-02168-4); and Varn 2022, where oligodendrocytes rise at IDHwt recurrence (P = 5e-6) alongside the AC-like fall (https://doi.org/10.1016/j.cell.2022.04.038). Graft % is the xenograft's purity. Every human analysis below therefore runs with and without a purity covariate. Graft % may also be a mediator (LITT causing smaller regrowth), in which case adjusting for it removes real effect.

## 2. Recommended human cohorts, ranked

Counts are IDH-wildtype GBM. P = primary, R = recurrence.

| # | Cohort | IDH-wt P / R (pairs) | Access | On disk (`raw/`) | Why this rank |
|---|---|---|---|---|---|
| 1 | **GLASS, Synapse release** (Varn 2022, https://doi.org/10.1016/j.cell.2022.04.038; https://www.synapse.org/glass, syn17038081) | 128 pairs (of 168 RNA pairs, 304 patients) | Free Synapse account plus self-sign click-through, AR 9604985 (DOWNLOAD; no certification or validation needed). The terms text could not be read without a login. | **No.** Needed: syn69961520 (gene TPM), syn71730743 (transcript counts), syn69917798 (clinical_surgeries), syn69917782 (rnaseq_pairs), syn69917776 (estimate_purity) | Largest paired set, linear-scale TPM and a purity table. **Needs Greg to log in and accept.** The treatment column list is unverified without a login. |
| 2 | **GSE222515** (Dhawan 2023, https://doi.org/10.3390/cancers15030670) | 40 pairs; IDH1 wt 13, mutant 1, unknown 26 | Open | Yes: TPM, 36,679 genes × 80 | Best open paired TPM. Treatment is given at cohort level only (TMZ 34/40, bevacizumab 13/40, lomustine 13/40), so bevacizumab cases cannot be excluded. |
| 3 | **GSE174554** snRNA-seq (Wang L 2022, Nature Cancer, https://doi.org/10.1038/s43018-022-00475-x) | 30 pair IDs with both time points (81 GSMs: 33 P, 48 R) | Open (processed); raw at EGA EGAS00001004909 | Yes: RAW.tar 1.65 GB, with tumour/normal labels | Malignant-cell pseudobulk is **purity-free, like the xengsort human stream**. The per-sample malignant fraction is the direct analogue of graft %. Recurrences followed standard TMZ + RT only. |
| 4 | **GLASS figshare subset** (Ho 2024, https://doi.org/10.1186/s40478-024-01790-3; https://doi.org/10.6084/m9.figshare.25051553.v3) | 87 pairs, 21,475 genes | Open, CC BY 4.0 | Yes, 3 files, md5 match | Usable **now**, for rank-based paired GSEA/GSVA only: the values are log2(TPM+1) residuals, so no CIBERSORT. It is a subset of #1, not independent. **Trap:** the gene column is "NA" in `GLASS_87p_..._source.txt`. Take symbols from the `_caseID` file (same row order; paired differences agree to 4e-13). |
| 5 | **CGGA mRNAseq_693 + _325** (https://doi.org/10.1016/j.gpb.2020.10.005; https://doi.org/10.1038/sdata.2017.24) | 183 P / 96 R, **unpaired** (693: 109/81; 325: 74/15) | Open | Yes: 6 zips (FPKM, counts, clinical) | Largest open contrast, and independent of GLASS. No patient key, two batches. Radio_status and Chemo_status flags are available. |
| 6 | Kim KH 2024 (https://doi.org/10.1016/j.ccell.2023.12.015; PRJNA1051047) | 84 pairs by SRA sample names (the paper says 86); IDH mixed | Open FASTQ, 991 GB | No | Largest open paired raw set, but it needs grid re-alignment and per-patient IDH from the paper's tables. Samsung and SNU cases may overlap GLASS. Run only if #1–#5 disagree. |
| 7 | TCGA-GBM (Xena/GDC) | 13 R; 6 pairs (Xena) or 9 (GDC, different batch) | Open | Yes, md5 30/30 | Supportive only. Its pairs are in GLASS. G-CIMP proxy: 11 of 13 non-G-CIMP. |
| 8 | CPTAC-3 brain (GDC) | 15 GBM R of 19 cases; 16 with a primary | Open | Yes, md5 238/238 | No IDH, no treatment, study of origin unknown. |

Secondary sets (downloaded; too small, off-design or overlapping):
- GSE339484: 11 IDH1-wt pairs. Labelled "glioma", count unit unverified.
- GSE155434: 11 IDH-wt patients. Exon-level counts, overlaps #3.
- GSE139533: 13 patients, 81 multi-region samples.
- GSE62153: 15 pairs, microarray.
- GSE212067: 7 pairs, array.
- GSE4271: 23 pairs, pre-IDH era.
- GSE79671: a bevacizumab contrast, not primary vs recurrence.
- GSE254873/4: NanoString, 784 genes.
- GSE228497: primary/recurrent labels not found.

### Which of our analyses to run on each

| Analysis | #1 GLASS | #2 GSE222515 | #3 GSE174554 | #4 figshare | #5 CGGA |
|---|---|---|---|---|---|
| RP-block paired GSEA (the 81-gene set, plus the other five) on the recurrent − primary ranking | yes | yes | yes, malignant pseudobulk | yes | unpaired, batch covariate |
| GSVA: Neftel AC/MES1/MES2/NPC/OPC, Garofano | yes | yes | yes | yes, rank-valid | yes |
| CIBERSORT v1.04, Neftel signature: in human bulk use `neftel4_confident_plus_nonmalignant` (c = 7), so AC-like is not confused with astrocytes | yes (TPM) | yes (TPM) | not needed: states are called per cell | **no** (residuals) | yes (FPKM → TPM) |
| Ciclopirox/iron axis: rerun the DSigDB reversal screen on each cohort's paired ranking and record the ciclopirox rank; paired change in TFRC, FTH1, FTL, DOHH, DHPS, EIF5A | yes | yes | yes | rank screen only | unpaired |
| Purity covariate | GLASS ESTIMATE table | ESTIMATE | malignant fraction | none (already regressed) | ESTIMATE |
| CIBERSORTx group mode (per class needs k ≥ 16–20 for c = 4, or 28–35 for c = 7) | yes (128 per class) | c = 4 only (40) | not needed | no | yes (P 183; R 96) |

Replication criteria, fixed before running:
- **RP block:** the paired median change is negative, with paired Wilcoxon p < 0.05, **and** it stays negative after purity adjustment. RP share must be computed with short ncRNA, rRNA and mitochondrial reads out of the denominator.
- **AC:** the same test on Neftel AC.
- **MES:** report whichever direction appears. Here the xenograft and the human literature predict opposite directions.
- **Ciclopirox:** report its rank in the screen. No threshold is set in advance.

Order of work: #4 (today, local R 4.5.2; pandas is broken in the local Python), then #2, then #3, then #5. #1 follows once AR 9604985 is accepted, and #6 on the grid only if needed.

## 3. What the literature says about each xenograft signal at human recurrence

None of these human cohorts received LITT. They are recurrences after surgery, RT and TMZ (and bevacizumab in some), so agreement in direction is the most that can be claimed.

| Xenograft signal | Human recurrence | Verdict | Sources |
|---|---|---|---|
| RP / "translation" fall | No paired primary–recurrent study reports a fall. "ribosom" returns zero hits in the Varn 2022 and Kim 2024 full texts, and none of Wang L 2022, Spitzer 2025 or Hoogstrate 2023 (abstract only) reports one. Indirect evidence: in Kim 2024, PPR fell 14.2 → 10.5 %, MTC 23.8 → 16.8 % and NEU rose 36.2 → 49.0 % (Fig 5F). Wang L 2022 found no change in cycling cells (Fig 3b). The same rat model reports cell cycle **up** in recurrence (Nagaraja 2026, 4 vs 4). | **Unknown.** In tension with the same-model study. | https://doi.org/10.1016/j.ccell.2023.12.015 ; https://doi.org/10.1038/s43018-022-00475-x ; https://doi.org/10.3171/2026.1.JNS241395 ; https://doi.org/10.1093/neuonc/noae165.1291 |
| Neftel AC-like fall | Varn 2022, IDHwt: "a significant decrease at recurrence in the astrocyte-like neoplastic cell state" (P = 7e-3, paired t, Fig S1G). It comes with more oligodendrocytes (P = 5e-6), i.e. the same admixture confound as graft % here, and remained significant after accounting for extent of resection. The effect size is not in the PMC text. Wang Q 2017 (91 IDHwt pairs): classical retained in 51 % and "least frequently found" at recurrence. Kim 2024 single-cell (Fig 5F) and Spitzer 2025 (ED Fig 3a,b) found no significant change. | **Agrees with bulk; not reproduced in single-cell.** | https://doi.org/10.1016/j.cell.2022.04.038 (PMC9189056) ; https://doi.org/10.1016/j.ccell.2017.06.003 ; https://doi.org/10.1038/s41588-025-02168-4 |
| MES fall | Wang L 2022: MES rose per patient (Fig 3a, 31 pairs, one-sided Wilcoxon P = 0.0397), with 19 of 31 patients transitioning to MES. Hoogstrate 2023: "preferential mesenchymal progression". Wang Q 2017: MES was the most stable subtype (65 %). Spitzer 2025: no consistent state enrichment, and MES-like lost in MGMT-low patients (Fig 4f). Ho 2024: about 50 % of pairs switch subtype; non-MES → MES was the most common switch but not significantly more frequent than the reverse. Proposed drivers are RT (Bhat 2013; Halliday 2014; Minata 2019) and myeloid OSM (Hara 2021). The athymic rat has myeloid cells but no T cells. | **Disagrees** with most cohorts; consistent only with Spitzer 2025. | https://doi.org/10.1038/s43018-022-00475-x ; https://doi.org/10.1016/j.ccell.2023.02.019 ; https://doi.org/10.1186/s40478-024-01790-3 ; https://doi.org/10.1016/j.ccr.2013.08.001 ; https://doi.org/10.1073/pnas.1321014111 ; https://doi.org/10.1016/j.celrep.2019.01.076 ; https://doi.org/10.1016/j.ccell.2021.05.002 |
| Stem-like / proliferative | Varn 2022: the proliferating stem-like rise is significant only in IDHmut (P = 1e-3). In IDHwt, differentiated-like cells fell (P = 4e-3), and proneural-to-mesenchymal transitions came with fewer stem-like cells (P = 3e-4, Fig 2E). | **Not comparable.** The Varn reference calls U251 about 100 % proliferating stem-like (not re-verified). | https://doi.org/10.1016/j.cell.2022.04.038 |
| Ciclopirox / iron | No human recurrence data. Preclinical evidence: U251 in vitro and subcutaneous in nude mice (Su 2021; IC50 34–74 µM by MTT); an independent signature-reversal screen that ranked CPX in its top five (Sun 2025; IC50 0.76–3.74 µM by ATP luminescence); CPX inhibits DOHH, blocking eIF5A hypusination, a translation link (Ofek 2023). GBM stem cells depend on iron (Schonberg 2015). Deferasirox loses activity on U251 at 3 % O₂ (Legendre 2016). ClinicalTrials.gov on 2026-09-26: no glioma trial of CPX, fosciclopirox, deferoxamine, deferiprone or deferasirox. Triapine (NCT06410248, NCT06860594) and gallium maltolate (NCT04319276, NCT07515924) are in trials. Brain PK of CPX has not been measured (Weir 2019 reports plasma and urine only). In Enrichr's DSigDB copy, the five ciclopirox CMap sets hold 0–1 RP genes of 61–177, so the rank is unlikely to rest on the RP block. The pipeline's own `hs.dsigdb.v1.0.gmt` is not on this machine, so this was not checked against it. | **Unknown in humans**; preclinical support only. | https://doi.org/10.1038/s41419-021-03535-9 ; https://doi.org/10.1186/s12967-024-06046-1 ; https://doi.org/10.1002/ijc.34545 ; https://doi.org/10.1016/j.ccell.2015.09.002 ; https://doi.org/10.1186/s12885-016-2074-y ; https://doi.org/10.1124/jpet.119.257972 ; https://clinicaltrials.gov/study/NCT06410248 |
| LITT bed | After LITT the BBB is open for weeks (Leuthardt 2016: back to baseline by about 6 weeks; Bartlett 2023: by 8 weeks). Elder 2019 describes a 1.3 ± 0.3 mm granulation rim with CD68 macrophages. No study measured translation or iron programmes in tumour regrowing after LITT. | **Unknown.** | https://doi.org/10.1371/journal.pone.0148613 ; https://doi.org/10.7759/cureus.37397 ; https://doi.org/10.1186/s13000-019-0794-4 |

## 4. Sample size for more xenografts

Assumptions:
- Two-sided α 0.05, 80 % power, equal groups.
- Exact Welch noncentral t, with each group keeping its own observed SD.
- The SDs come from 3 animals each (95 % CI 0.52–6.3×), so **every n below could be off by 0.27× to 39.5×**.
- No allowance for graft failure. IL64B failed to graft; the failure rate is unknown.
- About 25 M human reads per library.

Sources: `../power/power_table.csv`, `effects.json`, `dispersion.json`, `power_rnaseq.json`, `proper_quick.json`.

| Endpoint | Observed P vs R (mean ± SD) | d | n/group at d | at 0.8 d | at d/2 | Use as endpoint? |
|---|---|---|---|---|---|---|
| GSVA Neftel AC | 0.261 ± 0.120 vs −0.269 ± 0.139 | −4.08 | 3 | 3 | 5 | Yes, but GSVA d does not carry over to a new sample set, so plan on d/2 or less. |
| GSVA Neftel MES1 | 0.240 ± 0.257 vs −0.236 ± 0.265 | −1.82 | 6 | 9 | 20 | Only with graft % as a covariate. |
| RP-block score (80 genes, centred rlog) | 0.128 ± 0.280 vs −0.128 ± 0.029 | −1.28 | 12 | 17 | 40 | Yes. The primary SD is almost all IL68B. |
| CIBERSORT MES | 0.514 ± 0.055 vs 0.390 ± 0.102 | −1.52 | 9 | 12 | 29 | **No:** non-specific (N168B). |
| CIBERSORT AC | 0.221 ± 0.030 vs 0.271 ± 0.044 | +1.33 | 11 | 16 | 37 | **No:** points the opposite way to GSVA AC. |

| Genome-wide DE, 1.5-fold, FDR 0.05 (measured common BCV 0.269; DESeq2 median BCV 0.194) | n/group |
|---|---|
| Hart 2013 + Jung 2005, depth 463 (π0 0.99 / 0.95) | 18 / 14 |
| ssizeRNA_vary (π0 0.99 / 0.95 / 0.90) | 21 / 18 / 14 |
| PROPER marginal power (5 simulations) | n 3: 0.11 (realised FDR 0.22) · 8: 0.57 · 20: 0.80 · 25: 0.84 |

| Other constraints | Number |
|---|---|
| Exact permutation p floor, two-sided (2 / C(2n, n)) | n 3: 0.10 · 4: 0.029 · 5: 0.0079 · 6: 0.0022. **At least 4 per group** is needed before any exact p can fall below 0.05. |
| CIBERSORTx group mode, **per class** (Newman 2019; Steen 2020) | c = 4: minimum 5, recommended 16–20. c = 7: minimum 8, recommended 28–35. Run per class, 3 vs 3 gives k = 3 < c: underdetermined. |
| CIBERSORTx high-resolution | c = 4: minimum 10, 32–40 at the authors' window. |
| Separating recurrence from burden (covariate r ≈ 0.7 with group, 1.5-SD effect) | ≈ 15 per group. Normal approximation: 15.70 / 1.5² × 1/(1 − 0.49) = 13.7, before the t correction. |

**Recommendation:**
- 12 animals per group covers the targeted endpoints (RP block, AC) at the observed effects. 18–25 per group covers genome-wide DE.
- Harvest primaries across a range of MRI volumes, so that graft % overlaps the recurrences. Scores can then be compared at matched graft % instead of adjusted.
- A non-thermal cytoreduction arm separates regrowth after cytoreduction from thermal injury (design in `../graft_relation/SUMMARY.md`).

Sources: Welch 1947 https://doi.org/10.1093/biomet/34.1-2.28 · Hart 2013 https://doi.org/10.1089/cmb.2012.0283 · Jung 2005 https://doi.org/10.1093/bioinformatics/bti456 · Bi & Liu 2016 https://doi.org/10.1186/s12859-016-0994-9 · Wu 2015 https://doi.org/10.1093/bioinformatics/btu640 · Newman 2019 https://doi.org/10.1038/s41587-019-0114-2 · Steen 2020 https://doi.org/10.1007/978-1-0716-0301-7_7

## 5. What cannot be done, or was not found

- **No public human post-LITT tumour transcriptome exists** (searched GEO, Europe PMC and the web; not exhaustive).
  - Chandar 2023 and Nielsen 2026 are single patients studied by immunofluorescence or IHC (https://doi.org/10.1097/cji.0000000000000485 ; https://doi.org/10.1111/nan.70071).
  - Campian 2026 profiled biopsies taken *at* LITT (https://doi.org/10.1038/s41467-026-69522-w). Its EGA dataset EGAD50000001639 lists no bulk RNA-seq, so whether those 5 profiles were deposited is not established.
  - As a result, no human data can test the LITT-specific part of any signal.
- **The GLASS Synapse matrices** need a login and acceptance of AR 9604985. The clinical_surgeries column list is unverified.
- **Controlled access (EGA), not downloadable:**
  - G-SAM: EGAS00001005436 / EGAD00001007860
  - Körber 2019: EGAD00001004564 (16 pairs)
  - Kim J 2015: EGAD00001001424
  - Kim H 2015: EGAD00001001396
  - Wang J 2016: EGAS00001000579
  - Klughammer 2018: EGAD00001004076
- **Not established:**
  - IDH per patient in Kim 2024.
  - The study of origin and IDH status of the CPTAC-3 recurrences.
  - GSE228497's primary/recurrent labels.
  - The GSE339484 count unit.
  - The IDH status of the 12 Wang 2016 SRA pairs.
  - The Varn 2022 AC-like effect size (the Cell PDF returned 403).
  - The full texts of Hoogstrate 2023 and Nagaraja 2026 (paywalled).
- **Nagaraja 2026** (same U251N/Visualase model, 4 vs 4, immune-compromised rats): whether its animals overlap ours, and how host reads were removed, is unknown. Ask the authors.
- **CIBERSORTx:**
  - The official manual is behind a login.
  - Whether the v104 Neftel signature came from Smart-seq2 or 10x data is unknown, so the choice of B- or S-mode is open.
  - The v104 fits used no batch correction.
  - With a single cell line, purified state profiles would stay hard to interpret at any n.
- **Local-data gaps:**
  - The xengsort logs are not on this machine; graft % was checked against `metadata_full.csv` only.
  - Whether the rat both/host ratio carries over from control to tumour libraries has not been tested. It can be: spike control rat reads into a tumour library at known fractions and re-run xengsort.
  - What the IL/NL sample prefix denotes is undocumented, and it explains more of PC1 than group does (R² 0.84 vs 0.68; SUMMARY C8).
  - The cause of a ~0.05-point difference between two RP-share recounts is unresolved.
