# Single-cell references for CIBERSORT on the human reads from U251 rat xenografts (primary vs LITT-recurrent)

**Bottom line.**
- **Use first:** the Neftel Smart-seq2 cells you already have on the grid. Keep only the 4,916 adult malignant cells and collapse them to four states. Garofano 2021 ships a second label set for these same cells.
- **For the recurrence question:** GBM-CARE is the best public reference. It is single-nucleus 10x data, so it needs a different batch mode (S-mode).
- **Limits to state up front:**
  - No human single-cell data exist from a glioblastoma that recurred after LITT.
  - No paper assigns U251 to a Neftel state.
  - The one outside validation of CIBERSORTx on Neftel states found it barely tracked imaging ground truth (mean r = 0.02).

Treat the fractions as a hypothesis until marker scoring and a U251 culture baseline agree with them. Every candidate the checkers tested turned out to exist. The few leads that were never checked are listed in section 4.

## 1. References to use now, ranked

### 1a. Neftel SS2 (on hand), with the Garofano per-cell labels

| Item | Fact | Source |
|---|---|---|
| Cells | The adult SS2 set is 5,742 cells: 4,916 malignant, 536 macrophage, 210 oligodendrocyte, 80 T-cell. Drop the non-malignant cells, because in your data those cells are rat. Drop the paediatric cells too; two of them come from recurrences (BT1187, BT786). | https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP393/annotations |
| How to get states | Give each malignant cell the highest of MESlike1, MESlike2, AClike, OPClike, NPClike1 and NPClike2 (cycling scores excluded). MES = max(MES1, MES2) and NPC = max(NPC1, NPC2), as in `scalop::as_four_state_gbm`. About 15% of cells are hybrids. | https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/ ; https://github.com/jlaffy/scalop |
| Outside check (passed) | I rebuilt the butterfly-plot X and Y coordinates from the six score columns. They match SCP393's own hierarchy file to within 1e-14 for 6,576 cells, so the columns are being read correctly. 287 cells do not join by name: 244 from BT85 and 43 named in lower case "mgh102". | https://singlecell.broadinstitute.org/single_cell/api/v1/studies/SCP393/clusters |
| Second scheme on the same cells | Supp. Table 4a (tab 14) gives GPM/MTC/NEU/PPR plus a Neftel label for 17,367 cells. 4,873 of them are adult Neftel SS2 malignant cells and join to your matrix by cell name. Counts across all three datasets: MTC 5,358, PPR 4,147, GPM 4,013, NEU 3,849. | https://static-content.springer.com/esm/art%3A10.1038%2Fs43018-020-00159-4/MediaObjects/43018_2020_159_MOESM3_ESM.xlsx |
| Gene lists | Table S2 downloads openly from Cell: MES2 50, MES1 50, AC 39, OPC 50, NPC1 50, NPC2 50, G1/S 29, G2/M 45 genes. | https://ars.els-cdn.com/content/image/1-s2.0-S0092867419306877-mmc2.xlsx |
| Build recipe | GLASSx `neftel2019_signature_matrix.r`: CIBERSORTx signature with quantile normalisation off, k.max 999, q 0.01, 300–500 genes, B-mode. The matrix it produced was never published. | https://github.com/fsvarn/GLASSx/blob/master/R/single_cell/signature_matrix/neftel2019_signature_matrix.r |
| Published precedent | In IDH-wt patients, AC-like cells decreased at recurrence (P = 7e-3, paired t-test). That was standard therapy, not LITT, and the mixtures included host cells. | https://pmc.ncbi.nlm.nih.gov/articles/PMC9189056/ |
| Fit | Adult IDH-wt matches U251 (75-year-old male, no IDH1 mutation reported). Full-length SS2 is the closest platform to bulk poly-A RNA-seq. Limit: the tumours are untreated primaries, so this reference cannot show recurrence-specific states. | https://api.cellosaurus.org/cell-line/CVCL_0021?format=txt |
| Cells per state | Not counted in this search. Count after you assign states. | — |

**CIBERSORTx inputs** (web tool or the `cibersortx/fractions` container):
- **Single-cell reference file:**
  - tab-delimited, genes as rows, one column per cell;
  - row 1 holds each cell's label: MESlike, AClike, OPClike, NPClike (no hyphens or periods);
  - values in linear TPM, computed as 10·(2^E − 1) from the log values;
  - at least 3 cells per label, and no unlabelled cells;
  - 4,916 cells is under the recommended 5,000 (the hard limit is 10,000).
- **Mixture file:** genes × 6 samples in linear TPM, normalised the same way as the reference.
- **Settings:** `--single_cell TRUE`, B-mode, quantile normalisation off.
- Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/

**Classic `CIBERSORT.R` inputs:**
- The script cannot build a signature matrix itself ("use java version").
- Give it a linear genes × 4 matrix. Either export one from CIBERSORTx, or build it with the Newman 2015 recipe:
  1. For each state against the rest, run a two-sided unequal-variance t-test and keep genes at q < 0.3.
  2. Rank the kept genes by fold change and take the top G per state, for G from 50 to 200.
  3. Keep the matrix with the lowest `kappa()`, and report that condition number.
- Run with `QN = FALSE` and `perm` of at least 100. The script has no batch correction.
- Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC4739640/

### 1b. GBM-CARE (Nomura 2025, Spitzer 2025): the primary vs recurrent reference

| Item | Fact | Source |
|---|---|---|
| Files | `GSE274546_RAW.tar` (4.27 GB) holds 121 files named `GSM*_P#T#_umi_counts.RDS.gz`. They are gzipped twice: gunzip once, then `readRDS`. Each is a dense UMI matrix; P32T1 is 33,538 genes × 372 nuclei. Labels are in `malignant_meta_data_2025_01_08.RDS` (columns CellID, Sample, ID, State, isCC). Clinical data are in `spitzer_supptable1.xlsx`. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546 ; https://github.com/dravishays/GBM-CARE-WT |
| Cohort | 59 patients and 121 samples, all IDH-wt. Take status from the "Primary or recurrent" column, stripping trailing spaces first: 56 primary, 55 first recurrence, 10 second recurrence. Do not use the T1/T2 labels for status: for patients FR13, FR15 and NL08, T1 is already a first recurrence. | https://raw.githubusercontent.com/dravishays/GBM-CARE-WT/HEAD/data/spitzer_supptable1.xlsx |
| Nuclei per state (246,408 malignant) | AC 49,669; OPC 45,983; GPC 34,991; Unresolved 25,763; MES 23,248; Hypoxia 20,496; NEU 16,388; NPC 13,593; Stress 12,289; Cilia 3,988. At T1: AC 30,110, OPC 26,088, MES 13,558, NPC 8,697. At T2: AC 19,212, OPC 19,202, MES 7,873, NPC 4,860. | https://raw.githubusercontent.com/dravishays/GBM-CARE-WT/HEAD/data/malignant_meta_data_2025_01_08.RDS |
| States | The labels are shipped, not derived. Either keep all 9 states or collapse to Neftel's four; drop "Unresolved". | https://pmc.ncbi.nlm.nih.gov/articles/PMC12081307/ |
| Limits | (1) This is nuclear RNA and your bulk is whole-cell. Nomura reports GPC-like and NEU-like cells are enriched in nuclear data. (2) Recurrence followed standard care, not LITT, and 6 of the 55 first recurrences had no radiation. (3) Sample P58T3 has no malignant nuclei. (4) The CellID-to-matrix join was checked on P32T1 only (89 of 89 matched). | https://pmc.ncbi.nlm.nih.gov/articles/PMC12081298/ |

**Inputs:**
- **CIBERSORTx:** S-mode, because this is UMI data. S-mode also requires the single-cell reference file. Downsample to 5,000 nuclei or fewer, balanced by patient and state. Build per-state profiles in linear space.
- **Classic CIBERSORT:** build the signature yourself. Classic CIBERSORT cannot correct the gap between nuclear reference and whole-cell bulk.

### 1c. GBmap core and extended: whole cells, many studies, labels shipped

| Item | Fact | Source |
|---|---|---|
| Core | 338,564 cells from 110 donors and 16 studies, of which 127,521 are neoplastic. Level-3 states: AC 50,847; MES 33,167; NPC 22,117; OPC 21,390. Almost all recurrent malignant cells come from one donor (Johnson SM011: 4,563 cells, 4,517 of them NPC-like). | https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad |
| Extended | 1,135,677 cells from 240 donors, of which 403,575 are neoplastic: AC 132,265; OPC 118,334; MES 105,398; NPC 47,578. There are 38,993 recurrent malignant cells, and 79.5% of them come from Abdelfattah (5 donors: OPC 18,893, AC 6,254, MES 5,392, NPC 464). Labels for the added datasets were transferred by scArches, not curated, and there is no level 4. | https://datasets.cellxgene.cziscience.com/230c8701-7291-4dd5-bf36-688490c681ee.h5ad |
| Filtering | Keep `annotation_level_1 == "Neoplastic"` and `suspension_type == "cell"`. In `raw.X`, 10x cells hold UMI counts but Smart-seq2 cells hold rounded TPM, so use one assay family at a time. Do not concatenate core and extended: they overlap. | https://pmc.ncbi.nlm.nih.gov/articles/PMC12526130/ |
| Inputs | CIBERSORTx S-mode, with 5,000 or fewer cells balanced by donor and state. Classic CIBERSORT: build the matrix yourself. Label-transfer models are on Zenodo. | https://zenodo.org/records/6962901 |

## 2. Secondary references

| Reference | What it adds, and its limits | Source |
|---|---|---|
| Neftel 10x arm (GSM3828673) | 16,201 cells from 9 adult tumours. The GEO file has 30,314 genes on a linear, CPM-like scale. Malignant calls come from 3CA (11,172 cells); states from GBmap core (7,832 neoplastic: AC 2,969, MES 2,101, NPC 1,814, OPC 948). Useful as a platform-sensitivity check. It shares donors with the SS2 set, needs S-mode, and is primary tumour only. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3828673 ; https://www.weizmann.ac.il/sites/3CA/brain |
| Wang L 2022 (GSE174554) | Paired primary/recurrent snRNA: 81 GSMs (33 primary, 48 recurrent), with 31 pairs having both timepoints. Only Tumor/Normal is shipped per nucleus. The cohort is not in GBmap, so assign states with the GBmap Azimuth/scArches models. 7 of the 81 GSMs are labelled "Brain tumor" rather than GBM. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554 |
| Abdelfattah 2022 (GSE182109) | Recurrent tumours as whole-cell 10x. Drop LGG-03/04 and ndGBM-10/11. Primary and recurrent tumours are from different patients, with no IDH or treatment data. The GBmap barcode set (207,377) differs from the authors' (201,986), so check the overlap before transferring labels. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182109 |
| GSE211376 | 11 primary IDH-wt snRNA samples with labels: OPC 5,201, MES 2,790, AC 2,752, NPC 2,687. GBmap's labels for the same nuclei disagree (11,727 vs 13,430 malignant). Useful to test a nuclear reference against a whole-cell one. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211376 |
| Richards 2021 (Developmental vs Injury Response) | Per-cell AUC scores on SCP503 (sign-in needed). Table 7 lists 3,968 Developmental and 4,399 Injury Response genes. It is a gradient, so score it (GSVA/ssGSEA) rather than deconvolve. | https://www.nature.com/articles/s43018-020-00154-9 |
| Scoring cross-checks | Garofano Table 6j (4 × 50 genes). GBM-CARE Table S2 (15 × 50 genes; drop MP_12_LQ and MP_11_MIC). GBMdeconvoluteR marker lists, scored with MCP-counter (MCP_GBM mean r = 0.43 against imaging mass cytometry). | https://github.com/ajxa/GBMDeconvoluteR ; https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02167-5/MediaObjects/41588_2025_2167_MOESM2_ESM.xlsx |
| Ready-made matrices for `CIBERSORT.R` | Ivy GAP Table S14: 293 genes × CT/CTmvp/CTpan/LE, linear. These are anatomic compartments, so exploratory only; map old symbols (e.g. AGXT2L1 → ETNPPL). Varn Table S3: 4,130 genes × 12 classes, of which only 3 are tumour states; 9 columns are host cell types absent from your reads. Repair 7 symbols Excel turned into dates and 256 dot-for-dash symbols. It was built for S-mode with a reference file that was never released. | https://www.science.org/doi/suppl/10.1126/science.aaf2666/suppl_file/aaf2666_table-s14.xlsx ; https://ars.els-cdn.com/content/image/1-s2.0-S0092867422005360-mmc3.xlsx |
| U251 in vitro, single-cell (GSE178114) | U251 on 10x, but the matrices are raw and unfiltered (33,538 genes × 6,794,880 barcodes), so you must call cells yourself. Score the empty-vector arm (GSM5380412) to see which states cultured U251 occupies. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178114 |
| U251 culture vs rodent brain, bulk (GSE134470) | HuGene 1.0 ST arrays: U251_C1–C3 in culture and U251_X1–X2 from mouse brain, FACS-sorted to human cells. Run both through your signature as an outside check of the culture-to-brain shift. Caveat: the paper gives U251's source as ATCC HTB-14, which is U-87MG's ATCC number. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134470 |
| Human Protein Atlas, U-251MG | In culture: GFAP 1,739.8, SOX2 166.5, S100B 161.8, NES 136.0 nTPM. Check these genes in your xenograft reads. | https://www.proteinatlas.org/ENSG00000131095-GFAP/cell+line |
| LeBlanc 2022 (GSE173278) | Tissue, explant and gliomasphere line from the same patients (malignant cells: 23,556 / 34,006 / 10,137). Run your signature on each as pseudo-bulk to see how composition collapses in a line. States are available from GBmap extended. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173278 |
| Manoharan 2024 | Human tumour programmes learned from rodent xenografts. They are cNMF spectra, not expression, so use them as markers only. | https://zenodo.org/records/10437583 |

## 3. Traps that give a finite but wrong answer

| Trap | What goes wrong | Fix | Source |
|---|---|---|---|
| Log scale | SCP393 SS2 values are log2(TPM/10+1), and the GEO "TPM" SS2 file is the same log matrix. CIBERSORTx treats any matrix with max < 50 as log2 and anti-logs it, producing TPM/10+1: every zero becomes 1 and marker contrast shrinks. `CIBERSORT.R` never anti-logs the signature. | Upload 10·(2^E − 1) and do all averaging in linear space. | https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/ |
| Mixed scales | The GEO 10x file is linear. GBmap `raw.X` mixes UMI counts with rounded TPM. | Use one assay family per reference. | https://datasets.cellxgene.cziscience.com/861acfd8-25f0-418b-a445-aa96da232827.h5ad |
| Label strings | CIBERSORTx drops text after a period. R's `check.names` turns "MES-like" into "MES.like" (truncated to "MES") and "T-cell" into "T". The GLASSx script does this, and ties in `which(max)` can misalign labels. Its run scripts also contain another user's hard-coded CIBERSORTx token. | Use labels with no hyphens or periods. Check for ties. Use your own token. | https://github.com/fsvarn/GLASSx/blob/master/bin/cibersortx/cibersortx_fractions_run_neftel.sh |
| Batch mode | Smart-seq2 references take B-mode; 10x/UMI references take S-mode. Both need at least 3 mixtures and 10 are recommended; you have 6. Expression-purification modes (group and high-resolution) need 4–5 times as many mixtures as cell types. | B-mode for SS2. Run with and without B-mode and compare. Skip purification. | https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/ |
| Gene symbols | The Neftel matrix is hg19 RSEM (AARS → AARS1, AAED1 → PRXL2C). The marker lists use old names (C8orf4 → TCIM, HN1 → JPT1, SEPT3 → SEPTIN3, KIAA0101 → PCLAF). Unmapped genes drop silently when reference and bulk are intersected. | Map through HGNC previous symbols, then count matched genes per state. | https://rest.genenames.org/search/prev_symbol/AARS |
| Wrong status labels | GBM-CARE T1 is not always primary, and the status column has trailing spaces. GBmap stage is "unknown" for about 40% of cells, and core's recurrent block is effectively one patient. | Use the stripped clinical column. Do not treat GBmap core as a recurrence reference. | https://raw.githubusercontent.com/dravishays/GBM-CARE-WT/HEAD/data/spitzer_supptable1.xlsx |
| Diagnostics that cannot fail | The CIBERSORT p-value tests "no signature cell type present", which is always rejected when every read is tumour. Absolute mode has no meaning here. Script defaults are `perm = 0` (prints 9999) and `QN = TRUE`. Fractions are forced to sum to 1, so a programme the reference lacks, such as injury after ablation, gets redistributed across the four states. | Relative mode, `QN = FALSE`, `perm` ≥ 100. Compare per-sample correlation and RMSE between arms (my inference, not a published rule). | https://pmc.ncbi.nlm.nih.gov/articles/PMC4739640/ |
| Cell line vs patient states | Gliosphere lines sit in a more uniform state than their parent tumours (LeBlanc 2022). In Pine 2020, xenografts were enriched for AC-like and OPC-like cells and less NPC/OPC-like than primaries. U251 and U87 converge when grown intracerebrally (Camphausen 2005). If U251 occupies essentially one state, the four columns are nearly collinear, and a primary-vs-LITT difference could be noise shifted between them. | Report the condition number. Run in-vitro U251 (GSE178114, GSE134470, DepMap ACH-000232) through the same signature as a baseline. | https://pubmed.ncbi.nlm.nih.gov/35303420/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10256258/ ; https://pubmed.ncbi.nlm.nih.gov/15928080/ |
| No ground truth for U251 | No Neftel-state call for U251 has been published. Bulk subtype calls conflict: classical (Schnöller 2023) vs proneural (Schulz 2022). | Establish the state yourself; there is no literature value to check against. | https://pmc.ncbi.nlm.nih.gov/articles/PMC10007763/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9347152/ |
| Proliferation | U251 doubles in about 23–24 h, and Neftel's cycling cells are enriched in OPC/NPC-like states. If cell-cycle genes enter the signature, a fast line reads as OPC/NPC-high. | Untested suggestion: drop G1/S and G2/M genes, or add a cycling class, and report both runs. | https://pmc.ncbi.nlm.nih.gov/articles/PMC6703186/ |
| Hypoxia | MES2 is the hypoxia-dependent programme (HILPDA, DDIT3, ENO2, LDHA, ADM). Necrosis after ablation could raise MES through MES2. This link is my inference. | Report MES1 and MES2 separately before collapsing them. | https://ars.els-cdn.com/content/image/1-s2.0-S0092867419306877-mmc2.xlsx |
| Known weak performance | GBMdeconvoluteR's CIBERSORTx on Neftel states correlated with imaging mass cytometry at mean r = 0.02 (−0.24 in recurrent tumours, n = 5). Marker-based MCP_GBM reached 0.43. | Always pair CIBERSORT with marker scoring. | https://pmc.ncbi.nlm.nih.gov/articles/PMC10326489/ |
| Zebrafish U251 cells (GSE147526) | These 21 cells show zero GFAP, SOX2 and S100B reads, contrary to the Human Protein Atlas profile for U-251MG, so their identity is doubtful. Do not use them as evidence that U251 lacks AC/OPC/NPC markers. | Check those genes in your own reads. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147526 |
| Immune references | Rat immune cells were removed with the host reads. Any non-zero LM22 fraction means collinearity with tumour genes, not immune cells. | Use LM22 as a negative control only. | https://pubmed.ncbi.nlm.nih.gov/25822800/ |
| Species sorting (inference) | xengsort drops reads it classes as "both" or "ambiguous", so genes conserved between human and rat may lose reads unevenly. | Put human-only in-vitro U251 reads through the same pipeline and compare per-gene TPM with salmon on all reads. Exclude genes with large losses. | https://pmc.ncbi.nlm.nih.gov/articles/PMC8017614/ |
| Stock identity | SNB-19 and U-373MG (ATCC) are U251 derivatives; U-373MG (Uppsala) is a different line. Engineered U251 IDH1-R132 knock-ins also exist. | STR-type your stock. | https://api.cellosaurus.org/cell-line/CVCL_0021?format=txt |

## 4. What could not be established

- **LITT data:** no human single-cell data from glioblastoma recurring after LITT were found. Two related items were not verified and neither provides usable data:
  - Nagaraja 2024 SNO abstract: U251N rat xenografts ablated with Visualase (https://pmc.ncbi.nlm.nih.gov/articles/PMC11552993/).
  - Tao 2026 J Immunol: scRNA of rat cells after sublethal LITT (PMID 41764740); its data availability was not checked.
- **U251 in a rodent brain:** no single-cell data were found. Knudsen 2022 (nude rat, patient-derived spheroids, primary vs recurrent scRNA) is the same design, but no accession was found (https://pmc.ncbi.nlm.nih.gov/articles/PMC9248408/).
- **A published Neftel four-state CIBERSORT(x) matrix:** none exists. The GLASSx output file and GBMdeconvoluteR's 5,075-cell matrix are both unpublished, so you will have to build one.
- **Per-state counts** of the 4,916 adult SS2 malignant cells were not tabulated.
- **Score centring:** whether the SCP393 scores were centred per tumour or across all cells is not stated.
- **Method validity:**
  - CIBERSORTx has not been validated on plastic, continuous malignant states; Newman 2019 validated distinct cell types.
  - Whether six mixtures give a stable B-mode correction is unknown.
- **Joins checked only partly:**
  - GBM-CARE CellIDs, beyond sample P32T1.
  - GBmap cell IDs to GEO for the Neftel 10x cells.
  - Garofano's other 12,494 cells to GSE117891 and GSE103224.
- **Rat macrophages:** whether they drive human U251 toward MES through OSM (the Hara 2021 mechanism) was not established.
- **Other gaps:**
  - The GLASS-NL (Hoogstrate 2023) snRNA accession was not found.
  - Al-Dalahmah 2023: GEO lists 10 primaries and the paper 8, and the Google Drive objects could not be opened.
  - Couturier 2022 recurrent data were never deposited: EGA holds only the 2020 dataset and the GitHub repository is empty.
  - The DepMap U-251MG file was not checked. MIX-seq keys U251MG to ACH-000978, which Cellosaurus assigns to an endometrial line; use ACH-000232.
- **Checked and not useful as references:**
  - U251 is absent from Kinker 2020, MIX-seq (listed but never profiled) and Tahoe-100M.
  - STAMP U-373 MG is a targeted panel of 500–5,001 genes.
  - GSE225775 has no U251.
  - Wang R 2020 cells were CD133/GLAST-sorted, which biases state proportions.
  - Lee 2021 recurrences are confounded by anti-PD-1.
  - Wu 2023 is controlled access, with only 3 newly sequenced recurrences.
  - Johnson 2021 has one IDH-wt recurrence, and it is 99% NPC-like.