# invasion_niche: U251 in the contralateral hemisphere, its state, invasion programmes by arm, host niche

Item S5 of the U251N LITT multi-omics build (SELECTION.md). This is custom R on released-pipeline outputs:
- the nf-core/rnaseq **3.22.2** graft-only run (`M/rnaseq_graft`);
- the existing host run (nf-core/rnaseq 3.22.2 on host-only reads, `H/results_rnaseq`), and its relabelling fits
  (`H/verify`).

No Nextflow pipeline covers these questions, so none is launched.

This folder holds only code, parameters and a sample sheet with library names and grid paths. Results live on the grid
at `M/invasion_niche/results/`; `M/invasion` is a symlink to `M/invasion_niche`, because SELECTION.md names the item
`M/invasion`. `M` = `/rs/rs_grp_oschome/go2432/u251_multiomics`, `U` = the grid clone
`/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study`, `H` = `/rs/rs_grp_oschome/go2432/u251_host`.

## Claim ceiling (copied into REPORT.md)

- U251 RNA is present in the contralateral hemisphere of rats 68 and 69, both in the primary arm.
- Migration cannot be separated from carry-over during dissection.
- Nothing can be said about LITT and invasion: no ablated rat has contralateral tissue.
- With 3 primaries against 3 recurrences there are exactly 10 balanced relabellings. No primary-versus-recurrent result
  can go below p = 0.10.

## Questions and methods

### A. How much U251 is in each sample?

The readout is human chrY dosage. It does not depend on xengsort's call in the controls, because the rats are female
and U251 is male.

- **Genes.** The 11 X-degenerate Y genes: RPS4Y1, DDX3Y, KDM5D, UTY, USP9Y, EIF1AY, ZFY, TXLNGY, NLGN4Y, PRKY and
  TMSB4Y.
  - Their ids are resolved from the GENCODE v44 GTF by `gene_name`, on chrY, excluding `_PAR_Y` ids. The ids are never
    typed in.
  - TXLNGY resolves to two overlapping ids, an lncRNA and a transcribed pseudogene, and both are summed.
  - The ids used are listed in `A_y_dosage/gene_ids_from_gtf.tsv`.
- **Calibration.** Y counts per million graft pairs, in each of the six tumours. The range is reported, and C2B is
  shown beside it but is not part of the calibration.
- **Application.** The calibration is applied to IL64B, N168B and N269B.
  - Each gets an exact Poisson 95 % interval on its Y count, combined with the range of the tumour calibrations.
  - The result is U251 graft-bin-equivalent pairs as a share of all the library's pairs.
  - This is done on the published graft+both counts (`U/ANALYSIS/results_human_final`) and on the graft-only counts,
    and both are reported.
- **Y/X gametolog ratios** (RPS4Y1/RPS4X, DDX3Y/DDX3X, KDM5D/KDM5C, UTY/KDM6A, EIF1AY/EIF1AX, ZFY/ZFX, USP9Y/USP9X),
  plus a pooled ratio:
  - computed per library and per count source;
  - reported with the rank of the arm difference (Welch t, recurrent minus primary) among the 10 splits.
  - Cellosaurus records Y-negative U-251 subclones, so a shift in these ratios would be a clonal marker.
  - Written to `results/y_x_ratios.tsv`, which the rnavar variant integration reads.

### B. Are the contralateral U251 cells in an in-vivo or a culture-like state?

N269B is compared with IL69B and C2B on graft-only counts.

- **Universe.** Protein-coding genes (chrM excluded) whose mean CPM is at least 10 across C2B and the six tumours.
  - Genes are dropped where the graft-only background of IL64B or N168B, scaled to N269B's host read depth, exceeds 10 %
    of N269B's count. These are background-flagged genes.
- **Depth-matched null.** IL69B, C2B and the other five tumours are binomially thinned to N269B's deduplicated depth,
  in 200 draws.
  - The depth is D = N269B's universe count × (1 − PERCENT_DUPLICATION), where PERCENT_DUPLICATION comes from Picard
    MarkDuplicates in the graft-only run.
  - Thinning to the raw count would understate N269B's sampling noise, which is set by unique molecules.
- **Axis.** A culture-to-in-vivo axis position:
  - The axis genes are those where the five tumours other than IL69B all sit on the same side of C2B, by at least
    2-fold on average.
  - The position is the median over those genes of (log2 x − log2 C2B)/(log2 tumours − log2 C2B). It is 0 for C2B and
    1 for the tumours.
- **Verdict rule, fixed before the run.**
  - **In-vivo U251 state, not culture-like:** N269B's position lies inside the 95 % interval of the thinned in-vivo
    profiles (thinned IL69B, plus each other tumour thinned and placed on an axis built without it), and outside the
    95 % interval of thinned C2B.
  - **Culture-like:** the reverse.
  - **Not resolved:** anything else.
- **Also reported:**
  - Spearman correlations against every reference, with the thinned references' ranges.
  - NNLS of N269B on C2B plus the six tumours, leaving each donor out. The rows are scaled by the square root of the
    reference mean, and the same fit is run on thinned IL69B.
  - A gene-by-gene comparison with IL69B. A difference counts only if it lies outside the thinned-IL69B minimum and
    maximum **and** outside the six-tumour spread. The rule's own null rate is measured on the draws.
  - The libraries that share an i7 or i5 index with N269B, read from `M/reads/index_seq`.
- The claim stops at "in-vivo U251 state, not culture-like" (or its negation). Carry-over looks the same as migration,
  so the result never says "migrated".

### C. Do LITT recurrences carry more or less of published glioma infiltration programmes?

Scored by GSVA and fgsea, each ranked against the 10 relabellings.

**GSVA**
- Gaussian kernel (GSVA's default for continuous data), on DESeq2 vst of the graft-only counts, never TPM.
- All 10 libraries are scored. Two sensitivity versions are also written: the six tumours scored alone, and log2(CPM+1).
- The primary-versus-recurrent statistic is Welch t on the six tumours, ranked among the 10 splits.

**fgsea**
- Run on each split's Wald statistic, 10,000 permutations, sets of 10-500 genes.
- The fits are human_de_null pass B's (`M/de_null/results/passB`) when they are present and hold exactly the 10
  balanced splits. Otherwise the job fits its own DESeq2 models (listed below), and the report says which fits were used.
- On graft+both counts the job uses pass A in the same way.

**Gene sets**, from files only. Each file is copied into `results/genesets/`, and `SOURCES.tsv` gives its URL and md5.

| set | file | rule |
|---|---|---|
| `ivygap_LE` (primary), `ivygap_CT`, `ivygap_CTmvp`, `ivygap_CTpan` | Puchalski 2018 Table S14 (PMID 29748285), `U/ANALYSIS/cibersort/signatures_online/ivygap2018_anatomic4.txt` | a gene goes to the compartment where its value is largest; checked identical to the frozen `ANALYSIS/spread/signatures` sets |
| `neftel2019_MES1`, `MES2`, `AC`, `OPC`, `NPC1`, `NPC2`, `G1S`, `G2M` | Neftel 2019 Cell Table S2, "Meta-module gene lists" (the published lists), `mmc2.xlsx` | as published |
| `neftel_repo_MES`, `AC`, `OPC`, `NPC` | the repository's CIBERSORT signature built from the Neftel SS2 cells in `U/ANALYSIS/refs/neftel_2019` | largest column; derived, not the published lists |
| `venkataramani2022_invasivity_high` (primary), `_low` | Venkataramani 2022 Cell Table S2a (PMID 35914528), `mmc2.xlsx` | see below |
| `HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` (primary), `ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE` (primary) | msigdbr 26.1.0 (MSigDB 2026.1) in `U/ANALYSIS/refs/R_libs` | if msigdbr will not load in the image, the MSigDB 2023.2 GMT files are used; the report says which, and compares the two when both load |

How the Venkataramani groups were oriented:
- Table S2a labels its 172 genes "Correlated" (151) or "Anticorrelated" (21). The correlations are with pseudotime.
- Table S2e lists the 20 genes most correlated with the invasivity score, and all 20 are in the "Anticorrelated"
  group. So that group is taken as the invasive pole.
- This orientation comes from the file's own legends, not from the paper's text.

Optional sets that could not be obtained as files, and were not transcribed:
- Darmanis 2017: its Table S6 is an embedded image.
- Hoelzinger 2005 and Gill 2014: the tables are in the article body only; PMC does not allow full-text download, and
  pnas.org returns 403.
- Mariani 2001: no free full text.

The supplements were fetched on 2026-09-27 and converted cell by cell by `convert_supp.py`. The xlsx files and the
converted TSVs are on the grid in `M/invasion_niche/genesets_src/` and are not in this repository.

**Leak checks.** Every set is scored on graft-only and on graft+both counts. A set is flagged as measuring leak if any
of these holds:
- it scores higher in IL64B or N168B than in every tumour (graft-only);
- it tracks m = both/(graft+both) across the 10 libraries (Spearman, exact p < 0.05);
- its true-split rank is 1 or 2 on graft+both but 4 or worse, or of the opposite sign, on graft-only. That is, it moves
  only with the both bin.

The Ivy GAP LE set is the known risk: CIBERSORT on graft+both called IL64B 100 % leading edge.

**Reporting.**
- The scores are reported beside tumour content (graft %). The DNA tumour fraction is joined later and is not
  adjusted away.
- The four primary sets were declared before the run.
- With 21 or so sets, about 2 are expected to reach rank 1 of 10 by chance, and the report prints that expectation.
- Expected outcome: nothing reportable.

### D. Host niche: is host Mmp13 the one host change, and what tissue is the N269B piece?

- **Mmp13 rank.** The rank of the true split's |Wald| among the 10 host relabellings, from the existing host fits.
- **Mmp13 against cell-type markers.** Mmp13 is also taken per unit of Col1a1, Postn and Pdgfrb and of their geometric
  mean, and per unit of Cd68, Aif1 and Ptprc and of their geometric mean. Each of these is ranked among the 10 splits,
  because Mmp13 may simply follow stroma or myeloid content.
- **Host origin.** Human MMP13 in the graft-only and graft+both counts is shown beside rat Mmp13, to show the transcript
  is host.
- **Marker panels** (rat ids resolved from the rat GTF by name): myeloid, reactive glia, fibroblast and collagen,
  meninges, neuron, white matter. They compare N269B with N168B, IL64B and the tumours.
- Mmp13 stays "exploratory; needs ISH".

## Checks against something outside the computation

| check | against |
|---|---|
| Part A denominators: the 10 xengsort totals equal the 10 SRA spot counts of SRR39547363-72 | NCBI SRA runinfo (read 2026-09-27; `sra_spots` in params.yaml) |
| Part A on the published counts reproduces the preliminary shares (N269B 3.4-4.2 %, N168B about 0.31 %, IL64B at most 0.03 %) and the calibration (157-193 per million graft pairs; C2B 211) | the earlier tally; the local fixture test reproduced 3.903 % (3.401-4.189), 0.307 % and 0.025 % |
| Female host: rat chrY gene CPM in the host reads is about 0 | the host transcriptome, independent of the human reads |
| The human GTF header says "version 44 (Ensembl 110)" | GENCODE's own header |
| Ivy GAP argmax sets are identical to the frozen spread/signatures sets | sets frozen, with sha256, before any score |
| msigdbr membership against the MSigDB 2023.2 GMT files | a second MSigDB release |
| Supplement md5 at run time against the md5 at fetch | the fetch record |
| Gate: rnaseq_graft checks 1 (clip once), 2 (read counts from xengsort to STAR), 3 (reverse strand), 10 (no empty outputs) | the xengsort logs, the delivered read lengths, and the published strandedness |

## Settings that are not defaults, and why

| setting | value | why |
|---|---|---|
| counts | graft-only for B, C, arm statistics; graft+both beside them | the recurrent arm carries about twice the rat `both` contamination of the primary arm |
| thinning target | N269B universe counts × (1 − PERCENT_DUPLICATION) | sampling noise is set by unique molecules; duplicates were 65-79 % in the tumours |
| draws | 200, seed 20260927 | as specified; the min-max rule's null rate is measured on the draws themselves |
| universe | protein-coding, mean CPM ≥ 10, background share ≤ 10 % | keeps genes N269B's depth can measure and drops genes the near-empty controls fill |
| axis genes | all five other tumours on one side of C2B, mean |log2| ≥ 1 | a culture-versus-brain difference shared by every tumour, not one animal's |
| NNLS | exact, by enumerating supports; rows / sqrt(reference mean CPM) | exact for 7 columns; Pearson scaling for counting noise |
| GSVA | kcdf Gaussian on vst; minSize 10, maxSize 500 | vst is continuous; TPM inflated the translation result before |
| own DESeq2 fits (only if human_de_null's are absent) | six tumours, `~ group`, Wald, parametric, minReplicatesForReplace 99, genes with ≥ 10 counts in ≥ 3 of 6 | the published run's settings, where they apply; no length offsets (the fits only rank genes for fgsea) |
| fgsea | fgseaSimple, 10,000 permutations, sets 10-500, seed | as host_verify, with a smaller minimum because Neftel AC and G1/S have 29-39 genes |
| leak flag | Spearman with m, exact p < 0.05 (n = 10) | stated before the run |

## Run

- `run.sbatch` requests reqp with `--qos=requeue`, 4 CPU, 16 GB and 2 h, with `--output` under `M/logs/`. It runs:
  1. the gate;
  2. the GTF gene tables, built once into `M/refs/gene_tables/` (idempotent, keyed on the source md5);
  3. gene-set staging;
  4. step A (`bin/inv_stepA.R`, DESeq2 1.34.0 image: parts A, B and D, then vst and fits);
  5. step C (`bin/inv_stepC.R`, `go2432-bioconductor.sif` with `/rs/rs_grp_oschome/go2432/Rlib/R4.3`, installing GSVA
     or fgsea there only if missing);
  6. the report (`bin/inv_report.R`).
- `grid.config` holds the Singularity settings. Images are looked for in `M/singularity`, the shared cache, and then in
  the existing cache `/rs/rs_grp_oschome/go2432/singularity_cache`, and their sha256 is recorded. Sandboxes go under
  `/dev/shm` and are removed on exit.
- The gate: while the rnaseq_graft head or check job is live, the job resubmits itself `afterany` on the check job and
  exits 0.
- `submit.sh`, run on the login node, submits with `afterany` on the rnaseq_graft check job.
- There is no `-profile test` run, because no pipeline is launched. The code was run end to end locally on a fixture:
  - the published counts, with a synthetic graft-only table made by thinning each library by its graft share;
  - the real host tables;
  - R 4.5.2 / DESeq2 1.50.2 / GSVA 2.4.9 / fgsea 1.36.2.

  On the grid the versions are those of the images, and all are recorded.
- Expected wall time is 25-45 min.

## Outputs (`M/invasion_niche/results/`)

`REPORT.md` states, for each part, what it can and cannot show. The other outputs:
- `STATUS_A.tsv`, `EXIT_CODES.tsv`, `summary_A.tsv` and `summary_C.tsv`;
- `y_x_ratios.tsv`;
- `A_y_dosage/`, `B_contralateral_state/`, `C_programmes/` (expr, gsva, fits_own, fgsea, leak checks) and
  `D_host_niche/`;
- `genesets/`;
- `inputs/`: the gate record, input md5s, image sha256s, session info and step logs.

## What the result can and cannot show

**A**
- It can show how much human Y RNA, which a female host cannot make, sits in each contralateral sample, as a U251 share
  with an interval.
- It cannot show the source of that RNA (migration, carry-over or index hopping) or cell numbers.

**B**
- It can show whether N269B's U251 RNA looks grown-in-brain or culture-like at N269B's depth.
- It cannot show migration.
- Index-hopped reads from the libraries that share an index with N269B would carry those libraries' state.

**C**
- It can show whether any published infiltration programme separates the arms more than a relabelling does, and
  whether a separation survives removal of the both bin.
- Nothing can go below p = 0.10.
- Arm is confounded completely with tumour content.

**D**
- It can show whether the Mmp13 fall outranks all relabellings and survives division by stromal or myeloid markers,
  that the transcript is rat, and what compartments N269B's piece is enriched for.
- It cannot show where Mmp13 is expressed; that needs ISH.
