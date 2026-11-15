# human_de_null: does the published human recurrence result stand out among the relabellings?

Item 4 of the multi-omics plan (SELECTION.md S4). Custom R in the published pipeline's own DESeq2 container, plus
fgsea. No nf-core pipeline tests a result against its own relabellings, so this is a small documented step, not a
pipeline run.

## Question

The published human result is nf-core/differentialabundance 1.5.0 (`ANALYSIS/results_therapy_v3`), Recurrent_U2
against Primary_U2, 3 v 3: 19,351 genes tested, 107 at padj < 0.05, 82 in the filtered table. Its gene-set result is
a fall of the ribosomal-protein block in recurrence (KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION NES -1.99, FDR
q 0.022 at seed 1234, and 0.02-0.31 across seeds), whose size leans on one primary, IL68B.

With three tumours per arm there are exactly 10 ways to split the six into two groups of three. If the true split
(primary against recurrent) gives no more genes, and no stronger gene-set signal, than the nine meaningless splits,
the published result cannot be told apart from how these six tumours differ anyway. Without this test no human
recurrence claim is reportable. The host (rat) side has already failed it: the true split ranked 2 of 10 on genes
and 9 of 10 on gene sets (`/rs/rs_grp_oschome/go2432/u251_host/verify/results/summary.json`).

**The floor.** Rank 1 of 10 is p = 0.10. No result from this design can be smaller, whatever the tool.

## What runs

Each pass is one Slurm job running `run.sbatch`, with two steps.

1. `bin/human_de_null_de.R` runs in
   `depot.galaxyproject.org-singularity-bioconductor-deseq2-1.34.0--r41hc247a5b_3.img` (DESeq2 1.34.0, ashr 2.2-54).
   This is the image the published run's DESEQ2_DIFFERENTIAL task used. The script stops if another DESeq2 version
   loads (`expect_deseq2`). It runs 21 fits:
   - `full`, the published contrast, and `full_design_intercept`, a check that `~ Classification` gives the same calls;
   - R1: six leave-one-out fits, 2 v 3;
   - R2: the ten relabellings, of which split 1 is the true labelling;
   - R3: `batch_NL_vs_IL` (NL70B+NL71B against the other four; the IL/NL prefix explained PC1 better than the arm did,
     graft_relation C8) and `IL66B_vs_primaries`.
2. `bin/human_de_null_gsea.R` runs in `go2432-bioconductor.sif` (R 4.3.3, fgsea 1.28.0 in the image). It runs
   fgseaSimple on every fit's Wald statistic, for each collection and for the whole combined file.

Both scripts are copies of `ANALYSIS/host/verify/host_verify_de.R` and `host_verify_gsea.R`, adapted here; the host
originals are not edited. `bin/hv_json.R` is an unchanged copy (md5 6c2c2ff2...). The host scripts reproduced the rat
pipeline table to 2.3e-6.

**Pass A** is the **analysis of record** (pre-registered). It uses the published graft+both counts,
`U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_{counts,lengths}.tsv`. R0 has to pass first, or nothing
else runs.

**Pass B** is the pre-registered **rat-leak sensitivity analysis**. It uses the graft-only counts that rnaseq_graft
writes (`M/rnaseq_graft/results/star_salmon/salmon.merged.gene_{counts,lengths}.tsv`, columns `<lib>_graft`). It runs
the same fits and gene-set tests with pass A's universe rule, and has no R0, since no published table exists for
these counts. It exists because the recurrent arm carries about twice the rat `both` contamination of the primary
arm, so a "recurrence effect" in graft+both counts can be rat leak. Report A and B side by side; B does not replace A.

## The relabellings

IL67B stays in the reference group, so no split also appears as its mirror. In every fit the group without IL67B is
the target, so lfc = target / reference. A relabelling's sign is arbitrary, so ranks use counts and absolute values.

| split | reference (with IL67B) | target | note |
|---|---|---|---|
| 1 | IL67B IL68B IL69B | IL66B NL70B NL71B | **true labelling** (`split_IL66B_NL70B_NL71B`, identical to `full`) |
| 2 | IL67B IL68B IL66B | IL69B NL70B NL71B | |
| 3 | IL67B IL68B NL70B | IL69B IL66B NL71B | |
| 4 | IL67B IL68B NL71B | IL69B IL66B NL70B | |
| 5 | IL67B IL69B IL66B | IL68B NL70B NL71B | |
| 6 | IL67B IL69B NL70B | IL68B IL66B NL71B | |
| 7 | IL67B IL69B NL71B | IL68B IL66B NL70B | |
| 8 | IL67B IL66B NL70B | IL68B IL69B NL71B | |
| 9 | IL67B IL66B NL71B | IL68B IL69B NL70B | |
| 10 | IL67B NL70B NL71B | IL68B IL69B IL66B | |

Ranks are reported as a range when splits tie: rank_min = 1 + (splits strictly above the true split), rank_max =
(splits at or above it). The conservative p is rank_max / 10.

## Settings (every one that is not a tool default, and why)

| setting | value | why |
|---|---|---|
| samples | IL67B IL68B IL69B (Primary_U2), IL66B NL70B NL71B (Recurrent_U2), in `metadata_therapy.csv` order | the published sample sheet; the script refuses a samplesheet whose order differs from it |
| universe | >= 10 raw salmon reads (not rounded) in >= 1 of the **six** tumours | the published sample sheet held only the six, so the pipeline's filter saw six columns. Checked locally on the published tables: this rule gives exactly the 19,351 published ids. Alternatives give 17,186 (>= 2 samples), 16,154 (>= 3), 19,385 (rounded counts) and 20,989 (all ten libraries). The same universe is kept for every fit, so leave-one-out fits are not re-filtered |
| design | `~ 0 + Classification`, salmon gene lengths as `avgTxLength` | the pipeline's DESeq2 module with `transcript_length_matrix` set (published params JSON) |
| DESeq | Wald, parametric, sfType ratio, `minReplicatesForReplace` 99, useT FALSE | `therapy_v3_params.yaml` (99 turns off outlier replacement at n = 3) and the params JSON |
| results | alpha 0.1, independent filtering, BH, minmu 0.5, lfcThreshold 0, greaterAbs | params JSON |
| lfc | `lfcShrink(type = "ashr")`; p and padj from the unshrunk Wald test | the pipeline's `deseq2_shrink_lfc true`; for the rat table R0 found the log2FoldChange column is the ashr one |
| filtered table rule | padj < 0.05 and abs(log2FoldChange) >= log2(1.5), on the table's (ashr) log2FoldChange | derived from the published tables: this rule on the full table gives exactly the 82 filtered rows (58 up, 24 down). `<` against `<=` and `>` against `>=` make no difference there |
| RP block | the members of KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION, read from `combined_human.gmt` | the published "translation" fall is this block of about 80 RPL/RPS genes (memory note of 2026-09-26; `gsea_leave_one_out/SUMMARY.md`). Symbols are read from the file, never typed |
| gene sets | `U/ANALYSIS/refs/pathways/human/combined_human.gmt` split by collection (next section) | the file the published GSEA used |
| fgsea | fgseaSimple, 10,000 permutations, sets of 15-500 genes, seed 1234 before every call | as the host verification; 15-500 are the published `gsea_set_min`/`gsea_set_max` |
| ranking statistic | DESeq2 Wald statistic (unshrunk); one row per symbol, the id with the largest abs(stat) | as the host verification |

### Gene-set collections

The collections are split by set name, using the component files the combined file was built from. Membership comes
from the combined file itself. Counts are from a local read of the grid copies (`combined_human.gmt` md5
6410868154a7688033bbef8b6d1f67fd):

- **combined**: 15,650 names, none repeated. 33 lines carry no gene; fgsea drops them with the other sets outside
  15-500 genes.
- **hallmark**: 50 sets.
- **kegg**: 805 (legacy KEGG and KEGG MEDICUS).
- **reactome**: 1,692 (the REACTOME_ names of hs.c2_curated).
- **wikipathways**: 791 (the WP_ names).
- **c2_other**: 3,945 (BIOCARTA, PID and the chemical and genetic perturbation author sets).
- **gobp**: 7,647.
- **brain**: 878 of hs.brain_gmt_v2's 902.

How many sets of each fall inside 15-500 genes is written to `fgsea_collections.tsv`.

hs.dsigdb is not in the combined file and is not used. Every name in the combined file falls in at least one
collection. Each collection gets its own BH correction; "combined" is the whole file under one BH, as the published
run tested it.

The **pre-registered ribosomal-protein / translation sets** are the six leading sets of the published GSEA:
- KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION
- REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION
- REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY
- KEGG_RIBOSOME
- REACTOME_SELENOAMINO_ACID_METABOLISM
- REACTOME_CELLULAR_RESPONSE_TO_STARVATION

Four more are reported, not ranked as headline: WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS, GOBP_CYTOPLASMIC_TRANSLATION,
REACTOME_TRANSLATION and REACTOME_EUKARYOTIC_TRANSLATION_INITIATION.

## R0: the gate (pass A)

The job stops with exit 10 before any other fit, and nothing from this item is used, unless all of these hold:
- the universe is the same 19,351 gene ids as the published table;
- the padj NA pattern is identical;
- max abs(Δpadj) < 1e-4;
- max abs(Δlog2FC) < 1e-4, comparing the table's log2FoldChange with the ashr estimate;
- the filter rule reproduces the 82 published filtered rows from the published table;
- the fit's own filtered set is the same 82 genes.

`r0_check.json` also records the pvalue, lfcSE and baseMean differences, and which lfc column the table matches. It
also sets our size-factor estimates beside the pipeline's `sizefactors.tsv`; this is reported, not gated, because with
gene lengths DESeq2 uses per-gene normalisation factors.

## Checks against something outside this computation

- **R0.** The published table is the pipeline's own output, produced by the same container on the same counts.
  Reproducing it to 1e-4 shows the fits here are the published fits, not a plausible neighbour.
- **Size factors.** The pipeline's reported size factors are set beside several computations. A local check found
  which one it is: the median-of-ratios of the universe counts after dividing each gene by its lengths scaled to their
  geometric mean. Run locally in R, that reproduces all six published values to 2.4e-15, where plain median-of-ratios
  misses by about 0.03. `r0_check.json` `sizefactor_candidates_max_abs_diff.mor_on_length_adjusted_counts` should
  therefore be at machine precision on the grid too.
- **fgsea against the published Broad GSEA.** `fgsea_vs_broad_combined.tsv` and `summary.json`
  `gsea.vs_broad`. The Broad run used Diff_of_Classes on normalised counts with gene-set permutation, a different
  method, so this checks that the preranked tests see the same top of the list. It gives the Spearman correlation of
  NES, the top-20 overlaps, and the six leading sets side by side.
- **RP block against an independent estimator.** An earlier median-of-ratios computation
  (`SLIDES/07_contamination_magnitude.py`) put the fall of these sets at -0.38 to -0.40 log2 with all six tumours, and
  -0.10 to -0.11 without IL68B. The DESeq2 MLE mean of the block in `full` and `loo_IL68B` (`r1_loo_counts.tsv`,
  `rp_mean_lfc_mle`) should land near those numbers. A large gap means one of the two is measuring something else.
- **Pass B against pass A.** The same code on counts with the rat `both` reads removed. A result present in A and
  gone in B is rat leak.

## Pass B gate

`bin/gate_passb.py` opens the gate only when all of these hold:
1. rnaseq_graft's `NEXTFLOW_EXIT` is 0;
2. `checks.json` is newer than that exit;
3. checks 0, 1, 2, 3, 4 and 10 passed: versions, clip applied once, input counts, reverse strand, GENCODE v44, and no
   empty outputs. These are the checks the counts depend on. Check 9, the saved index, is reported but not required;
4. pass A's R0 passed;
5. both tables carry the six `<lib>_graft` columns.

While the rnaseq_graft check job or pass A is still queued or running, pass B resubmits itself afterany on it, at most
100 times. It does this because the rnaseq_graft chain resubmits its own head and check jobs, so a dependency on an old
id can release pass B too early.

## What the result can and cannot show

**It can show** whether the true labelling is unusual among the ten, on each of these:
- the number of genes at padj < 0.05;
- the number passing the published filter;
- the size of the RP-block fall;
- sets at padj < 0.05, per collection;
- the |NES| of each RP/translation set.

A rank of 1 of 10 is the best possible outcome and means p = 0.10. That is suggestive, and no more.

**It cannot**:
- reach significance at 0.05, by construction;
- separate LITT from anything else that differs between the arms. The arms also differ in tumour content (graft
  fraction, and DNA tumour fraction 0.49-0.66 against 0.71-0.83), in the NL/IL naming batch, and in whatever culture
  batch or passage each arm came from, which is not recorded. The relabelling null asks only whether the grouping is
  special, not why;
- say anything about invasion or the contralateral cells, which are not in these fits.

A rank near the middle of the ten means the published gene list and gene-set result are what any 3-v-3 split of these
tumours produces. A top rank in pass A that disappears in pass B means the result follows the rat leak. Leave-one-out
(R1) shows which single tumour carries each gene and the RP fall. `mean_NES_IL68B_in_target` and `..._in_reference`
in `fgsea_rp_rank.tsv` show whether the RP sets follow IL68B's group rather than the arm.

## Inputs (md5 of the local copies, checked equal to the grid copies on 2026-09-27)

| file | md5 |
|---|---|
| salmon.merged.gene_counts.tsv (published human) | a8d551fcfdb22e737f7c667602842c09 |
| salmon.merged.gene_lengths.tsv (published human) | 0c60637e5e3ded1a6bb912f5aaf3af4d |
| therapy_impact.deseq2.results.tsv | 1823dd29dc96bf3157830a2964e30a39 |
| therapy_impact.deseq2.results_filtered.tsv | 21dd8e7cb4bb12de0967b3923a6a4869 |
| therapy_impact.deseq2.sizefactors.tsv | 6c31a109f273bfd8b74ce60c38599b7a |
| metadata_therapy.csv | 6059de94e2ed0168845bc1efa6b089d2 |
| params_2026-05-24_19-23-04.json | 0c4faa5a7e35be57a9d24361c80c74ec |
| combined_human.gmt | 6410868154a7688033bbef8b6d1f67fd |

Each job also writes `inputs_md5.tsv`, with every input, script and GMT hashed at run time, and `containers.tsv`, with
the sha256 of both images.

## Outputs

`M/human_de_null/results/pass{A,B}/`. The job also links `M/de_null` to `M/human_de_null`, so the spec's path works.

| file | what |
|---|---|
| `r0_check.json` | R0 (pass A) |
| `fits/<fit>.tsv.gz`, `fits.tsv` | every fit: MLE and ashr lfc, Wald stat, p, padj; manifest |
| `ranks/<fit>.tsv` | each fit's Wald statistic (fgsea input) |
| `splits_summary.tsv` | one row per relabelling. DE counts, RP block, sets at padj < 0.05 per collection, NES of the six RP sets |
| `true_split_rank.tsv` | the true split's rank range on every metric (DE, gene sets, RP sets) |
| `r1_loo_counts.tsv`, `r1_loo_genes.tsv` | leave-one-out: counts, and per gene which tumour it needs |
| `r2_split_top_genes.tsv` | what drives each relabelling |
| `r3_counts.tsv`, `r3_il66b_and_nl_genes.tsv` | NL pair and IL66B alone |
| `rp_block_by_sample.tsv` | RP block level per tumour against the six-tumour mean |
| `fgsea_counts.tsv`, `fgsea_relabel_null.tsv`, `fgsea_rp_sets.tsv`, `fgsea_rp_rank.tsv`, `fgsea_top_full.tsv`, `fgsea_top_by_fit.tsv`, `fgsea_full_<coll>.tsv`, `fgsea_collections.tsv`, `fgsea_vs_broad_combined.tsv` | gene sets |
| `summary.json` | verdict lines, every rank, input md5s, container sha256s, R0, the DE and gene-set summaries |
| `STATUS` | running / done / r0_failed / gate_closed / deferred_to_<id> / ... |
| `../side_by_side_splits.tsv`, `../side_by_side_ranks.tsv` | written by pass B: A and B beside each other |

Logs: `/rs/rs_grp_oschome/go2432/u251_multiomics/logs/human_de_null_u251_de_null_{A,B}_<jobid>.out`.

## Resources and timing

Each pass is reqp `--qos=requeue --requeue`, 4 CPU, 24 GB, 3 h limit. The host run of the same code took about 10
minutes (job 40475691). The human side has more genes and about three times the gene sets, so expect **30-60 min per
pass**; most of that is converting the two images into a sandbox and running fgsea on 8 collections × 20 fits.

- Pass A starts when reqp has room.
- Pass B waits for the rnaseq_graft checks: prep_reads (1.5-3 h) → rnaseq 3.22.2 (4-6 h) → checks, then about 45 min.

No reference is downloaded or built. There is no Nextflow and no `-profile test` run.

## Runs

Submitted 2026-09-27 01:43 EDT with `submit.sh`:
- **Pass A: job 40481672.** It started 01:54 on rgd5.
  - The DE step finished at 01:57. R0 **passed**: 19,351 identical ids, padj NA pattern identical, max abs(Δpadj)
    4.1e-5, max abs(Δlog2FC) against ashr 3.2e-5, and 82 of 82 filtered genes. The margin on padj is only about 2.4
    times; the rat verification reached 2.6e-6.
  - Loaded versions: DESeq2 1.34.0, ashr 2.2.54, libPaths `/usr/local/lib/R/library` only.
  - Image sha256: deseq2 `425e77b97dbe465ad4b0b399c0b33b7c02c7ee60b082e32fd8a848e86832bfa3`; go2432-bioconductor.sif
    `69485ffa6dfa7cd884ead560227e05725fce05cc7e5ca28faa024ea7c4ae9fe0`.
  - fgsea (1.28.0 from the image) finished about 02:02; STATUS `done`. Its only warnings were rank ties
    (0.18-0.31 % of genes).
  - Size factors: our length-adjusted median-of-ratios matches the pipeline's file to 2.4e-15.
  - The largest padj differences sit on one Benjamini-Hochberg plateau at padj 0.674 (37 genes above 1e-5), driven by
    p-value differences of at most 3.8e-5. Why the p-values differ at that level, when the rat run matched to 2.3e-6,
    is not known. A different number of worker cores (4 here, 8 in the pipeline) is a guess, not a tested cause.
  - **Pass A result (record, conservative p = rank_max / 10).**
    - Genes: the true split ranks 1 of 10 on padj < 0.05 (107; next 70, median of the other nine 16) and on the
      published filter (82; next 49, median 14).
    - RP block: 1 of 10 on |mean lfc| (0.42; next 0.38).
    - Gene sets at padj < 0.05, by collection: combined 1 (1,630 against a median of 1,190), gobp 1, brain 1,
      c2_other 1, kegg 2, wikipathways 2, hallmark 2-3, reactome 6.
    - All six leading RP/translation sets rank 1 of 10 on |NES|. But KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION
      reaches padj < 0.05 in its own collection in **all ten** splits. Its sign follows IL68B's group: mean NES +3.13
      when IL68B is in the target, -3.32 when it is in the reference. The RP result is the strongest version of a
      signal that every relabelling produces, not one specific to the arm.
    - Every rank-1 result is p = 0.10, the floor.
    - Without IL68B the RP fall shrinks from -0.42 to -0.15. Leave-one-out keeps 40-76 of the 82 filtered genes, the
      fewest without NL70B.
- **Pass B: job 40481673**, afterany on 40481672 and on the rnaseq_graft check job 40480726. That check job is afterany
  on head 40480725, which is afterany on collect_prep 40480710, which is afterany on the prep array 40480709. If the
  chain resubmits its jobs, pass B follows the new ids itself (see Pass B gate). Expected about 45 min after the
  rnaseq_graft checks, which is roughly 7-10 h after 01:43 plus queue time.

**Do not edit `bin/*.R` on the grid while a pass is running.** Rscript reads the file as it executes. In the local
dry run, an edit made mid-run shifted what the parser read, and the run died on an error that is not in the code. To
change a script, wait until STATUS is final, or copy it to a new name.

A local dry run was done before submission, in scratchpad `multiomics/de_null_probe/`. It ran R 4.5.2, DESeq2 1.50.2,
ashr 2.2.63 and fgsea 1.36.2 with 1,000 permutations, so its numbers are not the record. It exercised both passes, the
R0 stop path (exit 10), the pass B gate (six cases) and the side-by-side tables. Two of its numbers were checked
against something outside the code:
- The RP block's MLE mean was -0.42 with all six tumours, -0.15 without IL68B, and -0.40 to -0.56 without each of
  the others. The earlier median-of-ratios estimate was -0.38 to -0.40, -0.10 to -0.11 and -0.37 to -0.53.
- 8,869 combined sets fell in 15-500, the same count as the published Broad run.

## Files

| file | role |
|---|---|
| `run.sbatch` | the job (both passes) |
| `submit.sh` | submits pass A now, and pass B afterany on pass A and the rnaseq_graft check job |
| `params.yaml` | every path and setting the job reads |
| `samplesheet.csv` | the six tumours, with their count column in each pass |
| `grid.config` | container and resource record, in the same shape as the other items (no Nextflow reads it) |
| `bin/human_de_null_de.R`, `bin/human_de_null_gsea.R`, `bin/hv_json.R`, `bin/gate_passb.py` | the code |
