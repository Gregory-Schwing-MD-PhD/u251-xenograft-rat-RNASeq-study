# rnasplice: splicing after LITT, calibrated by the 10-split relabelling null

Part of the U251N LITT multi-omics build (SELECTION.md, item S9). This folder holds only configuration, the sample
sheet (sample names and grid paths), the check code and the null code. No reads, BAMs or results are here; they live on
the grid under `/rs/rs_grp_oschome/go2432/u251_multiomics/` (`M` below), in `M/rnasplice/`.

## Question

Does recurrence after LITT change splicing in the U251N tumour cells? Three readouts, recurrent (IL66B, NL70B, NL71B)
against primary (IL67B, IL68B, IL69B), 3 v 3, unpaired:

- **exon usage** (DEU): DEXSeq on exon bins;
- **transcript usage** (DTU): DRIMSeq filter, DEXSeq, stageR confirmation;
- **splicing events**: rMATS (SE, MXE, A3SS, A5SS, RI) and SUPPA2 (local events and isoforms).

The pipeline's nominal FDR is **not trusted and not reported on its own**. Every count is ranked against the 10 balanced
relabellings of the same six tumours, computed outside the pipeline with the pipeline's own scripts and containers
(`null/`). With 3 v 3 there are exactly 10 such splits, so the smallest attainable p is 0.10.

**Expected result: nothing reportable.** At 4 residual degrees of freedom, 80 % power at alpha 0.001 needs a
difference of about 7.8 within-arm SDs; duplicates (65-79 % in these tumours, not removed by rnasplice) make rMATS's
nominal FDR optimistic. The host DE (38 genes at FDR < 0.05) already ranked 2 of 10 under the same kind of null.

## Pre-declared before any result existed (2026-09-27, before submission)

1. **Analysis of record: STAR + Salmon** (`results/star_salmon/`). `pseudo_aligner: null` is set so that the salmon
   pseudo-alignment branch does not run. A local `-preview` of 1.0.4 on Nextflow 23.10.0 with this exact params.yaml
   (paths swapped to local fixtures) showed nf-validation 1.1.3 accepts the null (it drops falsy params before
   validating) and the DAG holds no SALMON_INDEX, SALMON_QUANT_SALMON or `*_SALMON` DTU/SUPPA branch. If the grid run
   nevertheless produces `results/salmon/`, it is never used for a claim, and the null uses only `star_salmon`.
2. **Primary count per method** (secondary counts are reported beside them and never promoted afterwards):

   | method | primary count | secondary |
   |---|---|---|
   | DEU | genes with DEXSeq `perGeneQValue` < 0.05 | exon bins with padj < 0.05 |
   | DTU | genes with DEXSeq (DTU) `perGeneQValue` < 0.05 | transcripts confirmed by stageR at 0.05 |
   | rMATS | events with FDR < 0.05 and \|IncLevelDifference\| >= 0.1 in the **JC** files, summed over the 5 types | same in JCEC |
   | SUPPA2 | local events (ioe) with gene-corrected p < 0.05 | same with \|dPSI\| >= 0.1; isoforms (ioi) p < 0.05 |

   JC, not JCEC, for rMATS: this is total RNA (rRNA-depleted), and the exon- and intron-body reads JCEC adds include
   pre-mRNA, which moves with transcription, not splicing.
3. **Reportable** means: the true split's count ranks 1 of 10 (p = 0.10) **and** task 1 of the null reproduces the
   pipeline's table for that method. Anything else is "not distinguishable from an arbitrary regrouping".

## Inputs (grid)

| sample | condition | reads |
|---|---|---|
| IL67B_graft, IL68B_graft, IL69B_graft | primary | `M/reads/graft_clip/<lib>_graft_R{1,2}.fq.gz` |
| IL66B_graft, NL70B_graft, NL71B_graft | recurrent | same |

The graft-only records of each human file (xengsort `graft` bin; the rat-ambiguous `both` bin is excluded because the
recurrent arm carries about twice the rat `both` leak of the primary arm), with the Takara SMARTer Pico v2 3-nt R2
artefact clipped **once, upstream**, by prep_reads (gate A found the artefact in all ten libraries). Strandedness
`reverse` (ISR). **Six tumours only**: the DRIMSeq filter counts samples across the whole sheet, so the near-empty
controls or C2B would change which transcripts are tested.

Reference: GENCODE v44 (Ensembl 110) GRCh38 primary assembly, unzipped, md5-checked against GENCODE's MD5SUMS and
faidx'd by stage_refs (`M/refs/gencode_v44/`, marker `M/refs/.done/gencode`). The STAR index is built inside the run.

## Release

nf-core/rnasplice **1.0.4** (`-r 1.0.4`; GitHub refs/tags/1.0.4 = commit `1d0494ae3402d1a46e0adadad24f81a0ff855c77`,
checked by the head job after the run), Nextflow **23.10.0** via `NXF_VER`. The manifest allows `!>=23.04.0`, and its
plugin nf-validation 1.1.3 predates Nextflow 25, so a Nextflow contemporary with the release is used. Java 17 from
the grid's nextflow env runs it. The head job refuses to start if `nextflow -version` does not report 23.10.0.
Tools inside the release: STAR 2.7.10a, salmon 1.10.1, TrimGalore 0.6.7, HTSeq 2.0.2 (dexseq_count), DEXSeq 1.36.0,
DRIMSeq 1.18.0, stageR 1.12.0, rMATS 4.1.2, SUPPA 2.3.

## Every non-default setting and why

| parameter | value | default | why |
|---|---|---|---|
| `pseudo_aligner` | null | salmon | one result set (see pre-declaration 1) |
| `gencode` | true | false | GENCODE annotation; inert here (it only changes the salmon index, which is not built) |
| `igenomes_ignore` | true | false | never fall back to an iGenomes reference |
| `rmats_read_len` | 100 | 40 | usage.md: give the median for variable-length reads (R1 35-101, clipped R2 32-98); `--variable-read-length` is always on |
| `rmats_paired_stats` | false | true | the default pairs tumour i of one arm with tumour i of the other (PAIRADISE); the design is unpaired |
| `edger_exon` | false | true | same counts, another model: not independent evidence, one more count to choose from |
| `save_dexseq_annotation` | true | false | publishes the flattened GFF, which the null re-uses |
| `min_samps_gene_expr` | 6 | 4 (config) / 6 (schema) | the config and the schema disagree, so all six filter values are set |
| `min_samps_feature_expr`, `min_samps_feature_prop` | 3, 3 | 2 (config) / 0 (schema) | 3 = the smaller group (Love, Soneson & Patro 2018, F1000Res 7:952) |
| `min_gene_expr`, `min_feature_expr`, `min_feature_prop` | 10, 10, 0.1 | same | written out |
| `diffsplice_paired` | false | true | the default adds `-pa` (paired replicates); the design is unpaired |
| `diffsplice_method`, `diffsplice_gene_correction` | empirical, true | same | written out |
| `sashimi_plot` | false | true | MISO plots only |
| `max_cpus`, `max_memory`, `max_time` | 40, 240.GB, 96.h | 16, 128.GB, 240.h | caps for reqp nodes (Nextflow 23.10 has no `resourceLimits`) |
| `clip_r2` | **not set** | null | TrimGalore's ext.args hold only `--fastqc`/`--nextseq` (conf/modules.config), so it is not wired to anything; the clip was applied upstream |

`grid.config` (loaded with `-c`): slurm on `reqp`, `--qos=requeue`, retry on kill/pre-emption/OOM codes (3 retries),
`STAR_GENOMEGENERATE` 80 GB (process_high is 60 GB; the published human genomeGenerate peaked at 64.8 GB);
`DEXSEQ_EXON`/`DEXSEQ_DTU` 2 CPUs, 64 GB, 48 h (the R scripts use DEXSeq's default SerialParam, so process_high's 12
cores would sit idle); `DEXSEQ_COUNT` and `DIFFSPLICE_*` 2 CPUs (single-threaded tools). Resource settings only; no
analysis setting is changed by the config. Singularity cache: `M/singularity`, shared by the nf-core pipelines of this
round (rnafusion uses it too; the earlier items used `/rs/rs_grp_oschome/go2432/singularity_cache`, which holds only 4
of the 31 depot.galaxyproject.org images rnasplice 1.0.4's modules name: fastqc 0.12.1, python 3.9, ubuntu 20.04,
ucsc-bedclip). The strand check reads the samtools 1.21 image from that older cache. Nextflow 23.10 gets its
own `NXF_HOME` (`M/nxf_home_nf23`) so its plugin and asset checkout never share a directory with the Nextflow 25 runs.

## Jobs and gate

`submit.sh` submits three jobs: the head (`run.sbatch`, afterany on collect_prep while it is queued), the checks
(`check_rnasplice.sbatch`, afterany head) and the null launcher (`null/null_launch.sbatch`, afterany checks).

- **Gate** (head): Nextflow starts only when `M/reads/READY` **and** `M/refs/.done/gencode` exist. A missing marker
  whose producer is still queued makes the head resubmit itself afterany on it (at most 72 times); a producer that
  finished without its marker closes the gate (`NEXTFLOW_EXIT=gate_closed`).
- The head is **not** chained on stage_refs: its gencode step had already passed at submission
  (`.done/gencode`, 00:46), and the steps still running (known sites, VEP, the rnafusion bundle) are not used here.
  Waiting for them would add 3-6 h for nothing; the gate still requires the gencode marker.
- `-resume` is always on (reqp is PreemptMode=REQUEUE; a requeued head restarts from the top). One launch directory.

## Checks (`check_rnasplice.py`, writes `CHECKS.md` and `checks.json`)

Hard checks; the null launcher refuses to run unless all pass.

| # | check | against what |
|---|---|---|
| 0 | Nextflow exit 0; Nextflow 23.10.0; commit = tag 1.0.4; the params actually applied (`pipeline_info/params_*.json`) equal the intended values | the release tag on GitHub; params.yaml |
| 1 | every expected process ran for all six samples; no edgeR/MISO/bigWig; salmon pseudo-alignment branch absent (reported, not failed, with the pre-declared rule if present) | the DAG of the local preview |
| 2 | task scripts: rMATS prep and post `--libType fr-firststrand --readLength 100 --variable-read-length -t paired`, no `--paired-stats`, no `--novelSS`; dexseq_count `-s reverse -p yes -r pos -a 10`; diffSplice `-m empirical -gc`, no `-pa`; salmon `--libType=ISR`; TrimGalore without any clip; DRIMSeq filter called with `6 3 3 10 0.1 10`; STAR two-pass with TranscriptomeSAM | the intended settings, read from what actually ran |
| 3 | TrimGalore saw exactly N_graft pairs per library; STAR's input = pairs TrimGalore kept | the xengsort counts (prep_reads check b) |
| 4 | FastQC on the input: R1 max 101, R2 max 98 | the upstream clip: 101 - 3 = 98, applied once |
| 5 | salmon expected format ISR and compatible-fragment ratio; **independently**, R1 orientation on the exons of 20 named protein-coding genes (strand and exons from the GTF), samtools 1.21 in its image, MAPQ >= 30 | the published run's RSeQC antisense share, 97.4-98.4 % |
| 6 | six samples in every per-sample output; the bamlists hold the right three BAMs each; the rMATS prep lists every BAM exactly once; all 10 MATS files exist | the sample sheet |
| 7 | every file and container the null re-uses exists; `pipeline_bin/` md5 equals the 1.0.4 tarball's scripts; the exact diffSplice and rMATS option strings are taken from the pipeline's own task scripts | the GitHub 1.0.4 tarball |

Check 5 has a WARN band: a strand error reads near 0.5 (unstranded) or near 0.03 (reversed), so below 0.90 fails and
blocks; 0.90-0.97 is short of the published 97.4-98.4 % and is flagged for a look without blocking the null.

## The null (`null/`)

| task | label | recurrent (treatment) | primary (control) |
|--:|---|---|---|
| 1 | split01 (true) | IL66B NL70B NL71B | IL67B IL68B IL69B |
| 2 | split02 | IL69B NL70B NL71B | IL67B IL68B IL66B |
| 3 | split03 | IL69B IL66B NL71B | IL67B IL68B NL70B |
| 4 | split04 | IL69B IL66B NL70B | IL67B IL68B NL71B |
| 5 | split05 | IL68B NL70B NL71B | IL67B IL69B IL66B |
| 6 | split06 | IL68B IL66B NL71B | IL67B IL69B NL70B |
| 7 | split07 | IL68B IL66B NL70B | IL67B IL69B NL71B |
| 8 | split08 | IL68B IL69B NL71B | IL67B IL66B NL70B |
| 9 | split09 | IL68B IL69B NL70B | IL67B IL66B NL71B |
| 10 | split10 | IL68B IL69B IL66B | IL67B NL70B NL71B |
| 11-16 | loo_<lib> | true labels with one tumour dropped (IL66B, NL70B, NL71B, IL67B, IL68B, IL69B) | |

IL67B is always in the control group, so no split appears with its mirror. The group without IL67B is always labelled
`recurrent` and the other `primary`, so file names, contrast names and code paths are those of the pipeline run.

What each task re-runs (the pipeline's own scripts from `pipeline_bin/`, md5-checked; the pipeline's own containers,
read from the tasks' `.command.run`; the pipeline's published intermediates):

- **DEU**: `run_dexseq_exon.R` on `dexseq_exon/counts/<sample>.clean.count.txt` and `DEXSeq.gff`.
- **DTU**: the DRIMSeq-filtered counts rebuilt from the published `dmDSdata.rds` exactly as `run_drimseq_filter.R`
  writes them (only the .rds is published), then `run_dexseq_dtu.R` and `run_stager.R`. The DRIMSeq filter counts
  samples, never groups, so one filter serves every split; the leave-one-out rows reuse the six-sample filter.
- **SUPPA2**: `suppa_split_file.R` on the published TPM and PSI files with the relabelled sheet, then `suppa.py
  diffSplice` with the option string taken from the pipeline's DIFFSPLICE_IOE task (no `-pa`), on events.ioe and events.ioi.
- **rMATS**: `rmats.py --task post` on the pipeline's prep (`.rmats`) with regrouped `--b1/--b2` and the option string
  taken from the pipeline's RMATS_POST task. rMATS's own `--task inte` first confirms every BAM is named once (the 10
  splits). Post, not the cheaper `--task stat` on the published counts: rMATS's count filter (`filter_countfile`) keeps
  an event only if each group has reads, and it rewrites the published `*.raw.input.*` files in place under the true
  grouping, so regrouping those files would test the true split's event set in every split. Post re-applies the filter
  per grouping. Fallback, only if post fails: `--task both` on the published BAMs with `--libType fr-firststrand`
  (the route is recorded per task). **Never** `--source genome_bam`/`salmon_results`: `create_genome_bam_channel`
  drops strandedness, so rMATS would run fr-unstranded and dexseq_count would count sense.

**Task 1 is the null's own check.** It re-derives every table the pipeline published for the true split; null_collect
requires identical keys and |difference| <= 1e-6 (DEXSeq gene and bin padj, DTU gene padj, stageR, both SUPPA dpsi
files, all 10 MATS files). A method whose task 1 does not reproduce the pipeline gets no rank.

Two arrays of 16 (DEU single-threaded; DTU+SUPPA+rMATS on 4 CPUs), sized by `null_size.py` from the pipeline's own
trace (twice the pipeline's DEXSEQ_EXON realtime, 1.5 x its peak RSS, and so on), then `null_collect`
(afterany both). Each method writes `runs/<label>/<method>/DONE`; a requeued task skips finished methods.

Outputs: `null/results/NULL_REPORT.md`, `null_ranks.tsv` (per count: true value, the other nine, rank range, p range,
reproduced?), `null_loo.tsv` (per dropped tumour: the count and how many of the true split's hits survive),
`null_metrics.tsv`, `null_summary.json`. `null_collect` deletes `M/rnasplice/work` only when Nextflow exited 0, every
hard check passed, all 16 tasks finished and task 1 reproduced the pipeline.

## Tested locally before submission (WSL, nothing on the grid)

- `nextflow run nf-core/rnasplice -r 1.0.4 -preview` on Nextflow 23.10.0 with this params.yaml and grid.config (paths
  swapped): validation passed, `pseudo_aligner: null` applied, `params_*.json` shows every intended value, and the DAG
  holds one DTU/SUPPA branch (STAR_SALMON) and no salmon index or pseudo-alignment.
- null_task.sbatch (16 tasks x both method sets) and null_collect.py with a fake `singularity` that writes the tools'
  output files as a deterministic function of the grouping: every task finished, task 1 reproduced the fake pipeline
  tables, split01 used the pipeline's bamlist order, the rMATS fallback route ran when post was made to fail, the
  leave-one-out sheets held five samples, a rerun skipped finished methods, ranks and ties came out as computed by hand.
- check_rnasplice.py on a fake results/work tree shaped like 1.0.4 output (good, `--paired-stats`, sense strand,
  unclipped R2, unstranded reads).

## What the result can and cannot show

- **Can**: whether the number of splicing changes between the arms is larger than for any regrouping of the same six
  tumours (at best p = 0.10), and whether that number leans on one tumour (leave-one-out).
- **Cannot**: significance at 0.05 for anything; a claim about a single event unless its own statistic also ranks
  1 of 10; separate LITT from what else differs between the arms (tumour content, the IL/NL naming batch, the rats'
  lesion sizes, culture passage, which is not recorded).
- Duplicates are not removed (rnasplice has no duplicate marking), so every nominal FDR is optimistic; every split
  carries the same duplicates, which is why the rank, not the FDR, is the calibration.

## Not known

- Wall time on reqp under the current load; the only measured baseline is the published rnaseq trace (STAR index
  33 min, 65 GB).
- How long DEXSeq DEU takes single-threaded on the GENCODE v44 flattened annotation (hours; the null sizes itself from
  the pipeline's own run).
- Whether rMATS post on the published `.rmats` reproduces the pipeline's MATS files bit for bit; task 1 tests exactly
  this, and the fallback exists if post fails outright.
- Whether SUPPA2's empirical method behaves with two replicates in a group (leave-one-out rows only).

## Estimated cost

Pipeline: 10-18 h wall after the gate opens (STAR index about 1 h, alignment about 1 h, dexseq_count 2-5 h, DEXSeq
DEU 3-10 h single-threaded on the critical path), about 150-300 CPU-h. Null: 4-12 h wall, about 130-320 CPU-h
allocated (16 DEU tasks on 1 CPU each; 16 tasks on 4 CPUs). Results about 60 GB; work about 200 GB, deleted as above.
Queue waits on reqp are not included.
