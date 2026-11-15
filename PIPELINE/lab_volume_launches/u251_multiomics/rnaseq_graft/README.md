# rnaseq_graft: nf-core/rnaseq 3.22.2 on graft-only U251 reads

Part of the U251N LITT multi-omics build (SELECTION.md, item S3). This folder holds only configuration, the sample
sheet (sample names and grid paths) and the check code. No reads, BAMs or results are here; they live on the grid under
`/rs/rs_grp_oschome/go2432/u251_multiomics/` (`M` below). `U` is the grid clone
`/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study`.

## Questions

1. **What are the U251 profiles with the rat `both` leak removed?** The published human counts
   (`U/ANALYSIS/results_human_final`) were made from xengsort's `graft` + `both` bins. `both` reads match rat and human
   equally; they are 10.9-24.6 % of each tumour's human file (recurrent 14.6-24.6 %, primary 10.9-14.7 %), so they can turn rat
   contamination into an apparent recurrence effect. This run uses the `graft` bin only. It feeds the invasion and
   microenvironment readouts (S5 parts B, C), a graft-only sensitivity DE (S4 pass B), and later arms.
2. **How much does the Takara 3-nt Read 2 clip change the counts?** SMARTer Stranded Total RNA-Seq v2 Pico (634411):
   the first 3 nt of Read 2 come from the template-switching oligo. Neither published run clipped them (SCOUT s3).
3. **What are IL64B's 314 k (and N168B's 435 k) "graft" reads?** Rat mis-sorted as human, non-U251 human, or U251 cells?
4. **A count-based estimate of the rat leak** in the published tumour counts, independent of the xengsort rate model.
5. **Save the human STAR and salmon indices**, which the published run did not keep.

## Inputs (grid)

| sample | reads | what they are |
|---|---|---|
| `<lib>_graft` x 10 | `M/reads/graft_clip/<lib>_graft_R{1,2}.fq.gz` | the first N_graft records of each human file (graft bin), R2 clipped by 3 nt upstream by prep_reads |
| `IL69B_fullclip` | R1 `U/ANALYSIS/sorted_fastqs/IL69B_human_R1.fq.gz`, R2 `M/reads/human_clip/IL69B_human_R2.fq.gz` | exactly the records of the published IL69B (graft + both), R2 now clipped |

Libraries: C2B (culture), IL64B (primary, no tumour in the sample), IL66B, IL67B, IL68B, IL69B, N168B and N269B
(contralateral hemispheres of rats 68 and 69), NL70B, NL71B. Authoritative key: `ANALYSIS/SAMPLE_KEY.md`.

Reference: GENCODE v44 (Ensembl 110) GRCh38 primary assembly, the same two files the published run used
(`U/ANALYSIS/refs/human/GRCh38.primary_assembly.{genome.fa,annotation.gtf}.gz`). The indices are built inside the
run and saved to `M/rnaseq_graft/results/genome/index/`.

## Release

nf-core/rnaseq **3.22.2** (`-r 3.22.2`), Nextflow **25.10.2** (`NXF_VER`; the manifest requires `!>=25.04.0`). This is
the release, Nextflow and tool set of both published runs (v3.22.2-g3816d48: STAR 2.7.11b, salmon 1.10.3, TrimGalore
0.6.10 / cutadapt 4.9), so only the input reads and the clip differ. 3.27.0 was rejected: it moves to salmon 2.7.0,
which changes gene-level quantification by its own changelog.

## Settings that are not pipeline defaults, and why

`params.yaml`

| setting | value | why |
|---|---|---|
| `input`, `outdir` | `M/rnaseq_graft/{samplesheet.csv,results}` | CephFS; home is at 97 % of its quota |
| `fasta`, `gtf` | GENCODE v44 primary | identical reference to the published run, so tables join gene for gene |
| `featurecounts_group_type` | `gene_type` | GENCODE GTFs carry `gene_type`, not `gene_biotype` (as published) |
| `save_reference` | `true` | keep the human STAR and salmon indices this time |
| `aligner`, `skip_bbsplit`, `extra_trimgalore_args` | `star_salmon`, `true`, `null` | defaults, written out. `extra_trimgalore_args` stays null because the clip was applied upstream; `--clip_r2 3` here would remove 3 template bases a second time |

Left at defaults, as in the published run: `gencode` (false), `min_mapped_reads` (5), markduplicates on, trimmer
TrimGalore with adapter auto-detection, `min_trimmed_reads` (10000). `min_mapped_reads` removes a sample only from
genome-BAM steps (`workflows/rnaseq/main.nf` l.384-430); salmon quantification runs before it, so IL64B and N168B keep
their salmon counts even if STAR maps under 5 % of their reads.

`samplesheet.csv`

- `strandedness` is `reverse` for every row, not `auto`, so the control libraries' few, possibly non-U251 reads
  cannot drive the inference. The pipeline's own strand check then becomes a test (check 3).
- Sample names carry `_graft` / `_fullclip`, so any script that joins these tables to the published ones by the
  bare library name fails instead of returning plausible numbers.

`grid.config` (replaces `U/nextflow.config`, which sends tasks to `-q primary` and sets a xengsort default container)

- Tasks on `reqp` with `--qos=requeue`. reqp is `PreemptMode=REQUEUE`, so exit codes 130-145, 104 and 175-177
  (kill, pre-emption, out of memory) are retried up to 3 times; anything else finishes the run.
- `resourceLimits` 40 CPU / 240 GB / 96 h (largest reqp node types; the old `max_*` params are not used).
- `STAR_GENOMEGENERATE` 80 GB (published peak 64.8 GB) and `PICARD_MARKDUPLICATES` 100 GB (published peak 81.8 GB).
- Singularity with the shared CephFS image cache `/rs/rs_grp_oschome/go2432/singularity_cache`, which already holds the
  3.22.2 images from the host run. Every container gets `M/tmp` as `/tmp`, because node-local `/tmp` is swept
  (FastQC died there, job 40468032).

`run.sbatch` (head job: reqp requeue, 2 CPU, 12 GB, 4 days)

- `NXF_HOME`, `NXF_TEMP`, `TMPDIR` and the Singularity cache and temp directories all on CephFS under `M`.
- One launch directory for this run only (a shared one caused a resume-lock clash, job 40467893). `-resume` always,
  because a requeued head restarts from the top.
- **Gate.** Nextflow starts only if `M/reads/READY` exists; collect_prep writes it only when all ten prep checks
  passed (boundary, pairing, byte identity, clip shape). If prep has not finished, the head does not fail and does not
  sit on a node: it resubmits itself (afterany on collect_prep if that job is queued, else a 30-min `--begin`), at most
  72 times, and records the live id in `M/rnaseq_graft/HEAD_JOB`. If prep finished without READY, the gate closes
  (`NEXTFLOW_EXIT` = `gate_closed`). It also refuses to start if any FASTQ in the sample sheet is missing or empty.
- No `-profile test` run first: this exact release ran twice on this grid (human May, host September).

## Checks (`check_rnaseq_graft.py`, run by `check_rnaseq_graft.sbatch` afterany on the head)

Writes `M/rnaseq_graft/CHECKS.md`, `checks.json` and `tables/`. S4 pass B and invasion_niche read `all_hard_pass`.
Hard checks must all pass; the reports carry numbers and verdicts. The references each check uses are outside this
computation where one exists.

| # | check | kind | compared against |
|--:|---|---|---|
| 0 | same versions (pipeline commit, STAR, salmon, TrimGalore, cutadapt, FastQC, RSeQC); `extra_trimgalore_args` null | hard | the published run's versions file |
| 1 | clip exactly once: raw R1 max 101 nt, raw R2 max 98 nt in all 11; no TrimGalore report shows a 5' clip | hard | the 35-101 nt length range of the delivered reads |
| 2 | FastQC and TrimGalore read counts = N_graft (N_graft + N_both for IL69B_fullclip); STAR input = pairs TrimGalore kept | hard | xengsort classify logs (`U/ANALYSIS/xengsort_out`), re-parsed |
| 3 | RSeQC antisense share of assignable reads >= 0.95 for C2B and the six tumours; strand table not failed | hard | published run 97.4-98.4 %; Takara kit is reverse-stranded |
| 4 | a result GTF carries the "version 44 (Ensembl 110)" header; gene ids identical to the published table; 11 samples | hard | GENCODE's own header; published table |
| 5 | clip effect: IL69B_fullclip vs published IL69B, genes at CPM > 10 in both | report + decision | published IL69B |
| 6 | both-bin effect: IL69B_graft vs IL69B_fullclip; RP-block share by library, published vs graft | report | - |
| 7 | rat leak from counts (below) | report | prior estimate: primary 1.9-5.1 %, recurrent 5.1-8.0 % |
| 8 | IL64B and N168B graft reads: STAR metrics, top 50 genes, MT-/rRNA/histone/conserved small RNA shares, chrY and XIST; verdict | report | C2B_graft (pure U251) |
| 9 | STAR index (Genome, SA, SAindex, genomeParameters.txt, ...) and salmon index (full salmon 1.10 file set) saved, and identical file by file to the task's own output in `work/` | hard | the task output; the rat run published only 2 salmon files |
| 10 | `quant.genes.sf` has > 10,000 genes with NumReads > 0 for the 8 libraries with >= 3.8 M pairs; no trimmed FASTQ empty | hard | the empty-TRIMGALORE failure of job 40468266 |

Fixed before the run:

- **Clip decision (check 5).** The clip is *material* if more than 1 % of genes at CPM > 10 in both move by
  |log2 ratio| > 0.2. Material triggers the host rerun (SELECTION s4). Run-to-run differences are folded in, so this
  is an upper bound on the clip effect. Also reported: Spearman, the shares beyond 0.1 and 0.2, the ribosomal-protein
  block (protein-coding cytosolic RPL*/RPS*, no pseudogenes or RPS6K kinases), and mean log2(CPM+1) of the Neftel AC
  set (and MES, NPC, OPC for context). The Neftel sets are `genesets/neftel_*.txt`, copied from
  `ANALYSIS/spread/signatures/` (derived from the repository's Neftel CIBERSORT signature, not the published
  meta-module lists; see that folder's MANIFEST.md).
- **Rat leak (check 7).** C2B is pure human, so r_g = count(C2B_graft) / count(published C2B) is the share of gene g's
  human reads that xengsort keeps as graft; classification depends only on the read sequence. For genes with
  r_g >= 0.2 and published C2B >= 50, the rat part of tumour t's published count is P_gt - G_gt / r_g. Reported per
  tumour as a share of the eligible genes' published counts, signed and with genes clipped at 0.
- **Control verdict (check 8).** With U = STAR unique %, Q = share of counts in MT-, rRNA, histone and conserved
  small-RNA genes, Y = non-PAR chrY reads, Y_exp = C2B_graft's chrY CPM x the library's counts,
  rho = Spearman with C2B_graft: *rat mis-sorted as graft* if (U < 0.5 U_C2B and too short + too many mismatches > 20 %)
  or Q > max(0.5, 2 Q_C2B); *non-U251 human* if U >= 0.8 U_C2B, Q < max(0.3, 1.5 Q_C2B), >= 3,000 genes with >= 5
  counts and Y <= max(4, 0.2 Y_exp) (then every floor based on the controls must be revisited); *U251 human (tumour
  cells in the control)* if the same U and Q bounds hold, Y >= 0.5 Y_exp and rho >= 0.7 (a fourth outcome, added
  because U251 is male and the rats female); otherwise *unclassified* (triggers the competitive alignment).

The check job deletes `work/` (150-250 GB) only when Nextflow exited 0 and every hard check passed. If the published
index differs from the task output in `work/`, the missing files are copied from `work/` and this is recorded.

## What the result can and cannot show

- It **can** give U251 counts free of the `both` bin, measure the clip effect on one library, and say what the control
  libraries' human-classified reads are.
- It **cannot** remove rat reads that xengsort classified as `graft`; check 8 only characterises them in the two
  controls. The graft-only counts lose the human reads that xengsort called `both` (conserved genes lose more), so
  graft-only values are lower than the true human expression for conserved genes, by a gene-specific factor that
  check 7's r_g measures in C2B. Comparisons between arms share that factor; comparisons with external human
  data do not.
- The clip effect is measured on IL69B only. Check 5 bounds it; it does not show that every library moves the same way.
- The rat-leak estimate assumes that r_g is the same in culture and in the brain, i.e. that the isoform mix per gene
  does not change enough to change the fraction of reads with rat-matching k-mers. That is not tested here.
- The ambiguous and neither bins (0.75-2.49 % of pairs) were never saved and are not in any input.
- Nothing here is a primary-versus-recurrent result. With 3 v 3, any group claim needs the 10-split relabelling null
  (S4), and no result can go below p = 0.10.

## Cost

About 4-6 h after prep_reads if the queue is free (the published run took 3 h 27 min with index builds);
150-250 CPU-h. Peak memory: STAR genomeGenerate about 65 GB, STAR_ALIGN about 36 GB, MarkDuplicates up to about 82 GB.
Results about 60 GB plus about 50 GB of indices; work 150-250 GB until the checks pass.

## Files

| file | role |
|---|---|
| `samplesheet.csv`, `params.yaml`, `grid.config` | the run |
| `run.sbatch` | head job with the gate |
| `check_rnaseq_graft.py`, `check_rnaseq_graft.sbatch` | checks; installed to `M/bin/` by `submit.sh` |
| `submit.sh` | run on the login node after copying this folder to `M/rnaseq_graft/` |
| `genesets/neftel_{AC,MES,NPC,OPC}.txt` | gene sets used in checks 5 and 6 |
