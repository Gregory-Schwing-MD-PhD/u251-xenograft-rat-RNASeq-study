# rnavar: nf-core/rnavar 1.3.0 on graft-only U251 reads, plus variant integration

Part of the U251N LITT multi-omics build (SELECTION.md, item S6). This folder holds only configuration, the sample
sheet (sample names and grid paths) and code. No reads, BAMs, VCFs or results are here; they live on the grid under
`/rs/rs_grp_oschome/go2432/u251_multiomics/` (`M` below). `U` is the grid clone
`/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study`; `H` is `/rs/rs_grp_oschome/go2432/u251_host`.

## Questions

1. **Is the model U251?** Authentication against Cellosaurus CVCL_0021 (v52): TP53 p.Arg273His (c.818G>A) homozygous;
   PTEN p.Glu242ValfsTer15 (c.723_724dupTG) present; IDH1 codon 132 reference (none reported, and engineered
   R132C/R132H children exist); HLA class I A\*02:01, B\*18:01, C\*05:01, each homozygous (PubMed 25960936).
2. **Does the clonal make-up move after LITT?** Per-site VAF shifts and arm-level allelic imbalance, with C2B (the
   culture) as the parental baseline, every arm statistic ranked among the 10 balanced relabellings of the six tumours.

## Jobs

| job | script | resources | depends on |
|---|---|---|---|
| test | `run_test.sbatch` (`-profile test,singularity`, separate outdir, work and NXF_HOME) | head 2 CPU / 8 GB / 6 h; tasks capped at 4 CPU / 15 GB | nothing; nothing depends on it |
| head | `run.sbatch` (Nextflow 25.10.7, `nf-core/rnavar -r 1.3.0`) | reqp requeue, 2 CPU / 12 GB / 4 days | afterany collect_prep + stage_refs, then its own gate |
| integration | `integration/integrate.sbatch` + `integration/rnavar_integrate.py` | reqp requeue, 8 CPU / 48 GB / 8 h | afterany head |

`submit.sh` submits all three from the login node (squeue and sbatch only). Logs: `M/rnavar/logs/`. Status files:
`M/rnavar/{HEAD_JOB,INTEG_JOB,TEST_JOB,NEXTFLOW_EXIT,INTEGRATION_EXIT}`, `M/rnavar/test/TEST_EXIT`.

**Gate (head).** Nextflow starts only when `M/reads/READY` and `M/refs/.done/{gencode,known_sites,vep}` exist. If a
producer is still queued, the head resubmits itself afterany on it and exits 0 (at most 72 deferrals). If prep_reads
finished without READY, or stage_refs recorded one of those steps as FAILED, it writes `NEXTFLOW_EXIT=gate_closed`. If
stage_refs hit its time limit, the head re-checks every 30 min, so a resubmitted stage_refs
(`sbatch -D M/refs M/stage_refs/run.sbatch`) is picked up. The integration job follows `HEAD_JOB` the same way.

## Inputs

`samplesheet.csv`: ten rows `<lib>_graft`, `M/reads/graft_clip/<lib>_graft_R{1,2}.fq.gz` (the xengsort graft bin; R2
clipped by 3 nt once, upstream, by prep_reads). All ten are kept: IL64B and N168B (0.31 M and 0.44 M pairs) are cheap and
descriptive; their seq2HLA output is ignored. N269B against IL69B is the only within-animal pair (descriptive, n = 1).
rnavar has no trimming step (FASTQ goes straight to STAR), so the upstream clip is the only preprocessing; the human run
found 121 adapter hits in total, and GATK's RNA practice has no trimming step either.

References (stage_refs, each checked against its publisher): GENCODE v44 GRCh38 primary assembly (FASTA, fai, GTF),
Broad dbSNP 138 and Mills/1000G indels cut to chr1-22,X,Y,M with `##contig` lines rebuilt from the GENCODE fai, the
Ensembl VEP 115 indexed cache, UCSC `hg38ToRn7` chain, REDIportal hg38 v3 (optional).

## Release

nf-core/rnavar **1.3.0** (`-r 1.3.0`; tag commit a0e4641409eb7aed4b325e543d20577f9c08e437). The manifest requires
Nextflow `!>=25.10.4`; the grid env's launcher is 25.10.2, so the job sets `NXF_VER=25.10.7` and puts the jar on CephFS
(`NXF_BIN=M/nxf_dist/25.10.7/...`; the launcher would otherwise write it into the conda env in home).

## Settings that are not pipeline defaults, and why

`params.yaml`

| setting | value | why |
|---|---|---|
| `genome`, `igenomes_ignore` | `null`, `true` | the default `GRCh38` is an iGenomes key with no dbSNP, and main.nf fills dbsnp/known_indels/vep_* from it through `getGenomeAttribute()`; everything is given explicitly instead |
| `fasta`, `fasta_fai`, `gtf` | GENCODE v44 primary (stage_refs) | the release of the published human run, so genes and coordinates join; `fasta_fai` passed so the pipeline does not rebuild it |
| `read_length` | `101` | sjdbOverhang 100 for the index and STAR_ALIGN (R1 35-101 nt, clipped R2 32-98 nt). The default 150 would build an index for reads this data does not have |
| `save_reference` | `true` | keep the STAR index, dict and exon BED |
| `star_twopass` | `true` | default, written out |
| `dbsnp(_tbi)`, `known_indels(_tbi)` | Broad bundle, cut and reheadered | BQSR known sites and HaplotypeCaller rsIDs; header identical to the reference dictionary |
| `tools` | `seq2hla,vep` | ADDS to the defaults (`defineToolsList`), so BQSR, interval scattering, removeunknownregions and VariantFiltration stay on (confirmed in a local stub run: "Tools selected to be run: baserecalibrator,intervallisttools,multiqc,removeunknownregions,seq2hla,variantfiltration,vep") |
| `vep_cache`, `_version`, `_genome`, `_species` | `M/refs/vep_cache`, 115, GRCh38, homo_sapiens | local cache, no download (layout `<cache>/homo_sapiens/115_GRCh38`) |
| `vep_include_fasta` | `true` | offline HGVS needs the reference sequence |
| `vep_custom_args` | `--everything --total_length --offline --format vcf --allele_number` | the default adds `--filter_common` (drops every variant common in gnomAD, i.e. every germline heterozygous SNP the arm-imbalance step needs) and `--per_gene` (one transcript per gene, so the TP53 match on ENSP00000269305 could miss). `--allele_number` puts ALLELE_NUM in CSQ, so each consequence is tied to its ALT without re-deriving VEP's allele trimming |
| `snpeff_db` | `GRCh38.99` | snpEff is **not** run. `VCF_ANNOTATE_ALL` wires `VCF_ANNOTATE_SNPEFF` even when snpEff is off, and with `igenomes_ignore` the value is null, which stops the run at DAG construction ("A process input channel evaluates to null -- Invalid declaration `val db`"). Found by a local `-stub-run` of 1.3.0; GRCh38.99 is the pipeline's own iGenomes value |
| `generate_gvcf`, `seq_platform` | `false`, `illumina` | defaults, written out (per-sample VCFs; RG PL tag GATK needs) |

`grid.config`

- Tasks on `reqp` with `--qos=requeue`; exit codes 130-145, 104, 175-177 (kill, pre-emption, out of memory) retried
  up to 3 times; anything else finishes the run, except (a) any task of IL64B or N168B, and (b) SEQ2HLA for any
  library, which are *ignored* so that two descriptive libraries, or HLA typing, cannot stop variant calling for the
  other eight. The closure reads `meta.id` and falls through to `finish` for processes without a meta; tested on
  Nextflow 25.10.7 (IL64B failure ignored, C2B failure finishes, no-meta failure finishes, SEQ2HLA ignored). The
  checks list every IGNORED task and mark the library partial.
- `resourceLimits` 40 CPU / 240 GB / 96 h. `STAR_GENOMEGENERATE` 80 GB (published human peak 64.8 GB);
  `PICARD_MARKDUPLICATES` 72 GB (the published 81.8 GB peak measured the 200 GB allowance, since Picard's heap follows
  it); `GATK4_SPLITNCIGARREADS` 16 GB / 8 h (default 6 GB); `SEQ2HLA` 12 h. Memory and time scale with the attempt.
- Singularity cache `/rs/rs_grp_oschome/go2432/singularity_cache`: the task template named
  `M/singularity`, but the other items of this build (rnaseq_graft, stage_refs) use this existing CephFS cache, and it
  already holds the STAR image the rat index was built with; one cache is what makes it shared. Every container gets
  `M/tmp` as `/tmp` (node-local /tmp is swept).

`run.sbatch`: NXF_HOME, jar, temp and Singularity scratch on CephFS; one launch directory for this run; `-resume`
always (a requeued head restarts from the top); every FASTQ and reference checked non-empty before launch; Nextflow
must report 25.10.7 or the job stops. It also fetches the two pinned images the integration needs (pysam 0.24.1,
UCSC liftOver 482) into the cache; the integration job fetches them again if this failed.

## Pipeline checks (integration step `checks`; the hard ones decide whether anything else runs)

Each is compared with something the pipeline did not compute itself.

| check | expected | against |
|---|---|---|
| NEXTFLOW_EXIT | 0 | |
| arm labels | primary IL67B/IL68B/IL69B, recurrent IL66B/NL70B/NL71B | `U/ANALYSIS/metadata_therapy.csv` (the published DE design) |
| prep_reads status | OK for all ten | prep_reads checks a-g |
| BAM read group SM | `<lib>_graft` | the samplesheet |
| STAR "Number of input reads" | = N_graft exactly | prep_reads graft record count (xengsort graft bin) |
| flagstat **primary** | = 2 x (unique + multi-mapped pairs) within 1 % | STAR Log.final.out. The task's literal form (flagstat *total* within 1 % of 2 x the STAR input) cannot hold: STAR writes no unmapped reads and adds secondary records. It is reported as INFO |
| VCF records | > 0 for the eight libraries with >= 3.8 M pairs | |
| VEP | command line has `--allele_number`, lacks `--filter_common`/`--per_gene`; rows with gnomAD AF > 0.01 exist; ALLELE_NUM in CSQ; VEP records = VCF records | gnomAD frequencies in the Ensembl cache |
| trace | no `*DOWNLOAD` task; exactly one STAR_GENOMEGENERATE; every HaplotypeCaller interval complete per library | |
| `.command.sh` | STAR_ALIGN `--sjdbOverhang 100 --twopassMode Basic PL:illumina`; STAR_GENOMEGENERATE `--sjdbOverhang 100 --sjdbGTFfile`; VEP flags; HaplotypeCaller `--dont-use-soft-clipped-bases --dbsnp`; BQSR both known-site files | the commands that actually ran |
| saved STAR index | sjdbOverhang 100; contigs = GENCODE fai | stage_refs fai |
| R2 5' artefact | clip + mismatch rate at R2 cycles 1-3 <= 1.5 x cycles 4-10 (+0.005) in every new BAM; also samtools stats `-r` MPC on 300 k R2 records of chr1 | **positive control**: the old graft+both BAMs of the same libraries were never clipped and must show the artefact; if none does, the check is reported as having no power |

## Variant integration (`integration/`)

1. **Union** of PASS sites over the ten `*.haplotypecaller.filtered.vcf.gz` (read from the VEP VCFs, which carry the
   same records plus CSQ), every allele normalised (trim, left-align against GENCODE). Each library's FILTER
   (SnpCluster, FS, QD), GT and HaplotypeCaller AD kept as annotations.
2. **Forced counts** at every union site in every new `md.bam` and every old graft+both BAM:
   `bcftools mpileup -a AD,ADF,ADR,DP -q 60 -Q 20 -d 1000000 -L 1000000 -B` (bcftools 1.21). VAF comes only from these
   counts. Added to the specified command: `-L 1000000`, because the default `--max-idepth 250` silently skips indel
   calling above 250 reads (it would drop the PTEN duplication in a well-expressed gene); `-B`, because htslib skips
   BAQ for any read with an N operation (`realn.c` l.190 in 1.21, "do nothing if there is a reference skip"), so BAQ
   on would recalibrate only unspliced reads and make VAF depend on splicing. Duplicates, secondary, QC-fail and
   unmapped reads are excluded (mpileup default); overlapping mates count once. Indel sites are counted over a window
   (anchor to anchor + allele length + 10) and matched after normalisation, because STAR places gaps anywhere in a
   repeat.
3. **Artefact flags**, each tested against something outside the calls:
   - *rat leak, both bin*: alt-supporting reads in the old BAMs looked up in `M/reads/both_names/<lib>.names.gz`;
     flagged when > 50 % of the pooled alt support (core libraries) is both-bin. The lookup itself is checked: the
     share of both-bin names among a chr1-5 read sample must be the same order as xengsort's both share of the human
     file, or the both-bin filters are not applied;
   - *rat leak, rat base*: lifted to rn7 with the UCSC chain (liftOver 482), alt compared with the mRatBN7.2 base on the
     human strand. Checked by the human-reference-equals-rat-base identity, which must be similar for '+' and '-'
     chain blocks (a strand bug shows as a low '-' identity, and the filter is then switched off);
   - *rat leak, realignment*: alt-supporting pairs (<= 30 per site and BAM) realigned with STAR to the host run's
     mRatBN7.2 index, using the image that built that index; a pair "fits rat better" when its rat AS exceeds its
     human AS; flagged when most pooled alt pairs do;
   - *host fraction*: binomial logistic regression of alt/depth on the xengsort host fraction over C2B and the six
     tumours, flagged when the slope is positive with LRT p < 0.001 and C2B VAF <= 0.05 (the rat-leak shape). This
     replaces a rank correlation of six VAFs, which in the synthetic test missed the designed host-tracking site and
     flagged random ones. **It is confounded with the arm**: the recurrent tumours carry more host reads, so a real
     recurrence-gained allele tracks host fraction too (the synthetic test flagged all three designed
     recurrence-gained sites). It is therefore only in the strict set;
   - *position in read*: alt support within 6 nt of the R1 5' end (new BAMs) or 3 nt of the R2 5' end (old BAMs),
     against the same fractions for reference-supporting reads;
   - *editing*: REDIportal hg38 v3 if stage_refs staged it, otherwise A>G on + strand genes and T>C on - strand genes
     (GENCODE gene strand); the report names the rule used, and reports the strand-matched vs mirrored A>G excess;
   - *allelic sorting*: alt share among graft-bin vs both-bin reads of the old BAMs (Fisher p < 0.01, >= 0.2 apart).
     Whether xengsort routes alt reads differently is not known; this measures it where it can.
   - **Never** "absent in IL64B/N168B": they hold 0.3-0.5 % real U251 reads, so that rule would delete TP53 R273H.
   Filter sets: `raw`; `main` (both-bin, realignment, read position, editing, sorting, and rat base unless the allele is
   common in gnomAD); `strict` (main + every rat-base and host-fraction flag).
4. **Categories** on C2B + the six tumours, >= 20 unique reads in each, Wilson 95 % interval for every VAF;
   present = VAF >= 0.10 with >= 3 alt reads and lower bound > 0.05; absent = VAF <= 0.05: shared clonal (>= 0.3
   everywhere); present in C2B but lost in vivo; gained in one animal; lost in one animal; gained in all three
   recurrences and absent from all primaries and C2B (and the mirror).
5. **Arm-level allelic imbalance** at SNVs with gnomAD AF > 1 % (VEP) heterozygous in C2B (0.1-0.9, >= 20 reads):
   per library and arm, median |VAF - 0.5| (arms from UCSC hg38 cytoBand, whose chromosome lengths must equal the
   GENCODE fai; otherwise whole chromosomes). Also: SNPs homozygous in C2B but heterozygous in a tumour (non-U251 human
   reads), and heterozygous common SNPs on chrX outside the PARs (U251 is male; PARs from the X copies of GENCODE's
   `_PAR_Y` genes).
6. **Y/X**: fragments (R1, MAPQ >= 60, not duplicate) on the exons of RPS4Y1/RPS4X, DDX3Y/DDX3X, KDM5D/KDM5C, UTY/KDM6A,
   EIF1AY/EIF1AX, ZFY/ZFX, USP9Y/USP9X, gene models from the GTF by name (never typed ids). The spec takes Y/X from
   invasion_niche, which is not built; it is counted here from the rnavar BAMs, and an invasion_niche table is
   mentioned in the report if one exists when this runs.
7. **10-split ranks** (split 1 = true labels, IL67B always in A) for: sites gained in all of B, in all of A, in either;
   sites with |mean VAF B - mean VAF A| >= 0.2; summed arm-imbalance difference and arms differing by >= 0.1; |mean
   log2 Y/X difference|. Ties give a rank range; p = worst rank / 10, never below 0.10.
8. **Authentication table**: TP53 key derived twice (GENCODE MANE CDS codon 273, checked to read Arg; VEP match on
   unversioned ENSP00000269305 + p.Arg273His) and VAF >= 0.9 required in C2B and every tumour with >= 20 reads (0.5-0.8
   means non-U251 human or rat reads); PTEN key derived twice (c.723_724 checked to read TG, codon 242 Glu; VEP
   p.Glu242ValfsTer15 or c.723_724dup) and present; IDH1 codon 132 (checked to read Arg) reference fraction >= 0.95 at
   depth and no VEP p.Arg132 call; seq2HLA A/B/C 4-digit calls against A\*02:01, B\*18:01, C\*05:01 homozygous (raw
   file lines kept, because the parse of seq2HLA's table layout has not been checked on real output).
9. Also written: `private_allele_probe.tsv` (private alleles of each donor library in N269B, N168B and IL64B, with the
   shared i7/i5 from `M/reads/index_seq`; a library from another rat that shares an index and shows up there is
   evidence of index hopping), `n269b_vs_il69b.tsv`, `REPORT.md`.

The Nextflow work directory (300-400 GB) is deleted at the end only when the hard checks passed without
`FORCE_INTEGRATION`. `CLEAN_WORK=0` keeps it.

## Tested locally before submission

- `-stub-run` of nf-core/rnavar 1.3.0 on Nextflow 25.10.7 with this params.yaml (paths swapped) and grid.config's
  process block: parameters pass the schema; found the null `snpeff_db` DAG failure (fixed above); confirmed the
  published names the integration reads (`preprocessing/<id>/<id>.md.bam`,
  `variant_calling/<id>/<id>.haplotypecaller.filtered.vcf.gz`,
  `annotation/<id>/<id>.haplotypecaller.filtered_VEP.ann.vcf.gz`, `reports/star/<id>/<id>.Log.final.out`,
  `reports/samtools/<id>/<id>.md.flagstat`, `seq2hla/<id>/<id>-ClassI-class.HLAgenotype4digits`). Three module stub
  blocks needed local-only patches to run (STAR's transcriptome BAM, bgziptabix's same-name output, VEP's `args`);
  none affects the real run, whose script blocks differ.
- errorStrategy closure on Nextflow 25.10.7 (see grid.config above).
- The whole integration end to end on synthetic data (small genome, reads aligned with real STAR 2.7.11b, rat genome
  with one rat-specific base and a minus-strand chain block): the designed recurrence-gained sites rank 1 of 10; each
  designed artefact (both-bin leak, R1 5' artefact, rat-realignment site, editing site, host-tracking site) is flagged by
  its own filter; TP53/PTEN keys agree between the GTF derivation and VEP; IDH1 reads reference; a designed HLA
  mismatch is caught; minus-strand liftover identity 1.0.

## What the result can and cannot show

- It can show whether the U-251MG donor lineage's TP53, PTEN, IDH1 and HLA class I states are present at depth in
  the tumours, the culture and the contralateral hemisphere. That identifies the lineage (which includes SNB-19 and
  U-373MG ATCC), not this lab's stock. A disagreement only with the 1977 serology (A\*01,02; B\*08,14) is not
  misidentification.
- It can show whether site counts, arm imbalance or Y/X differ between the arms more than in the other nine balanced
  relabellings. With 3 v 3 the floor is p = 0.10; nothing can reach p < 0.05.
- It cannot separate LITT selection from founder bottlenecks: whether the arms and C2B share a culture batch or passage
  is not recorded.
- RNA allele fractions carry allele-specific expression and nonsense-mediated decay as well as copy number; only
  expressed genes are seen. Comparing with C2B at the same SNP removes the constant part, not all of it.
- xengsort's ambiguous and neither bins (0.75-2.49 % of pairs) were never saved; whether alleles are routed into them
  differently is not known.
- HaplotypeCaller's detection floor for subclonal alleles at 65-79 % duplication is not known; categories use forced
  counts, but a site no library called PASS is not in the union.

## Cost and timing (estimates)

rnavar about 8 h wall (5-14) once it starts, about 300 CPU-h, largest task about 80 GB (STAR index build); results
about 80 GB plus a 30 GB STAR index; work 300-400 GB, deleted after the checks. Integration about 2-3 h (limit 8 h).
Start: after stage_refs finishes (its own ETA 04:00-07:00 on 27 Sep, limit 12:33) and queue wait on reqp.

## Files

`samplesheet.csv`, `params.yaml`, `grid.config`, `run.sbatch`, `run_test.sbatch`, `test_limits.config`, `submit.sh`,
`integration/integrate.sbatch`, `integration/rnavar_integrate.py` (stdlib + pysam).
