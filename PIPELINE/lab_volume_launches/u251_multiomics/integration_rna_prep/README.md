# integration_rna_prep: the RNA half of the methylation integration

SELECTION S8 (EVAL_integration P0). A custom step rather than an nf-core pipeline, because no nf-core pipeline takes
bulk RNA-seq together with an Illumina methylation array. It needs nothing from the methylation run and reads nothing
from `/rs/rs_grp_oschome/go2432/u251_meth`, `ANALYSIS/methylation/` or any methylation job. It prepares every RNA-side
table the integration (I1-I5) needs, so that step can run as soon as the methylation run finishes.

## Questions

1. What share of each library's RNA is U251? This is the RNA-side partner of the DNA tumour fraction (I1).
2. How many of each gene's "human" reads are rat reads leaking through xengsort's `both` bin? The answer is per gene
   and per library, and it is used to drop genes from 5 and from I3/I4.
3. Is MGMT silent at the locus level, in culture, in the primaries and after LITT? What is the detection limit? And is
   CDKN2A/CDKN2B absent? This is the therapy question, answered later against STP27 in I2.
4. Stage the external anchors: DepMap 24Q4, CCLE RRBS, the EPIC hg38 manifest and the UCSC hg38 cytoBand.
5. What is the expression dosage per chromosome arm? This is the RNA partner of the array copy number (I3).
6. Is C2B the U-251 lineage at all? If it is not, the job stops.

## Run

```
wsl -e bash -lc 'scp -r ANALYSIS/multiomics/integration_rna_prep go2432@grid.wayne.edu:/rs/rs_grp_oschome/go2432/u251_multiomics/'
ssh go2432@grid.wayne.edu bash /rs/rs_grp_oschome/go2432/u251_multiomics/integration_rna_prep/submit.sh
```

`submit.sh` copies the scripts to `I/bin` and submits two jobs, each with `sbatch -D I`, where
`I = /rs/rs_grp_oschome/go2432/u251_multiomics/integration`:

- **`u251_irp`**, with `IRP_MODE=published` and no dependency. It runs items 1 to 6 on the published human run.
- **`u251_irp_graft`**, with `IRP_MODE=graft`. It depends `afterany` on the rnaseq_graft check job and on `u251_irp`,
  and reruns items 1, 2, 5 and 6 on the graft-only counts as a sensitivity analysis.
  - Its first step refuses to run unless `rnaseq_graft/checks.json` says `all_hard_pass: true`.
  - If an rnaseq_graft job is still queued, it re-queues itself behind that job, up to 8 times, instead of waiting.

Both jobs use `-p reqp --qos=requeue --requeue`, 4 CPU, 32 GB and 1 h. `TMPDIR=I/tmp`, and the containers' `/tmp` is
bound to it. The job is idempotent: a requeue redoes everything except the reference downloads whose md5 already
checks out.

| path | what |
|---|---|
| `I/results/p0/` | published run: every `*.tsv`, `summary.json` (status, every check, key numbers, md5 of every input including the ten BAMs and both images) and `REPORT.md` |
| `I/results/p0_graft/` | the graft-only rerun, with the same layout |
| `I/refs/` | references, each with a `.md5` file, plus `refs_status.tsv` |
| `I/tmp/p0_<mode>/` | intermediates: GTF extract, loci, raw locus counts, `env.tsv` |
| `I/logs/u251_irp_<job>.out` | log |

`summary.json` `status` is one of:
- `complete`;
- `complete_with_failed_checks` (the failed non-fatal checks are listed);
- `failed:<check>` (a fatal check failed);
- `STOP_identity`;
- `error: ...`.

The downstream I-job must gate on `status` and must read `failed_checks`.

## Inputs (read only)

`U = /wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study`

- xengsort classify logs `U/ANALYSIS/xengsort_out/<lib>.txt`: the five bins, parsed by their header row.
- The published human nf-core/rnaseq 3.22.2 run, `U/ANALYSIS/results_human_final/`:
  - `star_salmon/salmon.merged.gene_{counts,tpm}.tsv`;
  - `star_salmon/<lib>.markdup.sorted.bam`;
  - `star_salmon/log/<lib>.Log.final.out`;
  - `star_salmon/picard_metrics/*.MarkDuplicates.metrics.txt`;
  - `multiqc/.../multiqc_fastqc_fastqc_raw.txt`.
- The GENCODE v44 primary-assembly GTF `U/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz`, and
  GENCODE's own `MD5SUMS` beside it.
- The graft-only rerun reads `M/rnaseq_graft/results/star_salmon/salmon.merged.gene_{counts,tpm}.tsv` (columns
  `<lib>_graft`) and `M/rnaseq_graft/checks.json`.

The samplesheet holds only sample names and grid paths. `params.yaml` holds every setting, and the job reads it at run
time.

## What each item does, and every non-default setting

### 1. Calibrated human share (`h_share.tsv`, `calibrator_rates.tsv`)

- **Rates.** Fractions are taken of all five bins (host, graft, ambiguous, both, neither).
  - Pure human = C2B.
  - Rat = the mean of IL64B and N168B.
  - The job computes the rates from the logs, and refuses to run unless each rounds to the value stated in the spec:
    C2B graft 0.8619, both 0.1124, host 0.0046; IL64B/N168B host 0.9359/0.9357, graft 0.0033/0.0055, both 0.047/0.051.
- **Estimates.**
  - `h_graft = (G - g_R)/(g_H - g_R)`, from the graft bin.
  - `h_host = (r_R - H)/(r_R - r_H)`, from the host bin.
  - Both are carried forward as an interval `[h_lo, h_hi]`. Nothing downstream picks whichever suits it.
- **Fatal asserts.**
  - Ten libraries.
  - The five bins summed over the ten libraries give 822,628,276 pairs.
  - graft + both equals the FastQC Total Sequences of the human file, for R1 and R2 of every library.
  - The GTF's md5 equals GENCODE's `MD5SUMS` entry, and its header says `version 44 (Ensembl 110)`.
- **`h_graft` and `h_host` must agree within 0.03.** This check is recorded, not fatal (deviation 1 below).
- **A third-bin closure check.** Neither estimate uses the `both` bin, so `both_obs_minus_pred` (observed `both`
  fraction minus `h_host·b_H + (1-h_host)·b_R`) tests whether the two-rate model holds.
- **A chrY estimate, `h_chrY`.** The rats are female and U251 is male. `h_chrY` is chrY reads per total pair, relative
  to C2B, over the non-PAR protein-coding chrY genes at C2B CPM > 10. The genes are taken from the GTF, never typed
  in. This estimate uses no xengsort rate at all.

### 2. Rat reads per gene (`rat_leak_genes.tsv`, `rat_leak_by_library.tsv`)

- **Formula:** `expected_rat(g,t) = (1 - h_host,t) × total_t × (g_R + b_R) × p_ctrl(g)`, where
  `p_ctrl(g)` = the mean over IL64B and N168B of (human-file count of g) / (human-file read pairs).
- **Output:** `rat_leak_frac = expected / observed`. A gene with frac > 0.1 is marked in that library.
- In graft mode the multiplier is `g_R`, and `p_ctrl` is taken per graft pair from `IL64B_graft` and `N168B_graft`.
- `h_host` is used because the host bin is the direct measure of rat reads.
- **Sensitivity columns:**
  - `h_graft` in place of `h_host`;
  - a both-bin residual variant: rat both-reads = observed both − `h·N·b_H`;
  - each control's profile alone, which gives the spread;
  - how many marks each variant changes.
- **Self-check:** the calibrators run through their own model should give about 1.

### 3. MGMT, CDKN2A and CDKN2B loci (`mgmt_cdkn2a_loci.tsv`, `locus_nm_hist.tsv`, `locus_definitions.tsv`, `library_denominators.tsv`)

- **Coordinates and exons.** Gene spans and the union of exons of all transcripts come from the v44 GTF, by
  `gene_name`, and must be a unique protein-coding record.
- **Positive controls.** Two genes on the same chromosome with C2B TPM 1-5, not overlapping the target, not histone,
  nearest the target. The dry run picked C10orf143 and JAKMIP3 for MGMT (chr10), and DMRTA1 and ADAMTSL1 for
  CDKN2A/CDKN2B (chr9).
- **Reads.** `samtools 1.21 view -q 255 -F 0x904`: STAR-unique primary alignments only.
  - Counted per fragment (read name). A fragment is exonic if any aligned CIGAR block of either mate overlaps an exon,
    and intronic otherwise.
  - A read whose `N` gap only jumps over the gene is not counted.
  - Counts are given with and without Picard duplicates, and split sense/antisense. The library is reverse-stranded
    (ISR), so Read 1 is antisense.
  - The NM distribution is given per locus and library. Rat reads forced onto a human gene carry mismatches.
- **Detection limit.** A one-sided Poisson 95 % upper bound on MGMT CPM, which is 3.0 reads at k = 0 (the rule of
  three). It is computed over the salmon assigned counts, and for the BAM over deduplicated pairs.
- **DepMap anchors** (`depmap_anchor_genes.tsv`): MGMT, CDKN2A, CDKN2B and GAPDH in U251MG, U251MGDM and T98G, the
  MGMT-expressing positive control looked up in Model.csv, plus the pan-cancer median.

### 4. References (`I/refs/`, `refs_status.tsv`)

- Downloaded with GET only (`curl -L`), because figshare refuses HEAD.
- Seven DepMap 24Q4 files from figshare, and CCLE_RRBS_TSS1kb_20181022, each asserted against the md5 in the spec.
- The EPIC hg38 manifest (zhou-lab InfiniumAnnotationV1) and UCSC hg38 cytoBand have no published md5, so their md5
  is recorded on first download and checked on every rerun.
- A mismatch is kept as `<name>.md5_mismatch` and never used. The five files items 5 and 6 need are fatal if missing.
- `ccle_rrbs_loci_cns.tsv` extracts the MGMT/CDKN2A/CDKN2B/GAPDH TSS windows for the 44 CNS lines, for I2 and I4.

### 5. Arm dosage (`arm_dosage.tsv`, `arm_dosage_ngenes.tsv`, `arm_dosage_mark_sensitivity.tsv`, `depmap_arm_cn_family.tsv`, `arm_dosage_vs_depmap_cn.tsv`)

- **Score.** log2(TPM+1) minus the DepMap pan-cancer median of the same gene, over all 1,673 models.
- **Genes.**
  - GENCODE v44 protein-coding genes matched to DepMap's symbol by name. Names duplicated within either source are
    dropped, and every count is reported (`summary.json` `gene_matching`).
  - Histones are removed with the same regex as rnaseq_graft check 8, because DepMap is poly(A) and Pico v2 is total
    RNA.
  - Per sample, only genes at TPM ≥ 1 and not marked by item 2 are kept.
- **Arms.**
  - An arm is p or q by the gene midpoint against the hg38 centromere, which is the end of the p-arm `acen` band.
  - An arm is scored only with ≥ 20 genes.
  - The score is the median over the arm, centred on the sample's median arm.
  - Arm names (`1p` ... `Yq`) match the methylation run's `cnv_arm_log2.tsv`.
- **The same code runs on DepMap U251MG (ACH-000232) and U251MGDM (ACH-001172).** SNB19 (ACH-001198) is in Model.csv but
  has no omics profile at all in 24Q4 (no row in OmicsProfiles.csv), so it cannot be run. The job says so rather than
  substituting another line.
- **Positive control for the method.**
  - DepMap U251MG's own arm dosage against its own arm copy number, from `OmicsCNSegmentsProfile.csv` through
    `OmicsDefaultModelProfiles.csv`: log2 of the length-weighted mean relative copy ratio, median-centred.
  - C2B against DepMap CN of every CNS line (`identity_c2b_dosage_vs_cns_cn.tsv`) is a karyotype-level identity check,
    reported but not gating.

### 6. Identity (`depmap_identity.tsv`, `depmap_identity_by_library.tsv`)

- **Metric 1.** Spearman correlation between C2B and each DepMap CNS/Brain line with expression, over the 5,000 genes
  most variable across those lines.
- **Metric 2.** Pearson correlation of profiles centred on the pan-cancer median, over all 18,858 matched genes.
- **STOP** (exit 5, no arm dosage written) unless the nearest line is U-251 MG or U-251 MG DM on **both** metrics
  (deviation 3).
- The other nine libraries are reported, not gated.

## Deviations from the spec, and why

1. **`|h_graft - h_host| ≤ 0.03` is recorded, not asserted.**
   - With the spec's own rates, IL66B gives 0.3383 against 0.3696, a difference of 0.0313. The preview's
     "within 0.03 in every library" came from rounding (0.339 against 0.370).
   - A fatal assert would stop the whole job on a threshold that the preview already crossed.
   - The disagreement is itself a finding. IL66B also carries 2.3 points more `both` reads than the two-rate model
     predicts, while the other five tumours carry 0.5-1.7 points fewer.
   - So the check fails loudly in `summary.json`, and both estimates travel as an interval.
2. **The leak multiplier is `g_R + b_R`, not `b_R`, and `p_ctrl` is per read pair, not a profile that sums to 1.**
   - The human file holds rat reads from the graft bin (rate g_R) as well as from the both bin.
   - `p_ctrl` is measured per read of the controls' human file.
   - With `b_R` alone, the expected rat counts are about 8 % low.
   - With `p_ctrl` normalised to sum 1 over counted reads, they are 2.6-3.6 times too high: salmon assigns only 28 %
     (IL64B, 1.36 M of 4.86 M pairs) and 38 % (N168B, 1.70 M of 4.52 M) of the controls' human-file pairs. That would
     be a wrong but plausible number.
   - The self-check (the calibrators run through their own model should give about 1) is what shows which
     normalisation is right.
3. **The identity STOP requires both metrics.** In the dry run, metric 1 also put the two mostly-rat controls, IL64B and
   N168B, nearest U-251 MG, by margins of 0.01 and 0.04. Alone it is not a strict enough gate.
4. **References go to `I/refs`,** as the spec's WHERE line says, not to `M/refs`.
5. **Added outputs beyond the spec.** Each is an outside or closure check on a number the spec asks for:
   - `h_chrY` and the both-bin closure;
   - the per-control leak bounds;
   - the DepMap CN positive control;
   - the RRBS extract.

## Checks against something outside the computation

| number | checked against |
|---|---|
| xengsort parse | FastQC Total Sequences of the published run (a different tool on the same files), and the grand total 822,628,276 from SCOUT |
| calibrator rates | the values stated in the spec, to their printed precision |
| GTF | GENCODE's own MD5SUMS list and the header version line |
| human share h | a second bin (h_graft against h_host), a third bin (the both closure) and a species marker (chrY, female host) |
| rat leak | the calibrators through their own model; later the count-based estimate from rnaseq_graft check 7, which the graft rerun recomputes and compares (`rat_leak_two_estimates.tsv`) |
| MGMT zeros | same-chromosome positive controls counted by the same code in the same BAMs; DepMap U251MG MGMT TPM; CCLE RRBS U251MG MGMT promoter methylation |
| arm dosage | DepMap U251MG and U251MGDM through identical code, and each against its own DepMap copy number |
| C2B identity | the DepMap CNS panel (Model.csv lineage), with two metrics |
| family ids | Model.csv names (ACH-000232 = U251MG, ACH-001172 = U251MGDM, ACH-001198 = SNB19) |
| references | spec md5s (DepMap, CCLE) |

## Dry run (local, 2026-09-27; the grid job's numbers are the record)

The R part was run locally in R 4.3.2 on the same published tables, xengsort logs, GTF and references. The BAM counts
were synthetic, so none of the locus numbers below come from it. Everything else is real:

- **Human share.**
  - h_graft/h_host: IL66B 0.338/0.370 (the 0.031 flag); IL67B 0.491/0.485; IL68B 0.524/0.518; IL69B 0.746/0.730;
    NL70B 0.493/0.487; NL71B 0.380/0.369; N269B 0.052/0.071.
  - h_chrY runs below both in every in-vivo library: IL66B 0.23, IL67B 0.36, IL68B 0.41, IL69B 0.61, NL70B 0.46,
    NL71B 0.32, N269B 0.044. Either the xengsort rates overstate h, or U251 in vivo expresses less chrY per unit RNA
    than in culture (Cellosaurus records Y-negative subclones). This step cannot tell which.
- **Rat leak.**
  - Rat share of counts, by the model: IL66B 4.7 %, IL67B 2.7 %, IL68B 2.4 %, IL69B 0.9 %, NL70B 2.6 %, NL71B 3.9 %.
  - The bounds from each control's profile alone are about ±15-20 % of these values.
  - The both-residual variant gives IL66B 8.0 % and IL69B 0.1 %.
  - Genes marked at frac > 0.1: 160-455 of 11.8k-13.2k expressed genes per tumour.
  - Self-check: the calibrators' totals come out at 1.25 (IL64B) and 0.82 (N168B), not 1. The two control profiles
    differ per gene (IQR of log2 ratio −0.86 to 0.34), so the per-gene leak is uncertain by about a factor of 2 and
    the flag fails.
- **Identity.**
  - C2B's nearest CNS line is U-251 MG on both metrics: ρ 0.867 against 0.765 for the next line, DBTRG-05MG, a margin
    of 0.10; centred r margin 0.13.
- **Arm dosage.**
  - The method's positive control, DepMap U251MG against its own CN, gives Spearman 0.63 over 41 arms, and U251MGDM
    0.72.
  - C2B against DepMap U251MG CN gives only 0.22, and U-251 MG ranks 3rd of 93 CNS lines. Either the lab's U251N
    karyotype has drifted from the DepMap stock (Torsvik et al. 2014 documents drift among U-251 stocks), or the
    total-RNA library adds arm-level offsets. I3 separates the two with C2B's own array copy number.
  - Swapping the leak marks moves any tumour arm by at most 0.035 (N269B 0.08).
- **DepMap anchors.**
  - MGMT TPM: U251MG 0.16, U251MGDM 3.92, T98G 13.4.
  - CDKN2A TPM: U251MG 0.05, U251MGDM 0.
- **CCLE RRBS.**
  - U251MG MGMT TSS window 0.300.
  - No U251MG value at any CDKN2A/CDKN2B window.

## What the result can and cannot show

**Can:**
- A tumour-content denominator for every RNA library, from two xengsort bins and from chrY. Where these disagree, the
  size of the disagreement.
- Which genes' human counts are more than 10 % rat in which library, so that I3 and I4 can drop them explicitly.
- Whether MGMT is silent across the whole locus (exonic and intronic), in C2B and in every tumour, with a stated
  detection limit and positive controls that show the zero is not a counting artefact.
- Whether C2B is the U-251 lineage by expression.
- Arm-level expression dosage on the same arm names as the array copy number.

**Cannot:**
- Which absolute tumour share is right. The ratio of RNA per cell between U251 and rat brain is unknown.
- The rat share of any single gene to better than about a factor of 2. The two control profiles disagree.
- Anything about temozolomide response, or any primary-against-recurrent difference: every group comparison here is
  descriptive, and n = 3 against 3 has a floor of p = 0.10.
- Absolute copy number, or events below arm scale.
- Whether U251N's karyotype matches DepMap's U251MG stock. The expression dosage is only weakly correlated with
  DepMap's CN (0.22), and only the array (I3) can settle this.
- That C2B is this lab's stock specifically. DepMap identity places it in the U-251 lineage, which includes SNB-19
  and U-373MG (ATCC).

## Not known

- Why h_chrY runs 5-37 % below the xengsort estimates in vivo: lower chrY expression, Y-negative subclones, or rates
  that do not transfer from culture. The Y/X gametolog ratios in invasion_niche part A bear on this.
- Why IL66B carries excess `both` reads (+2.3 points), and whether they are rat or human.
- Whether the rat reads in a tumour-bearing hemisphere leak like those of the two tumour-free controls. The leak
  profile comes from the controls. Tumour-associated rat cells (myeloid, stromal) may express different conserved
  genes.
- Whether IL64B's and N168B's graft reads are rat mis-sorted as graft or real human reads. rnaseq_graft check 8
  decides this. It matters for the graft-mode leak model, which assumes they are rat.
