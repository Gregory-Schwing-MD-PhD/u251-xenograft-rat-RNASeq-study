# rnafusion: fusion transcripts as an identity check, and a search for recurrence-only fusions

nf-core/rnafusion **4.1.3** (commit `76ad76e7c3`) under Nextflow **25.10.7**, followed by a custom audit and a second
pass. `M` = `/rs/rs_grp_oschome/go2432/u251_multiomics`, `U` = the grid clone of this repository,
`H` = `/rs/rs_grp_oschome/go2432/u251_host`.

## Questions

1. **Is the line U251, checked against a truth set from outside this study?** CCLE's 2019 STAR-Fusion calls
   (`CCLE_Fusions_20181130.txt`) contain eight U251MG fusions with annotated-splice calls. Each one occurs in U251MG
   and in no other of the file's 1,017 cell lines (checked 2026-09-27): DHX38--AMPD1, CENPU (MLF1IP)--ENPP6,
   PSD3--RBPMS, MICALL2--MAD1L1 (the top four), and CNPY3--SPACA1, ZKSCAN7--KIF15, MRPS18B--XPO5, PCMTD1--TCEANC.
2. **Is any fusion confined to the LITT recurrences?** The test is present in all three recurrences, assessably
   absent in all three primaries, and ranked against the 10 balanced relabellings.
3. **Are N269B's human cells U251N, at clonal level?** If U251N-private fusion reads are present in the contralateral
   library, they came from U251-lineage transcripts.

## How it runs

`bash M/rnafusion/submit.sh` submits four jobs, each `afterany` on the one before. Each job opens only if the marker
files it depends on say the previous step passed:

| job name | script | opens when | does |
|---|---|---|---|
| `u251_rf_p1` | `run.sbatch` (RF_PASS=1) | `M/reads/READY`, `M/refs/.done/rnafusion`, and no prep, stage_refs or refs-recheck job still queued | Nextflow pass 1 |
| `u251_rf_audit` | `audit.sbatch` → `audit_rnafusion.py` | pass-1 exit 0 and `REFCHECK_PASS1` = OK | recall, artefact rules, identity verdict, `whitelist_pass2.tsv`, `AUDIT_OK` |
| `u251_rf_p2` | `run.sbatch` (RF_PASS=2) | pass 1 OK, `audit/AUDIT_OK`, `whitelist_pass2.tsv` present | same launch dir, `-resume`, `--whitelist whitelist_pass2.tsv` |
| `u251_rf_evid` | `evidence.sbatch` → `evidence_rnafusion.py` | pass-2 exit 0 and `REFCHECK_PASS2` = OK | per-library evidence, 10-split group claim |

- **Waiting on upstream jobs.** If a job starts while its upstream job is still queued (for example, a head that
  deferred itself), it resubmits itself `afterany` on that job, by job name, and exits 0. New ids go to `M/rnafusion/JOBS.log`.
- **Pass 1 when a marker is missing.** If nothing upstream is queued and a marker is still missing, pass 1 closes its
  gate at once when stage_refs recorded the rnafusion step as FAILED or the second-node preflight failed. Otherwise
  (for example, stage_refs killed at its 12 h limit writes no STATUS.tsv) it re-checks every 2 h for 24 h. A stage_refs
  that someone resubmits is therefore picked up without resubmitting this chain.
- **Resources.** Head jobs: reqp `--qos=requeue`, 2 CPU, 12 GB, 4 days. Logs go to `M/logs/rnafusion/`.
  Nextflow's own output goes to `M/logs/rnafusion/nextflow_pass<P>_<job>_r<restart>.out`.

**Before `nextflow run`, the head job does the following, on its own node:**
- Checks that the samplesheet is exactly the ten `<lib>_humclip` rows, with reverse strandedness and the expected
  paths, and that every FASTQ exists and is non-empty.
- Checks that every library has a prep_reads check row.
- Reruns `M/stage_refs/refs_util.py preflight`. This emulates `exists_not_empty()`: rnafusion tests one file per
  reference directory, chosen by readdir order on the node where Nextflow runs.
- Checks that every whitelist line is `GENE--GENE`, and that every gene exists both in the GENCODE v46 GTF and in the
  CTAT library's `ref_annot.gtf.gene_spans`.
- Cancels rnafusion tasks orphaned by an earlier attempt of the head. A resumed task reuses the same work directory, and
  an orphan would still be writing into it.

**While Nextflow runs, a watchdog kills it** in three cases:
- Any `BUILD_REFERENCES` task is submitted: a silent rebuild or download. STAR-Fusion's build takes 2 days at 150 GB.
- The launched revision is not `76ad76e7c3 [4.1.3]`.
- In pass 2 only, STAR, a caller, FastQC or MarkDuplicates is submitted, which means `-resume` did not reuse pass 1.

**After the run:**
- Every execution trace is searched for
  `BUILD|GENOMEGENERATE|DOWNLOAD|SALMON_INDEX|GENCODE|HGNC|FAIDX|CREATESEQUENCEDICTIONARY|BIOAWK`. Any hit writes
  `reference_rebuild` to `NEXTFLOW_EXIT_PASS<P>`, and the run is not accepted.
- For pass 2, `PASS2_RERUN_CHECK.tsv` lists every step that did not come from the cache. The expected ones are CAT_CAT,
  FusionInspector, AGAT, VCF_COLLECT, Arriba visualisation and MultiQC.
- Before pass 2 starts, pass 1's FusionInspector-dependent outputs are moved to `results/_pass1/`, so that a pass-2
  table can never be a stale pass-1 file. The evidence script also refuses any table older than `.pass2_start`.

## Inputs

- **Samplesheet** (`samplesheet.csv`, `sample,fastq_1,fastq_2,strandedness`). Ten rows `<lib>_humclip`:
  - R1 = `U/ANALYSIS/sorted_fastqs/<lib>_human_R1.fq.gz`, unchanged;
  - R2 = `M/reads/human_clip/<lib>_human_R2.fq.gz`, the same records with the Takara Pico v2 3-nt 5′ artefact removed
    once, by prep_reads;
  - strandedness `reverse`.
- **Why the full human file (graft + both).** The file holds the xengsort graft reads plus the `both` reads. A junction
  read that spans two conserved exons can land in `both`, so graft-only reads could miss real fusions. Every call is
  then audited read by read against `M/reads/both_names/<lib>.names.gz`.
- **Controls.** IL64B and N168B, whose human files are 93.5 % and 90.4 % `both`, are the artefact filter.
- **The `_humclip` suffix** makes any script that reads the published tables by mistake fail on the join.
- **Why the published rnaseq BAMs are not reused.** They were aligned without chimeric output.
- **References.** The S3 bundle `s3://nf-core-awsmegatests/rnafusion/references/`, mirrored to `M/refs/rnafusion` and
  verified by stage_refs:
  - GENCODE v46 FASTA and GTF;
  - STAR index (sjdbOverhang 99);
  - CTAT library with `cancer_splicing_lib`;
  - FusionCatcher database;
  - Arriba 2.5.0 files;
  - HGNC;
  - fusion-report databases, with `cosmic.db` (0 bytes on S3) removed.

## Settings that are not the defaults, and why

| setting | value | why |
|---|---|---|
| `tools` | arriba, starfusion, fusioncatcher, fusionreport, fusioninspector, ctatsplicing | three independent callers for the 2-of-3 rule; no `salmon` or `stringtie` (they would re-quantify on another GENCODE than the published run); no `fastp` (the reads were adapter-trimmed before delivery and R2 is clipped upstream; a second trim would change which reads the callers see relative to the audit's FASTQ scan) |
| `tools_cutoff` | 2 (default 1) | fusion-report's filtered list, and so FusionInspector's input, needs 2 of the 3 callers |
| `no_cosmic` | true | COSMIC needs credentials, and the S3 `cosmic.db` is empty; FUSIONREPORT_DETECT gets `--no-cosmic` |
| `genomes_base`, `genome`, `genome_gencode_version` | `M/refs/rnafusion`, GRCh38, 46 | the staged bundle; with these every reference path resolves to a checked file and nothing is built |
| `read_length` | 100 | matches the S3 STAR index's sjdbOverhang 99; it would only matter if an index were rebuilt |
| `skip_qc` | false | STAR and duplication metrics are compared with the published human run |
| `whitelist` | `whitelist_pass1.tsv` (the eight CCLE fusions), then `whitelist_pass2.tsv` | FusionInspector evidence for the markers in every library; pass 2 extends it to every candidate. **The schema requires a `.tsv` name** (pattern `\S+\.tsv$`), so the plan's `whitelist_pass1.txt` would have failed parameter validation |
| `seq_platform` | unset | when it is set, `conf/modules.config` l.290 writes `'SM:<id>'' PL:<platform>'`, which fuses SM and PL into one read-group token |
| ARRIBA_ARRIBA `ext.args` | `-X` | read identifiers (and transcript/peptide) also in `fusions.discarded.tsv`, so a missed marker can be traced to its reads and the filter that removed it; the main `fusions.tsv` has read identifiers by default |
| FUSIONCATCHER_FUSIONCATCHER | 16 CPU, 64 GB × attempt, 48 h | C2B has 95.5 M pairs; the module passes 80 % of memory as `--Xmx` |
| STAR_ALIGN | 80 GB × attempt | full GRCh38 index plus on-the-fly sjdb insertion and coordinate sort |
| GATK4_MARKDUPLICATES | 4 CPU, 48 GB × attempt, 12 h × attempt | the pipeline's process_low (12 GB, 4 h) is too small for 95 M pairs with up to 50 multimappers kept; the published Picard run peaked at 81.8 GB |
| executor | slurm, reqp, `--qos=requeue`, retries on 130-145/104/175-177, limits 40 CPU / 240 GB / 96 h | the grid's CPU partition; reqp pre-empts, so kill codes are retried |
| singularity | cache `M/singularity` (shared across the multi-omics pipelines), `-B M/tmp:/tmp`, pull timeout 90 min | node-local /tmp is swept (FastQC died there, job 40468032) |

## The audit (`audit/`, after pass 1)

**Recall of the CCLE markers.** Scored from each caller's own output:
- Arriba `fusions.tsv`;
- STAR-Fusion `abridged.tsv`;
- the FusionCatcher final list.

FusionInspector does not count, because whitelisting puts a fusion into FusionInspector whether or not a caller found
it. Gene names are normalised through HGNC previous symbols, so FusionCatcher's MLF1IP counts as CENPU. The CCLE file's
coordinates are hg19, so matching is by gene name and orientation. Reciprocal calls are listed separately. Arriba
discards of a marker are shown with their filters. Output: `recall_ccle.tsv`.

**Identity rules** (`IDENTITY.md`). A marker is "present" in a library when at least one caller reports it in the CCLE
orientation.
- Two or more of the top four in C2B confirm both the pipeline and the line.
- The markers should be present in C2B and all six tumours. A tumour lacking all four while C2B has them is a
  **mix-up until shown otherwise**.
- One missing marker is not a failure: U251N is not CCLE's U251MG (Torsvik 2014).
- For scale only, if expression matched CCLE (poly(A) there, total RNA here, counted before deduplication):

  | fusion | C2B | tumours |
  |---|---|---|
  | DHX38--AMPD1 | ~316 | ~112-183 |
  | CENPU--ENPP6 | ~48 | ~17-28 |
  | PSD3--RBPMS | ~36 | ~13-21 |
  | MICALL2--MAD1L1 | ~20 | ~7-11 |

  N269B: about 13-27 DHX38--AMPD1 reads. Arriba's count is shown both deduplicated and with its `duplicates(N)` added
  back, for this comparison.

**Artefact rules**, applied to every call in every library (`scored_calls.tsv`):
- **(a) Seen in a control.** Called by any caller in IL64B or N168B. For the CCLE markers this is reported, not used:
  a marker in a control means index hopping or residual tumour.
- **(b) Supported by `both` reads.** More than 50 % of the supporting reads are in `both_names`. Supporting reads are
  Arriba's `read_identifiers` plus STAR-Fusion's JunctionReads and SpanningFrags.
- **(c) Explained by rat.** More than 50 % of the supporting pairs align co-linearly to rat. They are realigned with
  STAR to the host run's mRatBN7.2 index (`H/results_rnaseq/genome/index/star`):
  - no chimeric output;
  - at least 90 % of the pair aligned;
  - at most 4 % mismatches;
  - STAR from the image rnafusion's own STAR_ALIGN used.
- **(d) Rat neighbours.** The rat 1:1 orthologs (`U/ANALYSIS/host/refs_online/rat_human_1to1.tsv`) lie within 200 kb
  on the same strand in `H/results_rnaseq/genome/Rattus_norvegicus.mRatBN7.2.110.filtered.gtf`. A rat read-through
  of neighbouring genes that leaks into the human file would look like a human fusion of the orthologs.

Rules (b) and (c) are marked **not assessable**, never 0 %, in any of these cases:
- the library's supporting read names were not found in its FASTQs (95 % minimum);
- the `both` names file does not have xengsort's N_both lines;
- the library could not be realigned;
- the SAM does not account for every extracted pair.

**A new fusion is reportable only if every one of these holds:**
- 2 of the 3 callers call it;
- it has 3 or more deduplicated fragments at annotated splice sites (Arriba rows with `splice-site` on both sides, or
  STAR-Fusion `ONLY_REF_SPLICE`). FusionCatcher counts are not deduplicated and never count toward this;
- the best breakpoint of every Arriba or STAR-Fusion call of it is on annotated splice sites. This excludes SMART
  template-switch chimeras;
- it is not in a control, and it is not flagged by (b), (c) or (d), with (b) and (c) assessable.

Output: `reportable_new.tsv`.

**`whitelist_pass2.tsv`** contains the eight CCLE markers, the union of every library's fusion-report filtered list and
every Arriba high-confidence call. Entries are dropped when they are intergenic, the same gene on both sides, or have a
gene missing from the CTAT library; the dropped entries and their reasons are in `whitelist_pass2_excluded.tsv`.
`AUDIT_OK` is written only when the hard checks pass:
- every output is present;
- STAR's input pairs equal xengsort's graft + both for all ten libraries;
- the pass-2 whitelist holds all eight markers.

## Pass 2 and the group claim (`evidence/`)

Pass 2 gives FusionInspector evidence for every candidate in every library. Each library then gets one status per
fusion:
- **present:** FusionInspector reports it.
- **absent:** FusionInspector does not report it, **and** the 5′ partner's expression in that library predicts 5 or
  more fragments at the carrier arm's mean FFPM. The prediction is expected = mean(FFPM) / mean(5′-partner CPM) × the
  5′-partner count here, with partner counts from STAR `ReadsPerGene` column 4 (reverse-stranded).
- **not assessable:** anything else.

For each of the 10 balanced splits, the script counts fusions present in all three of one arm and assessably absent in
all three of the other, in both directions. It reports:
- the rank of the true split among the 10 (floor p = 0.10, ties as a range), for all candidates and for the
  artefact-filtered set;
- whether each recurrence-only candidate is already in C2B, which would mean it existed in culture before LITT.

## Checks against something outside this computation

| check | outside reference |
|---|---|
| the eight markers are private to U251MG | CCLE_Fusions_20181130.txt: each of the eight occurs in 1 of 1,017 lines |
| STAR input pairs = graft + both, per library (hard) | xengsort classify logs, asserted by prep_reads |
| STAR unique-mapping % within 8 points; duplication within 0.08 | the published human rnaseq run's MultiQC (`multiqc_star.txt`, `multiqc_picard_dups.txt`; GENCODE v44, trimmed, R2 unclipped, so not identical) |
| `both` names file has N_both lines | xengsort classify logs |
| supporting read names found in the input FASTQs | the FASTQ headers themselves (a name-format mismatch would otherwise read as 0 % `both`) |
| partner counts come from the antisense column | the library was verified reverse-stranded (ISR) by salmon and RSeQC in the published runs |
| launched code = tag 4.1.3 | GitHub `refs/tags/4.1.3` → `76ad76e7c39b2ba9edc35aa3602e3dc454d842ec` |
| no reference rebuilt | the execution traces, and the preflight emulation of the pipeline's own directory test |

## What the result can and cannot show

**It can show:**
- Whether this lab's culture and tumours carry fusions that CCLE found only in U251MG. If they do, the cells are from
  the U-251 donor lineage, and the sample labels agree across C2B and the six tumours.
- Which fusions, if any, are called by 2 of 3 callers and survive every xenograft artefact rule, per library.
- For any fusion confined to one arm: its rank against the 10 relabellings, and whether C2B already carried it.
- Whether N269B holds U251-lineage fusion reads, and at what rate per million graft pairs relative to the tumours and
  to the two controls.

**It cannot show:**
- **That the cells are this lab's U251N rather than another U-251 derivative.** SNB-19 and U-373MG (ATCC) come from the
  same donor, and neither is in the CCLE fusion table, so "private to U251MG" means private among the lines CCLE sequenced.
- **That a missing marker means the wrong line.** U251N differs from U251MG (Torsvik 2014), and total-RNA libraries at
  65-79 % duplication give lower deduplicated counts than CCLE's poly(A) data.
- **A primary-versus-recurrent difference below p = 0.10.** With three against three there are only 10 balanced splits.
  A rank-1 result is the strongest statement this design allows, and it cannot separate LITT selection from a founder
  bottleneck: whether the arms and C2B came from the same culture batch is not recorded.
- **Whether N269B's human reads come from resident cells or from index hopping.**
  - N269B shares its i7 with IL66B and IL69B and its i5 with NL71B and C2B (`M/reads/INDEX_HOPPING.txt`).
  - Comparing N269B's rate with IL64B's and N168B's is the only handle on this. It assumes similar hopping exposure,
    which is not known.
- **FusionCatcher-only fusions cannot be read-audited.** Its final list gives no read names, so rules (b) and (c) are
  not assessable for them. No such fusion can be reportable, because 2 callers are required.
- **A single fragment.** "Not reported by FusionInspector" can still hide one: its default floor is 2 fragments. The
  5-fragment expectation keeps a false "absent" unlikely (P(≤1 | Poisson 5) ≈ 0.04), not impossible.
- **Whether the S3 bundle (CTAT library, FusionCatcher database) is scientifically right.** stage_refs checked that it
  arrived intact and was built on the GENCODE v46 FASTA, and nothing more.

## Cost and timing (estimates; queue waits on reqp not included)

- **Pass 1.** 12-20 h wall and 700-1,200 CPU-h. FusionCatcher takes 2-10 h per library, and container pulls take
  20-40 min on first use.
- **Audit.** 1-2 h, 8 CPU.
- **Pass 2.** 2-4 h and 100-200 CPU-h.
- **Evidence.** Under 15 min.
- **Disk.** Results 150-250 GB, including the STAR, sorted and MarkDuplicates BAMs. Work directory 300-500 GB, to be
  deleted after the evidence step passes.

## Files here

- `samplesheet.csv`, `params.yaml`, `grid.config`, `whitelist_pass1.tsv`: the pipeline's inputs.
- `common.sh`: paths, pins and the upstream-following helper.
- `run.sbatch`: the head job for both passes.
- `audit.sbatch` and `audit_rnafusion.py`: the audit.
- `evidence.sbatch` and `evidence_rnafusion.py`: the pass-2 evidence and the group claim.
- `submit.sh`: the submission.

The two Python scripts use the standard library only. They were run on a synthetic ten-library fixture on
2026-09-27, which exercised a CCLE marker, a `both` plus rat-adjacent artefact, a control artefact and a clean
recurrence-only fusion (rank 1 of 10, p = 0.10). The fixture itself is not kept here.
