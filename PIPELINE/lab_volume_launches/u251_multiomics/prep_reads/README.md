# prep_reads: the one canonical U251 read set

A plain SLURM array (bash, awk, pigz, gzip, md5sum). No Nextflow and no container: the work is cutting and renaming
FASTQ records, and every step is a line of awk that can be read in full. Plan: `SELECTION.md` S1 (scratch
`multiomics/`).

## Question

Every RNA arm of the multi-omics plan (rnaseq on graft-only reads, rnavar, rnafusion, rnasplice, circrna) needs the
same input: the U251 reads with the rat reads out and the Takara Pico v2 artefact off the 5' end of Read 2, removed
**exactly once**. rnavar and rnafusion have no 5' clip of their own. The `clip_r2` parameters of rnasplice and circrna
are not wired to any process. A second clip in rnaseq's TrimGalore would cut three template bases. So the clip is
done here, upstream, and every downstream arm is told not to clip.

Alongside it, the job answers three smaller questions that need a pass over the same files:

1. Does C2B's graft/both boundary sit where the xengsort counts say? The nine in-vivo libraries passed this check in
   `ANALYSIS/host/extract_host.sbatch`; C2B's task was never run.
2. Which reads are in xengsort's `both` bin (they are in both species' files)? Their names let later arms ask whether a
   variant, fusion or junction rests on reads that also match rat.
3. What index pair does each library carry? If N269B shares an i7 or an i5 with a human-rich library, index-specific
   hopping into N269B cannot be excluded as a source of its 4.9 % human reads.

## Inputs (read only, grid)

`U = /wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study`. Listed per library in `samplesheet.csv`.

| file | what it holds |
|---|---|
| `U/ANALYSIS/sorted_fastqs/<lib>_human_R{1,2}.fq.gz` | xengsort `graft` records, then `both` records (`U/main.nf`, SORT_READS) |
| `U/ANALYSIS/sorted_fastqs/<lib>_rat_R1.fq.gz` | `host` records, then `both` records |
| `U/ANALYSIS/xengsort_out/<lib>.txt` | classify counts, parsed exactly as `extract_host.sbatch` does |

The raw sequencer FASTQs are no longer on the grid. The XVAL inventory (RESUME, 2026-09-26 23:45) names SRA runs
SRR39547363-72 whose spot counts equal the xengsort totals; that was not checked here. The xengsort `ambiguous` and
`neither` bins (0.75-2.49 % of each library) were never saved, so nothing here can include them.

## Outputs (`M = /rs/rs_grp_oschome/go2432/u251_multiomics`)

| path | content |
|---|---|
| `M/reads/graft_clip/<lib>_graft_R1.fq.gz` | first `4*N_graft` lines of human R1, byte for byte |
| `M/reads/graft_clip/<lib>_graft_R2.fq.gz` | first `4*N_graft` lines of human R2, 3 nt off the 5' end of every sequence and quality line |
| `M/reads/human_clip/<lib>_human_R2.fq.gz` | every human R2 record (graft + both), clipped the same way; pairs with the unchanged `U/.../<lib>_human_R1.fq.gz` |
| `M/reads/both_names/<lib>.names.gz` | read names (first header token, no `@`) of records `N_graft+1..end` of human R1 |
| `M/reads/index_seq/<lib>.txt` | header comment of the first and the last human R1 record, the modal index field over all records, reads per lane |
| `M/reads/index_seq/<lib>.index_tally.tsv` | every distinct index field in the human R1 headers and its count |
| `M/reads/checks/<lib>.tsv` | one row of every checked quantity; written only when all checks pass |
| `M/reads/checks/<lib>.r2_pos1-10_{before,after}_clip.tsv` | R2 base composition at positions 1-10, first 1,000,000 records |
| `M/reads/CHECKS.tsv`, `M/reads/READY` | written by `collect_prep.sbatch`; READY only when all ten rows say OK and every output exists |
| `M/reads/index_pairs.tsv`, `M/reads/INDEX_HOPPING.txt` | the ten index pairs and the N269B verdict |

**Downstream arms read nothing here unless `M/reads/READY` exists.** A task preempted mid-run can leave some
final-named files from its first steps; READY is the only signal that a library passed.

## Settings that are not tool defaults, and why

| setting | value | why |
|---|---|---|
| R2 5' clip | 3 nt, sequence and quality, R2 only: `awk 'NR%4==2\|\|NR%4==0{print substr($0,4);next}{print}'` | Takara SMARTer Stranded Total RNA-Seq Kit v2 - Pico Input Mammalian (634411) user manual: the first three nucleotides of Read 2 derive from the Pico v2 SMART Adapter and must be trimmed before mapping. Neither published rnaseq run did it (`extra_trimgalore_args: null`; cutadapt ran only `-q 20 -O 1 -a AGATCGGAAGAGC`). |
| gate A | on the first 1,000,000 raw R2 records: mean C at positions 1-2 < 0.08, mean A at positions 1-2 < 0.10, mean C at positions 5-9 > 0.15; otherwise exit 3 and write nothing | the clip is only right if the artefact is there. FastQC of the sorted human files of C2B, IL64B, IL68B and NL70B: positions 1-2 C 2.7-5.8 %, A 4.4-9.0 %; positions 5-9 C 18.9-33.9 %. Reads already clipped, or a library without the artefact, fail the gate. |
| check g | clipped R2 positions 1-2 mean C > 0.10 | after the clip, positions 1-2 are the old positions 4-5, which read 19-28 % C. A clip that did not happen, or happened twice, moves this. |
| max clipped length | 98 | raw reads are 35-101 nt (already adapter-trimmed before delivery), so clipped R2 is 32-98 nt |
| expected xengsort counts | `params.yaml` `expected_<lib>` | the job refuses to cut unless the log parses to the numbers the published human run's FastQC also counted (graft + both = FastQC total for every library) |
| pigz | `/wsu/el7/pre-compiled/python/3.7/bin/pigz`, 6 threads for single streams; 4 + 4 for the two R2 writers that run at once | one decompress-clip-compress pass writes both R2 outputs |
| read-back decompressor | gzip's `zcat`, not pigz | the checks read every output with a second implementation |
| `LC_ALL=C` | | `length()` counts bytes; quality strings are ASCII |
| human-rich | graft / (host + graft + both) >= 0.20 | the lowest tumour (IL66B) is 0.30; N269B is 0.05, C2B 0.86, IL64B and N168B below 0.01 |
| SLURM | reqp, `--qos=requeue`, `--requeue`, `--open-mode=append`, 8 CPU, 8 GB, 6 h per task | CPU-only; memory is a few awk arrays. A requeue restarts the library from scratch (every file is written as `.part` and renamed); a library whose checks row says OK is skipped on a rerun (delete the row to redo) |
| TMPDIR | `M/tmp/prep_reads/<lib>` | node-local `/tmp` is swept |

## Checks (any failure: exit 2, outputs renamed `*.failed`)

- **a.** Human R1 and human R2 each hold `N_graft + N_both` records.
- **b.** Graft R1 and graft R2 hold `N_graft` each; human_clip R2 holds `N_graft + N_both`; the names file has `N_both`
  lines.
- **c.** The name of record `N_graft+1` of human R1 equals the name of record `N_host+1` of rat R1 (the first `both`
  read, seen from both species' files). Also: the first and last names in the names file are that record and the last
  human R1 record.
- **d.** The last names of graft R1 and graft R2 agree; the last names of human R1 and human_clip R2 agree; and every
  human_clip R2 record's name equals its R1 mate's, read in step (stricter than the last-name test).
- **e.** md5 of `zcat human_R1 | head -n 4*N_graft` equals md5 of `zcat graft_R1`; the same for the head of human_clip
  R2 against graft R2.
- **f.** Every clipped R2 sequence is 1-98 nt and its quality line is the same length; no raw R2 was shorter than 4 nt;
  clipped minimum, maximum and total bases equal raw minus 3 (minimum and maximum expected 32 and 98).
- **g.** As above.

## Checks against something outside the computation

- The xengsort counts are asserted against numbers FastQC counted independently in the published human run and the
  host extraction logs, not against themselves.
- The boundary (c) is tested against a different file, the rat R1, which xengsort wrote separately. An off-by-one in
  either direction changes the name.
- Gate A and check g are set from FastQC, a different program reading the same files, and from the kit manual.
- The byte identity (e) reads the outputs back with gzip, not the pigz that wrote them.
- A local dry run on synthetic reads (`scratch multiomics/build/drytest.sh`) gave the intended exit for each failure
  it was built to catch: no artefact (3), rat boundary shifted by one record (2, c1), two R2 mates swapped (2, d3),
  xengsort counts off (2), R2 truncated by four records (2, a4 b3 d2 d4 f6); the clean case passed and a rerun skipped.

## What this can and cannot show

- It shows that the graft-only pairs are exactly the xengsort `graft` records, correctly paired, clipped once.
- It does not show that 3 nt is right for every read. The manual says 3; the gate and check g test the population.
- For inserts shorter than the read, R1's 3' end may carry the complement of the same three adapter bases once the
  Illumina adapter is gone. Nothing here removes them (the kit manual asks only for the R2 clip); STAR soft-clips.
  Not measured.
- It does not measure how much the clip changes gene counts. `rnaseq_graft` check 5 does (IL69B_fullclip against the
  published IL69B).
- The index pair comes from the read headers. If xengsort dropped the header comment, `index_seq` says `NA` and the
  hopping question cannot be answered from these files; the sequencing core's sample sheet would settle it.
- A shared i7 or i5 with a human-rich library means index-specific hopping into N269B cannot be excluded. Unique i7 and
  i5 mean a read reaches N269B only if both indexes switch. Neither outcome rules out carry-over or contamination
  before indexing.

## Run

```bash
# from WSL, local repo -> grid
ssh go2432@grid.wayne.edu mkdir -p /rs/rs_grp_oschome/go2432/u251_multiomics/prep_reads
scp ANALYSIS/multiomics/prep_reads/* go2432@grid.wayne.edu:/rs/rs_grp_oschome/go2432/u251_multiomics/prep_reads/
ssh go2432@grid.wayne.edu bash /rs/rs_grp_oschome/go2432/u251_multiomics/prep_reads/submit.sh
```

`submit.sh` only makes directories, copies the two scripts to `M/bin` and calls `sbatch -D M/reads`. Logs:
`M/reads/logs/u251_prep_<array>_<task>.out`, `M/reads/logs/u251_prep_collect_<job>.out`.

Cost: C2B about 2-3 h, the others 0.5-1.5 h; 20-30 CPU-h; about 55 GB on CephFS.

## Submitted

2026-09-26 ~23:55 EDT: array **40480709** (tasks 1-10 in the order C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B
NL71B), collect **40480710** (afterany). Each task log ends `CHECK OK <lib>`, `CHECK FAILED` (exit 2), `GATE A FAILED`
(exit 3) or `SETUP ERROR` (exit 4). Results: `M/reads/CHECKS.tsv`, `M/reads/READY`, `M/reads/INDEX_HOPPING.txt`.
