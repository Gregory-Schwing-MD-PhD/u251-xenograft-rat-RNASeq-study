# Working checkpoints: U251 LITT recurrence study

Newest first. Everything for this study lives in this directory and this repository (Greg, 2026-09-26: "this project
should be contained in the u251 dir and the respective git repo"). The checkpoints below were first written into
spinesurg-ct-nnunet/docs/RESUME_2026-09-22_PM.md by mistake and were moved here the same day.

## 2026-09-27 00:36 EDT (grid clock): pipeline_v2 (rat host-DNA benchmark + downstream) SUBMITTED: 40481316 -> 40481317 -> 40481327 (held)

pipeline_v2 = every published host-DNA method with rat as the host, judged by PREREG.md (sha256 `67e9e250...5832`,
checked at the start of every job). Code: local `ANALYSIS/methylation_v2/` (README.md explains it; nothing committed),
grid `/rs/rs_grp_oschome/go2432/u251_meth/pipeline_v2/` (results in `pipeline_v2/results/`, work `u251_meth/work_v2`,
launch dir `u251_meth/launch_v2`, logs `u251_meth/logs/pipeline_v2_<id>.out`). Bespoke pipeline dir untouched.

| job | what | depends on | ETA (EDT, 27 Sep) |
|---|---|---|---|
| **40481316** | bench, `WITH_E5=false`: SETUP_V2, CLASSES, STAGE_C2B (C2B only), PREP (GSE174568 13 + C2B + GSE310817 37), PREP_V2 (GSE299969), E1-E4, F1/F2/F4/F5, C0, C1 (a/b/c), C2 fraction (400 mixtures) | afterany:40479710 (satisfied) | running since 00:32; done ~01:15-01:45 |
| **40481317** | bench, `WITH_E5=true`, `-resume`: E5, C2 exclusion, SELECT -> `results/selection.json` | afterany:40481316 and afterany:40481260 (ref collect, writes `ref/e5/DONE`) | ~15-30 min after 40481260; the bowtie-1 index builds (40481256/7) started 00:21, so roughly 03:00-05:30 |
| **40481327** | apply: STAGE_ALL, C3, QC_NORMALISE, ADAPT, CNV, DMP_DMR, REPORT | afterany:40481317, **HELD** | ~1 h after release |

**Manual gate before releasing 40481327 (PREREG section 2, order step 4):** when 40481317 has written
`results/selection.json`, put its sha256 (printed at the end of `logs/pipeline_v2_40481317.out`) as a new block at the
top of this file, then `echo <sha256> > .../pipeline_v2/results/selection.recorded` and `scontrol release 40481327`.
The loader refuses every tumour array until that file matches. If 40481316 or 40481317 fails, fix the code in
`ANALYSIS/methylation_v2/`, rsync it to the grid (excluding `ref_jobs`, `PREREG.*`, `results`), and resubmit the same
phase. `-resume` reuses every finished task.

Already seen at 00:36: CLASSES reproduced every manifest count that was fixed before the freeze. That covers:
- switch probes: 24 rat at AS 50 (16 G->R, 8 R->G), 74 at AS >= 45; mouse 23 / 71;
- E4a 25,672; E4b 58,599;
- R1: 24,936 unmasked cg, 0 disagreements;
- 46 chrY probes at rat AS >= 30.

These are manifest-only facts, not results.

## 2026-09-27 00:21 EDT (grid clock): BENCH reference data + E5 rat alignment chain SUBMITTED (40481253-40481260)

Rat host-DNA benchmark (BENCH workflow, PREREG `meth_bench/PREREG.md` sha256 67e9e250...). Public data only; nothing
from our chip. Scripts: local `ANALYSIS/methylation_v2/ref_jobs/` (not committed), grid copy
`u251_meth/pipeline_v2/ref_jobs/`; job ids also in `pipeline_v2/ref_jobs/JOBS.tsv`. Output names follow
`meth_bench/INTERFACE.md` (the pipeline_v2 reading contract). All reqp/requeue, logs `u251_meth/logs/ref_*_<id>.out`.

| job | what | ETA (EDT, 27 Sep) | output |
|---|---|---|---|
| 40481253 u251_ref_geo | GSE174568 (28), GSE310817 (37), GSE299969 (12, EPIC v2) RAW tars unpacked; Needhamsen files 1 and 4; Illumina B5 + Zhou manifests; sha256 + gzip + IDAT address-count table | ~01:30 | `ref/geo/<GSE>/`, `ref/geo/<GSE>/CHECKSUMS.tsv`, `ref/geo/GEO_FETCH_STATUS.txt` |
| 40481254 u251_ref_sesame | names(EPIC.addressSpecies$species) + rat flag (expect TRUE, Rnor_6.0) | ~00:45 | `ref/sesame_species.txt`, `ref/sesame_species_rat_flag.txt`, `ref/sesame_species_detail.tsv` |
| 40481255 u251_ref_tools | Bismark 0.14.5 + Bowtie 1.1.2 (published pair), fallback Bismark 0.20.0; smoke test picks the pair | ~00:45 | `u251_meth/envs/ref_tools.env` |
| 40481256 / 40481257 u251_ref_prep_{Rnor6,mR7} | Rnor_6.0 (Ensembl 101) / mRatBN7.2 (Ensembl 110) + `bismark_genome_preparation --bowtie1` (single-threaded bowtie-build: the long step) | ~04:00-07:00 | `ref/rat_genome/<asm>/genome/Bisulfite_Genome/PREP_DONE` |
| 40481258 / 40481259 u251_ref_align_{Rnor6,mR7} | Needhamsen's exact `bismark --bowtie1 -n 1 -l 28 -f` on their FASTA (+ `--ambiguous --un`); mR7 also the R-expanded deviation; spiked-mismatch test; Rnor6 vs Zhou's rat mapping | ~06:00-10:00 | `ref/e5/E5_Rnor6.tsv.gz`, `E5_mR7.tsv.gz`, `E5_mR7_Rexp.tsv.gz`; `ref/rat_alignment/checks/` |
| 40481260 u251_ref_collect | writes `ref/e5/DONE` (built / not built + reason), merged hits, PREREG lists | ~10:00 | `ref/e5/DONE`, `ref/rat_alignment/epic_v1_rat_hits.tsv`, `ref/rat_alignment/lists/E5-*.txt` |

Found while writing the jobs: Needhamsen's Additional file 1 is the reverse complement of AlleleA with every probe `R`
already set to `A` (checked on the first 46,868 records of the file; the full-file count is in
`ref/rat_alignment/checks/fasta_check_*.tsv` once the align jobs run), so Bowtie never saw an `R` in the published run, and after
Bismark's C->T read conversion every R expansion collapses to the same read. The METHODS.md E5 step 3 premise ("Bowtie
treats R as a mismatch") does not hold for their input; E5_mR7_Rexp is still built as INTERFACE.md specifies.

## 2026-09-26 23:58 EDT (grid clock): nf-core/methylarray arm FIXED and RESUBMITTED as 40480724 (40480122 failed)

- **Why 40480122 failed (21:51-21:56)**: SeSAMe sample QC dropped all six tumours ("Samples: initial = 6 | kept = 0"),
  so DENSITY_PLOTS, SEX_QC and SNP_HEATMAP failed on an empty matrix. PR #41 counts a probe as failed if it is NA in the
  SeSAMe betas, and the Q step (quality mask) alone makes ~11-12 % of EPIC v1 probes NA in every array, so the default
  `fail_sample_fr` 0.10 cannot be met on EPIC v1. Rates: IL66B 0.149, IL67B 0.142, IL68B 0.139, IL69B 0.136, IL70B
  0.142, IL71B 0.161, each 0.109-0.118 above the bespoke pOOBAH fraction (0.037/0.026/0.021/0.018/0.028/0.052). The mask
  share is inferred from that constant gap, not counted; the new run counts it (`N_masked_without_P`). The earlier
  README sentence "any tumour dropped would be a recurrent one" was a prediction that the run contradicted; rewritten.
- **Fix (patch P6, `bin/sesame_poobah_qc.R`)**: new `fail_sample_scope` = `unmasked`: SeSAMe also runs with P removed
  (`QCDB`) and the per-sample rate is the share of probes that pass every other step but are NA after `QCDPB` (the
  pOOBAH failures the parameter's help text describes). `fail_sample_fr` stays 0.10. Mock-tested locally (sesame is not
  installed locally); first real execution is 40480724. Chosen over raising `fail_sample_fr` to ~0.22, which would
  have given the threshold a different meaning per array.
- **Also fixed**: P3's DMRcate "no regions" pattern now matches DMRcate's misspelled message ("No signficant regions
  found"); `bin/relabel_null.py` ranks ties AGAINST the true split and prints permutation p = #(count >= true)/10 and
  the tie count (all-zero counts read "10 of 10, p = 1.00", not "1 of 10"); the rat list is content-addressed
  (`inputs/rat_exclusion_list.6f4a850ed0c60e0b.txt`, named in params.yaml, name checked against sha256 every run) so a
  changed list is a changed path and `-resume` re-runs SPLIT_COLLAPSE. Bespoke run 40480154 (done 22:10) did not change
  the list: sha256 `6f4a850e...` (85,764 IDs), same as 40480027's.
- **New patch**: `patches/pr41_3a7ff00_epicv1.patch` sha256 `0ef3365e6ef4d5a3409c8006f1c9f5d343d284175d82a4587aabe7cdbc35317e`
  (applies cleanly to pristine 3a7ff00 with git apply and patch -p1; results byte-identical to the patched tree).
  Submitted with `U251_MOVE_ASIDE=1`: the job moves 40480122's work/, results/, relabel/, launch/ to
  `nfcore_methylarray/superseded/patch_52e2e99892398c5b_moved_by_40480724/` and starts clean.
- **Job 40480724** (reqp, requeue, 8 CPU, 64 GB, 8 h): RUNNING on amx2 since 23:59:12. Checked at 23:59:34: steps
  1-3 passed (zip sha256, IDAT md5s, rat list matches its name, 40480122's dirs moved aside, patch 0ef3365e applied).
  Estimate, since P6 and the DMP/DMR steps have never run: main-run results ~00:20-00:45 EDT 27 Sep; `results/relabel_null.tsv`
  ~01:15-02:15 EDT. Read `run_status_40480724.txt`, `results/u251_checks.txt`,
  `results/sesame_poobah_qc/masking_failrate_per_sample.csv` (all six tumours kept? mask count per array).
- **Bespoke pipeline** is now 40480154 (COMPLETED 22:10); 40477975 and 40480027 are history. Not touched here.
- **`-profile test`** at the pinned commit failed at the container pull (quay.io 401) before reaching its test data;
  the test-data samplesheet URL also returns 404. Recorded only (`logs/test_profile_result.txt`), not re-run.
- Housekeeping: RESUME.md back to LF (the working tree had become CRLF; HEAD is LF); `MBR/` added to
  `.git/info/exclude` (untracked, no commit); the lab files copied into the session scratch
  (`scratchpad/methylarray/labfiles/`) deleted. Nothing committed or pushed.

## STATE 2026-09-26 23:45 EDT — SAVED ON GREG'S ORDER ("Save the resume file asap!"). Read this block, then the ones below it (newest first). All workflows run inside Claude session f1bdbd78 (Fable driving); a new session reads their journals (`C:\Users\grego\.claude\projects\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\subagents\workflows\<run>\journal.jsonl`), the scratch files under `...\f1bdbd78-...\scratchpad\`, and `sacct`. Nothing committed by any workflow; MBR/ stays untracked.

**TIMING DONE (wmem2cm7r, 2026-09-27):** no per-animal date for rats 64-71 in any figure, legend, table or supplement of the lab's papers (JNS 2026 published + submitted, Acta 2021 + video, six abstracts, eleven same-lab U251 papers) or in any lab file; the RNA-seq rats appear only in aggregate (JNS Fig. 8, S6). Hard bound: rats 66-71 harvested before 7 Apr 2022 09:22 (chip scan). Group level: SNO 2023 MODL-14 says recurrences were taken 2 weeks after LITT (about 4 weeks after implantation); no source gives the primaries' harvest age (about 2-4 weeks). "Raj" = Tavarekere N. Nagaraja. Applied: SAMPLE_KEY.md "What is not known" + dates section (MODL-14, 3-4 week survival, 5x10^4 vs 5x10^5 inoculum discrepancy); CONTRALATERAL_CONTROL_PRIOR_ART.md Acta pages 3455-3463; SLIDES/03_build_deck.py JNS 2026;145:364-377. Answer: scratch `u251_timing/ANSWER.md`, `TIMING.csv`, `PAPERS.md`, `FIGURE_READ.md`. Per-animal dates would come only from the lab: MRI folders `YYYYMMDD_HHMMSS_IL<rat>`, Visualase logs, IACUC #1509 log.
**nf-core/methylarray (wkj7t0b8b) DONE as a workflow:** first grid run 40480122 dropped all six tumours at SeSAMe QC (EPIC v1 quality mask counted as failures); patch P6 counts only unmasked probes (patch sha256 0ef3365e...); resubmitted as **40480724** (started 23:59); results ~00:20-00:45, relabel_null ~01:15-02:15.

**RESUMED 23:45 EDT (machine clock) on Opus 5.5; agents work again.** The limit actually hit at ~22:11. New task ids (same run ids): SPREAD v2 wjy1l2zie, XVAL wbyvmmgyq, SCIENCE whvuo0fak, CURVE wc2ib5g82, POWER wmvpy6nwg, BENCH wsg4cdgp5, multi-omics wi5bjtzih, nf-core/methylarray wkj7t0b8b, TIMING wmem2cm7r.

**SESSION LIMIT HIT ~23:55 EDT (resets 03:00 ET).** Agents started after that failed. Resume each run after 03:00 in the same Claude session with `Workflow({scriptPath: "<workflows\scripts\<name>-<run>.js>", resumeFromRunId: "<run>"})` (scripts under `C:\Users\grego\.claude\projects\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\workflows\scripts\` or the `C--...` variants; finished agents replay from cache). Known state at the limit:
- **SPREAD v2 (wf_e44b06bf-0e6)**: research, EVIDENCE_existing.md, DESIGN_IMPLICATIONS.md, CLONAL_LIT.md and PREREG_spread.md DONE (111 of 115 agents); the two analysis jobs, the write-up and the verify FAILED at the limit. **PREREG_spread.md hashes, recorded here as the document's section 8 requires: frozen base (sections 0-7) sha256 `454a24f3b3d6c438e2be6f0f498f5b11e4d90f2b7424204e9f6e40bae78be47d` (41,412 B, the digest DEVIATIONS.md and both RESULTS files cite); current file with the one appended section-8 entry of 23:15 (DiG unobtainable -> Ivy GAP leading-edge contrast is the primary directional statistic) sha256 `96131616a574eb4040edab18e934ee1fc73c95b277cae266a2dc8d70da7f2dd9` (43,113 B), byte-identical in `ANALYSIS/spread/PREREG_spread.md` and on the grid at `/rs/rs_grp_oschome/go2432/u251_science/spread/PREREG_spread.md`; sections 0-7 verified byte-identical between the two by diff. Parts 2 and 3 were prospective at freezing; E2 and E3 are labelled ALREADY OBSERVED in section 0.1 and carry no pre-registration credit. E1, the only pre-specified gate that can retire the artefact hypothesis, is still in flight (40480691 RUNNING, 40480692/40480693 PENDING behind it), so its result has not been examined. Facts found: the human count matrix was quantified from graft+both reads (both = 53.6 % of N269B's human input; > 90 % of the floors'), so the floor correction must model the both class; N168B's floor RNA is itself Y-bearing (92 chrY per million graft pairs, 0.59x cores: hopping, handling or a few cells; unexplained) while IL64B's is Y-poor (true artefact); index hopping cannot be excluded (one lane with six Y-bearing cores + C2B); the contralateral piece was dissected from the same tumour-bearing slices with the same blades on the same slide, so carry-over is permitted by the method and nothing in hand separates it from spread; DNA/RNA fraction ratio 4.0 at N2 vs 1.3-2.0 in cores (CURVE to settle); the literature has no false-human-call floor, and the public tumour-free rat brain corpus for it is GEO **GSE53960** (rat BodyMap, 32 F344 brain total-RNA libraries).
- **XVAL (wf_0186e254-651)**: INVENTORY_xval.md, PREREG_xval.md (sha256 140fdd44c15143ddfc70d5852fe9aa7d69e7a9897b2047ebc19b2afaf8568c75; copies on the grid and in ANALYSIS/xval/), reference annotations hashed under `u251_science/xval/ref/`, reads step DONE (no SRA fetch needed); array-derived tables, validate, write and verify FAILED at the limit. Corrections to earlier notes: the human STAR BAMs exist (45.5 GB) AND the xengsort graft/host FASTQs survive in `~/u251-xenograft-murine-RNASeq-study/ANALYSIS/sorted_fastqs/` (106 GB) with the index (`xengsort_index_clean/`, 17.9 GB, k = 25) and both reference FASTAs; only the ambiguous+neither bins were discarded. SRA runs SRR39547363-72 (spot counts equal the xengsort totals exactly). HOME is at 96.8 % of quota: nothing new goes there; use /rs.
- Other runs (SCIENCE, POWER, CURVE, BENCH, multi-omics, nf-core/methylarray, TIMING): state unknown at the limit; read their journals first.

| task id (run) | what | where |
|---|---|---|
| wrfcuzrxn (wf_e44b06bf-0e6) SPREAD v2 | Greg's hypothesis: contralateral U251 cells distinct from their core (H1-A) and closer to the recurrences (H1-B); STATE level with the floor-aware carry-over null + 20-labelling null; CLONE level (N2 CNV events, RNA variants); depth-matched CIBERSORT/Neftel states; design implications; manuscript section | grid `/rs/rs_grp_oschome/go2432/u251_science/spread/`; scratch `spread/`; `MBR/nature_draft/SPREAD_RESULTS.md`, `SECTION_spread.md` |
| wzy3qk6h5 (wf_0186e254-651) XVAL | matched array + RNA cross-validation: rs-probe genotypes vs RNA alleles (error model), CNV -> expression dosage incl. N2-only losses, promoter methylation vs expression, DNA support for the RNA lists; human BAMs exist so no SRA fetch unless missing | grid `u251_science/xval/`; `MBR/nature_draft/XVAL_RESULTS.md`; code `ANALYSIS/xval/` |
| wngku5ojo (wf_fce2e1f1-dba) SCIENCE | nested deep research -> three proposals -> one PREREG -> RNA/DNA jobs -> manuscript draft -> hostile referee -> revision | `MBR/nature_draft/manuscript.*`, `REFEREE_RESPONSE.md`; scratch `science/`; grid `u251_science/{rna,dna}/` |
| we2kzevy2 (wf_6afc3bcc-1d6) POWER | pre-stated targets, pilot variances, rats-needed chart (SchwingNet style), SCALED_EXPERIMENT.md + closing section; add the barcode clonal-tracking target (GOAL section) | `MBR/nature_draft/SCALED_EXPERIMENT.md`, `figures/fig_rats_needed.*`; grid `u251_science/power/` |
| wzvcizag3 (wf_480455e9-000) CURVE | in-silico human:rat titration curve validated on real human:mouse mixtures (GSE310817, GSE299969); our eight arrays placed on it | grid `u251_meth/curve/results/` (three PNGs); scratch `titration_curve/RESULT.md` |
| w1cox6j4w (wf_d243847e-fb8) BENCH | every published host-DNA method with rat, judged by PREREG checks C0-C3; pipeline_v2 after the bespoke run | grid `u251_meth/pipeline_v2/`, `u251_meth/ref/`; scratch `meth_bench/` |
| wc03287jk (wf_d598cc75-610) multi-omics | nf-core candidates evaluated; at most three submitted; `ANALYSIS/MULTIOMICS_PLAN.md` | grid `u251_multiomics/` |
| w3fa7sf50 (wf_2e081152-c82) nf-core/methylarray | PR #41 SeSAMe/limma/DMRcate arm; first grid run 40480122 failed at 5.5 min (SeSAMe dropped all six tumours: EPIC v1 quality mask > 0.10); fixed (patch P6) and resubmitted as **40480724** at 23:58, see the 23:58 block | `u251_meth/nfcore_methylarray/`, log `u251_meth/logs/nfcore_methylarray_40480724.out` |
| wlw03fxg7 (wf_1ce53ed2-059) TIMING | per-animal implant/LITT/harvest timing from Raj's papers' figures, lab zips, IDAT headers, GEO | scratch `u251_timing/ANSWER.md`, `TIMING.csv` |
| DONE 40480027 bespoke methylation | REPORT.md, summary.json; N2 carries U251N DNA (genotype 56/59, male Y, U251 CNV pattern); MGMT all M; Q5 no DMPs after fraction; Q4 no shared CNV change | `/rs/rs_grp_oschome/go2432/u251_meth/results/` |

Greg-side: send the GEO update (`GEO/seq_template_u251_FILLED.xlsx`, edited 22:55; confirm the sequencing core first); Raj: leftover DNA (IL64B/N168B for an EPIC v2 rat-brain reference; any sample for Alu qPCR), MRI volumes, harvest dates, Visualase/IACUC logs.

## 2026-09-26 22:55 EDT: GEO template edited; SPREAD v2 (Greg's H1-A/H1-B) resumed; DNA-RNA cross-validation launched; raw FASTQs are gone from the grid but BAMs exist

- **GEO GSE338105 metadata**: `GEO/seq_template_u251_FILLED.xlsx` edited (backup `...BEFORE_SUMMARYFIX.bak.xlsx`): title U251 -> U251N (matches the live record); summary's last sentence (translational shutdown / ISR / "Garofano Mitochondrial subtype") replaced by the honest design statement (lower tumour fraction in recurrences; relabelling null; ribosomal-protein block; contralaterals and culture; N269B holds tumour cells); xengsort step names GRCh38/GENCODE 44 + mRatBN7.2 (Ensembl 110). Raw and processed files unchanged (MD5s in the sheet = `GEO/md5_checksums.txt`), so it is a metadata-only update. Still Greg's: confirm the sequencing core (template says USC Molecular Genomics Core); send to GEO (geo@ncbi.nlm.nih.gov with the accession, or the portal's update). Live record still says "matched", "procedural / failed-graft control", "tumor-free" for N168B/N269B: wrong; the template fixes all of it.
- **Raw FASTQs**: not on the grid (DATA/RNASEQ/RAW holds only the Google-Drive download scripts); rat-sorted reads at `/rs/rs_grp_oschome/go2432/u251_host/fastq/`; **human STAR BAMs exist for all ten libraries** (`ANALYSIS/results_human_final/star_salmon/<lib>.markdup.sorted.bam` on the grid clone) and rat BAMs at `u251_host/results_rnaseq/star_salmon/`. SRA holds the ten runs (GSE338105) if a re-sort is ever needed; the XVAL workflow fetches from SRA only if BAMs are missing.
- **SPREAD v2** (task wrfcuzrxn, run wf_e44b06bf-0e6 resumed with the rewritten script): Greg's hypothesis (22:45) pre-registered as H1-A (contralateral U251 population distinct from its own core) and H1-B (closer to the recurrences than the primaries), STATE level (floor-corrected profiles; NULL-A = IL69B subsampled to N269B's U251 share PLUS the assignment-floor profile, 1,000 draws, seed 20260926; NULL-L = all 20 labelled 3-v-3 assignments, minimum p 0.05; CIBERSORT/Neftel state composition via the repo's guarded runner) and CLONE level (N2-only CNV events vs recurrences; RNA-derived variant sharing from the human BAMs); STATE-level rejection is consistent with plasticity as well as selection — only DNA-level sharing supports "clonal evolution". Plus the scaled Gawad-style test: lineage-barcoded U251N through LITT (GOAL section).
- **XVAL** (task wzy3qk6h5, run wf_0186e254-651): matched array + RNA of the same eight samples: (1) array genotypes at the 59 rs probes and SNP-in-probe CpGs vs RNA allele counts from the BAMs -> concordance and the RNA error model (the null for RNA-only variants; germline sites validate the pipeline, not somatic subclones); (2) CNV -> expression dosage, incl. the N2-only losses (3p, 4p, 4q, 12p) against N269B's floor-corrected expression (real subclone shows dosage; noise does not); (3) promoter methylation vs expression, within sample and delta-vs-delta with fraction; (4) DNA support for the existing RNA lists. Outputs `/rs/rs_grp_oschome/go2432/u251_science/xval/`, code `ANALYSIS/xval/`, note `MBR/nature_draft/XVAL_RESULTS.md`, scratch `xval/`.
  - **XVAL reads step (23:23 EDT):** the ten human BAMs were verified on the grid (sizes = INVENTORY_xval.md), so **no SRA fetch, xengsort re-sort or STAR job was submitted**. One bookkeeping job **40480668 `xval_bamchk`** (reqp, ~20-40 min; `ANALYSIS/xval/00_bam_check.sbatch`, grid copy `u251_science/xval/code/`) writes `xval/rna/bam_manifest.tsv` (+ `rna/bam_check/<lib>.{header.sam,idxstats,flagstat}`, `rna/xengsort_logs/`), checks each BAM against the pipeline's own flagstat and against the xengsort `graft + both` pair count, and prepares `xval/ref/GRCh38.primary_assembly.genome.fa(.fai)` + GTF after GENCODE md5 checks; result in `xval/rna/BAM_CHECK_STATUS.txt` (`OVERALL PASS|FAIL`), log `xval/logs/bamchk_40480668.out`. Recorded in `u251_raw/JOBS.txt` and scratch `xval/DEVIATIONS.md` (procedural only). `MBR/nature_draft/` added to `.gitignore`.

## 2026-09-26 22:20 EDT (grid clock): bespoke methylation run DONE (40480027); contralateral-spread evidence in hand; SPREAD workflow launched

Greg (22:15, Fable now driving): "SAY SOMETHING SIGNIFICANT ABOUT SOMETHING HOT ... EXTRACT THE EXACT BEST SIGNAL FROM MY MIXED RAT/TUMOR DNA AND SEE IF YOU CAN USE THE SIGNAL AS AN INDEPENDENT VERIFICATION OF SPREAD OF TUMOR ACROSS THE CORPUS CALLOSUM ... UNDERSTAND THE IMPLICATIONS OF THE EXPERIMENTAL DESIGN".

**Results in hand** (`/rs/rs_grp_oschome/go2432/u251_meth/results/REPORT.md`, `tumour_fraction/tumour_fraction.json`; RNA from `ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv`):
- **N2 = N269B (rat 69, contralateral) carries U251N DNA**: rs-probe genotype concordance with C2B 56/59 (95 %; tumours 58/59; permutation null 37 %, p = 0); chrY DNA dosage 0.88 of C2B with 64 % of Y probes detected (host female); human CNV profile r 0.72 with C2B, z 8.4 against the null, and N2 reproduces U251's arm pattern (1p, 7p, 9p, 15q, 18p gain; 10p, 11q, 13q, 18q loss). DNA fraction: two-mode 19.8 % (17.7-22.1), other estimators 17-47 %; absolute calibration pending from the CURVE workflow (real mouse titrations). MGMT-STP27: all eight "M" (U251 is MGMT-promoter methylated in the literature: outside check passes).
- **RNA chrY per million human counts**: C2B 256; cores 177-209; **N269B 112**; N168B 23.5; IL64B 2.9. So xengsort's "human-only" calls in a tumour-free rat sample (IL64B, 1.36 M counts) carry no Y: they are the assignment floor, not tumour. Correcting by chrY ratio: about two thirds of N269B's 4.16 M human counts are U251 (true U251 RNA share ≈ 3 %, chrY ratio then ≈ 175/M, matching the cores); N168B ≈ 10 % of its 1.70 M (≈ 0.06 % of reads): a trace, 8x the floor. Cores: the floor is a few % of their human counts.
- **Q5 (primary v recurrent DMPs)**: 0 at FDR < 0.05; the true split ranks 1 of 10 for p < 1e-3 counts (3,975) under arm alone but 326 (rank 4) with U251 fraction as covariate: the arm "difference" is mostly tumour fraction. **Q4**: no arm-level CNV change shared by all recurrences and no primary; the raw relabelling rank 1 becomes rank 2 after amplitude scaling (fraction again). Rat-probe list 85,764.
- **What this settles / does not**: DNA + RNA independently show U251N cells (genotype, male Y, copy-number pattern) in the contralateral hemisphere of one untreated animal. It does NOT separate in-vivo spread from dissection carry-over; the only evidence that can is the transcriptional state of N269B's human reads against its own core subsampled to the same depth (pre-registered in the SPREAD workflow). Both contralaterals are untreated rats: nothing about LITT and spread. Open point: DNA ~20 % vs true RNA ~3 % for N2; the two-mode estimator may inflate at low fraction (noise floor); the CURVE will show.
- **SPREAD workflow wo0cec1iy (run wf_e44b06bf-0e6)**: nested deep research (midline crossing rates/routes in U251 and GBM xenografts; carry-over controls; infiltrating-cell signatures; recurrence beyond the treated volume; pseudoprogression after LITT; xengsort floors); PREREG_spread.md (sha256) with E1 xengsort floor on public tumour-free rat brain RNA-seq, E2 chrY, E3 DNA identity, E4 DNA-vs-RNA, E5 depth-matched carry-over test (IL69B human counts multinomially subsampled to N269B's size, 1,000 draws, seed 20260926; fixed edge/core signatures), E6 host side; outputs `/rs/rs_grp_oschome/go2432/u251_science/spread/`, code `ANALYSIS/spread/`, write-ups `MBR/nature_draft/SPREAD_RESULTS.md`, `SECTION_spread.md`; scratch `spread/` (EVIDENCE_existing.md, DESIGN_IMPLICATIONS.md, PREREG_spread.md).

## 2026-09-26 22:00 EDT (grid clock): HANDOFF — everything running or waiting (Greg: "WRITE EVERYTHING YOU ARE CURRENTLY WORKING ON/WAITYING ON TOA RESUME CHECKPOINT FILE NOW"; Fable takes over)

**Greg's standing asks, newest first**
1. "please do your best to extract real biological insight, using validated, tested, reproducible, transparent approaches, THAT A NATURE PAPER OR CELL PAPER WOULD ACCEPT ... IDENTFIY SOME RELEVANT PROBLEM IN NEUROONCOLOGY AND ANSWER IT A LITTLE WITH THIS DATA LIKE PSEUDOPROGRESSION OR MIGRATION OR RECURRENCE OR CLONAL DYNAMICS OR LIQUID METASTASIS ... /deep-research ... write it to a manuscript". Workflow **SCIENCE** below. Be honest in the text about what n = 3 v 3 can carry.
2. Methylation: apply the PUBLISHED host-DNA methods with rat substituted, judged against a pre-registered correctness check ("just try all these"); never present our own estimator as a new method ("i dont like this idea"). nf-core wherever it exists (nf-core/methylarray exists, unreleased).
3. "are you sure its not in one of the figures from raj's papers?" (per-animal harvest timing). Workflow **TIMING**.

**GOAL of the paper and of everything below (Greg, 22:30: "detailing an experiment that should scale and give results that close a gap in the literature with whatever you found, similar to how I calculated power calculations for detecting rare classes and then ended my schwingnet paper teasing the reader about how it would work with more data")**
- The paper reports what this pilot measures honestly, then ENDS on a pre-specified scaled experiment with a "rats-needed" chart in the SchwingNet form. SchwingNet Figure 5a was the Wilson 95 % lower bound against the number of test persons, with pre-stated target lines (`spinesurg-ct-nnunet/SchwingNet/main.tex` Figure 5; `results/power/power.json`; `tools/fig_power.py`). Here: for each pre-stated target, the rats per arm needed, from pilot variances, drawn as curves with the pilot's own uncertainty.
- Gaps the scaled experiment would close (from today's findings; the SCIENCE deep research will sharpen them):
  - **G1, lesion composition after LITT.** What fraction of a regrowing lesion is tumour and what is host tissue, measured exactly by species in RNA and DNA; which host programmes make it up. This bears on the treatment-effect / pseudoprogression problem. Pilot: human reads 51 % in primaries vs 35 % in recurrences; array signal lower too.
  - **G2, the cells beyond the ablation zone.** How often tumour cells reach the contralateral hemisphere, whether they carry the core's copy-number profile, and whether LITT changes this. Pilot: 1 of 2 sampled contralaterals (rat 69), both from the primary arm.
  - **G2/G3 sharpened (Greg, 22:45): "the contralateral tumor cells (presumed migratory) are a distinct population from the primary biopsied cells, and the migratory cells have more in common with the recurrent tumors than the primaries ... reproducible clonal evolution ... prove it a priori with a defined null hypothesis".** Pilot test: SPREAD workflow v2 (H1-A distinct from own core; H1-B closer to recurrences; STATE level with a depth-matched-plus-floor null and a 20-labelling null; CLONE level from N2's CNV events and RNA-derived variant sharing; STATE-level rejection cannot separate clonal selection from plasticity). Scaled test, the Gawad-style design: lineage-barcoded U251N -> implant -> LITT -> barcode sequencing of core, contralateral and recurrence in every rat; null = the recurrence's barcodes are a random sample of the primary's; the hypothesis predicts enrichment of contralateral (invasive) barcodes in recurrences. Thousands of clones per animal give within-animal power; the rats-needed curve for this target goes on the chart (POWER workflow to add it, or Fable adds it when POWER finishes).
  - **G3, clonal selection by thermal ablation.** Copy-number and methylation evolution of recurrences against primaries and the parental culture. Pilot: bespoke CNV and methylation outputs, pending.
  - **G4, the methods gap.** No human:rat DNA mixture has ever been run on a human methylation array, and no RNU rat DNA is on one. A real titration series (U251N DNA into RNU brain DNA at 0, 1, 2, 5, 10, 20, 50, 100 %) on EPIC v2, plus pure RNU brain, is a citable resource in itself. It also calibrates every fraction estimate.
  - **G5, controls.** The contralateral hemisphere is not naive. A naive brain and an ablation-only (sham implant + LITT) arm separate injury from recurrence biology.
- Design skeleton, to be sized by the power workflow below:
  - **Arms:** (i) tumour, no LITT, time-matched harvest; (ii) tumour + LITT, harvest at MRI regrowth; (iii) tumour + LITT, early harvest in the treatment-effect window; (iv) sham implant + LITT; (v) naive.
  - **Per rat:** lesion core, rim and contralateral hemisphere; MRI volumes and harvest dates recorded.
  - **Per sample:** AllPrep DNA+RNA from the same piece; adjacent sections stained with a human-specific marker as the tissue-level ground truth for tumour share and invasion; human Alu qPCR; EPIC v2; total RNA-seq.
  - **Calibration:** the titration series (G4).
- Every target is pre-stated and fingerprinted before the pilot numbers feed it. Pilot effect sizes from n = 3 are inflated, so design for the conservative end of the pilot interval and show the whole curve.

**Background workflows (Claude session f1bdbd78; transcripts `C:\Users\grego\.claude\projects\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\subagents\workflows\wf_*\journal.jsonl`, one result line per finished agent; resume works only inside that session, so a new session reads the journals and scratch files instead). Scratch root: `C:\Users\grego\AppData\Local\Temp\claude\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\scratchpad\`.**

| task id | what | state at 22:00 | where results land |
|---|---|---|---|
| wpk9ifkvj | bespoke methylation pipeline (Nextflow, `ANALYSIS/methylation/`, grid `u251_meth/pipeline/`) | its run is grid **40480027** (u251_meth, RUNNING since 21:35; the earlier 40479710 was cancelled by the workflow after a fix). Published so far: qc (20:57), species_crosshyb, mgmt, cnv, tumour_fraction (21:36-21:42); DMP_DMR and REPORT pending | `/rs/rs_grp_oschome/go2432/u251_meth/results/` (REPORT.md, summary.json), log `u251_meth/logs/pipeline_40480027.out` |
| w3fa7sf50 | nf-core/methylarray standard arm (block below) | grid **40480122 FAILED after 5.5 min**; read `u251_meth/logs/nfcore_methylarray_40480122.out`; the workflow may already be fixing it | `u251_meth/nfcore_methylarray/`; README correction of the old "nf-core has no pipeline" claim is part of this workflow |
| wc03287jk | multi-omics plan: evaluate rnavar/rnadnavar/oncoanalyser, rnafusion/rnasplice/circrna, RNA x methylation integration, invasion/microenvironment; submit at most 3 released nf-core pipelines | running | `ANALYSIS/MULTIOMICS_PLAN.md`, `ANALYSIS/multiomics/<pipeline>/`, grid `/rs/rs_grp_oschome/go2432/u251_multiomics/`; scratch `multiomics/` |
| w1cox6j4w | rat host-DNA methods benchmark: probe removal E1 Guilhamon, E2 Ebata, E3 SeSAMe pOOBAH, E4 SeSAMe EPIC.addressSpecies rat mask, E5 EPIC-probe alignment to mRatBN7.2; fraction F1 Zhou M2, F2 Zhou M1 colour-switch rebuilt for EPIC, F4 chrY, F5 InfiniumPurify (negative control); bespoke E6/F3 comparison only. Checks C0-C3 and the selection rule fixed in PREREG.md (sha256 recorded) BEFORE results | running | scratch `meth_bench/` (METHODS.md, DATA.md, PREREG.md, INTERFACE.md); grid `u251_meth/pipeline_v2/`, `u251_meth/ref/` (GSE174568 IDATs, rat alignment); local `ANALYSIS/methylation_v2/` |
| wzvcizag3 (run wf_480455e9-000) CURVE | Greg: "synthetically mix pure rat and pure u251n and plot a curve ... and see where my samples lie" (after pointing to GSE299969). Pre-registered (PREREG_curve.md + sha256): Zhou 2022 Method-2-style within-array statistic (human-only vs shared probe signal, probe sets from sesameData EPIC.addressSpecies); step 1 validates in-silico mixing against REAL human:mouse mixtures GSE310817 (Zhou 2022, EPIC v1, two batches, pure ends in each) and GSE299969 (Ebata 2026, EPIC v2 despite GEO's label); step 2 rat curve = C2B (pure U251N) mixed in silico with each of 10 GSE174568 rat arrays; step 3 reads our eight arrays' human DNA fraction off it, next to RNA human reads. Independent of the w1cox6j4w benchmark by design (separate code and dirs) | launched 22:20 | grid `/rs/rs_grp_oschome/go2432/u251_meth/curve/` (ref/, results/ with three PNGs); code `ANALYSIS/methylation_curve/`; scratch `titration_curve/` (PREREG_curve.md, RESULT.md) |
| we2kzevy2 (run wf_6afc3bcc-1d6) POWER | the GOAL section's scaled experiment: pre-stated targets T1-T7 fingerprinted before pilot numbers (TARGETS.md), pilot variances from this data, rats-needed curves (pwr/Welch, Wilson/exact, PROPER or ssizeRNA with pilot dispersions, pwrEWAS), sized on the conservative end of the pilot interval, independent recompute of key n, then SCALED_EXPERIMENT.md + the manuscript's closing section + an update of the GOAL section here | launched 22:35 | `MBR/nature_draft/SCALED_EXPERIMENT.md`, `CLOSING_SECTION.md`, `figures/fig_rats_needed.{png,pdf}`; grid `/rs/rs_grp_oschome/go2432/u251_science/power/` (TARGETS.md, pilot/pilot.json, power.json); code `ANALYSIS/power/`; scratch `scaled_power/` |
| wlw03fxg7 | TIMING: every Nagaraja-lab U251/LITT paper, figure by figure, plus lab zips, IDAT headers, GEO GSE338105, MBR draft, for per-animal implant/LITT/harvest timing of rats 64-71 | running | scratch `u251_timing/ANSWER.md`, `TIMING.csv`, `FIGURE_READ.md`, `papers/` (PDFs: never into git) |
| wngku5ojo (run wf_fce2e1f1-dba) SCIENCE | nested deep research on candidate questions (treatment effect/pseudoprogression, invading cells beyond the ablation zone, clonal/CNV evolution, remote host response, MGMT, cell-state shift) + data inventory -> three proposals (lesion composition by species; rat-69 contralateral invaders vs core; recurrence evolution with tumour share modelled) -> editor merges into one PREREG (sha256; grid copy `/rs/rs_grp_oschome/go2432/u251_science/PREREG.md`) -> RNA and DNA analysis jobs (code `ANALYSIS/science/`, outputs `/rs/.../u251_science/rna|dna/`) -> manuscript draft with \PENDING{} placeholders -> hostile referee -> revision + its own RESUME block | launched 22:05 | `MBR/nature_draft/` (manuscript, figures/, references.bib, REFEREE_RESPONSE.md; MBR is untracked, never commit); scratch `science/` (INVENTORY.md, PROPOSAL_*.md, PREREG.md, DEVIATIONS.md) |

**Facts established today that the write-up needs**
- Rat DNA on EPIC v1 exists publicly only as GEO **GSE174568** (Arneson et al., Nat Commun 2022, doi 10.1038/s41467-022-28355-z): ten EpigenDx rat DNA standards (0-100 % methylated), strain not stated, plus three human and fifteen mouse standards. GEO has no RNU/nude-rat methylation data and no rat-host xenograft methylation data on any platform (searched 2026-09-26). The study rats are RNU/RNU (ABSTRACT.md).
- SeSAMe ships `EPIC.addressSpecies` with `inferSpecies()` / `updateSigDF(species=)`; whether rat is in its list is being checked by w1cox6j4w.
- Cited-tools deep research (`tasks\wx7k31mge.output` in the session temp dir): xengsort is the benchmark-recommended read splitter (Bhandari 2025, npj Precis Oncol) but validated only with mouse hosts; host-probe removal is published only for mouse (Guilhamon 2014 Genome Med; Ebata 2026 Epigenomics, EPIC v2); human-fraction by intensity: Zhou 2022 Cell Genomics (Method 2 transfers; Method 1 needs the mouse array); InfiniumPurify invalid here (human-normal model, needs >= 20 samples). Refuted, do not cite: Needhamsen 2017 numbers; "Guilhamon recommends >= 15 PDX per group". Unanswered (not "none exists"): integration-method minimum n, Heidelberg classifier / MGMT on U251, non-nf-core methylation pipelines.
- Methylation QC (bespoke run): every array < 0.2 % probes with detP > 0.01; median total signal C2B 9,811, tumours 6,247-8,597, N2 3,335; SeSAMe pOOBAH failures C2B 1.5 %, tumours 1.8-5.2 %, N2 16 %; minfi sex M for all eight. Rat-probe exclusion list 85,764 of 866,238 (30,308 by intensity, 11,843 by sequence, 43,613 not assessable); flagged probes align to rat (Rnor_6.0) at 41 % vs 1.5 % of unflagged (odds ratio 46).
- The design figure and slide 3 now show the methylation arrays (`SLIDES/12_design_figure.py` reads `ANALYSIS/methylation/assets/samplesheet.csv`; deck rebuilt, only slide 3 changed; previous deck backed up in scratch `deck_check/`). Uncommitted.

**Limits to state in every write-up** (told to Greg 21:50): 3 v 3, unpaired; recurrences hold less tumour (human reads mean 51 % primaries vs 35 % recurrences; array signal lower too) — the largest confound; per-animal timing unknown (TIMING workflow checking); planned 4 v 4 became 3 v 3 (IL64B no tumour, rat 65 no library); no naive rat brain and no ablation-only control; contralateral tissue only from rats 68 and 69, both primary arm, only 69's holds tumour; bulk tissue, no spatial data; two species in every sample (xengsort rat-host unvalidated; 10 % of array probes lost; no RNU rat reference on the array); no DNase step in RNA extraction; no replicates; EPIC v1 discontinued; culture passage unrecorded.

**Waiting on Greg / Raj**
- Raj: leftover DNA from IL64B or N168B (one EPIC v2 array = exact RNU rat-brain reference); any DNA for a human Alu qPCR (direct human fraction); MRI folders (tumour volumes), harvest dates, Visualase logs, IACUC log.
- Greg: GEO GSE338105 update (remove the "Garofano Mitochondrial" claim); confirm the sequencing core; whether to commit/push SLIDES and ANALYSIS changes (public repo: never MBR/, DATA/raw_lab/, IDATs, papers).
- Older open items: RUVSeq 2-control rerun (40476123; check sacct) then rerun 06/07 and replace "rerun pending" on B1; host verify 40475691 COMPLETED, B7 exact numbers still to be written in.

## 2026-09-26 21:50 EDT (grid clock): PENDING nf-core/methylarray standard arm (Greg: "i find it hard to believe that no nfcore pipeline exists for this data?")

> SUPERSEDED: 40480122 below FAILED at 21:56 (SeSAMe sample QC dropped all six tumours); the timings below never happened. Fixed and resubmitted as 40480724, see the 23:58 block.

| job | what | ETA | output |
|---|---|---|---|
| 40480122 | `run_methylarray.sbatch` (reqp/requeue, 8 cpu, 64G, 8h): zip sha256 + 16 IDATs (md5 checked against the bespoke STAGE table), frozen rat list, clone + verify PR #41 commit, apply patch, R library + hub cache, `-profile test` (recorded only), the real run on the six tumours, then the nine other 3-v-3 relabellings with `-resume` | wall ~2.5-3.5 h: main results ~23:00-23:30, relabelling ~00:30-01:30 | `/rs/rs_grp_oschome/go2432/u251_meth/nfcore_methylarray/` (`run_status_40480122.txt`, `results/`, `results/u251_checks.txt`, `results/relabel_null.tsv`), log `u251_meth/logs/nfcore_methylarray_40480122.out` |

- Greg was right. **nf-core/methylarray exists** (in development since Aug 2024, no release; `dev` HEAD `d34f531`,
  2026-07-27). The bespoke README section and the `main.nf` header claiming otherwise are corrected (local and grid
  copies; comments only, no logic). `run_pipeline.sbatch`'s nf-core check was left as is: it prints every match,
  methylarray included.
- Pinned: PR #41 (`schumz:feat/sesame-limma-dmrcate`, SeSAMe + limma + DMRcate) head
  `3a7ff00dcd71aab25efd57eb659dbd40370a3ffa`, plus `patches/pr41_3a7ff00_epicv1.patch` (sha256 `52e2e998...`). The patch
  covers EPIC v1 annotation, the v1 annotation package plus a guard at 90 % of probes, `cpg.annotate` EPICv1, a failed
  contrast written as NA instead of 0, and a rat `--exclude_probes`. Container: the bespoke Bioc 3.22 sandbox plus a
  separate R library (the PR's `quay.io/nf-core/methylarray:1.0.0dev` was never published). `dev` was not used: it
  hard-codes quantile normalisation, uses ChAMP `450K` on EPIC, and skips DMP/DMR without saying so.
- Non-default settings: `annotation ilm10b4.hg19` and `genome_build hg19` (EPIC v1); `exclude_probes` = frozen copy
  of the bespoke rat list; `keep_groups PRIMARY,RECURRENT`; `cofactors ","` (empty design `~0+group`);
  cell composition, covariate PCA / ChAMP SVD and coMethDMR off; six tumours only (N2 and C2B out: N2's 16 % pOOBAH
  failure would drop it or strip probes). Read only the `norm_only` outputs, because one chip gives one Plate.
- The rat list changed under us: 73,921 probes (bespoke 21:06 run) became 85,764 after the bespoke rerun 40480027
  (started 21:35 from another session; it supersedes 40477975 in the block below). The job freezes the list it uses,
  logs its sha256 and says at the end whether the bespoke list still matches. Compare Q5 between the arms only if it does.
- `-profile test` cannot pass at this commit: its test data branch does not exist in nf-core/test-datasets (404). The
  job records the failure and carries on.
- When it finishes, read in this order: `run_status_40480122.txt` (any `STOP`), `results/u251_checks.txt` (detected
  array, probes in annotation, rat probes removed, every DMP/DMR `ERROR`; errors in `plate_in_model` and `combat` are
  expected), `results/dmp_limma/DMP_limma__norm_only/`,
  `results/dmr_dmrcate/DMR_DMRcate_EPICv2__norm_only/`, `results/relabel_null.tsv`, then compare with the bespoke
  `results/dmp_dmr/`. Read Q5 as "recurrent vs primary, confounded with U251 DNA fraction".
- Local copy: `ANALYSIS/methylation/nfcore_methylarray/` (untracked; no IDATs or lab files in it).

## 2026-09-26 20:37 EDT (grid clock): PENDING methylation pipeline (EPIC chip 205648300021)

| job | what | ETA | output |
|---|---|---|---|
| 40477399 | `run_setup.sbatch`: pinned Bioconductor 3.22 / R 4.5.2 SIF (pulled 20:21) + R library into `rlib/`; then sesameData, CopyNeutralIMA, DMRcate caches and mgmtstp27 0.8 (tarball URL checked: HTTP 200). At 20:36: 169 packages in, installing the EPICv2 annotation, CPU busy | ~21:00-21:15 | `logs/setup_40477399.out`, `logs/setup_versions/` |
| 40477975 | `run_pipeline.sbatch` (afterany 40477399; reqp/requeue, 8 cpu, 64G, 6h): STAGE, QC_NORMALISE, SPECIES_CROSSHYB (Q6), CNV (Q4), TUMOUR_FRACTION (Q1, Q2), MGMT (Q3), DMP_DMR (Q5), REPORT | start after setup + queue; ~45-90 min run; results ~22:00-23:00 | `/rs/rs_grp_oschome/go2432/u251_meth/results/` (`REPORT.md`, `summary.json`), `logs/pipeline_40477975.out` |

- 40477200 (first setup) failed at the image pull: `mksquashfs` not on PATH. Fixed in the resubmission.
- Pipeline copy lives at `/rs/.../u251_meth/pipeline/` (local: `ANALYSIS/methylation/`, untracked); raw zip sha 1107d481... in `u251_meth/raw/`.
- The chip's eight: IL66B, IL67B, IL68B, IL69B, IL70B, IL71B, N2 (= N269B, contralateral of rat 69, carries tumour) and
  C2B (U251N culture). No tumour-free brain (no IL64B, no N168B) is on the chip.
- When it finishes: read `results/REPORT.md` and `summary.json`; check `pipeline_info/` for a failed process first
  (afterany starts the pipeline even if setup failed).

## 2026-09-26 18:44 EDT (grid clock): PENDING host verification for deck backup B7 (rat host DE, Recurrent vs Primary)

| job | what | ETA | output |
|---|---|---|---|
| 40475691 | `run_host_verify.sbatch` (reqp/requeue, 4 cpu, 24G, 6h). Step A in the pipeline's own DESeq2 1.34.0 + ashr container: R0 reproduce the 38-gene table; R1 leave one tumour out (6 fits, the 38 genes per fit, Mmp13); R2 all 10 balanced 3-v-3 relabellings; R3 NL70B+NL71B vs the other four, IL66B alone vs primaries, host marker panel; R5 human counts beside rat counts for the 38. Step B: fgsea on every fit's Wald stat (10,000 perms, 15-500, hallmark/KEGG/GO BP/brain separately), relabelling null, fgsea vs the Broad run | ~20:15 | `/rs/rs_grp_oschome/go2432/u251_host/verify/results/` (`summary.json`, `r0_*`..`r5_*`, `r4_gsea_*`), log `verify/logs/verify_40475691.out` |
| 40475693 | `run_de_therapy.sbatch`: nf-core/differentialabundance 1.5.0, host_therapy only, GSEA on; six-tumour sheet + the original 17,451-gene universe so every input equals the crashed run's | ~19:45 | `/rs/rs_grp_oschome/go2432/u251_host/results_de_therapy/`, `verify/de_therapy_compare.txt` (byte-compare with results_de and the unpublished GSEA) |

- Scripts: local copies in `ANALYSIS/host/verify/` (uncommitted), run from `/rs/.../u251_host/verify/`; grid clone untouched.
- Found while building: the two-contrast run's Broad GSEA for host_therapy DID finish (exit 0, never published because
  PLOT_EXPLORATORY crashed): `work_de/9a/6e148a...`, copied to `verify/broad_gsea_existing/`. 4,843 sets, **none at
  FDR < 0.25 either way** (min FDR 0.33); top nominal: immune/antigen-presentation up, ribosome/translation down.
- Pipeline lfc column is ashr-shrunk (lfcShrink type ashr, DESeq2 1.34.0); p-values are the unshrunk Wald test's. Model
  `~ 0 + Classification` with salmon lengths as avgTxLength; filter on all nine samples before the subset.
- Local dry run (DESeq2 1.50.2, no ashr; grid numbers will differ slightly): same 38 genes (Jaccard 1). Leave-one-out
  keeps 8-19 of the 38 when a primary or an NL tumour is dropped (37 without IL66B, which gives 398 genes). Relabelling
  null: the true split ranks 2nd of 10 (38 genes); IL69B+NL70B+NL71B vs the rest gives 2,752, along a myeloid-infiltrate
  (Ptprc, Cd68, Aif1 high) vs white-matter/neuron axis. NL70B+NL71B vs the other four: 502. The 23 down genes move in
  all three recurrences (IL66B ~100 % of the NL shift); the 15 up genes are NL-specific (IL66B median ~29 %). fgsea
  gene sets: every relabelling gives as many or more sets at padj < 0.05 as the true labels.
- When they finish: read `verify/results/summary.json` and `verify/de_therapy_compare.txt`; B7's footer and notes then
  need the gene-set line and the relabelling-null line.

## 2026-09-26 evening: host deconvolution (DONE 18:33, nothing pending)

- Jobs: 40473363 tasks 1 (LM22, requeued after preemption) and 2 (Zhang 2014) COMPLETED; task 3 (Bowman 2016) FAILED.
  Rerun 40475433 (prep) + array 40475461_1-4 all COMPLETED 18:33. Grid work dir
  `/rs/rs_grp_oschome/go2432/u251_host/deconv_fix/` (scripts/, results/, logs/); grid clone untouched.
- Bowman failure: permutation draw 169 (seed 42) had no positive weight at any nu, so CoreAlg divided 0 by 0,
  which.min() of three NaN RMSEs was empty and out[[mn]] stopped the run. About 1 draw in 1,000-1,500 with 4 columns.
  The signature itself is sound (no NA, no zero-variance rows, kappa 11.6); the deeper problem is that four myeloid
  columns cannot describe whole brain: real samples fit at r 0.01-0.08 (P 0.21-0.42), controls below 0.
- Fix (CIBERSORT.R untouched): `ANALYSIS/cibersort/run_v104_guarded.R` + `.sbatch` draws the same null outside
  CIBERSORT() and records degenerate draws instead of stopping. Reproduces 40473363_1/_2 exactly (fractions, r, RMSE
  and P, max diff 0). Plus `ANALYSIS/host/build_combined_signature.py`: Zhang non-myeloid (6) + Bowman myeloid (4),
  1,510 genes, kappa 102.
- Results: `ANALYSIS/cibersort/results/v104/fractions_v104_host_*.tsv` (4 signatures), 3 v 3 tables in
  `ANALYSIS/host/host_deconv_summary.tsv` (`summarise_host_deconv.py`). Primaries v recurrences: nothing reaches
  BH q < 0.1 in any signature. Combined signature: TAM_BMDM 0.16 v 0.33 (Welch p 0.32), driven by NL70B/NL71B (0.43,
  0.50; IL66B 0.07). Not interpretable: Bowman alone (r <= 0.08), and LM22, which fits the controls at r < 0 (P > 0.9)
  and puts 7-20 % T cells into tumours in athymic nude rats.

## 2026-09-26 ~15:45 EDT: drug rematch DONE (40468316), deck slide 9 + B2 updated

`ANALYSIS/drug_rematch/REMATCH.md`: clinical compounds 54 -> 62, both barrier models 13 -> 22 (salt records had inflated TPSA:
escitalopram, paroxetine, metoprolol, trimipramine, propantheline now pass; lysergide, mefloquine, tolonium, vandetanib enter).
Ciclopirox weighted rank unchanged in all three arms (1 / 1 / 3); unweighted 2 / 2 / 16 behind deferoxamine (ADMET-AI 0.53,
BOILED-Egg out). Top 5 weighted unchanged. S12 and the prior-art table (S15) of the manuscript need the same update.
Host rnaseq 40468361 aligning at 15:35 (N269B trimmed this time); DE / LOO / deconv chained afterok. GSC 40468077 still running.

## 2026-09-26 ~15:00 EDT (grid clock): PENDING GRID JOBS: drug rematch, host chain (Greg: "Yes" to the rerun; "add the host stuff")

| job | what | ETA | output |
|---|---|---|---|
| 40468316 | drug chain with exact-name ChEMBL matching, published + IL68B + IL66B, then compare | ~15:20 | `ANALYSIS/drug_rematch/REMATCH.md`, `rematch.json`, `runs/<arm>/subtypes/drug_ranking_final.csv` |
| 40468361 | host nf-core/rnaseq, resumed (N269B TrimGalore task deleted: it had exited 0 with empty output) | ~16:30 | `/rs/.../u251_host/results_rnaseq/star_salmon/` |
| 40468362 | host DE (afterok) | +40 min | `/rs/.../u251_host/results_de/` |
| 40468363 | host leave-one-out separation (afterok) | +10 min | `ANALYSIS/holdout_separation/loo_separation_host.tsv` |
| 40468364 | host deconvolution prep, then queues CIBERSORT v1.04 per signature (LM22, Zhang 2014, Bowman 2016) | +1 h | `ANALYSIS/host/results/`, `ANALYSIS/cibersort/results/v104/fractions_v104_host_*.tsv` |
| 40468077 | stem-cell drug screen on published GSC sets (running since ~13:10) | ? | `ANALYSIS/gsc_drugs/results_published/` |

- Host failures so far: 40468032 FastQC N269B could not write node-local /tmp (fix: host.config binds a CephFS dir as /tmp);
  40468266 FQ_SUBSAMPLE N269B got 0 records because the cached TRIMGALORE N269B task had empty outputs with exit 0.
- Grid clone: untracked copies that collided with the pull were moved to ~/u251_pull_backup_20260926 (only diff: a comment).
- Published arm of the rematch done: weighted top = ciclopirox, pentetrazol, nilutamide, primidone, pyrantel;
  unweighted top = deferoxamine, ciclopirox, metformin. When the job ends: update slide 9 table + text, B2, notes.
- ChEMBL audit (wf_f591d5b1-c19): 71 names, 2 tie-breaks, 33 overrides; `ANALYSIS/drug_rematch/audit_decisions.json`.

## 2026-09-26 afternoon: talk v4 built (SLIDES/CNS2026_Schwing_Abstract418_CNStemplate_v4.pptx)

- Greg: "synthesize all of these results into a compelling story and update the presentation", then "make sure these
  plots are styled like typical plots in this field ... not just inventing random ways to portray data".
- 13 slides + 4 backups, `SLIDES/09_build_deck_v4.py` (plan: `SLIDES/STORY_v4.md`). Figures in the field's standard
  forms from `SLIDES/10_standard_figures.py`: Broad GSEA enrichment plot (running ES recomputed from the ranked list,
  matches the pipeline's column to 5e-8), leading-edge heatmap (rlog z-scores; IL68B highest on all 77 genes), leave-
  one-out FDR table, score vs human share with OLS fit and 95 % band, GSVA subtype heatmap (vst counts; AC P 0.002,
  q 0.021, adjusted 0.024), CIBERSORT stacked fractions, DepMap 24Q4 violins (U-251 more DOHH-dependent than 78 % of
  1,178 lines). Drug ranking and whole-pipeline hold-outs are plain tables. v3 (`04_build_cns_template_deck.py`) kept.
- DepMap 24Q4 `CRISPRGeneEffect.csv` + `Model.csv` downloaded to `depmap/` (figshare 27993248; gitignored, 429 MB).
- Still to fold in when Greg asks: host rnaseq/DE/LOO (40468032/33/34), host CIBERSORT, GSC drug screens
  (40467921, 40468077), human-cohort plan. Open questions for Greg: what the IL/NL prefixes denote; MRI tumour volumes
  or tissue weights; the printed abstract's numbers (NES -3.17, ciclopirox 3.40, "mitochondrial", DMOG/LY-294002)
  disagree with the pipeline.

## 2026-09-26 13:10 EDT (grid clock): status

- Host rnaseq resubmitted as 40468032 (the first, 40467893, died on a Nextflow resume lock: it was launched from the
  repository root while the IL66B hold-out's Nextflow held the root's session; host runs now have their own launch
  directories). Host DE 40468033 and host leave-one-out 40468034 chained (afterok). ETA ~16:30-17:30.
- IL66B hold-out DONE (`ANALYSIS/holdout_IL66B/COMPARISON.md`): 140 DE genes vs 35, but translation initiation q 0.31
  and no set at q < 0.25; ciclopirox rank 3 (NES -1.64, FDR 0.04); Neftel AC -0.50 (p 0.066).
- Graft-fraction relation DONE (`ANALYSIS/graft_relation/`, job 40468073): graft % partly collinear with group
  (r -0.70); translation, AC, MES1 per-sample scores rise with graft % (r 0.76-0.84); adding graft % to DESeq2 keeps the
  translation fall (-0.40 -> -0.57) but DE genes 102 -> 42; group and graft fold changes correlate -0.66. Adversarial
  review running (workflow wf_47e66eea-c0c) -> SUMMARY.md.
- Stem-cell drug screen 40467921 still running; 40468077 (published GSC sets, out results_published) queued after it.
- CIBERSORT v1.04 all 9 signatures DONE; LM22 on human reads is a clean negative (fit R < 0, p > 0.9).
- Background workflows resumed: host refs + GSC sets (wf_5b325374-7af), human cohorts + power (wf_49b86515-404).

## 2026-09-26 12:25 EDT (grid clock): PENDING GRID JOBS: host reads, IL66B hold-out, stem-cell drug screen

| job | what | ETA | output |
|---|---|---|---|
| 40467892 (array 2-10) | host-only reads: first N_host records of each `*_rat` file (checks: counts + first `both` read ID) | ~12:30 | `/rs/rs_grp_oschome/go2432/u251_host/fastq/` |
| 40467893 | nf-core/rnaseq 3.22.2 star_salmon on rat (Ensembl 110), 9 in-vivo samples | ~16:00-17:00 | `/rs/.../u251_host/results_rnaseq/` |
| 40467929 | nf-core/differentialabundance on rat: Recurrent vs Primary; Tumour vs Control | +1 h after rnaseq | `/rs/.../u251_host/results_de/` |
| 40467930 | leave-one-tumour-out separation on host counts | +10 min after rnaseq | `ANALYSIS/holdout_separation/loo_separation_host.tsv` |
| 40467923 / 24 / 25 | IL66B held out: DE -> R figure -> drug chain + compare.py (control arm copied from holdout_IL68B) | ~14:30 | `ANALYSIS/holdout_IL66B/COMPARISON.md` |
| 40467921 | stem-cell (Varn) targeted drug screen, S1 ORA + S2 recurrence GSEA, ChEMBL + ADMET-AI | ~13:30 | `ANALYSIS/gsc_drugs/results/` |
| 40467819_5 | CIBERSORT v1.04, Varn 12-class | ~13:00-14:30 | `ANALYSIS/cibersort/results/v104/` |

Not yet submitted: host CIBERSORT (host_deconv_prep.py -> run_v104.sbatch with JOBS=results/v104_host_jobs.txt) waits
for the rat TPM and for the reference-fetch workflow (Zhang 2014, Bowman 2016, seq-ImmuCC, rat/mouse->human 1:1).
Human-read leave-one-out (done): IL66B out separates best (349 genes at padj < 0.05 vs 102; silhouette 0.43 vs 0.22);
IL68B out loses the PC1 split. Post hoc.

## 2026-09-26 08:20 EDT (grid clock): IL68B held out, full pipeline DONE (`ANALYSIS/holdout_IL68B/COMPARISON.md`)

Jobs 40467314 (DE), 40467315 (R figure), 40467509 (drug chain; 40467316 failed building the Python environment:
scikit-learn 1.9 has no wheel for the grid's glibc 2.17, so the compiled stack now comes from conda-forge).
- Control arm reproduces the published run: DSigDB top 100 identical (Jaccard 1.00), ciclopirox NES -2.261 FDR 1.7e-5
  rank 1, GSVA and STRING identical. The Python stage drifts slightly (rebuilt ADMET-AI environment): top-20
  Jaccard 0.82, 13 against 15 compounds clear both barrier models. September's ADMET-AI version was never recorded.
- Held out: DE 43 genes (33 up) against 35, 21 shared (Jaccard 0.37), log2FC Spearman 0.837. Broad GSEA: the six
  translation sets all hold (NES -1.90 to -1.96), and 12 sets clear q < 0.05 against 1 published.
- Drugs: ciclopirox stays rank 1 (NES -2.146, FDR 6.1e-4; score 3.02 against 3.14); DSigDB top 100 overlap 0.41 with
  the control; final top 20 overlap 0.43. Pentetrazol, primidone, nilutamide, pyrantel, ifosfamide, magnesium,
  diazepam, progesterone stay in the top 20; ozone, d-penicillamine, l-citrulline, paricalcitol, hypochlorous acid
  and others enter (not audited against the prior-art table).
- Subtypes: Neftel AC falls in recurrence with or without IL68B (-0.53, p 0.007; held out -0.57, p 0.044).

## 2026-09-26 06:00 EDT (grid clock): PENDING GRID JOBS: the full pipeline with IL68B held out

Greg: "rerun the full drug discovery pipeline etc with that sample held out then and compare the differences" (taken as
IL68B). `ANALYSIS/holdout_IL68B/` (README there): nf-core DE → R figure script (DSigDB drug GSEA, ChEMBL, STRING) →
September Python chain (ADMET-AI, BOILED-Egg, final ranking; GSVA) → `compare.py`. A **control** arm re-runs the published
six-tumour inputs through the R and Python stages today, to separate container / web-service drift from the held-out effect.

| job | what | submitted | ETA | output |
|---|---|---|---|---|
| 40467314 | nf-core differentialabundance 1.5.0, five tumours (qos secondary head; tasks on slurm -q primary) | 06:00 | 30–90 min (containers re-pulled to the CephFS cache) | `results_therapy_v3_noIL68B/` |
| 40467315 | R figure script, control + holdout arms | afterany | ~20–60 min | `{control,holdout}/publication_figure/` |
| 40467316 | Python env build (CephFS), drug chain + GSVA both arms, compare.py | afterany | ~20–40 min | `COMPARISON.md`, `comparison.json` |

Status 06:53 EDT (grid clock): 40467314 on its last task (the HTML report); every other nf-core task COMPLETED exit 0.
40467315 and 40467316 PENDING on it; expected finish of the chain about 08:00–09:00 EDT.
Interim, read from the finished nf-core tables: DE 128 genes at padj < 0.05 (91 up, 37 down) of 18,948, against 107
(74 up, 33 down) of 19,351 published; at the pipeline's filter (padj < 0.05, fold >= 1.5) 80 against 82, 48 shared.
Broad GSEA (seed 1234): translation initiation NES -1.93, q 0.027 (published -1.99, 0.022); elongation q 0.009,
ribosome 0.010, selenoamino acid 0.021, GCN2 0.034, starvation 0.065, so five of the six leading sets clear q < 0.05
without IL68B. The leave-one-out run with IL68B dropped (six-tumour normalisation kept) gave -1.93 / 0.022: the two
agree to within the size-factor change.

**On "check again":** read `COMPARISON.md`; first check the control arm reproduces the published drug profiles and ranking.

## 2026-09-26 06:00 EDT: leave-one-tumour-out GSEA DONE (jobs 40464722/23/24, all COMPLETED)

Results committed in `ANALYSIS/gsea_leave_one_out/` (SUMMARY.md, loo_sets.tsv, loo_screen.tsv, reproduction.json).
The full six-tumour run at seed 1234 reproduces the published report exactly (8,869 sets, every NES, p and q).

- Robust in all 41 runs: translation initiation NES −1.88 to −2.02, nominal p < 0.001; all six leading sets p ≤ 0.012.
- NOT robust: its FDR q. All six tumours, five seeds: 0.022, 0.026, 0.067, 0.139, 0.312 (q < 0.05 at 2 of 5).
  Without IL68B: 0.001–0.032 (q < 0.05 at 5 of 5; 12 sets at q < 0.05). Without NL70B: 0.001–0.015. Without NL71B:
  0.007–0.335. Without IL67B, IL69B or IL66B: above 0.05 at every seed (IL66B: no set at q < 0.25 at any seed).
- So IL68B carries the SIZE of the gene-level fall (−0.40 → −0.10) but not the enrichment's rank; the fragile claim is
  "the one set that clears FDR (q = 0.022)", which the manuscript and the old deck make. The deck now says this
  (slides 11, 13, 15, 23, 24 and notes), reading every number from these files.

## 2026-09-26 05:30 EDT: leave-one-tumour-out GSEA running on the grid (superseded above)

Question: does the GSEA result survive when any one tumour is left out? Translation initiation (NES −1.99, FDR
q = 0.022) is the only set at q < 0.05; one primary, IL68B, is the highest of the six tumours on every
translation-initiation and ribosome gene, and without it the gene-level fall is −0.10 instead of −0.40
(`SLIDES/07_contamination_magnitude.py`).

Scripts: `ANALYSIS/gsea_leave_one_out/` (README there). The exact nf-core `GSEA_GSEA` command (GSEA 4.3.2, the same
Galaxy-depot image re-pulled to `/rs/rs_grp_oschome/go2432/singularity_cache`, the pipeline's own GCT / CLS / GMT / chip
copied from `work/85/373d2dd…`), one tumour removed and the gene filter (≥ 10 reads in ≥ 1 sample) re-applied;
41 runs = full + six drop-one conditions at seeds 1234 and 1–4, and six unfiltered drop-one runs at 1234. The grid clone
of this repository is behind (9b4a583, with local changes); the scripts were rsynced into it, not pulled.

| job | what | status at 05:30 |
|---|---|---|
| 40464722 | prep: image, inputs, 10-permutation smoke test, `make_inputs.py` | COMPLETED |
| 40464723_[0-40] | GSEA runs (reqp, requeue qos, 2 CPU, 16 GB) | 22 of 41 COMPLETED, the rest ~19 min in (each run ~18-20 min) |
| 40464724 | `collect.py` | PENDING on the array |

First submission (40464589/90/91) failed in 9 s: `~/.local/bin/singularity` has no starter on the compute nodes; the
scripts use `~/mambaforge/envs/nextflow/bin/singularity`, as `run_*.sh` do.

**On "check again":** read `ANALYSIS/gsea_leave_one_out/SUMMARY.md` on the grid; `reproduction.json` must say
`"exact": true` before anything else is read; then update the Abstract 418 deck (slide 15, conclusions, limitations).

## 2026-09-26 05:30 EDT: the Abstract 418 build is self-contained

`SLIDES/cns_template.py` vendors the slide plumbing that used to be imported from CTSpinoPelvic1K-1/paper/cns2026;
`SLIDES/u251_paths.py` points 05/06/07 at this repository's result folders (gitignored, the same paths on the grid;
`SLIDES/fetch_grid_outputs.sh` fills them in a local clone); the CNS template goes in `SLIDES/templates/`
(gitignored: the CNS's file, and this repository is public). The deck is written to
`SLIDES/CNS2026_Schwing_Abstract418_CNStemplate.pptx`. Rebuilt this way, all 47 figure and data files are
byte-identical and every slide's text, notes and images match the previous build.

## 2026-09-26 03:25 EDT: Abstract 418 deck v3 (25 slides)

Talk SSTU02, Sunrise Science Session Tumor 2, Monday 2 November 2026, 7:00–7:06 AM, Room 147B, Walter E. Washington
Convention Center. Commit 16aab51. Greg: "19 is not enough, use the MBR manuscript and the grid" → 25 slides
(culture-to-tumour PCA, ten-library read sorting, running-sum plots, growth / glycolysis / respiration / sterol chart,
contamination magnitude, drug scatter, prior-art table, DepMap). Two adversarial checks (448 claims on v2 → 104
confirmed; 376 on the rebuild → 34 confirmed), all applied. Load-bearing: the translation fall leans on IL68B (slide 15,
conclusions, limitations); DepMap DHPS and EIF5A are needed by most lines, only DOHH is uncommon; S12 names 92
compounds, not 93.

Manuscript errors found, NOT edited: nominal p < 0.001 for all six (selenoamino-acid metabolism is 0.002); OXPHOS
"fell" depends on the collection; the S7 and facts.json BH q column is mis-ordered (lowest correct q 0.075); S15 marks
thioridazine (rank 12, tier A) as unranked; the strongest reversals are dimethyloxalylglycine and deferoxamine, not
LY-294002; "graft reads only" while main.nf merges the shared reads; "has not been profiled" against Nagaraja 2026's
own RNA-seq (whose abstract says 4 per arm, our Methods 4 + 3); PERMANOVA ran on all genes, not the 500; the audit
covered the top 20, not "every ranked compound"; 93 compounds is 92.

## 2026-09-26 01:05 EDT: Abstract 418 deck v2 (19 slides)

Rebuilt on the CNS template from the 8 September co-author deck (commit 0f011e6). Greg: "a six-minute talk, more
results and better explanations"; the rat cartoon dropped for a drawn pipeline strip; every data slide says what the
figure shows, the number and what it means; three new results slides (what else moves, the drug funnel, the DepMap
check). Ablation video (64 px, 7.9 s, from `Ablation video_TNN_18Sept26.pptx`) re-encoded at 512 px on the model slide.
Disclosures in the title notes (no disclosure slide); Greg has not confirmed that choice. PowerPoint's COM instance
holds the file for ~10 s after rendering: build to a temporary name and copy with retries.
