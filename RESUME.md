# Working checkpoints: U251 LITT recurrence study

Newest first. Everything for this study lives in this directory and this repository (Greg, 2026-09-26: "this project
should be contained in the u251 dir and the respective git repo"). The checkpoints below were first written into
spinesurg-ct-nnunet/docs/RESUME_2026-09-22_PM.md by mistake and were moved here the same day.

## 2026-09-27 08:25 EDT: the two assays are joined sample by sample; the methylation deposit is staged

**Greg's directive (2026-09-27): do NOT rename anything. "Just make sure GEO and people know which samples
correspond to multiomics data of the same sample since there's different names."** A rename scoped through the
pipeline was started and abandoned on his instruction; nothing was renamed. What the abandoned pass established, and
is worth not rediscovering:

- `N1` and `N2` in `ANALYSIS/power/` and `ANALYSIS/spread/design/` are **null-hypothesis labels**, not samples
  ("Null N1 : contralateral barcode counts ~ Multinomial(D, p_core)"). A token sweep of `N2` would have corrupted
  them. `ANALYSIS/power/pilot/` uses `N2` in the sample sense, so it does not even split by directory.
- `PREREG*.md` and `TARGETS.md` must never be edited: their SHA-256 is recorded inside the result files
  (`prereg_sha256`) and quoted in the manuscript. The result JSONs carry `_sha256` stamps, so their keys cannot be
  hand-edited either - they can only be regenerated.
- `IL` is the **laboratory's own** animal prefix (its RNA key says IL-70/IL-71; "IL" is its MRI prefix). `NL` appears
  only in the sequencing core's FASTQ file names, and from there in GEO and SRA. The public record carries the odd
  label and is the one that cannot be changed.

**Two things verified against the outside record, both of which decided the design.**

1. Of the ten GEO Series carrying both a methylation-array and a sequencing assay type, **ten of ten are
   SuperSeries**. So a SuperSeries is the route; arrays cannot be added to GSE338105.
2. In GSE335256, the sequencing SubSeries carries one `!Sample_relation = BioSample` per sample and the array
   SubSeries carries **none**. GEO mints no BioSample for an array sample, so **there is no accession the two arms
   share** and NCBI will not link them for anyone. The link exists only in metadata the submitter writes.

**What was built (commit `61f07f5`, pushed to origin/research-letter-expansion).** The correspondence is stated four
times and inferred nowhere: in each array sample's title (`U251N recurrent DNA methylation [IL70B = RNA library
NL70B]`), in five characteristics fields per sample (matched RNA library, GEO sample, SRA experiment, BioSample,
array label), in the Series summary, and in `GEO/multiomics_sample_key.tsv`, deposited as a supplementary file.
`GEO/make_sample_key.py` generates the key and cross-checks the array core's sheet against the accession map in both
directions, refusing to write it if they disagree on any pairing, rat number or chip position. The three two-name
pieces: **NL70B = IL70B, NL71B = IL71B, N269B = N2**. IL64B and N168B have no array.

**Grid: job 40483632 (`u251_geo_meth`), submitted 08:24, ETA ~15-25 min** (unzip 16 IDATs, pull the pinned minfi
container, noob betas + detection P values for 8 EPIC arrays, checksums). Log:
`~/u251-xenograft-murine-RNASeq-study/GEO/methylation/u251_geo_meth_40483632.out`. Its predecessor **40483631 failed
at exit 127** because `singularity` is on the login node's PATH and not a compute node's; the job now resolves
`singularity` or `apptainer`, loading the module if needed, and refuses in step 0 instead of after staging.

**Greg's to do, outward-facing.** `GEO/methylation/COVER_NOTE.md` is a ready-to-send draft to geo@ncbi.nlm.nih.gov
bundling the SuperSeries request with four corrections the audit found in the already-public record: no pointer from
GSE338105 to the methylation data, an SRA study descriptor that still says the data are private until July 2027,
IL64B described as a "failed graft" that the laboratory's key does not support (and the three controls carrying one
identical string), and `salmon.merged.gene_lengths.tsv` missing from the deposit so the DE step cannot be run from it.
Item 4 should be cleared with the laboratory first; the others are facts about the record.

## 2026-09-27 07:45 EDT: MANUSCRIPT DRAFT 1 AND DECK v5 ARE WRITTEN. Greg's directive: the instrument finding is the headline

**Two deliverables exist on disk now.**

1. `MBR/nature_draft/manuscript.md` (untracked, 651 lines) - draft 1 of the Nature/Cell-form paper, written from the
   pre-registered analyses: PREREG.md (sha256 `1b823ad5...2383d08`), SPREAD/SECTION_spread.md, XVAL_RESULTS.md, the
   CURVE result, the RNA and DNA science blocks, TARGETS.md (sha256 `88b5ac10...`) and `power_closed_only.json`
   (job 40481525). Structure: abstract, introduction, five aims of results, the contralateral subsection inserted
   verbatim from SECTION_spread.md between its delimiters, discussion, "the experiment this pilot sizes" with the
   target table, condensed methods, the figure list, and the editorial list of sentences that must not enter it.
   24 `\PENDING{}` slots remain and each names the job or file that fills it. `CHANGES.md` beside it records every
   insertion and every slot closed.

2. `SLIDES/CNS2026_Schwing_Abstract418_CNStemplate_v5.pptx` - **46 slides**: v4's 21 unchanged, plus 5 main and 20
   backups. Built by `SLIDES/14_build_deck_v5.py`, which opens the v4 deck and inserts; figures by
   `SLIDES/13_v5_figures.py`, which imports the manuscript's own drawing code and re-saves each panel at slide type.
   Plan and cut order: `SLIDES/STORY_v5.md`. Every slide was rendered to PNG through PowerPoint and looked at; three
   defects found that way were fixed (text colliding with the footer on two slides, the dilution null quoted in the
   wrong direction on the mesenchymal slide, the nine-panel power figure too small as a top figure).

**Greg's directive, 2026-09-27 (verbatim in intent).** Even if the contralateral signal is carry-over on the
instrument, *that is the positive finding for a neurosurgical audience*: the cellularity of the tumour is detectable
molecularly on what the surgeon touched, and surgeons should be aware of seeding distant brain with their
instruments. Both deliverables now carry it that way:
- manuscript: a new results subsection, "What both readings carry: tumour on the instrument is measurable", plus a
  paragraph in the discussion. No rate or risk during a human resection is asserted.
- deck: main slide 14, "Tumor in the opposite hemisphere - and what it means in the OR".
- anchor, measured and checked against Europe PMC before it was written: malignant cells recovered from gloves or
  instruments in 14 of 47 canine cancer operations (30 %), 57 % with involved margins against 19 % with clear
  margins, P = 0.016 (Orjefelt et al., J Small Anim Pract 2026;67:528-533). Case reports that such cells can
  establish: Steinmetz 2001 (J Neurooncol 55:167-171), Perrin & Bernstein 1998 (36:243-246), Bouillot-Eimer 2005
  (Clin Neuropathol 24:247-251), Bekar 2010 (World Neurosurg 73:719-721). What is new here is that species identity
  MEASURES the transferred material at the distant site; cytology and case reports cannot.

**Results that closed slots since the 02:30 block.**
- Bench lists (40481545): the H-C2 criterion holds under **all fourteen** published cross-species exclusion lists -
  true split rank 1 without the fraction covariate, rank 4 with it, in every list. The methylation arm difference is
  tumour fraction, replicated.
- Figure 6 draws from the closed-form power file; panels e and f (T4, T5) stay blank until `u251_pwr_asm` lands.
- Kept-gene count for the decomposition: 10,192 of 12,377 one-to-one ortholog pairs.

**Still pending on the grid (nothing blocks the drafts):** E1 external floor (40481267 running, 40481268 queued);
Aim 3 (pub_count 40481516 running, pub_noiseq 40481517 queued - and by pre-registration Aim 3 goes to the lab before
a sentence about it is drafted); the BENCH -> f_DNA chain (40481546-49) that would replace every uncalibrated
fraction; `u251_pwr_asm` for T4/T5; the multi-omics runs (rnavar, rnasplice, de_null_B, invasion_niche, irp_graft;
rnafusion blocked on its reference).

**For Greg.**
- The v5 deck and its figures are committed to this PUBLIC repository, as v3 and v4 were. The manuscript draft stays
  untracked (`MBR/` is gitignored). Say if slide 14 should wait until the manuscript is submitted.
- The dissection record does not say whether a fresh blade was used between pieces. If Raj's records say it was, the
  "one pass of an instrument" sentence changes in both deliverables.
- Unverified references still marked pending: Woodworth 2013 (29 % of reoperated recurrences held no active tumour),
  Wang L 2022 and Drexler 2024 for the mesenchymal shift, the J Neurosurg 2026 page range, and Shorrocks 2013 for the
  Shapley attribution (kept in the speaker notes, not the footer, until checked).

## 2026-09-27 02:30 EDT (grid clock): pipeline_v2 VERIFIED (refutation pass); GSE310817 batch-1 defect wired into SELECT; jobs 40481813, 40481814

**Checked and holding:**
- PREREG.md sha256 `67e9e250...832` recomputed; it matches in the scratchpad, local and grid copies.
- The bench code follows PREREG 3-9 as written: f grid, pairings, cross-fitted anchors, C1a/C1b criteria, the Q/Q_null
  metric, and the 8.1/8.2 rules.
- Mixing uses raw SigDF M/U fields with S* = S-con, not betas. The exclusion check pairs rat 0 % and 100 % arrays at the
  same f and the same human end.
- GSM-to-species and GSM-to-methylation-level mappings were re-checked against GEO for all three series.
- INTERFACE file names agree between the align/collect jobs and REFCHECK/E5.
- The bespoke `u251_meth/pipeline/` is untouched; only its README changed, at 23:54, by the nf-core workflow.
  40479710 had been cancelled at 21:34, before this workflow started; its replacement 40480027 finished at 21:57.
- README citations resolve and support their sentences. "provenance" does not occur. Nothing in
  `ANALYSIS/methylation_v2` was ever committed, and no IDAT is tracked.

**Found (high): GSE310817 batch 1 has no red channel.** Every Red IDAT of GSM9309404-416 is a byte copy of its Grn IDAT;
the titration workflow found this at 01:02. pipeline_v2's 01:15 block had read it as a "weak red channel". It decides
F2's C1b: F2 is "not evaluable in batch 1" and therefore fails under PREREG 5.1. It also puts 10 of C1a's 28 mixtures
and the pooled F1 curve on those arrays.

PREREG was not edited. SELECT still applies it as frozen, with batch 1 included. What was added:
- `select.R` records the defect in `selection.json` (`input_defects`), plus two reported-only items, with a flag saying
  they were never used for selection:
  - C1 over batch 2 alone (`sensitivity_c1_batch2.R`);
  - a C2-exclusion tail diagnostic (`diag_c2e_tail.R`: Q over probes the human end detects, and an independent
    recomputation of Q).
- `report_v2.R` prints these items.
- README.md and a METHODS.md section 9 describe the defect, plus two E5 corrections: the command admits 2 mismatches,
  not 1, and Needhamsen's FASTA has R already set to A.

Everything is synced to the grid and all 30 R scripts parse.

**Jobs:**

| job | what | state and ETA |
|---|---|---|
| 40481813 | u251_v2_verify: IDAT identity (`ref/geo/IDAT_CHANNEL_IDENTITY.tsv`), the batch-2 C1 sensitivity and tail diagnostic (`pipeline_v2/results/sensitivity_not_preregistered/`), and a code test of C2E/SELECT/REPORT in tmp (deleted after) | PENDING 02:25; about 20-30 min once it starts |
| 40481814 | u251_v2_xls: which ID forms Needhamsen's Additional file 4 holds | DONE 02:26. The file holds 17,944 IDs: 17,916 cg and 28 ch. R2's regex had dropped the 28 ch; `method_e5.R` is fixed and synced. R2 now also reports a cg-only Jaccard and whether their CHR/Mapinfo match our Rnor_6.0 loci, which answers the "assembly not stated" unknown |
| 40481317 | dependency is now `afterany:40481260:40481813`, so SELECT sees the files | ETA unchanged |

Logs: `u251_meth/logs/v2_verify_40481813.out` and `v2_xls_40481814.out`.

**Check before recording the selection sha256 (the 00:36 gate):**
1. The log shows `RED_IS_BYTE_COPY_OF_GRN` for 13 GSE310817 arrays and no others.
2. `SELECT field ... PRESENT` three times.
3. `compare_q.py` prints AGREE.

**Greg's decision, not made here:** whether to adopt "C1 on batch 2 only" as a PREREG section 11 deviation. It would
be dated after results were seen. `sensitivity_not_preregistered/c1_batch2.json` states what it would change, and adds
this: every F2 candidate's C2 mean MAE (0.096-0.124) is above the 0.05 adequacy bar. So even if an F2 variant became
eligible, PREREG 8.1 would still keep the fraction out of downstream adjustment.

## 2026-09-27 02:12 EDT (grid clock): XVAL (matched DNA array x RNA-seq cross-validation) DONE — all four tests run, nothing pending that changes a number

Greg's question: can the matched arrays show that SNPs or transcripts are real rather than errors? The write-up is
**`MBR/nature_draft/XVAL_RESULTS.md`**. It covers:
- what was cross-validated, and each criterion's outcome;
- what each test validates and what it cannot;
- the methylation-expression concordance rates with their nulls;
- draft Methods and Results paragraphs;
- captions of 60 words or fewer.

Identical copies are at scratch `xval/RESULTS_xval.md` and grid `u251_science/xval/results/RESULTS_xval.md` (sha256
`8494f8e8...` after the 02:40 hostile check; it was `6765d0a0...`). This block supersedes the "Pending" line in the XVAL
sub-bullets of the 22:55 block.

**Hostile check, 02:15-02:40 EDT.** No number changed a decision. What was checked:
- the PREREG hash in all three copies;
- that every table number traces to `xval/`;
- eleven rs positions and alleles against dbSNP (all match);
- the bin liftover, re-done with pyliftover (25,730 of 25,730 identical);
- the decile-matched dosage null, and that the floor handling comes before the dosage statistic;
- no `sra_*`/`xs_*`/`star_*` jobs, and nothing committed.

What was corrected in the write-up:
- The Results draft had said the 0.25-0.40 % descriptive rate sets the RNA-only threshold. kmin uses the
  pre-registered e_s upper bound per depth bin.
- The five sites: array betas of 0.90-0.94 at four of them allow at most about 10 % DNA ALT, against 14-57 % in RNA.
  The pre-registered paralogue check (PREREG 2.3) was not run, so three readings stay open.
- Job states: 40480668 was FAILED in Slurm, and 40481699 has COMPLETED.
- Caption word counts.

`ANALYSIS/xval/` was added to `.git/info/exclude`. The recomputed tables are in `xval/results/hostile_check/`.

**PREREG.** `PREREG_xval.md`, sha256 `140fdd44c15143ddfc70d5852fe9aa7d69e7a9897b2047ebc19b2afaf8568c75`. It was
re-hashed at 02:06, locally and on the grid, and is unchanged.
- Deviations 1-30 are in scratch `xval/DEVIATIONS.md` and grid `xval/DEVIATIONS_xval.md`.
- Entries 26, 27, 29 and 30 were written after the result they concern was read, and each is labelled so.
- No deviation changes a pre-registered decision.

**Jobs.** All ids are in `/rs/rs_grp_oschome/go2432/u251_raw/JOBS.txt`, with a status line appended at 02:11.

| step | job | state |
|---|---|---|
| SRA fetch (`sra_*`) | none | not submitted: all ten human STAR BAMs exist; PREREG_xval 1.5 excludes it |
| xengsort re-sort (`xs_*`) | none | not submitted: no test needs it |
| STAR (`star_*`) | none | not submitted: the BAMs exist |
| bcftools 1.21 pile-up | 40480752, 40480754 FAILED (no `--regions-overlap`); 40480755 CANCELLED (too slow); **40481325** (40 chunks) | 40481325 COMPLETED; chunk 30 was last, at 02:03:39 |
| test 1 analysis | 40481651 (primary only) COMPLETED 01:48; **40481326** (primary + duplicates kept) COMPLETED 02:05:02, exit 0 | the primary files are byte-identical between the two runs |
| tests 2-4 | 40480794 xval_deseq, 40480795 xval_tests, 40481406 xval_desc, 40481374/40481408 xval_fig | COMPLETED |
| cosmetic | 40481699 xval_geno1 (redraw of the primary-only figure) | COMPLETED 02:20:47, exit 0; duplicates-excluded rows identical to `rna/` |

**Result paths** (grid, `G/u251_science/xval/`):
- `rna/`: test 1. It holds `setB_contradiction.tsv`, `error_model*.tsv`, `kmin_table.tsv`, `rna_vs_rna_agreement.tsv`,
  `genotype.json` and the merged `pileup_sites*.vcf.gz`. The primary-only copies are in `rna/nodup_only/`.
- `cnv/`: test 2.
- `meth/`: test 3.
- `lists/`: test 4.
- `results/`: `fig_xval_{genotype,dosage,meth}.png` and `RESULTS_xval.md`.

**Findings.**
- *Test 1, SNPs.*
  - Set A (the 59 rs probes) is **not testable**: the probes are never expressed, with at most two reads.
  - Set B (CpG-SNP sites where the array implies hom-REF): 2.95 % contradictions (16 of 543; Wilson 1.8-4.7 %). That is
    between the 2 % tolerance and the 5 % failure line, so it is reported. The 16 come from five common-SNP sites
    that recur across libraries (FKBP9, PGM1, EP400, ELL2, RBM12B). What they are is not known: a low-fraction DNA
    allele, an RNA-DNA difference or paralogue reads.
  - B + B\* gives 1.39 %. With duplicates kept: 2.31 % (B) and 1.08 % (B + B\*).
  - RNA-vs-RNA agreement is 99.86 % (pass).
  - The pre-registered e_s is 1.5-2.7 %, above the 1 % bar, because the five sites are counted as errors. With them
    removed (descriptive), it is 0.17-0.40 %.
  - The kmin table, built from the pre-registered e_s upper bound per depth bin, is the null for any RNA-only
    variant. This validates the pipeline, **not somatic subclones**: every
    sample shares the U251N germline, and the array does not sequence.
- *Test 2, copy number.*
  - The positive controls pass: CDKN2A/B is silent, and the dosage slope is 0.54 (p 7e-6).
  - The pre-registered rule passes N2's 3p (p 0.0083) and, borderline, 4p (0.0120; a second draw gave 0.0134).
  - But N269B's arm-level expression does not follow N2's copy number anywhere: slope 0.17 (p 0.45), against 0.97 for
    IL69B.
  - So **the N2-only losses are neither confirmed nor shown to be noise**. They stay array-only.
- *Test 3, methylation.*
  - Within sample, rho is -0.33 to -0.36 in all seven U251-rich samples; no shuffle out of 10,000 reached it (pass).
  - Delta-delta: rho -0.003, eighth of ten relabellings, shuffle p 0.37, no gene selected. There is no coordinated
    change.
- *Test 4, the RNA lists.*
  - The RP block has no DNA support, and its unadjusted promoter hypomethylation was a tumour-fraction effect.
  - The fraction-adjusted up-list is "DNA-supported" only by the letter (+0.013 log2, three genes, p 0.043).
  - Absence of DNA support is not refutation.

**Shared inputs for the SPREAD workflow's C2.**
- **BAMs:** `/rs/rs_grp_oschome/go2432/u251_science/xval/bam/<lib>.bam(.bai)` are symlinks to the home-clone
  nf-core/rnaseq 3.22.2 BAMs, `~/u251-xenograft-murine-RNASeq-study/ANALYSIS/results_human_final/star_salmon/<lib>.markdup.sorted.bam`.
  Read them only; home is at 96.8 % of quota.
- **Re-sorted reads do not exist.** `/rs/rs_grp_oschome/go2432/u251_raw/xengsort/` was never created, because no
  xengsort re-sort was run.
  - The original sorted reads are at `~/u251-xenograft-murine-RNASeq-study/ANALYSIS/sorted_fastqs/<lib>_{human,rat}_R{1,2}.fq.gz`:
    human = graft + both, rat = host + both; 85 GB, 12 May 2026. The index is in `ANALYSIS/xengsort_index_clean/`.
  - A re-sort at another threshold needs SRA GSE338105 (SRR39547363-72, about 200 GB, into `u251_raw/`). Check
    `JOBS.txt` and `squeue` for `sra_*` and `xs_*` first.
- **Variant null:** `xval/rna/kmin_table.tsv`.
  - PREREG_xval 2.4 requires SPREAD C2 to state which threshold it used.
  - C2 ran (about 01:00) before the table existed (01:48). Its verdict is negative, so a stricter threshold cannot
    reverse it, but the text must name the threshold.

**Loose ends for other workflows** (not edited by XVAL):
1. `MBR/nature_draft/SPREAD_RESULTS.md`, around line 189 and in its status-table row at about line 467, still says
   the XVAL dosage check "has not run" because pile-up 40480752 failed.
   - The dosage check never depended on the pile-up. It ran at 00:29 and found that N269B's expression does not track
     N2's copy number.
   - Suggested replacement wording is in `XVAL_RESULTS.md`, under "Other documents that are now stale".
2. SCIENCE T2 gate G2a (01:30 block, item 3) waits on `xval/rna/genotype_concordance_per_library.tsv`. The file now
   exists but **has no rows** (Set A has no RNA depth), so G2a cannot pass on it as written, and SCIENCE needs its own
   deviation. The per-library Set B rows are in `xval/rna/setB_contradiction_per_library.tsv`.
3. `fig_xval_genotype.png` panel b draws N269B's zero rate at the 1e-4 axis floor. The draft caption says so.

**Not done, and not possible with these data:** somatic-subclone confirmation. The next experiment is targeted DNA
sequencing of the leftover array DNA. The cause of N269B's copy-number-independent expression is also untested.

**Nothing is committed.** `ANALYSIS/xval/` is untracked code, `MBR/nature_draft/` is git-ignored, and this RESUME edit
is uncommitted.

## 2026-09-27 02:10 EDT (grid clock): pipeline_v2: C2-exclusion / SELECT / REPORT code paths tested; ETA of 40481317 corrected to ~06:30-07:00 (latest ~10:30)

Code tests 40481655 and 40481656 ran C2 exclusion, SELECT and REPORT into `u251_meth/tmp/v2test_*`, with E5 marked as
not built. All three completed with exit 0. The test outputs are deleted: they are not the selection. Two cosmetic
fixes were synced before 40481317 runs: REPORT number formatting, and the PREREG hash read in `select.R`. The parse
check still passes (28 R scripts).

**ETA.** 40481317 (E5, C2 exclusion, SELECT) follows the ref collect job 40481260, which the fetch agent estimates at
~06:30 EDT, latest ~10:00. That puts the selection at ~06:30-07:00, latest ~10:30. The apply job 40481327 stays HELD
until the gate in the 00:36 block, then takes about 1 h plus queue time. The SCIENCE follow-up chain (40481545-40481549)
waits on 40481317 and 40481327.

**Seen in the test.** This is benchmark data only, with no tumour read. Every exclusion list had Q - Q_null of about
0.27 at the 99th percentile, against the 0.02 criterion, so none met it:
- U: Q 0.483, Q_null 0.176;
- none: Q 0.554.

If 40481317 gives the same with the E5 lists, the PREREG 8.2 outcome is "no published list brings rat contamination
to within 0.02 of the replicate floor". Downstream then uses the pre-declared fallback U, labelled as not protected.
The fraction outcome in the test was also "none eligible", as the 01:15 block predicted.

**Before any commit of `ANALYSIS/methylation_v2/`:** PREREG.md and METHODS.md quote two bespoke numbers derived from
our unpublished arrays (N2 0.264; 56 of 939 iDMCs). Nothing is committed.

## 2026-09-27 02:05 EDT (grid clock): SCIENCE DNA side + RNA-DNA integration: all registered DNA analyses run; BENCH follow-up chain queued (40481545-40481549)

SCIENCE PREREG (`scratchpad/science/PREREG.md`, sha256 `1b823ad5...`; grid `u251_science/PREREG.md`): every
DNA-methylation analysis (M2, H-M1 DNA, M3, L1, L4, L5, N0, H-C1, H-C2, H-C3, T3, T4, the S1 bound, the D3-DNA spike-in
floor) and the integration pieces (M3, Holm across RNA and DNA, the homogeneity triad, the DNA figure panels). Code
`ANALYSIS/science/dna/` (not committed; grid copy `u251_science/code/dna/`, `CODE_SHA256_<job>.txt` per job in
`u251_science/dna/`). Outputs `u251_science/dna/{j9,fdna,dilution,nfcore_cofactor/bespoke,integrate}/`; figures
`u251_science/figures/dna/` (fig1b, fig1c, fig4a-c, fig5c-d, PNG + PDF, `captions_dna.md`, all <= 60 words; every panel
looked at); logs `u251_science/logs/sci_*_<id>.out`; ids `u251_science/dna/JOBS_DNA.txt`. Frozen inputs
`u251_science/inputs/dna/` (bespoke work-dir objects copied with sha256), references `u251_science/ref/dna/` (LUMP 44,
Gabbutt 2022 EPIC fCpGs 1,794, Gabbutt 2025 lymphoid fCpGs 978, CCLE U251MG segments; `REF_SHA256_dna.txt`).
Deviations D-DNA-1 to D-DNA-10 appended to `scratchpad/science/DEVIATIONS.md` and `u251_science/DEVIATIONS.md` (incl.
a disclosure: the nf-core true-split FDR count, 46, was seen at 00:31 while locating files). Nothing committed or pushed.

| job | what | state / ETA (EDT) | output |
|---|---|---|---|
| 40481410 sci_dna_refs | freeze inputs, lists, CCLE; diptest; EPIC.5.SigDF.normal | DONE 00:49 | `inputs/dna/`, `ref/dna/` |
| 40481447 sci_dna1 (J9) | L1, L4, L5, D4, T3 gates + statistic, T4 + S1 | DONE 01:09 | `dna/j9/` |
| 40481488 sci_fdna | f_DNA resolver, M2 + H-M1, M3, D3-DNA spike-in | DONE 01:03 | `dna/fdna/` |
| 40481489 sci_dilution (J10) | N0: 680 mixtures; H-C1, H-C3, T3, T4 nulls (1,000 draws) | DONE 01:42 | `dna/dilution/` |
| 40481490 sci_nfcofactor (J11a) | nf-core/methylarray, cofactors = f_DNA (bespoke), ten splits | DONE 01:51 | `dna/nfcore_cofactor/bespoke/` |
| 40481491, 40481719 sci_dna_integrate | H-C2, triad, sentences, figures (40481719 = label-only rerun, D-DNA-10) | DONE 02:00 | `dna/integrate/`, `figures/dna/` |
| 40481723 sci_benchlists (interim) | H-C2 under the lists BENCH has already published (E1-E4, from its run 40481316); T3 exclusion without E5 | PENDING (priority) at 02:10; about 30-45 min once it starts | `dna/benchlists_interim_E1-E4/` (never read by 06_integrate; superseded by 40481545) |
| 40481545 sci_benchlists (J11b) | H-C2 under every BENCH list; T3 primary exclusion | afterany BENCH 40481317 (u251_v2; after ref_collect 40481260, ~10:00) | `dna/benchlists/` |
| 40481546 sci_fdna (rerun) | re-resolve f_DNA with BENCH | afterany BENCH C3 40481327 (**HELD** until BENCH's selection sha256 is recorded) | `dna/fdna/` (old kept as `fdna_before_<id>/`) |
| 40481547 sci_dilution (rerun) | T3 primary null; all nulls on the new f_DNA | afterany 40481545:40481546 | `dna/dilution_benchprimary/` |
| 40481548 sci_nfcofactor (rerun) | cofactor run on the new registered f_DNA | afterany 40481546 | `dna/nfcore_cofactor/<variant>/` |
| 40481549 sci_dna_integrate (rerun) | everything re-read, T2 picked up | afterany 40481547:40481548 | `dna/integrate/` |

**Registered f_DNA = the bespoke estimate, uncalibrated: every f_DNA-based number is `\PENDING{}`, every N0 percentile
is an "unvalidated null".** CURVE step 1 FAILED (MAE 0.242, max |error| 0.531); CURVE Amendment 1's calibrated reading
failed its gate too (and, post-freeze, was never eligible; D-DNA-2). BENCH has published nothing. The N0 machinery
itself checks out: its k_c equals CURVE's to 1e-15 for all ten rat arrays; the N0 code path reproduces the bespoke
amplitudes (|diff| <= 3.3e-4), fCpG W and MGMT promoter means (<= 1e-8); the exact 90 % interval code matches a
brute-force inversion of the 20-labelling test.

**Results (descriptive; 3 v 3 floor p = 0.05).**
- **L1 LUMP: FAIL.** 39 of 44 CpGs usable; LUMP 0.84-0.87 in the six tumours, 0.88 N2, 0.92 C2B; Spearman with RNA
  share -0.77. By section 8 every DNA-side statement carries this failure and "rat-blind" is withdrawn from the text.
- **L4 CCLE: PASS** (directional agreement 0.89 over the 19 arms CCLE calls). **L5 SeSAMe vs conumee: FAIL** on its
  letter (arm calls agree in 0.926 of 312 cells; amplitude slopes not rank-identical, Spearman 0.75). Both kept.
- **H-M1 DNA (bespoke, PENDING):** recurrent minus primary -0.196 (exact 90 % -0.340 to -0.059; Welch 95 % -0.367 to
  -0.025), rank 1 of 20 (p = 0.05); CURVE f_hat and f_cal also rank 1. RNA (replication) rank 2 (p = 0.10). Holm over
  RNA and DNA: 0.10, 0.10.
- **M3:** Spearman(RNA share, DNA fraction) 0.886 over six tumours, 0.929 over seven arrays. N2's DNA/RNA ratio 4.0
  (bespoke), 6.5 (CURVE f_hat), 1.0 (CURVE f_cal): unexplained, CURVE having failed.
- **H-C1 (unvalidated null, uncalibrated f_DNA):** C2B diluted in silico with rat DNA keeps its amplitude (A_dil 0.97-0.99
  at the tumours' fractions). The primaries sit below that curve (A - A_dil -0.064 / -0.094 / -0.103), the recurrences on
  or above it (+0.020 / +0.003 / +0.063). T_C1 0.115 (exact 90 % 0.067 to 0.166), rank 1 of 20, above all 1,000
  pseudo-recurrence draws (null 95 % 0.004 to 0.009); same under both CURVE fractions. Section 8 wording: "unexplained
  by dilution"; no homogeneity claim (the triad fails, below) and, with the null unvalidated, no claim at all yet.
- **H-C2:** replicated in the independent pipeline: nf-core/methylarray true split rank 1 of 10 at p < 1e-3 (2,800
  probes) without the covariate, rank 4 (341) with f_DNA; bespoke A rank 1 (3,975) -> B rank 4 (326). Criterion met in
  both; the covariate is the uncalibrated bespoke fraction; BENCH-list reruns pending.
- **H-C3:** STP27 "M" in all eight; flip distance 0.45-0.48 beta per recurrence against observed shifts -0.032 (Welch
  -0.077 to 0.013) and +0.057 (-0.006 to 0.120); promoter-mean difference -0.014 lies below all 1,000 null draws
  (null 95 % -0.002 to 0.003; unvalidated). MGMT TPM 0 in C2B and all tumours; 0-7 STAR exon reads per library.
- **T3 fCpG: UNINFORMATIVE** (G3b fails: C2B already multimodal at the fCpGs, dip p < 0.005 for every list; G3a fails
  for the lymphoid list and for noob betas). T_T3 about 0, rank 14 of 20. **Triad: not distinguishable from tumour
  share (T3 failed; T2 pending).**
- **T4:** T_T4 0.081, rank 5 of 20 (p = 0.25; above the pseudo-recurrence null q95 0.036, but the relabelling criterion
  fails): no convergence. **S1:** a shared single-copy loss is flagged on >= 80 % of arms at s >= 0.55 (from 2 copies) /
  0.80 (from 3); a shared gain never reaches 80 % (maximum 77 %: the unchanged rule cannot flag arms the primaries already
  call gained).
- **D3-DNA spike-in (POWER T5):** FDR < 0.05 recovery of a shared 0.05 / 0.10 / 0.15 / 0.20 beta shift: design A 19 / 46
  / 68 / 86 %; with the fraction covariate 0 / 0 / 14 / 22 %.

**Next.**
1. BENCH: the chain above runs by itself once BENCH's 40481317 finishes and its C3 job 40481327 is released (see the
   00:36 block for the manual gate). If 40481327 is cancelled instead, 40481546 re-resolves to bespoke and the chain
   still completes. Then read `dna/fdna/fdna_resolved.json`: if `registered` became `bench`, H-M1 (DNA), H-C1, H-C3 and
   the cofactor run stop being `\PENDING{}`; the N0 null stays unvalidated whatever BENCH picks (section 6 ties it to
   CURVE step 1).
2. H-T2 is RNA-side (sci_ai 40481604 -> `rna/allelic/t2.json`); `06_integrate.R` reads it on its next run (40481549, or
   `sbatch -D /rs/rs_grp_oschome/go2432/u251_science /rs/rs_grp_oschome/go2432/u251_science/code/dna/06_integrate.sbatch`).
   It cannot rescue the triad (T3 already failed its gate).
3. Not run, by rule: G3c (no verified accession of an independent U251 methylation array; NCBI eutils returned 500
   at 01:03); X1 (optional). Torsvik 2014 aCGH dropped (figures only; D-DNA-5).

## 2026-09-27 01:15 EDT (grid clock): pipeline_v2 bench run 1 DONE (40481316, 26 min, all tasks exit 0); 40481317 waits on the E5 alignment; 40481327 (apply) HELD

This block updates the pipeline_v2 block of 00:36 further down. Job ids are unchanged:
- **40481317** runs E5, C2 exclusion and SELECT. It depends on afterany:40481260, the ref collect job, which is still
  waiting on the bowtie-1 index builds 40481256/7 (running 1 h 05 min at 01:10). ETA about 03:00-05:30 EDT, plus
  15-30 min.
- **40481327** (apply) stays **HELD** until the selection sha256 is recorded. The gate is described in the 00:36 block.
- Code test **40481567**: C2 exclusion and SELECT run into `u251_meth/tmp/v2test_*`, not into results. It was queued
  at 01:10. A first attempt failed only because the test harness did not bind `/.rs`; Nextflow's own run of the same
  code path (C2 fraction) read C2B fine.

Benchmark facts from run 1. They come from public arrays and C2B only; nothing from the tumours was read.
- **R3 not reproduced.** The rule nearest Ebata's 133,377 is minfi "all" at 173,314, 30 % off; the criterion is 20 %.
- **C1a failed for every scaling.** The floor is median |d'| 0.101 and median |delta beta| 0.012. S* = S-con, so
  every C2 number is conditional on additive hybridisation.
- **C1b failed** for F1-i, F1-lin and F1-curve (MAE 0.15-0.22). F2 was **not evaluable** in batch 1 of GSE310817, so
  under PREREG 5.1 it fails C1b. In batch 1 the red Extension controls read 145-230 against 10,000-15,000 green, and
  RGdistort is 6-7.5. The batch-1 human array reads the human-green switch probes at r_p about 0.3, so the QC leaves no
  probes. Batch 2 behaves: pure human 0.93, pure mouse 0.005.
- **Consequence under PREREG 8.1:** no fraction candidate is eligible. SELECT should report the fraction as unresolved,
  and no tumour fraction will adjust anything downstream.
- **C0** passes only F2-50-lin, F2-50-curve and F2-45-lin (worst |error| 0.017-0.038).
- **C2 fraction** mean MAE runs from 0.088 (F1-i) to 0.40 (F1-Y); none is at or below 0.05.
- **Also seen:**
  - GSM5319816, a rat standard, has a failed red Extension control (568);
  - F4 precondition passes (C2B chrY 8.0x female blood);
  - F5's reimplementation matches InfiniumPurify exactly (939 iDMCs on EPIC, 66 rat-functional);
  - R1 reproduced;
  - E1-any 485,739 and E1-all 114,787 probes, as PREREG predicted for detection over ten arrays.
- **Not decided here:** whether to run C1b on batch 2 only. That would be a deviation after seeing results, so it
  goes in PREREG section 11 as a dated deviation, and only if Greg asks for it.

## 2026-09-27 01:30 EDT (grid clock): SCIENCE, RNA side: every pre-registered RNA analysis coded and submitted; the core has finished; S1, T2 and the Aim 3 reproduction are queued or running

PREREG `scratchpad/science/PREREG.md`, sha256 `1b823ad50b2eb6299a56cb4b13ef2e99d7dfc46a62ed2bc2508b486ca2383d08`. The grid copy
`u251_science/PREREG.md` is read-only, and every job checks the hash before it prints a number. Code is in
`ANALYSIS/science/` (not committed), with a README per analysis. The grid copy is `u251_science/code/` with
`CODE_SHA256.txt`. Outputs go to `u251_science/rna/`, except Aim 3, which goes to `u251_science/pub/` as PREREG section 3
says. RNA-side deviations P1-P17 are appended, never overwritten, to `u251_science/DEVIATIONS.md` and
`scratchpad/science/DEVIATIONS.md`. Job ids are in `u251_science/JOBS_RNA.txt`. `rna/RUN_RECORD_RNA.json` lists the job
states and which outputs exist. G = `/rs/rs_grp_oschome/go2432`.

**Finished; each log was read after its job completed:**

- **40481434 `sci_rna_core`**
  - **H-A1 (the primary endpoint) FAILS.** C-share is 0.466, and the threshold was > 0.50. Mixture is not the main
    source of the species-blind recurrent-versus-primary fold change.
  - Of the ten splits, the true split has the highest C-share. Spearman of C-share with |difference in f_RNA| = 0.85.
  - Sensitivities, all on the same side of 0.50: S2 0.449, S6 0.453, S3 (without IL66B) 0.369. The fraction-matched
    pair gives 0.031.
  - The covariance residual is flagged: 39 % of genes exceed 10 % of |LFC|.
  - H-M1 (RNA) ranks 2 of 20. This is a replication.
  - T1: the DiG contrast ranks 1 of 20, both fraction-adjusted (p = 0.05) and unadjusted. The PREREG expected it to fail.
  - H-E0 gate **FAILED**: IL64B (0.0024) is not below N269B (0.0019), so **H-E1 is abandoned**. The p-value the job
    printed is not a result.
- **40481453 (array) + 40481454: species-blind DE**
  - The true split has 133 DE genes; the ten splits range from 0 to 303.
  - **H-A2 fails**: 12.8 % of the true split's DE genes are composition-dominant.
  - H-A3 passes. It is descriptive.
  - C4: cell motility and inflammatory response are Holm-significant among the species-blind DE genes, and none of the
    three categories is in the human stream.
  - The composition-only view is empty by construction: phi_f cannot exceed about 0.62 here, and the view needs >= 1.
- **40481466 `sci_rna_r`**
  - C6: corr(f, arm) is -0.70 and the variance inflation 1.98. The arm term calls 30 DE genes (133 without f) and
    recovers 2 of 48 human-stream DE genes.
  - C5 (H-A5): D ranks 5 of 20 for Neftel MES and 9 of 20 for Wang MES, so the result is "not shown". But the
    human-only MES drop is below all 1,000 dilution-null draws for Neftel MES, and below 99.3 % of them for Wang MES.
    The rat-only MES score tracks myeloid cells (Spearman 0.66 and 0.60).
  - L2 ESTIMATE fails: Spearman 0.12 over the ten libraries.
  - T1 robustness: GSVA and the singscore package both rank DiG 1 of 20.
  - The S4 step failed here and was rerun as `sci_s4` (P17).
- **40481616 `sci_s4`**: S4 (HCOP orthologs) C-share 0.471, the same side as the main run.
- **40481582 `sci_host` (M4)**
  - Check (i) passes for BRETIGEA, NeuroExpresso MGP and CIBERSORT.
  - Check (ii): all six cell types are interpretable.
  - Check (iii): LM22 is rejected (T-cell share over 5 %).
  - H-M2a ranks 7 of 20 and H-M2b 9 of 20 (Holm 0.7 each).
  - The Liddelow/Zamanian, Keren-Shaul DAM and Dorrier fibrotic-scar lists were dropped: none is retrievable verbatim.
- **J0 jobs:** 40481377, 40481465, 40481439 and 40481379 are done. All references, the lab files, the JNS supplement, the
  published Neftel modules, the R library and the Bioconductor 3.14 and 3.13 libraries are in place.

**Queued or running (ETAs are grid clock, 27 Sep):**

| job | what | ETA | read |
|---|---|---|---|
| 40481475 `sci_s1_rnaseq` | S1: nf-core/rnaseq 3.22.2 on graft-only reads, then C2 | 04:00-06:00 | `rna/composition/composition_S1.json` |
| 40481604 `sci_ai` | T2 (allelic imbalance) | ~30 min after it starts; ~02:30 | `rna/allelic/t2.json` |
| 40481514 `pub_index` | Aim 3: rebuild the lab's second-pass STAR 2.7.1a index | ~03:00 | `logs/pub_index_40481514.out` (identity check), `pub/INDEX_DONE` |
| 40481515[0-6] `pub_align` | Aim 3: the lab's STAR command, route (ii) | ~06:00 | `pub/bam/*.ALIGN_DONE` |
| 40481516 `pub_count` | Aim 3: species split, featureCounts, D0-s, D0-c | ~08:00 | `pub/d0s_strand.tsv`, `pub/d0c_counts.tsv` |
| 40481517 `pub_noiseq` | Aim 3: D0-d gate, D2 (H-B1), D3 (H-B2), D4; Bioc 3.13 check | ~11:00 | `pub/noiseq/pub_noiseq.json`, `pub/noiseq_bioc313/` |

**When the jobs finish:**

1. S1: read the C-share in `composition_S1.json`. H-A1 has to be judged in sign under S1 and S4 together.
2. Aim 3: first check the index identity in the `pub_index` log. Then read `pub/d0s_strand.tsv` **first**. The D0-s
   outcome and the IL64B grouping go to Greg and to the lab (T. Nagaraja, N. Morosini) before any sentence about the
   published analysis is drafted. After that, read D0-c, then D0-d; on FAIL, only the discrepancy is reported.
3. T2: gate G2a stays pending until XVAL Test 1 writes `xval/rna/genotype_concordance_per_library.tsv`. XVAL's
   `xval_geno` 40481240 was cancelled. When the file exists, rerun only the `stats` stage; the command is the last call
   in `j8_allelic/sci_ai.sbatch`.
4. When SPREAD's E1 floor lands, rerun S2 and T1 on the E1-corrected export. This is not coded yet.
5. Rerun `code/common/run_record.py --sci G/u251_science`.

X1, the optional exploratory analysis, was not run.

## 2026-09-27 01:02 EDT (grid clock): human:rat titration curve DONE. In-silico mixing FAILED its test; order and N2 human DNA hold

Greg's question: mix pure rat and pure U251N in silico, draw a curve, and see where our eight arrays lie. **Nothing is
pending.**
- Full job **40481271** COMPLETED 00:26-00:33 (mdt78); test 40481270 COMPLETED.
- Verification jobs 40481412 (reimplementation) and 40481467 (Grn/Red check) COMPLETED.

Write-up: `C:\Users\grego\AppData\Local\Temp\claude\c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1\f1bdbd78-151f-470b-b703-dd9af9b3fecc\scratchpad\titration_curve\RESULT.md`
(the four figures sit beside it).

- **Step 1 FAIL.** In-silico mixing does not reproduce real human:mouse DNA mixtures (GSE310817): MAE 0.242, max
  0.531, against a bar of 0.05 and 0.10. Every mixture reads high: signal is concave in DNA fraction.
- **Amendment 1 FAIL.** The corrected reading misses its held-out bar: MAE 0.057, max 0.283. It matches the
  independent EPIC v2 set within 0.041.
- **Corrected human DNA fraction (not validated at the pre-set bar).** The uncorrected in-silico reads are 0.24-0.41
  higher; do not use them.

  | Array | Corrected | RNA human % |
  |---|---|---|
  | IL69B | 0.62 | 64 |
  | IL68B | 0.60 | 45 |
  | IL67B | 0.51 | 43 |
  | IL70B | 0.37 | 43 |
  | IL66B | 0.34 | 29 |
  | IL71B | 0.22 | 33 |
  | N2 | 0.05 | 4.9 |
  | C2B | 1.00 | 86 |

- **What holds:**
  - Order: Spearman 0.93 over the seven in vivo arrays.
  - Primaries lie above recurrents.
  - **N2 contains human DNA**: S is 2.3 log2 above every pure rat array. Arrays are not index-multiplexed, so this
    is separate from the RNA index-hopping question. Infiltration cannot be told apart from carry-over.
- **Input defect: GSE310817 batch 1's Red IDATs are byte copies of its Grn IDATs** (13/13 pairs; batch 2 is fine).
  - With green only (a sensitivity check, not pre-registered), both verdicts still FAIL, the order is unchanged, and
    the corrected values drop by at most 0.06.
  - **pipeline_v2 also stages GSE310817**, so tell that benchmark.
- **Open, for Greg:**
  - Whether to add an Amendment 2 for batch 1 (suggested: keep the pre-registered run and disclose the green-only
    check).
  - `ANALYSIS/methylation_curve/` is untracked but not git-ignored.
  - Cosmetic label fixes in figures b and c.
- **Paths:**
  - Grid results: `/rs/rs_grp_oschome/go2432/u251_meth/curve/results/`; verification in `curve/verify/`.
  - Pre-registrations: `curve/PREREG_curve.md` (sha `00d0f72b...`) and `PREREG_curve.v1.md`.
  - Code: `ANALYSIS/methylation_curve/` (not committed).

## 2026-09-27 01:00 EDT (grid clock): SPREAD write-ups DONE — H1-A met on its letter only, H1-B refuted, clone level empty, E1 still the open gate

Greg's contralateral-spread question is written up. Two documents, both untracked, in
`MBR/nature_draft/`: **`SPREAD_RESULTS.md`** (the full internal note: question, H-a/H-b/H-c, H1-A/H1-B, every line
of evidence with its number and criterion outcome, the decision under the pre-registered rule, the claims allowed
and forbidden, the design implications, and the scaled barcode design with its power) and **`SECTION_spread.md`**
(the manuscript-ready Results subsection, its Discussion paragraph, figure list with captions of 60 words or fewer,
and the sentences that must not enter the manuscript). `MBR/nature_draft/manuscript.md` does not exist yet, so
nothing was inserted and no CHANGES note was needed; `SECTION_spread.md` carries the exact insertion marker and
says to record it in `CHANGES.md` when the SCIENCE workflow's manuscript lands.

**PREREG_spread.md hashes.** Frozen base, sections 0-7: `454a24f3b3d6c438e2be6f0f498f5b11e4d90f2b7424204e9f6e40bae78be47d`
(41,412 B). With the 23:15 section-8 entry: `96131616a574eb4040edab18e934ee1fc73c95b277cae266a2dc8d70da7f2dd9`.
Current file, with the section-8 entry naming E1's second corpus:
`575d713b2d8f7d17c5bded9fcb69f2932f04e9320e52bc3690941f56814a10e3`; the first 41,412 B are still the frozen base and
the grid copy is byte-identical.

**The answer, in the order the pre-registration asks for it.**
- **H-c (species-assignment artefact): formally unrejected, because E1 has produced no number.** E2 and E3 are met:
  chromosome Y at 112 counts per million human counts against 2.9 in the floor library (38-fold), donor genotype
  56/59 (null 37 %, p = 0), chrY DNA dosage 0.879, CNV pattern r 0.721 at z 8.4. The internal floor is large and
  internally inconsistent (a quarter of N269B's human counts by phi-scaling, two-fifths by the Y mixture).
- **H1-A: met on its letter, weak by weight.** (i) whole-profile d 0.153 against a 99.5th percentile of 0.028;
  (ii) state-composition TVD 0.102 against a 1,000-draw maximum of 0.069. Permitted wording: "not the core
  diluted". **NOT permitted: "shifts toward an infiltrating programme"** — (iii), the only statistic with a
  pre-named direction, sits inside its band (p 0.56) and its leading-edge component moved the opposite way.
- **H1-B: not supported**, in the direction the literature predicted. Delta -0.028 (nearer the primaries); the true
  labelling ranks 15th of 20 on both statistics, p 0.75.
- **C1 and C2: both negative.** No contralateral-only arm event (the core carries all four pre-named arms at three
  to six times the measured bound U99 = 0.0254); labelling p 0.35. C2 as pre-registered gave Fisher p = 9.4e-10 in
  Greg's direction, and that was the rat-read floor: sharing tracks each tumour's human share at Spearman -1.00 and
  removing floor-origin variants leaves 0.830 vs 0.828, p = 0.93.
- **H-a vs H-b: not separated by anything in this study, as the pre-registration said in advance.** Nothing about
  LITT (both contralaterals are untreated animals). One animal = a case observation; one of two positive gives an
  exact interval of 1.3-98.7 %, so no rate exists.

**New since the 00:36 blocks.** Report **40481262 COMPLETED 00:57**: the third registered matrix
(bayesprism_gbm_ct6, 40480792_12) is in `rna/report/criteria.json` and the figures were rewritten. It fails its
negative control too — a pure human culture fits as 0.288 tumour with 0.263 pericyte and 0.154 myeloid, i.e. 71 % to
cell types a human-only stream from a rat host cannot contain, at the best fit of the three (r 0.63). Gates nothing;
lowers confidence in all three deconvolutions. **CURVE has landed** (`u251_meth/curve/results/`, written 00:33):
N2 in-silico 0.318 (0.256-0.344), mouse-calibrated **0.049** (0.039-0.060), bespoke two-mode 0.198 unchanged —
but **both carry `validated: FALSE`** against that workflow's own criteria (step 1 max |error| 0.531 vs 0.10;
Amendment 1 leave-one-batch-out MAE 0.0566 vs 0.05, max 0.283 vs 0.10). What passes there: the C2B anchor, the rat
leave-one-out floor (max 0.032), and rank agreement with the RNA share (Spearman 0.93, p 0.007). So **the ordering
is trustworthy and the level is not**: the calibrated reading would close the sixfold DNA-vs-RNA gap exactly and the
in-silico reading would widen it. No cell percentage is quoted from either platform.

**Jobs and ETAs (grid clock 00:59).**

| job | what | state at 00:59 | ETA 27 Sep EDT |
|---|---|---|---|
| **40481264** | E1 public array 0-21 (16 GSE53960 + 6 PRJNA627944 controls at 2x100 and 2x150) | tasks 0-15 COMPLETED, 16-21 RUNNING (30 min) | ~01:15-01:45 |
| **40481265** | E1 study libraries re-sorted at 50 bases single-end (matched chemistry) | tasks 0-6 COMPLETED, 7-9 RUNNING | ~01:10-01:30 |
| **40481266** | E1 human STAR/salmon index rebuild + IL64B reproduction (nf-core/rnaseq 3.22.2) | RUNNING 33 min | ~01:45-02:30 |
| **40481267** | E1 quantification of the public human streams (afterany 40481264:40481266) | PENDING | ~02:30-03:15 |
| **40481268** | E1 collect, phi_g, apply the criterion | PENDING | **~03:00-03:30** |
| 40480752 | XVAL pileup (the dosage check on N2's deeper 3p/4p/4q/12p losses) | PENDING | other workflow |
| 40481316 → 40481317 → 40481327 (held) | pipeline_v2 methods benchmark: the other input to the absolute DNA fraction | see the 00:36 block | 40481317 ~03:00-05:30 |

**Result paths.** `MBR/nature_draft/SPREAD_RESULTS.md`, `SECTION_spread.md`, `DEVIATIONS_spread_copy.md`.
Grid: `/rs/rs_grp_oschome/go2432/u251_science/spread/` — `rna/report/` (criteria.json + figures a-c, rewritten
00:57), `rna/iii_singscore_R/`, `cibersort/`, `dna/c1/`, `dna/c2/`, `dna/dna_evidence_N2.tsv`, `floor/` (E1, in
flight; read `floor/REPORT.md` and `floor/tables/floor.json` when 40481268 ends). Curve:
`u251_meth/curve/results/estimates.tsv`, `validation.json`. Code: `ANALYSIS/spread/` plus new
`ANALYSIS/spread/design/barcode_power_sim{,2,3}.py` (seed 20260927). Session scratch `spread/`: `RESULTS_rna.md`,
`RESULTS_dna.md`, `DEVIATIONS.md` (now with **W-D1**, this write-up's two additions), `DESIGN_IMPLICATIONS.md`,
`CLONAL_LIT.md`, `barcode_power_sim*.json`. Nothing committed; `MBR/` stays untracked.

**The scaled design, for the POWER workflow.** Barcoded U251N through LITT; compartments core, ipsilateral rim
beyond the margin, contralateral hemisphere, recurrence, plus a pre-LITT sample of the same animal; the injectate
sequenced per batch as the input null. Nulls N1 (contralateral barcodes = multinomial draw from the core), N2 (the
recurrence = draw from the same animal's pre-LITT core), N3 (bottleneck, measured per compartment). **Simulated
per-animal power (stated inputs, not estimates; no published U251N bottleneck exists): statistical power is not the
binding constraint.** Five recovered clones at 200 reads detect a bottleneck with certainty, and a reweighting with
no clone lost still gives 0.76-0.94. What binds is carry-over and recovery: with the contralateral pool sequenced
to 10,000 reads the test tolerates 90-95 % carry-over (power 0.93-1.00) and collapses at 99 % (0.06-0.22), and pure
carry-over is exactly the null (realised size 0.04-0.06 at alpha 0.05). Across animals it is a binomial: at
per-animal power 0.95, rejecting in at least n-1 of n animals has probability 0.99 (n=4) to 0.91 (n=10); at 0.80,
0.82 to 0.38. So the rats-needed figure is driven by N3, the recovery fraction, not by the statistics.

**Next steps.** (1) When 40481268 lands, read `floor/REPORT.md`; if the criterion rejects H-c, the two write-ups
need only their E1 rows and the H-c paragraph changed — the wording is already staged as conditional. If it does
not reject, PREREG section 5 requires Parts 2 and 3 to be reported as uninterpretable for origin, and both
documents must be cut back to Part 1. (2) Re-run Part 2 with the external phi_g
(`spread_rna.py --phi-tsv floor/tables/phi_g.tsv`; the switch self-test passed as 40481252) and record whether any
verdict moves. (3) XVAL 40480752 settles whether N2's deeper losses are a subclone or low-signal distortion. (4)
Ask Raj for the one cheap decisive control: a cutting-order blank — a naive rat hemisphere processed with the same
blade immediately after a tumour-bearing one — plus the adjacent sections for a human-specific stain and the
tared-vial weights. (5) The DiG/LiG lists (Chanoch-Myers 2026) are still the single most informative addition to
the state level. (6) `SECTION_spread.md` lists four clinical claims for which this project's verification run
produced **no citable source** (recurrence geometry beyond the treated volume, pseudoprogression after LITT,
hyperthermia and invasion, CSF liquid biopsy) and two that were refuted and must not be cited; those citations
have to be found before the Discussion can carry a stronger clinical frame.

## 2026-09-27 00:47 EDT (grid clock): multi-omics build VERIFIED; three items running (prep_reads, stage_refs, rnaseq_graft); plan in `ANALYSIS/MULTIOMICS_PLAN.md`

Plan (what exists, what runs, triggers, rejected items, what each result can claim and its outside check):
`ANALYSIS/MULTIOMICS_PLAN.md`. Scripts: local `ANALYSIS/multiomics/{prep_reads,stage_refs,rnaseq_graft}/` (untracked,
not committed); grid copies and all outputs under `/rs/rs_grp_oschome/go2432/u251_multiomics/` (`M`). Nothing was
cancelled or resubmitted: the verification found no high or medium defect in the three submitted items.

| job | pipeline @ release | question | state 00:46 | ETA (EDT, 27 Sep) | read when it finishes |
|---|---|---|---|---|---|
| 40480709 (array 1-10) | prep_reads (custom awk/pigz) | graft-only pairs, R2 clipped 3 nt once; clipped human R2; `both` names; index pair per library | 9 of 10 COMPLETED, all checks OK; C2B at its final checks | C2B ~01:05-01:25 | `M/reads/checks/C2B.tsv` (status OK?), log `M/reads/logs/u251_prep_40480709_1.out` |
| 40480710 | collect_prep (afterany 40480709) | 10 of 10 OK? READY; index pairs; N269B hopping verdict | PENDING | < 1 min after C2B | `M/reads/READY` (exists?), `M/reads/CHECKS.tsv`, `M/reads/INDEX_HOPPING.txt` |
| 40481320 | stage_refs (custom, network) | GENCODE v44, known sites, VEP 115, rnafusion 4.1.3 bundle (v46), hg38ToRn7 chain, REDIportal | RUNNING (rhi2); containers, gencode, chain PASSED; known_sites started 00:46 | ~04:00-07:00 (limit 12:33; on TIMEOUT resubmit `sbatch -D M/refs M/stage_refs/run.sbatch`, done steps skip) | `M/refs/STATUS.tsv`, `M/refs/REFS_CHECKS.tsv` (any FAIL?), `M/refs/.done/`, `M/refs/.recheck_submitted` (recheck id) |
| 40480725 | nf-core/rnaseq 3.22.2 (NXF 25.10.2), gate on READY | U251 profiles without the rat `both` bin; clip effect; IL64B/N168B graft reads; human indices saved | PENDING (afterany 40480710) | ~05:30-07:30 | `M/rnaseq_graft/NEXTFLOW_EXIT` (0? `gate_closed`?), log `M/rnaseq_graft/logs/u251_rnaseq_graft_40480725.out` |
| 40480726 | rnaseq_graft checks (afterany 40480725) | hard checks 0-4, 9, 10; reports 5-8 | PENDING | ~06:00-08:15 | `M/rnaseq_graft/CHECKS.md`: check 5 (clip material? triggers the host rerun), check 7 (rat share vs 1.9-5.1 % / 5.1-8.0 %), check 8 (IL64B/N168B verdict; "unclassified" triggers the competitive alignment) |

Found in the verification (details in the plan, Section 4):
- **Index hopping into N269B cannot be excluded.** Index pairs from the read headers: N269B = CTGAAGCT+AGGATAGG shares
  its i7 with IL66B and IL69B and its i5 with NL71B and C2B (C2B = TAATGCGC+AGGATAGG), all human-rich; the ten libraries
  were indexed combinatorially (3 i7 x 4 i5). A rough ceiling from IL64B (every one of its 314 k graft reads counted as
  hopped) is about 0.5 M of N269B's 3.82 M graft reads; it assumes similar per-molecule hopping for the two recipients,
  which is not known. N168B's i7 is also shared with human-rich libraries (IL68B, NL71B); its i5 is unique.
- **Two lanes, not one.** Every library's reads are split about evenly between lanes 1 and 2 of AAAMFFWHV. The
  23:45 block's "one lane" and `ANALYSIS/SAMPLE_KEY.md`'s "lane 1" describe only the first read; neither file was
  edited here. The GEO template should say two lanes if it names one.
- **Raw FASTQs are on SRA**: SRR39547363-SRR39547372, spot counts public (N269B = SRR39547367 = 78,075,952 = the
  xengsort total). The later "raw FASTQ re-download" item no longer needs the Drive token; not submitted.
- Every param name in the submitted and planned param files exists in its release's schema (rnaseq 3.22.2 = 3816d48,
  rnavar 1.3.0, rnafusion 4.1.3, rnasplice 1.0.4; schemas fetched from GitHub at the tag). Strandedness reverse in all
  11 rnaseq rows; the R2 clip is applied once upstream (gate A passed in all ten) and TrimGalore does not clip again.
- **Not built** (selected; see the plan, Section 5): human_de_null pass A and integration_rna_prep (no dependency, can
  go now), invasion_niche and pass B (after the rnaseq_graft checks), rnavar 1.3.0, rnafusion 4.1.3, rnasplice 1.0.4.
  The methylation integration (I1-I5) needs integration_rna_prep and human_de_null pass A, so nothing queued will
  start it. METH_FINAL itself holds (bespoke `MT/results/summary.json` 22:10:52, REPORT COMPLETED, no `u251_meth` job);
  whether the integration should use the bespoke result or pipeline_v2's is not settled.
- The rnafusion launcher, when built, must rerun `python3 M/stage_refs/refs_util.py preflight --base M/refs/rnafusion
  --ver 46` on its own head node before `nextflow run` (rnafusion 4.1.3 checks the last file its directory walk
  returns, and CephFS listing order differs by node).
- `beta_values.rda` (the lab's beta matrix, untracked) was not git-ignored; added to `.git/info/exclude` (local only).
  Still untracked and not ignored, for Greg to decide: manuscript PDFs under `manuscript/`, `Submission_Package_ForRaj*`,
  `NONC/`, `REVIEW/`, `SLIDES/*.pptx` copies. Nothing committed by this workflow; `ANALYSIS/multiomics/` and the plan are
  untracked.

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
| 40481256 / 40481257 u251_ref_prep_{Rnor6,mR7} | Rnor_6.0 (Ensembl 101) / mRatBN7.2 (Ensembl 110) + `bismark_genome_preparation --bowtie1` (single-threaded bowtie-build: the long step) | ~03:00-04:30 (bmax 492 M: ~7 sort blocks per conversion; 2 done in 35 min) | `ref/rat_genome/<asm>/genome/Bisulfite_Genome/PREP_DONE` |
| 40481258 / 40481259 u251_ref_align_{Rnor6,mR7} | Needhamsen's exact `bismark --bowtie1 -n 1 -l 28 -f` on their FASTA (+ `--ambiguous --un`); mR7 also the R-expanded deviation; spiked-mismatch test; Rnor6 vs Zhou's rat mapping | ~04:00-06:30 | `ref/e5/E5_Rnor6.tsv.gz`, `E5_mR7.tsv.gz`, `E5_mR7_Rexp.tsv.gz`; `ref/rat_alignment/checks/` |
| 40481260 u251_ref_collect | writes `ref/e5/DONE` (built / not built + reason), merged hits, PREREG lists | ~06:30 (latest ~10:00 if nodes are slow) | `ref/e5/DONE`, `ref/rat_alignment/epic_v1_rat_hits.tsv`, `ref/rat_alignment/lists/E5-*.txt` |

Done by 00:23: GEO OVERALL OK (28/37/12 arrays, address counts 1,051,815 / 1,051,943 / 1,105,209); SeSAMe carries
`rattus_norvegicus` (Rnor_6.0, 25,672 unmasked of 866,553); tools = Bowtie 1.1.2 + Bismark tag 0.14.5 (its script says
v0.14.4), smoke PASS, no fallback.

Found while writing the jobs (details in `meth_bench/INTERFACE.md` section 7; a local end-to-end test on Rnor_6.0 chr20
backs the last three):
- Needhamsen's Additional file 1 (866,894 sequences) is rc(AlleleA) with every probe `R` already set to `A`; Bowtie
  never saw an `R`, and after Bismark's C->T read conversion every R expansion collapses to the same read. The
  METHODS.md E5 step 3 premise does not hold for their input; E5_mR7_Rexp is still built as INTERFACE.md specifies.
- The published command tolerates **2** mismatches per probe (at most 1 in the 28-nt seed), not 1: Bowtie 1 rounds
  qualities to 30 by default, so `-e 70` admits two. Spike test on chr20: 2 outside the seed map 113/113.
- Bismark 0.14.5's NM tag counts C->T conversions; E5 `mismatches` come from Bowtie's converted-space count instead.
- Coordinate/colour rules agree 100 % with Zhou's Rnor_6.0 mapping on the probes both place at the same locus.

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

## 2026-09-27 05:20 EDT: WEEKLY LIMIT HIT (resets Sep 30 22:00 ET); Fable driving solo; Greg's priority order: LITT first, then the real-time video inference demo, then the corpus screen

- **Subagent workflows are unavailable until Sep 30 22:00** ("weekly limit"); work continues on the main thread only. Also: an Opus 5.5 safety classifier withheld one reply (about the cross-validation results) and flagged two science-workflow verification agents; nothing to redo, just do not re-emit that reply.
- **Workflow end states (all U251 runs are now finished or dead):** SPREAD v2 DONE; XVAL DONE (`MBR/nature_draft/XVAL_RESULTS.md`; hostile-verified); CURVE DONE; TIMING DONE; BENCH DONE as a workflow (pipeline_v2: 40481316 DONE 00:58, 40481317 DONE 02:50, **40481327 still held/pending**; results dir has c0/c1/c2e/c2f/e1-e5/classes but no REPORT.md yet); multi-omics: 15/16 agents (final verify failed at the limit) with **nine builds submitted**; SCIENCE: research, inventory, three proposals, PREREG and both analysis submissions DONE, **manuscript/referee/revision NOT written** (Fable writes it); POWER: TARGETS.md (sha256 88b5ac10...) and pilot.json DONE, **curves/figure/SCALED_EXPERIMENT.md NOT done** (u251_pwr_asm job pending on the grid); nf-core/methylarray **DONE 00:41 (40480724)**.
- **nf-core/methylarray result (PR #41 + patch 0ef3365e...):** six tumours kept (unmasked pOOBAH fail 1.4-4.4 %); 719,866 probes after QC, 702,636 annotated; **true split: 46 DMPs at FDR < 0.05, 43 significant CpGs, 9 DMRs; relabelling rank 1 of 10 (p = 0.10, the floor of a 10-split null); lambda_GC 1.88 for the true split vs 0.66-1.49 for the others** (inflation = tumour fraction, which this arm does not model). The bespoke arm with fraction as covariate gives 0. Read Q5 as "arm confounded with fraction"; the nf-core arm confirms the bespoke arm's headline once fraction is ignored.
- **Grid jobs running at 05:17 (from the finished workflows):** E1 floor: 40481264/40481265 arrays DONE, 40481266 DONE 03:13, **40481267 floor_quant RUNNING (2 h), 40481268 floor_final PENDING** -> `/rs/rs_grp_oschome/go2432/u251_science/spread/floor/REPORT.md`; multi-omics: `u251_rnaseq_graft` 40480725 RUNNING 4.3 h + check 40480726; `u251_rnavar` RUNNING (seq2HLA, STAR index) + `u251_rnavar_integ`; `u251_rnasplice` RUNNING (DEXSeq) + check + null launch; `u251_de_null_B`, `u251_invasion_niche`, `u251_irp_graft` PENDING; rnafusion has NEXTFLOW_EXIT_PASS1/PASS2 files (check them); **stage_refs 40481320 FAILED at 05:10 on the rnafusion reference (122 GB bundle) and REDIportal unavailable**: rnafusion's reference is not staged; SCIENCE jobs `sci_ai` RUNNING, `sci_dilution`, `sci_dna_integrate`, `sci_fdna`, `sci_nfcofactor` PENDING, `pub_count` RUNNING, `pub_noiseq` PENDING (outputs `/rs/rs_grp_oschome/go2432/u251_science/{rna,dna}/`); POWER `u251_pwr_asm` PENDING; BENCH `u251_v2` 40481327 PENDING.
- **Fable's plan (main thread):** (1) write `MBR/nature_draft/manuscript.md` from the verified pieces (SPREAD_RESULTS, XVAL_RESULTS, CURVE RESULT, science PREREG/proposals, TARGETS) with \PENDING slots for the running jobs; (2) compute the rats-needed curves from TARGETS + pilot.json and write SCALED_EXPERIMENT.md + the closing section; (3) collect E1/multi-omics/BENCH/SCIENCE results as they land.

**SPREAD v2 DONE (wjy1l2zie / wf_e44b06bf-0e6, 2026-09-27 ~01:05): Greg's hypothesis tested against pre-registered nulls (PREREG_spread.md sha256 96131616...; frozen base 454a24f3...). Write-ups `MBR/nature_draft/SPREAD_RESULTS.md`, `SECTION_spread.md` (hostile-verified; nine disclosure issues fixed in place).**
- **H-c artefact:** rejected in substance (chrY 112 vs 2.9 per million in the floor library; genotype 56/59 vs null 37 %; chrY DNA dosage 0.879; CNV r 0.721, z 8.4), FORMALLY open until the external floor E1 lands: jobs 40481264-40481268 (GSE53960 F344 brain + PRJNA627944 SD control brain through the study's own xengsort image and index; final report ~03:00-03:30 at `/rs/rs_grp_oschome/go2432/u251_science/spread/floor/REPORT.md`).
- **H1-A (contralateral U251 distinct from its core): met on its letter only.** Whole-profile distance 0.153 vs carry-over-null 99.5th pct 0.028; CIBERSORT state TVD 0.102 vs null max 0.069. Allowed wording: "not the core diluted", nothing more: the only directional statistic (Ivy GAP leading edge minus cellular tumour) sits inside its band (p 0.56) and its leading-edge part moved the wrong way; another animal's tumour (IL68B) differs from IL69B by the same TVD (0.103); deconvolution fits r 0.13-0.30 and fails its negative controls; the null carries an anti-conservative asymmetry (the two floor estimates disagree 1.7-fold, ~23 % unsubtracted floor counts), and (ii) was extended to 1,000 draws after the 200-draw result was seen (both disclosed).
- **H1-B (closer to recurrences than primaries): NOT SUPPORTED** (delta -0.028; true labelling 15th of 20; p 0.75 on both statistics).
- **Clone level: C1 NEGATIVE** (no contralateral-only CNV event; IL69B and C2B carry the four named arms 3p/4p/4q/12p at 3-6x the noise bound U99 0.0254; labelling p 0.35). **C2 NEGATIVE after floor control:** the pre-registered variant-sharing test gave Fisher p = 9.4e-10 in Greg's direction, but sharing tracked each tumour's human share at Spearman -1.00 (rat reads aligned to GRCh38 look like variants and are commonest in the low-human recurrences); removing variants the floor-only IL64B carries gives 0.830 vs 0.828, p = 0.93. **Every "recurrence-like" signal found is tumour fraction or host-read content, collinear with arm.** That is itself a citable methods warning for xenograft clonal analyses.
- **Spread vs carry-over: not separable by this design** (same slices, blades and slide). Nothing about LITT (both contralaterals untreated). One animal: a case observation.
- Flag: `ANALYSIS/spread/` and `ANALYSIS/methylation_curve/` are untracked but not ignored in the public repo (code only, but one `git add -A` publishes them).

**CURVE DONE (wc2ib5g82 / wf_480455e9-000, 2026-09-27 ~01:00):** full run 40481271 (66 arrays). Result `scratchpad/titration_curve/RESULT.md`; grid `/rs/rs_grp_oschome/go2432/u251_meth/curve/results/` (validation.json, estimates.tsv, four figures); PREREG v1 sha256 0bfa0b2e..., Amendment 1 sha256 00d0f72b... (added after the test run showed the problem, before any tumour array or N2 was read).
- **The pre-registered curve FAILED its test:** in-silico mixing of pure arrays does not behave like real DNA mixing (GSE310817, 28 real human:mouse mixtures: MAE 0.242, max 0.531 vs bar 0.05/0.10; GSE299969 also fails). The signal is concave in DNA fraction: a real 50:50 mix carries 79-89 % of the pure-human signal, a software mix 53-62 %.
- **The corrected reading (Amendment 1) also FAILED, narrowly:** leave-one-batch-out MAE 0.057, max 0.283 (23 of 28 within ±0.10; misses depend on the host tissue); the independent EPIC v2 set (not in the fit) lands within 0.041 on all four.
- **What survives:** the ordering (DNA vs RNA human share Spearman 0.93 across seven in-vivo arrays, p = 0.007; all three primaries above all three recurrents, descriptive) and **N2 contains human DNA** (its statistic is 2.3 log2 above the highest of ten pure-rat arrays). Corrected readings (not evidence: the correction failed its own check): IL69B 0.62, IL68B 0.60, IL67B 0.51, IL70B 0.37, IL66B 0.34, IL71B 0.22, **N2 0.049**, C2B 1.00. This points to the earlier "N2 DNA 20 % vs RNA 3 %" discrepancy as the bespoke two-mode estimator's inflation at low fraction (the same concavity), not biology; unproven.
- **GEO deposit error found:** in GSE310817 batch 1, all 13 Red IDATs are byte copies of their Green IDATs (batch 2 is fine). Green-only rerun: both verdicts still FAIL, order unchanged, corrected values move <= 0.06. **The BENCH workflow (wsg4cdgp5; pipeline_v2 jobs 40481316 -> 40481317 -> 40481327) uses GSE310817 for its C1 check and must use batch 2 or the green channel only: check its C1 when it lands.** Worth reporting to GEO / the Zhou lab (Greg's call).

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
  - **XVAL tests step (2026-09-27 00:05-01:20 EDT; DEVIATIONS 14-27; PREREG sha256 unchanged):**
    - The bam-check FAIL is only its own over-strict `read1 == read2` line: singletons are 0.01 %; the bound, flagstat,
      md5 and contigs all pass.
    - **Tests 2-4 DONE** (40480794 `xval_deseq`, 40480795 `xval_tests`, 40481406 `xval_desc`, 40481374/40481408 `xval_fig`).
    - *Test 2:*
      - Both positive controls pass: CDKN2A/B CPM <= 0.015, and the 273-point dosage slope is 0.54 (p 7e-6).
      - The pre-registered rule passes N2's 3p (p 0.0083) and 4p (0.0120; the second draw gave 0.0134, borderline).
        Distal 4q and 12p are not confirmed.
      - But an outside check (descriptive) shows **N269B's arm-level expression does not follow N2's copy number
        anywhere**: slope 0.17 (p 0.45), against IL69B 0.97 (p 6e-7). So the passes are not dosage evidence, and the
        N2-only losses stay unconfirmed.
    - *Test 3:* promoter beta vs expression gives rho -0.33 to -0.36 in all seven (pass). Delta-delta gives rho -0.003,
      rank 8/10, 0 genes selected, so there is no coordinated change.
    - *Test 4:* the RP block has no DNA support, and its unadjusted promoter hypomethylation vanishes with fraction
      modelled. The fraction-adjusted up-list is "DNA-supported" by the letter (CNV +0.013 log2, p 0.043, three genes).
    - *Test 1, primary (duplicates excluded), DONE* from the primary-only run 40481651, outputs in `xval/rna/nodup_only/`:
      - **Set A is not testable:** the 59 rs probes have at most 2 RNA reads.
      - Set B: 543 pairs, 2.95 % contradictions (B + B* 1.39 %). All 16 come from five reproducible germline sites
        (FKBP9, PGM1, EP400, ELL2, RBM12B).
      - RNA-vs-RNA agreement 99.86 %.
      - The pre-registered e_s is 1.5-2.7 %, above the 1 % bar, and is driven by those five sites; the descriptive rate
        with them set aside is 0.25-0.40 %. The kmin table is written.
    - **Pending:** duplicates-kept chunk 30 of 40481325 (3-h limit, by about 03:45 EDT), then `xval_geno` 40481326
      (afterany) reruns everything into `xval/rna/` with the sensitivity; 40481699 re-draws the primary-only figure.
      - If chunk 30 times out, the duplicates-kept sensitivity is missing that chunk and 40481326 STOPs. The primary
        numbers stand either way.
    - Write-up: `MBR/nature_draft/XVAL_RESULTS.md` (= scratch `xval/RESULTS_xval.md`, grid `xval/results/RESULTS_xval.md`).
    - Deviations 14-30 are in scratch `xval/DEVIATIONS.md`, and the grid copy is `xval/DEVIATIONS_xval.md`.

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
