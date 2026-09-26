# Working checkpoints: U251 LITT recurrence study

Newest first. Everything for this study lives in this directory and this repository (Greg, 2026-09-26: "this project
should be contained in the u251 dir and the respective git repo"). The checkpoints below were first written into
spinesurg-ct-nnunet/docs/RESUME_2026-09-22_PM.md by mistake and were moved here the same day.

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
