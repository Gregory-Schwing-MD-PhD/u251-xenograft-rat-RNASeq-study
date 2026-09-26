# Working checkpoints: U251 LITT recurrence study

Newest first. Everything for this study lives in this directory and this repository (Greg, 2026-09-26: "this project
should be contained in the u251 dir and the respective git repo"). The checkpoints below were first written into
spinesurg-ct-nnunet/docs/RESUME_2026-09-22_PM.md by mistake and were moved here the same day.

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

Status 06:25 EDT (grid clock): 40467314 RUNNING, DESeq2 finished for the five tumours (GSEA, plots and report still to
run); 40467315 and 40467316 PENDING on it. Expected finish of the chain: about 08:00–09:00 EDT.

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
