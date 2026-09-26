# Abstract 418 talk, v4: one story, 13 slides, about 5:21

CNS 2026, SSTU02, Monday 2 November 2026, 7:00–7:06 AM (six minutes). This file is the plan for the v4 rebuild of
`04_build_cns_template_deck.py`. It starts from v3 (`deck_v3_text_dump.txt`, 25 slides) and the results added since.
Built 2026-09-26 as `SLIDES/09_build_deck_v4.py` → `CNS2026_Schwing_Abstract418_CNStemplate_v4.pptx` (13 slides, 4
backups); the departures from this plan are listed in section 7. All paths are relative to the repository root,
`u251-xenograft-murine-RNASeq-study/`.
Every number below was read from the file named beside it on 2026-09-26.

## 1. The three candidate stories, scored

Each design is scored from 1 to 5 on four criteria, for a total out of 20.

| design | compelling | honest (n = 3 vs 3, one-tumor dependences) | fits 6 min | numbers traceable | total |
|---|---|---|---|---|---|
| A. The confound is the reason for the next experiment (14 slides, 338 s) | 4 | 5 | 4 | 5 | **18** |
| B. One finding, four attempts to break it (16 slides, 347 s) | 3 | 5 | 3 | 5 | 16 |
| C. Program, then state, then drug, with the limits at the end (14 slides, 337 s) | 5 | 3 | 4 | 4 | 16 |

**A** answers a surgeon's questions in order: what comes back, can I trust it, and what do I do next. The tumor-content
confound is shown on stage as the reason for the next experiment, not hidden in a caveat. It reports the robust
quantity (the direction and the nominal p) rather than the fragile one (a single FDR q). Weak spot: its slide 9
(hypoxia and cell cycle) is off the main line, and it runs 14 slides.

**B** has the most memorable structure. But five slides, about two minutes, go to statistics in front of a 7 AM
neurosurgical room. Its eight-row scorecard cannot be read in 25 s, and the rat-read test gets a whole slide when the
tumor-content question matters more. It is also the longest of the three (16 slides, 13 s spare).

**C** is the most persuasive clinically, but its order oversells:
- The drug is ranked first two slides before the audience learns that it drops to third without IL66B.
- It keeps glycolysis as a finding, although its nominal p reaches 0.13 when a tumor is left out (`loo_sets.tsv`).
- It reads the Varn GLASS "stem-like" call as the line's baseline. The rat-brain controls read the same, 0.88–0.99
  (`ANALYSIS/cibersort/results/v104/fractions_v104_varn2022_*.tsv`), so the call says nothing about the tumor.
- Two of its numbers exist only in `RESUME.md` or in superseded `facts.json` fields.

**Final plan: A's arc, with grafts.**
- From B: the single-gene instability line on slide 4, the "does it survive?" handoff on slide 5, and the hold-out
  chart as backup B2.
- From C: the read-sorting line on the model slide, the program, state and drug beats in "What this means", greying
  out the mitochondrial line on the subtype chart, and the risk that the printed abstract contradicts the talk.
- Removed from A: its "what else moves" slide, which becomes backup B3. That cuts the talk from 14 slides to 13 and
  leaves 39 s of slack.

## 2. The story

LITT kills the core of a glioblastoma, but the tumor grows back from cells that survived the sublethal margin, and what
those cells regrow into has not been read. In an MRI-guided LITT model of a human glioblastoma (U251N) in rat brain,
with three tumors taken before ablation and three after regrowth, every read was sorted by species so that the tumor
could be read on its own. The regrown tumor had turned its protein-synthesis machinery down: the six most depleted of
8,869 gene sets are all translation or amino-acid stress, and all 80 translation-initiation genes sit in the bottom
6.6 % of the ranked transcriptome.

We then tried to break that result.
- The direction held in all 41 leave-one-tumor-out re-runs, and rat reads cannot produce it.
- Its FDR q moved from 0.001 to 1.0.
- The recurrent samples held less tumor, so six animals cannot separate recurrence from tumor content.
- The astrocyte-like fall is the change that best survives adjustment for tumor content.

The same signature ranks ciclopirox first of 54 clinically available compounds, and third when one recurrent tumor is
held out. It blocks the eIF5A-activating step that U-251 depends on, and LITT holds the barrier open around the
ablation for weeks. So the next experiment is built around the confound: more animals matched on tumor content,
translation measured directly, and ciclopirox given in the barrier window.

**The one sentence to leave with:** "After LITT, the glioblastoma that grew back in this model had turned its
protein-synthesis machinery down. The direction survived every test six animals allow, its significance did not, and
it points to ciclopirox, given in the weeks LITT holds the barrier open, as the experiment to run."

**Why it is compelling and honest:**
- It follows the three questions a surgeon asks, in order: what comes back, can I trust it, and what do I do next.
- The stress tests form the middle of the talk (slides 6–8), not a limitation on slide 24.
- It separates what held (the direction, the nominal p, the drug in the top three) from what did not (the FDR q, the
  size of the fall, the single-gene list). It says out loud that tumor content is unresolved.
- The drug stays a hypothesis with a delivery window.
- Nothing is cut to make the story cleaner. The cuts remove side results (the roadmap, the methods slide, the culture
  PCA, the hubs, the prior-art table), and every fragile claim is reworded on screen.

## 3. Running order

| n | title (≤ 60 chars) | from v3 | action | s |
|---|---|---|---|---|
| 1 | Title slide (Abstract 418 title, unchanged) | 1 | keep | 10 |
| 2 | LITT kills the core; glioblastoma returns from the margin | 3, 4 | merge | 24 |
| 3 | The model: a human tumor in rat brain, ablated under MRI | 5, 6 | merge | 20 |
| 4 | Few genes move alone, so we read whole programs | 9, 10, 7 | merge | 18 |
| 5 | After LITT, the regrown tumor turns protein synthesis down | 11, 13 | edit | 30 |
| 6 | Leave any tumor out: the direction holds, the FDR does not | 15 (half), 11 notes | new | 32 |
| 7 | Recurrence after LITT, or less tumor in the sample? | 15 (half), 6 | new | 38 |
| 8 | The astrocyte-like program falls, and survives adjustment | 16 | edit | 25 |
| 9 | The signature points to one approved drug: ciclopirox | 18, 19, 20 | merge | 30 |
| 10 | Why ciclopirox: it hits a translation step U-251 needs | 22 (+ 21 line) | edit | 28 |
| 11 | What this means | 23 | rewrite | 28 |
| 12 | What six animals cannot tell us, and the experiment that can | 24 | rewrite | 28 |
| 13 | Acknowledgments and data availability | 25 | keep | 10 |
| | **total** | | | **321** |
| B1 | Backup: rat reads cannot produce the fall | 15 | backup | 0 |
| B2 | Backup: the whole pipeline without IL68B or IL66B | none | backup | 0 |
| B3 | Backup: hypoxia down, cell cycle up, nominally | 12, 14 | backup | 0 |
| B4 | Backup: cell-state deconvolution (CIBERSORT v1.04) | none | backup | 0 |
| B5 | Backup: what is already known about the top twenty | 21 | backup | 0 |

**Cut:**
- v3 2 (roadmap): slide 2 now ends on the question.
- v3 7 (methods): one sentence moves to the slide 4 notes.
- v3 8 (dish to brain): a single culture library, outside every test.
- v3 13 (one-sided): its point is now one line on slide 5.
- v3 17 (hubs): unstable in both hold-outs, and not on the story line.

**Cut order if running long:**
1. Slide 4 down to its first and last line (saves about 6 s).
2. The third line of slide 9 (the hold-out ranks, which are on backup B2).
3. The third line of slide 10.

Never cut slides 6 or 7; they are where the talk earns trust.

## 4. Slide by slide

Slide text is set at 18–20 pt with 3–4 lines. Charts use 17–30 pt type and are drawn as dots, lines or slopes, never
sideways bars. Captions are 60 words or fewer.

### 1. Title slide (keep v3 1), 10 s
- Characterizing the Transcriptomic Recurrence Signature of Glioblastoma Following Laser Interstitial Thermal Therapy
- Authors, affiliations, session line, time and room exactly as in v3 slide 1.
- NOTES:
  - Do not read the title.
  - Give the disclosures in one breath: I.Y.L. consults for Medtronic (Visualase, the LITT system used here) and
    Monteris (NeuroBlate); there was no commercial support; the drug candidates are computational predictions, and none
    has been tested in this model.
  - Use the new cut order in section 3.
- Figure: none.
- Numbers: none. Wording is from `SLIDES/deck_v3_text_dump.txt` slide 1.

### 2. LITT kills the core; glioblastoma returns from the margin (v3 3 + 4), 24 s
- LITT ablates the core; beyond the ablative threshold lies a sublethal margin.¹
- Recurrence grows from cells that survive there; the margin's vessels and immune response have been described.²˒³
- What has not been read is what program the regrown tumor runs, and whether a drug opposes it.
- Footer: 1. Chen C, et al. J Neurooncol 2021;151:429–442. 2. Cleary RT, et al. Neuro Oncol 2026;28:1649–1661.
  3. Tao R, et al. J Immunol 2026;215:vkaf327.
- NOTES:
  - Nagaraja 2026 did a first-pass recurrent-versus-primary comparison on these libraries, so do not say "never
    profiled". What is new is reading the tumor's own reads gene-set-wide and linking them to a drug.
  - Cleary 2026 comes back on slide 10: the barrier stays open for weeks.
  - End on the question and pause.
- Figure: `figures_cns/fig_margin.png` (existing). Caption: "Schematic, not to scale: thermal dose falls with distance
  from the laser fiber; cells beyond the ablative threshold survive in a sublethal margin."
- Numbers: none. References come from `SLIDES/deck_v3_text_dump.txt` slides 3–4.

### 3. The model: a human tumor in rat brain, ablated under MRI (v3 5 + one line of 6), 20 s
- Human U251N glioblastoma grown in the athymic rat; LITT guided and watched by diffusion MRI.⁴˒⁵
- Three tumors taken before ablation, three after regrowth.
- Every read is sorted by species (89–92 % assigned to one), so the tumor is read apart from the brain.⁶
- Footer: 4. Nagaraja TN, et al. J Neurosurg 2026;145:364–377. 5. Nagaraja TN, et al. Acta Neurochir 2021;163:3455–3463.
  6. Zentgraf J, Rahmann S. Algorithms Mol Biol 2021;16:2. Sequencing data: GEO GSE338105.
- NOTES:
  - Move fast.
  - MR thermometry could not run on the animal scanner.
  - The culture sample and the three rat-brain controls are in no test.
  - Reconcile the J Neurosurg abstract ("four per arm") with the Methods (4 primary + 3 recurrent) before the talk.
- Figure: `figures_cns/ablation_512.mp4` (poster `video_poster.png`), with `overview_pipeline.png` as a strip beneath
  (both existing). Caption: "Ablation under diffusion-weighted MRI in this model (7.9 s)."
- Numbers:
  - three per arm: `ANALYSIS/graft_relation/results/per_sample.tsv` (group)
  - 89–92 % assigned: graft_pct + host_pct, 88.6–92.3 %: `ANALYSIS/graft_relation/results/per_sample.tsv`
  - 7.9 s: `SLIDES/deck_v3_text_dump.txt` slide 5

### 4. Few genes move alone, so we read whole programs (v3 9 + 10; v3 7 into the notes), 18 s
- Every recurrent tumor lies right of every primary on the first axis; 27.8 % of the variance falls between arms
  (p = 0.10, the floor for three per arm).
- 35 genes pass FDR; without one tumor the list becomes 43 or 140, keeping 21 or 30 of the 35.
- So we test whole gene sets, permuted on genes, because three per arm allow only ten splits.
- NOTES:
  - The PERMANOVA F is 1.54, against 20 % of the variance expected by chance.
  - Nominal p is shown beside FDR q because gene-set permutation is anti-conservative for correlated sets (Maleki
    2019).
  - Without IL68B, the first axis no longer splits the arms.
  - Keep the readthrough and U6 top hits from v3 slide 10 for Q&A.
- Figure: `figures_cns/chart_pca.png` (existing). Caption: "Six tumors on the first two principal components of the
  500 most variable genes (41.8 % and 26.2 % of the variance). PERMANOVA on all genes: 27.8 % between arms, p = 0.10,
  the smallest value the ten possible splits allow."
- Numbers:
  - PC1 41.8 %, PC2 26.2 %, R² 0.278, F 1.54, p 0.10, 35 genes (23 up, 12 down): `SLIDES/figures/facts.json` (pc1,
    pc2, permanova, de_total, de_up, de_down)
  - 43 genes, 21 shared: `ANALYSIS/holdout_IL68B/COMPARISON.md` §1
  - 140 genes, 30 shared: `ANALYSIS/holdout_IL66B/COMPARISON.md` §1
  - first axis no longer splits the arms without IL68B (notes): `ANALYSIS/holdout_separation/loo_separation_human.tsv`
    (pc1_splits_groups FALSE)

### 5. After LITT, the regrown tumor turns protein synthesis down (v3 11 + 13), 30 s
- Of 8,869 gene sets, the six most depleted are all translation or amino-acid stress (NES −1.93 to −1.99).
- All 80 translation-initiation genes sit in the bottom 6.6 % of 19,236 ranked genes.
- Nothing on the up side survives correction (best q 0.61). With three per arm, does it hold?
- NOTES:
  - Walk the curve once: the score steps down at every other gene and jumps at every member.
  - All six sets are nominally significant (p ≤ 0.002). Do not quote q = 0.022 as the result; it is one permutation
    seed (next slide).
  - The average fall per gene is about 8 % (ashr fold changes, manuscript S2).
  - Hypoxia signalling (HIF1 targets, Buffa metagene) falls nominally in all 41 re-runs but never reaches q < 0.29.
    It is backup B3, for questions only.
- Figure: NEW `chart_running_sum_ti.png`, the translation-initiation panel of `chart_running_sum.png` redrawn alone at
  slide size (data `figures_cns/chart_running_sum.json`), with the label "NES −1.99 · nominal p < 0.001 · 80 genes" in
  place of "q = 0.022". Caption: "Genes ranked from most up in recurrence (left) to most down (right). Ticks mark the
  80 translation-initiation genes; the line is the running enrichment score, which reaches its floor where they crowd.
  Gene-set permutation, 1,000 permutations."
- Numbers:
  - 8,869 sets: `SLIDES/figures/facts.json` gsea_n_sets; `ANALYSIS/gsea_leave_one_out/SUMMARY.md` (exact reproduction)
  - NES −1.93 to −1.99, nominal p ≤ 0.0021: `ANALYSIS/gsea_leave_one_out/loo_sets.tsv` (condition full, seed 1234,
    lead = True)
  - 80 genes, bottom share 0.0657, 19,236 ranked: `SLIDES/figures_cns/chart_running_sum.json`
  - best up q 0.608: `ANALYSIS/gsea_leave_one_out/SUMMARY.md` (shape of the screen); `facts.json` gsea_up_best_q 0.61
  - about 8 % (mean ashr log2FC −0.125): `SLIDES/figures_cns/contamination_check.json` (manuscript_s2)
  - hypoxia nominal p ≤ 0.049 in 41/41, q ≥ 0.295: `loo_sets.tsv` (BUFFA_HYPOXIA_METAGENE, SEMENZA_HIF1_TARGETS)

### 6. Leave any tumor out: the direction holds, the FDR does not (new), 32 s
- We re-ran the published analysis 41 times, leaving out each tumor in turn, at five permutation seeds.
- The direction never moves: NES −1.88 to −2.02, nominal p < 0.001 in every run.
- FDR q runs from 0.001 to 1.0 (below 0.05 in 19 of 41); with all six, the seed alone moves it from 0.02 to 0.31.
- A consistent direction, not a settled false-discovery rate.
- NOTES:
  - The published run is reproduced exactly.
  - Without IL68B, q stays at 0.001–0.032 at every seed, and the full pipeline keeps five of the six translation sets at
    q < 0.05. The 12 sets at q < 0.05 there include stromal and matrix sets, so never say "12 translation sets".
  - Without IL67B, IL69B or IL66B, q is above 0.05 at every seed.
  - Without IL66B (post hoc: the recurrent whose removal separates the arms best), no set reaches q < 0.25 at four of
    five seeds.
  - IL68B carries the size of the gene-level fall (−0.40 → −0.10), not its rank.
- Figure: NEW `chart_loo_q.png` (add to `SLIDES/08_robustness_figures.py`).
  - Data: `ANALYSIS/gsea_leave_one_out/loo_sets.tsv`, set == KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION.
  - x: seven conditions (full, drop_IL67B, drop_IL68B, drop_IL69B, drop_IL66B, drop_NL70B, drop_NL71B), labelled "all
    six", "−IL67B" and so on.
  - y: fdr_q on a log axis from 0.001 to 1.
  - Marks: one filled dot per seed (1234, 1–4), jittered; the six `*_nofilter` runs as open diamonds. Navy when a
    primary is dropped, red when a recurrent is dropped.
  - Dashed lines at 0.05 and 0.25; type 20–28 pt.
  - Caption: "Translation-initiation FDR q (log scale) in 41 re-runs of the published enrichment: each tumor left out,
    one dot per permutation seed; diamonds, gene filter not re-applied. Navy, a primary left out; red, a recurrent.
    Dashed lines at q = 0.05 and 0.25. Nominal p < 0.001 in every run."
- Numbers:
  - 41 runs, NES −2.020 to −1.883, nominal p 0 in all, q 0.001–1.000, 19/41 below 0.05: `loo_sets.tsv` (computed over
    every row of the set)
  - all six, five seeds 0.022–0.312; drop_IL68B 0.001–0.032; drop_IL67B 0.080–0.590; drop_IL69B 0.106–1.000;
    drop_IL66B 0.129–1.000: `ANALYSIS/gsea_leave_one_out/SUMMARY.md` (table 1)
  - exact reproduction: `ANALYSIS/gsea_leave_one_out/reproduction.json`
  - drop_IL66B: no set at q < 0.25 at seeds 1234, 1, 2, 3; one (translation initiation, q 0.129) at seed 4:
    `ANALYSIS/gsea_leave_one_out/loo_screen.tsv`
  - five of six translation sets at q < 0.05 without IL68B (starvation 0.065); 12 sets overall:
    `ANALYSIS/holdout_IL68B/COMPARISON.md` §2
  - −0.40 → −0.10: `SLIDES/figures_cns/contamination_magnitude.json` (leave_one_out.sets.initiation)
  - IL66B separates the arms best: `ANALYSIS/holdout_separation/loo_separation_human.tsv`

### 7. Recurrence after LITT, or less tumor in the sample? (new), 38 s
- Recurrent samples held fewer human reads: 29–43 %, against 43–64 % before ablation.
- Not rat contamination: leaking rat reads could shift translation genes −0.01 at most; they fall −0.40.
- But each tumor's translation score tracks its human share (r = 0.75); adjusted for it, the per-tumor difference
  halves.
- Six animals cannot separate recurrence after LITT from less tumor in the sample.
- NOTES:
  - Arm and share correlate at r −0.70.
  - The per-tumor effect goes from −0.85 (p 0.09) to −0.49 (p 0.43). It was never significant, so say "halves", not
    "loses significance".
  - Gene by gene, the adjusted fall grows (−0.40 → −0.57), because share and arm are collinear and the size factors
    track share (r 0.92). The two adjustments disagree; if asked, give both.
  - Point at IL67B and NL70B: the same share (42.5 and 42.7 %) and the same translation score (−0.43 and −0.38). That
    is one pair, an illustration and not a test, and it points against the headline.
  - A lower share may mean a smaller tumor, more infiltration, or dissection. Tumor volume was not measured.
  - Genes at FDR < 0.05 with no fold threshold go from 102 to 42. They are not comparable to the 35, so keep that off
    the slide.
  - The contamination chart is backup B1.
- Figure: EDIT `chart_graft_scatter.png` into `chart_graft_translation.png`, the translation panel alone (08_robustness_
  figures.py). Fix the IL66B label, which is clipped under the x axis, and pull IL67B out from behind NL70B. Ring that
  pair, draw no fit line, and use 20–28 pt type. Keep navy circles for primaries and red squares for recurrents.
  Caption: "Each point is one tumor: its share of reads assigned to human (x) against its translation-initiation score
  (GSVA, y). Navy, before ablation; red, after regrowth. IL67B and NL70B carry the same share and score alike. Six
  tumors, r = 0.75."
- Numbers:
  - 29.4, 42.7, 33.1 % against 42.5, 45.4, 64.4 %: `ANALYSIS/graft_relation/results/per_sample.tsv` (graft_pct)
  - −0.0093 bound against −0.3996 observed: `SLIDES/figures_cns/contamination_magnitude.json` (sets.initiation)
  - r 0.755 (p 0.083; Spearman 0.94): `ANALYSIS/graft_relation/results/score_correlations.tsv` (translation
    initiation, r_graft)
  - per-tumor −0.845 (p 0.086) → −0.488 (p 0.426): `score_correlations.tsv` (group_effect, group_p,
    group_effect_adj_graft, group_p_adj)
  - arm against share r −0.704: `SLIDES/figures_cns/robustness.json` (graft.r_group_graft)
  - −0.400 → −0.569 gene level: `ANALYSIS/graft_relation/results/translation_sets_lfc_by_model.tsv`
  - size factor against share r 0.92: `score_correlations.tsv` (size_factor)
  - IL67B −0.43 and NL70B −0.38 at 42.52 and 42.67 %: `per_sample.tsv`
  - 102 → 42 (notes only): `ANALYSIS/graft_relation/results/de_models.tsv`

### 8. The astrocyte-like program falls, and survives adjustment (v3 16 + new), 25 s
- Of ten published glioblastoma state signatures, only astrocyte-like falls: +0.26 → −0.27 (p = 0.007; q = 0.075).⁷˒⁸
- It holds without IL68B (p 0.044), nearly without IL66B (p 0.066), and adjusted for human share (p 0.043).
- Deconvolution fits this cell line weakly (R ≈ 0.2) and cannot resolve astrocyte-like cells.
- Footer: 7. Neftel C, et al. Cell 2019;178:835–849. 8. Garofano L, et al. Nat Cancer 2021;2:141–156.
- NOTES:
  - At the matched pair the astrocyte-like score does differ (IL67B +0.26, NL70B −0.34), unlike translation. That is
    one pair.
  - The adjusted p comes from a recomputed GSVA (unadjusted p 0.0046), not from the published scores.
  - CIBERSORT v1.04 with the Neftel reference puts mesenchymal-like lower, 0.51 → 0.39 (p 0.14, not significant), and
    its ranking agrees with GSVA. The reference confuses astrocyte-like with mesenchymal.
  - The LM22 immune reference is a clean negative on the human reads.
  - The mitochondrial signature FALLS (−0.58, p 0.14). The printed abstract says "toward mitochondrial"; do not repeat
    it.
  - Do not quote the GLASS stem-like fractions.
- Figure: EDIT `figures_cns/chart_subtypes.png` (existing slope chart): grey out Garofano mitochondrial so that
  astrocyte-like is the only highlight. Caption: "Mean GSVA score per arm for ten published glioblastoma state
  signatures (Neftel 2019, Garofano 2021). Only astrocyte-like falls at nominal p < 0.05; none clears FDR across the
  ten (lowest q = 0.075)."
- Numbers:
  - +0.261 → −0.269, p 0.0075, q 0.075: `SLIDES/figures/facts.json` gsva_by_name.Neftel_AC
  - −0.57, p 0.044 without IL68B: `ANALYSIS/holdout_IL68B/COMPARISON.md` §5
  - −0.50, p 0.066 without IL66B: `ANALYSIS/holdout_IL66B/COMPARISON.md` §5
  - adjusted p 0.0431 (unadjusted 0.0046): `ANALYSIS/graft_relation/results/score_correlations.tsv` (Neftel_AC)
  - fit R 0.189–0.258 in the six tumors: `ANALYSIS/cibersort/results/v104/fractions_v104_neftel4_confident.tsv`
    (Correlation); `SLIDES/figures_cns/robustness.json` deconv
  - MES 0.514 → 0.390: `robustness.json` deconv; p 0.137: `score_correlations.tsv` cs_MES
  - astrocyte-like over-called by about 9 points: `ANALYSIS/cibersort/SUMMARY.md`
  - IL67B +0.26, NL70B −0.34: `ANALYSIS/graft_relation/results/per_sample.tsv` (Neftel_AC)
  - LM22 correlation −0.034 to −0.048, p 0.92–0.99 in the tumors:
    `ANALYSIS/cibersort/results/v104/fractions_v104_lm22_newman2015.tsv`
  - mitochondrial −0.58, p 0.14: `SLIDES/figures_cns/chart_subtypes.json`

### 9. The signature points to one approved drug: ciclopirox (v3 18 + 19 + 20), 30 s
- Drug gene sets are ranked by how strongly they oppose the recurrence, then weighted by predicted brain
  penetration.⁹⁻¹¹
- 92 compounds → 54 in clinical use → 15 that both barrier models say cross. Ciclopirox ranks first (3.14; next 2.34),
  with or without the weight.
- Re-run without one tumor: still first without IL68B, third without IL66B.
- It is approved as a topical antifungal; an oral form completed phase 1¹²; another group's glioblastoma screen also
  nominated it.¹³
- Footer:
  - 9. Yoo M, et al. Bioinformatics 2015;31:3069–3071.
  - 10. Swanson K, et al. Bioinformatics 2024;40:btae416.
  - 11. Daina A, Zoete V. ChemMedChem 2016;11:1117–1121.
  - 12. Minden MD, et al. Am J Hematol 2014;89:363–368.
  - 13. Sun S, et al. J Transl Med 2025;23:25.
  - "Computational predictions; none tested in this model."
- NOTES:
  - The filter keeps oxygen, ozone and magnesium, and drops deferoxamine, which has no recorded phase. The two
    strongest reversals overall, dimethyloxalylglycine and deferoxamine (|NES| 2.39), have no phase.
  - The rebuilt barrier model clears 13, not 15. Ciclopirox is first in both.
  - Without IL66B the unweighted rank is 14, so "with or without the weight" holds for all six tumors only.
  - The per-tumor GSVA score of the ciclopirox gene sets does not differ by arm (p 0.84). The rank rests on the
    ranked-list enrichment.
  - The abstract's score of 3.40 (barrier 1.0) came from an earlier ranking.
- Figure: `figures_cns/chart_drug_scatter.png` (existing dot chart). Caption: "Each point is one of the 54 clinically
  available compounds: how strongly its gene set opposes recurrence (|NES|, x) against predicted barrier permeability
  (ADMET-AI, y). Filled, both barrier models predict crossing; size, combined score."
- Numbers:
  - 100 signatures, 92 compounds, 54 clinical, 15 both: `SLIDES/figures_cns/chart_funnel.json`
  - ciclopirox 3.14, pyrantel 2.34: `SLIDES/figures/facts.json` drug_tierA
  - rank 1/1 (all six), 1/1 (−IL68B, score 3.02), control 13 both: `ANALYSIS/holdout_IL68B/COMPARISON.md` §4
  - rank 3/14 (−IL66B, score 2.02): `ANALYSIS/holdout_IL66B/COMPARISON.md` §4
  - DMOG and deferoxamine |NES| 2.388 and 2.387 with no phase: `SLIDES/figures_cns/chart_drug_scatter.json`
    (strongest_overall)
  - per-tumor ciclopirox GSVA +0.05, p 0.84: `ANALYSIS/graft_relation/results/score_correlations.tsv`
    (DSigDB_ciclopirox)
  - abstract 3.40 and 1.0: `ABSTRACT.md`
  - references 12 and 13: `SLIDES/deck_v3_text_dump.txt` slides 20 and 23

### 10. Why ciclopirox: it hits a translation step U-251 needs (v3 22 + one line of 21), 28 s
- Ciclopirox blocks DOHH, the enzyme that activates eIF5A, which ribosomes need to keep elongating.
- In DepMap, U-251 depends on DOHH (−0.58), a step most lines can spare (median −0.39), and on DHPS and EIF5A.¹⁴
- It slows U251 in culture and as a subcutaneous xenograft.¹⁵
- LITT opens the barrier around the ablation for weeks:² a window to deliver it.
- Footer: 14. DepMap 24Q4, Chronos gene effect, U-251 MG. 15. Su Z, et al. Cell Death Dis 2021;12:251. "Public
  cell-line data and a prediction; nothing has been dosed in this model."
- NOTES:
  - Most lines need DHPS and EIF5A (medians −0.83 and −0.97), so they are not a selective weakness. The margin on DOHH
    is modest.
  - Ribonucleotide reductase is common-essential (RRM1 −3.14).
  - If asked why a translation blocker would help a tumor already translating less: the genes the drug induces are
    lower in recurrence, so it opposes the recurrent state. A tumor running lean may have less to spare, but that is a
    hypothesis.
- Figure: NEW `chart_depmap_dots.png` (redraw of `chart_depmap.png` in `04_make_cns_figures.py`).
  - Genes: DOHH, DHPS and EIF5A, optionally with RRM1 and RRM2.
  - Marks: the U-251 value as a filled dot, the 1,178-line median as a tick, and a dashed line at −0.5. Type 20–28 pt.
  - Fallback: crop the existing chart to the hypusination group.
  - Caption: "DepMap CRISPR gene effect in U-251 MG (dot) and the median of 1,178 cell lines (tick). Below −0.5
    (dashed), the line depends on the gene. Public screen of cells in culture."
- Numbers:
  - DOHH −0.584 (median −0.394), DHPS −1.091 (−0.834), EIF5A −1.247 (−0.969), RRM1 −3.145:
    `depmap/depmap_u251_targets.csv`
  - 1,178 lines: `SLIDES/04_make_cns_figures.py` (lines 737 and 753, the chart legend)
  - Su 2021 and Cleary 2026: `SLIDES/deck_v3_text_dump.txt` slides 21 and 22

### 11. What this means (v3 23, rewritten), 28 s
- After LITT, the regrown tumor in this model ran less protein-synthesis machinery and lost its astrocyte-like
  signature.
- The direction held whichever tumor was left out, and rat reads cannot make it. Its FDR and its size lean on single
  tumors.
- It is entangled with how much tumor each sample held, which three per arm cannot separate.
- It names a drug, ciclopirox, and a window, the weeks after ablation, to test.
- Footer: "Investigational: one cell line (U251N), one rat xenograft model, n = 3 per arm; drug candidates are
  computational predictions; GEO GSE338105."
- NOTES: This is the one sentence of the talk. Say it, pause, and move on.
- Figure: none (text slide).
- Numbers: none new; each is sourced on slides 5–10.

### 12. What six animals cannot tell us, and the experiment that can (v3 24, rewritten), 28 s
Two columns, limitation → experiment, in black text.
- Three per arm, and tumor content differs by arm → more animals matched on tumor content, with MRI tumor volume as a
  covariate (power analysis under way).
- Translation inferred from RNA → puromycin incorporation or polysome profiling of the regrown tumor.
- One cell line in an athymic rat → an immunocompetent model, peri-ablation patient tissue, and public recurrence
  cohorts (under way).
- A predicted drug → ciclopirox dosed in the post-ablation barrier window, with survival as the endpoint.
- NOTES:
  - Methylation (EPIC, 866,238 probes, same six tumors): no probe passes FDR (smallest q 0.092), and purity and arm are
    confounded.
  - Before the talk, ask the Henry Ford group whether MRI tumor volumes exist for these six animals.
  - Quote no number from the pending analyses: host (rat) response, stem-cell drug screen, cohorts, power.
- Figure: none.
- Numbers:
  - methylation (notes only): `SLIDES/deck_v3_text_dump.txt` slide 24 notes (`REVIEW/purity_confound.py`)
  - cohorts downloaded: `ANALYSIS/human_cohorts/raw`

### 13. Acknowledgments and data availability (keep v3 25), 10 s
- As in v3: Nagaraja, Datta and Lee (Hermelin Brain Tumor Center, Henry Ford Health); Award A20050; IACUC #1509,
  ARRIVE; GEO GSE338105; Questions: go2432@wayne.edu.
- NOTES: Leave this slide up during questions. Backups B1–B5 follow.

### Backups (after Questions, not presented)

**B1. Rat reads cannot produce the fall** (v3 15 chart)
- The bound is −0.009 on each set. Observed: −0.400, −0.385 and −0.395 (initiation, elongation, ribosome).
  Abundance-matched genes: −0.13, −0.12 and −0.13.
- Without IL68B the sets fall −0.11, −0.10 and −0.10 (matched −0.03).
- RUVSeq with k = 2 gives −0.28 to −0.30, but those factors also separate the arms.
- Figure: `chart_contamination_magnitude.png` (existing vertical columns).
- Sources: `SLIDES/figures_cns/contamination_magnitude.json`; `contamination_check.json` (k2).

**B2. The whole pipeline without IL68B or IL66B**
- Two-fold DE genes: 35, 43 and 140. Translation-initiation q at seed 1234: 0.022, 0.027 and 0.309. Sets at q < 0.05:
  1, 12 and 0.
- Ciclopirox rank (weighted / unweighted): 1/1, 1/1 and 3/14. Astrocyte-like: −0.53 (p 0.007), −0.57 (0.044) and
  −0.50 (0.066).
- IL66B was chosen post hoc.
- Figure: `figures_cns/chart_holdout.png` (existing). EDIT: move the "1" label off the line in panel 2, and title
  panel 1 "FDR q, seed 1234".
- Sources: `ANALYSIS/holdout_IL68B/COMPARISON.md` and `ANALYSIS/holdout_IL66B/COMPARISON.md` §1, 2, 4, 5;
  `robustness.json` holdout.

**B3. Hypoxia down, cell cycle up, nominally** (v3 12 + 14)
- HIF1 targets −1.79 (p 0.016) and the Buffa metagene −1.79 (p 0.011): nominal p ≤ 0.049 in all 41 re-runs, q ≥ 0.29.
- Mitotic spindle +1.76 and G2M +1.71; the spindle p reaches 0.054 in some re-runs.
- Glycolysis −1.57 (p 0.055), up to 0.13 when a tumor is left out.
- None of these was tested against human share.
- Figure: text only. Optionally, the dot-and-range chart from design A (NES with all six tumors and its range across
  the six drops).
- Source: `ANALYSIS/gsea_leave_one_out/loo_sets.tsv`.

**B4. Cell-state deconvolution (CIBERSORT v1.04)**
- Neftel four-state reference: fit R 0.19–0.26 in every tumor; mesenchymal-like 0.51 → 0.39 (p 0.14); astrocyte-like
  not recoverable.
- LM22 immune reference on the human reads: correlation below zero in every sample.
- Figure: `figures_cns/chart_states.png` (existing dot chart).
- Sources: `results/v104/fractions_v104_neftel4_confident.tsv`, `fractions_v104_lm22_newman2015.tsv`,
  `robustness.json` deconv, `score_correlations.tsv` cs_MES.

**B5. What is already known about the top twenty** (v3 21)
- Figure: `table_prior_art.png` (existing).
- Source: `SLIDES/figures_cns/table_prior_art.json`.
- This audit covers the published top 20 only. The compounds that enter the hold-out top-20 lists were not audited.

## 5. Figures to build

**New:**

| figure | slide | drawn in | from |
|---|---|---|---|
| `chart_running_sum_ti.png` | 5 | `05_make_more_figures.py` | `chart_running_sum.json` |
| `chart_loo_q.png` | 6 | `08_robustness_figures.py` | `loo_sets.tsv` |
| `chart_graft_translation.png` | 7 | `08_robustness_figures.py` | edit of `chart_graft_scatter.png` |
| `chart_depmap_dots.png` | 10 | `04_make_cns_figures.py` | `depmap/depmap_u251_targets.csv` |

**Edits:**
- `chart_subtypes.png`: grey out the mitochondrial line.
- `chart_holdout.png`: fix the panel 2 label.

**Kept as they are:** `fig_margin`, `overview_pipeline`, the video, `chart_pca`, `chart_drug_scatter`,
`chart_contamination_magnitude`, `chart_states` and `table_prior_art`.

`08_robustness_figures.py` must also write the leave-one-out summary into `robustness.json`: 41 runs, the NES range,
the q range and 19/41. That way `04_build_cns_template_deck.py` types no number. The script is untracked; commit it
with the rebuild.

## 6. Risks

1. **The printed abstract contradicts the talk.** `ABSTRACT.md` gives:
   - elongation NES −3.17 (FDR 1.6e-26);
   - ciclopirox at 3.40 with barrier 1.0;
   - a shift "toward Garofano Mitochondrial";
   - DMOG and LY-294002 as candidates.

   In the released pipeline:
   - the mitochondrial score falls (−0.58, p 0.14);
   - DMOG has no clinical phase;
   - ciclopirox scores 3.14 with barrier 0.92.

   No file here shows where −3.17 came from. Have one sentence ready: the abstract was written from an earlier
   analysis, and the talk reports the released pipeline, reproduced exactly.
2. **Tumor content (slide 7) is the weakest point and the likeliest question.**
   - The two adjustments disagree. Gene by gene the fall grows (−0.40 → −0.57); per tumor it halves (p 0.086 → 0.43).
   - Share, arm and size factors are collinear (r −0.70 and 0.92).
   - The IL67B/NL70B pair argues that translation may follow tumor content, while the astrocyte-like score does not.
     That is n = 1, and the talk says so rather than hiding it.
   - The quickest real test may already exist: the Henry Ford group's MRI tumor volumes for these six animals.
3. **Both hold-outs were chosen post hoc.**
   - IL68B was picked after it was seen to dominate the translation genes. IL66B was picked because its removal
     separates the arms best.
   - Only these two had full-pipeline re-runs, so "ciclopirox first or third" covers two re-runs.
   - Without IL66B the unweighted rank is 14.
4. **Statements in the checkpoint or task text that the files do not support.**
   - "IL66B: no set at q < 0.25 at any seed" (`RESUME.md`): seed 4 has translation initiation at q 0.129
     (`loo_screen.tsv`).
   - "r 0.76–0.84": translation initiation is 0.755, and only MES1 reaches p < 0.05.
   - "12 translation sets without IL68B": five of six translation sets clear q < 0.05. The 12 include stromal and matrix
     sets.
   - "The translation difference loses significance": it was never significant (p 0.086).
5. **Two CIBERSORT files are stale.**
   - `ANALYSIS/cibersort/SUMMARY.md` still describes the re-implementation.
   - `results/crosscheck.txt` also reads `fractions_nusvr_*` (written 09:14, before v1.04 at 09:48). Its numbers
     differ only in the third decimal.
   - Quote `results/v104/*.tsv`, `robustness.json` deconv and `score_correlations.tsv` cs_* instead.
   - The Varn GLASS call stays off every slide: the rat-brain controls read 0.88–0.99 "proliferating stem-like" too.
6. **`ANALYSIS/graft_relation/SUMMARY.md` has not landed.** Its adversarial review, wf_47e66eea-c0c, is still running.
   Every graft number here was read from the result TSVs; re-check slides 7 and 8 when it lands. The adjusted
   astrocyte-like p comes from a recomputed GSVA (unadjusted −0.66, p 0.0046), not the published scores. The hold-out
   p values are Welch tests on 2 vs 3.
7. **The ciclopirox per-tumor score does not separate the arms** (DSigDB union GSVA +0.05, p 0.84). The rank rests on
   the ranked-list enrichment, so never say "each recurrent tumor" about the drug.
8. **Drift and interpretation.**
   - The rebuilt ADMET-AI clears 13 compounds, not 15, and the September version was never recorded.
   - Why use a translation blocker on a tumor already translating less? The slide 10 notes carry the answer.
9. **Pending analyses** (host DE and states, the stem-cell screen, the CGGA/GLASS/TCGA cohorts, power) are not in hand.
   - Slide 12 names them as under way, with no numbers.
   - If cohort replication lands before 2 November, it is the one result worth a slide: after slide 8, paid for from
     the 39 s of slack.
   - Do not change the headline for any of them until they have been through the same leave-one-out and hold-out
     checks.
10. **Open items carried over from v3.**
    - Reconcile "four per arm" (J Neurosurg) with 4 + 3 (Methods).
    - Greg has not yet confirmed reading the disclosures aloud.
    - Slide 2 uses a drawn schematic, where the standing preference is a published figure (permission pending).
11. **Timing and density.**
    - The plan totals 321 s of 360.
    - Slide 7 is the densest (four lines, 38 s); rehearse it.
    - If a four-line slide overflows at 18 pt, move the parenthetical p values to the notes rather than shrink the type.

## 7. As built (2026-09-26), and where it departs from this plan

Greg, after the plan: "make sure these plots are styled like typical plots in this field of literature ... essential
plots that any study like this in nature or so would have, not just inventing random ways to portray data". So the
figure forms changed; the story and running order did not.
- Slide 5: the Broad GSEA enrichment plot (running score, hit barcode, ranked-list colour bar and metric) beside a
  leading-edge heatmap (77 core-enrichment genes, rlog row z-scores, annotation rows for arm and human share). The
  heatmap shows IL68B highest on all 77, so the IL68B dependence is on the slide, not only in the notes.
- Slide 6: the leave-one-out result as a condition x seed table of FDR q (q < 0.05 shaded) with the NES range.
- Slide 7: score against human share, OLS fit with 95 % band, r and P, both panels (translation, astrocyte-like).
- Slide 8: a GSVA heatmap of the ten signatures with P and BH q, on DESeq2 variance-stabilised counts (the review's
  fix for the TPM composition artefact): AC P 0.002, q 0.021, adjusted for human share 0.024. The manuscript's TPM run
  (p 0.007, q 0.075) is in the notes.
- Slide 9: the drug ranking as a plain table (top six of the 54 with a clinical phase) instead of the scatter.
- Slide 10: DepMap 24Q4 distributions over all 1,178 lines with U-251 MG marked, instead of bars with a median tick.
- Backups: B1 contamination bars (kept), B2 the whole-pipeline hold-outs as a table, B3 CIBERSORT stacked fractions,
  B4 the prior-art table. The hypoxia / cell-cycle backup was dropped; the drawn slope, dot and line charts are unused.

