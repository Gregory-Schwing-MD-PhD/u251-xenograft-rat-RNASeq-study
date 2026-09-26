# CNS Abstract 418: printed-abstract numbers against the talk, 2026-09-26

Slide numbers below are from the deck before the design slide was inserted as slide 3 (add 1 to slides 3 onward).

**CNS Abstract 418: printed abstract compared with the v4 talk**

The printed abstract (the CNS program or Neurosurgery supplement text) is not on disk, so I could not compare against it directly. What exists is two drafts from 26 May 2026 with the same numbers. Every number in both drafts comes from one earlier report, produced by the first-pass figure script on 25 May.

## 1. What the printed abstract is known to be

- `C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study/ABSTRACT.md`. Last commit f3936cc, 26 May 15:02 CDT. `STORY_v4.md` and `RESUME.md` call this "the printed abstract", but nothing on disk confirms that it is the text that was submitted.
- `C:/Users/grego/OneDrive/Desktop/CNS_SNO Abstract_tnn edits.docx`. Author Nagaraja, last saved by Greg on 26 May at 14:39 CDT, which is 23 minutes before the last ABSTRACT.md commit. It has the same numbers plus four extra phrases:
  - "phase 4"
  - "LY-294002 ranked 7th (2.40), consistent with the translational/mTOR axis"
  - "residual host contamination was tested and found negligible"
  - "robust to contamination (which would bias ribosomal genes upward, opposite the finding)"
- The two CNS portal PDFs on the Desktop are other abstracts, not Abstract 418:
  - `MANUSCRIPT/2026 CNS Abstract Submission Portal….pdf` is the DMG abstract.
  - `NextPaper/old/…pdf` is the OSC abstract.
- Nothing matching was found in Downloads or Documents.
- **Where the numbers came from:** `publication_figure/therapy_impact_LLM_Analysis_Report.txt`, generated 2026-05-25 (commit b3972c6) by `create_publication_figure_600_dpi.R`. Every abstract number appears in that report.
- **Action for Greg:** get the Abstract 418 text from the CNS portal ("Abstract Submission Detail") or the program. Until then, the one line to use is: "the abstract was written from our first-pass figure report of 25 May; the talk reports the released pipeline."

## 2. Claim by claim

Paths are relative to the repository root. "Agree" means the same value.

| # | Abstract | Current (file) | Agree? | Why they differ |
|---|---|---|---|---|
| 1 | n = 3 per arm; primary "(pre-LITT)" | 3 v 3 (slides 3/4; `ANALYSIS/graft_relation/results/per_sample.tsv`) | Yes on n | "Pre-LITT" means animals that were never ablated. The design is unpaired (slide 3 notes). |
| 2 | xengsort k=25; DESeq2 Wald test, ashr, FDR<0.05, abs(log2FC)>1 | Same (slide 4 footer) | Yes | — |
| 3 | Methods: "subtypes by GSVA with limma" | Slide 8: GSVA on DESeq2 vst counts with the published Neftel/Garofano signatures, t-test, BH correction (`SLIDES/figures_cns/standard.json` subtypes_input) | **No** | The May subtype numbers were not GSVA. They were mean row z-scores of 5–6 hand-picked marker genes, tested with arrayWeights limma (R script panel C at b3972c6; the report says "mean signature z"). |
| 4 | PC1 41.8 %, PC2 26.2 % | 35 % / 25 % (slide 4 axes; `standard.json` pca.var_pc1 0.350, var_pc2 0.247) | **No** | May: prcomp on the 500 most variable genes of the vst matrix. Talk: the same selection on rlog. Different transform of the same counts. |
| 5 | "Globally distinct by PCA" | "The arms separate modestly"; clustering does not recover the arms (backup B5; `standard.json` sample_dist.k2_recovers_arms = false) | Framing | The abstract leaves out the PERMANOVA p. |
| 6 | PERMANOVA R² = 0.278 | R² 0.28, p = 0.10, the floor for 3 v 3 (slide 4; `SLIDES/figures/facts.json` permanova) | Yes | The abstract omits p = 0.10. |
| 7 | 35 DE genes (23 up, 12 down) | 35 (23/12) (`facts.json` de_up/de_down; slide 4, B2) | Yes | The talk adds 43 or 140 genes when one tumour is left out (B2). |
| 8 | Elongation NES −3.17, FDR 1.6×10⁻²⁶ | NES −1.96, nominal p < 0.001, q 0.158 (`facts.json` gsea_down_top; `standard.json` gsea.six_q) | **No** | See the GSEA note below the table. |
| 9 | Ribosome, translation initiation, GCN2 named (report: −3.15, −3.13, −3.11; FDR ~1e-25) | Ribosome −1.94 (q 0.204); KEGG MEDICUS translation initiation −1.99 (q 0.022, the only set below 0.05, 0.022–0.31 across seeds); GCN2 −1.95 (q 0.146) (slides 5–6) | Direction yes, significance no | Same method difference as row 8. The lead set changes from elongation to initiation. |
| 10 | "Translational shutdown with integrated stress" | Talk: lower ribosomal-protein (RP) mRNA, "not yet lower translation" (slides 5, 12) | **No** | Of the genes in each "stress" set, 85 of 99 (GCN2) and 85 of 150 (starvation) are RP genes (`standard.json` gsea.rp_in_other_leading_sets). The six sets' leading edges share 76 genes; their union is 86, of which 78 are RP genes. The translation-initiation set is 80 of 80 RP genes. |
| 11 | Shift toward Garofano Mitochondrial | Falls: 0.26 → −0.32, P 0.13, q 0.43 (`standard.json` subtypes.Garofano_MTC, slide 8). TPM run −0.58, p 0.14 (`facts.json` gsva_by_name; `chart_subtypes.json`). Without IL66B −0.73, p 0.127 (`ANALYSIS/holdout_IL66B/COMPARISON.md` §5) | **No, opposite sign** | The May "mitochondrial" score was 6 TCA-cycle genes (CS, ACO2, IDH2, IDH3A, OGDH, SDHA) as a mean z-score, not the published signature. OXPHOS cannot break the tie: Hallmark −1.45, KEGG +1.20, WP +1.17, GO +1.16, none with q < 0.6 (`ANALYSIS/results_therapy_v3/report_gsea/`). |
| 12 | Away from Neftel AC-like | +0.28 → −0.29, P 0.002, q 0.021; adjusted for human share P 0.024 (slide 8; `standard.json`). TPM run p 0.0075 | Yes | — |
| 13 | Away from NPC-like | NPC1 +0.03 (P 0.90); NPC2 −0.15 (P 0.40) (`standard.json`) | **No** | The May score was a 5-gene list (DCX, DLL3, ASCL1, NEUROG2, STMN2). |
| 14 | Hubs BGN, COL1A1, IGFBP3, CALB1 | Network reproduces (`ANALYSIS/holdout_IL68B/COMPARISON.md` §6). Only 11 of the 35 genes connect. Degree: BGN 5, COL1A1 4, IGFBP3 4, CALB1 1 (`SLIDES/figures_cns/chart_ppi.json`). Not on any v4 slide | Partly | The script labels the top 15 by degree as hubs, so all 11 connected genes became "hubs". CALB1 has one interaction. |
| 15 | Ciclopirox Integrated Score 3.40, BBB 1.0 (docx: phase 4) | 3.27, ADMET-AI BBB 0.96, BOILED-Egg "in". Rank 1 of 62 with the barrier weight, 2 without it (deferoxamine 3.69) (`ANALYSIS/drug_rematch/results/published_drug_ranking_final.csv`, `rematch.json`; slide 9). NES −2.261, FDR 1.66e-5 and phase 4 unchanged | Score no; rank and phase yes | 3.40 = 2.261^1.5. The May barrier rule (points for MW, logP, PSA, HBD, HBA, divided by 4, capped at 1.0) gave ciclopirox the maximum 1.0, so the abstract's 3.40 equals today's unweighted score. The earlier v4 figure of 3.14 (BBB 0.92, `facts.json` drug_top) came from an earlier ADMET-AI build. |
| 16 | DMOG, 2nd, 3.04 | Not in the talk. CHEMBL92309 has no clinical phase (audited, `ANALYSIS/drug_rematch/overrides.tsv`), so it is outside the 62. Its abs(NES) of 2.388 is the largest of all 92 compounds, tied with deferoxamine at 2.387 (`chart_drug_scatter.json`) | **No** | The May ranking had no clinical-phase filter. It used barrier score 0.825, and 2.388^1.5 × 0.825 = 3.04. |
| 17 | LY-294002 ranked 7th (docx: 2.40) | Excluded. CHEMBL98350 has max_phase −1 (`ANALYSIS/drug_rematch/chembl_mapping.csv`). Its abs(NES) of 2.039 is 6th of 92 | **No** | Same filter as row 16. "7th" counted rows; ciclopirox filled rows 1 and 2, so LY-294002 was the 6th compound. |
| 18 | Conclusions: "mitochondrial-metabolic shift"; candidates ciclopirox, DMOG, LY-294002 | Astrocyte-like fall only; ciclopirox is the only clinically available candidate; deferoxamine noted (slide 9 notes) | **No** | Rows 10, 11, 16, 17. |
| 19 | Docx only: contamination negligible; "would bias ribosomal genes upward" | About 0.01 log2 against a 0.40 fall; worst case about 0.05. The sign depends on one control: −0.01 with all three, +0.001 without N269B (slide 7 notes, B1) | Size yes, direction no | The abstract never raises the confound between arm and human share (r −0.70, slide 7). |

**Note on the GSEA numbers (rows 8–9):**
- **Where −3.17 comes from:** report line 28. It is clusterProfiler GSEA (fgsea), with genes ranked by the DESeq2 Wald statistic, gene permutation, and `eps=1e-50`.
- **What the pipeline specifies:** Broad GSEA with gene-set permutation, Diff_of_Classes, 1,000 permutations, seed 1234 (`therapy_v3_params.yaml`). With 1,000 permutations nominal p cannot go below 1e-3, so a q of 1.6e-26 is an extrapolated tail estimate, not a permutation result.
- **The script itself switched:** on 2026-07-18 (commit 392157d) Figure panel D moved to the Broad numbers. The script's comment says the fgsea estimate "is not appropriate at n=3/arm".

**A problem this exposes in the deck itself:** slide 9 still uses that same clusterProfiler/fgsea test. Its drug NES of −2.26 and the text "(36 at FDR < 0.05)" both come from it (36 matches `publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv`). Present those values as a ranking, and drop or qualify the "36 at FDR < 0.05". Otherwise the talk rejects the fgsea FDRs on slide 5 and relies on them on slide 9.

## 3(a) Sentences to say on stage

Slides 5, 8, 9 and 11 are the priority. With about 27 s per slide there may only be room for the slide 11 one; the others can be held for questions.

- **Slide 5, "After LITT, the regrown tumor turns its ribosomal-protein genes down":** "The printed abstract gives an NES of minus 3.17 and an FDR near ten to the minus 26 for this block. That came from a first-pass test that treats every gene as independent. The permutation test the pipeline specifies gives the same direction, an NES near minus 2, and one set at q = 0.022."
- **Slide 8, "Of ten published subtype signatures, only the astrocyte-like score separates the arms":** "The abstract says the tumors shifted toward the mitochondrial subtype. That came from a six-gene TCA-cycle score in our first pass. With the published signatures the mitochondrial score falls, not significantly, and the change that holds is the astrocyte-like fall the abstract also reported."
- **Slide 9, "A candidate drug against the recurrent state: ciclopirox":** "The abstract's 3.40 for ciclopirox is this same enrichment with a rule-of-thumb barrier score of 1.0. With a trained barrier model it is 3.27 and still first. DMOG and LY-294002, named in the abstract, drop out because neither has a clinical record."
- **Slide 11, "What this means":** "Three things differ from the printed abstract: I would now say lower ribosomal-protein genes rather than translational shutdown, the mitochondrial shift did not survive the published signatures, and of its three drugs only ciclopirox is clinically available."
- **Slide 4, "The arms separate modestly, and few single genes move" (only if asked):** "The abstract's 42 % on PC1 used a different transform of the same counts. The PERMANOVA value, 0.28, is unchanged, and its p of 0.10 is the smallest three against three allows."

## 3(b) Speaker-note lines, ready to paste

- **Slide 3, "A human tumor in a rat brain…":** Abstract "primary (pre-LITT)" means animals that were never ablated, not biopsies taken before ablation. The design is unpaired.
- **Slide 4, "The arms separate modestly, and few single genes move":**
  - Abstract PC1 41.8 % / PC2 26.2 % came from the May figure script (vst, 500 most variable genes). This plot uses rlog: 35 % / 25 %.
  - PERMANOVA R² 0.278 and 35 genes (23 up, 12 down) are unchanged. The abstract omits p = 0.10.
  - Abstract hubs (BGN, COL1A1, IGFBP3, CALB1): the network reproduces, but only 11 of the 35 genes connect. Degree: BGN 5, COL1A1 4, IGFBP3 4, CALB1 1.
- **Slide 5, "After LITT, the regrown tumor turns its ribosomal-protein genes down":**
  - Abstract: elongation NES −3.17, FDR 1.6e-26 (therapy_impact_LLM_Analysis_Report.txt, clusterProfiler/fgsea, Wald-stat ranking, eps 1e-50).
  - The pipeline's Broad GSEA, gene-set permutation, 1,000 permutations: elongation −1.96, q 0.158; initiation −1.99, q 0.022. It cannot reach q below 0.001.
  - Figure panel D was switched to these numbers on 2026-07-18 (392157d).
  - "Integrated stress": 85 of the 99 GCN2-set genes and 85 of the 150 starvation-set genes are RP genes, so these sets are the same block.
- **Slide 6, "Leave any tumor out…":** The abstract's FDR is not a Broad GSEA q. Here q runs 0.022–0.31 with all six tumours, depending on the seed alone.
- **Slide 8, "Of ten published subtype signatures…":**
  - Abstract "toward Garofano Mitochondrial, away from AC/NPC" came from May panel C: mean z-scores of hand-picked 5–6-gene lists (mitochondrial = CS, ACO2, IDH2, IDH3A, OGDH, SDHA), arrayWeights limma. It was not GSVA, although the abstract's Methods say GSVA.
  - Published signatures: mitochondrial 0.26 → −0.32, P 0.13 (TPM run p 0.14). AC-like agrees. NPC1 +0.03 (P 0.90), NPC2 −0.15 (P 0.40).
  - OXPHOS depends on the gene-set collection: Hallmark −1.45, KEGG +1.20, WP +1.17, GO +1.16.
- **Slide 9, "A candidate drug against the recurrent state: ciclopirox":**
  - Abstract 3.40 / BBB 1.0 = 2.261^1.5 × the May rule-of-thumb barrier score, which was capped at 1.0. ADMET-AI 0.96 gives 3.27: rank 1 of 62, 2 without the weight.
  - DMOG (CHEMBL92309, no phase; abs(NES) 2.388, largest of all 92) and LY-294002 (CHEMBL98350, phase −1; it was the 6th compound, since ciclopirox filled rows 1 and 2) are outside the clinical list.
  - The drug NES and "36 at FDR < 0.05" come from the same fgsea test as the abstract's −3.17: quote them as a ranking only.
- **Slide 11, "What this means":** Abstract conclusions:
  - "translational shutdown / integrated stress" → now lower RP mRNA, translation not measured.
  - "mitochondrial-metabolic shift" → not supported.
  - "ciclopirox, DMOG, LY-294002" → ciclopirox only.
- **Slide 12, "Limitations, and the experiments that answer them":** The abstract's "translational shutdown" is what the puromycin/polysome experiment tests. It has not been shown.
- **Slide 14, "Backup B1. Is the fall rat contamination?":** Use this only if the printed text is the docx version. Its claim that contamination "would bias ribosomal genes upward" does not match the bound here: about 0.01 log2, with the sign set by one control (−0.01 with all three controls, +0.001 without N269B). The size argument holds; the direction argument does not.

If the printed abstract includes the nine-panel figure, panels C, D and H carry the same May numbers.

Nothing was edited, committed or run on the grid. The committed `SLIDES/deck_v4_text_dump.txt` (15:36) is out of date. I regenerated the text from the pptx at HEAD (da7800f).

Files are in C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study:
- publication_figure/therapy_impact_LLM_Analysis_Report.txt
- ABSTRACT.md
- create_publication_figure_600_dpi.R
- SLIDES/figures/facts.json
- SLIDES/figures_cns/standard.json
- ANALYSIS/drug_rematch/results/published_drug_ranking_final.csv

Also:
- C:/Users/grego/OneDrive/Desktop/CNS_SNO Abstract_tnn edits.docx
- C:/Users/grego/AppData/Local/Temp/claude/c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1/f1bdbd78-151f-470b-b703-dd9af9b3fecc/scratchpad/deck_v4_now.txt (regenerated deck text)
