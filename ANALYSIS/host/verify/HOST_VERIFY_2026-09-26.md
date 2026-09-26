# Host (rat) verification and deconvolution, 2026-09-26

## Differential expression checks

Both grid jobs are submitted and were still pending in the queue at 18:45 EDT; neither had started, so there are no grid results yet. A local run of the same R code already shows the slide needs changing (preliminary, see "Changes how B7 should read").

**Jobs** (grid clock; ETA assumes each job starts within about 15 minutes):

| job | what | ETA | output |
|---|---|---|---|
| 40475691 | `run_host_verify.sbatch` (reqp, requeue QoS, 4 cpu, 24G, 6h). Step A runs in the pipeline's own DESeq2 1.34.0 + ashr container: R0 reproduce, R1 leave-one-out, R2 the 10 relabellings, R3 NL pair vs the other four and IL66B alone, R5 human counts beside rat. Step B: R4 fgsea on every fit's Wald statistic. | ~20:15 EDT | `/rs/rs_grp_oschome/go2432/u251_host/verify/results/` (`summary.json`, `r0_*`…`r5_*`, `r4_gsea_*`); log `verify/logs/verify_40475691.out` |
| 40475693 | `run_de_therapy.sbatch`: nf-core/differentialabundance 1.5.0, host_therapy only, GSEA on, launched like `run_host_de.sbatch` | ~19:45 EDT | `/rs/rs_grp_oschome/go2432/u251_host/results_de_therapy/`; `verify/de_therapy_compare.txt` compares it with the old table and old GSEA |

- **Step B fallback:** if fgsea is missing from the image, the job installs it into `/rs/rs_grp_oschome/go2432/Rlib/R<x.y>`. If step B fails anyway, it falls back to Broad GSEAPreranked on the full and leave-one-out fits.
- **Job 2 inputs:** the sample sheet has only the six tumours. The counts are the original run's 17,451 genes, cut to those six columns. Filtering is set to 0 so that gene list stays the same. This removes the crash, and the new table and GSEA should be byte-identical to the old run's.

**`loo_separation_host.tsv` (existing grid file)**
```
reads held_out n genes_tested padj05 padj05_fc1.5 pc1_var pc2_var silhouette_pc12 pc1_splits_groups ward_k2_recovers_groups
host  none   6 16856  38  38 57.9 20   0.089 FALSE FALSE
host  IL67B  5 16777  33  33 53.6 24.2 0.045 FALSE FALSE
host  IL68B  5 16756  27  27 59.8 21.1 -0.114 FALSE FALSE
host  IL69B  5 16812  29  29 66.6 20.6 0.122 FALSE FALSE
host  IL66B  5 16632 424 424 56.2 22.8 0.45 TRUE FALSE
host  NL70B  5 16727  48  48 64   26   0.061 FALSE FALSE
host  NL71B  5 16490  25  25 65.8 23.2 0.032 FALSE FALSE
```

**Changes how B7 should read**
1. **A gene-set test has already been run, so the footer's "no gene-set test has been run" is wrong.** The two-contrast pipeline's Broad GSEA for host_therapy finished (exit 0); it just never reached the results folder because a later plotting step crashed. It is in `work_de/9a/6e148a…`, and I copied the reports to `verify/broad_gsea_existing/`.
   - 4,843 sets were tested; **none reaches FDR < 0.25 in either direction** (lowest FDR 0.33).
   - Top by nominal p, higher in recurrences: antigen presentation, allograft rejection, interferon-gamma (NES about 1.85–1.92).
   - Lower in recurrences: KEGG_RIBOSOME (NES −2.04, p 0.042) and translation initiation.
2. **How the pipeline table was built:**
   - The log2FoldChange column is ashr-shrunk; pvalue and padj come from the unshrunk Wald test.
   - The model is `~ 0 + Classification`, with salmon gene lengths used as offsets.
   - Genes were filtered on all nine samples before the six tumours were subset out.
3. **Local run of the same R code (DESeq2 1.50.2 without ashr) — preliminary; the grid job will give the exact numbers:**
   - It finds the same 38 genes as the pipeline (identical set). p-values differ by up to 3% relative, which is expected across DESeq2 versions.
   - **The 38-gene list is unstable.** Dropping one primary or one NL tumour keeps only 8–19 of the 38. Dropping IL66B keeps 37 and raises the total to 398.
   - Only Mmp13, Ceacam4, Hbb and Ocstamp survive all six leave-one-out fits. Mmp13 stays 35–50-fold down in every fit (padj ≤ 1e-18).
   - **The true labels are not special.** They rank 2nd of the 10 possible 3-v-3 splits (38 genes; median of the others is 9).
   - The largest split is IL69B+NL70B+NL71B against the other three, at 2,752 genes. It follows immune infiltration (Ptprc, Cd68 and Aif1 high) against white matter and neurons (Mbp, Plp1, Pcp4 low). That axis cuts across the therapy groups.
   - NL70B+NL71B against the other four gives 502 genes. 11 of the 15 "up" genes are significant in that comparison.
   - **Only the down genes are shared by all three recurrences.** On the 23 down genes, IL66B moves about 100% as far as the NL pair. On the 15 up genes (Cxcl10, Batf2, Il2ra, Ighm) it moves a median of about 29%, so the interferon/lymphocyte signal comes from NL70B and NL71B.
   - **fgsea looks strong but does not beat chance.** Against the true labels it gives 14 hallmark sets at padj < 0.05 (interferon, allograft and TNF sets up, EMT down), all stable across leave-one-out. But every one of the 10 relabellings gives as many or more significant sets; the true split ranks 9th of 10 for hallmark, KEGG and the brain sets. The local GO BP run used a stand-in file, so only the grid job has real GO BP numbers.
   - **Possible human-read leakage into the rat counts.** The human ortholog count is at least the rat count in some tumour for Tspan6, Itgb4, Ass1, Batf2, Lrrc15 and Mmp13; for Mmp13 this is trivial (7 human against 5 rat reads in one tumour). Tspan6 is the most suspect: human is higher in all six tumours and falls with recurrence on both sides. 15 of the 38 genes (including Hbb, Hba-a1 and RT1-S3) have no 1:1 human ortholog in the map, so they are unchecked. Across the six tumours, the rat ribosomal-protein share tracks the human read share (Spearman 0.83), which fits leakage but does not prove it at n = 6.
   - One of the 38, LOC120093067, carries 0.66–5.1 million reads per tumour. It may be an rRNA locus, but I have not checked.
   - **Suggested B7 wording:** the host shows no gene set that survives either FDR or the relabelling null. The consistent host change is a set of genes lower in all three recurrences (Mmp13, Ceacam4, Hbb). The immune "up" genes are NL70B/NL71B-specific and line up with the immune-infiltration axis.

**Files** (nothing committed or pushed; the grid clone was not touched):
- Local scripts in `C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study/ANALYSIS/host/verify/`:
  - `hv_json.R`
  - `host_verify_de.R`
  - `host_verify_gsea.R`
  - `run_host_verify.sbatch`
  - `run_de_therapy.sbatch`
  - `metadata_therapy.csv`
  - `contrasts_therapy.csv`
  - `host_de_therapy_params.yaml`
- The same files are on the grid in `/rs/rs_grp_oschome/go2432/u251_host/verify/`.
- Checkpoint block added as the first section of `C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study/RESUME.md`, under the title and intro.
- Local dry-run outputs: `C:/Users/grego/AppData/Local/Temp/claude/c--Users-grego-OneDrive-Desktop-CTSpinoPelvic1K-1/f1bdbd78-151f-470b-b703-dd9af9b3fecc/scratchpad/localtest/`

## Deconvolution

Nothing is still running. The Bowman failure is diagnosed and fixed without editing CIBERSORT.R, and all four host signatures have outputs. No cell type differs between primaries and recurrences in any signature.

**Job status**
- **40473363_1 (LM22):** it was requeued after the 17:01 preemption and finished at 17:21, so there was nothing left to recover.
- **40473363_2 (Zhang 2014):** completed.
- **40473363_3 (Bowman 2016):** failed.
- **Rerun:** prep job 40475433 queued array 40475461, tasks 1–4. All completed by 18:33. The work directory is `/rs/rs_grp_oschome/go2432/u251_host/deconv_fix/` (scripts/, results/, logs/). The grid clone was not touched.

**Why Bowman failed**
- **The crash:** with seed 42, permutation draw 169 gave no positive weight at any of the three nu values. CoreAlg then computes 0/0, all three RMSEs are NaN, `which.min()` returns nothing, and `out[[mn]]` stops the whole run. In a local libsvm replica this happens about once per 1,000–1,500 draws with 4 columns. With 22 columns (LM22) it essentially never happens.
- **The signature file is clean:** no NA, no negative values, no zero-variance rows, condition number 11.6, 189/200 genes in the mixture. 11 genes are zero in every host sample; that is harmless.
- **The real problem:** four myeloid columns cannot describe whole brain. On these genes every sample correlates with every column at |r| < 0.08, and the three controls correlate negatively with all four. The sorted-microglia columns carry transcripts from neighbouring cells (KIF5A, NEFL, TPPP, S100B, COL3A1), which in bulk tissue come from neurons and glia.
- **The mixture's low retained TPM (34–44%) is not a mapping fault.** Most of what is lost is rat mitochondrial genes (AY172581.*, about 23%), rRNA and small ncRNAs.

**What I changed**
1. `ANALYSIS/cibersort/run_v104_guarded.R` and `.sbatch`: CIBERSORT.R is sourced unchanged. The permutation null uses the same seed and the same draws, and every fit is still CoreAlg. A draw that cannot be fitted is left out of the null and counted, which is the conservative choice. Check: repeating LM22 and Zhang through it reproduces the original fractions, r, RMSE and P exactly (difference 0).
2. `ANALYSIS/host/build_combined_signature.py`: six non-myeloid Zhang cell types plus the four Bowman myeloid states (1,510 genes, 1,456 in the mixture, condition number 102). Zhang's own microglia column is left out. The selection rule removes most of the carried-over transcripts. The two references come from different platforms, so some batch difference is built in.
3. `ANALYSIS/host/run_host_deconv_fix.sbatch` (prep, then the 4-task array), `summarise_host_deconv.py`, and a build log `signature_host_zhang6_bowman4_log.json`. The combined signature matrix stays on the grid only, because it is derived from third-party data.

**Fit per signature** (r, RMSE and P are CIBERSORT's; RMSE is in standardised units)

| Signature | Tumours r | Tumours RMSE | Tumours P | Controls r | Controls P | Usable? |
|---|---|---|---|---|---|---|
| Bowman alone | 0.01–0.08 | 1.01–1.17 | 0.21–0.42 | about −0.01 | 0.57–0.62 | No |
| LM22 | 0.29–0.44 | 0.90–0.99 | 0.004–0.011 | −0.05 to −0.03 | 0.92–0.98 | No |
| Zhang 2014 | 0.22–0.46 | 0.93–0.98 | 0.005–0.014 | 0.31–0.78 | 0.001–0.007 | Weak |
| Combined | IL 0.25–0.36, NL 0.77 | 0.87–0.98 | ≤ 0.005 | 0.32–0.75 | ≤ 0.004 | Best of the four, still weak for IL |

- **Bowman alone is too poor to interpret.**
- **LM22 fails its own check:** nude rats have almost no T cells, yet it assigns 7–20% T cells in tumours and 26–38% in controls. Its fractions (monocytes 28–50%, M2 macrophages 10–42%) are shares within a leukocyte-only model, not tissue fractions.

**Fractions per sample**, in the order IL67B, IL68B, IL69B | IL66B, NL70B, NL71B | IL64B, N168B, N269B:
- **Combined signature:**
  - TAM_BMDM: .115 .080 .272 | .065 .433 .498 | .003 0 .032
  - TAM_MG: .050 .093 .121 | .044 .018 .066 | 0 .006 .024
  - normal_MG: .039 .070 .019 | .040 .015 0 | .014 .010 .004
  - monocyte: .040 .048 .075 | .037 .070 .052 | 0 0 .005
  - endothelial: .101 .107 .095 | .092 .067 .075 | .066 .044 .089
  - neuron: .145 .142 .137 | .210 .135 .052 | .444 .382 .346
  - astrocyte: .150 .113 .100 | .187 .077 .111 | .146 .212 .160
- **Zhang 2014:**
  - microglia: .177 .247 .415 | .202 .273 .444 | .034 .042 .058
  - endothelial: .151 .151 .148 | .131 .151 .155 | .076 .054 .100
- **Cross-check between references:** total myeloid in the combined signature tracks Zhang's microglia column (Spearman 0.95 over 9 samples, 0.89 over the 6 tumours). Both are higher in tumours than in controls (Welch p 0.003).

**Primaries (3) v recurrences (3)**
I used a two-sided Welch t-test, with Benjamini–Hochberg correction across cell types within each signature. I also ran the exact Mann-Whitney test, but with 3 v 3 its smallest possible p is 0.10, so it cannot reach 0.05.
- **Combined:** nothing reaches significance. TAM_BMDM 0.16 v 0.33 (p 0.32), TAM_MG 0.09 v 0.04 (p 0.15), endothelial 0.10 v 0.08 (p 0.08); all q ≥ 0.58.
- **Zhang 2014:** every p ≥ 0.42.
- **LM22:** smallest is M2 0.35 v 0.18 (p 0.13, q 0.64).
- **Bowman alone:** TAM_MG gives p 0.036 (q 0.11). This is not evidence, because the fit behind it is r ≤ 0.08.

**Pattern found after looking, so exploratory only**
The recurrence rise in TAM_BMDM comes from NL70B and NL71B (0.43, 0.50), not IL66B (0.07). Comparing the four IL tumours with the two NL tumours gives 0.13 v 0.47, p 0.005. NL70B and NL71B are also the only tumours the combined signature fits well (r 0.77). This matches the host differential-expression finding: the host difference is NL70B/NL71B against the IL tumours. It is an IL/NL split, not a primary/recurrence one.

`RESUME.md` has the checkpoint under "## 2026-09-26 evening: host deconvolution". Nothing was committed.

Files are in `C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study/`:
- `ANALYSIS/cibersort/run_v104_guarded.R`
- `ANALYSIS/cibersort/run_v104_guarded.sbatch`
- `ANALYSIS/host/build_combined_signature.py`
- `ANALYSIS/host/run_host_deconv_fix.sbatch`
- `ANALYSIS/host/summarise_host_deconv.py`
- `ANALYSIS/host/host_deconv_summary.tsv`
- `ANALYSIS/host/signature_host_zhang6_bowman4_log.json`
- `ANALYSIS/cibersort/results/v104/fractions_v104_host_lm22.tsv`
- `ANALYSIS/cibersort/results/v104/fractions_v104_host_zhang2014_brain_mean7_tpm.tsv`
- `ANALYSIS/cibersort/results/v104/fractions_v104_host_bowman2016_myeloid_mean4.tsv`
- `ANALYSIS/cibersort/results/v104/fractions_v104_host_zhang6_bowman4.tsv`
- `RESUME.md`
