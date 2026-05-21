# RUVSeq contamination-adjustment results

This directory holds the output of
[`ANALYSIS/ruvseq_contamination_adjustment.R`](../ruvseq_contamination_adjustment.R),
run on the cluster via [`run_ruvseq.sh`](../../run_ruvseq.sh).

## Why this analysis exists

The xenograft RNA-seq pipeline (xengsort PDX split → Salmon on the human
transcriptome) leaves residual rat reads from highly **sequence-conserved**
genes (ribosomal proteins, translation factors, histones, tubulins) misassigned
to their human orthologs. This contamination scales with the rat tissue
fraction, which differs systematically between **Recurrent** (lower graft
fraction) and **Primary** (higher graft fraction) tumours — so a naive
Primary-vs-Recurrent contrast can manufacture artifactual DE in conserved gene
families (e.g. a spurious "translation/ribosome downregulation in recurrence").

The three **Control** samples (`IL64B`, `N168B`, `N269B`; rat brain only,
~0.3–4.9% graft) went through the identical pipeline, so any human-gene signal
they carry **is** contamination. `RUVSeq::RUVs` (Risso et al., *Nat Biotechnol*
2014) estimates factors of unwanted variation **W** from these Control
replicates; **W** is then added as covariate(s) to the DESeq2 model for the
Primary-vs-Recurrent contrast. Controls anchor **W** but are excluded from the
contrast itself.

## Files

| File | Description |
|------|-------------|
| `ruvseq_baseline_de.tsv` | Unadjusted DESeq2 (`~ Classification`), Recurrent vs Primary. Reference. |
| `ruvseq_adjusted_de_k1.tsv` / `_k2.tsv` / `_k3.tsv` | DESeq2 with design `~ W_1..W_k + Classification`, Recurrent vs Primary. log2FC > 0 = up in Recurrent. |
| `ruvseq_estimated_W_factors.csv` | Estimated W per sample, per k (long format: `sample, Classification, k_setting, factor, value`). |
| `ruvseq_concordance_summary.csv` | Per-setting metrics vs baseline: Spearman(log2FC), Jaccard(sig sets), n sig up/down, top-20 overlap, Primary-vs-Recurrent silhouette. `recommended` column flags the chosen k. |
| `ruvseq_pca_before_after.pdf` | 2×2 PCA grid (baseline + k=1,2,3) coloured by Classification; adjusted panels have the W effect residualised out. |
| `ruvseq_recommended_k.txt` | Single integer (1, 2, or 3): the k with the highest Primary-vs-Recurrent silhouette. |

## Method notes

- **Pre-filter:** keep genes with count ≥ 10 in ≥ 3 samples.
- **RUVs:** all retained genes as negative controls (`cIdx = rownames(set)`);
  `scIdx` = a single replicate group of the 3 Control samples. A single
  3-replicate group spans at most 2 effective factors, so **k=3 is reported for
  completeness but may be near-degenerate** — treat it with caution.
- **Significance threshold:** `padj < 0.05` and `|log2FC| > 1`.
- **k selection:** maximum Primary-vs-Recurrent silhouette in PC1–PC2 over
  k ∈ {1,2,3}. The silhouette heuristic is **noisy at this sample size (n=3 per
  group)**; k=2 is the conventional RUVSeq default and a reasonable manual
  override if the automatic pick is k=1 or k=3.
- **Significance is not interpreted here.** Adopt a k, then point the
  publication-figure script at `ruvseq_adjusted_de_k{k}.tsv` if desired.

## Metadata provenance

The graft/host/both fractions in `ANALYSIS/metadata_full.csv` are **not**
hand-typed. They are derived from the per-sample xengsort 2.1.0 classify logs
(`ANALYSIS/xengsort_out/<sample>.txt`) by
[`ANALYSIS/build_metadata_from_xengsort.py`](../build_metadata_from_xengsort.py),
which parses the `## Classification Statistics` counts block and computes each
percentage as `count / total * 100`. The design assignment (sample →
Control/Primary/Recurrent) lives in `ANALYSIS/metadata_base.csv`;
`run_ruvseq.sh` regenerates `metadata_full.csv` from it on every run.

## Reproduce

```bash
sbatch run_ruvseq.sh          # on Warrior; rebuilds metadata, logs to slurm_logs/

# the two steps it runs, manually:
python3 ANALYSIS/build_metadata_from_xengsort.py \
    --base ANALYSIS/metadata_base.csv \
    --xengsort-dir ANALYSIS/xengsort_out \
    --out ANALYSIS/metadata_full.csv

Rscript ANALYSIS/ruvseq_contamination_adjustment.R \
    ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv \
    ANALYSIS/metadata_full.csv \
    ANALYSIS/results_ruvseq
```
