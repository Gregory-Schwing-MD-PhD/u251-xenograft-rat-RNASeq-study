#!/usr/bin/env bash
# Copy the pipeline outputs the slide scripts read (SLIDES/u251_paths.py) from the grid u251 directory into the same
# gitignored paths of a local clone. Run from the repository root; network through WSL on Windows:
#   wsl -e bash -lc 'cd /mnt/c/.../u251-xenograft-murine-RNASeq-study && bash SLIDES/fetch_grid_outputs.sh'
set -euo pipefail
REMOTE=go2432@grid.wayne.edu:/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
G=ANALYSIS/results_therapy_v3/report/gsea/therapy_impact/combined_human
P=therapy_impact.combined_human
files=(
  ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv
  ANALYSIS/results_therapy_v3/tables/processed_abundance/all.rlog.tsv
  ANALYSIS/metadata_full.csv
  ANALYSIS/results_ruvseq/ruvseq_baseline_de.tsv
  ANALYSIS/results_ruvseq/ruvseq_adjusted_de_k1.tsv
  ANALYSIS/results_ruvseq/ruvseq_adjusted_de_k2.tsv
  ANALYSIS/results_ruvseq/ruvseq_concordance_summary.csv
  ANALYSIS/results_ruvseq/ruvseq_estimated_W_factors.csv
  ANALYSIS/results_decontamination/contamination_bias_per_gene.csv
  ANALYSIS/results_decontamination/contaminated_genes.tsv
  $G/$P.gsea_report_for_Primary_U2.tsv
  $G/$P.gsea_report_for_Recurrent_U2.tsv
  $G/$P.ranked_gene_list_Recurrent_U2_versus_Primary_U2.tsv
)
for s in KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION KEGG_RIBOSOME \
         REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY REACTOME_SELENOAMINO_ACID_METABOLISM \
         REACTOME_CELLULAR_RESPONSE_TO_STARVATION REACTOME_EUKARYOTIC_TRANSLATION_INITIATION; do
  files+=("$G/$P.$s.tsv")
done
for f in "${files[@]}"; do
  mkdir -p "$(dirname "$f")"
  rsync -t -e "ssh -o BatchMode=yes" "$REMOTE/$f" "$f"
done
# the pipeline's own running-sum plots of the six leading sets (enplots_composite.png)
rsync -rt -e "ssh -o BatchMode=yes" --include="$P.enplot_*.png" --exclude='*' "$REMOTE/$G/" "$G/"
echo "fetched ${#files[@]} files and the enrichment plots"
