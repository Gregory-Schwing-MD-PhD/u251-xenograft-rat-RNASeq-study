#!/bin/bash
# nf-core/differentialabundance 1.5.0 for host_therapy ONLY (Recurrent_U2 vs Primary_U2, six tumours), GSEA on,
# launched as ANALYSIS/host/run_host_de.sbatch launches the two-contrast run, with the same settings
# (host_de_therapy_params.yaml = host_de_params.yaml with a new outdir and filtering_min_abundance 0).
#
# Why the inputs change: the two-contrast run (results_de) crashed after DESeq2 because DESEQ2_NORM subset the
# normalised matrix to the first contrast's six tumours, and PLOT_EXPLORATORY / SHINYNGS_APP then could not find the
# three control brains named in the nine-sample sheet. Its host_therapy GSEA had finished (exit 0) in
# work_de/9a/6e148a... but was never published. Here the sample sheet is the six tumours, and the counts and lengths
# are the six tumour columns of the 17,451 genes that the original CUSTOM_MATRIXFILTER kept on all nine samples
# (>= 10 reads in >= 1 sample); filtering_min_abundance 0 keeps that universe unchanged. Every DESeq2 and GSEA input
# therefore equals the original run's, so this run is expected to reproduce results_de's host_therapy table and the
# unpublished GSEA exactly, and to finish the exploratory plots and report the crashed run never produced.
#   sbatch -D /rs/rs_grp_oschome/go2432/u251_host/verify /rs/rs_grp_oschome/go2432/u251_host/verify/run_de_therapy.sbatch
#SBATCH --job-name=u251_hostde_th
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/de_therapy_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=/rs/rs_grp_oschome/go2432/u251_host
V=$H/verify
R=$H/results_rnaseq/star_salmon
N=$H/results_de_therapy
ORIG=$H/results_de/tables/differential/host_therapy.deseq2.results.tsv
OLDG=$H/work_de/9a/6e148abc127c567ca55315e3cf8f64
IN=$V/de_therapy_inputs; mkdir -p "$IN"

# six tumour columns of the original run's gene universe (the gene_id column of its host_therapy table)
subset6() {
  awk -F'\t' -v ids="$ORIG" '
    BEGIN { while ((getline l < ids) > 0) { split(l, f, "\t"); if (f[1] != "gene_id") keep[f[1]] = 1 }
            n = split("IL67B IL68B IL69B IL66B NL70B NL71B", want, " ") }
    NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i
              printf "gene_id\tgene_name"; for (j = 1; j <= n; j++) printf "\t%s", want[j]; printf "\n"; next }
    ($1 in keep) { printf "%s\t%s", $(col["gene_id"]), $(col["gene_name"])
                   for (j = 1; j <= n; j++) printf "\t%s", $(col[want[j]]); printf "\n" }' "$1" > "$2"
}
subset6 "$R/salmon.merged.gene_counts.tsv" "$IN/gene_counts_therapy.tsv"
subset6 "$R/salmon.merged.gene_lengths.tsv" "$IN/gene_lengths_therapy.tsv"
echo "counts rows (with header): $(wc -l < "$IN/gene_counts_therapy.tsv"); lengths: $(wc -l < "$IN/gene_lengths_therapy.tsv"); original table: $(wc -l < "$ORIG")"
head -2 "$IN/gene_counts_therapy.tsv"

source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
mkdir -p "$H/tmp"
L=$H/launch_de_therapy; mkdir -p "$L"; cd "$L"   # own launch directory: no Nextflow lock shared with other runs
set +e
nextflow run nf-core/differentialabundance -r 1.5.0 -profile singularity \
  -c "$U/nextflow.config" -c "$U/ANALYSIS/host/host.config" \
  --input "$V/metadata_therapy.csv" --contrasts "$V/contrasts_therapy.csv" \
  --matrix "$IN/gene_counts_therapy.tsv" --transcript_length_matrix "$IN/gene_lengths_therapy.tsv" \
  -params-file "$V/host_de_therapy_params.yaml" -w "$H/work_de_therapy" -ansi-log false -resume
rc=$?
set -e

# compare with the original run: the DE table and the unpublished GSEA
C=$V/de_therapy_compare.txt
{
  echo "nextflow exit status: $rc   ($(date))"
  NEW=$N/tables/differential/host_therapy.deseq2.results.tsv
  [ -s "$NEW" ] || NEW=$(find "$H/work_de_therapy" -name host_therapy.deseq2.results.tsv -size +1k 2>/dev/null | head -1)
  if [ -n "$NEW" ] && [ -s "$NEW" ]; then
    if cmp -s "$NEW" "$ORIG"; then echo "DE table: byte-identical to results_de ($NEW)"
    else echo "DE table DIFFERS from results_de ($NEW); first differences:"; diff "$NEW" "$ORIG" | head -20; fi
    echo "padj < 0.05: $(awk -F'\t' 'NR > 1 && $6 != "NA" && $6 < 0.05' "$NEW" | wc -l)"
  else echo "no host_therapy DE table found"; fi
  for side in Recurrent_U2 Primary_U2; do
    G=$(find "$N" "$H/work_de_therapy" -name "host_therapy.combined_rat.gsea_report_for_${side}.tsv" 2>/dev/null | head -1)
    O=$OLDG/host_therapy.combined_rat.gsea_report_for_${side}.tsv
    if [ -n "$G" ]; then
      if diff -q <(cut -f1,4-9 "$G") <(cut -f1,4-9 "$O") > /dev/null; then echo "GSEA $side: identical to the unpublished run ($G)"
      else echo "GSEA $side DIFFERS from the unpublished run ($G)"; diff <(cut -f1,4-9 "$G") <(cut -f1,4-9 "$O") | head -10; fi
      echo "GSEA $side: sets $(($(wc -l < "$G") - 1)); FDR < 0.25: $(awk -F'\t' 'NR > 1 && $8 < 0.25' "$G" | wc -l); nominal p < 0.01: $(awk -F'\t' 'NR > 1 && $7 < 0.01' "$G" | wc -l)"
    else echo "GSEA $side: no report found"; fi
  done
  echo "published files:"; find "$N" -maxdepth 3 -type d 2>/dev/null | head -40
} > "$C" 2>&1
cat "$C"
exit $rc
