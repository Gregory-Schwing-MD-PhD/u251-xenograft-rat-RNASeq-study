#!/bin/bash
# Step 1 of the IL68B-held-out rerun: nf-core/differentialabundance 1.5.0 exactly as run_de_pdx_v3.sh ran it for the
# published analysis (same params file, counts, gene lengths and contrast), with IL68B removed from the sample sheet.
# Produces DESeq2 results, rlog matrix and the Broad GSEA report for the five-tumour contrast.
#   sbatch -D ANALYSIS/holdout_IL68B ANALYSIS/holdout_IL68B/run_de.sbatch
#SBATCH --job-name=u251_ho_de
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/de_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=$U/ANALYSIS/holdout_IL68B
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
mkdir -p "$XDG_RUNTIME_DIR" "$NXF_SINGULARITY_CACHEDIR"
cd "$U"                                   # the repository's nextflow.config (slurm executor, -q primary) applies, as in May
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    -c "$H/holdout.config" \
    --input "$H/ANALYSIS/metadata_therapy.csv" \
    --contrasts "$U/ANALYSIS/contrasts_therapy.csv" \
    --matrix "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv" \
    -params-file "$U/therapy_v3_params.yaml" \
    --outdir "$H/results_therapy_v3_noIL68B" \
    -w "$H/work" \
    -ansi-log false
ls "$H/results_therapy_v3_noIL68B/tables/differential/"
echo "DE done"
