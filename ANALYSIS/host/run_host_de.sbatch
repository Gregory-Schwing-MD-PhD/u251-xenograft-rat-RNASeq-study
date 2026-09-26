#!/bin/bash
# nf-core/differentialabundance 1.5.0 on the host (rat) gene counts, the published human settings with rat
# references (host_de_params.yaml): Recurrent vs Primary (3 v 3) and Tumour vs Control brain (6 v 3), each contrast
# fitted on its own samples (differential_subset_to_contrast_samples).
#   sbatch -D ANALYSIS/host --dependency=afterok:<rnaseq> ANALYSIS/host/run_host_de.sbatch
#SBATCH --job-name=u251_hostde
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=logs/de_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
R=/rs/rs_grp_oschome/go2432/u251_host/results_rnaseq/star_salmon
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
mkdir -p /rs/rs_grp_oschome/go2432/u251_host/tmp
L=/rs/rs_grp_oschome/go2432/u251_host/launch_de; mkdir -p "$L"; cd "$L"   # own launch directory (see run_host_rnaseq)
nextflow run nf-core/differentialabundance -r 1.5.0 -profile singularity \
  -c "$U/nextflow.config" -c "$U/ANALYSIS/host/host.config" \
  --input "$U/ANALYSIS/host/metadata_host.csv" --contrasts "$U/ANALYSIS/host/contrasts_host.csv" \
  --matrix "$R/salmon.merged.gene_counts.tsv" --transcript_length_matrix "$R/salmon.merged.gene_lengths.tsv" \
  -params-file "$U/ANALYSIS/host/host_de_params.yaml" -w /rs/rs_grp_oschome/go2432/u251_host/work_de -ansi-log false
ls /rs/rs_grp_oschome/go2432/u251_host/results_de/tables/differential/
echo "host DE done"
