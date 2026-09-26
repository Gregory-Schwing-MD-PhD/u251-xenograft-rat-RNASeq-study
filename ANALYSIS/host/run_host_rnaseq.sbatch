#!/bin/bash
# nf-core/rnaseq on the host-only rat reads (head job; tasks go to slurm through the repository's nextflow.config).
#   sbatch -D ANALYSIS/host --dependency=afterok:<extract> ANALYSIS/host/run_host_rnaseq.sbatch
#SBATCH --job-name=u251_hostrna
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/rnaseq_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
mkdir -p "$XDG_RUNTIME_DIR"
# own launch directory: from the repository root, -resume picked up (and could not lock) a concurrent Nextflow run
L=/rs/rs_grp_oschome/go2432/u251_host/launch_rnaseq; mkdir -p "$L"; cd "$L"
nextflow run nf-core/rnaseq -r 3.22.2 -profile singularity \
  -c "$U/nextflow.config" -c "$U/rnaseq_custom.config" -c "$U/ANALYSIS/host/host.config" \
  -params-file "$U/ANALYSIS/host/host_rnaseq_params.yaml" \
  -w /rs/rs_grp_oschome/go2432/u251_host/work -resume -ansi-log false
ls /rs/rs_grp_oschome/go2432/u251_host/results_rnaseq/star_salmon/ | grep merged
echo "host rnaseq done"
