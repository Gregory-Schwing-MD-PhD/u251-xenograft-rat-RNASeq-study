#!/bin/bash
# Leave-one-tumour-out separation table (loo_separation.R) on the human reads; rerun on the host reads with
# COUNTS=<host gene counts> LABEL=host once the host quantification exists.
#   sbatch -D ANALYSIS/holdout_separation ANALYSIS/holdout_separation/run_loo.sbatch
#SBATCH --job-name=u251_loo
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/loo_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
SIF=/rs/rs_grp_oschome/go2432/singularity_cache/go2432-bioconductor.sif
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID}; mkdir -p "$XDG_RUNTIME_DIR"
COUNTS=${COUNTS:-$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv}
LABEL=${LABEL:-human}
cd "$U/ANALYSIS/holdout_separation"
"$SING" exec --bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432 "$SIF" \
  Rscript loo_separation.R "$COUNTS" loo_separation "$LABEL"
