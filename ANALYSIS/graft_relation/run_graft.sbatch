#!/bin/bash
# Graft/host fraction against every result (graft_relation.R), human reads, six tumours.
#   sbatch -D ANALYSIS/graft_relation ANALYSIS/graft_relation/run_graft.sbatch
#SBATCH --job-name=u251_graft
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=logs/graft_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
SIF=/rs/rs_grp_oschome/go2432/singularity_cache/go2432-bioconductor.sif
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID}; mkdir -p "$XDG_RUNTIME_DIR"
cd "$U/ANALYSIS/graft_relation"
"$SING" exec --bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432 "$SIF" \
  Rscript graft_relation.R "$U" results
