#!/bin/bash
# CIBERSORT v1.04 through run_v104_guarded.R (a permutation draw no nu-SVR can fit is recorded, not fatal), one signature
# per array task. Lines of $JOBS: name<TAB>signature<TAB>perm<TAB>mixture. Same image and R library as run_v104.sbatch.
#   sbatch -D <dir> --array=1-N --export=ALL,JOBS=<list>,OUT=<outdir>,SCRIPTS=<dir with run_v104_guarded.R>,CIBERSORT_R=<vendor/CIBERSORT.R> run_v104_guarded.sbatch
#SBATCH --job-name=u251_cs104g
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/cs104g_%A_%a.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
C=$U/ANALYSIS/cibersort
S=${SCRIPTS:-$C}
OUT=${OUT:-$C/results/v104}
CS=${CIBERSORT_R:-$C/vendor/CIBERSORT.R}
SIF=/rs/rs_grp_oschome/go2432/singularity_cache/go2432-bioconductor.sif
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
RLIB=/rs/rs_grp_oschome/go2432/Rlib_cibersort
mkdir -p "$OUT"
unset PYTHONPATH LD_LIBRARY_PATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID}; mkdir -p "$XDG_RUNTIME_DIR"
RTMP=/dev/shm/${USER}_rtmp_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID:-0}; mkdir -p "$RTMP"
trap 'rm -rf "$RTMP" "$XDG_RUNTIME_DIR"' EXIT
export SINGULARITYENV_R_LIBS_USER=$RLIB SINGULARITYENV_TMPDIR=$RTMP SINGULARITYENV_CIBERSORT_R=$CS CIBERSORT_R=$CS
R() { "$SING" exec --bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432 "$SIF" "$@"; }
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOBS:?JOBS list required}")
name=$(echo "$line" | cut -f1); sig=$(echo "$line" | cut -f2); perm=$(echo "$line" | cut -f3); perm=${perm:-1000}
mix=$(echo "$line" | cut -f4); mix=${mix:-$C/results/mixture_u251_tpm.txt}
echo "signature $name: $sig (perm $perm); mixture $mix; CIBERSORT $CS"
cd "$S"
R Rscript run_v104_guarded.R "$sig" "$mix" "$OUT/fractions_v104_${name}.tsv" "$perm"
