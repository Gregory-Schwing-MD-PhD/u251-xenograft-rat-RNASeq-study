#!/bin/bash
# CIBERSORT v1.04 on the bulk mixture, one signature per array task (lines of v104_jobs.txt: name<TAB>signature path).
# R from the lab's Bioconductor image; e1071 and preprocessCore installed once into a CephFS library (SETUP=1).
#   sbatch -D ANALYSIS/cibersort --export=ALL,SETUP=1 run_v104.sbatch            # install only
#   sbatch -D ANALYSIS/cibersort --array=1-N run_v104.sbatch                     # run
#SBATCH --job-name=u251_cs104
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=3
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/cs104_%A_%a.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
C=$U/ANALYSIS/cibersort
SIF=/rs/rs_grp_oschome/go2432/singularity_cache/go2432-bioconductor.sif
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
RLIB=/rs/rs_grp_oschome/go2432/Rlib_cibersort
mkdir -p "$RLIB" "$C/results/v104"
unset PYTHONPATH LD_LIBRARY_PATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID}; mkdir -p "$XDG_RUNTIME_DIR"
export SINGULARITYENV_R_LIBS_USER=$RLIB
R() { "$SING" exec --bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432 "$SIF" "$@"; }
if [ "${SETUP:-0}" = 1 ]; then
  R Rscript -e 'lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()));
    if (!requireNamespace("e1071", quietly=TRUE)) install.packages("e1071", lib=lib, repos="https://cloud.r-project.org");
    if (!requireNamespace("preprocessCore", quietly=TRUE)) BiocManager::install("preprocessCore", lib=lib, update=FALSE, ask=FALSE, configure.args="--disable-threading");
    library(e1071); library(preprocessCore); cat("ok", as.character(packageVersion("e1071")), as.character(packageVersion("preprocessCore")), "\n")'
  exit 0
fi
# JOBS=<list> selects another job list; an optional 4th field names the mixture (default: the human-read mixture)
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${JOBS:-$C/v104_jobs.txt}")
name=$(echo "$line" | cut -f1); sig=$(echo "$line" | cut -f2); perm=$(echo "$line" | cut -f3); perm=${perm:-1000}
mix=$(echo "$line" | cut -f4); mix=${mix:-$C/results/mixture_u251_tpm.txt}
echo "signature $name: $sig (perm $perm); mixture $mix"
cd "$C"
R Rscript run_v104.R "$sig" "$mix" "$C/results/v104/fractions_v104_${name}.tsv" "$perm"
