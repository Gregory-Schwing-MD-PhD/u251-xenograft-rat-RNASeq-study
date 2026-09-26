#!/bin/bash
# Bulk deconvolution only, from the signature matrices and mixture already in results/ (deconv_only.py).
#   sbatch -D ANALYSIS/cibersort ANALYSIS/cibersort/run_deconv.sbatch
#SBATCH --job-name=u251_deconv
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/deconv_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
PY=$HOME/mambaforge/envs/spineps/bin/python
unset PYTHONPATH
cd "$U/ANALYSIS/cibersort"
"$PY" -u deconv_only.py --results results --cohorts "$U/ANALYSIS/sample_cohorts.csv" --perm 1000 \
  neftel4_confident neftel4_allcells neftel4_confident_plus_nonmalignant
