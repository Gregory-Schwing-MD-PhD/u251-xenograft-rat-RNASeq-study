#!/bin/bash
# Neftel 2019 reference, CIBERSORT-format signature matrices, the CIBERSORTx single-cell reference, the
# leave-one-tumour-out check and a CIBERSORT-core deconvolution of the ten bulk profiles (cibersort_neftel.py).
#   sbatch -D ANALYSIS/cibersort ANALYSIS/cibersort/run.sbatch
#SBATCH --job-name=u251_cibersort
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/cibersort_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
PY=$HOME/mambaforge/envs/spineps/bin/python
unset PYTHONPATH
cd "$U/ANALYSIS/cibersort"
"$PY" -u cibersort_neftel.py \
  --neftel "$U/ANALYSIS/refs/neftel_2019" \
  --tpm "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv" \
  --cohorts "$U/ANALYSIS/sample_cohorts.csv" \
  --out results
ls -la results
