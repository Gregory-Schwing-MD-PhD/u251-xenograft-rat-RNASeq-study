#!/bin/bash
# Drug discovery aimed at the tumour stem-cell program (gsc_screen.py). EXTRA="<list files>" adds published GSC sets.
#   sbatch -D ANALYSIS/gsc_drugs ANALYSIS/gsc_drugs/run_gsc.sbatch
#SBATCH --job-name=u251_gsc
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/gsc_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
PY=/rs/rs_grp_oschome/go2432/envs/u251_drug/bin/python
unset PYTHONPATH LD_LIBRARY_PATH
cd "$U/ANALYSIS/gsc_drugs"
"$PY" -u gsc_screen.py --dsigdb "$U/ANALYSIS/refs/pathways/human/hs.dsigdb.v1.0.gmt" \
  --de "$U/ANALYSIS/results_therapy_v3/tables/differential/therapy_impact.deseq2.results.tsv" \
  --tpm "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_tpm.tsv" \
  --varn "$U/ANALYSIS/cibersort/signatures_online/varn2022_glass_sc12.txt" \
  --extra ${EXTRA:-} --out results
