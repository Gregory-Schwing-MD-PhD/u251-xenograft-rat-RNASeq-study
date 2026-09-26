#!/bin/bash
# Host (rat) cell-state deconvolution: host_deconv_prep.py builds the rat->human 1:1 mixture and the signatures
# (LM22 as is; Zhang 2014 brain cell types and Bowman 2016 glioma myeloid, mouse->human 1:1, TPM), then queues one
# CIBERSORT v1.04 task per signature (../cibersort/run_v104.sbatch with JOBS=results/v104_host_jobs.txt).
# seq-ImmuCC is left out: its values are percentile-scaled read counts, not TPM (refs_online/SOURCES.tsv).
#   sbatch -D ANALYSIS/host --dependency=afterany:<rnaseq> ANALYSIS/host/run_host_deconv.sbatch
#SBATCH --job-name=u251_hostcs
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/deconv_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=$U/ANALYSIS/host
TPM=/rs/rs_grp_oschome/go2432/u251_host/results_rnaseq/star_salmon/salmon.merged.gene_tpm.tsv
PY=/rs/rs_grp_oschome/go2432/envs/u251_drug/bin/python
unset PYTHONPATH LD_LIBRARY_PATH
[ -s "$TPM" ] || { echo "host TPM missing: $TPM (rnaseq not finished)"; exit 3; }
R=$H/refs_online
cd "$H"
"$PY" host_deconv_prep.py --tpm "$TPM" --rat2human "$R/rat_human_1to1.tsv" --mouse2human "$R/mouse_human_1to1.tsv" \
  --lm22 "$U/ANALYSIS/cibersort/signatures_online/lm22_newman2015.txt" \
  --mouse-mean "$R/zhang2014_brain_mean7_tpm.txt" "$R/bowman2016_myeloid_mean4.txt" --out results
N=$(grep -c . results/v104_host_jobs.txt)
mkdir -p "$U/ANALYSIS/cibersort/logs"
sbatch -D "$U/ANALYSIS/cibersort" --array=1-"$N" --export=ALL,JOBS="$H/results/v104_host_jobs.txt" "$U/ANALYSIS/cibersort/run_v104.sbatch"
echo "host deconvolution queued: $N signatures"
