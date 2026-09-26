#!/bin/bash
# Rerun of the host deconvolution after the Bowman 2016 task (40473363_3) stopped on a degenerate permutation draw.
# Builds the tissue-complete signature (build_combined_signature.py: Zhang 2014 non-myeloid + Bowman 2016 myeloid), then
# queues four CIBERSORT v1.04 tasks through ../cibersort/run_v104_guarded.R, 1,000 permutations each:
#   1 Bowman 2016 as built by host_deconv_prep.py (the failed task)   2 the combined signature
#   3 Zhang 2014 and 4 LM22 again: the guarded wrapper must reproduce 40473363_2 and _1 (fractions, r, RMSE and P).
# W is a working directory on CephFS holding copies of the four scripts under W/scripts (the grid clone is not touched;
# its host mixture, signatures, references and vendor/CIBERSORT.R are read in place).
#   sbatch -D $W --export=ALL,W=$W $W/scripts/run_host_deconv_fix.sbatch
#SBATCH --job-name=u251_hostfix
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/fix_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=$U/ANALYSIS/host
W=${W:?set W to the working directory}
PY=/rs/rs_grp_oschome/go2432/envs/u251_drug/bin/python
unset PYTHONPATH LD_LIBRARY_PATH
R=$H/refs_online
MIX=$H/results/mixture_host_hs.txt
mkdir -p "$W/results" "$W/logs"
cd "$W/scripts"
"$PY" build_combined_signature.py --zhang "$R/zhang2014_brain_mean7_tpm.txt" --bowman "$R/bowman2016_myeloid_mean4.txt" \
  --mouse2human "$R/mouse_human_1to1.tsv" --mixture "$MIX" --out "$W/results"
{
  printf 'host_bowman2016_myeloid_mean4\t%s\t1000\t%s\n' "$H/results/signature_host_bowman2016_myeloid_mean4.txt" "$MIX"
  printf 'host_zhang6_bowman4\t%s\t1000\t%s\n' "$W/results/signature_host_zhang6_bowman4.txt" "$MIX"
  printf 'host_zhang2014_brain_mean7_tpm\t%s\t1000\t%s\n' "$H/results/signature_host_zhang2014_brain_mean7_tpm.txt" "$MIX"
  printf 'host_lm22\t%s\t1000\t%s\n' "$H/results/signature_host_lm22.txt" "$MIX"
} > "$W/results/jobs.txt"
N=$(grep -c . "$W/results/jobs.txt")
sbatch -D "$W" --array=1-"$N" \
  --export=ALL,JOBS="$W/results/jobs.txt",OUT="$W/results",SCRIPTS="$W/scripts",CIBERSORT_R="$U/ANALYSIS/cibersort/vendor/CIBERSORT.R" \
  "$W/scripts/run_v104_guarded.sbatch"
echo "queued $N guarded CIBERSORT tasks"
