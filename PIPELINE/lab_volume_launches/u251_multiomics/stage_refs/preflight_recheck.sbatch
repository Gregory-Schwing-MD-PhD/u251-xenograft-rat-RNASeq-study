#!/bin/bash
# Repeat rnafusion's reference-directory check (refs_util.py preflight) from a node other than the one that wrote the
# files. exists_not_empty() in rnafusion 4.1.3 tests one file per directory, chosen by readdir order; the writing
# node's cache might list a directory in another order than a fresh client. If this node disagrees, .done/rnafusion
# is removed so the rnafusion gate stays closed. Submitted by run.sbatch with --exclude=<writer node>.
# The rnafusion launcher should run the same command on its own head node before `nextflow run`:
#   python3 /rs/rs_grp_oschome/go2432/u251_multiomics/stage_refs/refs_util.py preflight \
#     --base /rs/rs_grp_oschome/go2432/u251_multiomics/refs/rnafusion --ver 46      (exit 1 = do not launch)
#SBATCH --job-name=u251_refs_recheck
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/refs/logs/%x_%j.out
set -u
M=/rs/rs_grp_oschome/go2432/u251_multiomics
R=$M/refs
NODE=${SLURMD_NODENAME:-$(hostname -s)}
OUT=$R/.checks/rnafusion_preflight_recheck_$NODE.tsv
echo "recheck on $NODE (files written on ${WRITER_NODE:-unknown}), job ${SLURM_JOB_ID:-none}, $(date -Is)"
$HOME/mambaforge/envs/nextflow/bin/python3 "$M/stage_refs/refs_util.py" preflight --base "$R/rnafusion" --ver 46 |
  awk -v n="$NODE" 'BEGIN { FS = OFS = "\t" } { $1 = $1 "@" n; print "rnafusion_recheck", $0 }' > "$OUT"
cat "$OUT"
if awk -F'\t' '$5 == "FAIL" { f = 1 } END { exit !f }' "$OUT"; then
  rm -f "$R/.done/rnafusion"
  printf 'readdir-order check failed on %s (written on %s): see %s\n' "$NODE" "${WRITER_NODE:-unknown}" "$OUT" \
    > "$R/rnafusion/PREFLIGHT_FAILED_$NODE"
  echo "FAILED on $NODE: .done/rnafusion removed"
  exit 1
fi
echo "recheck passed on $NODE"
