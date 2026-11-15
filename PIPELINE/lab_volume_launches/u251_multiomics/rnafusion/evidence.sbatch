#!/bin/bash
# rnafusion pass-2 evidence and the recurrence group claim (afterany on the pass-2 head). Runs only when pass 2 exited
# 0 and its reference check says OK. Writes M/rnafusion/evidence/{EVIDENCE.md,evidence_matrix.tsv,
# group_claim_splits.tsv,true_split_candidates.tsv,summary.json,EVIDENCE_OK}.
#SBATCH --job-name=u251_rf_evid
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/logs/rnafusion/%x_%j.out
set -uo pipefail
source /rs/rs_grp_oschome/go2432/u251_multiomics/rnafusion/common.sh
E=$D/evidence
mkdir -p "$E"
echo "${SLURM_JOB_ID:-none}" > "$D/EVIDENCE_JOB"
log "rnafusion evidence: host $(hostname -s), job ${SLURM_JOB_ID:-none}, follow ${RF_FOLLOW:-0}"
follow "$D/evidence.sbatch" RF_FOLLOW 40 u251_rf_p1 u251_rf_audit u251_rf_p2; fr=$?
[ "$fr" = 0 ] && exit 0
skip() { log "EVIDENCE NOT RUN: $*"; echo "$*" > "$E/EVIDENCE_SKIPPED"; exit 3; }
[ "$fr" = 2 ] && skip "follow limit reached"
[ "$(head -n 1 "$D/NEXTFLOW_EXIT_PASS2" 2>/dev/null)" = 0 ] || skip "pass 2 exit is '$(head -n 1 "$D/NEXTFLOW_EXIT_PASS2" 2>/dev/null)'"
[ "$(head -n 1 "$D/REFCHECK_PASS2" 2>/dev/null)" = OK ]    || skip "REFCHECK_PASS2 is '$(head -n 1 "$D/REFCHECK_PASS2" 2>/dev/null)'"
rm -f "$E/EVIDENCE_OK" "$E/EVIDENCE_FAILED" "$E/EVIDENCE_SKIPPED"
[ -f "$D/PASS2_RERUN_CHECK.tsv" ] && head -n 1 "$D/PASS2_RERUN_CHECK.tsv"
"$PY" "$D/evidence_rnafusion.py" --d "$D" --gene-info "$G/star/geneInfo.tab"
rc=$?
log "evidence exit $rc"
exit $rc
