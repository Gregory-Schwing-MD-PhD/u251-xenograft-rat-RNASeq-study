# Sourced by every rnafusion job (run.sbatch, audit.sbatch, evidence.sbatch). Paths, pins and two helpers.
# Nothing here writes; the helpers only read squeue and call sbatch.

M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnafusion                                   # this item: configs, launch/, work/, results/, audit/, evidence/
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=/rs/rs_grp_oschome/go2432/u251_host
REFS=$M/refs/rnafusion                           # staged by stage_refs (S3 mirror of nf-core-awsmegatests)
G=$REFS/GRCh38/gencode_v46
RESULTS=$D/results
LOGS=$M/logs/rnafusion
PY=$HOME/mambaforge/envs/nextflow/bin/python3    # the same interpreter stage_refs' recheck uses
LIBS="C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B NL71B"
RF_TAG=4.1.3
RF_TAG_SHA=76ad76e7c39b2ba9edc35aa3602e3dc454d842ec   # GitHub refs/tags/4.1.3 (lightweight tag -> this commit), read 2026-09-27
RF_NXF_VER=25.10.7                                   # rnafusion 4.1.3 manifest: nextflowVersion '!>=25.10.4'

log() { echo "$(date -Is) $*"; }

# upstream_ids NAME...  ->  colon-joined array/job ids (%A) of this user's jobs with those names still in squeue
upstream_ids() {
  local n
  for n in "$@"; do squeue -h -u "$USER" -n "$n" -o %A 2>/dev/null; done | sort -u | paste -sd: -
}

# follow SCRIPT COUNTVAR MAX NAME...
#   If any job named NAME... is still queued or running, resubmit SCRIPT (same job name, same environment) with
#   afterany on those jobs, append the new id to $D/JOBS.log and return 0 (the caller exits 0).
#   Returns 1 when nothing upstream is queued (the caller then decides: gate open or gate closed).
#   Returns 2 when the follow limit is reached (the caller closes its gate).
follow() {
  local script=$1 cvar=$2 max=$3; shift 3
  local ids; ids=$(upstream_ids "$@")
  [ -z "$ids" ] && return 1
  local c=${!cvar:-0}
  if [ "$c" -ge "$max" ]; then log "FOLLOW LIMIT reached ($cvar=$c >= $max) while waiting on: $*"; return 2; fi
  local new
  new=$(sbatch --parsable -J "${SLURM_JOB_NAME:-u251_rf}" -D "$D/launch" --dependency=afterany:"$ids" \
          --export=ALL,"$cvar"=$((c + 1)) "$script")
  new=${new%%;*}
  if [ -z "$new" ]; then log "FOLLOW: resubmission of $script failed"; return 2; fi
  log "follow: upstream still queued ($* -> $ids); resubmitted $script as $new (afterany:$ids, $cvar=$((c + 1)))" |
    tee -a "$D/JOBS.log"
  echo "$new" > "$D/$(basename "$script" .sbatch)_${RF_PASS:-x}.jobid"
  return 0
}
