#!/bin/bash
# Checks for rnaseq_graft. Submitted by submit.sh with --dependency=afterany:<head>. Installed as
# $M/bin/check_rnaseq_graft.sbatch next to $M/bin/check_rnaseq_graft.py.
# If the head job deferred itself (READY was missing, see run.sbatch), NEXTFLOW_EXIT does not exist yet and the
# current head id is in $D/HEAD_JOB: this job then resubmits itself afterany on that id and exits 0.
# Writes $D/CHECKS.md, $D/checks.json and $D/tables/. Deletes $D/work only when Nextflow exited 0 AND every hard check
# passed (CLEAN_WORK=0 keeps it).
#SBATCH --job-name=u251_rnaseq_graft_check
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnaseq_graft/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnaseq_graft
CLEAN_WORK=${CLEAN_WORK:-1}
echo "rnaseq_graft checks: host $(hostname), job ${SLURM_JOB_ID:-none}, start $(date -Is)"

if [ ! -f "$D/NEXTFLOW_EXIT" ]; then
  H=$(cat "$D/HEAD_JOB" 2>/dev/null || true)
  if [ -n "$H" ] && [ "$H" != "${SLURM_JOB_ID:-}" ] && squeue -h -j "$H" -o %i 2>/dev/null | grep -q .; then
    NEW=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H" --export=ALL,CLEAN_WORK="$CLEAN_WORK" "$M/bin/check_rnaseq_graft.sbatch")
    NEW=${NEW%%;*}
    echo "$(date -Is) head $H still queued/running (no NEXTFLOW_EXIT yet); checks resubmitted as $NEW" | tee -a "$D/HEAD_JOBS.log"
    echo "$NEW" > "$D/CHECK_JOB"
    exit 0
  fi
  echo "no NEXTFLOW_EXIT and no live head job (HEAD_JOB=${H:-none}): the checks will report what is there"
fi

PY=""
for p in "$HOME/mambaforge/bin/python3" "$HOME/mambaforge/envs/spineps/bin/python" "$(command -v python3 || true)"; do
  if [ -n "$p" ] && [ -x "$p" ]; then PY=$p; break; fi
done
echo "python: $PY ($("$PY" --version 2>&1))"
"$PY" "$M/bin/check_rnaseq_graft.py" --d "$D" --genesets "$D/genesets" > "$D/logs/check_stdout_${SLURM_JOB_ID:-x}.md"
rc=$?
echo "check exit $rc (0 = all hard checks pass)"
sed -n '1,12p' "$D/CHECKS.md" 2>/dev/null || true

nfx=$(cat "$D/NEXTFLOW_EXIT" 2>/dev/null || echo missing)
if [ "$rc" -eq 0 ] && [ "$nfx" = "0" ] && [ "$CLEAN_WORK" = "1" ] && [ -d "$D/work" ]; then
  sz=$(getfattr --only-values -n ceph.dir.rbytes "$D/work" 2>/dev/null || echo unknown)
  echo "$(date -Is) all hard checks passed and Nextflow exited 0: deleting $D/work ($sz bytes)"
  rm -rf "$D/work" && echo "work deleted $(date -Is)" | tee -a "$D/CHECKS.md"
else
  echo "work kept (check rc $rc, nextflow exit $nfx, CLEAN_WORK $CLEAN_WORK)"
fi
exit $rc
