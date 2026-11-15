#!/bin/bash
# Checks for rnasplice. Submitted by submit.sh with --dependency=afterany:<head>. Installed as
# $M/bin/check_rnasplice.sbatch next to $M/bin/check_rnasplice.py.
# If the head job deferred itself (a gate marker was missing, see run.sbatch), NEXTFLOW_EXIT does not exist yet and
# the current head id is in $D/HEAD_JOB: this job then resubmits itself afterany on that id and exits 0.
# Writes $D/CHECKS.md, $D/checks.json and, only when every hard check passed, $D/null/inputs.env.
# It never deletes the work directory: null_collect does that after the null has reproduced the pipeline.
#SBATCH --job-name=u251_rnasplice_check
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnasplice
echo "rnasplice checks: host $(hostname), job ${SLURM_JOB_ID:-none}, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$D/CHECK_JOB"

if [ ! -f "$D/NEXTFLOW_EXIT" ]; then
  H=$(cat "$D/HEAD_JOB" 2>/dev/null || true)
  if [ -n "$H" ] && [ "$H" != "${SLURM_JOB_ID:-}" ] && squeue -h -j "$H" -o %i 2>/dev/null | grep -q .; then
    NEW=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H" "$M/bin/check_rnasplice.sbatch")
    NEW=${NEW%%;*}
    echo "$(date -Is) head $H still queued/running (no NEXTFLOW_EXIT yet); checks resubmitted as $NEW" | tee -a "$D/HEAD_JOBS.log"
    echo "$NEW" > "$D/CHECK_JOB"
    exit 0
  fi
  echo "no NEXTFLOW_EXIT and no live head job (HEAD_JOB=${H:-none}): the checks will report what is there"
fi

export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH      # singularity 3.8.6 for the samtools strand count
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache SINGULARITY_TMPDIR=$M/tmp/sing_tmp TMPDIR=$M/tmp
unset PYTHONPATH
PY=""
for p in "$HOME/mambaforge/bin/python3" "$HOME/mambaforge/envs/spineps/bin/python" "$(command -v python3 || true)"; do
  if [ -n "$p" ] && [ -x "$p" ]; then PY=$p; break; fi
done
echo "python: $PY ($("$PY" --version 2>&1))"
rm -f "$D/null/inputs.env"
"$PY" "$M/bin/check_rnasplice.py" --d "$D"
rc=$?
echo "check exit $rc (0 = all hard checks pass)"
sed -n '1,8p' "$D/CHECKS.md" 2>/dev/null || true
exit $rc
