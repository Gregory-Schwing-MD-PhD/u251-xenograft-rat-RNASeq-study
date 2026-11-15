#!/bin/bash
# rnaseq_graft head job: nf-core/rnaseq 3.22.2 on the graft-only U251 reads. Tasks go to slurm through grid.config.
# Submitted by submit.sh (sbatch -D $D/launch --dependency=afterany:<collect_prep>); the checks run in a separate job
# (check_rnaseq_graft.sbatch, afterany on this one).
#
# GATE. Nextflow starts only when $M/reads/READY exists (collect_prep writes it only when all ten prep checks passed).
# If READY is missing because prep_reads has not finished (or not been submitted yet), this job does not fail and does
# not wait on a node: it resubmits itself (afterany on collect_prep when that job is queued, otherwise a 30-min
# --begin), at most RG_MAX_DEFER times, records the new id in $D/HEAD_JOB and exits 0. If prep_reads has finished
# without READY (CHECKS.tsv present, no prep job queued), the gate closes: NEXTFLOW_EXIT=gate_closed, exit 3.
#SBATCH --job-name=u251_rnaseq_graft
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnaseq_graft/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnaseq_graft
RG_DEFER=${RG_DEFER:-0}
RG_MAX_DEFER=${RG_MAX_DEFER:-72}          # 72 x 30 min = 36 h of following the prep chain, then give up
echo "rnaseq_graft head: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, defer $RG_DEFER, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$D/HEAD_JOB"

# ---- gate
if [ ! -f "$M/reads/READY" ]; then
  PC=$(squeue -h -u "$USER" -n u251_prep_collect -o %i 2>/dev/null | head -n 1)
  PA=$(squeue -h -u "$USER" -n u251_prep -o %i 2>/dev/null | head -n 1)
  if [ -z "$PC" ] && [ -z "$PA" ] && [ -f "$M/reads/CHECKS.tsv" ]; then
    echo "GATE CLOSED: prep_reads finished without READY. CHECKS.tsv:"; cat "$M/reads/CHECKS.tsv"
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"
    exit 3
  fi
  if [ "$RG_DEFER" -ge "$RG_MAX_DEFER" ]; then
    echo "GATE CLOSED: READY still missing after $RG_DEFER deferrals (prep_reads never finished)."
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"
    exit 3
  fi
  if [ -n "$PC" ]; then dep="--dependency=afterany:$PC"; why="collect_prep $PC is queued"
  else dep="--begin=now+30minutes"; why="no collect_prep queued yet (prep array: ${PA:-none})"; fi
  NEW=$(sbatch --parsable -D "$D/launch" $dep --export=ALL,RG_DEFER=$((RG_DEFER + 1)) "$D/run.sbatch")
  NEW=${NEW%%;*}
  if [ -z "$NEW" ]; then echo "GATE: READY missing and the resubmission failed"; echo "gate_resubmit_failed" > "$D/NEXTFLOW_EXIT"; exit 5; fi
  echo "$NEW" > "$D/HEAD_JOB"
  echo "$(date -Is) deferred: READY missing, $why; head resubmitted as $NEW ($dep)" | tee -a "$D/HEAD_JOBS.log"
  exit 0
fi
rm -f "$D/NEXTFLOW_EXIT"
echo "gate open: $M/reads/READY ($(stat -c %y "$M/reads/READY"))"
cat "$M/reads/READY"
[ -f "$M/reads/CHECKS.tsv" ] && cat "$M/reads/CHECKS.tsv"

# ---- every FASTQ in the samplesheet exists and is not empty (nf-schema checks existence only)
bad=0
while IFS=, read -r s f1 f2 st; do
  [ "$s" = "sample" ] && continue
  for f in "$f1" "$f2"; do
    if [ -s "$f" ]; then echo "input ok  $s  $(stat -c %s "$f") bytes  $f"; else echo "INPUT MISSING OR EMPTY  $s  $f"; bad=1; fi
  done
done < "$D/samplesheet.csv"
if [ "$bad" -ne 0 ]; then echo "input_missing" > "$D/NEXTFLOW_EXIT"; exit 4; fi

# ---- environment (home is at 97 % of its quota: Nextflow home, temp and container scratch all on CephFS)
export NXF_VER=25.10.2
export NXF_HOME=$M/nxf_home
export NXF_SINGULARITY_CACHEDIR=/rs/rs_grp_oschome/go2432/singularity_cache
export NXF_SINGULARITY_HOME_MOUNT=true
export NXF_TEMP=$M/tmp
export TMPDIR=$M/tmp
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export NXF_OPTS='-Xms1g -Xmx6g'
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH
mkdir -p "$NXF_HOME" "$NXF_TEMP" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR" "$D/launch" "$D/logs"
echo "nextflow: $(command -v nextflow)   singularity: $(command -v singularity)"
nextflow -version 2>&1 | grep -i version || true
singularity --version || true

# ---- run (one launch directory for this run only; -resume always, because a requeued head restarts from the top)
cd "$D/launch"
nextflow run nf-core/rnaseq -r 3.22.2 \
  -profile singularity \
  -c ../grid.config \
  -params-file ../params.yaml \
  -w ../work \
  -resume \
  -ansi-log false
rc=$?
echo "$rc" > "$D/NEXTFLOW_EXIT"
echo "nextflow exit $rc at $(date -Is)"
# the commit that ran (the published run was v3.22.2-g3816d48; check 0 compares the versions files)
git -C "$NXF_HOME/assets/nf-core/rnaseq" rev-parse HEAD 2>/dev/null | sed 's/^/pipeline commit: /' || true
exit $rc
