#!/bin/bash
# rnavar head job: nf-core/rnavar 1.3.0 (Nextflow 25.10.7) on the graft-only U251 reads. Tasks go to slurm through
# grid.config. Submitted by submit.sh (sbatch -D $D/launch --dependency=afterany:<collect_prep>:<stage_refs>); the
# variant integration runs in a separate job (integration/integrate.sbatch, afterany on this one).
#
# GATE. Nextflow starts only when $M/reads/READY exists (collect_prep writes it only when all ten prep checks passed)
# AND $M/refs/.done/{gencode,known_sites,vep} exist (stage_refs writes each only when that step's checks passed).
# If something is missing because its producer is still queued or running, this job does not fail and does not wait on
# a node: it resubmits itself afterany on those jobs, records the new id in $D/HEAD_JOB and exits 0. If a producer
# finished and recorded a failure (prep CHECKS.tsv without READY; stage_refs STATUS.tsv with the step FAILED), the
# gate closes: NEXTFLOW_EXIT=gate_closed, exit 3. Anything else (for example stage_refs hit its time limit before
# writing STATUS.tsv) is re-checked every 30 min, at most RV_MAX_DEFER times, so that a resubmitted stage_refs is
# picked up; resubmit it with  sbatch -D $M/refs $M/stage_refs/run.sbatch  (finished steps are skipped).
#SBATCH --job-name=u251_rnavar
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnavar/logs/%x_%j.out
set -uo pipefail
umask 0002
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnavar
C=/rs/rs_grp_oschome/go2432/singularity_cache
RV_DEFER=${RV_DEFER:-0}
RV_MAX_DEFER=${RV_MAX_DEFER:-72}          # 72 x 30 min = 36 h of following the producers, then give up
NEEDED_REFS="gencode known_sites vep"
echo "rnavar head: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, defer $RV_DEFER, start $(date -Is)"
mkdir -p "$D/launch" "$D/logs"
echo "${SLURM_JOB_ID:-none}" > "$D/HEAD_JOB"

# ---- gate
missing=()
[ -f "$M/reads/READY" ] || missing+=("reads/READY")
for s in $NEEDED_REFS; do [ -f "$M/refs/.done/$s" ] || missing+=("refs/.done/$s"); done
if [ ${#missing[@]} -gt 0 ]; then
  echo "gate: missing ${missing[*]}"
  # producers still queued or running (array jobs by their base id)
  live=$(squeue -h -u "$USER" -n u251_prep,u251_prep_collect,u251_stage_refs -o %A 2>/dev/null | sort -u | paste -sd: -)
  closed=""
  if [ -z "$live" ]; then
    if [ ! -f "$M/reads/READY" ] && [ -f "$M/reads/CHECKS.tsv" ]; then
      closed="prep_reads finished without READY; CHECKS.tsv: $(tr '\n' ' ' < "$M/reads/CHECKS.tsv" | cut -c1-2000)"
    fi
    if [ -f "$M/refs/STATUS.tsv" ]; then
      for s in $NEEDED_REFS; do
        if [ ! -f "$M/refs/.done/$s" ] && awk -F'\t' -v s="$s" '$1 == s && $2 == "FAILED" { f = 1 } END { exit !f }' "$M/refs/STATUS.tsv"; then
          closed="$closed stage_refs finished with step $s FAILED (see $M/refs/.checks/$s.tsv);"
        fi
      done
    fi
  fi
  if [ -n "$closed" ]; then
    echo "GATE CLOSED: $closed"
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"
    exit 3
  fi
  if [ "$RV_DEFER" -ge "$RV_MAX_DEFER" ]; then
    echo "GATE CLOSED: ${missing[*]} still missing after $RV_DEFER deferrals"
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"
    exit 3
  fi
  if [ -n "$live" ]; then dep="--dependency=afterany:$live"; why="producers queued or running: $live"
  else dep="--begin=now+30minutes"; why="no producer queued (resubmit stage_refs if it timed out: sbatch -D $M/refs $M/stage_refs/run.sbatch)"; fi
  NEW=$(sbatch --parsable -D "$D/launch" $dep --export=ALL,RV_DEFER=$((RV_DEFER + 1)) "$D/run.sbatch")
  NEW=${NEW%%;*}
  if [ -z "$NEW" ]; then echo "GATE: inputs missing and the resubmission failed"; echo "gate_resubmit_failed" > "$D/NEXTFLOW_EXIT"; exit 5; fi
  echo "$NEW" > "$D/HEAD_JOB"
  echo "$(date -Is) deferred: ${missing[*]} missing, $why; head resubmitted as $NEW ($dep)" | tee -a "$D/HEAD_JOBS.log"
  exit 0
fi
rm -f "$D/NEXTFLOW_EXIT"
echo "gate open: $M/reads/READY ($(stat -c %y "$M/reads/READY"))"
for s in $NEEDED_REFS; do echo "  refs/.done/$s: $(cat "$M/refs/.done/$s")"; done

# ---- every input the run reads exists and is not empty (nf-schema checks existence only)
bad=0
while IFS=, read -r s f1 f2; do
  [ "$s" = "sample" ] && continue
  for f in "$f1" "$f2"; do
    if [ -s "$f" ]; then echo "input ok  $s  $(stat -c %s "$f") bytes  $f"; else echo "INPUT MISSING OR EMPTY  $s  $f"; bad=1; fi
  done
done < "$D/samplesheet.csv"
for f in $(awk -F': *' '/^(fasta|fasta_fai|gtf|dbsnp|dbsnp_tbi|known_indels|known_indels_tbi):/ { gsub(/"/, "", $2); print $2 }' "$D/params.yaml"); do
  if [ -s "$f" ]; then echo "ref ok    $(stat -c %s "$f") bytes  $f"; else echo "REFERENCE MISSING OR EMPTY  $f"; bad=1; fi
done
if [ ! -f "$M/refs/vep_cache/homo_sapiens/115_GRCh38/info.txt" ]; then echo "VEP CACHE LAYOUT MISSING: $M/refs/vep_cache/homo_sapiens/115_GRCh38/info.txt"; bad=1; fi
if [ "$bad" -ne 0 ]; then echo "input_missing" > "$D/NEXTFLOW_EXIT"; exit 4; fi

# ---- environment (home is at 97 % of its quota: Nextflow home, jar, temp and container scratch all on CephFS)
export NXF_VER=25.10.7                     # rnavar 1.3.0 manifest: nextflowVersion '!>=25.10.4'
export NXF_BIN=$M/nxf_dist/$NXF_VER/nextflow-$NXF_VER-one.jar   # the launcher downloads the jar here (not into home)
export NXF_HOME=$M/nxf_home
export NXF_SINGULARITY_CACHEDIR=$C
export NXF_SINGULARITY_HOME_MOUNT=true
export NXF_TEMP=$M/tmp
export TMPDIR=$M/tmp
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export NXF_OPTS='-Xms1g -Xmx6g'
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH   # the env's launcher script (honours NXF_VER), java 17, singularity 3.8.6
mkdir -p "$NXF_HOME" "$NXF_TEMP" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR" "$(dirname "$NXF_BIN")"
echo "nextflow launcher: $(command -v nextflow)   singularity: $(command -v singularity) ($(singularity --version 2>&1))"
v=$(nextflow -version 2>&1 | grep -io 'version [0-9.]*' | head -1)
echo "nextflow reports: $v"
case "$v" in
  *25.10.7*) : ;;
  *) echo "NEXTFLOW VERSION: expected 25.10.7, got '$v' (jar $NXF_BIN)"; echo "nextflow_version" > "$D/NEXTFLOW_EXIT"; exit 6 ;;
esac

# ---- containers the integration job needs that no pipeline pulls (pinned tags; non-fatal here, fatal there)
getimg() { # name url
  local out=$C/$1 url=$2
  [ -s "$out" ] && { echo "image present: $out"; return 0; }
  if curl --fail --location --silent --show-error --retry 5 --retry-delay 20 --connect-timeout 60 -o "$out.part.$$" "$url"; then
    mv -f "$out.part.$$" "$out" && echo "image fetched: $out ($(stat -c %s "$out") bytes)"
  else
    rm -f "$out.part.$$"; echo "image NOT fetched (the integration job retries): $url"
  fi
}
getimg depot.galaxyproject.org-singularity-pysam-0.24.1--py312hf5ad864_0.img https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0
getimg depot.galaxyproject.org-singularity-ucsc-liftover-482--h0b57e2e_0.img https://depot.galaxyproject.org/singularity/ucsc-liftover:482--h0b57e2e_0

# ---- run (one launch directory for this run only; -resume always, because a requeued head restarts from the top)
cd "$D/launch"
nextflow run nf-core/rnavar -r 1.3.0 \
  -profile singularity \
  -c ../grid.config \
  -params-file ../params.yaml \
  -w ../work \
  -resume \
  -ansi-log false
rc=$?
echo "$rc" > "$D/NEXTFLOW_EXIT"
echo "nextflow exit $rc at $(date -Is)"
git -C "$NXF_HOME/assets/nf-core/rnavar" rev-parse HEAD 2>/dev/null | sed 's/^/pipeline commit: /' | tee "$D/PIPELINE_COMMIT" || true
exit $rc
