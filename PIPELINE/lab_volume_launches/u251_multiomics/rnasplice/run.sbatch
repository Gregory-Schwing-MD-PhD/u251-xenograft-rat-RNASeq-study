#!/bin/bash
# rnasplice head job: nf-core/rnasplice 1.0.4 (Nextflow 23.10.0) on the six graft-only tumours. Tasks go to slurm
# through grid.config. Submitted by submit.sh (sbatch -D $D/launch). The checks run in a separate job
# (check_rnasplice.sbatch, afterany on this one); the relabelling null is launched after the checks (null/).
#
# GATE. Nextflow starts only when BOTH $M/reads/READY (collect_prep wrote it: all ten prep checks passed) and
# $M/refs/.done/gencode (stage_refs' gencode step passed its md5/header/faidx checks) exist. If one is missing because
# its producer is still queued, this job does not fail and does not hold a node: it resubmits itself afterany on that
# producer (or with a 30-min --begin), at most RS_MAX_DEFER times, records the new id in $D/HEAD_JOB and exits 0. If a
# producer finished without writing its marker, the gate closes: NEXTFLOW_EXIT=gate_closed, exit 3.
#SBATCH --job-name=u251_rnasplice
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnasplice
RS_DEFER=${RS_DEFER:-0}
RS_MAX_DEFER=${RS_MAX_DEFER:-72}          # 72 x 30 min = 36 h of following the prep/refs chain, then give up
PIPE_TAG=1.0.4
PIPE_SHA=1d0494ae3402d1a46e0adadad24f81a0ff855c77   # GitHub refs/tags/1.0.4 (checked 2026-09-27)
NF_WANT=23.10.0
echo "rnasplice head: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, defer $RS_DEFER, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$D/HEAD_JOB"

defer() {   # defer <why> <dependency option>
  if [ "$RS_DEFER" -ge "$RS_MAX_DEFER" ]; then
    echo "GATE CLOSED: $1, still after $RS_DEFER deferrals"; echo "gate_closed" > "$D/NEXTFLOW_EXIT"; exit 3
  fi
  NEW=$(sbatch --parsable -D "$D/launch" $2 --export=ALL,RS_DEFER=$((RS_DEFER + 1)) "$D/run.sbatch")
  NEW=${NEW%%;*}
  if [ -z "$NEW" ]; then echo "GATE: $1 and the resubmission failed"; echo "gate_resubmit_failed" > "$D/NEXTFLOW_EXIT"; exit 5; fi
  echo "$NEW" > "$D/HEAD_JOB"
  echo "$(date -Is) deferred: $1; head resubmitted as $NEW ($2)" | tee -a "$D/HEAD_JOBS.log"
  exit 0
}

# ---- gate 1: prep_reads READY
if [ ! -f "$M/reads/READY" ]; then
  PC=$(squeue -h -u "$USER" -n u251_prep_collect -o %i 2>/dev/null | head -n 1)
  PA=$(squeue -h -u "$USER" -n u251_prep -o %i 2>/dev/null | head -n 1)
  if [ -z "$PC" ] && [ -z "$PA" ] && [ -f "$M/reads/CHECKS.tsv" ]; then
    echo "GATE CLOSED: prep_reads finished without READY. CHECKS.tsv:"; cat "$M/reads/CHECKS.tsv"
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"; exit 3
  fi
  if [ -n "$PC" ]; then defer "READY missing, collect_prep $PC queued" "--dependency=afterany:$PC"
  else defer "READY missing, no collect_prep queued (prep array: ${PA:-none})" "--begin=now+30minutes"; fi
fi
# ---- gate 2: GENCODE v44 staged and checked
if [ ! -f "$M/refs/.done/gencode" ]; then
  SR=$(squeue -h -u "$USER" -n u251_stage_refs -o %i 2>/dev/null | head -n 1)
  if [ -z "$SR" ]; then
    echo "GATE CLOSED: refs/.done/gencode missing and no stage_refs job queued"; ls -la "$M/refs/.done" 2>&1
    echo "gate_closed" > "$D/NEXTFLOW_EXIT"; exit 3
  fi
  defer "refs/.done/gencode missing, stage_refs $SR queued" "--dependency=afterany:$SR"
fi
rm -f "$D/NEXTFLOW_EXIT"
echo "gate open: $M/reads/READY ($(stat -c %y "$M/reads/READY")), $M/refs/.done/gencode ($(cat "$M/refs/.done/gencode"))"

# ---- every input exists and is not empty (nf-validation checks the sample sheet's format, not the files)
bad=0
while IFS=, read -r s f1 f2 st cond; do
  [ "$s" = "sample" ] && continue
  for f in "$f1" "$f2"; do
    if [ -s "$f" ]; then echo "input ok  $s ($cond, $st)  $(stat -c %s "$f") bytes  $f"; else echo "INPUT MISSING OR EMPTY  $s  $f"; bad=1; fi
  done
done < "$D/samplesheet.csv"
for f in "$M/refs/gencode_v44/GRCh38.primary_assembly.genome.fa" "$M/refs/gencode_v44/GRCh38.primary_assembly.annotation.gtf"; do
  if [ -s "$f" ]; then echo "reference ok  $(stat -c %s "$f") bytes  $f"; else echo "REFERENCE MISSING OR EMPTY  $f"; bad=1; fi
done
if [ "$bad" -ne 0 ]; then echo "input_missing" > "$D/NEXTFLOW_EXIT"; exit 4; fi

# ---- environment (home is at 97 % of its quota: Nextflow home, temp and container scratch all on CephFS).
# Nextflow 23.10 gets its own NXF_HOME so its plugin (nf-validation 1.1.3) and asset checkout never share a directory
# with the 25.x runs of this build.
export NXF_VER=$NF_WANT
export NXF_HOME=$M/nxf_home_nf23
export NXF_SINGULARITY_CACHEDIR=$M/singularity      # shared by the nf-core pipelines of this round
export NXF_SINGULARITY_HOME_MOUNT=true
export NXF_TEMP=$M/tmp
export TMPDIR=$M/tmp
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export NXF_OPTS='-Xms1g -Xmx6g'
export NXF_DISABLE_CHECK_LATEST=true
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH
mkdir -p "$NXF_HOME" "$NXF_SINGULARITY_CACHEDIR" "$NXF_TEMP" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR" "$D/launch" "$D/logs"
echo "nextflow launcher: $(command -v nextflow)   java: $(command -v java)   singularity: $(command -v singularity)"
java -version 2>&1 | head -n 1
NFV=$(nextflow -version 2>&1 | grep -o 'version [0-9][0-9.]*' | head -n 1)
echo "nextflow reports: $NFV"
if [ "$NFV" != "version $NF_WANT" ]; then
  echo "WRONG NEXTFLOW: wanted $NF_WANT (NXF_VER), got '$NFV'. Not running."; echo "wrong_nextflow" > "$D/NEXTFLOW_EXIT"; exit 6
fi
singularity --version || true

# ---- run (one launch directory for this run only; -resume always, because a requeued head restarts from the top)
cd "$D/launch"
nextflow run nf-core/rnasplice -r $PIPE_TAG \
  -profile singularity \
  -c ../grid.config \
  -params-file ../params.yaml \
  -w ../work \
  -resume \
  -ansi-log false
rc=$?
echo "$rc" > "$D/NEXTFLOW_EXIT"
echo "nextflow exit $rc at $(date -Is)"

# ---- record the commit that ran, and copy the pipeline's own bin/ for the null (md5-checked against the 1.0.4 tarball)
ASSET=$NXF_HOME/assets/nf-core/rnasplice
SHA=$(git -C "$ASSET" rev-parse HEAD 2>/dev/null || true)
if [ -z "$SHA" ] && [ -f "$ASSET/.git/HEAD" ]; then   # no git on the node: read the checkout (a tag checkout is detached)
  SHA=$(cat "$ASSET/.git/HEAD")
  case $SHA in "ref: "*) ref=${SHA#ref: }; SHA=$(cat "$ASSET/.git/$ref" 2>/dev/null || awk -v r="$ref" '$2 == r {print $1}' "$ASSET/.git/packed-refs") ;; esac
fi
SHA=${SHA:-unknown}
echo "$SHA" > "$D/PIPELINE_COMMIT"
if [ "$SHA" = "$PIPE_SHA" ]; then echo "pipeline commit $SHA = tag $PIPE_TAG"; else echo "PIPELINE COMMIT MISMATCH: $SHA, expected $PIPE_SHA"; fi
mkdir -p "$D/pipeline_bin"
cp -f "$ASSET"/bin/*.R "$ASSET"/bin/*.py "$D/pipeline_bin/" 2>/dev/null || true
( cd "$D/pipeline_bin" && md5sum -c "$D/null/pipeline_bin.md5" ) && echo "pipeline_bin md5 OK" || echo "PIPELINE_BIN MD5 MISMATCH (the null will refuse to run)"
exit $rc
