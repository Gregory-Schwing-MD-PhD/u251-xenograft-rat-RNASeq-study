#!/bin/bash
# null_launch: afterany on the rnasplice checks. Launches the 10-split relabelling null (plus the leave-one-out of the
# true split) only if every hard check passed and check_rnasplice.py wrote null/inputs.env. It submits
#   u251_rs_null_deu   array 1-16, DEXSeq DEU (single-threaded R), sized from the pipeline's own DEXSEQ_EXON trace row
#   u251_rs_null_fast  array 1-16, DEXSeq DTU + stageR, SUPPA2 diffSplice, rMATS post; sized from the trace
#   u251_rs_null_collect  afterany on both arrays: metrics, reproduction of the pipeline, ranks, report
# and writes their ids to null/JOBS. If the checks have not run yet (the head deferred itself), it follows CHECK_JOB.
#SBATCH --job-name=u251_rnasplice_nulllaunch
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnasplice
N=$D/null
echo "null_launch: host $(hostname), job ${SLURM_JOB_ID:-none}, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$N/LAUNCH_JOB"
PY=""
for p in "$HOME/mambaforge/bin/python3" "$HOME/mambaforge/envs/spineps/bin/python" "$(command -v python3 || true)"; do
  if [ -n "$p" ] && [ -x "$p" ]; then PY=$p; break; fi
done

if [ -f "$N/JOBS" ] && grep -q '^collect ' "$N/JOBS"; then
  echo "null already launched:"; cat "$N/JOBS"; exit 0     # a requeued launcher must not submit twice
fi
if [ ! -f "$D/checks.json" ]; then
  K=$(cat "$D/CHECK_JOB" 2>/dev/null || true)
  if [ -n "$K" ] && squeue -h -j "$K" -o %i 2>/dev/null | grep -q .; then
    NEW=$(sbatch --parsable -D "$N" --dependency=afterany:"$K" "$N/null_launch.sbatch"); NEW=${NEW%%;*}
    echo "$(date -Is) checks $K not finished; null_launch resubmitted as $NEW" | tee -a "$D/HEAD_JOBS.log"
    echo "$NEW" > "$N/LAUNCH_JOB"
    exit 0
  fi
  echo "NULL GATE CLOSED: no checks.json and no live check job (CHECK_JOB=${K:-none})" | tee "$N/NULL_GATE"
  exit 3
fi
HP=$("$PY" -c 'import json,sys; print(json.load(open(sys.argv[1])).get("hard_pass"))' "$D/checks.json")
if [ "$HP" != "True" ] || [ ! -s "$N/inputs.env" ]; then
  echo "NULL GATE CLOSED: hard_pass=$HP, inputs.env $( [ -s "$N/inputs.env" ] && echo present || echo missing ). See $D/CHECKS.md" | tee "$N/NULL_GATE"
  exit 3
fi
echo "gate open: hard_pass True, inputs.env present" | tee "$N/NULL_GATE"

# size both arrays from what the same code took inside the pipeline (trace realtime and peak_rss)
read -r DEU_T DEU_M FAST_T FAST_M < <("$PY" "$N/null_size.py" "$D/results/pipeline_info")
echo "sizing: deu --time=$DEU_T --mem=$DEU_M; fast --time=$FAST_T --mem=$FAST_M (4 CPUs)"
A1=$(sbatch --parsable -D "$N" -J u251_rs_null_deu --array=1-16 --cpus-per-task=1 --mem="$DEU_M" --time="$DEU_T" \
     --export=ALL,NULL_METHODS=deu "$N/null_task.sbatch"); A1=${A1%%;*}
A2=$(sbatch --parsable -D "$N" -J u251_rs_null_fast --array=1-16 --cpus-per-task=4 --mem="$FAST_M" --time="$FAST_T" \
     --export=ALL,NULL_METHODS=dtu,suppa,rmats "$N/null_task.sbatch"); A2=${A2%%;*}
if [ -z "$A1" ] || [ -z "$A2" ]; then echo "SUBMISSION FAILED: deu '$A1' fast '$A2'"; exit 5; fi
K2=$(sbatch --parsable -D "$N" --dependency=afterany:"$A1":"$A2" "$N/null_collect.sbatch"); K2=${K2%%;*}
{ echo "deu $A1 array 1-16 time $DEU_T mem $DEU_M"
  echo "fast $A2 array 1-16 time $FAST_T mem $FAST_M"
  echo "collect $K2 afterany:$A1:$A2"
  echo "submitted $(date -Is) by null_launch ${SLURM_JOB_ID:-none}"; } > "$N/JOBS"
cat "$N/JOBS" | tee -a "$D/HEAD_JOBS.log"
