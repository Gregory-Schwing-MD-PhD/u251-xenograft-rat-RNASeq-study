#!/bin/bash
# Submit rnavar: the -profile test job (no dependency, never blocks anything), the Nextflow head job (afterany on
# collect_prep and stage_refs when they are still queued; its own gate checks READY and the refs .done markers) and the
# variant integration (afterany on the head). Runs on the grid login node and does only mkdir, chmod, squeue and sbatch.
# The folder is expected at $M/rnavar (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/rnavar/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnavar
mkdir -p "$D/launch" "$D/logs" "$D/test/launch" "$D/integration/work" "$M/tmp/sing_cache" "$M/tmp/sing_tmp" \
         "$M/nxf_home" "$M/nxf_dist"
chmod +x "$D/run.sbatch" "$D/run_test.sbatch" "$D/integration/integrate.sbatch" "$D/integration/rnavar_integrate.py"

if squeue -h -u "$USER" -n u251_rnavar -o %i | grep -q .; then
  echo "a u251_rnavar head job is already queued: $(squeue -h -u "$USER" -n u251_rnavar -o '%i %T')"; exit 1
fi

# producers still queued or running: collect_prep (writes READY) and stage_refs (writes refs/.done/*)
deps=$(squeue -h -u "$USER" -n u251_prep,u251_prep_collect,u251_stage_refs -o %A | sort -u | paste -sd: -)
if [ -n "$deps" ]; then DEP="--dependency=afterany:$deps"; why="afterany on $deps"
else DEP=""; why="no producer queued (the head's gate checks READY and refs/.done itself)"; fi

rm -f "$D/NEXTFLOW_EXIT" "$D/INTEGRATION_EXIT"
TJ=$(sbatch --parsable -D "$D/test/launch" "$D/run_test.sbatch"); TJ=${TJ%%;*}
echo "$TJ" > "$D/TEST_JOB"
HJ=$(sbatch --parsable -D "$D/launch" $DEP "$D/run.sbatch"); HJ=${HJ%%;*}
echo "$HJ" > "$D/HEAD_JOB"
IJ=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$HJ" "$D/integration/integrate.sbatch"); IJ=${IJ%%;*}
echo "$IJ" > "$D/INTEG_JOB"
echo "$(date -Is) submitted test $TJ, head $HJ ($why), integration $IJ (afterany:$HJ)" | tee -a "$D/HEAD_JOBS.log"
echo "test $TJ"
echo "head $HJ"
echo "integration $IJ"
echo "logs: $D/logs/u251_rnavar_$HJ.out, u251_rnavar_integ_$IJ.out, u251_rnavar_test_$TJ.out"
