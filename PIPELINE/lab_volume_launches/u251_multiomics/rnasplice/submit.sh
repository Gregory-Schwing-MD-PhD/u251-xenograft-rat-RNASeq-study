#!/bin/bash
# Submit rnasplice: the Nextflow head job, its checks (afterany head) and the null launcher (afterany checks).
# Runs on the grid login node and does only mkdir, cp, chmod, squeue and sbatch. The folder is expected at
# $M/rnasplice (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/submit.sh
# Dependency of the head: afterany on collect_prep while it is queued (the head's own gate then checks READY and
# refs/.done/gencode). stage_refs is not a dependency: its gencode step had already passed (.done/gencode 00:46) and
# the steps still running (known sites, VEP, rnafusion bundle) are not used here; the gate still requires the marker.
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnasplice
mkdir -p "$D/launch" "$D/logs" "$D/tables" "$D/null/runs" "$D/null/results" "$M/bin" "$M/tmp/sing_cache" "$M/tmp/sing_tmp" "$M/nxf_home_nf23"
cp -f "$D/check_rnasplice.py" "$M/bin/check_rnasplice.py"
cp -f "$D/check_rnasplice.sbatch" "$M/bin/check_rnasplice.sbatch"
chmod +x "$D/run.sbatch" "$M/bin/check_rnasplice.sbatch" "$M/bin/check_rnasplice.py" "$D"/null/*.sbatch "$D"/null/*.py "$D"/null/*.R

if squeue -h -u "$USER" -n u251_rnasplice -o %i | grep -q .; then
  echo "a u251_rnasplice head job is already queued: $(squeue -h -u "$USER" -n u251_rnasplice -o '%i %T')"; exit 1
fi
PC=$(squeue -h -u "$USER" -n u251_prep_collect -o %i | head -n 1)
if [ -f "$M/reads/READY" ]; then
  DEP=""; why="READY already present"
elif [ -n "$PC" ]; then
  DEP="--dependency=afterany:$PC"; why="afterany on collect_prep $PC"
else
  DEP=""; why="no collect_prep queued and no READY: the head's gate follows the prep chain (30-min deferrals)"
fi
rm -f "$D/NEXTFLOW_EXIT" "$D/checks.json" "$D/null/inputs.env" "$D/null/JOBS" "$D/null/NULL_GATE"
H=$(sbatch --parsable -D "$D/launch" $DEP "$D/run.sbatch"); H=${H%%;*}
echo "$H" > "$D/HEAD_JOB"
K=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H" "$M/bin/check_rnasplice.sbatch"); K=${K%%;*}
echo "$K" > "$D/CHECK_JOB"
L=$(sbatch --parsable -D "$D/null" --dependency=afterany:"$K" "$D/null/null_launch.sbatch"); L=${L%%;*}
echo "$L" > "$D/null/LAUNCH_JOB"
{ echo "$(date -Is) submitted head $H ($why), checks $K (afterany:$H), null_launch $L (afterany:$K)"; } | tee -a "$D/HEAD_JOBS.log"
echo "head $H"
echo "check $K"
echo "null_launch $L"
echo "logs: $D/logs/u251_rnasplice_$H.out, u251_rnasplice_check_$K.out, u251_rnasplice_nulllaunch_$L.out"
