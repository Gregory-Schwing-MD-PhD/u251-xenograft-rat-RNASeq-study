#!/bin/bash
# Submit rnaseq_graft: the Nextflow head job (afterany on collect_prep when it is queued) and the checks (afterany on
# the head). Runs on the grid login node and does only mkdir, cp, squeue and sbatch. The folder is expected at
# $M/rnaseq_graft (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/rnaseq_graft/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnaseq_graft
mkdir -p "$D/launch" "$D/logs" "$D/tables" "$M/bin" "$M/tmp/sing_cache" "$M/tmp/sing_tmp" "$M/nxf_home"
cp -f "$D/check_rnaseq_graft.py" "$M/bin/check_rnaseq_graft.py"
cp -f "$D/check_rnaseq_graft.sbatch" "$M/bin/check_rnaseq_graft.sbatch"
chmod +x "$D/run.sbatch" "$M/bin/check_rnaseq_graft.sbatch" "$M/bin/check_rnaseq_graft.py"

if squeue -h -u "$USER" -n u251_rnaseq_graft -o %i | grep -q .; then
  echo "a u251_rnaseq_graft head job is already queued: $(squeue -h -u "$USER" -n u251_rnaseq_graft -o '%i %T')"; exit 1
fi

PC=$(squeue -h -u "$USER" -n u251_prep_collect -o %i | head -n 1)
if [ -f "$M/reads/READY" ]; then
  DEP=""; why="READY already present"
elif [ -n "$PC" ]; then
  DEP="--dependency=afterany:$PC"; why="afterany on collect_prep $PC"
else
  DEP=""; why="no collect_prep queued and no READY: the head's gate follows the prep chain (30-min deferrals)"
fi
rm -f "$D/NEXTFLOW_EXIT"
H=$(sbatch --parsable -D "$D/launch" $DEP "$D/run.sbatch"); H=${H%%;*}
echo "$H" > "$D/HEAD_JOB"
K=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H" "$M/bin/check_rnaseq_graft.sbatch"); K=${K%%;*}
echo "$K" > "$D/CHECK_JOB"
{ echo "$(date -Is) submitted head $H ($why), checks $K (afterany:$H)"; } | tee -a "$D/HEAD_JOBS.log"
echo "head $H"
echo "check $K"
echo "logs: $D/logs/u251_rnaseq_graft_$H.out and u251_rnaseq_graft_check_$K.out"
