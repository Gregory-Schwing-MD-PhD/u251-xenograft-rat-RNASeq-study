#!/bin/bash
# Submit human_de_null: pass A now (no dependency), pass B afterany on pass A and on the rnaseq_graft check job.
# Runs on the grid login node and does only mkdir, squeue and sbatch. The folder is expected at $M/human_de_null
# (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/human_de_null/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/human_de_null
mkdir -p "$M/logs" "$D/results"
for q in u251_de_null_A u251_de_null_B; do
  if squeue -h -u "$USER" -n "$q" -o %i | grep -q .; then
    echo "a $q job is already queued: $(squeue -h -u "$USER" -n "$q" -o '%i %T')"; exit 1
  fi
done
A=$(sbatch --parsable -D "$D" -J u251_de_null_A --export=ALL,PASS=A "$D/run.sbatch"); A=${A%%;*}
echo "$A" > "$D/PASSA_JOB"
DEP="afterany:$A"
K=$(cat "$M/rnaseq_graft/CHECK_JOB" 2>/dev/null || true)
if [ -n "$K" ] && squeue -h -j "$K" -o %i 2>/dev/null | grep -q .; then DEP="$DEP:$K"; fi
B=$(sbatch --parsable -D "$D" -J u251_de_null_B --dependency="$DEP" --export=ALL,PASS=B "$D/run.sbatch"); B=${B%%;*}
echo "$B" > "$D/PASSB_JOB"
echo "$(date -Is) submitted pass A $A (no dependency); pass B $B ($DEP)" | tee -a "$D/JOBS.log"
echo "pass A $A"
echo "pass B $B"
echo "logs: $M/logs/human_de_null_u251_de_null_A_$A.out and human_de_null_u251_de_null_B_$B.out"
