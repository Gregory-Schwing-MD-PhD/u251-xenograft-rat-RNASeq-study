#!/bin/bash
# Submit prep_reads: the 10-task array, then collect_prep afterany on it. Runs on the grid login node and does only
# mkdir, cp and sbatch. The folder is expected at $M/prep_reads (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/prep_reads/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
mkdir -p "$M/reads/logs" "$M/tmp/prep_reads" "$M/bin"
cp -f "$M/prep_reads/run.sbatch" "$M/bin/prep_reads_run.sbatch"
cp -f "$M/prep_reads/collect_prep.sbatch" "$M/bin/collect_prep.sbatch"
A=$(sbatch --parsable -D "$M/reads" "$M/bin/prep_reads_run.sbatch")
A=${A%%;*}
C=$(sbatch --parsable -D "$M/reads" --dependency=afterany:"$A" --export=ALL,PREP_ARRAY="$A" "$M/bin/collect_prep.sbatch")
C=${C%%;*}
echo "prep_reads array $A (tasks 1-10: C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B NL71B)"
echo "collect_prep $C (afterany:$A)"
echo "logs: $M/reads/logs/u251_prep_${A}_<task>.out and u251_prep_collect_${C}.out"
