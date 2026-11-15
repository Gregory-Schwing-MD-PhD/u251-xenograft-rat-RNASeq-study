#!/bin/bash
# Submit integration_rna_prep. Runs on the grid login node and does only mkdir, cp, cat, squeue and sbatch.
# The folder is expected at $M/integration_rna_prep (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/integration_rna_prep/submit.sh
# Two jobs:
#   u251_irp        IRP_MODE=published, no dependency (needs nothing from any other arm)
#   u251_irp_graft  IRP_MODE=graft, afterany on the rnaseq_graft check job (if queued) and on u251_irp; its first step
#                   refuses to run unless rnaseq_graft/checks.json says all_hard_pass, and re-queues itself behind a
#                   still-queued rnaseq_graft job instead of waiting
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
C=$M/integration_rna_prep
I=$M/integration
mkdir -p "$I/bin" "$I/logs" "$I/refs" "$I/results" "$I/tmp"
for f in run.sbatch prep_rna_refs.R irp_json.R gtf_extract.awk locus_count.awk; do cp -f "$C/$f" "$I/bin/$f"; done
chmod +x "$I/bin/run.sbatch"

if squeue -h -u "$USER" -n u251_irp,u251_irp_graft -o %i | grep -q .; then
  echo "a u251_irp job is already queued: $(squeue -h -u "$USER" -n u251_irp,u251_irp_graft -o '%i %j %T' | tr '\n' ' ')"; exit 1
fi

A=$(sbatch --parsable -D "$I" --job-name=u251_irp --export=ALL,IRP_MODE=published "$I/bin/run.sbatch"); A=${A%%;*}

CK=$(cat "$M/rnaseq_graft/CHECK_JOB" 2>/dev/null || true)
if [ -n "$CK" ] && squeue -h -j "$CK" -o %i 2>/dev/null | grep -q .; then
  DEP="--dependency=afterany:$CK:$A"; why="afterany on rnaseq_graft check $CK and u251_irp $A"
else
  DEP="--dependency=afterany:$A"; why="afterany on u251_irp $A (rnaseq_graft check job not queued; the gate decides)"
fi
G=$(sbatch --parsable -D "$I" $DEP --job-name=u251_irp_graft --export=ALL,IRP_MODE=graft "$I/bin/run.sbatch"); G=${G%%;*}

echo "$A" > "$I/P0_JOB"
echo "$G" > "$I/P0_GRAFT_JOB"
echo "$(date -Is) u251_irp $A (no dependency); u251_irp_graft $G ($why)" | tee -a "$I/JOBS.log"
echo "logs: $I/logs/u251_irp_$A.out and $I/logs/u251_irp_graft_$G.out"
