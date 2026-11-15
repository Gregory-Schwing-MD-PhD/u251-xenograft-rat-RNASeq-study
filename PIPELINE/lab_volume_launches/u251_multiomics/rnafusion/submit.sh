#!/bin/bash
# Submit the rnafusion chain: head pass 1 -> audit -> head pass 2 -> evidence, each afterany on the one before and each
# gated on marker files (see README.md). Runs on the grid login node and does only mkdir, chmod, squeue and sbatch.
# The folder is expected at /rs/rs_grp_oschome/go2432/u251_multiomics/rnafusion (copied there from the local repo).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/rnafusion/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/rnafusion
mkdir -p "$D/launch" "$D/checks" "$D/audit" "$D/evidence" "$M/logs/rnafusion" "$M/singularity" \
         "$M/tmp/sing_cache" "$M/tmp/sing_tmp" "$M/nxf_home"
chmod +x "$D/run.sbatch" "$D/audit.sbatch" "$D/evidence.sbatch" "$D/audit_rnafusion.py" "$D/evidence_rnafusion.py"

for n in u251_rf_p1 u251_rf_audit u251_rf_p2 u251_rf_evid; do
  if squeue -h -u "$USER" -n "$n" -o %i | grep -q .; then
    echo "a $n job is already queued: $(squeue -h -u "$USER" -n "$n" -o '%i %T')"; exit 1
  fi
done

# pass 1 waits for whatever of prep_reads / stage_refs / the refs recheck is still queued; its gate re-checks markers
UP=$(for n in u251_prep u251_prep_collect u251_stage_refs u251_refs_recheck; do squeue -h -u "$USER" -n "$n" -o %A; done |
     sort -u | paste -sd: -)
DEP=""
[ -n "$UP" ] && DEP="--dependency=afterany:$UP"
H1=$(sbatch --parsable -D "$D/launch" $DEP --export=ALL,RF_PASS=1 "$D/run.sbatch"); H1=${H1%%;*}
AU=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H1" "$D/audit.sbatch"); AU=${AU%%;*}
H2=$(sbatch --parsable -J u251_rf_p2 -D "$D/launch" --dependency=afterany:"$AU" --export=ALL,RF_PASS=2 "$D/run.sbatch"); H2=${H2%%;*}
EV=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$H2" "$D/evidence.sbatch"); EV=${EV%%;*}
echo "$(date -Is) submitted: head pass 1 $H1 (${DEP:-no dependency}), audit $AU (afterany:$H1), head pass 2 $H2" \
     "(afterany:$AU), evidence $EV (afterany:$H2)" | tee -a "$D/JOBS.log"
echo "head1 $H1"
echo "audit $AU"
echo "head2 $H2"
echo "evidence $EV"
echo "logs: $M/logs/rnafusion/"
