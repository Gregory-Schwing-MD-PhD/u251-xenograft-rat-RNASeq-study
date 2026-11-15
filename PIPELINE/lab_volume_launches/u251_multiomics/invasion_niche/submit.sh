#!/bin/bash
# Submit invasion_niche. Runs on the grid login node and does only mkdir, ln, cat, squeue and sbatch. The folder is
# expected at $M/invasion_niche (copied from the local repository with scp), with the fetched gene-set supplements in
# $M/invasion_niche/genesets_src (copied from the session scratchpad; never in the repository).
#   bash /rs/rs_grp_oschome/go2432/u251_multiomics/invasion_niche/submit.sh
set -eu
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/invasion_niche
RG=$M/rnaseq_graft
mkdir -p "$D/logs" "$D/results" "$M/logs" "$M/refs/gene_tables" "$M/singularity"
# SELECTION.md names this item's directory M/invasion; the build uses M/invasion_niche. Both resolve to one place.
[ -e "$M/invasion" ] || ln -s invasion_niche "$M/invasion"
chmod +x "$D/run.sbatch" "$D/bin/"*.sh
for f in SOURCES.tsv neftel2019_tableS2_metamodules.tsv venkataramani2022_tableS2a_invasivity.tsv convert_supp.py; do
  [ -s "$D/genesets_src/$f" ] || { echo "missing $D/genesets_src/$f"; exit 1; }
done
if squeue -h -u "$USER" -n u251_invasion_niche -o %i | grep -q .; then
  echo "a u251_invasion_niche job is already queued: $(squeue -h -u "$USER" -n u251_invasion_niche -o '%i %T')"; exit 1
fi
CK=$(cat "$RG/CHECK_JOB" 2>/dev/null || true)
if [ -n "$CK" ] && squeue -h -j "$CK" -o %i 2>/dev/null | grep -q .; then
  DEP="--dependency=afterany:$CK"; why="afterany on rnaseq_graft checks $CK"
else
  DEP=""; why="no live rnaseq_graft check job: the gate reads checks.json now or follows the chain"
fi
J=$(sbatch --parsable -D "$D" $DEP "$D/run.sbatch"); J=${J%%;*}
echo "$J" > "$D/JOB"
echo "$(date -Is) submitted $J ($why)" | tee -a "$D/JOBS.log"
echo "job $J"
echo "log: $M/logs/u251_invasion_niche_$J.out; results: $D/results/REPORT.md"
