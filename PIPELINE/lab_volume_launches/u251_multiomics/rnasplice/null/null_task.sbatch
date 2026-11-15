#!/bin/bash
# null_task: one relabelling (tasks 1-10: the 10 balanced 3 v 3 splits, task 1 = the true labels) or one
# leave-one-out of the true labels (tasks 11-16), for the methods in NULL_METHODS (deu | dtu,suppa,rmats).
# Submitted by null_launch.sbatch as two arrays (-J, --time, --mem, --cpus-per-task set there).
#
# Everything is the pipeline's own code at 1.0.4 in the pipeline's own containers, on the pipeline's published
# intermediate files; only the sample-to-group assignment changes. The group that does not contain IL67B is always
# labelled "recurrent" (treatment) and the other "primary" (control), so the file names, contrast names and code
# paths are exactly those of the pipeline run; task 1 must reproduce the pipeline's tables (null_collect checks it).
#   deu    run_dexseq_exon.R on dexseq_exon/counts/<sample>.clean.count.txt + DEXSeq.gff        (bioconductor-dexseq 1.36.0)
#   dtu    counts(dmDSdata.rds) -> run_dexseq_dtu.R -> run_stager.R                           (drimseq 1.18.0, dexseq 1.36.0, stager 1.12.0)
#          (the DRIMSeq filter counts samples only, never groups, so its output is shared by every split)
#   suppa  suppa_split_file.R (r-base 3.4.2) -> suppa.py diffSplice with the pipeline's own options, no -pa (suppa 2.3)
#   rmats  rmats.py --task post on the pipeline's prep (.rmats) with regrouped --b1/--b2 and the pipeline's own options
#          (the rMATS README: "The assignment of BAMs to prep steps does not restrict the choice of --b1 and --b2 for a
#          later post step"). Fallback, only if post fails: --task both on the published BAMs, --libType fr-firststrand.
# Each method writes runs/<label>/<method>/DONE only when its outputs exist; a requeued task skips finished methods.
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnasplice/logs/%x_%A_%a.out
set -uo pipefail
M=${RS_M:-/rs/rs_grp_oschome/go2432/u251_multiomics}   # RS_M only for the local dry run
D=$M/rnasplice
N=$D/null
T=${SLURM_ARRAY_TASK_ID:?array task id}
CPUS=${SLURM_CPUS_PER_TASK:-1}
METHODS=${NULL_METHODS:?NULL_METHODS}
[ -s "$N/inputs.env" ] || { echo "no inputs.env: the checks did not pass"; exit 3; }
# shellcheck disable=SC1091
. "$N/inputs.env"
( cd "$PB" && md5sum -c --quiet "$N/pipeline_bin.md5" ) || { echo "pipeline_bin does not match the 1.0.4 scripts"; exit 3; }

row=$(awk -F'\t' -v t="$T" '$1 == t' "$N/splits.tsv")
[ -n "$row" ] || { echo "no row $T in splits.tsv"; exit 3; }
IFS=$'\t' read -r _ LABEL KIND TRT CTL DROP <<< "$row"
ORDER="IL67B_graft IL68B_graft IL69B_graft IL66B_graft NL70B_graft NL71B_graft"   # the pipeline's sample-sheet order
echo "null task $T: $LABEL ($KIND) recurrent=[$TRT] primary=[$CTL] dropped=$DROP methods=$METHODS cpus=$CPUS host $(hostname) start $(date -Is)"

export TMPDIR=$M/tmp SINGULARITY_CACHEDIR=$M/tmp/sing_cache SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export SINGULARITYENV_TMPDIR=/tmp SINGULARITYENV_PYTHONNOUSERSITE=1
export SINGULARITYENV_R_PROFILE_USER=/.Rprofile SINGULARITYENV_R_ENVIRON_USER=/.Renviron
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH
unset PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
sx() { local img=$1; shift; singularity exec --cleanenv -B /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432 -B "$M/tmp:/tmp" "$img" "$@"; }

cond_of() { case ",$TRT," in *",$1,"*) echo recurrent ;; *) case ",$CTL," in *",$1,"*) echo primary ;; *) echo none ;; esac ;; esac; }
in_task=""
for s in $ORDER; do [ "$(cond_of "$s")" != none ] && in_task="$in_task $s"; done
in_task=${in_task# }
n_in=$(wc -w <<< "$in_task")
echo "samples in this task ($n_in): $in_task"
write_sheet() {   # write_sheet <file> : sample,condition in the pipeline's order
  { echo "sample,condition"; for s in $in_task; do echo "$s,$(cond_of "$s")"; done; } > "$1"
}
write_contrasts() { printf 'contrast,treatment,control\nrecurrent_vs_primary,recurrent,primary\n' > "$1"; }

R=$N/runs/$LABEL
mkdir -p "$R"
fail=0

# ------------------------------------------------------------------ DEXSeq DEU
run_deu() {
  local W=$R/deu
  [ -f "$W/DONE" ] && { echo "deu: already done"; return 0; }
  rm -rf "$W"; mkdir -p "$W/counts"
  for s in $in_task; do ln -s "$DEU_COUNTS/$s.clean.count.txt" "$W/counts/$s.clean.count.txt"; done
  write_sheet "$W/samplesheet.csv"; write_contrasts "$W/contrasts.csv"
  ( cd "$W" && sx "$IMG_DEXSEQ" Rscript "$PB/run_dexseq_exon.R" counts "$DEU_GFF" samplesheet.csv contrasts.csv 10 ) > "$W/run.log" 2>&1
  local rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$W/perGeneQValue.recurrent_vs_primary.csv" ] && [ -s "$W/DEXSeqResults.recurrent_vs_primary.csv" ]; then
    echo "$(date -Is) rc 0" > "$W/DONE"; echo "deu: done"
  else echo "deu: FAILED rc $rc (see $W/run.log)"; tail -n 20 "$W/run.log"; fail=1; fi
}

# ------------------------------------------------------------------ DEXSeq DTU + stageR
run_dtu() {
  local W=$R/dtu
  [ -f "$W/DONE" ] && { echo "dtu: already done"; return 0; }
  rm -rf "$W"; mkdir -p "$W"
  { printf 'sample_id\tcondition\n'; for s in $in_task; do printf '%s\t%s\n' "$s" "$(cond_of "$s")"; done; } > "$W/samples.tsv"
  write_contrasts "$W/contrasts.csv"
  ( cd "$W" && sx "$IMG_DRIMSEQ" Rscript "$N/drimseq_counts.R" "$DMDS" "$(tr ' ' , <<< "$in_task")" counts.tsv \
      && sx "$IMG_DTU" Rscript "$PB/run_dexseq_dtu.R" samples.tsv contrasts.csv counts.tsv 10 \
      && sx "$IMG_STAGER" Rscript "$PB/run_stager.R" recurrent-primary DEXSeqResults.recurrent-primary.tsv perGeneQValue.recurrent-primary.tsv dexseq ) > "$W/run.log" 2>&1
  local rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$W/perGeneQValue.recurrent-primary.tsv" ] && [ -s "$W/getAdjustedPValues.recurrent-primary.tsv" ]; then
    echo "$(date -Is) rc 0" > "$W/DONE"; echo "dtu: done"
  else echo "dtu: FAILED rc $rc (see $W/run.log)"; tail -n 20 "$W/run.log"; fail=1; fi
}

# ------------------------------------------------------------------ SUPPA2 diffSplice (local events and isoforms)
run_suppa() {
  local W=$R/suppa
  [ -f "$W/DONE" ] && { echo "suppa: already done"; return 0; }
  case " $SUPPA_OPTS " in *" -pa "*) echo "suppa: REFUSED, options carry -pa: $SUPPA_OPTS"; fail=1; return 1 ;; esac
  rm -rf "$W"; mkdir -p "$W"
  write_sheet "$W/samplesheet.csv"
  ( cd "$W" \
    && sx "$IMG_RBASE" Rscript "$PB/suppa_split_file.R" "$TPM" samplesheet.csv .tpm false \
    && sx "$IMG_RBASE" Rscript "$PB/suppa_split_file.R" "$PSI_L" samplesheet.csv .psi true local \
    && sx "$IMG_RBASE" Rscript "$PB/suppa_split_file.R" "$PSI_I" samplesheet.csv .psi true transcript \
    && sx "$IMG_SUPPA" suppa.py diffSplice $SUPPA_OPTS -i "$IOE" -p local_recurrent.psi local_primary.psi \
         -e recurrent.tpm primary.tpm -o recurrent-primary_local_diffsplice \
    && sx "$IMG_SUPPA" suppa.py diffSplice $SUPPA_OPTS -i "$IOI" -p transcript_recurrent.psi transcript_primary.psi \
         -e recurrent.tpm primary.tpm -o recurrent-primary_transcript_diffsplice ) > "$W/run.log" 2>&1
  local rc=$?
  if [ "$rc" -eq 0 ] && [ -s "$W/recurrent-primary_local_diffsplice.dpsi" ] && [ -s "$W/recurrent-primary_transcript_diffsplice.dpsi" ]; then
    echo "$(date -Is) rc 0 options: $SUPPA_OPTS" > "$W/DONE"; echo "suppa: done"
  else echo "suppa: FAILED rc $rc (see $W/run.log)"; tail -n 20 "$W/run.log"; fail=1; fi
}

# ------------------------------------------------------------------ rMATS post (fallback: both)
run_rmats() {
  local W=$R/rmats
  [ -f "$W/DONE" ] && { echo "rmats: already done"; return 0; }
  case " $RMATS_OPTS " in *" --paired-stats "*) echo "rmats: REFUSED, options carry --paired-stats"; fail=1; return 1 ;; esac
  case " $RMATS_OPTS " in *" --libType fr-firststrand "*) ;; *) echo "rmats: REFUSED, options lack --libType fr-firststrand"; fail=1; return 1 ;; esac
  rm -rf "$W"; mkdir -p "$W/tmp_links"
  for f in "$RMATS_TEMP"/*.rmats; do ln -s "$f" "$W/tmp_links/$(basename "$f")"; done
  # group lists use the BAM names recorded in the prep; task 1 keeps the pipeline's own list order
  local b1="" b2="" s
  if [ "$LABEL" = split01 ]; then
    b1=$(cat "$BAMLIST_T"); b2=$(cat "$BAMLIST_C")
  else
    for s in $in_task; do
      if [ "$(cond_of "$s")" = recurrent ]; then b1="$b1,${s}_sorted.bam"; else b2="$b2,${s}_sorted.bam"; fi
    done
    b1=${b1#,}; b2=${b2#,}
  fi
  echo "$b1" > "$W/b1.txt"; echo "$b2" > "$W/b2.txt"
  local rc=1 route=post
  if [ "$KIND" = split ]; then   # rMATS's own integrity check: every BAM named once in the .rmats files (full six only)
    ( cd "$W" && sx "$IMG_RMATS" rmats.py --gtf "$GTF" --b1 b1.txt --b2 b2.txt --od inte_od --tmp tmp_links $RMATS_OPTS --task inte ) > "$W/inte.log" 2>&1
    grep -q '^Ok\.' "$W/inte.log" && echo "rmats: integrity Ok" || { echo "rmats: integrity check did not say Ok"; cat "$W/inte.log"; }
  fi
  ( cd "$W" && sx "$IMG_RMATS" rmats.py --gtf "$GTF" --b1 b1.txt --b2 b2.txt --od rmats_post --tmp tmp_links \
      $RMATS_OPTS --nthread "$CPUS" --tstat "$CPUS" --task post ) > "$W/rmats_post.log" 2>&1
  rc=$?
  if [ "$rc" -ne 0 ] || [ ! -s "$W/rmats_post/SE.MATS.JC.txt" ]; then
    echo "rmats: post failed (rc $rc); fallback --task both on the published BAMs"; tail -n 20 "$W/rmats_post.log"
    route=both
    mkdir -p "$W/bams" "$W/tmp_both"
    local l1="" l2=""
    for s in $in_task; do
      ln -sf "$BAM_DIR/${s}_sorted.bam" "$W/bams/${s}_sorted.bam"; ln -sf "$BAM_DIR/${s}_sorted.bam.bai" "$W/bams/${s}_sorted.bam.bai"
      if [ "$(cond_of "$s")" = recurrent ]; then l1="$l1,$W/bams/${s}_sorted.bam"; else l2="$l2,$W/bams/${s}_sorted.bam"; fi
    done
    echo "${l1#,}" > "$W/b1_both.txt"; echo "${l2#,}" > "$W/b2_both.txt"
    ( cd "$W" && sx "$IMG_RMATS" rmats.py --gtf "$GTF" --b1 b1_both.txt --b2 b2_both.txt --od rmats_both --tmp tmp_both \
        $RMATS_OPTS --nthread "$CPUS" --tstat "$CPUS" --task both ) > "$W/rmats_both.log" 2>&1
    rc=$?
    [ "$rc" -eq 0 ] && [ -s "$W/rmats_both/SE.MATS.JC.txt" ] && ln -sfn rmats_both "$W/rmats_out"
  else
    ln -sfn rmats_post "$W/rmats_out"
  fi
  if [ "$rc" -eq 0 ] && [ -s "$W/rmats_out/SE.MATS.JC.txt" ]; then
    echo "$route" > "$W/ROUTE"; echo "$(date -Is) rc 0 route $route options: $RMATS_OPTS" > "$W/DONE"; echo "rmats: done ($route)"
  else echo "rmats: FAILED rc $rc"; fail=1; fi
}

for m in ${METHODS//,/ }; do
  case $m in
    deu) run_deu ;;
    dtu) run_dtu ;;
    suppa) run_suppa ;;
    rmats) run_rmats ;;
    *) echo "unknown method $m"; fail=1 ;;
  esac
done
echo "null task $T ($LABEL) finished $(date -Is), fail=$fail"
exit $fail
