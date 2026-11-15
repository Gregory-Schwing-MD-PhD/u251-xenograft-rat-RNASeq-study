#!/bin/bash
# human_de_null: does the published human primary-vs-recurrent DE (nf-core/differentialabundance 1.5.0, 3 v 3; 19,351
# genes tested, 82 in the filtered table) and its gene-set result stand out among the 10 balanced 3-v-3 relabellings of
# the six tumours? No Nextflow: two R steps in two containers, run by `singularity exec`. One script for both passes:
#   PASS=A  the published graft+both counts; the analysis of record. R0 must reproduce the published table, or the job
#           stops with exit 10 before any other fit and nothing from this item is used.
#   PASS=B  graft-only counts from rnaseq_graft; the pre-registered rat-leak sensitivity analysis. Gated on rnaseq_graft's
#           checks.json (checks 0-4 and 10 passed, written after a Nextflow exit 0) and on pass A's R0. While the
#           rnaseq_graft check job or pass A is still queued or running, this job resubmits itself afterany on it and
#           exits 0 (at most DN_MAX_DEFER times): the rnaseq_graft chain resubmits its own jobs, so a dependency on an
#           old id can release this job before the checks exist.
# Submitted by submit.sh; by hand:
#   sbatch -D $D -J u251_de_null_A --export=ALL,PASS=A $D/run.sbatch
#   sbatch -D $D -J u251_de_null_B --dependency=afterany:<pass A>:<rnaseq_graft check> --export=ALL,PASS=B $D/run.sbatch
# Requeue-safe: every fit and every gene-set run is written as .part and renamed; a restart reads what exists.
#SBATCH --job-name=u251_de_null
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/logs/human_de_null_%x_%j.out
set -uo pipefail
PASS=${PASS:-A}
case "$PASS" in A|B) ;; *) echo "PASS must be A or B, got '$PASS'"; exit 2 ;; esac
DN_DEFER=${DN_DEFER:-0}
DN_MAX_DEFER=${DN_MAX_DEFER:-100}
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/human_de_null
RG=$M/rnaseq_graft
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
BIN=$D/bin
RES=$D/results
OUT=$RES/pass$PASS
PY=$HOME/mambaforge/bin/python3
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
NCPU=${SLURM_CPUS_PER_TASK:-4}
mkdir -p "$OUT" "$M/logs"
echo "human_de_null pass $PASS: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, defer $DN_DEFER, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$D/PASS${PASS}_JOB"
# the item spec calls this directory M/de_null; the link lets either name resolve
[ -e "$M/de_null" ] || ln -s human_de_null "$M/de_null" 2>/dev/null || true
status() { echo "$1" > "$OUT/STATUS.part" && mv -f "$OUT/STATUS.part" "$OUT/STATUS"; echo "STATUS $1 ($(date -Is))"; }

# ---- params.yaml (flat key: value)
p() {
  awk -v k="$1" '{ l = $0; sub(/[ \t]+#.*$/, "", l) }
                 index(l, k ":") == 1 { v = substr(l, length(k) + 2); gsub(/^[ \t]+|[ \t]+$/, "", v); print v; exit }' \
      "$D/params.yaml"
}
for k in passA_counts passA_lengths passB_counts passB_lengths passB_required_checks metadata published_table \
         published_filtered published_sizefactors published_params published_design_params broad_gsea_dir gmt_dir \
         deseq2_image expect_deseq2 bioconductor_image rlib filtering_min_abundance filtering_min_samples deseq2_alpha \
         padj_cut differential_min_fold_change lfc_shrink r0_tolerance rp_block_set gsea_nperm gsea_set_min \
         gsea_set_max gsea_seed; do
  v=$(p "$k")
  if [ -z "$v" ]; then echo "params.yaml: $k is missing"; status params_missing; exit 2; fi
  printf -v "P_$k" '%s' "$v"
done
if [ "$PASS" = A ]; then
  COUNTS=$P_passA_counts; LENGTHS=$P_passA_lengths; R0=enforce; CMP=none
else
  COUNTS=$P_passB_counts; LENGTHS=$P_passB_lengths; R0=skip; CMP=$RES/passA
fi
ASHR=no; [ "$P_lfc_shrink" = ashr ] && ASHR=yes

# ---- pass B: follow the rnaseq_graft chain and pass A, then gate
live() { [ -n "${1:-}" ] && [ "$1" != "${SLURM_JOB_ID:-}" ] && squeue -h -j "$1" -o %i 2>/dev/null | grep -q .; }
if [ "$PASS" = B ]; then
  K=$(cat "$RG/CHECK_JOB" 2>/dev/null || true)
  AJ=$(cat "$D/PASSA_JOB" 2>/dev/null || true)
  WAIT=""
  if live "$K"; then WAIT=$K; fi
  if live "$AJ"; then WAIT=${WAIT:+$WAIT:}$AJ; fi
  if [ -n "$WAIT" ]; then
    if [ "$DN_DEFER" -ge "$DN_MAX_DEFER" ]; then
      echo "GATE CLOSED: still waiting on $WAIT after $DN_DEFER deferrals"; status gate_closed_max_deferrals; exit 3
    fi
    NEW=$(sbatch --parsable -D "$D" -J u251_de_null_B --dependency=afterany:"$WAIT" \
            --export=ALL,PASS=B,DN_DEFER=$((DN_DEFER + 1)) "$D/run.sbatch")
    NEW=${NEW%%;*}
    if [ -z "$NEW" ]; then echo "GATE: resubmission failed"; status gate_resubmit_failed; exit 5; fi
    echo "$NEW" > "$D/PASSB_JOB"
    echo "$(date -Is) pass B deferred (#$((DN_DEFER + 1))): waiting on $WAIT; resubmitted as $NEW" | tee -a "$D/JOBS.log"
    status "deferred_to_$NEW"
    exit 0
  fi
  "$PY" "$BIN/gate_passb.py" --checks "$RG/checks.json" --nextflow-exit "$RG/NEXTFLOW_EXIT" \
      --required "$P_passB_required_checks" --r0 "$RES/passA/r0_check.json" --counts "$COUNTS" --lengths "$LENGTHS" \
      --samplesheet "$D/samplesheet.csv" --out "$OUT/gate.json"
  if [ $? -ne 0 ]; then echo "GATE CLOSED (pass B): reasons in $OUT/gate.json"; status gate_closed; exit 3; fi
fi

# ---- inputs present and not empty
for f in "$COUNTS" "$LENGTHS" "$P_metadata" "$P_published_table" "$P_published_filtered" "$P_gmt_dir/combined_human.gmt" \
         "$P_deseq2_image" "$P_bioconductor_image" "$BIN/human_de_null_de.R" "$BIN/human_de_null_gsea.R" "$BIN/hv_json.R" \
         "$D/samplesheet.csv"; do
  if [ ! -s "$f" ]; then echo "INPUT MISSING OR EMPTY: $f"; status input_missing; exit 4; fi
done
rm -f "$OUT/STATUS"
status running

# ---- md5 of every input; sha256 of both images (go2432-bioconductor.sif is not digest-pinned)
h() {
  if [ -f "$2" ]; then printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$(stat -c %s "$2")" "$(md5sum < "$2" | cut -d' ' -f1)"
  else printf '%s\t%s\tNA\tmissing\n' "$1" "$2"; fi
}
{
  printf 'role\tpath\tbytes\tmd5\n'
  h counts "$COUNTS"; h lengths "$LENGTHS"; h metadata "$P_metadata"; h samplesheet "$D/samplesheet.csv"
  h params "$D/params.yaml"; h published_table "$P_published_table"; h published_filtered "$P_published_filtered"
  h published_sizefactors "$P_published_sizefactors"; h published_params "$P_published_params"
  h published_design_params "$P_published_design_params"
  for g in combined_human.gmt hs.hallmark.v2023.2.gmt hs.kegg.v2023.2.gmt hs.c2_curated.v2023.2.gmt \
           hs.gobp.v2023.2.gmt hs.brain_gmt_v2.gmt; do h "gmt:$g" "$P_gmt_dir/$g"; done
  for r in Recurrent_U2 Primary_U2; do h "broad_report:$r" "$P_broad_gsea_dir/therapy_impact.combined_human.gsea_report_for_$r.tsv"; done
  for s in human_de_null_de.R human_de_null_gsea.R hv_json.R gate_passb.py; do h "script:$s" "$BIN/$s"; done
  h "script:run.sbatch" "$D/run.sbatch"
  if [ "$PASS" = B ]; then h rnaseq_graft_checks "$RG/checks.json"; h passA_r0_check "$RES/passA/r0_check.json"; fi
} > "$OUT/inputs_md5.tsv.part" && mv -f "$OUT/inputs_md5.tsv.part" "$OUT/inputs_md5.tsv"
{
  printf 'role\tpath\tbytes\tsha256\n'
  for kv in "deseq2_image|$P_deseq2_image" "bioconductor_image|$P_bioconductor_image"; do
    f=${kv#*|}; printf '%s\t%s\t%s\t%s\n' "${kv%%|*}" "$f" "$(stat -c %s "$f")" "$(sha256sum < "$f" | cut -d' ' -f1)"
  done
  printf 'singularity\t%s\tNA\t%s\n' "$SING" "$("$SING" --version 2>&1 | head -n 1)"
} > "$OUT/containers.tsv.part" && mv -f "$OUT/containers.tsv.part" "$OUT/containers.tsv"
cat "$OUT/containers.tsv"

# ---- container environment: temp and runtime dirs in /dev/shm (node /tmp is swept), no R library from home
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH JAVA_HOME
SHM=/dev/shm/${USER}_${SLURM_JOB_ID:-$$}
export XDG_RUNTIME_DIR=$SHM/xdg
export TMPDIR=$SHM/tmp
mkdir -p "$XDG_RUNTIME_DIR" "$TMPDIR"
export SINGULARITYENV_TMPDIR=$TMPDIR
export SINGULARITYENV_R_LIBS_USER=$SHM/no_user_library
trap 'rm -rf "$SHM"' EXIT
BND=(--bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432)

echo "== parse both R scripts in the DESeq2 container, $(date -Is)"
"$SING" exec "${BND[@]}" "$P_deseq2_image" Rscript -e \
  "for (f in c('$BIN/human_de_null_de.R', '$BIN/human_de_null_gsea.R')) invisible(parse(f)); cat('parse ok\n')"
if [ $? -ne 0 ]; then status parse_failed; exit 6; fi

echo "== step 1: DESeq2 fits in the published pipeline's container, $(date -Is)"
"$SING" exec "${BND[@]}" "$P_deseq2_image" Rscript "$BIN/human_de_null_de.R" \
  out="$OUT" pass="$PASS" samplesheet="$D/samplesheet.csv" metadata="$P_metadata" \
  counts="$COUNTS" lengths="$LENGTHS" gmt="$P_gmt_dir/combined_human.gmt" \
  pub="$P_published_table" pub_filtered="$P_published_filtered" pub_sf="$P_published_sizefactors" \
  r0="$R0" cores="$NCPU" ashr="$ASHR" expect_deseq2="$P_expect_deseq2" \
  min_abundance="$P_filtering_min_abundance" min_samples="$P_filtering_min_samples" alpha="$P_deseq2_alpha" \
  padj_cut="$P_padj_cut" min_fc="$P_differential_min_fold_change" r0_tol="$P_r0_tolerance" rp_set="$P_rp_block_set"
rc=$?
if [ "$rc" -eq 10 ]; then
  echo "R0 FAILED: the published table was not reproduced; see $OUT/r0_check.json. Nothing else from this item is used."
  status r0_failed; exit 10
fi
if [ "$rc" -ne 0 ]; then status "de_failed_$rc"; exit "$rc"; fi

echo "== step 2: fgsea in the Bioconductor image, $(date -Is)"
"$SING" exec "${BND[@]}" "$P_bioconductor_image" Rscript "$BIN/human_de_null_gsea.R" \
  out="$OUT" pass="$PASS" gmt_dir="$P_gmt_dir" rlib="$P_rlib" cores="$NCPU" broad_dir="$P_broad_gsea_dir" \
  inputs_md5="$OUT/inputs_md5.tsv" containers="$OUT/containers.tsv" compare="$CMP" \
  nperm="$P_gsea_nperm" set_min="$P_gsea_set_min" set_max="$P_gsea_set_max" seed="$P_gsea_seed" padj_cut="$P_padj_cut"
rc=$?
if [ "$rc" -ne 0 ]; then status "gsea_failed_$rc"; exit "$rc"; fi
status done
ls -la "$OUT"
echo "done $(date -Is)"
