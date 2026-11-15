#!/bin/bash
# invasion_niche: U251 in the contralateral hemisphere (Y dosage), the state of those cells, published infiltration
# programmes by arm (10-split rank), and the host niche (Mmp13, marker panels). Custom R on the outputs of
# nf-core/rnaseq 3.22.2 (graft-only run, M/rnaseq_graft) and the existing host run (H). No Nextflow.
#
# GATE. Runs only when rnaseq_graft's checks.json exists, is newer than its NEXTFLOW_EXIT, and shows checks 1, 2, 3 and
# 10 passed. While the rnaseq_graft head or its check job is still queued or running (both defer themselves while
# prep_reads runs), this job resubmits itself afterany on the live check job (or a 30-min --begin) and exits 0, at
# most IN_MAX_DEFER times; the live id is kept in $D/JOB. If rnaseq_graft finished and the checks failed, the gate
# closes: results/GATE_CLOSED, exit 3.
#
#   submitted by submit.sh:  sbatch -D $D --dependency=afterany:<rnaseq_graft check job> $D/run.sbatch
#SBATCH --job-name=u251_invasion_niche
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/logs/%x_%j.out
set -uo pipefail
M=/rs/rs_grp_oschome/go2432/u251_multiomics
D=$M/invasion_niche
RG=$M/rnaseq_graft
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=/rs/rs_grp_oschome/go2432/u251_host
OUT=$D/results
IN_DEFER=${IN_DEFER:-0}
IN_MAX_DEFER=${IN_MAX_DEFER:-60}
source "$D/grid.config"
echo "invasion_niche: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, defer $IN_DEFER, start $(date -Is)"
echo "${SLURM_JOB_ID:-none}" > "$D/JOB"

## ---------------------------------------------------------------------------------------------------- gate
live() { [ -n "${1:-}" ] && squeue -h -j "$1" -o %i 2>/dev/null | grep -q . ; }
CK=$(cat "$RG/CHECK_JOB" 2>/dev/null || true)
HD=$(cat "$RG/HEAD_JOB" 2>/dev/null || true)
if live "$CK" || live "$HD"; then
  if [ "$IN_DEFER" -ge "$IN_MAX_DEFER" ]; then
    echo "GATE CLOSED: rnaseq_graft still running after $IN_DEFER deferrals"; mkdir -p "$OUT"; echo "rnaseq_graft still running after $IN_DEFER deferrals" > "$OUT/GATE_CLOSED"; exit 3
  fi
  if live "$CK"; then dep="--dependency=afterany:$CK"; why="rnaseq_graft check job $CK is live"
  else dep="--begin=now+30minutes"; why="rnaseq_graft head $HD is live, no live check job"; fi
  NEW=$(sbatch --parsable -D "$D" $dep --export=ALL,IN_DEFER=$((IN_DEFER + 1)) "$D/run.sbatch"); NEW=${NEW%%;*}
  if [ -z "$NEW" ]; then echo "GATE: resubmission failed"; exit 5; fi
  echo "$NEW" > "$D/JOB"
  echo "$(date -Is) deferred ($why); resubmitted as $NEW ($dep)" | tee -a "$D/JOBS.log"
  exit 0
fi
PY=""
for p in "$HOME/mambaforge/bin/python3" "$HOME/mambaforge/envs/spineps/bin/python" "$(command -v python3 || true)"; do
  if [ -n "$p" ] && [ -x "$p" ]; then PY=$p; break; fi
done
mkdir -p "$OUT/inputs" "$D/logs"
rm -f "$OUT/GATE_CLOSED"
if [ ! -s "$RG/checks.json" ] || { [ -f "$RG/NEXTFLOW_EXIT" ] && [ "$RG/checks.json" -ot "$RG/NEXTFLOW_EXIT" ]; }; then
  echo "GATE CLOSED: $RG/checks.json missing or older than NEXTFLOW_EXIT, and no rnaseq_graft job is live"
  echo "checks.json missing or stale; no rnaseq_graft job live (HEAD_JOB=${HD:-none}, CHECK_JOB=${CK:-none})" > "$OUT/GATE_CLOSED"; exit 3
fi
GATE=$("$PY" - "$RG/checks.json" <<'EOF'
import json, sys
c = json.load(open(sys.argv[1]))
ch = c.get("checks", {})
need = ["1", "2", "3", "10"]
res = {k: ch.get(k, {}).get("pass") for k in need}
ok = all(v is True for v in res.values())
print(("OPEN" if ok else "CLOSED") + " " + " ".join("check%s=%s" % (k, res[k]) for k in need)
      + " all_hard_pass=%s failed_hard=%s nextflow_exit=%s finished=%s" % (c.get("all_hard_pass"), c.get("failed_hard"),
                                                                          c.get("nextflow_exit"), c.get("finished")))
EOF
)
echo "gate: $GATE"
echo "rnaseq_graft checks.json ($(stat -c %y "$RG/checks.json")): $GATE" > "$OUT/inputs/gate.txt"
case "$GATE" in
  OPEN*) ;;
  *) echo "GATE CLOSED: checks 1-3 and 10 of rnaseq_graft did not all pass"; cp "$OUT/inputs/gate.txt" "$OUT/GATE_CLOSED"; exit 3 ;;
esac
if [ ! -s "$RG/results/star_salmon/salmon.merged.gene_counts.tsv" ]; then
  echo "GATE CLOSED: graft-only count table missing"; echo "input missing: graft-only salmon.merged.gene_counts.tsv" > "$OUT/GATE_CLOSED"; exit 3
fi
# part B alone needs N269B's Picard metrics; without them part B fails and says so, the other parts run
[ -s "$RG/results/star_salmon/picard_metrics/N269B_graft.markdup.sorted.MarkDuplicates.metrics.txt" ] || echo "WARNING: N269B_graft Picard metrics missing; part B will not run"

## ---------------------------------------------------------------------------------------------------- environment
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH JAVA_HOME
JOBTAG=${SLURM_JOB_ID:-$$}
export SINGULARITY_TMPDIR=$SHM_BASE/${USER}_inv_sing_$JOBTAG
export TMPDIR=$SHM_BASE/${USER}_inv_tmp_$JOBTAG
export XDG_RUNTIME_DIR=$SHM_BASE/${USER}_inv_xdg_$JOBTAG
export SINGULARITYENV_TMPDIR=$TMPDIR
mkdir -p "$SINGULARITY_TMPDIR" "$TMPDIR" "$XDG_RUNTIME_DIR"
trap 'rm -rf "$SINGULARITY_TMPDIR" "$TMPDIR" "$XDG_RUNTIME_DIR"' EXIT
img() { if [ -f "$SING_CACHE/$1" ]; then echo "$SING_CACHE/$1"; elif [ -f "$SING_CACHE_LEGACY/$1" ]; then echo "$SING_CACHE_LEGACY/$1"; fi; }
DESEQ_IMG=$(img "$DESEQ_IMG_NAME"); SIF=$(img "$SIF_NAME")
if [ -z "$DESEQ_IMG" ] || [ -z "$SIF" ]; then echo "image missing (DESeq2: ${DESEQ_IMG:-none}; sif: ${SIF:-none})"; exit 4; fi
{ printf 'role\tpath\tbytes\tsha256\n'
  for x in "stepA_report:$DESEQ_IMG" "stepC:$SIF"; do r=${x%%:*}; p=${x#*:}
    printf '%s\t%s\t%s\t%s\n' "$r" "$p" "$(stat -c %s "$p")" "$(sha256sum "$p" | cut -d' ' -f1)"; done
} > "$OUT/inputs/images.tsv"
cat "$OUT/inputs/images.tsv"
B=(--bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432)
"$SING" --version

## ---------------------------------------------------------------------------------------------------- references (idempotent)
bash "$D/bin/gtf_gene_table.sh" "$U/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" "$M/refs/gene_tables/gencode_v44_genes.tsv" || { echo "human gene table failed"; exit 6; }
grep -q "version 44 (Ensembl 110)" "$M/refs/gene_tables/gencode_v44_genes.tsv.header" || { echo "human GTF header is not GENCODE v44 (Ensembl 110)"; exit 6; }
bash "$D/bin/gtf_gene_table.sh" "$H/results_rnaseq/genome/Rattus_norvegicus.mRatBN7.2.110.filtered.gtf" "$M/refs/gene_tables/mRatBN7.2_110_genes.tsv" || { echo "rat gene table failed"; exit 6; }
bash "$D/bin/stage_genesets.sh" "$U/ANALYSIS/cibersort/signatures_online/ivygap2018_anatomic4.txt" "$U/ANALYSIS/spread/signatures" \
  "$D/genesets_src" "$U/ANALYSIS/refs/pathways/human/hs.hallmark.v2023.2.gmt" "$U/ANALYSIS/refs/pathways/human/hs.c2_curated.v2023.2.gmt" \
  "$OUT/genesets" || { echo "gene-set staging failed"; exit 6; }
{ printf 'file\tbytes\tmd5\n'
  for f in "$RG/results/star_salmon/salmon.merged.gene_counts.tsv" "$RG/results/star_salmon/picard_metrics/N269B_graft.markdup.sorted.MarkDuplicates.metrics.txt" \
           "$RG/checks.json" "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
           "$H/results_rnaseq/star_salmon/salmon.merged.gene_counts.tsv" "$H/verify/results/fits.tsv" "$H/verify/results/r1_mmp13.tsv" \
           "$H/verify/results/r2_relabel_splits.tsv" "$M/refs/gene_tables/gencode_v44_genes.tsv" "$M/refs/gene_tables/mRatBN7.2_110_genes.tsv" \
           "$U"/ANALYSIS/xengsort_out/*.txt "$D/params.yaml"; do
    printf '%s\t%s\t%s\n' "$f" "$(stat -c %s "$f" 2>/dev/null || echo NA)" "$(md5sum "$f" 2>/dev/null | cut -d' ' -f1)"; done
  for f in "$M/refs/gene_tables/gencode_v44_genes.tsv" "$M/refs/gene_tables/mRatBN7.2_110_genes.tsv"; do
    printf '%s (source GTF md5)\tNA\t%s\n' "$f" "$(cat "$f.src_md5")"; done
} > "$OUT/inputs/inputs_md5.tsv"
for d in "$M/de_null/results/passA" "$M/de_null/results/passB"; do echo "human_de_null fits $d: $([ -f "$d/fits.tsv" ] && echo present || echo absent)"; done

## ---------------------------------------------------------------------------------------------------- steps
cd "$D"
echo "== step A (parts A, B, D; vst and fits) in $DESEQ_IMG, $(date -Is)"
"$SING" exec "${B[@]}" "$DESEQ_IMG" Rscript "$D/bin/inv_stepA.R" "$D/params.yaml" > "$D/logs/stepA_$JOBTAG.log" 2>&1
rcA=$?; echo "step A exit $rcA"; tail -n 12 "$D/logs/stepA_$JOBTAG.log"
rcC=9
if [ -f "$OUT/C_programmes/expr/vst_graft_only_all10.tsv.gz" ]; then
  echo "== step C (part C: GSVA, fgsea) in $SIF, $(date -Is)"
  "$SING" exec "${B[@]}" "$SIF" Rscript "$D/bin/inv_stepC.R" "$D/params.yaml" > "$D/logs/stepC_$JOBTAG.log" 2>&1
  rcC=$?; echo "step C exit $rcC"; tail -n 12 "$D/logs/stepC_$JOBTAG.log"
else
  echo "step C skipped: step A wrote no vst matrix"
fi
echo "== report, $(date -Is)"
"$SING" exec "${B[@]}" "$DESEQ_IMG" Rscript "$D/bin/inv_report.R" "$D/params.yaml" > "$D/logs/report_$JOBTAG.log" 2>&1
rcR=$?; echo "report exit $rcR"
cp -f "$D/logs/stepA_$JOBTAG.log" "$D/logs/stepC_$JOBTAG.log" "$OUT/inputs/" 2>/dev/null || true
printf 'stepA\t%s\nstepC\t%s\nreport\t%s\nfinished\t%s\n' "$rcA" "$rcC" "$rcR" "$(date -Is)" > "$OUT/EXIT_CODES.tsv"
sed -n '1,40p' "$OUT/REPORT.md" 2>/dev/null || true
if [ "$rcA" -eq 0 ] && [ "$rcC" -eq 0 ] && [ "$rcR" -eq 0 ]; then echo "invasion_niche done $(date -Is)"; exit 0; fi
echo "invasion_niche finished with failures (A $rcA, C $rcC, report $rcR) $(date -Is)"; exit 2
