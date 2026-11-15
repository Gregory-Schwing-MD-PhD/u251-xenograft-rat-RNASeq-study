#!/bin/bash
# integration_rna_prep: the RNA half of the U251N methylation integration (SELECTION S8). README.md in this folder
# says what each step answers and why every setting is what it is. Reads nothing from the methylation run
# (/rs/rs_grp_oschome/go2432/u251_meth) or from ANALYSIS/methylation/.
#   IRP_MODE=published (default)  items 1-6 on the published human run (graft+both counts, markdup BAMs)
#   IRP_MODE=graft                items 1, 2, 5 and 6 again on the rnaseq_graft counts (sensitivity analysis); gated
#                                 on rnaseq_graft/checks.json all_hard_pass, and defers itself behind a queued
#                                 rnaseq_graft job (at most graft_max_deferrals times) instead of waiting
# Submitted by submit.sh as $I/bin/run.sbatch with sbatch -D $I. Idempotent: a requeue redoes everything except the
# reference downloads whose md5 already checks out.
# Exit: 0 complete (non-fatal failed checks are listed in summary.json); 2 fatal check; 3 graft gate not met;
#       4 setup error; 5 identity STOP (C2B not nearest to the U-251 family); 1 unexpected R error.
#SBATCH --job-name=u251_irp
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/integration/logs/%x_%j.out
set -euo pipefail
export LC_ALL=C

MODE=${IRP_MODE:-published}
CONF=${IRP_CONF:-/rs/rs_grp_oschome/go2432/u251_multiomics/integration_rna_prep}
PF=$CONF/params.yaml
SS=$CONF/samplesheet.csv
ts()  { date '+%F %T'; }
log() { echo "[$(ts)] $*"; }
die() { echo "[$(ts)] SETUP ERROR: $*"; exit 4; }
pv()  { awk -v k="$1" '{ sub(/[ \t]+#.*$/, "") } index($0, k ":") == 1 { sub(/^[^:]*:[ \t]*/, ""); print; f = 1; exit } END { if (!f) exit 1 }' "$PF"; }

[ -r "$PF" ] && [ -r "$SS" ] || die "params.yaml or samplesheet.csv missing in $CONF"
case "$MODE" in published|graft) ;; *) die "IRP_MODE must be published or graft, got $MODE" ;; esac
I=$(pv I) || die "I not in params"
M=$(pv M) || die "M not in params"
U=$(pv U) || die "U not in params"
REFS=$(pv refs_dir) || die "refs_dir not in params"
SING=$(pv singularity) || die "singularity not in params"
BIN=$I/bin
OUT=$I/results/p0; [ "$MODE" = graft ] && OUT=$I/results/p0_graft
WD=$I/tmp/p0_$MODE
mkdir -p "$I/logs" "$REFS" "$I/results" "$I/tmp" "$WD" "$OUT"
export TMPDIR=$I/tmp
export SINGULARITYENV_TMPDIR=/tmp
export SINGULARITY_CACHEDIR=$I/tmp/sing_cache SINGULARITY_TMPDIR=$I/tmp/sing_tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR"
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID:-$$}
mkdir -p "$XDG_RUNTIME_DIR"
trap 'rm -rf "$XDG_RUNTIME_DIR"' EXIT
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH JAVA_HOME
log "integration_rna_prep mode=$MODE job=${SLURM_JOB_ID:-none} node=$(hostname) conf=$CONF out=$OUT"

# ---------------------------------------------------------------- graft gate (sensitivity rerun only)
if [ "$MODE" = graft ]; then
  CJ=$(pv graft_checks); CG=$(pv counts_graft); TG=$(pv tpm_graft)
  if ! { [ -s "$CJ" ] && grep -Eq '"all_hard_pass": *true' "$CJ" && [ -s "$CG" ] && [ -s "$TG" ]; }; then
    n=${IRP_DEFER:-0}; maxd=$(pv graft_max_deferrals)
    pend=$(squeue -h -u "$USER" -n u251_rnaseq_graft_check,u251_rnaseq_graft -o %i 2>/dev/null | head -n 1 || true)
    if [ -n "$pend" ] && [ "$n" -lt "$maxd" ]; then
      J=$(sbatch --parsable -D "$I" --dependency=afterany:"$pend" --job-name=u251_irp_graft \
            --export=ALL,IRP_MODE=graft,IRP_DEFER=$((n + 1)) "$BIN/run.sbatch")
      log "GATE: $CJ missing or all_hard_pass not true; deferred behind rnaseq_graft job $pend as job ${J%%;*} (deferral $((n + 1)) of $maxd)"
      exit 0
    fi
    log "GATE NOT MET: $CJ missing or all_hard_pass not true, and no rnaseq_graft job queued (or $maxd deferrals used). Nothing run."
    exit 3
  fi
  log "gate ok: $CJ all_hard_pass true"
fi

# ---------------------------------------------------------------- containers (pinned tags; md5 recorded)
find_img() {   # name url -> path; looks in img_dirs, pulls into the first one if absent
  local name=$1 url=$2 d first=""
  IFS=, read -ra DIRS <<< "$(pv img_dirs)"
  for d in "${DIRS[@]}"; do
    [ -z "$first" ] && first=$d
    if [ -s "$d/$name" ]; then echo "$d/$name"; return 0; fi
  done
  mkdir -p "$first"
  "$SING" pull "$first/$name" "$url" >&2 || return 1
  echo "$first/$name"
}
IMG_R=$(find_img "$(pv img_r)" "$(pv img_r_url)") || die "R image unavailable"
IMG_S=$(find_img "$(pv img_samtools)" "$(pv img_samtools_url)") || die "samtools image unavailable"
BIND=(-B "$U:$U" -B "$M:$M" -B "$I/tmp:/tmp")
RSCRIPT=("$SING" exec "${BIND[@]}" "$IMG_R" Rscript --vanilla)
SAMTOOLS=("$SING" exec "${BIND[@]}" "$IMG_S" samtools)
sv=$("$SING" --version)
stv=$("${SAMTOOLS[@]}" --version); stv=${stv%%$'\n'*}
rv=$("${RSCRIPT[@]}" -e 'cat(R.version.string)')
{
  printf 'singularity\t%s\n' "$sv"
  printf 'img_r\t%s\n' "$IMG_R"; printf 'img_r_md5\t%s\n' "$(md5sum < "$IMG_R" | cut -d' ' -f1)"; printf 'img_r_R\t%s\n' "$rv"
  printf 'img_samtools\t%s\n' "$IMG_S"; printf 'img_samtools_md5\t%s\n' "$(md5sum < "$IMG_S" | cut -d' ' -f1)"
  printf 'samtools_version\t%s\n' "$stv"
  printf 'slurm_job\t%s\n' "${SLURM_JOB_ID:-none}"; printf 'node\t%s\n' "$(hostname)"; printf 'conf_dir\t%s\n' "$CONF"
} > "$WD/env.tsv"
log "R image $IMG_R ($rv); samtools image $IMG_S ($stv); $sv"
"${RSCRIPT[@]}" -e "invisible(parse('$BIN/prep_rna_refs.R')); invisible(parse('$BIN/irp_json.R')); cat('R parse ok\n')"

# ---------------------------------------------------------------- md5 of the ten BAMs, in the background
MD5PID=""
if [ "$MODE" = published ]; then
  rm -f "$WD/bam_md5.txt"
  awk -F, 'NR > 1 { print $3 }' "$SS" | xargs -P 4 -I{} md5sum {} > "$WD/bam_md5.txt.part" &
  MD5PID=$!
fi

# ---------------------------------------------------------------- item 4: references (GET only; md5 asserted or recorded)
fetch_ref() {   # name url want -> "status md5 how"
  local name=$1 url=$2 want=$3 f=$REFS/$1 got
  if [ -s "$f" ]; then
    got=$(md5sum < "$f" | cut -d' ' -f1)
    if [ "$want" = record ]; then
      if [ -s "$f.md5" ] && [ "$got" = "$(cut -d' ' -f1 < "$f.md5")" ]; then echo "ok $got kept"; return 0; fi
    elif [ "$got" = "$want" ]; then echo "ok $got kept"; return 0; fi
  fi
  rm -f "$f.part"
  if ! curl -sSL --fail --retry 5 --retry-delay 15 --connect-timeout 60 -o "$f.part" "$url"; then
    rm -f "$f.part"; echo "download_failed NA none"; return 0
  fi
  got=$(md5sum < "$f.part" | cut -d' ' -f1)
  if [ "$want" != record ] && [ "$got" != "$want" ]; then mv -f "$f.part" "$f.md5_mismatch"; echo "md5_mismatch $got downloaded"; return 0; fi
  mv -f "$f.part" "$f"; echo "$got  $name" > "$f.md5"; echo "ok $got downloaded"
}
printf 'name\turl\texpected_md5\tmd5\tbytes\tstatus\thow\n' > "$REFS/refs_status.tsv.part"
for k in $(grep -oE '^ref_[0-9]+' "$PF"); do
  line=$(pv "$k")
  read -r name url want <<< "$line"
  res=$(fetch_ref "$name" "$url" "$want")
  read -r st md how <<< "$res"
  bytes=$(stat -c %s "$REFS/$name" 2>/dev/null || echo NA)
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$url" "$want" "$md" "$bytes" "$st" "$how" >> "$REFS/refs_status.tsv.part"
  log "ref $name: $st $md ($how, $bytes bytes)"
done
mv -f "$REFS/refs_status.tsv.part" "$REFS/refs_status.tsv"

# ---------------------------------------------------------------- GENCODE v44 gene table + chr9/chr10 exons
GTF=$(pv gtf)
gzip -cd "$GTF" | awk -v G="$WD/genes.tsv" -v E="$WD/exons_chr9_chr10.tsv" -v H="$WD/gtf_header.txt" -f "$BIN/gtf_extract.awk"
log "GTF extracted: $(($(wc -l < "$WD/genes.tsv") - 1)) genes"

# ---------------------------------------------------------------- item 3: MGMT, CDKN2A/CDKN2B and controls in the BAMs
if [ "$MODE" = published ]; then
  "${RSCRIPT[@]}" "$BIN/prep_rna_refs.R" loci "$PF" "$SS" "$WD" < /dev/null
  # stdin of every command inside the two read loops comes from /dev/null, or singularity would eat the samplesheet
  MAPQ=$(pv mapq_unique); XF=$(pv sam_exclude_flags)
  : > "$WD/locus_counts.tsv.part"; : > "$WD/locus_nm.tsv"
  while IFS=$'\t' read -r gname gidv chr gs ge gstr role rest; do
    [ "$gname" = gene_name ] && continue
    while IFS=, read -r lib xlog bam slog pic; do
      [ "$lib" = sample ] && continue
      [ -s "$bam" ] && [ -s "$bam.bai" ] || die "BAM or index missing: $bam"
      "${SAMTOOLS[@]}" view -q "$MAPQ" -F "$XF" "$bam" "$chr:$gs-$ge" < /dev/null \
        | awk -v EX="$WD/exons_$gname.txt" -v GS="$gs" -v GE="$ge" -v GSTR="$gstr" -v LIB="$lib" -v LOC="$gname" \
              -v NMF="$WD/locus_nm.tsv" -f "$BIN/locus_count.awk" >> "$WD/locus_counts.tsv.part"
    done < "$SS"
    log "locus $gname ($role, $chr:$gs-$ge $gstr) counted in all libraries"
  done < "$WD/loci.tsv"
  mv -f "$WD/locus_counts.tsv.part" "$WD/locus_counts.tsv"
  wait "$MD5PID"
  mv -f "$WD/bam_md5.txt.part" "$WD/bam_md5.txt"
  log "BAM md5 done ($(wc -l < "$WD/bam_md5.txt") files)"
fi

# ---------------------------------------------------------------- items 1, 2, 5, 6 (+ 3 and 4 summaries)
set +e
"${RSCRIPT[@]}" "$BIN/prep_rna_refs.R" main "$PF" "$SS" "$WD" "$OUT" "$MODE"
rc=$?
set -e
log "prep_rna_refs.R exit $rc"
if [ -s "$OUT/summary.json" ]; then
  grep -m 1 '"status"' "$OUT/summary.json" || true
  awk '/"failed_checks"/, /\]/' "$OUT/summary.json" || true
fi
ls -la "$OUT"
exit "$rc"
