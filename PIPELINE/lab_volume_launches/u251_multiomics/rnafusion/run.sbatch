#!/bin/bash
# rnafusion head job: nf-core/rnafusion 4.1.3 under Nextflow 25.10.7. RF_PASS=1 (callers + FusionInspector on the
# CCLE whitelist) or RF_PASS=2 (same launch dir, -resume, --whitelist whitelist_pass2.tsv: only CAT_CAT,
# FusionInspector, the VCF step, visualisation and MultiQC should rerun). Tasks go to slurm through grid.config.
# Submitted by submit.sh; pass 2 is submitted with -J u251_rf_p2 --export=ALL,RF_PASS=2.
#
# GATE (pass 1): M/reads/READY and M/refs/.done/rnafusion exist, and no prep/stage_refs/recheck job is still queued
#   (the second-node recheck can still delete .done/rnafusion). While anything upstream is queued this job resubmits
#   itself with afterany on it and exits 0; when nothing is queued and a marker is missing the gate closes (exit 3).
# GATE (pass 2): pass 1 exit 0, REFCHECK_PASS1 = OK, audit/AUDIT_OK, whitelist_pass2.tsv not empty.
# Before launch: the exists_not_empty() emulation is rerun on THIS node (readdir order decides which file rnafusion
#   tests), every whitelist gene must exist in the v46 GTF and in the CTAT library's gene spans, and the samplesheet
#   must be exactly the ten <lib>_humclip rows.
# WATCHDOG: if Nextflow submits any BUILD_REFERENCES task (a silently rebuilt or downloaded reference; STAR-Fusion's
#   build is 2 days at 150 GB) the head kills Nextflow. In pass 2 it also kills Nextflow if a caller or STAR reruns
#   (resume failed), and in both passes if the launched revision is not the 4.1.3 commit.
#SBATCH --job-name=u251_rf_p1
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=4-00:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/logs/rnafusion/%x_%j.out
set -uo pipefail
source /rs/rs_grp_oschome/go2432/u251_multiomics/rnafusion/common.sh
P=${RF_PASS:-1}
case "$P" in 1|2) ;; *) echo "RF_PASS must be 1 or 2 (got '$P')"; exit 9 ;; esac
EXITF=$D/NEXTFLOW_EXIT_PASS$P
WL=$D/whitelist_pass$P.tsv
mkdir -p "$D/launch" "$D/checks" "$LOGS"
echo "${SLURM_JOB_ID:-none}" > "$D/HEAD${P}_JOB"
log "rnafusion head pass $P: host $(hostname -s), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, follow ${RF_FOLLOW:-0}"

close_gate() { log "GATE CLOSED (pass $P): $*"; echo "gate_closed: $*" > "$EXITF"; exit 3; }
stop()       { log "STOP (pass $P): $2"; echo "$1: $2" > "$EXITF"; exit "${3:-4}"; }

# ------------------------------------------------------------------ gate
why=""
if [ "$P" = 1 ]; then
  [ -f "$M/reads/READY" ]         || why="$why M/reads/READY missing;"
  [ -f "$M/refs/.done/rnafusion" ] || why="$why M/refs/.done/rnafusion missing;"
  UP="u251_prep u251_prep_collect u251_stage_refs u251_refs_recheck"
else
  [ "$(head -n 1 "$D/NEXTFLOW_EXIT_PASS1" 2>/dev/null)" = 0 ] || why="$why pass 1 exit is '$(head -n 1 "$D/NEXTFLOW_EXIT_PASS1" 2>/dev/null)';"
  [ "$(head -n 1 "$D/REFCHECK_PASS1" 2>/dev/null)" = OK ]    || why="$why REFCHECK_PASS1 is not OK;"
  [ -f "$D/audit/AUDIT_OK" ] || why="$why audit/AUDIT_OK missing;"
  [ -s "$WL" ]               || why="$why $WL missing or empty;"
  UP="u251_rf_p1 u251_rf_audit"
fi
# follow anything upstream that is still queued, even when the markers exist (the refs recheck may still remove one)
follow "$D/run.sbatch" RF_FOLLOW 40 $UP; fr=$?
[ "$fr" = 0 ] && exit 0
[ "$fr" = 2 ] && close_gate "follow limit reached;$why"
if [ "$P" = 1 ] && [ -n "$why" ]; then
  # nothing upstream is queued and a marker is missing: a recorded failure closes the gate; an unfinished upstream
  # (stage_refs killed at its 12 h limit writes no STATUS.tsv) is re-checked every 2 h for 24 h, so a resubmitted
  # stage_refs or prep_reads is picked up without resubmitting this chain
  if [ ! -f "$M/refs/.done/rnafusion" ] && grep -qE '^rnafusion[[:space:]]+FAILED' "$M/refs/STATUS.tsv" 2>/dev/null; then
    close_gate "stage_refs recorded the rnafusion step as FAILED ($M/refs/STATUS.tsv);$why"
  fi
  if ls "$M/refs/rnafusion"/PREFLIGHT_FAILED_* > /dev/null 2>&1; then
    close_gate "the second-node preflight failed: $(ls "$M/refs/rnafusion"/PREFLIGHT_FAILED_* | tr '\n' ' ');$why"
  fi
  if [ ! -f "$M/reads/READY" ] && [ -f "$M/reads/CHECKS.tsv" ]; then
    close_gate "prep_reads finished without READY ($M/reads/CHECKS.tsv);$why"
  fi
  w=${RF_WAIT:-0}
  if [ "$w" -lt 12 ]; then
    new=$(sbatch --parsable -J "${SLURM_JOB_NAME:-u251_rf_p1}" -D "$D/launch" --begin=now+2hours \
            --export=ALL,RF_WAIT=$((w + 1)),RF_FOLLOW=0 "$D/run.sbatch"); new=${new%%;*}
    log "deferred 2 h (RF_WAIT=$((w + 1)) of 12): $why -> $new" | tee -a "$D/JOBS.log"
    echo "$new" > "$D/HEAD1_JOB"
    exit 0
  fi
  close_gate "markers still missing after 12 x 2 h;$why"
fi
[ -z "$why" ] || close_gate "$why"
rm -f "$EXITF"
log "gate open"
[ "$P" = 1 ] && { cat "$M/reads/READY"; cat "$M/refs/.done/rnafusion"; }

# ------------------------------------------------------------------ 1. samplesheet and inputs
SS=$D/samplesheet.csv
[ "$(head -n 1 "$SS" | tr -d '\r')" = "sample,fastq_1,fastq_2,strandedness" ] || stop input_check "samplesheet header differs" 4
bad=""; n=0
while IFS=, read -r s f1 f2 st; do
  st=${st%$'\r'}
  [ "$s" = sample ] && continue
  n=$((n + 1)); lib=${s%_humclip}
  case " $LIBS " in *" $lib "*) ;; *) bad="$bad $s:unknown_library" ;; esac
  [ "$s" = "${lib}_humclip" ] || bad="$bad $s:name"
  [ "$st" = reverse ] || bad="$bad $s:strandedness=$st"
  [ "$f1" = "$U/ANALYSIS/sorted_fastqs/${lib}_human_R1.fq.gz" ] || bad="$bad $s:fastq_1_path"
  [ "$f2" = "$M/reads/human_clip/${lib}_human_R2.fq.gz" ]       || bad="$bad $s:fastq_2_path"
  for f in "$f1" "$f2"; do
    if [ -s "$f" ]; then echo "input ok  $s  $(stat -c %s "$f") bytes  $f"; else bad="$bad $s:missing:$f"; fi
  done
  [ -s "$M/reads/checks/${lib}.tsv" ] || bad="$bad $s:no_prep_check_row"
done < "$SS"
[ "$n" = 10 ] || bad="$bad rows=$n"
[ -z "$bad" ] || stop input_check "samplesheet/input problems:$bad" 4

# ------------------------------------------------------------------ 2. exists_not_empty() emulation on this node
PF=$D/checks/preflight_pass${P}_$(hostname -s)_${SLURM_JOB_ID:-none}.tsv
"$PY" "$M/stage_refs/refs_util.py" preflight --base "$REFS" --ver 46 > "$PF"; rc=$?
cat "$PF"
[ "$rc" = 0 ] || stop preflight_failed "rnafusion would rebuild or download a reference from this node: $PF" 5

# ------------------------------------------------------------------ 3. whitelist: format, and every gene annotated
WLG=$D/checks/whitelist_genes_pass$P.tsv
badl=$(tr -d '\r' < "$WL" | awk -F'--' 'NF && (NF != 2 || $1 == "" || $2 == "" || $0 ~ /[[:space:]]/)' | head -n 5)
[ -z "$badl" ] || stop whitelist_format "lines that are not GENE--GENE in $WL: $badl" 6
genes=$(tr -d '\r' < "$WL" | awk -F'--' 'NF == 2 { print $1; print $2 }' | sort -u | tr '\n' ' ')
log "whitelist $WL: $(tr -d '\r' < "$WL" | awk 'NF' | wc -l) fusions, $(echo $genes | wc -w) genes"
{
  # GENCODE v46 GTF (Arriba's and the STAR index's annotation): gene_name attribute
  awk -v want="$genes" 'BEGIN { n = split(want, a, " "); for (i = 1; i <= n; i++) w[a[i]] = 1 }
       $3 == "gene" && match($0, /gene_name "[^"]*"/) { g = substr($0, RSTART + 11, RLENGTH - 12); if (g in w) f[g] = 1 }
       END { for (g in w) print "gencode_v46_gtf\t" g "\t" ((g in f) ? "found" : "MISSING") }' \
       "$G/gencode/Homo_sapiens_GRCh38_46.gtf"
  # CTAT genome library (FusionInspector's and STAR-Fusion's annotation): any field of ref_annot.gtf.gene_spans
  awk -F'\t' -v want="$genes" 'BEGIN { n = split(want, a, " "); for (i = 1; i <= n; i++) w[a[i]] = 1 }
       { for (i = 1; i <= NF; i++) if ($i in w) f[$i] = 1 }
       END { for (g in w) print "ctat_gene_spans\t" g "\t" ((g in f) ? "found" : "MISSING") }' \
       "$G/starfusion/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans"
} | sort > "$WLG"
nmiss=$(grep -c MISSING "$WLG"); nfound=$(grep -c found "$WLG")
log "whitelist genes: $nfound found, $nmiss missing ($WLG)"
[ "$nfound" -gt 0 ] || stop whitelist_genes "no whitelist gene found at all: the annotation files are not what this check expects" 6
[ "$nmiss" = 0 ]    || stop whitelist_genes "genes missing from the annotation: $(grep MISSING "$WLG" | cut -f1,2 | tr '\t\n' ': ')" 6

# ------------------------------------------------------------------ 4. pass 2: set pass 1's FusionInspector-dependent outputs aside (once)
if [ "$P" = 2 ] && [ ! -f "$D/.pass1_outputs_moved" ]; then
  mkdir -p "$RESULTS/_pass1"
  for d in fusioninspector cat vcf arriba_visualisation agat multiqc; do
    if [ -e "$RESULTS/$d" ]; then mv "$RESULTS/$d" "$RESULTS/_pass1/$d" && log "moved results/$d -> results/_pass1/$d"; fi
  done
  date -Is > "$D/.pass2_start"
  touch "$D/.pass1_outputs_moved"
fi

# ------------------------------------------------------------------ 5. environment (home is at 97 %: everything on CephFS)
export NXF_VER=$RF_NXF_VER
export NXF_HOME=$M/nxf_home
export NXF_SINGULARITY_CACHEDIR=$M/singularity
export NXF_SINGULARITY_HOME_MOUNT=true
export NXF_TEMP=$M/tmp
export TMPDIR=$M/tmp
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export NXF_OPTS='-Xms1g -Xmx6g'
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH
mkdir -p "$NXF_HOME" "$NXF_TEMP" "$NXF_SINGULARITY_CACHEDIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR"
log "nextflow launcher: $(command -v nextflow); singularity: $(command -v singularity) ($(singularity --version 2>&1))"
NV=$(nextflow -version 2>&1)
echo "$NV" | grep -iE "version" | head -n 2
echo "$NV" | grep -q "version $RF_NXF_VER" || stop nextflow_version "NXF_VER=$RF_NXF_VER did not give Nextflow $RF_NXF_VER" 8

# ------------------------------------------------------------------ 6. rnafusion tasks left by an earlier attempt of this head
# Only one rnafusion head runs at a time, so any rnafusion task still in the queue now was orphaned by a requeue of
# this head. A resumed run would give the same task the same work directory while the orphan still writes into it.
orph=$(squeue -h -u "$USER" -o "%i %j" | awk '$2 ~ /^nf-NFCORE_RNAFUSION/ { print $1 }' | tr '\n' ' ')
if [ -n "${orph// /}" ]; then
  log "cancelling rnafusion tasks orphaned by an earlier attempt: $orph"
  scancel $orph
  sleep 60
fi

# ------------------------------------------------------------------ 7. run, with the watchdog
cd "$D/launch"
NFOUT=$LOGS/nextflow_pass${P}_${SLURM_JOB_ID:-none}_r${SLURM_RESTART_COUNT:-0}.out
log "nextflow output: $NFOUT"
nextflow run nf-core/rnafusion -r "$RF_TAG" \
  -profile singularity \
  -c ../grid.config \
  -params-file ../params.yaml \
  --whitelist "$WL" \
  -w ../work \
  -resume \
  -ansi-log false > "$NFOUT" 2>&1 &
NFPID=$!

REBUILD_RE='(Submitted|Cached) process > [^ ]*(BUILD_REFERENCES|STARFUSION_BUILD|STAR_GENOMEGENERATE|GENCODE_DOWNLOAD|HGNC_DOWNLOAD|FUSIONREPORT_DOWNLOAD|SALMON_INDEX|SAMTOOLS_FAIDX|GATK4_CREATESEQUENCEDICTIONARY|BIOAWK|CTATSPLICING_PREPGENOMELIB|GFFREAD|UCSC_GTFTOGENEPRED)'
HEAVY_RE='Submitted process > [^ ]*(STAR_ALIGN|FUSIONCATCHER_FUSIONCATCHER|ARRIBA_ARRIBA|STARFUSION_DETECT|GATK4_MARKDUPLICATES|CTATSPLICING_STARTOCANCERINTRONS|FASTQC)'
killed=""
while kill -0 "$NFPID" 2>/dev/null; do
  if grep -qE "$REBUILD_RE" "$NFOUT"; then
    killed=reference_rebuild
    log "WATCHDOG: a reference-building task was submitted: $(grep -m 3 -E "$REBUILD_RE" "$NFOUT" | tr '\n' ' ')"
  elif [ "$P" = 2 ] && grep -qE "$HEAVY_RE" "$NFOUT"; then
    killed=resume_failed
    log "WATCHDOG: pass 2 reran a pass-1 step (resume did not reuse the cache): $(grep -m 3 -E "$HEAVY_RE" "$NFOUT" | tr '\n' ' ')"
  elif grep -q "revision:" "$NFOUT" && ! grep -q "revision: ${RF_TAG_SHA:0:10} \[$RF_TAG\]" "$NFOUT"; then
    killed=wrong_revision
    log "WATCHDOG: launched revision is not ${RF_TAG_SHA:0:10} [$RF_TAG]: $(grep -m 1 'revision:' "$NFOUT")"
  fi
  if [ -n "$killed" ]; then
    kill -TERM "$NFPID" 2>/dev/null
    for _ in $(seq 1 30); do kill -0 "$NFPID" 2>/dev/null || break; sleep 10; done
    kill -KILL "$NFPID" 2>/dev/null
    left=$(squeue -h -u "$USER" -o "%i %j" | awk '$2 ~ /^nf-NFCORE_RNAFUSION/ { print $1 }' | tr '\n' ' ')
    [ -n "${left// /}" ] && { log "cancelling remaining rnafusion tasks: $left"; scancel $left; }
    break
  fi
  sleep 60
done
wait "$NFPID"; rc=$?
tail -n 40 "$NFOUT"
log "nextflow exit $rc${killed:+ (killed by the watchdog: $killed)}"
grep -m 1 "revision:" "$NFOUT" | sed 's/^/launched: /'
if [ -n "$killed" ]; then echo "$killed" > "$EXITF"; else echo "$rc" > "$EXITF"; fi

# ------------------------------------------------------------------ 8. reference check: every execution trace of this run directory
REFC=$D/REFCHECK_PASS$P
traces=$(ls "$RESULTS"/pipeline_info/execution_trace_*.txt 2>/dev/null)
if [ -z "$traces" ]; then
  echo "NO_TRACE" > "$REFC"
else
  hits=$(grep -hE 'BUILD|GENOMEGENERATE|DOWNLOAD|SALMON_INDEX|GENCODE|HGNC|FAIDX|CREATESEQUENCEDICTIONARY|BIOAWK' $traces)
  if [ -n "$hits" ]; then { echo FAIL; echo "$hits"; } > "$REFC"
  else { echo OK; echo "$(echo "$traces" | wc -l) execution trace(s) searched; no reference-building or download task"; } > "$REFC"; fi
fi
log "reference check: $(head -n 1 "$REFC")"
if [ "$(head -n 1 "$REFC")" = FAIL ] && [ "$(head -n 1 "$EXITF")" = 0 ]; then
  echo "reference_rebuild" > "$EXITF"
  log "the run is NOT accepted: a reference was rebuilt or downloaded (see $REFC)"
fi

# ------------------------------------------------------------------ 9. pass 2: which steps actually reran
if [ "$P" = 2 ] && [ -f "$D/.pass2_start" ]; then
  T2=$(ls -t "$RESULTS"/pipeline_info/execution_trace_*.txt 2>/dev/null | head -n 1)
  RR=$D/PASS2_RERUN_CHECK.tsv
  if [ -n "$T2" ]; then
    awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
         { nm = $(c["name"]); sub(/ \(.*$/, "", nm); n = split(nm, a, ":"); p = a[n]; s = $(c["status"])
           k = p "\t" s; cnt[k]++ }
         END { for (k in cnt) print k "\t" cnt[k] }' "$T2" | sort > "$RR.tmp"
    allowed='^(CAT_CAT|FUSIONINSPECTOR|AGAT_CONVERTSPGFF2TSV|VCF_COLLECT|ARRIBA_VISUALISATION|MULTIQC)$'
    unexpected=$(awk -F'\t' -v re="$allowed" '$2 != "CACHED" && $1 !~ re' "$RR.tmp")
    { if [ -z "$unexpected" ]; then echo -e "verdict\tOK\tonly FusionInspector-dependent steps reran ($T2)"
      else echo -e "verdict\tWARN\tsteps other than the FusionInspector-dependent ones reran ($T2)"; fi
      echo -e "process\tstatus\tn"; cat "$RR.tmp"; } > "$RR"
    rm -f "$RR.tmp"
    cat "$RR"
  fi
fi

# ------------------------------------------------------------------ 10. per-sample outputs the audit reads (sizes only; the audit checks content)
for lib in $LIBS; do
  s=${lib}_humclip
  printf '%s' "$s"
  for f in "arriba/$s.arriba.fusions.tsv" "starfusion/$s.starfusion.fusion_predictions.tsv" \
           "fusioncatcher/$s.fusion-genes.txt" "fusionreport/$s/$s.fusionreport_filtered.tsv" \
           "fusioninspector/$s/$s.FusionInspector.fusions.abridged.tsv" "star/$s.Log.final.out"; do
    printf '  %s=%s' "${f%%/*}" "$(stat -c %s "$RESULTS/$f" 2>/dev/null || echo MISSING)"
  done
  echo
done
log "head pass $P done: exit file says '$(head -n 1 "$EXITF")'"
[ "$(head -n 1 "$EXITF")" = 0 ] && exit 0 || exit 1
