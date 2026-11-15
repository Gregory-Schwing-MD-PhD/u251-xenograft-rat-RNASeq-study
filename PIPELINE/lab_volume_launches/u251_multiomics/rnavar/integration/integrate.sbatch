#!/bin/bash
# rnavar variant integration (afterany on the rnavar head job). Every step writes $W/.done/<step> when it finished, so a
# requeue or a resubmission resumes at the first unfinished step. Tools run in pinned containers:
#   pysam 0.24.1 (python), bcftools 1.21, samtools 1.21, UCSC liftOver 482, and STAR from the image that built the rat
#   index (htslib_samtools_star_gawk, STAR 2.7.11b, the rnaseq/rnavar STAR_GENOMEGENERATE image).
# If the head job deferred itself (inputs not ready), NEXTFLOW_EXIT does not exist yet and the current head id is in
# $D/HEAD_JOB: this job then resubmits itself afterany on that id and exits 0.
# Exit: 0 done; 2 hard pipeline checks failed (FORCE_INTEGRATION=1 runs the rest anyway); 3 no rnavar results;
#       4 a step failed (see the log; resubmit to resume).
#SBATCH --job-name=u251_rnavar_integ
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/rnavar/logs/%x_%j.out
set -u
umask 0002
export LC_ALL=C
# INTEG_* overrides exist only for the local synthetic test (tools on PATH instead of images); none is set on the grid
M=${INTEG_M:-/rs/rs_grp_oschome/go2432/u251_multiomics}
D=$M/rnavar
R=$D/results
I=$D/integration
W=$I/work
O=$I/results
C=${INTEG_C:-/rs/rs_grp_oschome/go2432/singularity_cache}
U=${INTEG_U:-/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study}
H=${INTEG_H:-/rs/rs_grp_oschome/go2432/u251_host}
FASTA=$M/refs/gencode_v44/GRCh38.primary_assembly.genome.fa
GTF=$M/refs/gencode_v44/GRCh38.primary_assembly.annotation.gtf
CHAIN=$M/refs/chains/hg38ToRn7.over.chain.gz
RATFA=$H/results_rnaseq/genome/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa
RATIDX=$H/results_rnaseq/genome/index/star
CYTO=$M/refs/ucsc/cytoBand.hg38.txt.gz
CYTO_URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
PYS=$C/depot.galaxyproject.org-singularity-pysam-0.24.1--py312hf5ad864_0.img
PYS_URL=https://depot.galaxyproject.org/singularity/pysam:0.24.1--py312hf5ad864_0
LIFT=$C/depot.galaxyproject.org-singularity-ucsc-liftover-482--h0b57e2e_0.img
LIFT_URL=https://depot.galaxyproject.org/singularity/ucsc-liftover:482--h0b57e2e_0
BCF=$C/depot.galaxyproject.org-singularity-bcftools-1.21--h8b25389_0.img
SAMI=$C/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img
STARI=$C/community-cr-prod.seqera.io-docker-registry-v2-blobs-sha256-26-268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d-data.img
LIBS="C2B IL64B IL66B IL67B IL68B IL69B N168B N269B NL70B NL71B"
LARGE="C2B IL66B IL67B IL68B IL69B N269B NL70B NL71B"
NP=${SLURM_CPUS_PER_TASK:-8}
FORCE_INTEGRATION=${FORCE_INTEGRATION:-0}
CLEAN_WORK=${CLEAN_WORK:-1}

export TMPDIR=$M/tmp/rnavar_integ
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH     # singularity 3.8.6, curl
mkdir -p "$W/.done" "$W/mp" "$W/mpc" "$W/rl" "$W/yx" "$W/rat" "$O" "$TMPDIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" \
         "$XDG_RUNTIME_DIR" "$(dirname "$CYTO")"
log() { echo "[$(date '+%F %T')] $*"; }
# sxrun IMAGE CMD...: every containerised command goes through this one wrapper (the parallel task lines too)
cat > "$W/sxrun" <<EOF
#!/bin/bash
img=\$1; shift
[ -n "${INTEG_LOCAL:-}" ] && exec "\$@"
exec singularity exec -B /rs/rs_grp_oschome/go2432 -B $U -B $TMPDIR:/tmp "\$img" "\$@"
EOF
chmod +x "$W/sxrun"
SX=$W/sxrun
sx() { "$SX" "$@"; }
PY() { sx "$PYS" python3 "$I/rnavar_integrate.py" "$@"; }
log "rnavar integration: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, cpus $NP"

# ---------------------------------------------------------------- gate: follow the head job until it has finished
if [ ! -f "$D/NEXTFLOW_EXIT" ]; then
  HJ=$(cat "$D/HEAD_JOB" 2>/dev/null || true)
  if [ -n "$HJ" ] && [ "$HJ" != "${SLURM_JOB_ID:-}" ] && squeue -h -j "$HJ" -o %i 2>/dev/null | grep -q .; then
    NEW=$(sbatch --parsable -D "$D/launch" --dependency=afterany:"$HJ" \
          --export=ALL,FORCE_INTEGRATION="$FORCE_INTEGRATION",CLEAN_WORK="$CLEAN_WORK" "$I/integrate.sbatch")
    NEW=${NEW%%;*}
    echo "$NEW" > "$D/INTEG_JOB"
    echo "$(date -Is) head $HJ still queued/running (no NEXTFLOW_EXIT); integration resubmitted as $NEW" | tee -a "$D/HEAD_JOBS.log"
    exit 0
  fi
  log "no NEXTFLOW_EXIT and no live head job (HEAD_JOB=${HJ:-none}): nothing to integrate"
  echo "no_nextflow_exit" > "$D/INTEGRATION_EXIT"; exit 3
fi
NFX=$(cat "$D/NEXTFLOW_EXIT")
case "$NFX" in
  gate_closed|input_missing|nextflow_version|gate_resubmit_failed)
    log "rnavar did not run (NEXTFLOW_EXIT=$NFX): nothing to integrate"; echo "no_results_$NFX" > "$D/INTEGRATION_EXIT"; exit 3 ;;
esac
log "rnavar NEXTFLOW_EXIT=$NFX"
rm -f "$D/INTEGRATION_EXIT"

FAILED_STEP=""
run() { # step function
  local s=$1 t0 rc
  if [ -f "$W/.done/$s" ]; then log "== $s: already done ($(cat "$W/.done/$s"))"; return 0; fi
  t0=$(date +%s); log "== $s: start"
  "step_$s"; rc=$?
  if [ $rc -eq 0 ]; then
    echo "$(date -Is) job ${SLURM_JOB_ID:-none} $(( $(date +%s) - t0 )) s" > "$W/.done/$s"; log "== $s: done in $(( $(date +%s) - t0 )) s"
  else
    log "== $s: FAILED (rc $rc) after $(( $(date +%s) - t0 )) s"; FAILED_STEP=$s
  fi
  return $rc
}
getimg() { # out url
  [ -s "$1" ] && return 0
  curl --fail --location --silent --show-error --retry 5 --retry-delay 20 --connect-timeout 60 -o "$1.part.$$" "$2" &&
    mv -f "$1.part.$$" "$1"
}
par() { # run the command lines in file $1, $NP at a time; fail if any failed
  local f=$1
  tr '\n' '\0' < "$f" | xargs -0 -P "$NP" -n 1 bash -c 'eval "$0" || { echo "TASK FAILED: $0" >&2; exit 1; }'
}

# ---------------------------------------------------------------- steps
step_containers() {
  [ -n "${INTEG_LOCAL:-}" ] && { log "local test: tools on PATH"; return 0; }
  getimg "$PYS" "$PYS_URL" || { log "pysam image not fetched"; return 1; }
  getimg "$LIFT" "$LIFT_URL" || { log "liftOver image not fetched"; return 1; }
  local v
  v=$(sx "$PYS" python3 -c 'import pysam, sys; print(pysam.__version__, sys.version.split()[0])' 2>/dev/null); log "pysam image: $v"
  [[ "$v" == 0.24.1* ]] || return 1
  v=$(sx "$BCF" bcftools --version 2>/dev/null | head -1); log "bcftools image: $v"; [[ "$v" == "bcftools 1.21"* ]] || return 1
  v=$(sx "$SAMI" samtools --version 2>/dev/null | head -1); log "samtools image: $v"; [[ "$v" == "samtools 1.21"* ]] || return 1
  v=$(sx "$STARI" STAR --version 2>/dev/null | head -1); log "STAR image: $v"; [ -n "$v" ] || return 1
  sx "$LIFT" liftOver 2>&1 | head -2 | sed 's/^/  liftOver: /'
  sx "$LIFT" liftOver 2>&1 | grep -q "liftOver - Move annotations" || return 1
  log "rat index: $(grep -E '^versionGenome' "$RATIDX/genomeParameters.txt")"
  for f in "$FASTA" "$FASTA.fai" "$GTF" "$CHAIN" "$RATFA" "$RATFA.fai" "$RATIDX/SA"; do
    [ -s "$f" ] || { log "missing $f"; return 1; }
  done
  for l in $LIBS; do [ -s "$U/ANALYSIS/results_human_final/star_salmon/$l.markdup.sorted.bam.bai" ] || log "old BAM index missing for $l"; done
}

step_mpc() { # samtools stats -r on 300 k R2 records of chr1, new and old BAM of each large library (spec: MPC check)
  local t=$W/mpc/tasks.txt l
  : > "$t"
  for l in $LARGE; do
    for side in new old; do
      if [ $side = new ]; then b=$R/preprocessing/${l}_graft/${l}_graft.md.bam; else b=$U/ANALYSIS/results_human_final/star_salmon/$l.markdup.sorted.bam; fi
      [ -s "$b" ] || continue
      echo "$SX $SAMI sh -c 'samtools view -h -f 0x80 -F 0xF04 -q 60 $b chr1 | head -n 300000 | samtools stats -r $FASTA - > $W/mpc/$l.$side.r2.stats.part' && mv -f $W/mpc/$l.$side.r2.stats.part $W/mpc/$l.$side.r2.stats" >> "$t"
    done
  done
  par "$t"
}

step_checks() {
  PY checks --results "$R" --work "$D/work" --prep-checks "$M/reads/checks" --nextflow-exit "$D/NEXTFLOW_EXIT" \
     --fasta "$FASTA" --u "$U" --metadata "$U/ANALYSIS/metadata_therapy.csv" --mpc-dir "$W/mpc" --out-dir "$W"
  local rc=$?
  cp -f "$W/CHECKS.tsv" "$O/CHECKS.tsv" 2>/dev/null
  awk -F'\t' 'NR > 1 && ($5 == "FAIL" || $5 == "WARN") { print "  " $5 ": " $1 " " $2 ": " $4 }' "$W/CHECKS.tsv" 2>/dev/null | head -60
  if [ $rc -ne 0 ]; then
    if [ "$FORCE_INTEGRATION" = 1 ]; then log "hard checks FAILED; FORCE_INTEGRATION=1, continuing (results are not trustworthy)"; return 0; fi
    return 2
  fi
}

step_sites() {
  PY sites --results "$R" --fasta "$FASTA" --gtf "$GTF" --work "$W"
  local rc=$?
  [ $rc -eq 3 ] && { log "authentication loci: problems recorded in $W/auth.json (reported as FAIL)"; return 0; }
  return $rc
}

step_mpileup() { # forced counts at every union site: new graft-only BAMs and old graft+both BAMs
  local t=$W/mp/tasks.txt l b
  : > "$t"
  for l in $LIBS; do
    for side in new old; do
      if [ $side = new ]; then b=$R/preprocessing/${l}_graft/${l}_graft.md.bam; else b=$U/ANALYSIS/results_human_final/star_salmon/$l.markdup.sorted.bam; fi
      [ -s "$b" ] || { log "no BAM for $l.$side"; continue; }
      [ -s "$W/mp/$l.$side.vcf.gz" ] && continue
      echo "$SX $BCF bcftools mpileup -f $FASTA -R $W/regions.tsv -a AD,ADF,ADR,DP -q 60 -Q 20 -d 1000000 -L 1000000 -B -O z -o $W/mp/$l.$side.vcf.gz.part $b 2> $W/mp/$l.$side.log && $SX $BCF bcftools index -t -f $W/mp/$l.$side.vcf.gz.part && mv -f $W/mp/$l.$side.vcf.gz.part.tbi $W/mp/$l.$side.vcf.gz.tbi && mv -f $W/mp/$l.$side.vcf.gz.part $W/mp/$l.$side.vcf.gz" >> "$t"
    done
  done
  par "$t" || return 1
  for l in $LARGE; do [ -s "$W/mp/$l.new.vcf.gz" ] || { log "mpileup output missing for $l.new"; return 1; }; done
}

step_counts() { PY counts --fasta "$FASTA" --work "$W" --mp-dir "$W/mp"; }
step_select() { PY select --work "$W"; }

step_liftover() {
  sx "$LIFT" liftOver "$W/liftover_in.bed" "$CHAIN" "$W/lifted.bed.part" "$W/unlifted.bed" && mv -f "$W/lifted.bed.part" "$W/lifted.bed"
  log "liftover: $(wc -l < "$W/liftover_in.bed") in, $(wc -l < "$W/lifted.bed") lifted to rn7"
}

step_readlevel() {
  local t=$W/rl/tasks.txt l b nm
  : > "$t"
  for l in $LIBS; do
    b=$R/preprocessing/${l}_graft/${l}_graft.md.bam
    [ -s "$b" ] && [ ! -s "$W/rl/$l.new.sites.tsv" ] && \
      echo "$SX $PYS python3 $I/rnavar_integrate.py readlevel --bam $b --tag $l.new --candidates $W/candidates.tsv --out-prefix $W/rl/$l.new 2> $W/rl/$l.new.log" >> "$t"
    b=$U/ANALYSIS/results_human_final/star_salmon/$l.markdup.sorted.bam; nm=$M/reads/both_names/$l.names.gz
    [ -s "$b" ] && [ ! -s "$W/rl/$l.old.sites.tsv" ] && \
      echo "$SX $PYS python3 $I/rnavar_integrate.py readlevel --bam $b --tag $l.old --candidates $W/candidates.tsv --both-names $nm --out-prefix $W/rl/$l.old 2> $W/rl/$l.old.log" >> "$t"
    b=$R/preprocessing/${l}_graft/${l}_graft.md.bam
    [ -s "$b" ] && [ ! -s "$W/yx/$l.new.yx.tsv" ] && \
      echo "$SX $PYS python3 $I/rnavar_integrate.py yx --bam $b --tag $l.new --work $W --out $W/yx/$l.new.yx.tsv 2> $W/yx/$l.new.log" >> "$t"
  done
  par "$t"
}

step_rat() { # realign the alt-supporting pairs to mRatBN7.2 with the host run's STAR index
  local r1=$W/rat/all_R1.fq r2=$W/rat/all_R2.fq f
  : > "$r1"; : > "$r2"
  for f in "$W"/rl/*.R1.fq; do cat "$f" >> "$r1"; cat "${f%.R1.fq}.R2.fq" >> "$r2"; done
  log "rat realignment: $(( $(wc -l < "$r1") / 4 )) pairs"
  if [ ! -s "$r1" ]; then : > "$W/rat/Aligned.out.sam"; return 0; fi
  sx "$STARI" STAR --runThreadN "$NP" --genomeDir "$RATIDX" --readFilesIn "$r1" "$r2" --outSAMtype SAM \
     --outSAMattributes NH HI AS nM --outSAMunmapped Within KeepPairs --outFilterMultimapNmax 20 \
     --outFileNamePrefix "$W/rat/" > "$W/rat/star.stdout" 2>&1 || { tail -20 "$W/rat/star.stdout"; return 1; }
  grep -E "Uniquely mapped reads %|% of reads mapped to multiple loci|Number of input reads" "$W/rat/Log.final.out" | sed 's/^/  rat STAR: /'
}

step_cytoband() { # UCSC hg38 cytoBand for chromosome arms (optional: without it the imbalance is per chromosome)
  if [ ! -s "$CYTO" ]; then
    curl --fail --location --silent --show-error --retry 5 --connect-timeout 60 -o "$CYTO.part" "$CYTO_URL" && mv -f "$CYTO.part" "$CYTO" ||
      { log "cytoBand not fetched: arms fall back to whole chromosomes"; rm -f "$CYTO.part"; return 0; }
  fi
  log "cytoBand: $(md5sum "$CYTO" | cut -d' ' -f1) $(stat -c %s "$CYTO") bytes ($CYTO_URL)"
}

step_final() {
  local redi=""
  [ -f "$M/refs/.done/rediportal" ] && redi=$M/refs/rediportal/TABLE1_hg38_v3.txt.gz
  PY final --results "$R" --work "$W" --fasta "$FASTA" --u "$U" --prep-checks "$M/reads/checks" --rl-dir "$W/rl" \
     --yx-dir "$W/yx" --index-dir "$M/reads/index_seq" --out-dir "$O" --rat-fasta "$RATFA" \
     --rat-sam "$W/rat/Aligned.out.sam" --lifted "$W/lifted.bed" --rediportal "$redi" --cytoband "$CYTO" \
     --invasion-dir "$M/invasion/results"
}

for s in containers mpc checks sites mpileup counts select liftover readlevel rat cytoband final; do
  run "$s"; rc=$?
  if [ $rc -ne 0 ]; then
    if [ "$s" = checks ]; then
      log "hard pipeline checks FAILED: $O/CHECKS.tsv. Nothing downstream ran. FORCE_INTEGRATION=1 overrides."
      echo "checks_failed" > "$D/INTEGRATION_EXIT"; exit 2
    fi
    echo "step_failed_$s" > "$D/INTEGRATION_EXIT"; exit 4
  fi
done
cp -f "$W/CHECKS.tsv" "$O/CHECKS.tsv" 2>/dev/null
cp -f "$W/r2_profiles.json" "$W/auth.json" "$O/" 2>/dev/null
echo "0" > "$D/INTEGRATION_EXIT"
log "integration done: $O/REPORT.md"
sed -n '1,40p' "$O/REPORT.md"

# the Nextflow work directory (300-400 GB) goes only when the hard checks passed without FORCE and everything ran
if [ "$CLEAN_WORK" = 1 ] && [ "$FORCE_INTEGRATION" != 1 ] && [ -f "$W/.done/checks" ] && [ -d "$D/work" ]; then
  sz=$(getfattr --only-values -n ceph.dir.rbytes "$D/work" 2>/dev/null || echo unknown)
  log "deleting the Nextflow work directory $D/work ($sz bytes)"
  rm -rf "$D/work" && log "work deleted"
fi
exit 0
