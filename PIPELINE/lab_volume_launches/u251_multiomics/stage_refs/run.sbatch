#!/bin/bash
# stage_refs: stage and verify every reference the released pipelines read, each against a source outside this build.
#   sbatch -D /rs/rs_grp_oschome/go2432/u251_multiomics/refs /rs/rs_grp_oschome/go2432/u251_multiomics/stage_refs/run.sbatch
# Steps (in order): containers gencode chain known_sites vep rnafusion rediportal(optional).
# A step writes $R/.done/<step> only when every one of its checks passed, so a requeue resumes at the first
# unfinished step. Check rows: $R/.checks/<step>.tsv, gathered into $R/REFS_CHECKS.tsv at the end.
# A failed step does not stop the others (only known_sites and chain need gencode). README.md says why each check.
# Exit: 0 every required step done; 2 at least one required step failed.
#SBATCH --job-name=u251_stage_refs
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/refs/logs/%x_%j.out
set -u
export LC_ALL=C
umask 0002

# STAGE_REFS_* overrides exist only for the local dry run (dryrun.sh); on the grid none is set
M=${STAGE_REFS_M:-/rs/rs_grp_oschome/go2432/u251_multiomics}
D=$M/stage_refs                     # this folder: scripts, expected S3 listing, bin/crc32c
R=$M/refs
DONE=$R/.done
CHK=$R/.checks
C=${STAGE_REFS_C:-/rs/rs_grp_oschome/go2432/singularity_cache}
U=${STAGE_REFS_U:-/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study}
H=${STAGE_REFS_H:-/rs/rs_grp_oschome/go2432/u251_host}
ENVB=$HOME/mambaforge/envs/nextflow/bin
PY=${STAGE_REFS_PY:-$ENVB/python3}
CURL=$ENVB/curl; [ -x "$CURL" ] || CURL=$(command -v curl)
SAM_IMG=$C/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img
BCF_TAG=1.21--h8b25389_0            # the tag nf-core/modules bcftools uses for 1.21; pinned, never :latest
BCF_URL=https://depot.galaxyproject.org/singularity/bcftools:$BCF_TAG
BCF_IMG=$C/depot.galaxyproject.org-singularity-bcftools-$BCF_TAG.img   # Nextflow's own name for that URL

export TMPDIR=$M/tmp/stage_refs     # node-local /tmp is swept on this grid
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache
export SINGULARITY_TMPDIR=$M/tmp/sing_tmp
export XDG_RUNTIME_DIR=$HOME/xdr
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export PATH=$ENVB:$PATH             # python3, curl, singularity 3.8.6
mkdir -p "$R/logs" "$DONE" "$CHK" "$TMPDIR" "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$XDG_RUNTIME_DIR" \
         "$M/stage_refs/results"

GFA=$R/gencode_v44/GRCh38.primary_assembly.genome.fa
GFAI=$GFA.fai
GGTF=$R/gencode_v44/GRCh38.primary_assembly.annotation.gtf
VEP_BASE=${STAGE_REFS_VEP_BASE:-https://ftp.ensembl.org/pub/release-115/variation/indexed_vep_cache}
RP_BASE=${STAGE_REFS_RP_BASE:-http://rediportal.cloud.ba.infn.it/download}
# rnafusion S3 subset (122.1 GB): NOT salmon, NOT gffread. Trailing slashes keep star/ from matching starfusion/
RF_PREFIXES=${STAGE_REFS_RF_PREFIXES:-"GRCh38/gencode_v46/starfusion/ctat_genome_lib_build_dir/ GRCh38/gencode_v46/star/ GRCh38/gencode_v46/fusioncatcher/ GRCh38/gencode_v46/gencode/ GRCh38/arriba/2.5.0/ GRCh38/hgnc/ GRCh38/fusion_report_db/"}
PRIMARY="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"

log() { echo "[$(date '+%F %T')] $*"; }
log "stage_refs: host $(hostname), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, cpus ${SLURM_CPUS_PER_TASK:-?}"
log "python $($PY --version 2>&1), $($CURL --version 2>/dev/null | head -1), singularity $(singularity --version 2>&1)"

# ------------------------------------------------------------------ helpers
STEP=none; CKF=/dev/null; STEP_FAIL=0
clean() { printf '%s' "$*" | tr '\t\n\r' '   '; }
ck() {   # name expected observed PASS|FAIL|INFO [outside_reference]
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$STEP" "$(clean "$1")" "$(clean "$2")" "$(clean "$3")" "$4" "$(clean "${5:-}")" >> "$CKF"
  echo "  [$4] $1: $(clean "$3")"
  [ "$4" = FAIL ] && STEP_FAIL=$((STEP_FAIL + 1))
  return 0
}
pf() { if [ "$1" = "$2" ] && [ -n "$1" ]; then echo PASS; else echo FAIL; fi; }
pyck() { # refs_util.py subcommand; its rows (5 columns) go to the step's check file with the step name in front
  local out rc
  out=$($PY "$D/refs_util.py" "$@"); rc=$?
  if [ -n "$out" ]; then
    printf '%s\n' "$out" | awk -v s="$STEP" 'BEGIN { FS = OFS = "\t" } { print s, $0 }' >> "$CKF"
    printf '%s\n' "$out" | awk -F'\t' '{ printf "  [%s] %s: %s\n", $4, $1, $3 }'
  fi
  if [ $rc -ne 0 ]; then
    if printf '%s\n' "$out" | awk -F'\t' '$4 == "FAIL" { f = 1 } END { exit !f }'; then
      STEP_FAIL=$((STEP_FAIL + 1))
    else
      ck "refs_util:$1" "runs" "exit $rc without a FAIL row (python error; see log)" FAIL
    fi
  fi
  return $rc
}
sx() { local img=$1; shift; singularity exec -B /rs/rs_grp_oschome/go2432 -B "$TMPDIR":/tmp "$img" "$@"; }
[ -n "${STAGE_REFS_LOCAL:-}" ] && sx() { shift; "$@"; }   # dry run: tools on PATH instead of the images
bx() { sx "$BCF_IMG" "$@"; }
clen() { # final Content-Length after redirects
  $CURL -sIL --fail --retry 3 --connect-timeout 60 "$1" | tr -d '\r' |
    awk 'tolower($1) == "content-length:" { v = $2 } END { print v }'
}
getfile() { # URL FINAL SIZE : resumable download to FINAL.part, moved into place only at the listed size
  local url=$1 out=$2 want=$3 i rc have now t0 dt
  if [ -f "$out" ] && [ "$(stat -c %s "$out")" = "$want" ]; then log "present at listed size: $out"; return 0; fi
  rm -f "$out"
  mkdir -p "$(dirname "$out")"
  for i in 1 2 3 4 5 6; do
    have=0; [ -f "$out.part" ] && have=$(stat -c %s "$out.part")
    if [ "$have" -gt "$want" ]; then rm -f "$out.part"; have=0; fi
    if [ "$have" -lt "$want" ] || [ "$want" -eq 0 ]; then
      t0=$(date +%s)
      if [ "$have" -gt 0 ]; then
        $CURL --fail --location --silent --show-error --retry 5 --retry-delay 20 --connect-timeout 60 \
              --speed-limit 10240 --speed-time 300 -C - -o "$out.part" "$url"; rc=$?
      else
        $CURL --fail --location --silent --show-error --retry 5 --retry-delay 20 --connect-timeout 60 \
              --speed-limit 10240 --speed-time 300 -o "$out.part" "$url"; rc=$?
      fi
      [ "$want" -eq 0 ] && [ ! -f "$out.part" ] && : > "$out.part"
      dt=$(( $(date +%s) - t0 )); now=$(stat -c %s "$out.part" 2>/dev/null || echo 0)
      log "fetch $(basename "$out") attempt $i: curl rc=$rc, $(( (now - have) / 1000000 )) MB in ${dt} s ($(( (now - have) / 1000000 / (dt > 0 ? dt : 1) )) MB/s)"
      [ "$rc" = 33 ] && rm -f "$out.part"   # server refused the range: start again
    fi
    if [ -f "$out.part" ] && [ "$(stat -c %s "$out.part")" = "$want" ]; then mv -f "$out.part" "$out"; return 0; fi
    sleep $((30 * i))
  done
  return 1
}
CRC=""
pick_crc32c() { # the static binary shipped in bin/, else compile; either must pass the CRC-32C check vector
  local c
  for c in "$D/bin/crc32c" "$R/tools/crc32c"; do
    if [ -x "$c" ] && [ "$(printf 123456789 | "$c" - 2>/dev/null)" = "e3069283" ]; then CRC=$c; return 0; fi
  done
  mkdir -p "$R/tools"
  for g in gcc cc /opt/ohpc/pub/compiler/gcc/7.3.0/bin/gcc; do
    if command -v "$g" > /dev/null 2>&1 && "$g" -O2 -o "$R/tools/crc32c" "$D/crc32c.c" 2>/dev/null &&
       [ "$(printf 123456789 | "$R/tools/crc32c" -)" = "e3069283" ]; then CRC=$R/tools/crc32c; return 0; fi
  done
  return 1
}
FAILED=(); DONE_NOW=()
run_step() { # name [optional]
  local s=$1 t0 rc
  if [ -f "$DONE/$s" ]; then log "== $s: already done ($(cat "$DONE/$s"))"; return 0; fi
  STEP=$s; CKF=$CHK/$s.tsv; : > "$CKF"; STEP_FAIL=0; t0=$(date +%s)
  log "== $s: start"
  "step_$s"; rc=$?
  if [ $rc -eq 0 ] && [ $STEP_FAIL -eq 0 ]; then
    echo "$(date -Is) job ${SLURM_JOB_ID:-none} host $(hostname) $(( $(date +%s) - t0 )) s" > "$DONE/$s"
    log "== $s: PASS in $(( $(date +%s) - t0 )) s"; DONE_NOW+=("$s")
  else
    log "== $s: FAILED (rc=$rc, failed checks=$STEP_FAIL) after $(( $(date +%s) - t0 )) s"
    [ "${2:-}" = optional ] || FAILED+=("$s")
  fi
  STEP=none; CKF=/dev/null
}

# ------------------------------------------------------------------ step: containers
step_containers() {
  local v v2 want
  [ -s "$SAM_IMG" ] || { ck samtools_image "exists" "missing $SAM_IMG" FAIL; return 1; }
  # stdout only, matched by line: singularity writes its own warnings to stderr
  v=$(sx "$SAM_IMG" samtools --version 2> /dev/null | grep -m1 '^samtools ')
  ck samtools_version "samtools 1.21" "$v" "$([[ "$v" == "samtools 1.21"* ]] && echo PASS || echo FAIL)" "image tag 1.21--h50ea8bc_0"
  if [ ! -s "$BCF_IMG" ]; then
    want=$(clen "$BCF_URL")
    [ -n "$want" ] || { ck bcftools_image "Content-Length from depot" "no answer from $BCF_URL" FAIL; return 1; }
    getfile "$BCF_URL" "$BCF_IMG" "$want" || { ck bcftools_image "downloaded" "failed" FAIL; return 1; }
    ck bcftools_image_size "$want (depot Content-Length)" "$(stat -c %s "$BCF_IMG")" "$(pf "$want" "$(stat -c %s "$BCF_IMG")")" "$BCF_URL"
  fi
  v=$(bx bcftools --version 2> /dev/null | grep -E -m2 '^(bcftools|Using htslib) ' | tr '\n' ' ')
  ck bcftools_version "bcftools 1.21 Using htslib 1.21" "$v" "$([[ "$v" == "bcftools 1.21 Using htslib 1.21"* ]] && echo PASS || echo FAIL)" "$BCF_URL"
  v=$(bx bgzip --version 2>&1 | grep -m1 '^bgzip '); v2=$(bx tabix --version 2>&1 | grep -m1 '^tabix ')
  ck htslib_tools "bgzip (htslib) 1.21 + tabix (htslib) 1.21 in the bcftools image" "$v; $v2" \
     "$([[ "$v" == "bgzip (htslib) 1.21"* && "$v2" == "tabix (htslib) 1.21"* ]] && echo PASS || echo FAIL)"
  ck bcftools_image_path "info" "$BCF_IMG" INFO "pinned tag $BCF_TAG"
  if pick_crc32c; then
    ck crc32c_tool "CRC-32C('123456789') = e3069283" "$CRC ok" PASS "RFC 3720 check value"
  else
    ck crc32c_tool "a working CRC-32C program" "none (bin/crc32c failed and no compiler worked)" FAIL
  fi
}

# ------------------------------------------------------------------ step: gencode
step_gencode() {
  local G=$R/gencode_v44 S=$U/ANALYSIS/refs/human name path got listed live hdr
  local fa_gz=$S/GRCh38.primary_assembly.genome.fa.gz gtf_gz=$S/GRCh38.primary_assembly.annotation.gtf.gz
  mkdir -p "$G"
  if ! $CURL -s --fail --retry 3 --connect-timeout 60 -o "$G/MD5SUMS.ebi_live" \
       https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/MD5SUMS; then rm -f "$G/MD5SUMS.ebi_live"; fi
  cp -f "$S/MD5SUMS" "$G/MD5SUMS.u_copy"
  for pair in "GRCh38.primary_assembly.genome.fa.gz|$fa_gz" \
              "$(basename "$(readlink -f "$gtf_gz")")|$gtf_gz"; do
    name=${pair%%|*}; path=${pair#*|}
    [ -r "$path" ] || { ck "md5:$name" "readable" "missing $path" FAIL; return 1; }
    got=$(md5sum < "$path" | awk '{ print $1 }')
    listed=$(awk -v n="$name" '$2 == n { print $1 }' "$S/MD5SUMS")
    ck "md5:$name:vs_U_MD5SUMS" "$listed" "$got" "$(pf "$listed" "$got")" "GENCODE release_44 MD5SUMS kept with the download"
    if [ -s "$G/MD5SUMS.ebi_live" ]; then
      live=$(awk -v n="$name" '$2 == n { print $1 }' "$G/MD5SUMS.ebi_live")
      ck "md5:$name:vs_EBI_live" "$live" "$got" "$(pf "$live" "$got")" "ftp.ebi.ac.uk GENCODE release_44 MD5SUMS, fetched now"
    else
      ck "md5:$name:vs_EBI_live" "EBI MD5SUMS" "not reachable from $(hostname)" INFO
    fi
  done
  gzip -dc "$fa_gz" > "$GFA.part" && mv -f "$GFA.part" "$GFA" || { ck gunzip_fasta "gzip CRC ok" "gzip failed" FAIL; return 1; }
  gzip -dc "$gtf_gz" > "$GGTF.part" && mv -f "$GGTF.part" "$GGTF" || { ck gunzip_gtf "gzip CRC ok" "gzip failed" FAIL; return 1; }
  ck gunzip "gzip exit 0 (CRC-32 and length verified)" "fasta $(stat -c %s "$GFA") bytes, gtf $(stat -c %s "$GGTF") bytes" PASS
  rm -f "$GFAI"
  sx "$SAM_IMG" samtools faidx "$GFA" && [ -s "$GFAI" ] || { ck faidx "samtools faidx exit 0" "failed" FAIL; return 1; }
  ck faidx "samtools faidx exit 0" "$(wc -l < "$GFAI") sequences" PASS
  pyck grch38 --fai "$GFAI"
  hdr=$(head -n 5 "$GGTF" | grep -m1 '^##description')
  ck gtf_header "contains 'version 44 (Ensembl 110)'" "$hdr" \
     "$([[ "$hdr" == *"version 44 (Ensembl 110)"* ]] && echo PASS || echo FAIL)" "GENCODE GTF header"
  return 0
}

# ------------------------------------------------------------------ step: chain
step_chain() {
  local CH=$R/chains f=hg38ToRn7.over.chain.gz url want live got
  url=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/$f
  mkdir -p "$CH"
  want=$(clen "$url"); [ -n "$want" ] || { ck chain_download "Content-Length" "no answer from UCSC" FAIL; return 1; }
  getfile "$url" "$CH/$f" "$want" || { ck chain_download "downloaded" "failed" FAIL; return 1; }
  if $CURL -s --fail --retry 3 -o "$CH/md5sum.txt.ucsc_live" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/md5sum.txt; then
    live=$(awk -v f="$f" '$2 == f { print $1 }' "$CH/md5sum.txt.ucsc_live")
  else
    live=""
  fi
  got=$(md5sum < "$CH/$f" | awk '{ print $1 }')
  if [ -n "$live" ]; then
    ck chain_md5 "$live (UCSC md5sum.txt, fetched now)" "$got" "$(pf "$live" "$got")" "hgdownload liftOver/md5sum.txt"
  else
    ck chain_md5 "4300ea13578990496b20fc94103e7562 (UCSC md5sum.txt on 2026-09-26)" "$got" \
       "$(pf 4300ea13578990496b20fc94103e7562 "$got")" "planning copy of UCSC md5sum.txt (live copy unreachable)"
  fi
  gzip -t "$CH/$f" && ck chain_gzip "gzip -t ok" ok PASS || ck chain_gzip "gzip -t ok" failed FAIL
  [ -f "$DONE/gencode" ] || { ck chain_sizes "GENCODE v44 fai (step gencode)" "gencode not done" FAIL; return 1; }
  pyck chain --chain "$CH/$f" --hfai "$GFAI" --rat "$H/results_rnaseq/genome/index/star/chrNameLength.txt"
  return 0
}

# ------------------------------------------------------------------ step: known_sites
gcs_get() { # URL OUT : download a Broad bundle object and check it against the hashes Google stores for it
  local url=$1 out=$2 hdr want crc md5 gen exp got try
  hdr=$($CURL -sI --fail --retry 3 --connect-timeout 60 "$url" | tr -d '\r')
  want=$(awk 'tolower($1) == "content-length:" { print $2 }' <<< "$hdr")
  crc=$(awk 'tolower($1) == "x-goog-hash:" && $2 ~ /^crc32c=/ { sub(/^crc32c=/, "", $2); print $2 }' <<< "$hdr")
  md5=$(awk 'tolower($1) == "x-goog-hash:" && $2 ~ /^md5=/ { sub(/^md5=/, "", $2); print $2 }' <<< "$hdr")
  gen=$(awk 'tolower($1) == "x-goog-generation:" { print $2 }' <<< "$hdr")
  [ -n "$want" ] && [ -n "$crc" ] || { ck "gcs_head:$(basename "$out")" "size and crc32c headers" "none" FAIL; return 1; }
  printf '%s\n' "$hdr" > "$out.headers"
  if [ -f "$out.gen" ] && [ "$(cat "$out.gen")" != "$gen" ]; then rm -f "$out" "$out.part"; fi   # object replaced
  echo "$gen" > "$out.gen"
  for try in 1 2; do
    getfile "$url" "$out" "$want" || { ck "download:$(basename "$out")" "$want bytes" "failed" FAIL; return 1; }
    exp=$($PY "$D/refs_util.py" b64hex "$crc"); got=$("$CRC" "$out")
    [ "$exp" = "$got" ] && break
    log "crc32c mismatch on $(basename "$out") ($got, expected $exp): fetching again"; rm -f "$out"
  done
  ck "crc32c:$(basename "$out")" "$exp (x-goog-hash, generation $gen)" "$got" "$(pf "$exp" "$got")" "Google Cloud Storage object hash"
  if [ -n "$md5" ]; then
    exp=$($PY "$D/refs_util.py" b64hex "$md5"); got=$(md5sum < "$out" | awk '{ print $1 }')
    ck "md5:$(basename "$out")" "$exp (x-goog-hash)" "$got" "$(pf "$exp" "$got")" "Google Cloud Storage object hash"
  fi
  [ "$exp" = "$got" ]
}
process_vcf() { # NAME SRC(bgzipped + tbi) OUT REF_CHECK_CONTIGS(comma list)
  local n=$1 src=$2 out=$3 w=$R/known_sites/work regs rr n_in
  bx bcftools index -s "$src" > "$w/$n.src.idxstats" || { ck "$n:index_stats" "bcftools index -s" "failed" FAIL; return 1; }
  bx bcftools view -h "$src" > "$w/$n.src.hdr" || { ck "$n:header" "bcftools view -h" "failed" FAIL; return 1; }
  pyck vcfhdr --header "$w/$n.src.hdr" --fai "$GFAI" --label "$n:source_header_lengths_vs_gencode" \
       --ref "Broad hg38 bundle header vs GENCODE v44 genome: two providers of GRCh38; a build mismatch differs here"
  regs=$(for c in $PRIMARY; do awk -v c="$c" '$1 == c { print c }' "$w/$n.src.idxstats"; done | paste -sd, -)
  ck "$n:regions_kept" "chr1-22,X,Y,M present in the source" "$regs" INFO
  bx bcftools view --threads 3 -r "$regs" -Oz -o "$w/$n.cut.vcf.gz.part" "$src" && mv -f "$w/$n.cut.vcf.gz.part" "$w/$n.cut.vcf.gz" ||
    { ck "$n:cut" "bcftools view -r" "failed" FAIL; return 1; }
  bx tabix -f -p vcf "$w/$n.cut.vcf.gz" || { ck "$n:cut_tabix" "sorted and indexable" "tabix failed" FAIL; return 1; }
  bx bcftools index -s "$w/$n.cut.vcf.gz" > "$w/$n.cut.idxstats"
  pyck idxsub --full "$w/$n.src.idxstats" --cut "$w/$n.cut.idxstats" --label "$n:cut_keeps_every_primary_record"
  bx bcftools reheader --fai "$GFAI" -o "$w/$n.rh.vcf.gz.part" "$w/$n.cut.vcf.gz" && mv -f "$w/$n.rh.vcf.gz.part" "$w/$n.rh.vcf.gz" ||
    { ck "$n:reheader" "bcftools reheader --fai" "failed" FAIL; return 1; }
  bx bcftools view -h "$w/$n.rh.vcf.gz" > "$w/$n.rh.hdr"
  if diff -q <(grep '^##contig=' "$w/$n.rh.hdr" | sed 's/^##contig=<ID=\([^,>]*\).*/\1/') <(cut -f1 "$GFAI") > /dev/null; then
    ck "$n:reheader_method" "info" "bcftools reheader --fai: contig lines are exactly the GENCODE fai, in its order" INFO
  else
    # --fai kept contig lines the GENCODE genome does not have: write the header explicitly from the fai instead
    awk -v FAI="$GFAI" 'BEGIN { while ((getline l < FAI) > 0) { split(l, a, "\t"); n++; c[n] = a[1]; L[n] = a[2] } }
         /^##contig=/ { next }
         /^#CHROM/ { for (i = 1; i <= n; i++) printf "##contig=<ID=%s,length=%s>\n", c[i], L[i] } { print }' \
         "$w/$n.rh.hdr" > "$w/$n.explicit.hdr"
    bx bcftools reheader -h "$w/$n.explicit.hdr" -o "$w/$n.rh.vcf.gz.part" "$w/$n.cut.vcf.gz" &&
      mv -f "$w/$n.rh.vcf.gz.part" "$w/$n.rh.vcf.gz" || { ck "$n:reheader_h" "bcftools reheader -h" "failed" FAIL; return 1; }
    ck "$n:reheader_method" "info" "--fai left contig lines outside GENCODE; header rebuilt from the fai with reheader -h" INFO
  fi
  mv -f "$w/$n.rh.vcf.gz" "$out"
  rm -f "$out.tbi"
  # tabix refuses anything that is not BGZF or not sorted, so success here is the format check
  bx tabix -f -p vcf "$out" || { ck "$n:final_tabix" "BGZF, sorted, indexable" "tabix failed" FAIL; return 1; }
  ck "$n:final_tabix" "BGZF, sorted, indexable" "tabix ok" PASS
  bx bcftools index -s "$out" > "$w/$n.final.idxstats"
  pyck idxcmp --a "$w/$n.cut.idxstats" --b "$w/$n.final.idxstats" --label "$n:counts_identical_across_reheader"
  bx bcftools view -h "$out" > "$w/$n.final.hdr"
  pyck vcfhdr --header "$w/$n.final.hdr" --fai "$GFAI" --label "$n:final_header_lengths_vs_gencode" --ref "GENCODE v44 fai"
  if diff -q <(grep '^##contig=' "$w/$n.final.hdr" | sed 's/^##contig=<ID=\([^,>]*\).*/\1/') <(cut -f1 "$GFAI") > /dev/null; then
    ck "$n:final_contig_dictionary" "GENCODE fai names in fai order (what GATK compares with the reference .dict)" "identical" PASS
  else
    ck "$n:final_contig_dictionary" "GENCODE fai names in fai order" "differs" FAIL
  fi
  # REF alleles against the GENCODE genome sequence itself: an evenly spaced sample of up to 20,000 records with plain
  # ACGT REF from the contigs in $4 that the file has, each REF compared with samtools faidx of the same span
  rr=$(for c in $(echo "$4" | tr ',' ' '); do awk -v c="$c" '$1 == c { print c }' "$w/$n.final.idxstats"; done | paste -sd, -)
  n_in=$(for c in $(echo "$rr" | tr ',' ' '); do awk -v c="$c" '$1 == c { print $3 }' "$w/$n.final.idxstats"; done |
         awk '{ s += $1 } END { print s + 0 }')
  [ -n "$rr" ] && [ "$n_in" -gt 0 ] || { ck "$n:REF_vs_GENCODE_genome" "records to test" "none in $4" FAIL; return 1; }
  bx bcftools query -r "$rr" -f '%CHROM\t%POS\t%REF\n' "$out" |
    $PY "$D/refs_util.py" refsample --out-regions "$w/$n.refsample.regions" --out-expect "$w/$n.refsample.expect" \
        --label "$n" --total "$n_in" --n 20000 > "$w/$n.refsample.rows"
  awk -v s="$STEP" 'BEGIN { FS = OFS = "\t" } { print s, $0 }' "$w/$n.refsample.rows" >> "$CKF"
  sx "$SAM_IMG" samtools faidx -r "$w/$n.refsample.regions" "$GFA" > "$w/$n.refsample.fa" 2> "$w/$n.refsample.err" ||
    { ck "$n:REF_vs_GENCODE_genome" "samtools faidx of the sample" "failed: $(tail -1 "$w/$n.refsample.err")" FAIL; return 1; }
  pyck redicmp --expect "$w/$n.refsample.expect" --faidx-out "$w/$n.refsample.fa" --min 0.99 \
       --label "$n:REF_vs_GENCODE_genome($rr)" \
       --ref "GENCODE v44 genome sequence; a build or coordinate mismatch makes most REF alleles disagree"
  return 0
}
step_known_sites() {
  local K=$R/known_sites B=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
  local mills=Mills_and_1000G_gold_standard.indels.hg38.vcf.gz dbs=Homo_sapiens_assembly38.dbsnp138.vcf f w
  [ -f "$DONE/gencode" ] && [ -f "$DONE/containers" ] || { ck gate "gencode and containers done" "not done" FAIL; return 1; }
  [ -n "$CRC" ] || pick_crc32c || { ck crc32c_tool "working" "missing" FAIL; return 1; }
  mkdir -p "$K/raw" "$K/work"; w=$K/work
  for f in "$mills" "$mills.tbi" "$dbs"; do gcs_get "$B/$f" "$K/raw/$f" || return 1; done
  # Mills: bgzipped by the Broad. Index a copy here (tabix refuses plain gzip and unsorted input); if tabix refuses
  # the copy, recompress it with bgzip and index again
  cp -f "$K/raw/$mills" "$w/mills.src.vcf.gz.part" && mv -f "$w/mills.src.vcf.gz.part" "$w/mills.src.vcf.gz"
  if ! bx tabix -f -p vcf "$w/mills.src.vcf.gz" 2> "$w/mills.tabix.err"; then
    ck mills:source_format "BGZF" "tabix refused the Broad file ($(tail -1 "$w/mills.tabix.err")); recompressed with bgzip" INFO
    gzip -dc "$K/raw/$mills" | bx bgzip -@ 4 -c > "$w/mills.src.vcf.gz.part" && mv -f "$w/mills.src.vcf.gz.part" "$w/mills.src.vcf.gz"
    bx tabix -f -p vcf "$w/mills.src.vcf.gz" || { ck mills:tabix "sorted and indexable" "tabix failed" FAIL; return 1; }
  fi
  process_vcf mills "$w/mills.src.vcf.gz" "$K/mills_1000g.primary.vcf.gz" "$(echo $PRIMARY | tr ' ' ',')" || return 1
  # dbSNP 138: uncompressed at the source
  if [ ! -f "$w/dbsnp138.src.ok" ]; then
    bx bgzip -@ 4 -c "$K/raw/$dbs" > "$w/dbsnp138.src.vcf.gz.part" && mv -f "$w/dbsnp138.src.vcf.gz.part" "$w/dbsnp138.src.vcf.gz" ||
      { ck dbsnp:bgzip "bgzip exit 0" "failed" FAIL; return 1; }
    bx tabix -f -p vcf "$w/dbsnp138.src.vcf.gz" || { ck dbsnp:tabix "sorted and indexable" "tabix failed" FAIL; return 1; }
    date -Is > "$w/dbsnp138.src.ok"
  fi
  process_vcf dbsnp "$w/dbsnp138.src.vcf.gz" "$K/dbsnp138.primary.vcf.gz" "chr21,chr22,chrY,chrM" || return 1
  if [ "$STEP_FAIL" -eq 0 ]; then   # keep the finals, hashes, headers and index stats; drop the 17 GB of intermediates
    rm -f "$K/raw/$dbs" "$w"/*.src.vcf.gz "$w"/*.src.vcf.gz.tbi "$w"/*.cut.vcf.gz "$w"/*.cut.vcf.gz.tbi
  fi
  ck final_files "info" "$(ls -l "$K"/*.primary.vcf.gz* | awk '{ printf "%s %s; ", $9, $5 }')" INFO
  return 0
}

# ------------------------------------------------------------------ step: vep
step_vep() {
  local V=$R/vep_cache base=$VEP_BASE
  local f=homo_sapiens_vep_115_GRCh38.tar.gz T want exp got try n_list n_disk miss c info
  T=$V/download/$f
  mkdir -p "$V/download"
  want=$(clen "$base/$f"); [ -n "$want" ] || { ck vep_download "Content-Length" "no answer from Ensembl" FAIL; return 1; }
  ck vep_listed_size "25272957721 bytes when the plan was made (2026-09-26)" "$want" \
     "$([ "$want" = 25272957721 ] && echo PASS || echo INFO)" "ftp.ensembl.org"
  $CURL -s --fail --retry 3 -o "$V/download/CHECKSUMS" "$base/CHECKSUMS" ||
    { ck vep_checksums "Ensembl CHECKSUMS" "not downloadable" FAIL; return 1; }
  exp=$(awk -v f="$f" '$3 == f { print $1 + 0, $2 + 0 }' "$V/download/CHECKSUMS")
  [ -n "$exp" ] || { ck vep_checksums "a line for $f" "absent" FAIL; return 1; }
  for try in 1 2; do
    getfile "$base/$f" "$T" "$want" || { ck vep_download "$want bytes" "failed" FAIL; return 1; }
    got=$(sum "$T" | awk '{ print $1 + 0, $2 + 0 }')
    [ "$got" = "$exp" ] && break
    log "VEP tarball sum $got, expected $exp: fetching again"; rm -f "$T"
  done
  ck vep_sum "$exp (Ensembl release-115 CHECKSUMS; BSD sum, 1 KiB blocks)" "$got" "$(pf "$exp" "$got")" \
     "$base/CHECKSUMS"
  [ "$got" = "$exp" ] || return 1
  rm -rf "$V/.extract"; mkdir -p "$V/.extract"
  tar -xzvf "$T" -C "$V/.extract" > "$V/download/tar_members.txt" 2> "$V/download/tar.err" ||
    { ck vep_extract "tar exit 0" "failed: $(tail -2 "$V/download/tar.err")" FAIL; return 1; }
  n_list=$(grep -vc '/$' "$V/download/tar_members.txt")
  n_disk=$(find "$V/.extract" \( -type f -o -type l \) | wc -l)
  ck vep_extract "every non-directory member on disk" "$n_list members, $n_disk on disk" "$(pf "$n_list" "$n_disk")" "the tarball's own listing"
  ck vep_top_level "homo_sapiens/115_GRCh38 only" "$(cut -d/ -f1-2 "$V/download/tar_members.txt" | sort -u | paste -sd' ' -)" INFO
  [ -f "$V/.extract/homo_sapiens/115_GRCh38/info.txt" ] || { ck vep_layout "homo_sapiens/115_GRCh38/info.txt" "absent" FAIL; return 1; }
  mkdir -p "$V/homo_sapiens"; rm -rf "$V/homo_sapiens/115_GRCh38"
  mv "$V/.extract/homo_sapiens/115_GRCh38" "$V/homo_sapiens/115_GRCh38" && rm -rf "$V/.extract"
  info=$V/homo_sapiens/115_GRCh38/info.txt
  ck vep_layout "\${vep_cache}/homo_sapiens/115_GRCh38/info.txt (rnavar usage.md l.236)" "$info" "$([ -s "$info" ] && echo PASS || echo FAIL)"
  ck vep_info_assembly "a line naming GRCh38" "$(grep -m1 -i 'assembly' "$info")" "$(grep -qi 'assembly.*GRCh38' "$info" && echo PASS || echo FAIL)" "info.txt written by Ensembl"
  ck vep_info_var_type "info" "$(grep -m1 -i '^var_type' "$info")" INFO "tabix = the indexed cache"
  miss=""; for c in $(seq 1 22) X Y MT; do [ -d "$V/homo_sapiens/115_GRCh38/$c" ] || miss="$miss $c"; done
  ck vep_chromosome_dirs "1-22, X, Y, MT (Ensembl names)" "${miss:-all present}" "$([ -z "$miss" ] && echo PASS || echo FAIL)"
  ck vep_size "info" "$(du -sh "$V/homo_sapiens/115_GRCh38" | cut -f1) extracted; tarball kept at $T" INFO
  return 0
}

# ------------------------------------------------------------------ step: rnafusion
step_rnafusion() {
  local RF=$R/rnafusion G46 try rc FR v
  G46=$RF/GRCh38/gencode_v46
  mkdir -p "$RF"
  for try in 1 2 3; do
    : > "$CHK/rnafusion.s3.tsv"
    $PY "$D/s3_mirror.py" --bucket nf-core-awsmegatests --root rnafusion/references/ \
      $(for p in $RF_PREFIXES; do printf -- '--prefix %s ' "$p"; done) \
      --dest "$RF" --expected "$D/expected_s3_listing.txt" --checks "$CHK/rnafusion.s3.tsv" --jobs 4 --curl "$CURL"
    rc=$?
    [ $rc -eq 0 ] && break
    log "s3_mirror attempt $try exit $rc"; sleep 60
  done
  awk -v s="$STEP" 'BEGIN { FS = OFS = "\t" } { print s, $0 }' "$CHK/rnafusion.s3.tsv" >> "$CKF"
  awk -F'\t' '{ printf "  [%s] %s: %s\n", $4, $1, $3 }' "$CHK/rnafusion.s3.tsv"
  rm -f "$CHK/rnafusion.s3.tsv"
  [ $rc -eq 0 ] || { STEP_FAIL=$((STEP_FAIL + 1)); return 1; }

  # 0-byte files in fusion_report_db. cosmic.db: nothing reads it (no_cosmic=true gives --no-cosmic, conf/modules.config
  # l.131). fusion_report.log: a log the download module moved there; nothing reads it. Either could be the file
  # exists_not_empty() lands on, which would send the pipeline to FUSIONREPORT_DOWNLOAD.
  FR=$RF/GRCh38/fusion_report_db
  mkdir -p "$RF/_set_aside"
  printf 'key\tsize\taction\treason\twhen\n' > "$RF/REMOVED.tsv.part"
  if [ -f "$FR/cosmic.db" ]; then
    if [ "$(stat -c %s "$FR/cosmic.db")" = 0 ]; then rm -f "$FR/cosmic.db"; else ck cosmic_db "0 bytes" "$(stat -c %s "$FR/cosmic.db") bytes" FAIL; fi
  fi
  printf 'GRCh38/fusion_report_db/cosmic.db\t0\tdeleted\tempty on S3 (MD5 d41d8cd98f00b204e9800998ecf8427e); no_cosmic=true runs FUSIONREPORT_DETECT with --no-cosmic\t%s\n' "$(date -Is)" >> "$RF/REMOVED.tsv.part"
  if [ -f "$FR/fusion_report.log" ] && [ "$(stat -c %s "$FR/fusion_report.log")" = 0 ]; then
    mv -f "$FR/fusion_report.log" "$RF/_set_aside/fusion_report.log"
  fi
  printf 'GRCh38/fusion_report_db/fusion_report.log\t0\tmoved to _set_aside/\tempty log file; nothing reads it; removes the second 0-byte file the directory check could land on\t%s\n' "$(date -Is)" >> "$RF/REMOVED.tsv.part"
  mv -f "$RF/REMOVED.tsv.part" "$RF/REMOVED.tsv"
  ck fusion_report_db_zero_byte_files "none left" "$(find "$FR" -type f -size 0 | wc -l) left; cosmic.db deleted, fusion_report.log set aside" \
     "$([ "$(find "$FR" -type f -size 0 | wc -l)" = 0 ] && echo PASS || echo FAIL)" "REMOVED.tsv"

  # what rnafusion 4.1.3 will decide when it looks at these paths
  pyck preflight --base "$RF" --ver 46

  # the S3 bundle against independent copies of the same assembly and against itself
  if [ -f "$DONE/gencode" ]; then
    pyck faicmp --a "$G46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" --b "$GFAI" --primary \
         --label v46_vs_v44_fai_primary --a-name "S3 GENCODE v46" --b-name "EBI GENCODE v44" \
         --ref "GRCh38 is one assembly: chr1-22,X,Y,M must agree across GENCODE releases"
    $PY "$D/refs_util.py" faicmp --a "$G46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" --b "$GFAI" \
         --label v46_vs_v44_fai_all_194 --a-name "S3 GENCODE v46" --b-name "EBI GENCODE v44" |
      awk -v s="$STEP" 'BEGIN { FS = OFS = "\t" } { $4 = ($4 == "FAIL") ? "INFO" : $4; print s, $0 }' | tee -a "$CKF"
  else
    ck v46_vs_v44_fai_primary "GENCODE v44 fai" "gencode step not done" INFO
  fi
  pyck faicmp --a "$G46/star/chrNameLength.txt" --b "$G46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" \
       --label star_index_built_on_the_v46_fasta --a-name "star/chrNameLength.txt" --b-name "v46 fai"
  pyck faicmp --a "$G46/starfusion/ctat_genome_lib_build_dir/ref_genome.fa.fai" \
       --b "$G46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" \
       --label ctat_lib_built_on_the_v46_fasta --a-name "ctat ref_genome.fa.fai" --b-name "v46 fai"
  pyck faicmp --a "$G46/starfusion/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt" \
       --b "$G46/gencode/Homo_sapiens_GRCh38_46_dna_primary_assembly.fa.fai" \
       --label ctat_star_index_on_the_v46_fasta --a-name "ctat star chrNameLength" --b-name "v46 fai"
  v=$(head -n 5 "$G46/gencode/Homo_sapiens_GRCh38_46.gtf" | grep -m1 '^##description')
  ck v46_gtf_header "contains 'version 46 (Ensembl 112)'" "$v" "$([[ "$v" == *"version 46 (Ensembl 112)"* ]] && echo PASS || echo FAIL)" "GENCODE GTF header"
  ck star_versionGenome "info (rnafusion 4.1.3 aligns with STAR 2.7.11b)" \
     "star: $(awk '$1 == "versionGenome" { print $2 }' "$G46/star/genomeParameters.txt"); ctat: $(awk '$1 == "versionGenome" { print $2 }' "$G46/starfusion/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt"); ctat sjdbOverhang $(awk '$1 == "sjdbOverhang" { print $2 }' "$G46/starfusion/ctat_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt")" INFO
  cp -f "$RF/MANIFEST.tsv" "$M/stage_refs/results/rnafusion_MANIFEST.tsv" 2> /dev/null
  return 0
}

# ------------------------------------------------------------------ step: rediportal (optional)
step_rediportal() {
  local RP=$R/rediportal f=TABLE1_hg38_v3.txt.gz url want why=""
  url=$RP_BASE/$f
  mkdir -p "$RP"
  want=$(clen "$url")
  if [ -z "$want" ]; then why="no Content-Length from $url"
  elif ! getfile "$url" "$RP/$f" "$want"; then why="download failed"
  elif ! gzip -t "$RP/$f" 2> /dev/null; then why="gzip -t failed"; rm -f "$RP/$f"
  fi
  if [ -z "$why" ]; then
    ck rediportal_size "$want (server Content-Length)" "$(stat -c %s "$RP/$f")" "$(pf "$want" "$(stat -c %s "$RP/$f")")" "$url"
    ck rediportal_gzip "gzip -t ok" ok PASS
    if [ -f "$DONE/gencode" ] && pyck redisample --table "$RP/$f" --out-regions "$RP/sample.regions" --out-expect "$RP/sample.expect"; then
      sx "$SAM_IMG" samtools faidx -r "$RP/sample.regions" "$GFA" > "$RP/sample.fa" 2> "$RP/sample.faidx.err"
      pyck redicmp --expect "$RP/sample.expect" --faidx-out "$RP/sample.fa" || why="REF sample check failed"
    else
      why="table checks failed or gencode not done"
    fi
  fi
  if [ -n "$why" ]; then
    printf 'REDIportal hg38 v3 not staged: %s (%s, job %s). S6 falls back to the strand rule.\n' "$why" "$(date -Is)" "${SLURM_JOB_ID:-none}" > "$RP/NOT_AVAILABLE"
    ck rediportal "staged" "NOT_AVAILABLE: $why" INFO
    return 1
  fi
  rm -f "$RP/NOT_AVAILABLE"
  return 0
}

# ------------------------------------------------------------------ run
[ -n "${STAGE_REFS_SOURCE_ONLY:-}" ] && return 0   # dry run sources the functions only
run_step containers
run_step gencode
run_step chain
run_step known_sites
run_step vep
run_step rnafusion
run_step rediportal optional

# ------------------------------------------------------------------ gather
{
  printf 'step\tcheck\texpected\tobserved\tresult\toutside_reference\n'
  for s in containers gencode chain known_sites vep rnafusion rediportal; do
    [ -f "$CHK/$s.tsv" ] && cat "$CHK/$s.tsv"
  done
  for f in "$CHK"/rnafusion_preflight_recheck_*.tsv; do [ -f "$f" ] && cat "$f"; done
} > "$R/REFS_CHECKS.tsv.part" && mv -f "$R/REFS_CHECKS.tsv.part" "$R/REFS_CHECKS.tsv"
{
  printf 'step\tstatus\tdetail\n'
  for s in containers gencode chain known_sites vep rnafusion rediportal; do
    if [ -f "$DONE/$s" ]; then printf '%s\tdone\t%s\n' "$s" "$(cat "$DONE/$s")"
    elif [ "$s" = rediportal ] && [ -f "$R/rediportal/NOT_AVAILABLE" ]; then printf '%s\tnot_available\t%s\n' "$s" "$(clean "$(cat "$R/rediportal/NOT_AVAILABLE")")"
    else printf '%s\tFAILED\tsee %s\n' "$s" "$CHK/$s.tsv"; fi
  done
} > "$R/STATUS.tsv"
cp -f "$R/REFS_CHECKS.tsv" "$R/STATUS.tsv" "$M/stage_refs/results/"
log "REFS_CHECKS: $(awk -F'\t' 'NR > 1 { n[$5]++ } END { for (k in n) printf "%s=%d ", k, n[k] }' "$R/REFS_CHECKS.tsv")"
cat "$R/STATUS.tsv"

# The directory check in rnafusion depends on readdir order, which the node that wrote the files may see differently
# from a fresh client. Repeat it once from another node; the recheck removes .done/rnafusion if it disagrees.
NODE=${SLURMD_NODENAME:-$(hostname -s)}
if [ -f "$DONE/rnafusion" ] && ! ls "$CHK"/rnafusion_preflight_recheck_*.tsv > /dev/null 2>&1 && [ ! -f "$R/.recheck_submitted" ]; then
  if command -v sbatch > /dev/null 2>&1 &&
     jid=$(sbatch --parsable --exclude="$NODE" --export=ALL,WRITER_NODE="$NODE" -D "$R" "$D/preflight_recheck.sbatch" 2>&1); then
    echo "$jid $(date -Is) from $NODE" > "$R/.recheck_submitted"; log "preflight recheck submitted: job $jid (excluding $NODE)"
  else
    log "preflight recheck NOT submitted from this node (${jid:-no sbatch}); the rnafusion gate must run it"
  fi
fi

if [ ${#FAILED[@]} -gt 0 ]; then log "required steps FAILED: ${FAILED[*]}"; exit 2; fi
log "every required step done"
exit 0
