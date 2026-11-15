#!/bin/bash
# rnafusion pass-1 audit (afterany on the pass-1 head). Follows the head if it deferred itself; runs only when pass 1
# exited 0 and the reference check says OK. Writes M/rnafusion/audit/{AUDIT.md,IDENTITY.md,recall_ccle.tsv,
# scored_calls.tsv,reportable_new.tsv,checks.tsv,whitelist_pass2*.tsv} and AUDIT_OK only when every hard check passes;
# copies whitelist_pass2.tsv to M/rnafusion/ for the pass-2 head.
#SBATCH --job-name=u251_rf_audit
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/rs/rs_grp_oschome/go2432/u251_multiomics/logs/rnafusion/%x_%j.out
set -uo pipefail
source /rs/rs_grp_oschome/go2432/u251_multiomics/rnafusion/common.sh
A=$D/audit
mkdir -p "$A"
echo "${SLURM_JOB_ID:-none}" > "$D/AUDIT_JOB"
log "rnafusion audit: host $(hostname -s), job ${SLURM_JOB_ID:-none}, restart ${SLURM_RESTART_COUNT:-0}, follow ${RF_FOLLOW:-0}"

follow "$D/audit.sbatch" RF_FOLLOW 40 u251_rf_p1; fr=$?
[ "$fr" = 0 ] && exit 0
skip() { log "AUDIT NOT RUN: $*"; echo "$*" > "$A/AUDIT_SKIPPED"; exit 3; }
[ "$fr" = 2 ] && skip "follow limit reached"
[ "$(head -n 1 "$D/NEXTFLOW_EXIT_PASS1" 2>/dev/null)" = 0 ] || skip "pass 1 exit is '$(head -n 1 "$D/NEXTFLOW_EXIT_PASS1" 2>/dev/null)'"
[ "$(head -n 1 "$D/REFCHECK_PASS1" 2>/dev/null)" = OK ]    || skip "REFCHECK_PASS1 is '$(head -n 1 "$D/REFCHECK_PASS1" 2>/dev/null)'"
rm -f "$A/AUDIT_OK" "$A/AUDIT_FAILED" "$A/AUDIT_SKIPPED"

export PATH=$HOME/mambaforge/envs/nextflow/bin:$PATH      # singularity 3.8.6 with a working starter on compute nodes
export SINGULARITY_CACHEDIR=$M/tmp/sing_cache SINGULARITY_TMPDIR=$M/tmp/sing_tmp TMPDIR=$M/tmp
mkdir -p "$SINGULARITY_CACHEDIR" "$SINGULARITY_TMPDIR" "$A/containers"

# STAR for the rat realignment: the image rnafusion's STAR_ALIGN used (seqera htslib_samtools_star_gawk), from the cache
STAR_SHA=268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d
IMG=$(ls "$M"/singularity/*"$STAR_SHA"* /rs/rs_grp_oschome/go2432/singularity_cache/*"$STAR_SHA"* "$A"/containers/*"$STAR_SHA"* 2>/dev/null | head -n 1)
if [ -z "$IMG" ]; then
  log "STAR image not in the caches; pulling it"
  singularity pull "$A/containers/star_$STAR_SHA.sif" \
    "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/$STAR_SHA/data" && IMG=$A/containers/star_$STAR_SHA.sif
fi
log "STAR image: ${IMG:-none}"
log "rat index: $(grep -E '^versionGenome|^genomeFastaFiles' "$H/results_rnaseq/genome/index/star/genomeParameters.txt" | tr '\n' ' ')"

"$PY" "$D/audit_rnafusion.py" all \
  --results "$RESULTS" --audit "$A" --reads "$M/reads" --u "$U" \
  --hgnc "$REFS/GRCh38/hgnc/hgnc_complete_set.txt" \
  --gtf "$G/gencode/Homo_sapiens_GRCh38_46.gtf" \
  --ctat-spans "$G/starfusion/ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans" \
  --rat-gtf "$H/results_rnaseq/genome/Rattus_norvegicus.mRatBN7.2.110.filtered.gtf" \
  --rat-index "$H/results_rnaseq/genome/index/star" \
  --orthologs "$U/ANALYSIS/host/refs_online/rat_human_1to1.tsv" \
  --star-img "${IMG:-}" --tmp "$M/tmp" --threads "${SLURM_CPUS_PER_TASK:-8}"
rc=$?
log "audit exit $rc"
if [ -f "$A/AUDIT_OK" ]; then
  cp -f "$A/whitelist_pass2.tsv" "$D/whitelist_pass2.tsv"
  log "whitelist_pass2.tsv: $(wc -l < "$D/whitelist_pass2.tsv") fusions -> $D/whitelist_pass2.tsv"
fi
cat "$A/AUDIT.md" 2>/dev/null
# the extracted supporting reads are working files (they are human-file reads: never copied anywhere else)
rm -rf "$A/reads/"*.fq
exit $rc
