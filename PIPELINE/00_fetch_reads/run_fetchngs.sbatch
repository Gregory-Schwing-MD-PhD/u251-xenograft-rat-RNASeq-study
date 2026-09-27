#!/bin/bash
#SBATCH --job-name=u251_fetchngs
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
# Step 00: get the reads from the public accession.
#
# This replaces the repository's previous route to the reads, which was a private Google Drive folder behind a
# personal OAuth token: nobody outside this laboratory could run it, and the accession GSE338105 appeared nowhere in
# the README. Everything here is a published tool driven by a public accession.
#
#   nf-core/fetchngs resolves GSE338105 -> its GSM samples -> the SRA runs -> the FASTQ files, downloads them,
#   checks each against the md5 the archive publishes, and writes a samplesheet already in nf-core/rnaseq's format.
#
# Pinned on purpose: -r 1.12.0 is a released revision, not a branch. Change it deliberately, never silently.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-/rs/rs_grp_oschome/go2432/u251_reads}
REV=1.12.0

module load singularity 2>/dev/null || true
export NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-/rs/rs_grp_oschome/go2432/containers}
export NXF_HOME=${NXF_HOME:-/rs/rs_grp_oschome/go2432/nxf_home}
mkdir -p "$OUT" "$NXF_SINGULARITY_CACHEDIR" "$NXF_HOME"

echo "START $(date) on $(hostname); fetchngs $REV -> $OUT"
nextflow run nf-core/fetchngs -r "$REV" \
    --input "$HERE/ids.csv" \
    --outdir "$OUT" \
    --nf_core_pipeline rnaseq \
    --download_method sratools \
    -profile singularity \
    -resume

echo "== what came back"
ls -la "$OUT"/fastq 2>/dev/null | head
echo "== the samplesheet nf-core/rnaseq will consume"
find "$OUT" -name "samplesheet*.csv" -maxdepth 3 -print -exec head -3 {} \;

echo "== check the run accessions against the map this project uses"
python3 "$HERE/check_library_map.py" --samplesheet "$(find "$OUT" -name 'samplesheet*.csv' | head -1)" \
                                     --map "$HERE/library_map.csv"
echo "DONE $(date)"
