#!/bin/bash
#SBATCH --job-name=u251_geo_meth
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
# Stage the DNA methylation submission for GEO: the sixteen raw IDATs and the processed matrices, with checksums.
#
# The raw IDATs live in the laboratory's archive, not in a pipeline working directory: a Nextflow work/ tree is
# swept, and the study's only other copy was in one. This job reads the archive, never writes to it.
#
#   sbatch -D ~/logs GEO/methylation/stage_submission.sbatch
#
# Output: /rs/rs_grp_oschome/go2432/u251_geo_meth/{idat,processed,md5_checksums.txt,FILES.txt}
set -euo pipefail
ZIP=/rs/rs_grp_oschome/go2432/u251_meth/raw/DNAmethylation-20260926T220706Z-1-001.zip
OUT=/rs/rs_grp_oschome/go2432/u251_geo_meth
CHIP=205648300021
mkdir -p "$OUT/idat" "$OUT/processed"
echo "START $(date) on $(hostname)"

# A container runtime is NOT on the compute node's default PATH -- the login node has one under ~/.local/bin, the
# compute nodes do not, and the first version of this job spent two minutes staging IDATs before dying at exit 127.
# Resolve it first, and refuse immediately if there is none, so the failure costs a second instead of a stage.
echo "== 0. the container runtime"
if ! command -v singularity >/dev/null 2>&1 && ! command -v apptainer >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh 2>/dev/null || true
  module load singularity/3.5.2 2>/dev/null || module load apptainer/1.2.2 2>/dev/null || true
fi
RUNTIME=$(command -v singularity || command -v apptainer || true)
[ -n "$RUNTIME" ] || { echo "REFUSING: no singularity or apptainer on this node, and neither module loaded"; exit 1; }
echo "   using $RUNTIME ($($RUNTIME --version 2>&1 | head -1))"

echo "== 1. the sixteen raw IDATs, straight out of the laboratory archive"
cd "$OUT/idat"
unzip -j -o "$ZIP" "DNAmethylation/IDAT/${CHIP}_IDAT FILES/*.idat"
ls -la | head -20
N=$(ls -1 *.idat 2>/dev/null | wc -l)
echo "   $N IDAT files (expected 16: eight arrays x two channels)"
[ "$N" -eq 16 ] || { echo "REFUSING: expected 16 IDATs, found $N"; exit 1; }
for p in R01C01 R02C01 R03C01 R04C01 R05C01 R06C01 R07C01 R08C01; do
  for ch in Grn Red; do
    [ -f "${CHIP}_${p}_${ch}.idat" ] || { echo "REFUSING: missing ${CHIP}_${p}_${ch}.idat"; exit 1; }
  done
done
gzip -f *.idat
echo "   gzipped, as GEO serves them"

echo "== 2. the processed matrices, from those same IDATs with a pinned public tool"
CACHE=${NXF_SINGULARITY_CACHEDIR:-/rs/rs_grp_oschome/go2432/containers}
mkdir -p "$CACHE"
SIF=$CACHE/bioconductor-minfi_1.48.0.sif
if [ ! -f "$SIF" ]; then
  echo "   pulling the pinned minfi container into $CACHE"
  "$RUNTIME" pull "$SIF" docker://quay.io/biocontainers/bioconductor-minfi:1.48.0--r43hdfd78af_0
fi
cat > "$OUT/make_processed.R" <<'RSCRIPT'
suppressMessages({library(minfi); library(IlluminaHumanMethylationEPICmanifest)})
idat <- Sys.getenv("IDAT_DIR"); out <- Sys.getenv("OUT_DIR")
sheet <- read.csv(Sys.getenv("SHEET"), stringsAsFactors = FALSE)
sheet$Basename <- file.path(idat, paste0(sheet$sentrix_id, "_", sheet$sentrix_position))
rg <- read.metharray.exp(targets = sheet, force = TRUE)
sampleNames(rg) <- paste0(sheet$sample, "_EPIC")
cat("arrays read:", ncol(rg), "\n")
det <- detectionP(rg)
ms  <- preprocessNoob(rg)
b   <- getBeta(ms)
stopifnot(identical(colnames(b), colnames(det)))
write.table(data.frame(probe = rownames(b), b, check.names = FALSE),
            gzfile(file.path(out, "U251N_EPIC_beta_noob.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(probe = rownames(det), det, check.names = FALSE),
            gzfile(file.path(out, "U251N_EPIC_detection_pvalues.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
cat("probes:", nrow(b), " arrays:", ncol(b), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_processed.txt"))
RSCRIPT
export IDAT_DIR="$OUT/idat" OUT_DIR="$OUT/processed" SHEET="$SLURM_SUBMIT_DIR/array_samplesheet.csv"
[ -f "$SHEET" ] || SHEET=$HOME/u251-xenograft-murine-RNASeq-study/GEO/methylation/array_samplesheet.csv
export SHEET
"$RUNTIME" exec -B /rs -B "$HOME" "$SIF" Rscript "$OUT/make_processed.R"

echo "== 3. the multi-omics sample key travels with the deposit"
# Three of the eight pieces are named differently in the two arms (IL70B/NL70B, IL71B/NL71B, N2/N269B), and GEO
# mints no BioSample for an array sample, so nothing but this file and the sample characteristics links the arms.
KEY=$(dirname "$SHEET")/../multiomics_sample_key.tsv
[ -f "$KEY" ] || { echo "REFUSING: the sample key is missing; run python GEO/make_sample_key.py"; exit 1; }
cp "$KEY" "$OUT/processed/multiomics_sample_key.tsv"
awk -F'\t' 'NR==1{next} $10!="" && $10!=$4 {print "   two names: array "$10" = RNA library "$4}' \
  "$OUT/processed/multiomics_sample_key.tsv"

echo "== 4. checksums and the file list GEO asks for"
cd "$OUT"
md5sum idat/*.idat.gz processed/*.gz processed/multiomics_sample_key.tsv > md5_checksums.txt
{ echo "# raw"; ls -1 idat; echo "# processed"; ls -1 processed; } > FILES.txt
cat md5_checksums.txt
echo "DONE $(date)"
echo "Upload with GEO's own instructions (FTP or Aspera to the submitter's personal upload space), then fill the"
echo "series/sample spreadsheet from GEO/methylation/metadata_methylation.xlsx."
