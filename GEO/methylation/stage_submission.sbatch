#!/bin/bash
#SBATCH --job-name=u251_geo_meth
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
# Stage the DNA methylation submission for GEO: the sixteen raw IDATs, the processed matrices, the multi-omics sample
# key, and checksums for all of it.
#
# The raw IDATs live in the laboratory's archive, not in a pipeline working directory: a Nextflow work/ tree is swept,
# and the study's only other copy was in one. This job reads the archive, never writes to it.
#
#   sbatch -D <repo>/GEO/methylation <repo>/GEO/methylation/stage_submission.sbatch
#
# Output: /rs/rs_grp_oschome/go2432/u251_geo_meth/{idat,processed,md5_checksums.txt,FILES.txt}
#
# THE CONTAINER ROUTE IS NOT NEGOTIABLE ON THIS CLUSTER, and two jobs were lost learning it:
#   * 40483631 died at exit 127 -- no container runtime on a compute node's default PATH.
#   * 40483632 found ~/.local/bin/singularity, which is a partial install: it pulls, then dies converting the SIF to
#     a sandbox with "/usr/libexec/singularity/bin/starter not found".
# The route that works, and that every other R step in this study already uses, is the conda singularity plus a
# PRE-UNPACKED sandbox (no SIF conversion at run time) plus the project's own R library, which already holds minfi
# and the EPIC manifest. Using it also means the deposited matrix is built with the same minfi as the analyses.
set -euo pipefail

M=/rs/rs_grp_oschome/go2432/u251_meth
OUT=/rs/rs_grp_oschome/go2432/u251_geo_meth
ZIP=$M/raw/DNAmethylation-20260926T220706Z-1-001.zip
ZIP_SHA=1107d481844de343f398f3ccd22414a4502c2877b13eda3341e17277390a9fc2
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
SANDBOX=$M/containers/bioconductor_docker_RELEASE_3_22-R-4.5.2.sandbox
RLIB=$M/rlib/bioc-3.22_R-4.5.2
CHIP=205648300021
TMP=$OUT/tmp
mkdir -p "$OUT/idat" "$OUT/processed" "$TMP"
echo "START $(date) on $(hostname)"

echo "== 0. the container route, checked before anything is staged"
# Refuse in seconds rather than after unpacking 16 IDATs, and never fall back to a runtime on PATH: on this cluster
# the one on PATH is the broken one.
for p in "$SING" "$SANDBOX" "$RLIB" "$ZIP"; do
  [ -e "$p" ] || { echo "REFUSING: missing $p"; exit 1; }
done
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH
export SINGULARITY_TMPDIR=$TMP APPTAINER_TMPDIR=$TMP TMPDIR=$TMP
# Bind the SPECIFIC group path, not /rs: `-B /rs` does not make the tree visible inside on this cluster, and job
# 40483633 died with "cannot open file .../tmp/probe.R" although the file was there on the host. Every other R step in
# this study binds /rs/rs_grp_oschome/go2432. /tmp must also be writable, or even Rscript -e cannot open its temp file.
SB=(-B /rs/rs_grp_oschome/go2432 -B "$HOME" -B "$TMP:/tmp")
echo 'cat(R.version.string, "\n"); .libPaths(c(Sys.getenv("U251_RLIB"), .libPaths()))
for (p in c("minfi", "IlluminaHumanMethylationEPICmanifest")) cat(" ", p, as.character(packageVersion(p)), "\n")' \
  > "$TMP/probe.R"
"$SING" exec "${SB[@]}" --env "U251_RLIB=$RLIB" "$SANDBOX" Rscript /tmp/probe.R
echo "   the route works"

echo "== 1. the source archive, checked against its recorded digest"
got=$(sha256sum "$ZIP" | cut -d' ' -f1)
[ "$got" = "$ZIP_SHA" ] || { echo "REFUSING: zip sha256 $got != $ZIP_SHA"; exit 1; }
echo "   sha256 ok: $ZIP_SHA"

echo "== 2. the sixteen raw IDATs, straight out of the laboratory archive"
cd "$OUT/idat"
rm -f ./*.idat ./*.idat.gz          # a requeue must not mix a previous attempt's files with this one's
unzip -j -o "$ZIP" "DNAmethylation/IDAT/${CHIP}_IDAT FILES/*.idat"
N=$(ls -1 ./*.idat 2>/dev/null | wc -l)
echo "   $N IDAT files (expected 16: eight arrays x two channels)"
[ "$N" -eq 16 ] || { echo "REFUSING: expected 16 IDATs, found $N"; exit 1; }
for p in R01C01 R02C01 R03C01 R04C01 R05C01 R06C01 R07C01 R08C01; do
  for ch in Grn Red; do
    [ -f "${CHIP}_${p}_${ch}.idat" ] || { echo "REFUSING: missing ${CHIP}_${p}_${ch}.idat"; exit 1; }
  done
done

echo "== 3. the processed matrices, from those same IDATs"
cat > "$TMP/make_processed.R" <<'RSCRIPT'
.libPaths(c(Sys.getenv("U251_RLIB"), .libPaths()))
suppressMessages({library(minfi); library(IlluminaHumanMethylationEPICmanifest)})
idat <- Sys.getenv("IDAT_DIR"); out <- Sys.getenv("OUT_DIR")
sheet <- read.csv(Sys.getenv("SHEET"), stringsAsFactors = FALSE)
sheet$Basename <- file.path(idat, paste0(sheet$sentrix_id, "_", sheet$sentrix_position))
rg <- read.metharray.exp(targets = sheet, force = TRUE)
# The column names are the array names, unchanged: IL66B .. IL71B, N2, C2B. Three of them differ from the matched RNA
# library's name (IL70B/NL70B, IL71B/NL71B, N2/N269B). That is deliberate and the pairing travels in the sample key.
sampleNames(rg) <- paste0(sheet$sample, "_EPIC")
cat("arrays read:", ncol(rg), "\n")
stopifnot(ncol(rg) == nrow(sheet))
det <- detectionP(rg)
ms  <- preprocessNoob(rg)
b   <- getBeta(ms)
stopifnot(identical(colnames(b), colnames(det)))
gz <- function(f) gzfile(file.path(out, f))
write.table(data.frame(probe = rownames(b), b, check.names = FALSE),
            gz("U251N_EPIC_beta_noob.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(probe = rownames(det), det, check.names = FALSE),
            gz("U251N_EPIC_detection_pvalues.tsv.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("probes:", nrow(b), " arrays:", ncol(b), "\n")
cat("columns:", paste(colnames(b), collapse = ", "), "\n")
cat("minfi:", as.character(packageVersion("minfi")), "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_processed.txt"))
RSCRIPT
SHEET=$SLURM_SUBMIT_DIR/array_samplesheet.csv
[ -f "$SHEET" ] || SHEET=$HOME/u251-xenograft-murine-RNASeq-study/GEO/methylation/array_samplesheet.csv
[ -f "$SHEET" ] || { echo "REFUSING: cannot find array_samplesheet.csv"; exit 1; }
cp "$SHEET" "$TMP/sheet.csv"          # the sheet may live under $HOME or the submit dir; /tmp is always bound
"$SING" exec "${SB[@]}" \
  --env "U251_RLIB=$RLIB,IDAT_DIR=$OUT/idat,OUT_DIR=$OUT/processed,SHEET=/tmp/sheet.csv" \
  "$SANDBOX" Rscript /tmp/make_processed.R

echo "== 4. the multi-omics sample key travels with the deposit"
# Three of the eight pieces are named differently in the two arms, and GEO mints no BioSample for an array sample, so
# nothing but this file and the sample characteristics links the two assays.
KEY=$(dirname "$SHEET")/../multiomics_sample_key.tsv
[ -f "$KEY" ] || { echo "REFUSING: the sample key is missing; run python GEO/make_sample_key.py"; exit 1; }
cp "$KEY" "$OUT/processed/multiomics_sample_key.tsv"
awk -F'\t' 'NR>1 && $10!="" && $10!=$4 {print "   two names: array "$10" = RNA library "$4}' \
  "$OUT/processed/multiomics_sample_key.tsv"

echo "== 5. gzip the IDATs, as GEO serves them, and checksum everything"
cd "$OUT/idat" && gzip -f ./*.idat
cd "$OUT"
md5sum idat/*.idat.gz processed/* > md5_checksums.txt
{ echo "# raw"; ls -1 idat; echo "# processed"; ls -1 processed; } > FILES.txt
cat md5_checksums.txt
echo "DONE $(date)"
echo "Upload with GEO's own instructions (FTP or Aspera to the submitter's personal upload space), then fill the"
echo "series/sample spreadsheet from GEO/methylation/metadata_methylation.xlsx and send GEO/methylation/COVER_NOTE.md."
