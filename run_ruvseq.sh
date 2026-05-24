#!/bin/bash
# ==============================================================================
# RUVSeq CONTAMINATION ADJUSTMENT - SLURM SUBMISSION (Wayne State Warrior HPC)
# ==============================================================================
# Estimates RUVSeq factors of unwanted variation from the Control (rat-brain)
# samples and re-runs the Primary-vs-Recurrent DESeq2 contrast adjusted for
# them. CPU/memory-bound (no GPU). Mirrors the R-via-Singularity invocation in
# run_publication_figure.sh.
# ==============================================================================
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --job-name=ruvseq_contam
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=slurm_logs/ruvseq_%j.out
#SBATCH --error=slurm_logs/ruvseq_%j.err

set -euo pipefail

echo "================================================================"
echo " RUVSeq contamination adjustment"
echo " Job ID : ${SLURM_JOB_ID:-local}   Node: $(hostname)"
echo " Start  : $(date)"
echo "================================================================"

# ------------------------------------------------------------------------------
# Environment (matches run_publication_figure.sh)
# ------------------------------------------------------------------------------
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
unset JAVA_HOME

export XDG_RUNTIME_DIR="${HOME}/xdr"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p "$XDG_RUNTIME_DIR" "$NXF_SINGULARITY_CACHEDIR"

export NXF_SINGULARITY_HOME_MOUNT=true
unset LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE

# ------------------------------------------------------------------------------
# Paths (relative to the repo root = $PWD, bound to /data inside the container)
# ------------------------------------------------------------------------------
R_SCRIPT="ANALYSIS/ruvseq_contamination_adjustment.R"
COUNTS="ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
GEN_SCRIPT="ANALYSIS/build_metadata_from_xengsort.py"
BASE_METADATA="ANALYSIS/metadata_base.csv"
XENGSORT_DIR="ANALYSIS/xengsort_out"
METADATA="ANALYSIS/metadata_full.csv"

# Optional positional args (to compare control sets):
#   $1 = comma-separated RUV anchor control IDs (e.g. "N168B,N269B" to drop the
#        failed-graft IL64B). Empty -> all Classification==Control samples.
#   $2 = output dir (use a new one to keep both runs, e.g. ANALYSIS/results_ruvseq_2ctrl)
RUV_CONTROLS="${1:-}"
OUT_DIR="${2:-ANALYSIS/results_ruvseq}"

mkdir -p "$OUT_DIR" slurm_logs

# ------------------------------------------------------------------------------
# Intermediate step: derive contamination fractions from the xengsort logs and
# append them to the design metadata (so they are never hand-maintained).
# ------------------------------------------------------------------------------
if [ -d "$XENGSORT_DIR" ]; then
    echo "Building $METADATA from xengsort logs in $XENGSORT_DIR ..."
    python3 "$GEN_SCRIPT" \
        --base "$BASE_METADATA" \
        --xengsort-dir "$XENGSORT_DIR" \
        --out "$METADATA"
elif [ -f "$METADATA" ]; then
    echo "WARNING: $XENGSORT_DIR not found; using existing $METADATA as-is."
else
    echo "ERROR: neither $XENGSORT_DIR nor $METADATA is present; cannot build metadata."
    exit 1
fi

echo "Verifying input files..."
for path in "$R_SCRIPT" "$COUNTS" "$METADATA"; do
    if [ ! -e "$path" ]; then
        echo "ERROR: Required path not found: $path"
        exit 1
    fi
done
echo "All required files found."

# ------------------------------------------------------------------------------
# Container (same Bioconductor image as run_publication_figure.sh)
# ------------------------------------------------------------------------------
IMG_PATH="${NXF_SINGULARITY_CACHEDIR}/go2432-bioconductor.sif"
if [[ ! -f "$IMG_PATH" ]]; then
    echo "Pulling Bioconductor container..."
    singularity pull "$IMG_PATH" docker://go2432/bioconductor:latest
fi
echo "Container ready: $IMG_PATH"

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
echo "Launching RUVSeq adjustment..."
echo "  R script : $R_SCRIPT"
echo "  Counts   : $COUNTS"
echo "  Metadata : $METADATA"
echo "  Out dir  : $OUT_DIR"
echo "  Controls : ${RUV_CONTROLS:-<all Classification==Control>}"

set +e
singularity exec --bind "$PWD:/data" --pwd /data "$IMG_PATH" \
    Rscript "/data/$R_SCRIPT" \
    "/data/$COUNTS" \
    "/data/$METADATA" \
    "/data/$OUT_DIR" \
    "$RUV_CONTROLS"
exit_code=$?
set -e

echo "================================================================"
echo " End      : $(date)"
echo " Exit code: $exit_code"
echo "================================================================"

if [ $exit_code -ne 0 ]; then
    echo "ERROR: RUVSeq adjustment failed (exit code: $exit_code)"
    exit $exit_code
fi

echo "Outputs in $OUT_DIR:"
ls -lh "$OUT_DIR"
echo "Recommended k: $(cat "$OUT_DIR/ruvseq_recommended_k.txt" 2>/dev/null || echo 'N/A')"
exit 0
