#!/bin/bash
# Step 2: create_publication_figure_600_dpi.R (DSigDB drug enrichment, ChEMBL lookups through the shared
# .drug_discovery_cache, STRING network, Supplementary_Data.xlsx) in two arms:
#   control  the published six-tumour inputs (ANALYSIS/results_therapy_v3) - must reproduce the published
#            publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv, which checks that today's container
#            and ChEMBL responses still give May's answer
#   holdout  the IL68B-held-out results of step 1
#   sbatch -D ANALYSIS/holdout_IL68B --dependency=afterany:<de> ANALYSIS/holdout_IL68B/run_figure.sbatch
#SBATCH --job-name=u251_ho_fig
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/figure_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=$U/ANALYSIS/holdout_IL68B
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
unset JAVA_HOME LD_LIBRARY_PATH PYTHONPATH R_LIBS R_LIBS_USER R_LIBS_SITE
export XDG_RUNTIME_DIR="${HOME}/xdr"; mkdir -p "$XDG_RUNTIME_DIR"
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
CACHE=/rs/rs_grp_oschome/go2432/singularity_cache
IMG=$CACHE/go2432-bioconductor.sif
[ -s "$H/results_therapy_v3_noIL68B/tables/differential/therapy_impact.deseq2.results.tsv" ] || { echo "step 1 output missing"; exit 3; }
if [ ! -s "$IMG" ]; then "$SING" pull "$IMG" docker://go2432/bioconductor:latest; fi
md5sum "$IMG"; "$SING" inspect "$IMG" 2>/dev/null | head -20 || true
cd "$U"                                   # the R script's ChEMBL cache is ./.drug_discovery_cache, shared with May's run
ls .drug_discovery_cache | wc -l

run_arm () {  # arm results_dir metadata
  local arm=$1 res=$2 meta=$3
  mkdir -p "$H/$arm/publication_figure"
  "$SING" exec --bind "$U:/data" --pwd /data "$IMG" \
    Rscript /data/create_publication_figure_600_dpi.R \
      "/data/${res}/tables/processed_abundance/all.rlog.tsv" \
      "/data/${res}" \
      "/data/ANALYSIS/refs/pathways/human" \
      "/data/ANALYSIS/refs/pathways/human/human_string" \
      "/data/ANALYSIS/holdout_IL68B/${arm}/publication_figure" \
      therapy_impact \
      "/data/${meta}" \
      ""
  ls -la "$H/$arm/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv" "$H/$arm/publication_figure/Supplementary_Data.xlsx"
}
run_arm control ANALYSIS/results_therapy_v3 ANALYSIS/metadata_therapy.csv
run_arm holdout ANALYSIS/holdout_IL68B/results_therapy_v3_noIL68B ANALYSIS/holdout_IL68B/ANALYSIS/metadata_therapy.csv
echo "--- control Drug_Profiles against the published one"
cmp "$H/control/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv" "$U/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv" && echo "IDENTICAL" || echo "DIFFERS (compare.py quantifies)"
echo "figure done"
