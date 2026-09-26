#!/bin/bash
# Step 3: the September drug and subtype chain, unmodified scripts from subtypes/ copied into each arm so their
# HERE/ROOT-relative paths resolve inside the arm (arm/subtypes, arm/manuscript, arm/ANALYSIS):
#   fetch_smiles.py (ChEMBL) -> rerun_bbb.py (ADMET-AI) -> bbb_consensus.py (BOILED-Egg, B3DB) -> rebuild_drug_ranking.py
#   rerun_subtypes.py (GSVA over the published Neftel / Garofano signatures)
# then compare.py (published vs control vs held-out). The Python environment is built once on CephFS (the local
# Windows venvs of September no longer exist); its package versions are written to logs.
#   sbatch -D ANALYSIS/holdout_IL68B --dependency=afterany:<figure> ANALYSIS/holdout_IL68B/run_drugs.sbatch
#SBATCH --job-name=u251_ho_drugs
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/drugs_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
H=$U/ANALYSIS/holdout_IL68B
ENV=/rs/rs_grp_oschome/go2432/envs/u251_drug
unset PYTHONPATH LD_LIBRARY_PATH
for arm in control holdout; do
  [ -s "$H/$arm/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv" ] || { echo "step 2 output missing for $arm"; exit 3; }
done

if [ ! -x "$ENV/bin/python" ]; then
  MAMBA=$HOME/mambaforge/bin/mamba; [ -x "$MAMBA" ] || MAMBA=$HOME/mambaforge/bin/conda
  "$MAMBA" create -y -p "$ENV" -c conda-forge python=3.11 rdkit pandas numpy openpyxl scipy pip
  "$ENV/bin/pip" install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu
  "$ENV/bin/pip" install --no-cache-dir admet-ai gseapy
fi
"$ENV/bin/pip" freeze | grep -i -E "admet|chemprop|gseapy|rdkit|pandas|numpy|torch" || true
PY=$ENV/bin/python

for arm in control holdout; do
  A=$H/$arm
  mkdir -p "$A/subtypes" "$A/manuscript" "$A/ANALYSIS"
  cp "$H"/subtypes_src/*.py "$H"/subtypes_src/*.gmt "$H"/subtypes_src/counts.tsv.gz "$A/subtypes/"
  cp "$A/publication_figure/therapy_impact_Drug_Profiles_Comprehensive.csv" "$A/manuscript/"
  if [ "$arm" = control ]; then cp "$U/ANALYSIS/metadata_therapy.csv" "$A/ANALYSIS/"; else cp "$H/ANALYSIS/metadata_therapy.csv" "$A/ANALYSIS/"; fi
  echo "=== $arm: $(wc -l < "$A/manuscript/therapy_impact_Drug_Profiles_Comprehensive.csv") drug rows; samples $(tail -n +2 "$A/ANALYSIS/metadata_therapy.csv" | cut -d, -f1 | tr '\n' ' ')"
  ( cd "$A/subtypes" && "$PY" fetch_smiles.py && "$PY" rerun_bbb.py && "$PY" bbb_consensus.py && "$PY" rebuild_drug_ranking.py \
      && "$PY" rerun_subtypes.py counts.tsv.gz )
  ls -la "$A/subtypes/drug_ranking_final.csv" "$A/subtypes/subtype_rerun_results.csv"
done

"$PY" "$H/compare.py"
echo "drugs done"
