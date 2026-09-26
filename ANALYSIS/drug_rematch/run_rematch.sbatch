#!/bin/bash
# The drug chain again with exact-name ChEMBL matching (resolve_chembl.py apply chembl_mapping.csv), for the published
# six tumours and both hold-outs. Everything after the matching is the September chain, unmodified:
#   fetch_smiles.py (SMILES by the corrected ChEMBL ID) -> rerun_bbb.py (ADMET-AI) -> bbb_consensus.py (BOILED-Egg, B3DB)
#   -> rebuild_drug_ranking.py (clinical phase 1-4 and a structure; |NES|^1.5 x BBB)
# then compare_rematch.py (old matching vs exact matching, per arm) -> REMATCH.md, rematch.json.
#   sbatch -D ANALYSIS/drug_rematch ANALYSIS/drug_rematch/run_rematch.sbatch
#SBATCH --job-name=u251_rematch
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=logs/rematch_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
R=$U/ANALYSIS/drug_rematch
PY=/rs/rs_grp_oschome/go2432/envs/u251_drug/bin/python
SRC=$U/ANALYSIS/holdout_IL68B/subtypes_src                 # the September subtypes/ scripts, as the hold-outs ran them
unset PYTHONPATH LD_LIBRARY_PATH
"$PY" -c "import admet_ai, pandas, rdkit; print('env ok')"
for arm in published IL68B IL66B; do
  A=$R/runs/$arm; rm -rf "$A"; mkdir -p "$A/subtypes" "$A/manuscript" "$A/ANALYSIS"
  cp "$SRC"/fetch_smiles.py "$SRC"/rerun_bbb.py "$SRC"/bbb_consensus.py "$SRC"/rebuild_drug_ranking.py "$A/subtypes/"
  if [ "$arm" = published ]; then IN=$R/inputs/published_drug_profiles.csv; else IN=$R/inputs/${arm}_holdout_drug_profiles.csv; fi
  "$PY" "$R/resolve_chembl.py" apply "$R/chembl_mapping.csv" "$IN" "$A/manuscript/therapy_impact_Drug_Profiles_Comprehensive.csv"
  echo "=== $arm"
  ( cd "$A/subtypes" && "$PY" fetch_smiles.py && "$PY" rerun_bbb.py && "$PY" bbb_consensus.py && "$PY" rebuild_drug_ranking.py ) \
    > "$A/chain.log" 2>&1 || { echo "chain failed for $arm"; tail -30 "$A/chain.log"; exit 4; }
  tail -25 "$A/chain.log"
done
"$PY" "$R/compare_rematch.py"
echo "rematch done"
