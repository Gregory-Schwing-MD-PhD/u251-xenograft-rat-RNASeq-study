#!/bin/bash
# Prep for the leave-one-tumour-out GSEA: the same GSEA 4.3.2 image the nf-core run used (re-pulled; the old
# singularity cache was deleted), the pipeline's own GSEA inputs copied from its work directory, and the per-condition
# matrices. Submit with:  sbatch -D ANALYSIS/gsea_leave_one_out ANALYSIS/gsea_leave_one_out/prep.sbatch
#SBATCH --job-name=u251_gsea_prep
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/prep_%j.out
set -euo pipefail
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
D=$U/ANALYSIS/gsea_leave_one_out
W=$U/work/85/373d2dd0f9132fb3c16154e6a6bded           # the nf-core GSEA_GSEA task of the therapy_v3 run
CACHE=/rs/rs_grp_oschome/go2432/singularity_cache
IMG=$CACHE/depot.galaxyproject.org-singularity-gsea-4.3.2--hdfd78af_0.img
SING=$(command -v singularity || echo /wsu/home/go/go24/go2432/.local/bin/singularity)
PY=/wsu/home/go/go24/go2432/mambaforge/bin/python

mkdir -p "$CACHE" "$D/inputs/original" "$D/logs"
if [ ! -s "$IMG" ]; then
  curl -fsSL --retry 5 -o "$IMG.part" https://depot.galaxyproject.org/singularity/gsea:4.3.2--hdfd78af_0
  mv "$IMG.part" "$IMG"
fi
md5sum "$IMG"
"$SING" exec "$IMG" bash -c 'command -v gsea-cli; ls "$(dirname "$(readlink -f "$(command -v gsea-cli)")")" | head -5'

cp -L "$W/therapy_impact.combined_human.therapy_impact.gct" "$D/inputs/therapy_impact.gct"
cp -L "$W/therapy_impact.combined_human.therapy_impact.cls" "$D/inputs/therapy_impact.cls"
cp -L "$W/therapy_impact.combined_human.combined_human.gmt" "$D/inputs/combined_human.gmt"
cp -L "$W/therapy_impact.combined_human.GRCh38.anno.feature_metadata.chip" "$D/inputs/GRCh38.anno.feature_metadata.chip"
cp -L "$W/therapy_impact.combined_human.gsea_report_for_Primary_U2.tsv" "$W/therapy_impact.combined_human.gsea_report_for_Recurrent_U2.tsv" \
      "$W/therapy_impact.combined_human.Gsea.rpt" "$W/.command.sh" "$D/inputs/original/"
md5sum "$D"/inputs/*.* "$D"/inputs/original/*

# smoke test: the full command on the pipeline's own matrix with 10 permutations. It runs BEFORE tasks.tsv is written,
# so if the container or the command fails the array tasks find no tasks.tsv and stop at once.
rm -rf "$D/smoke"; mkdir -p "$D/smoke/tmp"
export SINGULARITYENV_JAVA_TOOL_OPTIONS="-Duser.home=$D/smoke -Djava.io.tmpdir=$D/smoke/tmp"
( cd "$D/smoke" && "$SING" exec -B /wsu/home/go/go24/go2432 "$IMG" gsea-cli GSEA \
    -res "$D/inputs/therapy_impact.gct" -cls "$D/inputs/therapy_impact.cls#Recurrent_U2_versus_Primary_U2" \
    -gmx "$D/inputs/combined_human.gmt" -chip "$D/inputs/GRCh38.anno.feature_metadata.chip" -collapse true \
    -out "$D/smoke" --rpt_label therapy_impact.combined_human -nperm 10 -permute gene_set -scoring_scheme weighted \
    -metric Diff_of_Classes -sort real -order descending -set_max 500 -set_min 15 -norm meandiv -rnd_type no_balance \
    -make_sets false -median false -num 100 -plot_top_x 0 -rnd_seed 1234 -save_rnd_lists false -zip_report false )
ls "$D"/smoke/therapy_impact.combined_human.Gsea.*/gsea_report_for_Primary_U2_*.tsv
unset SINGULARITYENV_JAVA_TOOL_OPTIONS

"$PY" "$D/make_inputs.py"
echo "prep done"
