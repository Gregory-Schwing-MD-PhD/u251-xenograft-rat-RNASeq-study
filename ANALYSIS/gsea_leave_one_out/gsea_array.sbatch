#!/bin/bash
# One GSEA run per line of tasks.tsv, with the nf-core run's exact command (GSEA 4.3.2, gene-set permutation,
# Diff_of_Classes, 1,000 permutations, weighted scoring, meandiv, set size 15-500, chip collapse) and only the matrix,
# class file and seed changed. 'light' runs skip the per-set pages and plots (-make_sets false -plot_top_x 0), which
# does not change any statistic. Submit after prep with:
#   sbatch -D ANALYSIS/gsea_leave_one_out --dependency=afterany:<prep> ANALYSIS/gsea_leave_one_out/gsea_array.sbatch
#SBATCH --job-name=u251_gsea_loo
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --array=0-40
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/gsea_%A_%a.out
set -euo pipefail
D=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study/ANALYSIS/gsea_leave_one_out
IMG=/rs/rs_grp_oschome/go2432/singularity_cache/depot.galaxyproject.org-singularity-gsea-4.3.2--hdfd78af_0.img
SING=$(command -v singularity || echo /wsu/home/go/go24/go2432/.local/bin/singularity)
[ -s "$IMG" ] || { echo "no GSEA image at $IMG (prep failed?)"; exit 2; }
[ -s "$D/tasks.tsv" ] || { echo "no tasks.tsv (prep failed?)"; exit 3; }

line=$(awk -F'\t' -v t="$SLURM_ARRAY_TASK_ID" 'NR > 1 && $1 == t' "$D/tasks.tsv")
[ -n "$line" ] || { echo "no task $SLURM_ARRAY_TASK_ID"; exit 4; }
cond=$(cut -f2 <<< "$line"); seed=$(cut -f3 <<< "$line"); mode=$(cut -f4 <<< "$line")
echo "task $SLURM_ARRAY_TASK_ID: condition $cond seed $seed output $mode on $(hostname)"
OUT=$D/runs/${cond}_s${seed}
rm -rf "$OUT"                                  # a requeued task starts clean
mkdir -p "$OUT/tmp"
cd "$OUT"
if [ "$mode" = full ]; then EXTRA="-make_sets true -plot_top_x 20"; else EXTRA="-make_sets false -plot_top_x 0"; fi
export SINGULARITYENV_JAVA_TOOL_OPTIONS="-Duser.home=$OUT -Djava.io.tmpdir=$OUT/tmp"   # no shared ~/gsea_home between tasks

"$SING" exec -B /wsu/home/go/go24/go2432 "$IMG" gsea-cli GSEA \
    -res "$D/inputs/$cond/therapy_impact.gct" \
    -cls "$D/inputs/$cond/therapy_impact.cls#Recurrent_U2_versus_Primary_U2" \
    -gmx "$D/inputs/combined_human.gmt" \
    -chip "$D/inputs/GRCh38.anno.feature_metadata.chip" -collapse true \
    -out "$OUT" \
    --rpt_label therapy_impact.combined_human \
    -nperm 1000 -permute gene_set -scoring_scheme weighted -metric Diff_of_Classes -sort real -order descending \
    -set_max 500 -set_min 15 -norm meandiv -rnd_type no_balance $EXTRA -median false -num 100 \
    -rnd_seed "$seed" -save_rnd_lists false -zip_report false

ls "$OUT"/therapy_impact.combined_human.Gsea.*/gsea_report_for_Primary_U2_*.tsv
rm -rf "$OUT/tmp"
echo "task $SLURM_ARRAY_TASK_ID done"
