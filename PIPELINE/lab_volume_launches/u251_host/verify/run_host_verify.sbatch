#!/bin/bash
# Host (rat) DE verification, Recurrent_U2 vs Primary_U2 (deck backup B7). Runs from the grid verify directory, never
# from the repository clone.
#   step A  host_verify_de.R in the pipeline's own DESeq2 container (1.34.0 + ashr 2.2-54): R0 reproduce the table,
#           R1 leave one tumour out, R2 ten 3-v-3 relabellings, R3 NL-prefix contrast + IL66B alone, R5 human-read
#           leakage; writes each fit's Wald statistic to results/ranks/
#   step B  host_verify_gsea.R in the lab's Bioconductor image: fgsea (10,000 permutations, sets 15-500) per fit and
#           collection; installs fgsea into /rs/rs_grp_oschome/go2432/Rlib/R<x.y> if the image lacks it. If step B
#           fails, Broad GSEAPreranked (the GSEA 4.3.2 image) runs on the full and leave-one-out fits instead.
#   sbatch -D /rs/rs_grp_oschome/go2432/u251_host/verify /rs/rs_grp_oschome/go2432/u251_host/verify/run_host_verify.sbatch
#SBATCH --job-name=u251_hostverify
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/verify_%j.out
set -euo pipefail
V=/rs/rs_grp_oschome/go2432/u251_host/verify
H=/rs/rs_grp_oschome/go2432/u251_host
U=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study
C=/rs/rs_grp_oschome/go2432/singularity_cache
DESEQ_IMG=$C/depot.galaxyproject.org-singularity-bioconductor-deseq2-1.34.0--r41hc247a5b_3.img
SIF=$C/go2432-bioconductor.sif
GSEA_IMG=$C/depot.galaxyproject.org-singularity-gsea-4.3.2--hdfd78af_0.img
SING=$HOME/mambaforge/envs/nextflow/bin/singularity
OUT=$V/results
NCPU=${SLURM_CPUS_PER_TASK:-4}
unset R_LIBS R_LIBS_USER R_LIBS_SITE PYTHONPATH LD_LIBRARY_PATH JAVA_HOME
export XDG_RUNTIME_DIR=/dev/shm/${USER}_xdg_${SLURM_JOB_ID}; mkdir -p "$XDG_RUNTIME_DIR"
export TMPDIR=/dev/shm/${USER}_tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
export SINGULARITYENV_TMPDIR=$TMPDIR
trap 'rm -rf "$XDG_RUNTIME_DIR" "$TMPDIR"' EXIT
mkdir -p "$OUT"
cd "$V"
B=(--bind "$U:$U" --bind /rs/rs_grp_oschome/go2432:/rs/rs_grp_oschome/go2432)
echo "host $(hostname), job $SLURM_JOB_ID, $(date)"
"$SING" --version

echo "== step A: DESeq2 fits in the pipeline's container, $(date)"
"$SING" exec "${B[@]}" "$DESEQ_IMG" Rscript -e 'invisible(parse("host_verify_de.R")); cat("parse ok\n")'
"$SING" exec "${B[@]}" "$DESEQ_IMG" Rscript host_verify_de.R "$OUT" \
  "$H/results_rnaseq/star_salmon/salmon.merged.gene_counts.tsv" \
  "$H/results_rnaseq/star_salmon/salmon.merged.gene_lengths.tsv" \
  "$H/results_de/tables/differential/host_therapy.deseq2.results.tsv" \
  "$U/ANALYSIS/host/refs_online/rat_human_1to1.tsv" \
  "$U/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
  "$NCPU" ashr

echo "== step B: fgsea in the Bioconductor image, $(date)"
if "$SING" exec "${B[@]}" "$SIF" Rscript host_verify_gsea.R "$OUT" "$U/ANALYSIS/refs/pathways/rat" \
     /rs/rs_grp_oschome/go2432/Rlib "$NCPU" "$V/broad_gsea_existing"; then
  echo "step B ok"
else
  echo "step B (fgsea) failed: Broad GSEAPreranked fallback on the full and leave-one-out fits, $(date)"
  F=$OUT/gsea_preranked_fallback; mkdir -p "$F/tmp"
  export SINGULARITYENV_JAVA_TOOL_OPTIONS="-Duser.home=$F -Djava.io.tmpdir=$F/tmp"
  for fitname in full loo_IL67B loo_IL68B loo_IL69B loo_IL66B loo_NL70B loo_NL71B; do
    # one row per symbol, the Ensembl id with the largest |stat|
    awk -F'\t' 'NR > 1 && $2 != "" && $2 != "NA" { a = ($3 < 0 ? -$3 : $3); print a "\t" $2 "\t" $3 }' "$OUT/ranks/$fitname.tsv" \
      | sort -t$'\t' -k1,1gr | awk -F'\t' '!seen[$2]++ { print $2 "\t" $3 }' > "$F/$fitname.rnk"
    for co in hallmark:rn.hallmark.v2023.2.gmt kegg:rn.kegg.v2023.2.gmt gobp:rn.gobp.v2023.2.gmt brain:rn.brain_gmt_v2.gmt; do
      "$SING" exec "${B[@]}" "$GSEA_IMG" gsea-cli GSEAPreranked -rnk "$F/$fitname.rnk" \
        -gmx "$U/ANALYSIS/refs/pathways/rat/${co#*:}" -collapse No_Collapse -mode Max_probe -norm meandiv \
        -nperm 10000 -rnd_seed 1234 -scoring_scheme weighted -set_max 500 -set_min 15 -make_sets false \
        -plot_top_x 0 -create_svgs false -include_only_symbols true -zip_report false \
        -rpt_label "${fitname}.${co%%:*}" -out "$F" || echo "GSEAPreranked failed: $fitname ${co%%:*}"
    done
  done
  rm -rf "$F/tmp"
fi
ls -la "$OUT"
echo "done $(date)"
