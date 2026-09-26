#!/bin/bash
# Collect the leave-one-tumour-out GSEA runs. Submit with:
#   sbatch -D ANALYSIS/gsea_leave_one_out --dependency=afterany:<array> ANALYSIS/gsea_leave_one_out/collect.sbatch
#SBATCH --job-name=u251_gsea_collect
#SBATCH --partition=reqp
#SBATCH --qos=requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/collect_%j.out
set -euo pipefail
D=/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study/ANALYSIS/gsea_leave_one_out
/wsu/home/go/go24/go2432/mambaforge/bin/python "$D/collect.py"
du -sh "$D/runs"
