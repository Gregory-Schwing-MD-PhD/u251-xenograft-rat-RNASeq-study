# -*- coding: utf-8 -*-
"""Where the slide scripts read the pipeline's outputs: the repository's own gitignored result folders, the same paths
on the grid (/wsu/home/go/go24/go2432/u251-xenograft-murine-RNASeq-study) and in a local clone after
SLIDES/fetch_grid_outputs.sh has copied them down. Every file was traced to these paths by name and md5 on 2026-09-26.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "ANALYSIS"
COUNTS = A / "results_human_final" / "star_salmon" / "salmon.merged.gene_counts.tsv"          # nf-core/rnaseq, human + shared reads
GSEA_DIR = A / "results_therapy_v3" / "report" / "gsea" / "therapy_impact" / "combined_human"     # nf-core/differentialabundance GSEA
RLOG = A / "results_therapy_v3" / "tables" / "processed_abundance" / "all.rlog.tsv"
RUVSEQ_DIR = A / "results_ruvseq"                    # ANALYSIS/ruvseq_contamination_adjustment.R
DECONTAM_DIR = A / "results_decontamination"         # ANALYSIS/filter_contaminated_genes.R, test_contamination_bias.R
METADATA = A / "metadata_full.csv"                   # ANALYSIS/build_metadata_from_xengsort.py
GSEA_LOO = A / "gsea_leave_one_out"                  # leave-one-tumour-out GSEA results (committed: loo_sets.tsv, loo_screen.tsv, SUMMARY.md)


def gsea_file(name):
    """A file of the GSEA report folder, e.g. 'KEGG_RIBOSOME.tsv' or 'gsea_report_for_Primary_U2.tsv'."""
    return GSEA_DIR / f"therapy_impact.combined_human.{name}"


def check():
    missing = [p for p in (COUNTS, GSEA_DIR, RUVSEQ_DIR, DECONTAM_DIR, METADATA) if not p.exists()]
    if missing:
        raise SystemExit("pipeline outputs missing (run SLIDES/fetch_grid_outputs.sh in a local clone):\n  " + "\n  ".join(map(str, missing)))
