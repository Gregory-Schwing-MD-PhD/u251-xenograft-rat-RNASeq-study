#!/usr/bin/env Rscript
# ==============================================================================
# Compare DESeq2 DE result tables across decontamination approaches
# ==============================================================================
# Quantifies how much the Recurrent-vs-Primary call changes between approaches
# (baseline / hard contamination filter / RUVSeq adjustment). The FIRST table
# is the reference for the concordance columns.
#
# Each table is a TSV with at least: gene_id, log2FoldChange, padj
# (DESeq2 / nf-core differentialabundance results, or ruvseq_*_de.tsv).
#
# Usage:
#   Rscript compare_de_tables.R <out.csv> <label1>=<path1> <label2>=<path2> ...
#
# Example:
#   Rscript ANALYSIS/compare_de_tables.R ANALYSIS/de_comparison.csv \
#     baseline=ANALYSIS/results_therapy_v3_baseline/tables/differential/therapy_impact.deseq2.results.tsv \
#     hardfilter=ANALYSIS/results_therapy_v3/tables/differential/therapy_impact.deseq2.results.tsv \
#     ruvseq_k2=ANALYSIS/results_ruvseq/ruvseq_adjusted_de_k2.tsv
# ==============================================================================

PADJ <- 0.05
LFC  <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("Usage: Rscript compare_de_tables.R <out.csv> label1=path1 label2=path2 [...]\n")
    quit(status = 1)
}
OUT   <- args[1]
specs <- args[-1]
labels <- sub("=.*$", "", specs)
paths  <- sub("^[^=]*=", "", specs)

read_de <- function(path) {
    if (!file.exists(path)) stop("DE table not found: ", path)
    df <- read.delim(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    gid <- if ("gene_id" %in% colnames(df)) df$gene_id else as.character(df[[1]])
    if (!"log2FoldChange" %in% colnames(df) || !"padj" %in% colnames(df))
        stop(path, ": needs 'log2FoldChange' and 'padj' columns")
    data.frame(gene_id = as.character(gid),
               log2FoldChange = suppressWarnings(as.numeric(df$log2FoldChange)),
               padj = suppressWarnings(as.numeric(df$padj)),
               stringsAsFactors = FALSE)
}
sig_set <- function(d) d$gene_id[!is.na(d$padj) & d$padj < PADJ & abs(d$log2FoldChange) > LFC]
top_ids <- function(d, n = 20) {
    x <- d[!is.na(d$padj), ]
    if (!nrow(x)) return(character(0))
    head(x$gene_id[order(x$padj)], n)
}
jacc <- function(a, b) { u <- length(union(a, b)); if (!u) NA_real_ else length(intersect(a, b)) / u }

tabs <- lapply(paths, read_de)
names(tabs) <- labels
ref     <- tabs[[1]]
ref_sig <- sig_set(ref)
ref_top <- top_ids(ref)

rows <- lapply(seq_along(tabs), function(i) {
    d <- tabs[[i]]; s <- sig_set(d)
    if (i == 1) {
        sp <- 1; jc <- 1; ov <- length(ref_top)
    } else {
        m <- merge(ref[, c("gene_id", "log2FoldChange")],
                   d[, c("gene_id", "log2FoldChange")],
                   by = "gene_id", suffixes = c("_ref", "_x"))
        sp <- suppressWarnings(cor(m$log2FoldChange_ref, m$log2FoldChange_x,
                                   method = "spearman", use = "complete.obs"))
        jc <- jacc(ref_sig, s)
        ov <- length(intersect(ref_top, top_ids(d)))
    }
    data.frame(table = labels[i],
               n_genes_tested = sum(!is.na(d$padj)),
               n_sig = length(s),
               n_sig_up   = sum(!is.na(d$padj) & d$padj < PADJ & d$log2FoldChange >  LFC),
               n_sig_down = sum(!is.na(d$padj) & d$padj < PADJ & d$log2FoldChange < -LFC),
               spearman_lfc_vs_ref = round(sp, 4),
               jaccard_sig_vs_ref  = round(jc, 4),
               top20_overlap_vs_ref = ov,
               stringsAsFactors = FALSE)
})
res <- do.call(rbind, rows)
cat(sprintf("Reference = %s | sig: padj<%.2g & |log2FC|>%g\n\n", labels[1], PADJ, LFC))
print(res, row.names = FALSE)
write.csv(res, OUT, row.names = FALSE)
cat("\nWrote", OUT, "\n")
