#!/usr/bin/env Rscript
# Write the DRIMSeq-filtered counts exactly as the pipeline's run_drimseq_filter.R does (write.table of counts(d),
# tab-separated, no quotes, no row names), from the dmDSdata.rds the pipeline published, keeping the named samples.
# The pipeline publishes only dmDSdata.rds (conf/modules.config pattern "*.rds"), not counts.tsv.
# Usage: drimseq_counts.R <dmDSdata.rds> <sample1,sample2,...> <out.tsv>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
suppressPackageStartupMessages(library(DRIMSeq))
d <- readRDS(args[1])
keep <- strsplit(args[2], ",")[[1]]
cts <- counts(d)
stopifnot(all(c("gene_id", "feature_id") %in% colnames(cts)), all(keep %in% colnames(cts)))
cts <- cts[, c("gene_id", "feature_id", keep), drop = FALSE]
write.table(cts, args[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("drimseq_counts: %d features in %d genes, %d samples (%s); DRIMSeq %s\n", nrow(cts),
            length(unique(cts$gene_id)), length(keep), paste(keep, collapse = ","), as.character(packageVersion("DRIMSeq"))))
