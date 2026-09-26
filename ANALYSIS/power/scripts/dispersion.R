lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
suppressMessages({library(edgeR); library(DESeq2); library(RNASeqPower)})
args <- commandArgs(TRUE); outdir <- args[1]
f <- "C:/Users/grego/OneDrive/Desktop/u251-xenograft-murine-RNASeq-study/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv"
x <- read.delim(f, check.names = FALSE)
PRI <- c("IL67B","IL68B","IL69B"); REC <- c("IL66B","NL70B","NL71B")
res <- list()
fit_disp <- function(samples, tag) {
  cnt <- round(as.matrix(x[, samples])); rownames(cnt) <- x$gene_id
  grp <- factor(ifelse(samples %in% PRI, "P", "R"), levels = c("P","R"))
  y <- DGEList(cnt, group = grp)
  keep <- filterByExpr(y); y <- y[keep, , keep.lib.sizes = FALSE]; y <- normLibSizes(y)
  design <- model.matrix(~grp)
  y <- estimateDisp(y, design, robust = TRUE)
  cpmmean <- rowMeans(cpm(y))
  normcnt <- sweep(y$counts, 2, y$samples$lib.size * y$samples$norm.factors, "/") * mean(y$samples$lib.size * y$samples$norm.factors)
  mu <- rowMeans(normcnt)
  qt <- quantile(sqrt(y$tagwise.dispersion), c(.1,.25,.5,.75,.9))
  fitq <- glmQLFit(y, design, robust = TRUE)
  # DESeq2 for comparison
  dds <- DESeqDataSetFromMatrix(cnt[keep, ], data.frame(grp = grp), ~grp)
  dds <- estimateSizeFactors(dds); dds <- estimateDispersions(dds, quiet = TRUE)
  dd <- dispersions(dds); dg <- mcols(dds)$dispGeneEst
  out <- list(tag = tag, n_genes = sum(keep), lib_sizes_M = round(y$samples$lib.size/1e6, 2),
              common_disp = y$common.dispersion, common_BCV = sqrt(y$common.dispersion),
              tagwise_BCV_quantiles = as.list(qt),
              trended_BCV_median = median(sqrt(y$trended.dispersion)),
              QL_df_prior_median = median(fitq$df.prior),
              deseq2_final_disp_median = median(dd, na.rm = TRUE), deseq2_final_BCV_median = sqrt(median(dd, na.rm = TRUE)),
              deseq2_genewise_disp_median = median(dg, na.rm = TRUE),
              mean_norm_count_median = median(mu), mean_norm_count_quantiles = as.list(quantile(mu, c(.1,.25,.5,.75,.9))))
  saveRDS(list(mu = mu, disp = y$tagwise.dispersion, trended = y$trended.dispersion), file.path(outdir, paste0("disp_", tag, ".rds")))
  print(str(out)); out
}
res$all6 <- fit_disp(c(PRI, REC), "all6")
res$no_IL68B <- fit_disp(setdiff(c(PRI, REC), "IL68B"), "no_IL68B")
res$no_IL66B <- fit_disp(setdiff(c(PRI, REC), "IL66B"), "no_IL66B")
# within-group dispersion of the three primaries only and three recurrences only (null model, no group)
for (g in list(list("pri_only", PRI), list("rec_only", REC))) {
  cnt <- round(as.matrix(x[, g[[2]]])); y <- DGEList(cnt); keep <- filterByExpr(y, min.count = 10, group = rep(1, 3)); y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y); y <- estimateDisp(y, matrix(1, 3, 1))
  res[[g[[1]]]] <- list(common_BCV = sqrt(y$common.dispersion), tagwise_BCV_median = median(sqrt(y$tagwise.dispersion)), n_genes = sum(keep))
  print(g[[1]]); print(res[[g[[1]]]])
}
jsonlite::write_json(res, file.path(outdir, "dispersion.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
