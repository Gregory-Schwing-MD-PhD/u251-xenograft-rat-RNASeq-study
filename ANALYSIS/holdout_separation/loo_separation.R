# Which held-out tumour lets primary and recurrent separate best?  All six tumours, then each one left out.
#   Rscript loo_separation.R <gene_counts.tsv> <out_prefix> [label]
# Per subset: the pipeline's filter (>= 10 reads in >= 1 sample), DESeq2 ~ Classification (Recurrent vs Primary);
# genes at padj < 0.05, and at padj < 0.05 with |fold| >= 1.5; rlog, PCA on the 500 most variable genes (DESeq2's
# plotPCA default): silhouette width of the two groups in PC1-2 and whether PC1 alone splits them; Ward.D2 clustering
# on 1 - Spearman over the 5000 most variable genes (the pipeline's exploratory settings), cut at k = 2.
# Choosing the held-out sample by this table is post hoc: all six rows are reported, not only the best.
suppressPackageStartupMessages({ library(DESeq2); library(cluster) })
a <- commandArgs(trailingOnly = TRUE)
counts <- read.delim(a[1], check.names = FALSE)
out <- a[2]; label <- if (length(a) > 2) a[3] else "human"
meta <- data.frame(sample = c("IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B"),
                   Classification = c(rep("Primary_U2", 3), rep("Recurrent_U2", 3)))
id <- if ("gene_id" %in% names(counts)) "gene_id" else names(counts)[1]
m <- round(as.matrix(counts[, meta$sample])); rownames(m) <- counts[[id]]
run <- function(keep) {
  md <- meta[meta$sample %in% keep, ]; md$Classification <- factor(md$Classification, c("Primary_U2", "Recurrent_U2"))
  x <- m[, md$sample]; x <- x[rowSums(x >= 10) >= 1, ]
  dds <- DESeqDataSetFromMatrix(x, md, ~ Classification)
  dds <- DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
  r <- results(dds, contrast = c("Classification", "Recurrent_U2", "Primary_U2"))
  n05 <- sum(r$padj < 0.05, na.rm = TRUE)
  n05fc <- sum(r$padj < 0.05 & abs(r$log2FoldChange) >= log2(1.5), na.rm = TRUE)
  rl <- assay(rlog(dds, blind = TRUE))
  top <- head(order(apply(rl, 1, var), decreasing = TRUE), 500)
  pc <- prcomp(t(rl[top, ]))
  xy <- pc$x[, 1:2]; g <- as.integer(md$Classification)
  sil <- mean(silhouette(g, dist(xy))[, "sil_width"])
  pc1_split <- (max(xy[g == 1, 1]) < min(xy[g == 2, 1])) || (min(xy[g == 1, 1]) > max(xy[g == 2, 1]))
  top5k <- head(order(apply(rl, 1, var), decreasing = TRUE), 5000)
  hc <- hclust(as.dist(1 - cor(rl[top5k, ], method = "spearman")), method = "ward.D2")
  k2 <- cutree(hc, 2)
  clus_ok <- length(unique(k2[g == 1])) == 1 && length(unique(k2[g == 2])) == 1 && k2[g == 1][1] != k2[g == 2][1]
  data.frame(reads = label, held_out = setdiff(meta$sample, keep)[1], n = length(keep), genes_tested = nrow(x),
             padj05 = n05, padj05_fc1.5 = n05fc, pc1_var = round(100 * pc$sdev[1]^2 / sum(pc$sdev^2), 1),
             pc2_var = round(100 * pc$sdev[2]^2 / sum(pc$sdev^2), 1), silhouette_pc12 = round(sil, 3),
             pc1_splits_groups = pc1_split, ward_k2_recovers_groups = clus_ok)
}
res <- run(meta$sample); res$held_out <- "none"
for (s in meta$sample) res <- rbind(res, run(setdiff(meta$sample, s)))
write.table(res, paste0(out, "_", label, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
