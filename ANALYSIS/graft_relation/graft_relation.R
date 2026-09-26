# Does the xengsort graft/host fraction relate to the results?  Six tumours (human reads).
#   Rscript graft_relation.R <U> <outdir>
# Per sample: graft %, host %, both %, m = both / (graft + both) (the share of the human stream that is xengsort's
# conserved 'both' bin), DESeq2 size factor, rlog PC1/PC2, GSVA scores (Gaussian kernel on DESeq2 vst counts) of the six
# leading translation sets, the Neftel/Garofano signatures and the DSigDB ciclopirox set, and the CIBERSORT v1.04
# fractions.
# Tests (all with n = 6; graft % and group are partly collinear, so read the adjusted estimates as descriptive):
#   graft % by group (Welch); each score against graft % and against m (Pearson, Spearman); score ~ group + graft %
#   (group effect with and without graft %); DESeq2 genome-wide: ~ group; ~ graft_z + group; ~ graft_z alone; overlap
#   and fold-change agreement between the group and graft-% gene lists; the translation sets' mean log2FC under each.
suppressPackageStartupMessages({ library(DESeq2); library(GSVA); library(GSEABase) })
a <- commandArgs(trailingOnly = TRUE); U <- a[1]; out <- a[2]; dir.create(out, showWarnings = FALSE, recursive = TRUE)
A <- file.path(U, "ANALYSIS")
TUM <- c("IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B")
meta <- read.csv(file.path(A, "metadata_full.csv"))
meta <- meta[match(TUM, meta$sample), ]
meta$group <- factor(ifelse(meta$Classification == "Primary", "Primary", "Recurrent"), c("Primary", "Recurrent"))
meta$m <- meta$both_pct / (meta$graft_pct + meta$both_pct)
meta$graft_z <- as.numeric(scale(meta$graft_pct))
rownames(meta) <- meta$sample

cnt <- read.delim(file.path(A, "results_human_final/star_salmon/salmon.merged.gene_counts.tsv"), check.names = FALSE)
tpm <- read.delim(file.path(A, "results_human_final/star_salmon/salmon.merged.gene_tpm.tsv"), check.names = FALSE)
m <- round(as.matrix(cnt[, TUM])); rownames(m) <- cnt$gene_id
m <- m[rowSums(m >= 10) >= 1, ]
sym <- setNames(cnt$gene_name, cnt$gene_id)

# ---- DESeq2 models
fit <- function(design) {
  d <- DESeqDataSetFromMatrix(m, meta, design); DESeq(d, quiet = TRUE, minReplicatesForReplace = Inf)
}
d_g <- fit(~ group); r_g <- results(d_g, name = "group_Recurrent_vs_Primary")
d_gg <- fit(~ graft_z + group); r_gg <- results(d_gg, name = "group_Recurrent_vs_Primary"); r_gg_graft <- results(d_gg, name = "graft_z")
d_graft <- fit(~ graft_z); r_graft <- results(d_graft, name = "graft_z")
sf <- sizeFactors(d_g)
sig <- function(r, fc = 0) which(r$padj < 0.05 & abs(r$log2FoldChange) >= fc)
S_g <- rownames(r_g)[sig(r_g)]; S_gg <- rownames(r_gg)[sig(r_gg)]; S_graft <- rownames(r_graft)[sig(r_graft)]
jac <- function(x, y) if (length(union(x, y))) length(intersect(x, y)) / length(union(x, y)) else NA
de <- data.frame(
  model = c("~ group", "~ graft_z + group (group term)", "~ graft_z + group (graft term)", "~ graft_z (tumours only)"),
  padj05 = c(length(S_g), length(S_gg), sum(r_gg_graft$padj < 0.05, na.rm = TRUE), length(S_graft)),
  padj05_fc1.5 = c(length(sig(r_g, log2(1.5))), length(sig(r_gg, log2(1.5))), length(sig(r_gg_graft, log2(1.5))), length(sig(r_graft, log2(1.5)))))
agree <- data.frame(
  comparison = c("group genes vs graft-% genes", "group genes vs group|graft genes"),
  jaccard = c(jac(S_g, S_graft), jac(S_g, S_gg)),
  shared = c(length(intersect(S_g, S_graft)), length(intersect(S_g, S_gg))),
  spearman_lfc_all_genes = c(cor(r_g$log2FoldChange, r_graft$log2FoldChange, method = "spearman", use = "complete.obs"),
                             cor(r_g$log2FoldChange, r_gg$log2FoldChange, method = "spearman", use = "complete.obs")))
write.table(de, file.path(out, "de_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(agree, file.path(out, "de_agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
res_all <- data.frame(gene_id = rownames(r_g), symbol = sym[rownames(r_g)], baseMean = r_g$baseMean,
                      lfc_group = r_g$log2FoldChange, padj_group = r_g$padj,
                      lfc_group_adj = r_gg$log2FoldChange, padj_group_adj = r_gg$padj,
                      lfc_graft = r_graft$log2FoldChange, padj_graft = r_graft$padj)
write.table(res_all, file.path(out, "de_genes_three_models.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- gene sets
gmt <- function(p) { l <- strsplit(readLines(p), "\t"); setNames(lapply(l, function(x) unique(x[-(1:2)][x[-(1:2)] != ""])), sapply(l, `[`, 1)) }
hs <- gmt(file.path(A, "refs/pathways/human/combined_human.gmt"))
TR <- c("KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
        "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "KEGG_RIBOSOME",
        "REACTOME_SELENOAMINO_ACID_METABOLISM", "REACTOME_CELLULAR_RESPONSE_TO_STARVATION")
sets <- hs[intersect(TR, names(hs))]
subf <- c(file.path(U, "subtypes/subtype_signatures_remapped.gmt"), file.path(A, "holdout_IL68B/subtypes_src/subtype_signatures_remapped.gmt"))
sub <- gmt(subf[file.exists(subf)][1]); sets <- c(sets, sub)
dsig <- gmt(file.path(A, "refs/pathways/human/hs.dsigdb.v1.0.gmt"))
cic <- grep("^ciclopirox", names(dsig), ignore.case = TRUE, value = TRUE)
if (length(cic)) sets[["DSigDB_ciclopirox"]] <- unique(unlist(dsig[cic]))
# translation-set mean log2FC under each model (by symbol)
tr_lfc <- do.call(rbind, lapply(TR[TR %in% names(hs)], function(s) {
  k <- res_all$symbol %in% hs[[s]]
  data.frame(set = s, n = sum(k), mean_lfc_group = mean(res_all$lfc_group[k], na.rm = TRUE),
             mean_lfc_group_adj = mean(res_all$lfc_group_adj[k], na.rm = TRUE), mean_lfc_graft_per_sd = mean(res_all$lfc_graft[k], na.rm = TRUE))
}))
write.table(tr_lfc, file.path(out, "translation_sets_lfc_by_model.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# GSVA input: DESeq2 variance-stabilised counts (median-of-ratios), not TPM. The adversarial review (SUMMARY.md, C2/C4)
# found that un-normalised log2 TPM carries a per-sample composition constant (7SK/7SL/Y RNA and rRNA take 40-65 % of
# TPM) that correlates with graft %; count-based input removes it. Symbols summed after back-transforming.
vs <- assay(vst(d_g, blind = TRUE))
x <- 2^vs; rownames(x) <- sym[rownames(vs)]
x <- x[!is.na(rownames(x)) & rownames(x) != "", ]
x <- log2(rowsum(x, rownames(x)))
gs <- tryCatch(gsva(gsvaParam(x, sets, kcdf = "Gaussian"), verbose = FALSE),
               error = function(e) gsva(x, sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE))

# ---- per-sample table
rl <- assay(rlog(d_g, blind = TRUE)); top <- head(order(apply(rl, 1, var), decreasing = TRUE), 500)
pc <- prcomp(t(rl[top, ]))
ps <- data.frame(meta[, c("sample", "group", "graft_pct", "host_pct", "both_pct", "m")], size_factor = sf[TUM],
                 PC1 = pc$x[TUM, 1], PC2 = pc$x[TUM, 2], t(gs[, TUM]), check.names = FALSE)
v104 <- file.path(A, "cibersort/results/v104/fractions_v104_neftel4_confident.tsv")
if (file.exists(v104)) { f <- read.delim(v104, check.names = FALSE); rownames(f) <- f$sample
  ps$cs_MES <- f[TUM, "MES"]; ps$cs_AC <- f[TUM, "AC"]; ps$cs_NPC <- f[TUM, "NPC"]; ps$cs_fitR <- f[TUM, "Correlation"] }
write.table(ps, file.path(out, "per_sample.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- correlations and adjusted group effects
scores <- setdiff(names(ps), c("sample", "group", "graft_pct", "host_pct", "both_pct", "m"))
cr <- do.call(rbind, lapply(scores, function(s) {
  y <- ps[[s]]; g <- as.integer(ps$group == "Recurrent")
  f1 <- lm(y ~ g); f2 <- lm(y ~ g + ps$graft_pct)
  data.frame(score = s,
             r_graft = cor(y, ps$graft_pct), p_graft = cor.test(y, ps$graft_pct)$p.value,
             rho_graft = cor(y, ps$graft_pct, method = "spearman"),
             r_m = cor(y, ps$m), p_m = cor.test(y, ps$m)$p.value,
             group_effect = unname(coef(f1)[2]), group_p = summary(f1)$coefficients[2, 4],
             group_effect_adj_graft = unname(coef(f2)[2]), group_p_adj = summary(f2)$coefficients[2, 4],
             graft_coef_adj_group = unname(coef(f2)[3]), graft_p_adj = summary(f2)$coefficients[3, 4])
}))
write.table(cr, file.path(out, "score_correlations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
wt <- t.test(ps$graft_pct[ps$group == "Recurrent"], ps$graft_pct[ps$group == "Primary"])
cat("graft % primary", paste(ps$graft_pct[ps$group == "Primary"], collapse = ", "), "| recurrent",
    paste(ps$graft_pct[ps$group == "Recurrent"], collapse = ", "), "| Welch p", signif(wt$p.value, 3), "\n")
cat("point-biserial r(group, graft %):", round(cor(as.integer(ps$group == "Recurrent"), ps$graft_pct), 3), "\n")
print(de); print(agree); print(tr_lfc, digits = 3)
print(cr[order(-abs(cr$r_graft)), c("score", "r_graft", "p_graft", "r_m", "group_effect", "group_p", "group_effect_adj_graft", "group_p_adj")], digits = 3, row.names = FALSE)
