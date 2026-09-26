#!/usr/bin/env Rscript
# Step B of run_host_verify.sbatch: preranked gene-set tests (fgsea) on every host fit's DESeq2 Wald statistic from
# step A (host_verify_de.R, ranks/<fit>.tsv, manifest fits.tsv). 10,000 permutations, sets of 15-500 genes; hallmark,
# KEGG, GO BP and the brain collection tested separately (each with its own BH), and the combined collection for the
# full fit and the NL batch fit, to set beside the pipeline's Broad GSEA of the same contrast.
#   Rscript host_verify_gsea.R <outdir> <gmt dir> <R library root> [cores] [dir with the Broad GSEA reports]
# Genes are ranked by symbol; a symbol carried by several Ensembl ids keeps the id with the largest |stat|.
# The ten 3-v-3 relabellings give an empirical null for how many sets reach padj < 0.05 when the labels mean nothing.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "hv_json.R"))
a <- commandArgs(trailingOnly = TRUE)
out <- a[1]; gmt_dir <- a[2]
lib <- file.path(a[3], paste0("R", R.version$major, ".", sub("\\..*$", "", R.version$minor)))
cores <- if (length(a) >= 4) as.integer(a[4]) else 4L
broad <- if (length(a) >= 5) a[5] else ""
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
installed_now <- FALSE
if (!requireNamespace("fgsea", quietly = TRUE)) {
  message("fgsea is not in this image; installing it into ", lib)
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
  BiocManager::install("fgsea", lib = lib, update = FALSE, ask = FALSE)
  installed_now <- TRUE
}
suppressPackageStartupMessages({ library(fgsea); library(BiocParallel) })
BP <- if (cores > 1) MulticoreParam(cores) else SerialParam()
NPERM <- 10000L
dir.create(file.path(out, "gsea"), showWarnings = FALSE)
wt <- function(df, name) write.table(df, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)

gmts <- c(hallmark = "rn.hallmark.v2023.2.gmt", kegg = "rn.kegg.v2023.2.gmt", gobp = "rn.gobp.v2023.2.gmt",
          brain = "rn.brain_gmt_v2.gmt", combined = "combined_rat.gmt")
core <- c("hallmark", "kegg", "gobp", "brain")
PW <- lapply(gmts, function(f) {
  p <- gmtPathways(file.path(gmt_dir, f)); p <- p[!duplicated(names(p))]
  lapply(p, function(g) unique(g[!is.na(g) & nzchar(g)])) })

man <- read.delim(file.path(out, "fits.tsv"), stringsAsFactors = FALSE)
man$true_split <- as.logical(man$true_split); man$nl_together <- as.logical(man$nl_together)
true_lab <- man$label[man$type == "relabel_3v3" & man$true_split]

load_stats <- function(label) {
  r <- read.delim(file.path(out, "ranks", paste0(label, ".tsv")), stringsAsFactors = FALSE)
  r <- r[!is.na(r$stat) & !is.na(r$gene_name) & nzchar(r$gene_name), ]
  r <- r[order(-abs(r$stat)), ]; r <- r[!duplicated(r$gene_name), ]
  setNames(r$stat, r$gene_name)
}
run1 <- function(label, coll) {
  f <- file.path(out, "gsea", paste0(label, ".", coll, ".tsv.gz"))
  if (file.exists(f)) return(read.delim(gzfile(f), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  set.seed(1234)
  res <- as.data.frame(fgseaSimple(PW[[coll]], load_stats(label), nperm = NPERM, minSize = 15, maxSize = 500,
                                   BPPARAM = BP))
  res$leadingEdge <- vapply(res$leadingEdge, paste, "", collapse = ",")
  res <- res[order(res$pval, -abs(res$NES)), c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
  gz <- gzfile(f, "w"); write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  cat(sprintf("gsea %-28s %-9s %5d sets %6.1f s\n", label, coll, nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res
}

# the true relabelling is the full fit (step A checks the statistics are identical), so it reuses the full fit's sets
R <- list()
for (i in seq_len(nrow(man))) {
  lab <- man$label[i]
  if (lab == true_lab) next
  colls <- if (man$type[i] %in% c("full", "batch_2v4")) names(gmts) else core
  if (man$type[i] == "single_1v3") colls <- core
  for (co in colls) R[[paste(lab, co, sep = "|")]] <- run1(lab, co)
}
for (co in core) R[[paste(true_lab, co, sep = "|")]] <- R[[paste("full", co, sep = "|")]]
get <- function(lab, co) R[[paste(lab, co, sep = "|")]]

## counts per fit and collection -------------------------------------------------------------------------------------
cn <- do.call(rbind, lapply(names(R), function(k) {
  p <- strsplit(k, "|", fixed = TRUE)[[1]]; d <- R[[k]]; m <- man[man$label == p[1], ]
  data.frame(fit = p[1], type = m$type, true_split = m$true_split, nl_together = m$nl_together, collection = p[2],
             sets_tested = nrow(d), padj05 = sum(d$padj < 0.05, na.rm = TRUE),
             padj05_up = sum(d$padj < 0.05 & d$NES > 0, na.rm = TRUE),
             padj05_down = sum(d$padj < 0.05 & d$NES < 0, na.rm = TRUE),
             padj25 = sum(d$padj < 0.25, na.rm = TRUE), min_padj = suppressWarnings(min(d$padj, na.rm = TRUE)),
             stringsAsFactors = FALSE) }))
wt(cn, "r4_gsea_counts.tsv")

## empirical null from the ten relabellings ----------------------------------------------------------------------------
nul <- do.call(rbind, lapply(core, function(co) {
  s <- cn[cn$type == "relabel_3v3" & cn$collection == co, ]
  t5 <- s$padj05[s$true_split]; t25 <- s$padj25[s$true_split]
  data.frame(collection = co, true_padj05 = t5, true_rank_padj05 = min(rank(-s$padj05, ties.method = "min")[s$true_split]),
             splits_with_padj05_ge_true = sum(s$padj05 >= t5), max_split_padj05 = max(s$padj05),
             max_split = paste(s$fit[s$padj05 == max(s$padj05)], collapse = ";"),
             median_other_padj05 = median(s$padj05[!s$true_split]),
             true_padj25 = t25, true_rank_padj25 = min(rank(-s$padj25, ties.method = "min")[s$true_split]),
             median_other_padj25 = median(s$padj25[!s$true_split]), stringsAsFactors = FALSE) }))
wt(nul, "r4_gsea_relabel_null.tsv")

## the full fit's top sets, and how they fare without each tumour and under the relabellings ----------------------------
loo_labs <- man$label[man$type == "leave_one_out"]
null_labs <- man$label[man$type == "relabel_3v3" & !man$true_split]
top <- do.call(rbind, lapply(c(core, "combined"), function(co) {
  f <- get("full", co)
  up <- head(f[f$NES > 0, ], 10); dn <- head(f[f$NES < 0, ], 10)
  t <- rbind(cbind(direction = rep("up", nrow(up)), up), cbind(direction = rep("down", nrow(dn)), dn))
  if (!nrow(t)) return(NULL)
  lk <- function(lab, pw) { d <- get(lab, co); if (is.null(d)) return(matrix(NA_real_, length(pw), 2))
                            m <- d[match(pw, d$pathway), ]; cbind(m$NES, m$padj) }
  below <- function(x, th) !is.na(x) & x < th
  for (j in c("n_loo_same_sign_padj25", "n_loo_same_sign_padj05", "n_null_splits_padj05")) t[[j]] <- 0L
  for (l in loo_labs) t[[paste0("NES_wo_", sub("^loo_", "", l))]] <- NA_real_
  if (co %in% core) {
    for (l in loo_labs) { x <- lk(l, t$pathway); ss <- !is.na(x[, 1]) & sign(x[, 1]) == sign(t$NES)
      t$n_loo_same_sign_padj25 <- t$n_loo_same_sign_padj25 + (ss & below(x[, 2], 0.25))
      t$n_loo_same_sign_padj05 <- t$n_loo_same_sign_padj05 + (ss & below(x[, 2], 0.05))
      t[[paste0("NES_wo_", sub("^loo_", "", l))]] <- round(x[, 1], 3) }
    for (l in null_labs) { x <- lk(l, t$pathway); t$n_null_splits_padj05 <- t$n_null_splits_padj05 + below(x[, 2], 0.05) }
  } else { t$n_loo_same_sign_padj25 <- NA; t$n_loo_same_sign_padj05 <- NA; t$n_null_splits_padj05 <- NA }
  xb <- lk("batch_NL_vs_IL", t$pathway); t$batch_NL_vs_IL_NES <- round(xb[, 1], 3); t$batch_NL_vs_IL_padj <- xb[, 2]
  t$leadingEdge <- substr(t$leadingEdge, 1, 300)
  cbind(collection = co, t, stringsAsFactors = FALSE) }))
top_cols <- c("collection", "direction", "pathway", "size", "NES", "pval", "padj", "n_loo_same_sign_padj25",
              "n_loo_same_sign_padj05", "n_null_splits_padj05", "batch_NL_vs_IL_NES", "batch_NL_vs_IL_padj")
wt(top, "r4_gsea_top_full.tsv")
for (co in c(core, "combined")) wt(get("full", co), paste0("r4_gsea_full_", co, ".tsv"))
# the five strongest sets each way in every fit (what the relabellings and the NL contrast pick up)
tbf <- do.call(rbind, lapply(names(R), function(k) {
  p <- strsplit(k, "|", fixed = TRUE)[[1]]; d <- R[[k]]
  s <- rbind(head(d[d$NES > 0, ], 5), head(d[d$NES < 0, ], 5))
  if (!nrow(s)) return(NULL)
  data.frame(fit = p[1], collection = p[2], s[, c("pathway", "size", "NES", "pval", "padj")], stringsAsFactors = FALSE) }))
wt(tbf, "r4_gsea_top_by_fit.tsv")

## beside the pipeline's Broad GSEA (Diff_of_Classes on normalised counts, gene-set permutation, combined collection) -----
cmp <- NULL
rf <- file.path(broad, c("host_therapy.combined_rat.gsea_report_for_Recurrent_U2.tsv",
                         "host_therapy.combined_rat.gsea_report_for_Primary_U2.tsv"))
if (nzchar(broad) && all(file.exists(rf))) {
  b <- do.call(rbind, lapply(rf, function(x) { d <- read.delim(x, check.names = FALSE, stringsAsFactors = FALSE)
    data.frame(pathway = d$NAME, NES = as.numeric(d$NES), p = as.numeric(d[["NOM p-val"]]), fdr = as.numeric(d[["FDR q-val"]])) }))
  f <- get("full", "combined"); f$key <- toupper(f$pathway)
  m <- merge(f, b, by.x = "key", by.y = "pathway", suffixes = c("_fgsea", "_broad"))
  topk <- function(d, nes, pv, k = 20, s = 1) { d <- d[sign(d[[nes]]) == s, ]; head(d$key[order(d[[pv]], -abs(d[[nes]]))], k) }
  cmp <- list(sets_matched = nrow(m), sets_fgsea = nrow(f), sets_broad = nrow(b),
              spearman_NES = cor(m$NES_fgsea, m$NES_broad, method = "spearman"),
              top20_up_overlap = length(intersect(topk(m, "NES_fgsea", "pval"), topk(m, "NES_broad", "p"))),
              top20_down_overlap = length(intersect(topk(m, "NES_fgsea", "pval", s = -1), topk(m, "NES_broad", "p", s = -1))),
              broad_fdr25_up = sum(b$fdr < 0.25 & b$NES > 0, na.rm = TRUE),
              broad_fdr25_down = sum(b$fdr < 0.25 & b$NES < 0, na.rm = TRUE),
              broad_min_fdr = min(b$fdr, na.rm = TRUE), broad_min_nominal_p = min(b$p, na.rm = TRUE),
              broad_top10_up = head(b[b$NES > 0, ][order(b$p[b$NES > 0], -b$NES[b$NES > 0]), c("pathway", "NES", "p", "fdr")], 10),
              broad_top10_down = head(b[b$NES < 0, ][order(b$p[b$NES < 0], b$NES[b$NES < 0]), c("pathway", "NES", "p", "fdr")], 10))
  wt(m[, c("key", "size", "NES_fgsea", "pval", "padj", "NES_broad", "p", "fdr")], "r4_fgsea_vs_broad_combined.tsv")
}

gs <- list(method = list(fgsea = as.character(packageVersion("fgsea")), installed_this_run = installed_now,
                         library = lib, R = R.version.string, function_used = "fgseaSimple", nperm = NPERM,
                         set_size = "15-500", statistic = "DESeq2 Wald stat (unshrunk)", seed = 1234,
                         collections_tested_separately = core),
           counts = cn[cn$type %in% c("full", "leave_one_out", "batch_2v4", "single_1v3"), ],
           relabel_null = nul, top_full = top[, top_cols], vs_broad = cmp)
gs_txt <- to_json(gs, 1)
writeLines(gs_txt, file.path(out, "summary_gsea.json"))
de_f <- file.path(out, "summary_de.json")
de_txt <- if (file.exists(de_f)) paste(readLines(de_f), collapse = "\n") else "null"
writeLines(paste0("{\n\"de\": ", de_txt, ",\n\"gsea\": ", gs_txt, "\n}"), file.path(out, "summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_gsea.txt"))
cat("step B done\n")
