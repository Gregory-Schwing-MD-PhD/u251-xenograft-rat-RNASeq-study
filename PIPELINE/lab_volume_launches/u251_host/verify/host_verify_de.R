#!/usr/bin/env Rscript
# Host (rat) differential expression, Recurrent_U2 vs Primary_U2 (deck backup B7): verification of the
# nf-core/differentialabundance table. Step A of run_host_verify.sbatch. Runs in the pipeline's own DESeq2 container
# (bioconductor-deseq2 1.34.0, ashr 2.2-54, R 4.1.3) so that R0 can reproduce the table exactly.
#
#   Rscript host_verify_de.R <outdir> <rat gene counts> <rat gene lengths> <pipeline host_therapy table> \
#           <rat_human_1to1.tsv> <human gene counts> [cores] [ashr|no_ashr]
#
# Settings are those of the pipeline's DESEQ2_DIFFERENTIAL task (u251_host/work_de/00/9decfd.../.command.sh):
#   filter   >= 10 reads in >= 1 of all NINE samples (CUSTOM_MATRIXFILTER runs before the contrast subset)
#   subset   the six tumours (differential_subset_to_contrast_samples)
#   design   ~ 0 + Classification (cell means); salmon gene lengths as avgTxLength (tximport-style offsets)
#   DESeq    Wald, parametric, sfType ratio, minReplicatesForReplace 99 (no outlier replacement), useT FALSE
#   results  alpha 0.1, independent filtering, BH, minmu 0.5, lfcThreshold 0
#   lfc      lfcShrink(type = "ashr"); pvalue and padj stay those of the unshrunk Wald test
# Every other fit here (leave-one-out, relabellings, NL batch, IL66B alone) uses the same universe and settings.
#
# R0 reproduce the table; R1 leave one tumour out; R2 all 10 balanced 3-v-3 relabellings; R3 NL-prefix contrast and
# IL66B alone against the primaries; R5 human reads beside rat reads for the DE genes. Each fit's Wald statistic goes
# to ranks/ for step B (host_verify_gsea.R).
suppressPackageStartupMessages({ library(DESeq2); library(BiocParallel) })
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "hv_json.R"))

a <- commandArgs(trailingOnly = TRUE)
stopifnot(length(a) >= 6)
out <- a[1]; f_counts <- a[2]; f_len <- a[3]; f_pipe <- a[4]; f_orth <- a[5]; f_hum <- a[6]
cores <- if (length(a) >= 7) as.integer(a[7]) else 4L
use_ashr <- !(length(a) >= 8 && a[8] == "no_ashr")
if (use_ashr && !requireNamespace("ashr", quietly = TRUE)) stop("ashr is not installed in this container")
for (d in c("fits", "ranks")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
BP <- if (cores > 1) MulticoreParam(cores) else SerialParam()
set.seed(1234)
wt <- function(df, name) write.table(df, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE)

tum   <- c("IL67B", "IL68B", "IL69B", "IL66B", "NL70B", "NL71B")    # the pipeline's sample-sheet order
prim  <- tum[1:3]; recur <- tum[4:6]
ctrl  <- c("IL64B", "N168B", "N269B")

cnt <- read.delim(f_counts, check.names = FALSE, stringsAsFactors = FALSE)
rownames(cnt) <- cnt$gene_id
sym <- setNames(ifelse(is.na(cnt$gene_name) | cnt$gene_name == "", cnt$gene_id, cnt$gene_name), cnt$gene_id)
universe <- rownames(cnt)[rowSums(as.matrix(cnt[, c(tum, ctrl)]) >= 10) >= 1]
len <- read.delim(f_len, check.names = FALSE, stringsAsFactors = FALSE); rownames(len) <- len$gene_id
stopifnot(all(universe %in% rownames(len)))
pipe <- read.delim(f_pipe, stringsAsFactors = FALSE); rownames(pipe) <- pipe$gene_id
cat("universe", length(universe), "genes; pipeline table", nrow(pipe), "rows\n")

## one DESeq2 fit: `ref` samples are the reference level, `tgt` the target; lfc = target over reference ---------------
fit <- function(label, ref, tgt, design = ~ 0 + Classification, keep_dds = FALSE) {
  f <- file.path(out, "fits", paste0(label, ".tsv.gz"))
  if (file.exists(f) && !keep_dds) return(read.delim(gzfile(f), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  smp <- tum[tum %in% c(ref, tgt)]
  cls <- factor(ifelse(smp %in% tgt, "Recurrent_U2", "Primary_U2"), levels = c("Primary_U2", "Recurrent_U2"))
  md <- data.frame(sample = smp, Classification = cls, row.names = smp)
  dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt[universe, smp])), md, design)
  L <- as.matrix(len[universe, smp]); dimnames(L) <- dimnames(dds)
  assays(dds)[["avgTxLength"]] <- L
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = 99, useT = FALSE,
               sfType = "ratio", parallel = cores > 1, BPPARAM = BP, quiet = TRUE)
  ctr <- c("Classification", "Recurrent_U2", "Primary_U2")
  r <- results(dds, lfcThreshold = 0, altHypothesis = "greaterAbs", independentFiltering = TRUE, alpha = 0.1,
               pAdjustMethod = "BH", minmu = 0.5, contrast = ctr)
  df <- data.frame(gene_id = rownames(r), gene_name = unname(sym[rownames(r)]), baseMean = r$baseMean,
                   lfc_mle = r$log2FoldChange, lfcSE_mle = r$lfcSE, stat = r$stat, pvalue = r$pvalue,
                   padj = r$padj, stringsAsFactors = FALSE)
  if (use_ashr) {
    s <- lfcShrink(dds, type = "ashr", contrast = ctr, quiet = TRUE)
    stopifnot(identical(rownames(s), df$gene_id))
    df$lfc_ashr <- s$log2FoldChange; df$lfcSE_ashr <- s$lfcSE; df$padj_shrinkcall <- s$padj
  } else {
    df$lfc_ashr <- NA_real_; df$lfcSE_ashr <- NA_real_; df$padj_shrinkcall <- NA_real_
  }
  gz <- gzfile(f, "w"); write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  rk <- df[!is.na(df$stat), c("gene_id", "gene_name", "stat")]
  write.table(rk, file.path(out, "ranks", paste0(label, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("fit %-28s ref %-26s tgt %-26s %5.1f s\n", label, paste(ref, collapse = ","),
              paste(tgt, collapse = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (keep_dds) attr(df, "dds") <- dds
  df
}

sig <- function(df) !is.na(df$padj) & df$padj < 0.05
cnts <- function(df) {
  s <- sig(df)
  c(genes_with_padj = sum(!is.na(df$padj)), padj05 = sum(s), padj05_up = sum(s & df$lfc_mle > 0),
    padj05_down = sum(s & df$lfc_mle < 0),
    padj05_absLFC1_ashr = if (use_ashr) sum(s & abs(df$lfc_ashr) >= 1) else NA,
    padj05_absLFC1_mle = sum(s & abs(df$lfc_mle) >= 1),
    padj05_fc1.5_mle = sum(s & abs(df$lfc_mle) >= log2(1.5)))
}

## the fits ----------------------------------------------------------------------------------------------------------
manifest <- data.frame(label = character(), type = character(), reference = character(), target = character(),
                       true_split = logical(), nl_together = logical(), stringsAsFactors = FALSE)
addm <- function(label, type, ref, tgt, true_split = FALSE) {
  nl_tog <- (all(c("NL70B", "NL71B") %in% ref) || all(c("NL70B", "NL71B") %in% tgt))
  manifest[nrow(manifest) + 1, ] <<- list(label, type, paste(ref, collapse = ","), paste(tgt, collapse = ","),
                                          true_split, nl_tog)
}
full <- fit("full", prim, recur, keep_dds = TRUE); dds_full <- attr(full, "dds"); attr(full, "dds") <- NULL
addm("full", "full", prim, recur, TRUE)
full_int <- fit("full_design_intercept", prim, recur, design = ~ Classification)
loo <- list()
for (s in tum) { lab <- paste0("loo_", s); loo[[s]] <- fit(lab, setdiff(prim, s), setdiff(recur, s))
                 addm(lab, "leave_one_out", setdiff(prim, s), setdiff(recur, s)) }
# every partition of the six into two groups of three; the group holding IL67B is the reference
splits <- list()
for (p in combn(setdiff(tum, "IL67B"), 2, simplify = FALSE)) {
  ref <- tum[tum %in% c("IL67B", p)]; tgt <- setdiff(tum, ref)
  lab <- paste0("split_", paste(tgt, collapse = "_"))
  splits[[lab]] <- fit(lab, ref, tgt)
  addm(lab, "relabel_3v3", ref, tgt, true_split = identical(sort(tgt), sort(recur)))
}
batch <- fit("batch_NL_vs_IL", c(prim, "IL66B"), c("NL70B", "NL71B")); addm("batch_NL_vs_IL", "batch_2v4", c(prim, "IL66B"), c("NL70B", "NL71B"))
il66 <- fit("IL66B_vs_primaries", prim, "IL66B"); addm("IL66B_vs_primaries", "single_1v3", prim, "IL66B")
wt(manifest, "fits.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_de.txt"))

## R0: reproduce the pipeline's table ---------------------------------------------------------------------------------
p_sig <- pipe$gene_id[!is.na(pipe$padj) & pipe$padj < 0.05]
m_sig <- full$gene_id[sig(full)]
cm <- intersect(pipe$gene_id, full$gene_id)
P <- pipe[cm, ]; M <- full[match(cm, full$gene_id), ]
dif <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  if (!any(ok)) return(c(max_abs = NA, max_rel = NA))
  big <- ok & abs(y) > 1e-12
  c(max_abs = max(abs(x[ok] - y[ok])), max_rel = if (any(big)) max(abs(x[big] - y[big]) / abs(y[big])) else NA)
}
implied_p <- 2 * pnorm(-abs(P$log2FoldChange / P$lfcSE))
imp_ok <- !is.na(implied_p) & !is.na(P$pvalue) & P$pvalue > 0
r0 <- list(
  universe_n_here = length(universe), universe_n_pipeline = nrow(pipe),
  universe_identical = setequal(universe, pipe$gene_id),
  padj_nonNA_here = sum(!is.na(full$padj)), padj_nonNA_pipeline = sum(!is.na(pipe$padj)),
  padj_NA_pattern_mismatches = sum(is.na(P$padj) != is.na(M$padj)),
  padj05_here = length(m_sig), padj05_pipeline = length(p_sig),
  padj05_intersection = length(intersect(m_sig, p_sig)),
  jaccard = if (length(union(m_sig, p_sig))) length(intersect(m_sig, p_sig)) / length(union(m_sig, p_sig)) else NA,
  only_pipeline = unname(sym[setdiff(p_sig, m_sig)]), only_here = unname(sym[setdiff(m_sig, p_sig)]),
  pipeline_counts = c(padj05 = length(p_sig),
                      up = sum(pipe[p_sig, "log2FoldChange"] > 0), down = sum(pipe[p_sig, "log2FoldChange"] < 0),
                      padj05_absLFC1 = sum(abs(pipe[p_sig, "log2FoldChange"]) >= 1)),
  here_counts = cnts(full),
  baseMean_diff = dif(M$baseMean, P$baseMean), pvalue_diff = dif(M$pvalue, P$pvalue), padj_diff = dif(M$padj, P$padj),
  table_lfc_vs_ashr = dif(M$lfc_ashr, P$log2FoldChange), table_lfc_vs_mle = dif(M$lfc_mle, P$log2FoldChange),
  table_lfcSE_vs_ashr = dif(M$lfcSE_ashr, P$lfcSE), table_lfcSE_vs_mle = dif(M$lfcSE_mle, P$lfcSE),
  table_pvalue_implied_by_its_lfc_over_lfcSE_within_1pct = mean(abs(implied_p[imp_ok] / P$pvalue[imp_ok] - 1) < 0.01),
  shrinkcall_padj_equals_results_padj = if (use_ashr) isTRUE(all.equal(full$padj, full$padj_shrinkcall)) else NA,
  design_intercept_vs_cellmeans = c(stat_max_abs_diff = max(abs(full_int$stat - full$stat), na.rm = TRUE),
                                    padj05_same_set = setequal(full_int$gene_id[sig(full_int)], m_sig))
)
r0$lfc_column_is <- if (!use_ashr) "unknown (ashr not run)" else
  if (isTRUE(r0$table_lfc_vs_ashr[["max_abs"]] < 1e-4)) "ashr-shrunk (lfcShrink type ashr)" else
  if (isTRUE(r0$table_lfc_vs_mle[["max_abs"]] < 1e-4)) "unshrunk MLE" else "neither ashr nor MLE matches"
r0$reproduces <- isTRUE(r0$universe_identical) && isTRUE(r0$jaccard == 1) &&
  isTRUE(r0$padj_diff[["max_rel"]] < 1e-4) && r0$padj_NA_pattern_mismatches == 0
wt(data.frame(key = names(unlist(r0)), value = unname(sapply(unlist(r0), as.character))), "r0_reproduce.tsv")
cat("R0 reproduces:", r0$reproduces, "| lfc column:", r0$lfc_column_is, "\n")

## R1: leave one tumour out -----------------------------------------------------------------------------------------
g38 <- p_sig
dir38 <- sign(pipe[g38, "log2FoldChange"])
r1_counts <- data.frame(fit = c("full", paste0("loo_", tum)), held_out = c("none", tum),
                        do.call(rbind, c(list(cnts(full)), lapply(loo, cnts))), check.names = FALSE)
r1_counts$of_pipeline_genes_kept_same_sign <- sapply(c(list(full), loo), function(df) {
  d <- df[match(g38, df$gene_id), ]; sum(!is.na(d$padj) & d$padj < 0.05 & sign(d$lfc_mle) == dir38) })
wt(r1_counts, "r1_loo_counts.tsv")
gt <- data.frame(gene_id = g38, gene_name = unname(sym[g38]), direction = ifelse(dir38 > 0, "up", "down"),
                 pipeline_lfc = pipe[g38, "log2FoldChange"], pipeline_padj = pipe[g38, "padj"], stringsAsFactors = FALSE)
keep_n <- integer(length(g38)); lost_in <- character(length(g38))
for (s in tum) {
  d <- loo[[s]][match(g38, loo[[s]]$gene_id), ]
  gt[[paste0("lfc_ashr_wo_", s)]] <- d$lfc_ashr; gt[[paste0("lfc_mle_wo_", s)]] <- d$lfc_mle
  gt[[paste0("padj_wo_", s)]] <- d$padj
  k <- !is.na(d$padj) & d$padj < 0.05 & sign(d$lfc_mle) == dir38
  keep_n <- keep_n + k; lost_in[!k] <- paste0(lost_in[!k], ifelse(nzchar(lost_in[!k]), ",", ""), s)
}
gt$n_loo_kept <- keep_n; gt$lost_without <- lost_in
gt <- gt[order(gt$direction, -gt$n_loo_kept, gt$pipeline_padj), ]
wt(gt, "r1_loo_genes.tsv")
mmp <- names(sym)[sym == "Mmp13"][1]
r1_mmp13 <- do.call(rbind, lapply(c(full = list(full), setNames(loo, paste0("loo_", tum))), function(df) {
  d <- df[df$gene_id == mmp, ]
  data.frame(lfc_ashr = d$lfc_ashr, fold_down_ashr = 2^-d$lfc_ashr, lfc_mle = d$lfc_mle, fold_down_mle = 2^-d$lfc_mle,
             padj = d$padj) }))
r1_mmp13 <- data.frame(fit = rownames(r1_mmp13), r1_mmp13); wt(r1_mmp13, "r1_mmp13.tsv")

## R2: relabelling null -------------------------------------------------------------------------------------------------
sm <- manifest[manifest$type == "relabel_3v3", ]
r2 <- data.frame(sm[, c("label", "reference", "target", "true_split", "nl_together")],
                 do.call(rbind, lapply(splits[sm$label], cnts)), check.names = FALSE)
rank_of <- function(col) { r <- rank(-r2[[col]], ties.method = "min"); c(rank = r[r2$true_split], tied = sum(r2[[col]] == r2[[col]][r2$true_split])) }
r2_sum <- list(
  n_splits = nrow(r2),
  true_split_padj05 = r2$padj05[r2$true_split], true_split_rank_padj05 = rank_of("padj05"),
  true_split_padj05_absLFC1 = r2$padj05_absLFC1_ashr[r2$true_split],
  true_split_rank_padj05_absLFC1 = if (use_ashr) rank_of("padj05_absLFC1_ashr") else NA,
  largest_split_padj05 = r2$label[r2$padj05 == max(r2$padj05)], largest_padj05 = max(r2$padj05),
  median_other_splits_padj05 = median(r2$padj05[!r2$true_split]),
  true_split_equals_full_fit = isTRUE(all.equal(splits[[sm$label[sm$true_split]]]$stat, full$stat)),
  nl_together_splits = r2$label[r2$nl_together],
  nl_together_padj05 = r2$padj05[r2$nl_together], nl_apart_padj05 = r2$padj05[!r2$nl_together])
r2 <- r2[order(-r2$padj05), ]; wt(r2, "r2_relabel_splits.tsv")
# what drives each relabelling: its strongest genes each way, and how far its statistic agrees with the true split's
r2_top <- do.call(rbind, lapply(sm$label, function(l) {
  d <- splits[[l]]; d <- d[!is.na(d$stat), ]
  data.frame(label = l, padj05 = sum(sig(d)),
             spearman_stat_vs_true = cor(d$stat, full$stat[match(d$gene_id, full$gene_id)], method = "spearman",
                                         use = "complete.obs"),
             top_up_in_target = paste(head(d$gene_name[order(-d$stat)], 15), collapse = ","),
             top_down_in_target = paste(head(d$gene_name[order(d$stat)], 15), collapse = ","), stringsAsFactors = FALSE) }))
wt(r2_top[order(-r2_top$padj05), ], "r2_split_top_genes.tsv")

## R3: NL-prefix contrast; IL66B alone ----------------------------------------------------------------------------------
nc <- counts(dds_full, normalized = TRUE)
pm <- rowMeans(nc[g38, prim, drop = FALSE])
l2 <- function(x) log2((x + 1) / (pm + 1))
bt <- batch[match(g38, batch$gene_id), ]; i1 <- il66[match(g38, il66$gene_id), ]
r3g <- data.frame(gene_id = g38, gene_name = unname(sym[g38]), direction = ifelse(dir38 > 0, "up", "down"),
                  pipeline_lfc = pipe[g38, "log2FoldChange"], round(nc[g38, tum], 1),
                  l2_IL66B_vs_prim = l2(nc[g38, "IL66B"]), l2_NL70B_vs_prim = l2(nc[g38, "NL70B"]),
                  l2_NL71B_vs_prim = l2(nc[g38, "NL71B"]), check.names = FALSE, stringsAsFactors = FALSE)
pmin_ <- apply(nc[g38, prim, drop = FALSE], 1, min); pmax_ <- apply(nc[g38, prim, drop = FALSE], 1, max)
beyond <- function(x) ifelse(dir38 > 0, x > pmax_, x < pmin_)     # outside every primary, in the recurrence direction
r3g$IL66B_beyond_all_primaries <- beyond(nc[g38, "IL66B"])
r3g$NL70B_beyond_all_primaries <- beyond(nc[g38, "NL70B"])
r3g$NL71B_beyond_all_primaries <- beyond(nc[g38, "NL71B"])
r3g$IL66B_inside_primary_range <- nc[g38, "IL66B"] >= pmin_ & nc[g38, "IL66B"] <= pmax_
r3g$IL66B_share_of_NL_shift <- r3g$l2_IL66B_vs_prim / ((r3g$l2_NL70B_vs_prim + r3g$l2_NL71B_vs_prim) / 2)
r3g$batch_NL_vs_IL_lfc_mle <- bt$lfc_mle; r3g$batch_NL_vs_IL_padj <- bt$padj
r3g$IL66B_alone_lfc_mle <- i1$lfc_mle; r3g$IL66B_alone_padj <- i1$padj
wt(r3g, "r3_il66b_and_batch_genes.tsv")
bsig <- batch$gene_id[sig(batch)]; isig <- il66$gene_id[sig(il66)]
by_dir <- function(dn) {
  k <- r3g$direction == dn
  list(n = sum(k), IL66B_beyond_all_primaries = sum(r3g$IL66B_beyond_all_primaries[k]),
       NL70B_and_NL71B_beyond_all_primaries = sum(r3g$NL70B_beyond_all_primaries[k] & r3g$NL71B_beyond_all_primaries[k]),
       IL66B_inside_primary_range = sum(r3g$IL66B_inside_primary_range[k]),
       median_IL66B_share_of_NL_shift = median(r3g$IL66B_share_of_NL_shift[k], na.rm = TRUE),
       significant_in_batch_NL_vs_IL = sum(r3g$gene_id[k] %in% bsig),
       significant_in_IL66B_alone_same_sign = sum(k & r3g$gene_id %in% isig & sign(r3g$IL66B_alone_lfc_mle) == dir38))
}
r3 <- list(batch_NL_vs_IL = as.list(cnts(batch)), IL66B_alone_vs_primaries = as.list(cnts(il66)),
           batch_sig_overlap_with_pipeline_genes = sum(g38 %in% bsig),
           up = by_dir("up"), down = by_dir("down"))
wt(data.frame(fit = c("batch_NL_vs_IL", "IL66B_vs_primaries"), rbind(cnts(batch), cnts(il66)), check.names = FALSE),
   "r3_counts.tsv")
# host composition per sample (rat-bin CPM, all nine samples): parenchyma, glia, myeloid, lymphoid, vessel, blood
mk <- c("Snap25", "Syt1", "Rbfox3", "Slc17a7", "Gad1", "Gfap", "Aqp4", "Mbp", "Plp1", "Mog", "Pdgfra", "P2ry12",
        "Cx3cr1", "Aif1", "Cd68", "Ptprc", "Cd3e", "Cd19", "Ms4a1", "Nkg7", "Pecam1", "Cldn5", "Col1a1", "Hbb")
mid <- names(sym)[match(mk, sym)]; mid <- mid[!is.na(mid)]
cpm <- sweep(as.matrix(cnt[, c(tum, ctrl)]), 2, colSums(as.matrix(cnt[, c(tum, ctrl)])), "/") * 1e6
wt(data.frame(gene = unname(sym[mid]), round(cpm[mid, , drop = FALSE], 1), check.names = FALSE), "r3_marker_panel_cpm.tsv")

## R5: human reads beside rat reads (xengsort leakage) ------------------------------------------------------------------
orth <- read.delim(f_orth, stringsAsFactors = FALSE)
hum <- read.delim(f_hum, check.names = FALSE, stringsAsFactors = FALSE)
hum$gid <- sub("\\.[0-9]+$", "", hum$gene_id)
hrow <- function(rid) { h <- orth$human_ensembl_id[orth$rat_ensembl_id == rid]
                        if (!length(h)) return(NA_integer_); m <- match(h[1], hum$gid); m }
r5 <- do.call(rbind, lapply(g38, function(g) {
  m <- hrow(g)
  rr <- unlist(cnt[g, tum]); hh <- if (is.na(m)) rep(NA_real_, 6) else unlist(hum[m, tum])
  data.frame(gene_id = g, gene_name = sym[[g]], human_gene = if (is.na(m)) NA else hum$gene_name[m],
             direction = ifelse(pipe[g, "log2FoldChange"] > 0, "up", "down"),
             t(setNames(rr, paste0("rat_", tum))), t(setNames(hh, paste0("human_", tum))),
             human_ge_rat_any = if (is.na(m)) NA else any(hh >= rr),
             max_human_over_rat = if (is.na(m)) NA else max(hh / pmax(rr, 0.5)),
             check.names = FALSE, stringsAsFactors = FALSE)
}))
wt(r5, "r5_leakage_genes.tsv")
# ribosomal-protein block: does the rat-bin RP share follow the tumour's human share?
rp <- names(sym)[grepl("^Rp[ls][0-9]+[a-z]?[0-9]*$", sym) & names(sym) %in% universe]
lib_rat <- colSums(cnt[, tum]); lib_hum <- colSums(hum[, tum])
rp_h <- na.omit(sapply(rp, hrow))
r5rp <- data.frame(tumour = tum, rat_reads = lib_rat, human_reads = lib_hum,
                   human_share = lib_hum / (lib_hum + lib_rat),
                   rat_RP_fraction = colSums(cnt[rp, tum]) / lib_rat,
                   human_RP_fraction = colSums(hum[rp_h, tum]) / lib_hum, row.names = NULL)
wt(r5rp, "r5_rp_block.tsv")
r5_sum <- list(genes_checked = length(g38), genes_without_1to1_human_ortholog = r5$gene_name[is.na(r5$human_gene)],
               genes_human_ge_rat_in_any_tumour = r5$gene_name[which(r5$human_ge_rat_any)],
               n_rat_RP_genes = length(rp),
               spearman_rat_RP_fraction_vs_human_share = cor(r5rp$rat_RP_fraction, r5rp$human_share, method = "spearman"))

summ <- list(settings = list(deseq2 = as.character(packageVersion("DESeq2")),
                             ashr = if (use_ashr) as.character(packageVersion("ashr")) else NA,
                             R = R.version.string, universe_rule = ">= 10 reads in >= 1 of 9 samples",
                             design = "~ 0 + Classification, avgTxLength offsets",
                             minReplicatesForReplace = 99, alpha = 0.1, lfc_rule = "|ashr lfc| >= 1"),
             pipeline_genes = data.frame(gene = unname(sym[g38]), direction = ifelse(dir38 > 0, "up", "down"),
                                         lfc = round(pipe[g38, "log2FoldChange"], 3), padj = signif(pipe[g38, "padj"], 3)),
             R0_reproduce = r0,
             R1_leave_one_out = list(counts = r1_counts, mmp13 = r1_mmp13,
                                     genes_kept_in_all_6 = sum(gt$n_loo_kept == 6),
                                     genes_kept_in_5_or_more = sum(gt$n_loo_kept >= 5),
                                     genes_lost_only_without = table(gt$lost_without[gt$n_loo_kept == 5])),
             R2_relabel_null = list(summary = r2_sum, splits = r2[, c("label", "true_split", "nl_together", "padj05",
                                                                      "padj05_absLFC1_ashr")]),
             R3_batch_and_IL66B = r3, R5_leakage = list(summary = r5_sum, rp_block = r5rp))
summ$R1_leave_one_out$genes_lost_only_without <- as.list(summ$R1_leave_one_out$genes_lost_only_without)
writeLines(to_json(summ), file.path(out, "summary_de.json"))
cat("step A done\n")
