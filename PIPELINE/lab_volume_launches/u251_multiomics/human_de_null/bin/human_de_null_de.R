#!/usr/bin/env Rscript
# human_de_null_de.R: step 1 of human_de_null/run.sbatch. A copy of ANALYSIS/host/verify/host_verify_de.R adapted for
# the human (U251) side; the host original is not edited. Runs in the published pipeline's own DESeq2 container
# (depot.galaxyproject.org bioconductor-deseq2 1.34.0--r41hc247a5b_3: DESeq2 1.34.0, ashr 2.2-54, R 4.1) so that R0
# can reproduce the published nf-core/differentialabundance 1.5.0 table exactly.
#
#   Rscript human_de_null_de.R key=value ...      (run.sbatch passes every key; see the arg() calls below)
#
# Settings are the published human run's (therapy_v3_params.yaml; results_therapy_v3/pipeline_info/
# params_2026-05-24_19-23-04.json) as implemented by the pipeline's DESEQ2_DIFFERENTIAL module:
#   filter   >= 10 reads (raw salmon counts, not rounded) in >= 1 of the SIX tumours. The published sample sheet holds
#            only the six, so CUSTOM_MATRIXFILTER saw six columns; this rule gives exactly the published 19,351 ids
#            (the alternatives - all ten libraries, >= 2 samples, rounded counts - give 16,154-21,024).
#   design   ~ 0 + Classification (cell means); salmon gene lengths as avgTxLength (tximport-style offsets)
#   DESeq    Wald, parametric, sfType ratio, minReplicatesForReplace 99 (no outlier replacement), useT FALSE
#   results  alpha 0.1, independent filtering, BH, minmu 0.5, lfcThreshold 0, altHypothesis greaterAbs
#   lfc      lfcShrink(type = "ashr"); pvalue and padj stay those of the unshrunk Wald test
#   filtered padj < 0.05 and |log2FoldChange| >= log2(1.5), on the table's (ashr) log2FoldChange. Derived from the
#            published tables: this rule on the published full table gives exactly the 82 published filtered rows.
# Every other fit (leave-one-out, the ten relabellings, NL pair, IL66B alone) uses the same universe and settings.
#
# R0 reproduce the published table (pass A only; if it fails under r0=enforce the script stops with status 10 before
#    any other fit); R1 leave one tumour out; R2 the ten balanced 3-v-3 relabellings (IL67B fixed in the reference
#    group, so no split appears with its mirror); R3 NL70B+NL71B against the other four, and IL66B alone against the
#    primaries. Each fit's Wald statistic goes to ranks/<fit>.tsv for step 2 (human_de_null_gsea.R).
#
# Direction: in every fit the group WITHOUT IL67B is the target, so lfc = target over reference. For the true split the
# target is the recurrent arm. A relabelling's sign is arbitrary, so every rank below uses counts or absolute values.
suppressPackageStartupMessages({ library(DESeq2); library(BiocParallel) })
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "hv_json.R"))

A <- commandArgs(trailingOnly = TRUE)
arg <- function(k, default = NULL) {
  m <- grep(paste0("^", k, "="), A, value = TRUE)
  if (!length(m)) { if (is.null(default)) stop("missing argument ", k, "="); return(default) }
  sub(paste0("^", k, "="), "", m[1])
}
out      <- arg("out")
pass     <- arg("pass")                    # A (published graft+both counts) or B (graft-only counts)
f_sheet  <- arg("samplesheet")             # sample,passA_column,passB_column
f_meta   <- arg("metadata")                # ANALYSIS/metadata_therapy.csv (sample order and arm)
f_counts <- arg("counts")
f_len    <- arg("lengths")
f_gmt    <- arg("gmt")                     # combined_human.gmt (only the RP-block set is read here)
f_pub    <- arg("pub")                     # published full table (R0 in pass A; reference gene lists in both)
f_pubf   <- arg("pub_filtered")            # published filtered table (82 rows)
f_sf     <- arg("pub_sf", "none")          # published size factors (reported beside ours, not gated)
r0_mode  <- arg("r0", "enforce")           # enforce | report | skip
cores    <- as.integer(arg("cores", "4"))
use_ashr <- arg("ashr", "yes") == "yes"
MIN_AB   <- as.numeric(arg("min_abundance", "10"))
MIN_N    <- as.integer(arg("min_samples", "1"))
ALPHA    <- as.numeric(arg("alpha", "0.1"))
QCUT     <- as.numeric(arg("padj_cut", "0.05"))
MINFC    <- as.numeric(arg("min_fc", "1.5"))
TOL      <- as.numeric(arg("r0_tol", "1e-4"))
RP_SET   <- arg("rp_set", "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION")
EXPECT   <- arg("expect_deseq2", "none")    # the published run's DESeq2; a library from home would change it silently
L2FC     <- log2(MINFC)
if (use_ashr && !requireNamespace("ashr", quietly = TRUE)) stop("ashr is not installed in this container")
if (EXPECT != "none" && as.character(packageVersion("DESeq2")) != EXPECT)
  stop("DESeq2 ", as.character(packageVersion("DESeq2")), " loaded from ", dirname(find.package("DESeq2")),
       "; expected ", EXPECT, " (the published run's)")
cat("DESeq2", as.character(packageVersion("DESeq2")), "| ashr", if (use_ashr) as.character(packageVersion("ashr")) else "not used",
    "| libPaths", paste(.libPaths(), collapse = ":"), "\n")
for (d in c("fits", "ranks")) dir.create(file.path(out, d), recursive = TRUE, showWarnings = FALSE)
BP <- if (cores > 1) MulticoreParam(cores) else SerialParam()
set.seed(1234)
wt <- function(df, name) {
  f <- file.path(out, name)
  write.table(df, paste0(f, ".part"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file.rename(paste0(f, ".part"), f))
}
t_start <- Sys.time()

## samples ------------------------------------------------------------------------------------------------------------
sh <- read.csv(f_sheet, stringsAsFactors = FALSE, check.names = FALSE)
md <- read.csv(f_meta, stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(sh$sample, md$sample)) stop("samplesheet order differs from ", f_meta)
tum   <- sh$sample
prim  <- tum[md$Classification == "Primary_U2"]
recur <- tum[md$Classification == "Recurrent_U2"]
if (length(tum) != 6 || length(prim) != 3 || length(recur) != 3 || !("IL67B" %in% prim))
  stop("expected 3 Primary_U2 (with IL67B) and 3 Recurrent_U2 in ", f_meta)
col_of <- setNames(sh[[paste0("pass", pass, "_column")]], tum)
cat("pass", pass, "| primary", paste(prim, collapse = ","), "| recurrent", paste(recur, collapse = ","),
    "| count columns", paste(col_of, collapse = ","), "\n")

## counts, universe, lengths -------------------------------------------------------------------------------------------
cnt <- read.delim(f_counts, check.names = FALSE, stringsAsFactors = FALSE)
miss <- setdiff(col_of, colnames(cnt))
if (length(miss)) stop("count columns missing from ", f_counts, ": ", paste(miss, collapse = ","))
if (anyDuplicated(cnt$gene_id)) stop("duplicate gene_id in ", f_counts)
rownames(cnt) <- cnt$gene_id
X <- as.matrix(cnt[, col_of]); colnames(X) <- tum
if (any(!is.finite(X)) || any(X < 0)) stop("counts contain NA, Inf or negative values")
sym <- setNames(ifelse(is.na(cnt$gene_name) | cnt$gene_name == "", cnt$gene_id, cnt$gene_name), cnt$gene_id)
universe <- rownames(X)[rowSums(X >= MIN_AB) >= MIN_N]
univ_alt <- c(rule_used = length(universe),
              min_samples_2 = sum(rowSums(X >= MIN_AB) >= 2), min_samples_3 = sum(rowSums(X >= MIN_AB) >= 3),
              rounded_counts = sum(rowSums(round(X) >= MIN_AB) >= MIN_N))
len <- read.delim(f_len, check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(len$gene_id)) stop("duplicate gene_id in ", f_len)
rownames(len) <- len$gene_id
missL <- setdiff(universe, rownames(len))
if (length(missL)) stop(length(missL), " universe genes have no row in ", f_len)
missC <- setdiff(col_of, colnames(len))
if (length(missC)) stop("length columns missing from ", f_len, ": ", paste(missC, collapse = ","))
LEN <- as.matrix(len[universe, col_of]); colnames(LEN) <- tum
nbad <- sum(!is.finite(LEN) | LEN <= 0)
if (nbad) stop(nbad, " non-positive or missing gene lengths in the universe")
cat(sprintf("%d genes in the count table; universe %d (>= %g reads in >= %d of the six tumours)\n",
            nrow(X), length(universe), MIN_AB, MIN_N))

## the ribosomal-protein block (the published 'translation' fall): members of RP_SET, read from the GMT, never typed
gl <- grep(paste0("^", RP_SET, "\t"), readLines(f_gmt, warn = FALSE), value = TRUE, useBytes = TRUE)
if (!length(gl)) stop("set ", RP_SET, " is not in ", f_gmt)
rp_sym <- unique(strsplit(gl[1], "\t", fixed = TRUE)[[1]][-(1:2)]); rp_sym <- rp_sym[nzchar(rp_sym)]
rp_ids <- universe[unname(sym[universe]) %in% rp_sym]
cat(sprintf("RP block: %d symbols in %s, %d ids in the universe\n", length(rp_sym), RP_SET, length(rp_ids)))

## one DESeq2 fit: `ref` samples are the reference level, `tgt` the target; lfc = target over reference ---------------
fit <- function(label, ref, tgt, design = ~ 0 + Classification, keep_dds = FALSE) {
  f <- file.path(out, "fits", paste0(label, ".tsv.gz"))
  if (file.exists(f) && !keep_dds) return(read.delim(gzfile(f), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  smp <- tum[tum %in% c(ref, tgt)]
  cls <- factor(ifelse(smp %in% tgt, "Recurrent_U2", "Primary_U2"), levels = c("Primary_U2", "Recurrent_U2"))
  cd <- data.frame(sample = smp, Classification = cls, row.names = smp)
  dds <- DESeqDataSetFromMatrix(round(X[universe, smp, drop = FALSE]), cd, design)
  L <- LEN[, smp, drop = FALSE]; dimnames(L) <- dimnames(dds)
  assays(dds)[["avgTxLength"]] <- L
  dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = 99, useT = FALSE,
               sfType = "ratio", parallel = cores > 1, BPPARAM = BP, quiet = TRUE)
  ctr <- c("Classification", "Recurrent_U2", "Primary_U2")
  r <- results(dds, lfcThreshold = 0, altHypothesis = "greaterAbs", independentFiltering = TRUE, alpha = ALPHA,
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
  part <- paste0(f, ".part"); gz <- gzfile(part, "w")
  write.table(df, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz); file.rename(part, f)
  wt(df[!is.na(df$stat), c("gene_id", "gene_name", "stat")], file.path("ranks", paste0(label, ".tsv")))
  cat(sprintf("fit %-28s ref %-26s tgt %-26s %5.1f s\n", label, paste(ref, collapse = ","),
              paste(tgt, collapse = ","), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (keep_dds) attr(df, "dds") <- dds
  df
}

lfc_tab <- function(df) if (use_ashr) df$lfc_ashr else df$lfc_mle     # the published table's column is ashr (R0)
sig  <- function(df) !is.na(df$padj) & df$padj < QCUT
filt <- function(df) sig(df) & !is.na(lfc_tab(df)) & abs(lfc_tab(df)) >= L2FC
rpm <- function(df) {
  d <- df[df$gene_id %in% rp_ids, ]
  c(rp_ids = nrow(d), rp_mean_lfc_mle = mean(d$lfc_mle, na.rm = TRUE),
    rp_median_lfc_mle = median(d$lfc_mle, na.rm = TRUE),
    rp_mean_lfc_ashr = if (use_ashr) mean(d$lfc_ashr, na.rm = TRUE) else NA_real_,
    rp_median_stat = median(d$stat, na.rm = TRUE), rp_padj05 = sum(sig(d)), rp_padj05_down = sum(sig(d) & d$lfc_mle < 0))
}
cnts <- function(df) {
  s <- sig(df); f <- filt(df)
  c(genes_with_padj = sum(!is.na(df$padj)), padj05 = sum(s), padj05_up = sum(s & df$lfc_mle > 0),
    padj05_down = sum(s & df$lfc_mle < 0), filtered = sum(f), filtered_up = sum(f & lfc_tab(df) > 0),
    filtered_down = sum(f & lfc_tab(df) < 0), filtered_on_mle = sum(s & abs(df$lfc_mle) >= L2FC),
    padj_below_alpha = sum(!is.na(df$padj) & df$padj < ALPHA), rpm(df))
}

## the published tables: reference gene lists in both passes; R0 in pass A ------------------------------------------
pub  <- read.delim(f_pub, stringsAsFactors = FALSE); rownames(pub) <- pub$gene_id
pubf <- read.delim(f_pubf, stringsAsFactors = FALSE)
g_ref <- pubf$gene_id                                     # the 82 published filtered genes
dir_ref <- sign(pub[g_ref, "log2FoldChange"])
g_p05 <- pub$gene_id[!is.na(pub$padj) & pub$padj < QCUT]  # the published padj < 0.05 genes
dir_p05 <- sign(pub[g_p05, "log2FoldChange"])

manifest <- data.frame(label = character(), type = character(), reference = character(), target = character(),
                       true_split = logical(), nl_together = logical(), IL68B_in_target = logical(),
                       stringsAsFactors = FALSE)
addm <- function(label, type, ref, tgt, true_split = FALSE) {
  nl_tog <- (all(c("NL70B", "NL71B") %in% ref) || all(c("NL70B", "NL71B") %in% tgt))
  manifest[nrow(manifest) + 1, ] <<- list(label, type, paste(ref, collapse = ","), paste(tgt, collapse = ","),
                                          true_split, nl_tog, "IL68B" %in% tgt)
}

full <- fit("full", prim, recur, keep_dds = TRUE); dds_full <- attr(full, "dds"); attr(full, "dds") <- NULL
addm("full", "full", prim, recur, TRUE)

dmax <- function(x, y) { ok <- !is.na(x) & !is.na(y); if (!any(ok)) NA_real_ else max(abs(x[ok] - y[ok])) }
r0 <- list(mode = r0_mode, pass = NA)
if (r0_mode != "skip") {
  cm <- intersect(pub$gene_id, full$gene_id)
  P <- pub[cm, ]; M <- full[match(cm, full$gene_id), ]
  rule_pub <- pub$gene_id[!is.na(pub$padj) & pub$padj < QCUT & abs(pub$log2FoldChange) >= L2FC]
  here_f <- full$gene_id[filt(full)]
  big <- !is.na(P$baseMean) & P$baseMean > 1e-9
  r0 <- list(
    mode = r0_mode, tolerance = TOL,
    published_rows = nrow(pub), universe_here = length(universe),
    universe_rule = sprintf(">= %g raw salmon reads in >= %d of the six tumours", MIN_AB, MIN_N),
    universe_identical = setequal(universe, pub$gene_id) && length(universe) == nrow(pub),
    universe_alternatives = as.list(univ_alt),
    padj_NA_pattern_mismatches = sum(is.na(P$padj) != is.na(M$padj)),
    pvalue_NA_pattern_mismatches = sum(is.na(P$pvalue) != is.na(M$pvalue)),
    max_abs_dpadj = dmax(M$padj, P$padj), max_abs_dpvalue = dmax(M$pvalue, P$pvalue),
    max_abs_dlfc_vs_ashr = dmax(M$lfc_ashr, P$log2FoldChange), max_abs_dlfc_vs_mle = dmax(M$lfc_mle, P$log2FoldChange),
    max_abs_dlfcSE_vs_ashr = dmax(M$lfcSE_ashr, P$lfcSE), max_abs_dlfcSE_vs_mle = dmax(M$lfcSE_mle, P$lfcSE),
    max_rel_dbaseMean = if (any(big)) max(abs(M$baseMean[big] - P$baseMean[big]) / P$baseMean[big]) else NA_real_,
    filter_rule = sprintf("padj < %g and |log2FoldChange| >= log2(%g)", QCUT, MINFC),
    filter_rule_on_published_table = length(rule_pub),
    filter_rule_reproduces_published_filtered = setequal(rule_pub, pubf$gene_id),
    published_filtered_rows = nrow(pubf), filtered_here = length(here_f),
    filtered_identical = setequal(here_f, pubf$gene_id),
    filtered_only_published = unname(sym[setdiff(pubf$gene_id, here_f)]),
    filtered_only_here = unname(sym[setdiff(here_f, pubf$gene_id)]),
    padj05_published = length(g_p05), padj05_here = sum(sig(full)),
    shrinkcall_padj_equals_results_padj = if (use_ashr) isTRUE(all.equal(full$padj, full$padj_shrinkcall)) else NA)
  r0$lfc_column_is <- if (!use_ashr) "unknown (ashr not run)" else
    if (isTRUE(r0$max_abs_dlfc_vs_ashr < TOL)) "ashr-shrunk (lfcShrink type ashr)" else
    if (isTRUE(r0$max_abs_dlfc_vs_mle < TOL)) "unshrunk MLE" else "neither ashr nor MLE matches"
  # size factors, reported beside the pipeline's file (with avgTxLength DESeq2 uses per-gene normalisation factors)
  if (f_sf != "none" && file.exists(f_sf)) {
    sfp <- read.delim(f_sf, stringsAsFactors = FALSE); sfp <- setNames(sfp$sizeFactor, sfp$sample)[tum]
    # A local check found the published values are the median-of-ratios of the universe counts divided by each gene's
    # lengths scaled to their geometric mean (local R: max diff 2.4e-15), i.e. DESeq2's internal length-adjusted size factor
    lenN <- LEN / exp(rowMeans(log(LEN)))
    cand <- list(sizeFactors = sizeFactors(dds_full),
                 mor_on_length_adjusted_counts = estimateSizeFactorsForMatrix(round(X[universe, tum]) / lenN),
                 mor_on_rounded_counts = estimateSizeFactorsForMatrix(counts(dds_full)),
                 geomean_of_normalizationFactors = if (!is.null(normalizationFactors(dds_full)))
                   exp(colMeans(log(normalizationFactors(dds_full)))) else NULL)
    r0$sizefactor_candidates_max_abs_diff <- lapply(cand[!vapply(cand, is.null, TRUE)],
                                                    function(v) max(abs(v[tum] - sfp)))
    r0$published_sizefactors <- as.list(sfp)
  }
  r0$pass <- isTRUE(r0$universe_identical) && r0$padj_NA_pattern_mismatches == 0 &&
    isTRUE(r0$max_abs_dpadj < TOL) && isTRUE(r0$max_abs_dlfc_vs_ashr < TOL) &&
    isTRUE(r0$filtered_identical) && isTRUE(r0$filter_rule_reproduces_published_filtered)
  writeLines(to_json(r0), file.path(out, "r0_check.json"))
  cat("R0 pass:", r0$pass, "| lfc column:", r0$lfc_column_is, "| max |dpadj|", format(r0$max_abs_dpadj),
      "| max |dlfc| (ashr)", format(r0$max_abs_dlfc_vs_ashr), "| filtered", r0$filtered_here, "of",
      r0$published_filtered_rows, "\n")
  if (!r0$pass && r0_mode == "enforce") {
    writeLines(to_json(list(pass = pass, stopped = "R0 failed; no other fit was run and nothing else is used",
                            R0_reproduce = r0)), file.path(out, "summary_de.json"))
    cat("R0 FAILED: stopping before R1-R3 (r0=enforce)\n")
    quit(save = "no", status = 10)
  }
}

## R1 leave one out; R2 the ten relabellings; R3 NL pair and IL66B alone --------------------------------------------
full_int <- fit("full_design_intercept", prim, recur, design = ~ Classification)
loo <- list()
for (s in tum) {
  lab <- paste0("loo_", s); loo[[s]] <- fit(lab, setdiff(prim, s), setdiff(recur, s))
  addm(lab, "leave_one_out", setdiff(prim, s), setdiff(recur, s))
}
others <- setdiff(tum, "IL67B")
splits <- list()
for (p in combn(others, 2, simplify = FALSE)) {
  ref <- tum[tum %in% c("IL67B", p)]; tgt <- setdiff(tum, ref)
  lab <- paste0("split_", paste(tgt, collapse = "_"))
  splits[[lab]] <- fit(lab, ref, tgt)
  addm(lab, "relabel_3v3", ref, tgt, true_split = setequal(tgt, recur))
}
stopifnot(length(splits) == 10, sum(manifest$true_split[manifest$type == "relabel_3v3"]) == 1)
nl <- c("NL70B", "NL71B")
batch <- fit("batch_NL_vs_IL", setdiff(tum, nl), nl); addm("batch_NL_vs_IL", "batch_2v4", setdiff(tum, nl), nl)
il66 <- fit("IL66B_vs_primaries", prim, "IL66B"); addm("IL66B_vs_primaries", "single_1v3", prim, "IL66B")
wt(manifest, "fits.tsv")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_de.txt"))

## R1 -----------------------------------------------------------------------------------------------------------------
kept <- function(df, g, dirs) { d <- df[match(g, df$gene_id), ]; !is.na(d$padj) & d$padj < QCUT & sign(d$lfc_mle) == dirs }
r1_counts <- data.frame(fit = c("full", paste0("loo_", tum)), held_out = c("none", tum),
                        do.call(rbind, c(list(cnts(full)), lapply(loo, cnts))), check.names = FALSE)
r1_counts$of_published_filtered_kept_same_sign <- sapply(c(list(full), loo), function(df) sum(kept(df, g_ref, dir_ref)))
r1_counts$of_published_padj05_kept_same_sign <- sapply(c(list(full), loo), function(df) sum(kept(df, g_p05, dir_p05)))
wt(r1_counts, "r1_loo_counts.tsv")
gt <- data.frame(gene_id = g_ref, gene_name = unname(sym[g_ref]), direction = ifelse(dir_ref > 0, "up", "down"),
                 published_lfc = pub[g_ref, "log2FoldChange"], published_padj = pub[g_ref, "padj"],
                 RP_block = g_ref %in% rp_ids, stringsAsFactors = FALSE)
keep_n <- integer(length(g_ref)); lost_in <- character(length(g_ref))
for (s in tum) {
  d <- loo[[s]][match(g_ref, loo[[s]]$gene_id), ]
  gt[[paste0("lfc_ashr_wo_", s)]] <- d$lfc_ashr; gt[[paste0("lfc_mle_wo_", s)]] <- d$lfc_mle
  gt[[paste0("padj_wo_", s)]] <- d$padj
  k <- kept(loo[[s]], g_ref, dir_ref)
  keep_n <- keep_n + k; lost_in[!k] <- paste0(lost_in[!k], ifelse(nzchar(lost_in[!k]), ",", ""), s)
}
gt$n_loo_kept <- keep_n; gt$lost_without <- lost_in
gt <- gt[order(gt$direction, -gt$n_loo_kept, gt$published_padj), ]
wt(gt, "r1_loo_genes.tsv")

## R2 -----------------------------------------------------------------------------------------------------------------
sm <- manifest[manifest$type == "relabel_3v3", ]
r2 <- data.frame(sm[, c("label", "reference", "target", "true_split", "nl_together", "IL68B_in_target")],
                 do.call(rbind, lapply(splits[sm$label], cnts)), check.names = FALSE, stringsAsFactors = FALSE)
true_lab <- sm$label[sm$true_split]
r2$spearman_stat_vs_true <- sapply(sm$label, function(l) {
  d <- splits[[l]]; cor(d$stat, full$stat[match(d$gene_id, full$gene_id)], method = "spearman", use = "complete.obs") })
r2$of_published_filtered_same_sign <- sapply(sm$label, function(l) sum(kept(splits[[l]], g_ref, dir_ref)))
r2$abs_rp_mean_lfc_mle <- abs(r2$rp_mean_lfc_mle); r2$abs_rp_median_stat <- abs(r2$rp_median_stat)
rank_range <- function(col) {
  v <- r2[[col]]; t <- v[r2$true_split]
  data.frame(metric = col, value_true = t, rank_min = 1 + sum(v > t), rank_max = sum(v >= t), n_splits = length(v),
             p_conservative = sum(v >= t) / length(v), median_other = median(v[!r2$true_split]),
             max_other = max(v[!r2$true_split]), max_other_split = paste(r2$label[!r2$true_split][v[!r2$true_split] ==
                                                                          max(v[!r2$true_split])], collapse = ";"),
             stringsAsFactors = FALSE)
}
rk <- do.call(rbind, lapply(c("padj05", "filtered", "padj_below_alpha", "abs_rp_mean_lfc_mle", "abs_rp_median_stat",
                              "rp_padj05"), rank_range))
rk <- cbind(pass = pass, family = "DE", rk)
wt(rk, "true_split_rank_de.tsv")
r2_top <- do.call(rbind, lapply(sm$label, function(l) {
  d <- splits[[l]]; d <- d[!is.na(d$stat), ]
  data.frame(label = l, padj05 = sum(sig(d)), filtered = sum(filt(d)),
             top_up_in_target = paste(head(d$gene_name[order(-d$stat)], 15), collapse = ","),
             top_down_in_target = paste(head(d$gene_name[order(d$stat)], 15), collapse = ","), stringsAsFactors = FALSE) }))
wt(r2_top[order(-r2_top$padj05), ], "r2_split_top_genes.tsv")
wt(r2[order(-r2$padj05), ], "splits_summary.tsv")     # step 2 adds the gene-set columns to this file
r2_sum <- list(
  n_splits = nrow(r2), true_split = true_lab,
  true_split_equals_full_fit = isTRUE(all.equal(splits[[true_lab]]$stat, full$stat)),
  ranks = rk, nl_together_splits = r2$label[r2$nl_together],
  floor_note = "10 balanced splits: rank 1 of 10 is p = 0.10, the smallest p this design allows")

## R3 -----------------------------------------------------------------------------------------------------------------
nc <- counts(dds_full, normalized = TRUE)
g3 <- intersect(g_ref, rownames(nc)); d3 <- sign(pub[g3, "log2FoldChange"])
pm <- rowMeans(nc[g3, prim, drop = FALSE])
l2 <- function(x) log2((x + 1) / (pm + 1))
bt <- batch[match(g3, batch$gene_id), ]; i1 <- il66[match(g3, il66$gene_id), ]
r3g <- data.frame(gene_id = g3, gene_name = unname(sym[g3]), direction = ifelse(d3 > 0, "up", "down"),
                  published_lfc = pub[g3, "log2FoldChange"], round(nc[g3, tum, drop = FALSE], 1),
                  l2_IL66B_vs_prim = l2(nc[g3, "IL66B"]), l2_NL70B_vs_prim = l2(nc[g3, "NL70B"]),
                  l2_NL71B_vs_prim = l2(nc[g3, "NL71B"]), check.names = FALSE, stringsAsFactors = FALSE)
pmin_ <- apply(nc[g3, prim, drop = FALSE], 1, min); pmax_ <- apply(nc[g3, prim, drop = FALSE], 1, max)
beyond <- function(x) ifelse(d3 > 0, x > pmax_, x < pmin_)     # outside every primary, in the recurrence direction
for (s in recur) r3g[[paste0(s, "_beyond_all_primaries")]] <- beyond(nc[g3, s])
r3g$IL66B_inside_primary_range <- nc[g3, "IL66B"] >= pmin_ & nc[g3, "IL66B"] <= pmax_
r3g$IL66B_share_of_NL_shift <- r3g$l2_IL66B_vs_prim / ((r3g$l2_NL70B_vs_prim + r3g$l2_NL71B_vs_prim) / 2)
r3g$batch_NL_vs_IL_lfc_mle <- bt$lfc_mle; r3g$batch_NL_vs_IL_padj <- bt$padj
r3g$IL66B_alone_lfc_mle <- i1$lfc_mle; r3g$IL66B_alone_padj <- i1$padj
wt(r3g, "r3_il66b_and_nl_genes.tsv")
bsig <- batch$gene_id[sig(batch)]; isig <- il66$gene_id[sig(il66)]
by_dir <- function(dn) {
  k <- r3g$direction == dn
  list(n = sum(k), IL66B_beyond_all_primaries = sum(r3g$IL66B_beyond_all_primaries[k]),
       NL70B_and_NL71B_beyond_all_primaries = sum(r3g$NL70B_beyond_all_primaries[k] & r3g$NL71B_beyond_all_primaries[k]),
       IL66B_inside_primary_range = sum(r3g$IL66B_inside_primary_range[k]),
       median_IL66B_share_of_NL_shift = median(r3g$IL66B_share_of_NL_shift[k], na.rm = TRUE),
       significant_in_batch_NL_vs_IL = sum(r3g$gene_id[k] %in% bsig),
       significant_in_IL66B_alone_same_sign = sum(k & r3g$gene_id %in% isig & sign(r3g$IL66B_alone_lfc_mle) == d3))
}
r3 <- list(batch_NL_vs_IL = as.list(cnts(batch)), IL66B_alone_vs_primaries = as.list(cnts(il66)),
           true_split = as.list(cnts(full)),
           batch_sig_overlap_with_published_filtered = sum(g3 %in% bsig),
           published_filtered_up = by_dir("up"), published_filtered_down = by_dir("down"))
wt(data.frame(fit = c("full", "batch_NL_vs_IL", "IL66B_vs_primaries"), rbind(cnts(full), cnts(batch), cnts(il66)),
              check.names = FALSE), "r3_counts.tsv")
# the RP block per sample (normalised counts, log2 of the ratio to the six-tumour mean): whose level sets the fall
rpn <- nc[intersect(rp_ids, rownames(nc)), tum, drop = FALSE]
rp_sample <- data.frame(sample = tum, arm = ifelse(tum %in% prim, "Primary_U2", "Recurrent_U2"),
                        rp_median_log2_vs_six_mean = apply(log2((rpn + 1) / (rowMeans(rpn) + 1)), 2, median),
                        genome_median_log2_vs_six_mean = apply(log2((nc[, tum] + 1) / (rowMeans(nc[, tum]) + 1)), 2, median),
                        row.names = NULL)
wt(rp_sample, "rp_block_by_sample.tsv")

summ <- list(
  pass = pass,
  settings = list(deseq2 = as.character(packageVersion("DESeq2")),
                  ashr = if (use_ashr) as.character(packageVersion("ashr")) else NA, R = R.version.string,
                  universe_rule = sprintf(">= %g raw reads in >= %d of the six tumours", MIN_AB, MIN_N),
                  universe_n = length(universe), universe_alternatives = as.list(univ_alt),
                  design = "~ 0 + Classification, salmon gene lengths as avgTxLength",
                  test = "Wald, parametric, sfType ratio, minReplicatesForReplace 99, useT FALSE",
                  results = sprintf("alpha %g, independent filtering, BH, minmu 0.5", ALPHA),
                  lfc = if (use_ashr) "lfcShrink ashr; p and padj from the unshrunk Wald test" else "MLE (ashr not run)",
                  filtered_rule = sprintf("padj < %g and |lfc (ashr)| >= log2(%g)", QCUT, MINFC),
                  rp_block = list(set = RP_SET, symbols = length(rp_sym), ids_in_universe = length(rp_ids)),
                  count_columns = as.list(col_of), counts = f_counts, lengths = f_len),
  R0_reproduce = r0,
  R1_leave_one_out = list(counts = r1_counts,
                          published_filtered_kept_in_all_6 = sum(gt$n_loo_kept == 6),
                          published_filtered_kept_in_5_or_more = sum(gt$n_loo_kept >= 5),
                          lost_only_without = as.list(table(gt$lost_without[gt$n_loo_kept == 5]))),
  R2_relabel_null = list(summary = r2_sum,
                         splits = r2[, c("label", "true_split", "nl_together", "IL68B_in_target", "padj05", "filtered",
                                         "rp_mean_lfc_mle", "rp_median_stat", "spearman_stat_vs_true")]),
  R3_nl_pair_and_IL66B = r3,
  RP_block_by_sample = rp_sample,
  runtime_min = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2))
writeLines(to_json(summ), file.path(out, "summary_de.json"))
cat("step 1 (DE) done\n")
