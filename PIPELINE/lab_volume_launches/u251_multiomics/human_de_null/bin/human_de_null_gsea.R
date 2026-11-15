#!/usr/bin/env Rscript
# human_de_null_gsea.R: step 2 of human_de_null/run.sbatch. A copy of ANALYSIS/host/verify/host_verify_gsea.R adapted
# for the human side; the host original is not edited. Runs in go2432-bioconductor.sif (R 4.3.3, fgsea 1.28.0 in the
# image; installs fgsea into <rlib>/R<x.y> only if the image lacks it).
#
# Preranked gene-set tests (fgseaSimple, 10,000 permutations, sets of 15-500 genes, seed 1234 before every call) on
# every fit's DESeq2 Wald statistic from step 1 (ranks/<fit>.tsv, manifest fits.tsv). Genes are ranked by symbol; a
# symbol carried by several Ensembl ids keeps the id with the largest |stat|.
#
# Gene sets: U/ANALYSIS/refs/pathways/human/combined_human.gmt, the file the published GSEA used. Membership is taken
# from it (first occurrence of a duplicated name). It is split into collections by NAME, using the component files it
# was built from: hallmark (hs.hallmark), kegg (hs.kegg: legacy KEGG and KEGG MEDICUS), reactome and wikipathways (the
# REACTOME_ and WP_ names of hs.c2_curated), c2_other (the rest of hs.c2_curated: BIOCARTA, PID and the CGP author
# sets), gobp (hs.gobp) and brain (hs.brain_gmt_v2). Each collection gets its own BH; 'combined' is the whole file
# with one BH, as the published run tested it.
#
# The ten 3-v-3 relabellings give an empirical null for (a) how many sets reach padj < 0.05 per collection and (b) the
# |NES| of the pre-registered ribosomal-protein / translation sets: the six leading sets of the published GSEA
# (ANALYSIS/gsea_leave_one_out/SUMMARY.md) plus four more that name the same machinery.
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "hv_json.R"))
A <- commandArgs(trailingOnly = TRUE)
arg <- function(k, default = NULL) {
  m <- grep(paste0("^", k, "="), A, value = TRUE)
  if (!length(m)) { if (is.null(default)) stop("missing argument ", k, "="); return(default) }
  sub(paste0("^", k, "="), "", m[1])
}
out     <- arg("out")
pass    <- arg("pass")
gmt_dir <- arg("gmt_dir")
rlib    <- arg("rlib")
cores   <- as.integer(arg("cores", "4"))
broad   <- arg("broad_dir", "none")        # the published Broad GSEA reports (combined collection)
f_md5   <- arg("inputs_md5", "none")
f_cont  <- arg("containers", "none")
cmp     <- arg("compare", "none")          # pass A results directory (pass B only): side-by-side tables
NPERM   <- as.integer(arg("nperm", "10000"))
SMIN    <- as.integer(arg("set_min", "15"))
SMAX    <- as.integer(arg("set_max", "500"))
SEED    <- as.integer(arg("seed", "1234"))
QCUT    <- as.numeric(arg("padj_cut", "0.05"))
t_start <- Sys.time()

lib <- file.path(rlib, paste0("R", R.version$major, ".", sub("\\..*$", "", R.version$minor)))
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
dir.create(file.path(out, "gsea"), showWarnings = FALSE)
wt <- function(df, name) {
  f <- file.path(out, name)
  write.table(df, paste0(f, ".part"), sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(file.rename(paste0(f, ".part"), f))
}
rtsv <- function(f) read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
rr <- function(v, is_true) { t <- v[is_true]; c(rank_min = 1 + sum(v > t, na.rm = TRUE), rank_max = sum(v >= t, na.rm = TRUE)) }

## collections -------------------------------------------------------------------------------------------------------
GMT <- c(combined = "combined_human.gmt", hallmark = "hs.hallmark.v2023.2.gmt", kegg = "hs.kegg.v2023.2.gmt",
         c2 = "hs.c2_curated.v2023.2.gmt", gobp = "hs.gobp.v2023.2.gmt", brain = "hs.brain_gmt_v2.gmt")
raw <- lapply(GMT, function(f) suppressWarnings(gmtPathways(file.path(gmt_dir, f))))
comb <- raw$combined
dup_names <- sum(duplicated(names(comb)))
comb <- comb[!duplicated(names(comb))]
comb <- lapply(comb, function(g) unique(g[!is.na(g) & nzchar(g)]))
nm <- lapply(raw[-1], function(p) unique(names(p)))
reac <- nm$c2[startsWith(nm$c2, "REACTOME_")]; wp <- nm$c2[startsWith(nm$c2, "WP_")]
coll_names <- list(hallmark = nm$hallmark, kegg = nm$kegg, reactome = reac, wikipathways = wp,
                   c2_other = setdiff(nm$c2, c(nm$kegg, reac, wp)), gobp = nm$gobp, brain = nm$brain)
PW <- lapply(coll_names, function(n) comb[names(comb) %in% n])
PW$combined <- comb
colls <- names(PW); split_colls <- setdiff(colls, "combined")
unassigned <- setdiff(names(comb), unlist(coll_names))
coll_info <- data.frame(collection = colls, sets_in_file = vapply(PW, length, 1L),
                        sets_in_size_range = NA_integer_, stringsAsFactors = FALSE)
cat(sprintf("combined_human.gmt: %d names (%d duplicate lines dropped); unassigned to a collection: %d\n",
            length(comb), dup_names, length(unassigned)))

RP6 <- c("KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION", "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
         "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY", "KEGG_RIBOSOME",
         "REACTOME_SELENOAMINO_ACID_METABOLISM", "REACTOME_CELLULAR_RESPONSE_TO_STARVATION")
RPX <- c("WP_CYTOPLASMIC_RIBOSOMAL_PROTEINS", "GOBP_CYTOPLASMIC_TRANSLATION", "REACTOME_TRANSLATION",
         "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION")
rp_all <- c(RP6, RPX)
rp_missing <- setdiff(rp_all, names(comb))
rp_use <- setdiff(rp_all, rp_missing)
home_of <- vapply(rp_use, function(s) { h <- split_colls[vapply(split_colls, function(co) s %in% names(PW[[co]]), TRUE)]
                                        if (length(h)) h[1] else NA_character_ }, "")

## fits ----------------------------------------------------------------------------------------------------------------
man <- rtsv(file.path(out, "fits.tsv"))
for (k in c("true_split", "nl_together", "IL68B_in_target")) man[[k]] <- as.logical(man[[k]])
true_lab <- man$label[man$type == "relabel_3v3" & man$true_split]
load_stats <- function(label) {
  r <- rtsv(file.path(out, "ranks", paste0(label, ".tsv")))
  r <- r[!is.na(r$stat) & !is.na(r$gene_name) & nzchar(r$gene_name), ]
  r <- r[order(-abs(r$stat)), ]; r <- r[!duplicated(r$gene_name), ]
  setNames(r$stat, r$gene_name)
}
true_equals_full <- isTRUE(all.equal(load_stats("full"), load_stats(true_lab)))
WARN <- list()                                 # every fgsea warning, tallied into summary_gsea.json
run1 <- function(label, coll) {
  f <- file.path(out, "gsea", paste0(label, ".", coll, ".tsv.gz"))
  if (file.exists(f)) return(read.delim(gzfile(f), stringsAsFactors = FALSE))
  t0 <- Sys.time()
  set.seed(SEED)
  res <- withCallingHandlers(
    as.data.frame(fgseaSimple(PW[[coll]], load_stats(label), nperm = NPERM, minSize = SMIN, maxSize = SMAX, BPPARAM = BP)),
    warning = function(w) { WARN[[length(WARN) + 1]] <<- paste(coll, gsub("[[:space:]]+", " ", conditionMessage(w)), sep = ": ")
                            invokeRestart("muffleWarning") })
  res$leadingEdge <- vapply(res$leadingEdge, paste, "", collapse = ",")
  res <- res[order(res$pval, -abs(res$NES)), c("pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
  part <- paste0(f, ".part"); gz <- gzfile(part, "w")
  write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz); file.rename(part, f)
  cat(sprintf("gsea %-28s %-12s %5d sets %6.1f s\n", label, coll, nrow(res),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  res
}
R <- list()
for (i in seq_len(nrow(man))) {
  lab <- man$label[i]
  if (lab == true_lab && true_equals_full) next          # the true relabelling is the full fit (checked just above)
  for (co in colls) R[[paste(lab, co, sep = "|")]] <- run1(lab, co)
}
if (true_equals_full) for (co in colls) R[[paste(true_lab, co, sep = "|")]] <- R[[paste("full", co, sep = "|")]]
get <- function(lab, co) R[[paste(lab, co, sep = "|")]]
coll_info$sets_in_size_range <- vapply(colls, function(co) nrow(get("full", co)), 1L)
wt(coll_info, "fgsea_collections.tsv")

## counts per fit and collection; the relabelling null ------------------------------------------------------------
cn <- do.call(rbind, lapply(names(R), function(k) {
  p <- strsplit(k, "|", fixed = TRUE)[[1]]; d <- R[[k]]; m <- man[man$label == p[1], ]
  data.frame(fit = p[1], type = m$type, true_split = m$true_split, nl_together = m$nl_together,
             IL68B_in_target = m$IL68B_in_target, collection = p[2], sets_tested = nrow(d),
             padj05 = sum(d$padj < QCUT, na.rm = TRUE), padj05_up = sum(d$padj < QCUT & d$NES > 0, na.rm = TRUE),
             padj05_down = sum(d$padj < QCUT & d$NES < 0, na.rm = TRUE), padj25 = sum(d$padj < 0.25, na.rm = TRUE),
             min_padj = suppressWarnings(min(d$padj, na.rm = TRUE)), stringsAsFactors = FALSE) }))
wt(cn, "fgsea_counts.tsv")
nul <- do.call(rbind, lapply(colls, function(co) {
  s <- cn[cn$type == "relabel_3v3" & cn$collection == co, ]
  a <- rr(s$padj05, s$true_split); b <- rr(s$padj25, s$true_split); o <- !s$true_split
  data.frame(collection = co, sets_tested = s$sets_tested[s$true_split], true_padj05 = s$padj05[s$true_split],
             rank_min_padj05 = a[["rank_min"]], rank_max_padj05 = a[["rank_max"]],
             p_conservative_padj05 = a[["rank_max"]] / nrow(s), median_other_padj05 = median(s$padj05[o]),
             max_other_padj05 = max(s$padj05[o]),
             max_other_split = paste(s$fit[o][s$padj05[o] == max(s$padj05[o])], collapse = ";"),
             true_padj25 = s$padj25[s$true_split], rank_min_padj25 = b[["rank_min"]], rank_max_padj25 = b[["rank_max"]],
             median_other_padj25 = median(s$padj25[o]), stringsAsFactors = FALSE) }))
wt(nul, "fgsea_relabel_null.tsv")

## the ribosomal-protein / translation sets in every fit ------------------------------------------------------------
rp_rows <- do.call(rbind, lapply(man$label, function(lab) {
  dc <- get(lab, "combined"); mrow <- man[man$label == lab, ]
  do.call(rbind, lapply(rp_use, function(s) {
    m <- dc[match(s, dc$pathway), ]; h <- home_of[[s]]
    dh <- if (!is.na(h)) get(lab, h) else NULL
    ph <- if (!is.null(dh)) dh$padj[match(s, dh$pathway)] else NA_real_
    le <- if (is.na(m$leadingEdge)) 0L else length(strsplit(m$leadingEdge, ",", fixed = TRUE)[[1]])
    data.frame(fit = lab, type = mrow$type, true_split = mrow$true_split, IL68B_in_target = mrow$IL68B_in_target,
               set = s, leading_six = s %in% RP6, size = m$size, NES = m$NES, pval = m$pval, padj_combined = m$padj,
               home_collection = h, padj_home = ph, leading_edge_n = le, stringsAsFactors = FALSE) }))
}))
wt(rp_rows, "fgsea_rp_sets.tsv")
rp_rank_one <- function(s, d) {
  v <- abs(d$NES); a <- rr(v, d$true_split); o <- !d$true_split
  data.frame(set = s, leading_six = s %in% RP6, NES_true = d$NES[d$true_split], padj_home_true = d$padj_home[d$true_split],
             abs_NES_rank_min = a[["rank_min"]], abs_NES_rank_max = a[["rank_max"]],
             p_conservative = a[["rank_max"]] / nrow(d), median_other_abs_NES = median(v[o], na.rm = TRUE),
             max_other_abs_NES = max(v[o], na.rm = TRUE),
             max_other_split = paste(d$fit[o][which(v[o] == max(v[o], na.rm = TRUE))], collapse = ";"),
             splits_padj_home05 = sum(d$padj_home < QCUT, na.rm = TRUE),
             mean_NES_IL68B_in_target = mean(d$NES[d$IL68B_in_target], na.rm = TRUE),
             mean_NES_IL68B_in_reference = mean(d$NES[!d$IL68B_in_target], na.rm = TRUE), stringsAsFactors = FALSE)
}
rs <- rp_rows[rp_rows$type == "relabel_3v3", ]
rp_rank <- do.call(rbind, lapply(rp_use, function(s) rp_rank_one(s, rs[rs$set == s, ])))
six <- aggregate(NES ~ fit + true_split + IL68B_in_target, data = rs[rs$leading_six, ], FUN = mean)
six$padj_home <- NA_real_
rp_rank <- rbind(rp_rank, rp_rank_one("mean_NES_of_the_leading_six", six))
wt(rp_rank, "fgsea_rp_rank.tsv")

## the full fit's top sets: stability without each tumour and how often the relabellings reach them ------------------
loo_labs <- man$label[man$type == "leave_one_out"]
null_labs <- man$label[man$type == "relabel_3v3" & !man$true_split]
top <- do.call(rbind, lapply(colls, function(co) {
  f <- get("full", co)
  up <- head(f[f$NES > 0, ], 10); dn <- head(f[f$NES < 0, ], 10)
  t <- rbind(cbind(direction = rep("up", nrow(up)), up), cbind(direction = rep("down", nrow(dn)), dn))
  if (!nrow(t)) return(NULL)
  lk <- function(lab, pw) { d <- get(lab, co); m <- d[match(pw, d$pathway), ]; cbind(m$NES, m$padj) }
  below <- function(x, th) !is.na(x) & x < th
  t$n_loo_same_sign_padj05 <- 0L; t$n_loo_same_sign_padj25 <- 0L; t$n_null_splits_padj05 <- 0L
  for (l in loo_labs) {
    x <- lk(l, t$pathway); ss <- !is.na(x[, 1]) & sign(x[, 1]) == sign(t$NES)
    t$n_loo_same_sign_padj05 <- t$n_loo_same_sign_padj05 + (ss & below(x[, 2], QCUT))
    t$n_loo_same_sign_padj25 <- t$n_loo_same_sign_padj25 + (ss & below(x[, 2], 0.25))
    t[[paste0("NES_wo_", sub("^loo_", "", l))]] <- round(x[, 1], 3)
  }
  for (l in null_labs) { x <- lk(l, t$pathway); t$n_null_splits_padj05 <- t$n_null_splits_padj05 + below(x[, 2], QCUT) }
  xb <- lk("batch_NL_vs_IL", t$pathway); t$batch_NL_vs_IL_NES <- round(xb[, 1], 3); t$batch_NL_vs_IL_padj <- xb[, 2]
  t$leadingEdge <- substr(t$leadingEdge, 1, 300)
  cbind(collection = co, t, stringsAsFactors = FALSE) }))
wt(top, "fgsea_top_full.tsv")
for (co in colls) wt(get("full", co), paste0("fgsea_full_", co, ".tsv"))
tbf <- do.call(rbind, lapply(names(R), function(k) {
  p <- strsplit(k, "|", fixed = TRUE)[[1]]; d <- R[[k]]
  s <- rbind(head(d[d$NES > 0, ], 5), head(d[d$NES < 0, ], 5))
  if (!nrow(s)) return(NULL)
  data.frame(fit = p[1], collection = p[2], s[, c("pathway", "size", "NES", "pval", "padj")], stringsAsFactors = FALSE) }))
wt(tbf, "fgsea_top_by_fit.tsv")

## beside the published Broad GSEA (Diff_of_Classes on normalised counts, gene-set permutation, combined) ------------
vsb <- NULL
rf <- file.path(broad, c("therapy_impact.combined_human.gsea_report_for_Recurrent_U2.tsv",
                         "therapy_impact.combined_human.gsea_report_for_Primary_U2.tsv"))
if (broad != "none" && all(file.exists(rf))) {
  b <- do.call(rbind, lapply(rf, function(x) { d <- rtsv(x)
    data.frame(key = toupper(d$NAME), NES = suppressWarnings(as.numeric(d$NES)), p = suppressWarnings(as.numeric(d[["NOM p-val"]])),
               fdr = suppressWarnings(as.numeric(d[["FDR q-val"]])), stringsAsFactors = FALSE) }))
  f <- get("full", "combined"); f$key <- toupper(f$pathway)
  m <- merge(f, b, by = "key", suffixes = c("_fgsea", "_broad"))
  topk <- function(d, nes, pv, k = 20, s = 1) { d <- d[!is.na(d[[nes]]) & sign(d[[nes]]) == s, ]
                                                head(d$key[order(d[[pv]], -abs(d[[nes]]))], k) }
  rp_b <- m[m$key %in% toupper(RP6), c("key", "NES_fgsea", "pval", "padj", "NES_broad", "p", "fdr")]
  vsb <- list(sets_matched = nrow(m), sets_fgsea = nrow(f), sets_broad = nrow(b),
              spearman_NES = cor(m$NES_fgsea, m$NES_broad, method = "spearman", use = "complete.obs"),
              top20_up_overlap = length(intersect(topk(m, "NES_fgsea", "pval"), topk(m, "NES_broad", "p"))),
              top20_down_overlap = length(intersect(topk(m, "NES_fgsea", "pval", s = -1), topk(m, "NES_broad", "p", s = -1))),
              broad_fdr25_up = sum(b$fdr < 0.25 & b$NES > 0, na.rm = TRUE),
              broad_fdr25_down = sum(b$fdr < 0.25 & b$NES < 0, na.rm = TRUE),
              leading_six_side_by_side = rp_b)
  wt(m[, c("key", "size", "NES_fgsea", "pval", "padj", "NES_broad", "p", "fdr")], "fgsea_vs_broad_combined.tsv")
}

## splits_summary.tsv gets the gene-set columns; true_split_rank.tsv collects every rank ----------------------------
ss <- rtsv(file.path(out, "splits_summary.tsv"))
for (co in colls) { s <- cn[cn$type == "relabel_3v3" & cn$collection == co, ]; ss[[paste0("sets05_", co)]] <- s$padj05[match(ss$label, s$fit)] }
for (s in intersect(RP6, rp_use)) { d <- rs[rs$set == s, ]; ss[[paste0("NES_", s)]] <- round(d$NES[match(ss$label, d$fit)], 4) }
wt(ss, "splits_summary.tsv")
rk_de <- rtsv(file.path(out, "true_split_rank_de.tsv"))
rk_gs <- data.frame(pass = pass, family = "gene_sets", metric = paste0("sets_padj05_", nul$collection),
                    value_true = nul$true_padj05, rank_min = nul$rank_min_padj05, rank_max = nul$rank_max_padj05,
                    n_splits = 10L, p_conservative = nul$p_conservative_padj05, median_other = nul$median_other_padj05,
                    max_other = nul$max_other_padj05, max_other_split = nul$max_other_split, stringsAsFactors = FALSE)
rk_rp <- data.frame(pass = pass, family = "rp_sets_abs_NES", metric = rp_rank$set, value_true = abs(rp_rank$NES_true),
                    rank_min = rp_rank$abs_NES_rank_min, rank_max = rp_rank$abs_NES_rank_max, n_splits = 10L,
                    p_conservative = rp_rank$p_conservative, median_other = rp_rank$median_other_abs_NES,
                    max_other = rp_rank$max_other_abs_NES, max_other_split = rp_rank$max_other_split, stringsAsFactors = FALSE)
trk <- rbind(rk_de, rk_gs, rk_rp)
wt(trk, "true_split_rank.tsv")

## pass B: side by side with pass A (the analysis of record) --------------------------------------------------------
sbs_sum <- NULL
if (cmp != "none" && file.exists(file.path(cmp, "splits_summary.tsv"))) {
  keys <- c("label", "reference", "target", "true_split", "nl_together", "IL68B_in_target")
  a_ss <- rtsv(file.path(cmp, "splits_summary.tsv"))
  sbs <- merge(a_ss, ss, by = keys, suffixes = c("_passA", "_passB"))
  wt(sbs[order(!sbs$true_split, -sbs$padj05_passA), ], file.path("..", "side_by_side_splits.tsv"))
  if (file.exists(file.path(cmp, "true_split_rank.tsv"))) {
    a_rk <- rtsv(file.path(cmp, "true_split_rank.tsv"))
    srk <- merge(a_rk[, setdiff(names(a_rk), "pass")], trk[, setdiff(names(trk), "pass")], by = c("family", "metric"),
                 suffixes = c("_passA", "_passB"), all = TRUE)
    wt(srk, file.path("..", "side_by_side_ranks.tsv"))
  }
  fa <- read.delim(gzfile(file.path(cmp, "fits", "full.tsv.gz")), stringsAsFactors = FALSE)
  fb <- read.delim(gzfile(file.path(out, "fits", "full.tsv.gz")), stringsAsFactors = FALSE)
  mm <- merge(fa, fb, by = "gene_id", suffixes = c("_A", "_B"))
  sA <- fa$gene_id[which(fa$padj < QCUT)]; sB <- fb$gene_id[which(fb$padj < QCUT)]
  fA <- fa$gene_id[which(fa$padj < QCUT & abs(fa$lfc_ashr) >= log2(1.5))]
  fB <- fb$gene_id[which(fb$padj < QCUT & abs(fb$lfc_ashr) >= log2(1.5))]
  sbs_sum <- list(genes_in_both_universes = nrow(mm), universe_A = nrow(fa), universe_B = nrow(fb),
                  spearman_stat_full_fit = cor(mm$stat_A, mm$stat_B, method = "spearman", use = "complete.obs"),
                  padj05_A = length(sA), padj05_B = length(sB), padj05_both = length(intersect(sA, sB)),
                  filtered_A = length(fA), filtered_B = length(fB), filtered_both = length(intersect(fA, fB)))
}

## summaries ------------------------------------------------------------------------------------------------------
pick <- function(fam, met) trk[trk$family == fam & trk$metric == met, ]
say <- function(fam, met, what) {
  r <- pick(fam, met); if (!nrow(r)) return(NA_character_)
  sprintf("%s: true split %s, rank %d-%d of %d (conservative p = %.2f); median of the other nine %s",
          what, format(signif(r$value_true, 4)), r$rank_min, r$rank_max, r$n_splits, r$p_conservative,
          format(signif(r$median_other, 4)))
}
verdict <- list(
  genes_padj05 = say("DE", "padj05", "genes at padj < 0.05"),
  genes_published_filter = say("DE", "filtered", "genes passing the published filter"),
  rp_block_abs_mean_lfc = say("DE", "abs_rp_mean_lfc_mle", "|mean MLE log2FC| of the RP block"),
  sets_combined = say("gene_sets", "sets_padj05_combined", "sets at padj < 0.05, combined collection"),
  sets_hallmark = say("gene_sets", "sets_padj05_hallmark", "sets at padj < 0.05, hallmark"),
  sets_kegg = say("gene_sets", "sets_padj05_kegg", "sets at padj < 0.05, kegg"),
  sets_reactome = say("gene_sets", "sets_padj05_reactome", "sets at padj < 0.05, reactome"),
  sets_gobp = say("gene_sets", "sets_padj05_gobp", "sets at padj < 0.05, gobp"),
  rp_six_mean_abs_NES = say("rp_sets_abs_NES", "mean_NES_of_the_leading_six", "|mean NES| of the six leading RP/translation sets"),
  kegg_medicus_translation_initiation = say("rp_sets_abs_NES", "KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION",
                                            "|NES| of KEGG_MEDICUS_REFERENCE_TRANSLATION_INITIATION"),
  floor = "10 balanced relabellings: rank 1 of 10 is p = 0.10, the smallest p this design can give")
gs <- list(method = list(fgsea = as.character(packageVersion("fgsea")), installed_this_run = installed_now,
                         library = lib, R = R.version.string, function_used = "fgseaSimple", nperm = NPERM,
                         set_size = sprintf("%d-%d", SMIN, SMAX), statistic = "DESeq2 Wald stat (unshrunk)", seed = SEED,
                         collections_tested_separately = split_colls, combined_also_tested = TRUE,
                         duplicate_set_lines_dropped = dup_names, sets_unassigned_to_a_collection = length(unassigned),
                         rp_sets_missing = rp_missing, true_split_stats_equal_full_fit = true_equals_full,
                         fgsea_warnings = if (length(WARN)) as.list(table(unlist(WARN))) else list()),
           collections = coll_info, relabel_null = nul, rp_rank = rp_rank,
           counts_non_relabel = cn[cn$type != "relabel_3v3", ], vs_broad = vsb, side_by_side_full_fit = sbs_sum,
           runtime_min = round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2))
gs_txt <- to_json(gs, 1)
writeLines(gs_txt, file.path(out, "summary_gsea.json"))
rd_txt <- function(f) if (file.exists(f)) paste(readLines(f, warn = FALSE), collapse = "\n") else "null"
md5 <- if (f_md5 != "none" && file.exists(f_md5)) rtsv(f_md5) else NULL
cont <- if (f_cont != "none" && file.exists(f_cont)) rtsv(f_cont) else NULL
writeLines(paste0("{\n\"pass\": \"", pass, "\",\n\"written\": \"", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                  "\",\n\"verdict\": ", to_json(verdict, 1), ",\n\"true_split_rank\": ", to_json(trk, 1),
                  ",\n\"inputs\": ", to_json(md5, 1), ",\n\"containers\": ", to_json(cont, 1),
                  ",\n\"r0_check\": ", rd_txt(file.path(out, "r0_check.json")),
                  ",\n\"de\": ", rd_txt(file.path(out, "summary_de.json")), ",\n\"gsea\": ", gs_txt, "\n}"),
           file.path(out, "summary.json"))
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_gsea.txt"))
if (length(WARN)) { cat("fgsea warnings (message: count):\n"); print(table(unlist(WARN))) }
cat(paste(unlist(verdict), collapse = "\n"), "\n")
cat("step 2 (gene sets) done\n")
