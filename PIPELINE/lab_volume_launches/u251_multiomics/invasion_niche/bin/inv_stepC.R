#!/usr/bin/env Rscript
# invasion_niche step C (part C). Runs in go2432-bioconductor.sif (R 4.3.3) with /rs/rs_grp_oschome/go2432/Rlib/R<x.y>
# first on the library path; GSVA and fgsea are installed there only if the image lacks them.
#   Rscript inv_stepC.R <params.yaml>
# Gene sets come from files only (results/genesets/, each with its source URL and md5). Scores: GSVA (Gaussian kernel)
# on DESeq2 vst of the graft-only counts, never TPM; the same on the published graft+both counts, so that a set moving
# only with the both bin shows. Arm statistic: Welch t on the six tumours, ranked among the ten balanced relabellings.
# fgsea on each relabelling's Wald statistic: human_de_null's fits when present and well formed, else step A's own.
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "inv_common.R"))
P <- read_params(args[1])
OUT <- P$outdir; Cd <- file.path(OUT, "C_programmes"); GS <- file.path(OUT, "genesets")
dir.create(Cd, recursive = TRUE, showWarnings = FALSE)
kv_file <- file.path(OUT, "summary_C.tsv"); cat("part\tkey\tvalue\n", file = kv_file)
LIBS <- splitv(P$libraries); PRIM <- splitv(P$primary); RECUR <- splitv(P$recurrent); TUM <- c(PRIM, RECUR)
CTRL <- splitv(P$controls)
SP <- make_splits(PRIM, RECUR); TL <- true_label(SP)
SEED <- as.integer(P$seed)

lib <- file.path(P$Rlib, paste0("R", R.version$major, ".", sub("\\..*$", "", R.version$minor)))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
installed_now <- character(0)
for (p in c("GSVA", "fgsea")) if (!requireNamespace(p, quietly = TRUE)) {
  message(p, " is not in this image; installing it into ", lib)
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
  BiocManager::install(p, lib = lib, update = FALSE, ask = FALSE)
  installed_now <- c(installed_now, p)
}
suppressPackageStartupMessages({ library(GSVA); library(fgsea); library(BiocParallel) })
kv("C", "R_version_stepC", R.version.string)
kv("C", "GSVA_version", as.character(packageVersion("GSVA"))); kv("C", "fgsea_version", as.character(packageVersion("fgsea")))
kv("C", "installed_into_Rlib_now", if (length(installed_now)) paste(installed_now, collapse = ",") else "none")

## ------------------------------------------------------------------------------------------------ gene sets
sets <- list(); origin <- list()
add_set <- function(name, genes, src, rule) {
  genes <- unique(genes[!is.na(genes) & nzchar(genes)])
  sets[[name]] <<- genes
  origin[[name]] <<- data.frame(set = name, genes = length(genes), source_file = src, rule = rule, stringsAsFactors = FALSE)
}
# Ivy GAP (Puchalski 2018, Table S14 CIBERSORT signature): a gene belongs to the compartment where its value is largest
iv <- read.delim(file.path(GS, "ivygap2018_anatomic4.txt"), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(iv)[1] == "GeneSymbol", setequal(colnames(iv)[-1], c("CT", "CTmvp", "CTpan", "LE")))
am <- colnames(iv)[-1][apply(as.matrix(iv[, -1]), 1, which.max)]
for (cmp in c("LE", "CT", "CTmvp", "CTpan")) add_set(paste0("ivygap_", cmp), iv$GeneSymbol[am == cmp], "ivygap2018_anatomic4.txt",
                                                     paste0("argmax over CT, CTmvp, CTpan, LE == ", cmp))
# cross-check against the frozen spread/signatures copies (sha256 recorded there before any score was computed)
for (cmp in c("LE", "CT", "CTmvp", "CTpan")) {
  f <- file.path(GS, paste0("spread_ivygap_", cmp, ".txt"))
  if (file.exists(f)) kv("C", paste0("ivygap_", cmp, "_identical_to_spread_frozen_set"),
                         setequal(readLines(f, warn = FALSE)[nzchar(readLines(f, warn = FALSE))], sets[[paste0("ivygap_", cmp)]]))
}
# Neftel 2019 Table S2 meta-modules (the published lists)
nf <- read.delim(file.path(GS, "neftel2019_tableS2_metamodules.tsv"), stringsAsFactors = FALSE)
for (s in unique(nf$set)) add_set(paste0("neftel2019_", gsub("[^A-Za-z0-9]", "", s)), nf$gene[nf$set == s], "neftel2019_tableS2_metamodules.tsv",
                                  paste0("Table S2 column ", s))
# Neftel states as derived in this repository from the Neftel SS2 cells (ANALYSIS/refs/neftel_2019 -> CIBERSORT signature)
for (s in c("MES", "AC", "OPC", "NPC")) {
  f <- file.path(GS, paste0("spread_neftel_", s, ".txt"))
  if (file.exists(f)) add_set(paste0("neftel_repo_", s), readLines(f, warn = FALSE), basename(f),
                              "argmax over the repository's Neftel CIBERSORT signature (spread/signatures)")
}
# Venkataramani 2022 Table S2a invasivity-score genes. Table S2e lists 20 genes as the 20 most correlated with the
# invasivity score; all 20 are in the group the file labels "Anticorrelated" (its correlations are with pseudotime),
# so that group is the invasive pole. Stated as read from the file's own legends, not from the paper's text.
vk <- read.delim(file.path(GS, "venkataramani2022_tableS2a_invasivity.tsv"), stringsAsFactors = FALSE)
add_set("venkataramani2022_invasivity_high", vk$gene[vk$direction == "Anticorrelated"], "venkataramani2022_tableS2a_invasivity.tsv",
        "Table S2a genes labelled Anticorrelated (the group holding Table S2e's 20 genes most correlated with invasivity)")
add_set("venkataramani2022_invasivity_low", vk$gene[vk$direction == "Correlated"], "venkataramani2022_tableS2a_invasivity.tsv",
        "Table S2a genes labelled Correlated (the opposite pole)")
# MSigDB: msigdbr (version recorded) when it loads, the MSigDB 2023.2 GMT files otherwise; both compared when both exist
msig_names <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "ANASTASSIOU_MULTICANCER_INVASIVENESS_SIGNATURE")
gmt_read <- function(f) { l <- strsplit(readLines(f, warn = FALSE), "\t"); setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1)) }
gmt <- c(gmt_read(P$msigdb_gmt_hallmark), gmt_read(P$msigdb_gmt_c2))
from_msigdbr <- tryCatch({
  .libPaths(c(.libPaths(), P$msigdbr_lib))
  suppressPackageStartupMessages(library(msigdbr))
  kv("C", "msigdbr_version", as.character(packageVersion("msigdbr")))
  a <- msigdbr(species = "Homo sapiens", collection = "H")
  b <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CGP")
  ab <- rbind(a[, c("gs_name", "gene_symbol")], b[, c("gs_name", "gene_symbol")])
  kv("C", "msigdbr_db_version", tryCatch(paste(unique(a$db_version), collapse = ","), error = function(e) "not reported"))
  setNames(lapply(msig_names, function(n) unique(ab$gene_symbol[ab$gs_name == n])), msig_names)
}, error = function(e) { kv("C", "msigdbr_error", conditionMessage(e)); NULL })
for (n in msig_names) {
  g_gmt <- gmt[[n]]
  if (!is.null(from_msigdbr) && length(from_msigdbr[[n]])) {
    add_set(n, from_msigdbr[[n]], "msigdbr", "msigdbr collection H / C2:CGP")
    if (!is.null(g_gmt)) kv("C", paste0(n, "_msigdbr_vs_gmt2023.2"), sprintf("msigdbr %d, GMT %d, shared %d", length(from_msigdbr[[n]]), length(g_gmt), length(intersect(from_msigdbr[[n]], g_gmt))))
  } else if (!is.null(g_gmt)) {
    add_set(n, g_gmt, basename(if (grepl("^HALLMARK", n)) P$msigdb_gmt_hallmark else P$msigdb_gmt_c2), "MSigDB v2023.2 GMT (msigdbr unavailable)")
  } else kv("C", paste0("missing_set_", n), "not in msigdbr or the GMT files")
}
origin <- do.call(rbind, origin)
PRIMARY_SETS <- splitv(P$primary_sets)
origin$primary <- origin$set %in% PRIMARY_SETS
missing_primary <- setdiff(PRIMARY_SETS, names(sets))
kv("C", "primary_sets", paste(PRIMARY_SETS, collapse = ",")); kv("C", "primary_sets_missing", if (length(missing_primary)) paste(missing_primary, collapse = ",") else "none")
kv("C", "optional_sets_missing", P$optional_sets_missing)
memb <- do.call(rbind, lapply(names(sets), function(n) data.frame(set = n, gene = sets[[n]], stringsAsFactors = FALSE)))
wt(memb, file.path(GS, "sets_membership.tsv"))

## ------------------------------------------------------------------------------------------------ scoring
rd <- function(f) {
  d <- read.delim(gzfile(f), check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(d[, -(1:2)]); rownames(m) <- d$gene_id
  # rows keyed by gene symbol; a symbol carried by several ids keeps the id with the highest mean
  nm <- d$gene_name; o <- order(-rowMeans(m)); m <- m[o, , drop = FALSE]; nm <- nm[o]
  keep <- !duplicated(nm) & !is.na(nm) & nzchar(nm); m <- m[keep, , drop = FALSE]; rownames(m) <- nm[keep]; m
}
gsva_scores <- function(x) {
  x <- x[apply(x, 1, var) > 0, , drop = FALSE]
  mn <- as.integer(P$gsva_min_size); mx <- as.integer(P$gsva_max_size)
  s <- tryCatch(gsva(gsvaParam(x, sets, kcdf = "Gaussian", minSize = mn, maxSize = mx), verbose = FALSE),
                error = function(e) gsva(x, sets, method = "gsva", kcdf = "Gaussian", min.sz = mn, max.sz = mx, verbose = FALSE))
  as.matrix(s)
}
xe <- read.delim(file.path(OUT, "inputs", "xengsort_bins.tsv"), stringsAsFactors = FALSE); rownames(xe) <- xe$lib
m_both <- setNames(xe$m_both_share, xe$lib); graft_pct <- setNames(xe$graft_pct, xe$lib)
matched <- list(); scores <- list(); summ <- list()
for (src in c("graft_only", "graft_plus_both")) for (set in c("all10", "six")) for (kind in c("vst", "logcpm")) {
  f <- file.path(Cd, "expr", sprintf("%s_%s_%s.tsv.gz", kind, src, if (kind == "logcpm") "all10" else set))
  if (kind == "logcpm" && set == "six") next
  if (!file.exists(f)) { kv("C", paste("missing_input", kind, src, set), f); next }
  x <- rd(f)
  if (kind == "logcpm" && set == "all10") x <- x[, LIBS]
  sc <- gsva_scores(x)
  key <- paste(kind, src, set, sep = ":")
  scores[[key]] <- sc
  wt(data.frame(set = rownames(sc), sc, check.names = FALSE), file.path(Cd, "gsva", sprintf("gsva_%s_%s_%s.tsv", kind, src, set)))
  if (kind == "vst" && set == "all10") matched[[src]] <- vapply(names(sets), function(n) sum(sets[[n]] %in% rownames(x)), 0)
  for (n in rownames(sc)) {
    v <- sc[n, ]; st <- split_stats(v, SP); r <- rank_true(st, TL)
    row <- data.frame(scores = kind, counts = src, scored_libraries = set, set = n, primary_set = n %in% PRIMARY_SETS,
                      mean_primary = mean(v[PRIM]), mean_recurrent = mean(v[RECUR]),
                      diff_recurrent_minus_primary = mean(v[RECUR]) - mean(v[PRIM]), welch_t_true = st[TL],
                      rank_min = r[1], rank_max = r[2], rank = fmt_rank(r), p_floor = r[1] / 10,
                      spearman_with_graft_pct_six = suppressWarnings(cor(v[TUM], graft_pct[TUM], method = "spearman")),
                      stringsAsFactors = FALSE)
    if (set == "all10") {
      ct <- suppressWarnings(cor.test(v[LIBS], m_both[LIBS], method = "spearman", exact = TRUE))
      row$C2B <- v["C2B"]; row$IL64B <- v["IL64B"]; row$N168B <- v["N168B"]; row$N269B <- v["N269B"]
      row$tumour_min <- min(v[TUM]); row$tumour_max <- max(v[TUM])
      row$high_in_controls <- max(v[CTRL]) > max(v[TUM])
      row$spearman_with_m_all10 <- unname(ct$estimate); row$spearman_p_with_m <- ct$p.value
      row$tracks_m <- ct$p.value < 0.05
    }
    summ[[length(summ) + 1]] <- cbind(row, t(st))
  }
}
bind_fill <- function(l) { cols <- unique(unlist(lapply(l, colnames)))
  do.call(rbind, lapply(l, function(d) { for (c in setdiff(cols, colnames(d))) d[[c]] <- NA; d[, cols] })) }
summ <- bind_fill(summ); rownames(summ) <- NULL
# leak verdict per set (vst, all ten scored): high in the near-empty controls, tracks the both-bin share, or an arm
# difference present on graft+both that is gone on graft-only
lk <- do.call(rbind, lapply(names(sets), function(n) {
  g <- summ[summ$scores == "vst" & summ$counts == "graft_only" & summ$scored_libraries == "all10" & summ$set == n, ]
  b <- summ[summ$scores == "vst" & summ$counts == "graft_plus_both" & summ$scored_libraries == "all10" & summ$set == n, ]
  if (!nrow(g) || !nrow(b)) return(NULL)
  moves_only_with_both <- (b$rank_min <= 2) && (g$rank_min >= 4 || sign(g$welch_t_true) != sign(b$welch_t_true))
  data.frame(set = n, primary_set = n %in% PRIMARY_SETS, genes = length(sets[[n]]), matched_graft_only = matched$graft_only[n],
             rank_graft_only = g$rank, t_graft_only = g$welch_t_true, rank_graft_plus_both = b$rank, t_graft_plus_both = b$welch_t_true,
             high_in_controls_graft_only = g$high_in_controls, tracks_m_graft_only = g$tracks_m, rho_m_graft_only = g$spearman_with_m_all10,
             high_in_controls_graft_plus_both = b$high_in_controls, tracks_m_graft_plus_both = b$tracks_m,
             moves_only_with_both_bin = moves_only_with_both,
             leak_verdict = ifelse(g$high_in_controls || g$tracks_m || moves_only_with_both, "measures leak (see flags)", "no leak flag"),
             stringsAsFactors = FALSE)
}))
wt(summ, file.path(Cd, "gsva_split_summary.tsv"))
wt(lk, file.path(Cd, "leak_checks.tsv"))
origin$matched_graft_only <- matched$graft_only[origin$set]; origin$matched_graft_plus_both <- matched$graft_plus_both[origin$set]
wt(origin, file.path(GS, "sets_used.tsv"))
# every library's primary-set scores beside its tumour content (graft %, human %, both share); DNA fraction joined later
sc10 <- scores[["vst:graft_only:all10"]]
if (!is.null(sc10)) {
  arm <- ifelse(LIBS %in% PRIM, "primary", ifelse(LIBS %in% RECUR, "recurrent", ifelse(LIBS == "C2B", "culture", ifelse(LIBS %in% CTRL, "control", "contralateral"))))
  bt <- data.frame(library = LIBS, arm = arm, graft_pct = xe[LIBS, "graft_pct"], human_pct = xe[LIBS, "human_pct"],
                   both_share_m = xe[LIBS, "m_both_share"], dna_tumour_fraction = NA,
                   t(sc10[intersect(PRIMARY_SETS, rownames(sc10)), LIBS, drop = FALSE]), check.names = FALSE)
  wt(bt, file.path(Cd, "scores_beside_tumour_content.tsv"))
}
prim <- summ[summ$scores == "vst" & summ$counts == "graft_only" & summ$scored_libraries == "all10" & summ$primary_set, ]
for (i in seq_len(nrow(prim))) kv("C", paste0("gsva_rank_", prim$set[i]), sprintf("%s (Welch t %.2f, recurrent-primary %.3f)", prim$rank[i], prim$welch_t_true[i], prim$diff_recurrent_minus_primary[i]))
a10 <- summ[summ$scores == "vst" & summ$counts == "graft_only" & summ$scored_libraries == "all10", ]
kv("C", "gsva_sets_at_rank1_graft_only", sprintf("%d of %d (expected under the null: %.1f)", sum(a10$rank_min == 1, na.rm = TRUE), nrow(a10), nrow(a10) / 10))

## ------------------------------------------------------------------------------------------------ fgsea on the relabelling fits
load_fits <- function(dir, want_cols_strip = "_graft") {
  # human_de_null layout (as host_verify): fits.tsv manifest + ranks/<label>.tsv or fits/<label>.tsv.gz with gene_name, stat
  mf <- file.path(dir, "fits.tsv")
  if (!file.exists(mf)) return(list(ok = FALSE, why = paste("no manifest", mf)))
  man <- read.delim(mf, stringsAsFactors = FALSE)
  if (!all(c("label", "type", "target") %in% colnames(man))) return(list(ok = FALSE, why = "manifest columns differ"))
  man <- man[man$type == "relabel_3v3", ]
  norm <- function(s) paste(sort(sub(paste0(want_cols_strip, "$"), "", splitv(s))), collapse = ",")
  want <- vapply(SP, function(s) paste(sort(s$B), collapse = ","), "")
  if (nrow(man) != 10) return(list(ok = FALSE, why = sprintf("%d relabel_3v3 rows, not 10", nrow(man))))
  got <- vapply(man$target, norm, "")
  if (!setequal(got, want)) return(list(ok = FALSE, why = "the ten targets are not the ten balanced splits"))
  st <- list()
  for (i in seq_len(nrow(man))) {
    lab <- names(want)[match(got[i], want)]
    f1 <- file.path(dir, "ranks", paste0(man$label[i], ".tsv")); f2 <- file.path(dir, "fits", paste0(man$label[i], ".tsv.gz"))
    f3 <- file.path(dir, paste0(man$label[i], ".tsv.gz"))                     # this job's own fits (step A)
    d <- if (file.exists(f1)) read.delim(f1, stringsAsFactors = FALSE) else if (file.exists(f2)) read.delim(gzfile(f2), stringsAsFactors = FALSE) else
         if (file.exists(f3)) read.delim(gzfile(f3), stringsAsFactors = FALSE) else NULL
    if (is.null(d) || !all(c("gene_name", "stat") %in% colnames(d))) return(list(ok = FALSE, why = paste("no stat table for", man$label[i])))
    d <- d[!is.na(d$stat) & !is.na(d$gene_name) & nzchar(d$gene_name), ]
    d <- d[order(-abs(d$stat)), ]; d <- d[!duplicated(d$gene_name), ]
    # ranks use |NES|, so a fit written as A over B (sign flipped) ranks the same; NES signs are reported as B over A
    # only for this job's own fits (contrast B vs A)
    st[[lab]] <- setNames(d$stat, d$gene_name)
  }
  list(ok = TRUE, stats = st, why = "")
}
fit_sources <- list(
  graft_only = list(ext = P$passB_dir, own = file.path(Cd, "fits_own", "graft_only")),
  graft_plus_both = list(ext = P$passA_dir, own = file.path(Cd, "fits_own", "graft_plus_both")))
BP <- if (as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")) > 1) MulticoreParam(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))) else SerialParam()
fg_all <- list(); fg_rank <- list()
for (src in names(fit_sources)) {
  own <- load_fits(fit_sources[[src]]$own)
  ext <- load_fits(fit_sources[[src]]$ext)
  use <- if (ext$ok) ext else own
  kv("C", paste0("fgsea_fits_", src), if (ext$ok) paste("human_de_null:", fit_sources[[src]]$ext) else paste("this job's own DESeq2 fits (human_de_null not usable:", ext$why, ")"))
  if (ext$ok && own$ok) kv("C", paste0("fits_agreement_true_split_", src), {
    a <- ext$stats[[TL]]; b <- own$stats[[TL]]; g <- intersect(names(a), names(b))
    sprintf("Spearman of Wald statistics, human_de_null vs own, true split: %.3f over %d genes", cor(a[g], b[g], method = "spearman"), length(g)) })
  if (!use$ok) { kv("C", paste0("fgsea_not_run_", src), use$why); next }
  nes <- list()
  for (lab in names(SP)) {
    set.seed(SEED)
    r <- as.data.frame(fgseaSimple(sets, use$stats[[lab]], nperm = as.integer(P$fgsea_nperm), minSize = as.integer(P$fgsea_min_size),
                                   maxSize = as.integer(P$fgsea_max_size), BPPARAM = BP))
    r$leadingEdge <- vapply(r$leadingEdge, paste, "", collapse = ",")
    r$split <- lab; r$true_split <- lab == TL; r$counts <- src
    fg_all[[length(fg_all) + 1]] <- r[, c("counts", "split", "true_split", "pathway", "size", "ES", "NES", "pval", "padj", "leadingEdge")]
    nes[[lab]] <- setNames(r$NES, r$pathway)
  }
  for (n in names(sets)) {
    v <- vapply(names(SP), function(l) if (n %in% names(nes[[l]])) nes[[l]][[n]] else NA_real_, 0)
    rr <- rank_true(v, TL)
    fg_rank[[length(fg_rank) + 1]] <- data.frame(counts = src, set = n, primary_set = n %in% PRIMARY_SETS, NES_true = v[TL],
                                                 rank_min = rr[1], rank_max = rr[2], rank = fmt_rank(rr), p_floor = rr[1] / 10, t(v),
                                                 check.names = FALSE, stringsAsFactors = FALSE)
  }
}
if (length(fg_all)) {
  wt(do.call(rbind, fg_all), file.path(Cd, "fgsea_all_splits.tsv"))
  fr <- do.call(rbind, fg_rank); wt(fr, file.path(Cd, "fgsea_split_ranks.tsv"))
  pf <- fr[fr$primary_set & fr$counts == "graft_only", ]
  for (i in seq_len(nrow(pf))) kv("C", paste0("fgsea_rank_", pf$set[i]), sprintf("%s (NES %.2f)", pf$rank[i], pf$NES_true[i]))
}
writeLines(capture.output(sessionInfo()), file.path(OUT, "inputs", "sessionInfo_stepC.txt"))
cat("step C done\n")
