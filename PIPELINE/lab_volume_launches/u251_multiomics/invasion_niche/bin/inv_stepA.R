#!/usr/bin/env Rscript
# invasion_niche step A. Runs in the DESeq2 1.34.0 image (R 4.1.3). Parts A (Y dosage), B (state of the contralateral
# U251 cells), D (host niche) are base R; the DESeq2 part writes the vst matrices and the ten relabelling fits that
# part C (step C, other image) scores. Each part runs inside its own tryCatch: one failing part does not stop the
# others, and the failure is written to results/STATUS_A.tsv.
#   Rscript inv_stepA.R <params.yaml>
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "inv_common.R"))
P <- read_params(args[1])
OUT <- P$outdir
for (d in c("A_y_dosage", "B_contralateral_state", "C_programmes", "D_host_niche", "inputs"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)
kv_file <- file.path(OUT, "summary_A.tsv"); cat("part\tkey\tvalue\n", file = kv_file)
status <- data.frame(part = character(), status = character(), message = character(), stringsAsFactors = FALSE)
run_part <- function(name, expr) {
  t0 <- Sys.time()
  r <- tryCatch({ force(expr); "ok" }, error = function(e) { message("PART ", name, " FAILED: ", conditionMessage(e)); conditionMessage(e) })
  status[nrow(status) + 1, ] <<- list(name, if (identical(r, "ok")) "ok" else "failed", if (identical(r, "ok")) "" else r)
  cat(sprintf("part %s: %s (%.1f s)\n", name, if (identical(r, "ok")) "ok" else "FAILED", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  wt(status, file.path(OUT, "STATUS_A.tsv"))
}
SEED <- as.integer(P$seed)
LIBS <- splitv(P$libraries); PRIM <- splitv(P$primary); RECUR <- splitv(P$recurrent); TUM <- c(PRIM, RECUR)
CTRL <- splitv(P$controls)
SP <- make_splits(PRIM, RECUR); TL <- true_label(SP)
arm_of <- function(l) ifelse(l %in% PRIM, "primary", ifelse(l %in% RECUR, "recurrent", ifelse(l == "C2B", "culture",
                      ifelse(l %in% CTRL, "control", "contralateral"))))

## ---------------------------------------------------------------------------------------------- shared inputs
hg <- read_gene_table(file.path(P$gene_tables_dir, "gencode_v44_genes.tsv"))
rownames(hg) <- hg$gene_id
G <- read_counts(P$graft_counts)
PUB <- read_counts(P$pub_counts)
gcols <- paste0(LIBS, "_graft")
if (!all(gcols %in% colnames(G$m))) stop("graft-only table lacks: ", paste(setdiff(gcols, colnames(G$m)), collapse = ","))
if (!all(LIBS %in% colnames(PUB$m))) stop("published table lacks: ", paste(setdiff(LIBS, colnames(PUB$m)), collapse = ","))
gm <- G$m[, gcols]; colnames(gm) <- LIBS
pm <- PUB$m[, LIBS]
if (!setequal(rownames(gm), rownames(pm))) stop("graft-only and published tables do not carry the same gene ids")
pm <- pm[rownames(gm), ]
if (!all(rownames(gm) %in% hg$gene_id)) stop("count-table gene ids missing from the GENCODE v44 gene table")
gname <- setNames(hg$gene_name, hg$gene_id)
XE <- read_xeng(P$xengsort_dir, LIBS); rownames(XE) <- XE$lib
XE$m_both_share <- XE$both / (XE$graft + XE$both)
XE$graft_pct <- 100 * XE$graft / XE$total
XE$human_pct <- 100 * (XE$graft + XE$both) / XE$total
sra <- as.numeric(splitv(P$sra_spots))
xe_sra <- isTRUE(all.equal(sort(XE$total), sort(sra)))
wt(XE, file.path(OUT, "inputs", "xengsort_bins.tsv"))
kv("inputs", "xengsort_totals_equal_SRA_spot_counts", xe_sra)
kv("inputs", "genes_in_count_tables", nrow(gm))
cat("xengsort totals equal the SRA spot counts (multiset):", xe_sra, "\n")
if (!xe_sra) warning("xengsort totals do not equal the SRA spot counts; Part A denominators are suspect")

## ============================================================================================== PART A
run_part("A", {
  A <- file.path(OUT, "A_y_dosage")
  yn <- splitv(P$y_genes)
  yr <- hg[hg$chr == "chrY" & hg$gene_name %in% yn & !grepl("_PAR_Y$", hg$gene_id), ]
  miss <- setdiff(yn, yr$gene_name)
  if (length(miss)) stop("Y genes not found on chrY outside the PAR: ", paste(miss, collapse = ","))
  pr <- do.call(rbind, lapply(splitv(P$yx_pairs), function(p) { s <- strsplit(p, "/", fixed = TRUE)[[1]]
    data.frame(y = s[1], x = s[2], stringsAsFactors = FALSE) }))
  xr <- hg[hg$chr == "chrX" & hg$gene_name %in% pr$x, ]
  if (!setequal(unique(xr$gene_name), pr$x)) stop("X gametologs not all found on chrX: ", paste(setdiff(pr$x, xr$gene_name), collapse = ","))
  ids <- rbind(data.frame(role = "Y (X-degenerate, non-PAR)", yr[, c("gene_name", "gene_id", "chr", "start", "end", "gene_type")]),
               data.frame(role = "X gametolog", xr[, c("gene_name", "gene_id", "chr", "start", "end", "gene_type")]))
  ids$in_count_table <- ids$gene_id %in% rownames(gm)
  wt(ids, file.path(A, "gene_ids_from_gtf.tsv"))
  if (!all(ids$in_count_table)) stop("some Y/X gene ids are not in the count tables")
  kv("A", "y_gene_ids", paste(sprintf("%s=%s", yr$gene_name, yr$gene_id), collapse = ";"))
  yid <- yr$gene_id
  srcs <- list(graft_only = gm, graft_plus_both = pm)
  expd <- splitv(P$y_expected_pub, ";")
  expd <- setNames(lapply(expd, function(e) { s <- strsplit(e, "=", fixed = TRUE)[[1]]; as.numeric(strsplit(s[2], "-", fixed = TRUE)[[1]]) }),
                   vapply(expd, function(e) strsplit(e, "=", fixed = TRUE)[[1]][1], ""))
  rows <- list(); per_gene <- list()
  c2b_conv <- (XE["C2B", "graft"] + XE["C2B", "both"]) / XE["C2B", "graft"]
  for (s in names(srcs)) {
    M <- srcs[[s]]
    yt <- colSums(M[yid, , drop = FALSE])
    per_gene[[s]] <- data.frame(source = s, gene_name = gname[yid], gene_id = yid, M[yid, , drop = FALSE], check.names = FALSE)
    cal <- yt[TUM] / XE[TUM, "graft"] * 1e6
    cmed <- median(cal); cmin <- min(cal); cmax <- max(cal)
    kv("A", paste0(s, "_calibration_Y_per_million_graft_pairs_tumours"), paste(sprintf("%s=%.1f", TUM, cal), collapse = ";"))
    kv("A", paste0(s, "_calibration_range"), sprintf("%.1f-%.1f (median %.1f); C2B %.1f", cmin, cmax, cmed, yt["C2B"] / XE["C2B", "graft"] * 1e6))
    for (l in LIBS) {
      y <- unname(yt[l]); yi <- round(y); ci <- pois_ci(yi)
      tot <- XE[l, "total"]
      eq <- y / (cmed * 1e-6)
      r <- data.frame(source = s, library = l, arm = arm_of(l), y_count = y, y_count_rounded = yi,
                      poisson95_lo = ci[1, "lo"], poisson95_hi = ci[1, "hi"],
                      graft_pairs = XE[l, "graft"], both_pairs = XE[l, "both"], total_pairs = tot,
                      y_per_million_graft_pairs = y / XE[l, "graft"] * 1e6,
                      cal_median = cmed, cal_min = cmin, cal_max = cmax,
                      u251_graft_equiv_pairs = eq,
                      u251_pct_of_library = 100 * eq / tot,
                      u251_pct_cal_range_lo = 100 * y / (cmax * 1e-6) / tot,
                      u251_pct_cal_range_hi = 100 * y / (cmin * 1e-6) / tot,
                      u251_pct_lo_poisson_and_cal = 100 * ci[1, "lo"] / (cmax * 1e-6) / tot,
                      u251_pct_hi_poisson_and_cal = 100 * ci[1, "hi"] / (cmin * 1e-6) / tot,
                      u251_pct_of_graft_bin = 100 * eq / XE[l, "graft"],
                      u251_pct_all_pairs_c2b_conversion = 100 * eq * c2b_conv / tot,
                      stringsAsFactors = FALSE)
      r$preliminary_pct_range <- if (s == "graft_plus_both" && l %in% names(expd)) paste(expd[[l]], collapse = "-") else NA
      r$matches_preliminary <- if (s == "graft_plus_both" && l %in% names(expd))
        (r$u251_pct_of_library >= expd[[l]][1] - 1e-9 & r$u251_pct_of_library <= expd[[l]][2] + 1e-9) else NA
      rows[[length(rows) + 1]] <- r
    }
  }
  yd <- do.call(rbind, rows)
  wt(yd, file.path(A, "y_dosage.tsv"))
  wt(do.call(rbind, per_gene), file.path(A, "y_gene_counts.tsv"))
  for (l in splitv(P$y_targets)) for (s in names(srcs)) {
    r <- yd[yd$source == s & yd$library == l, ]
    kv("A", paste0(s, "_", l), sprintf("Y %.1f reads; U251 %.3f%% of library (cal range %.3f-%.3f; with Poisson %.3f-%.3f); %.1f%% of the graft bin",
                                      r$y_count, r$u251_pct_of_library, r$u251_pct_cal_range_lo, r$u251_pct_cal_range_hi,
                                      r$u251_pct_lo_poisson_and_cal, r$u251_pct_hi_poisson_and_cal, r$u251_pct_of_graft_bin))
    if (!is.na(r$matches_preliminary)) kv("A", paste0("preliminary_reproduced_", l), r$matches_preliminary)
  }

  # Y/X gametolog ratios (a clonal marker: Cellosaurus records Y-negative U-251 subclones)
  yx <- list(); rk <- list()
  for (s in names(srcs)) {
    M <- srcs[[s]]
    cnt <- function(nm, chr) { id <- hg$gene_id[hg$gene_name == nm & hg$chr == chr & !grepl("_PAR_Y$", hg$gene_id)]; colSums(M[id, , drop = FALSE]) }
    L2 <- list()
    for (i in seq_len(nrow(pr))) {
      yc <- cnt(pr$y[i], "chrY"); xc <- cnt(pr$x[i], "chrX")
      nm <- paste0(pr$y[i], "/", pr$x[i])
      yx[[length(yx) + 1]] <- data.frame(source = s, library = LIBS, arm = arm_of(LIBS), pair = nm, y_gene = pr$y[i], x_gene = pr$x[i],
                                         y_count = yc[LIBS], x_count = xc[LIBS], ratio = ifelse(xc[LIBS] > 0, yc[LIBS] / xc[LIBS], NA),
                                         log2_ratio_pc0.5 = log2((yc[LIBS] + 0.5) / (xc[LIBS] + 0.5)), stringsAsFactors = FALSE)
      L2[[nm]] <- log2((yc + 0.5) / (xc + 0.5))
    }
    ys <- Reduce(`+`, lapply(seq_len(nrow(pr)), function(i) cnt(pr$y[i], "chrY")))
    xs <- Reduce(`+`, lapply(seq_len(nrow(pr)), function(i) cnt(pr$x[i], "chrX")))
    yx[[length(yx) + 1]] <- data.frame(source = s, library = LIBS, arm = arm_of(LIBS), pair = "pooled_7_pairs", y_gene = "sum", x_gene = "sum",
                                       y_count = ys[LIBS], x_count = xs[LIBS], ratio = ifelse(xs[LIBS] > 0, ys[LIBS] / xs[LIBS], NA),
                                       log2_ratio_pc0.5 = log2((ys[LIBS] + 0.5) / (xs[LIBS] + 0.5)), stringsAsFactors = FALSE)
    L2[["pooled_7_pairs"]] <- log2((ys + 0.5) / (xs + 0.5))
    for (nm in names(L2)) {
      v <- L2[[nm]]; st <- split_stats(v, SP); r <- rank_true(st, TL)
      rk[[length(rk) + 1]] <- data.frame(source = s, pair = nm, mean_log2_primary = mean(v[PRIM]), mean_log2_recurrent = mean(v[RECUR]),
                                         log2_C2B = v["C2B"], diff_recurrent_minus_primary = mean(v[RECUR]) - mean(v[PRIM]),
                                         welch_t_true = st[TL], rank_min = r[1], rank_max = r[2], rank = fmt_rank(r),
                                         p_floor_rank_over_10 = r[1] / 10, t(st), check.names = FALSE, stringsAsFactors = FALSE)
    }
  }
  yx <- do.call(rbind, yx); rk <- do.call(rbind, rk)
  wt(yx, file.path(A, "y_x_ratios.tsv"))
  wt(yx, file.path(OUT, "y_x_ratios.tsv"))            # the copy the variant integration reads
  wt(rk, file.path(A, "y_x_split_ranks.tsv"))
  p <- rk[rk$pair == "pooled_7_pairs", ]
  for (i in seq_len(nrow(p))) kv("A", paste0("yx_pooled_rank_", p$source[i]), sprintf("%s (t %.2f; recurrent-primary %.3f log2)", p$rank[i], p$welch_t_true[i], p$diff_recurrent_minus_primary[i]))
})

## ============================================================================================== PART B
run_part("B", {
  Bd <- file.path(OUT, "B_contralateral_state")
  set.seed(SEED)
  ND <- as.integer(P$thin_draws)
  gt <- setNames(hg$gene_type, hg$gene_id); gc <- setNames(hg$chr, hg$gene_id)
  pc <- rownames(gm)[gt[rownames(gm)] == "protein_coding" & gc[rownames(gm)] != "chrM"]
  # background: what the two near-empty controls put in their graft bin, scaled to N269B's host depth
  s64 <- XE["N269B", "host"] / XE["IL64B", "host"]; s168 <- XE["N269B", "host"] / XE["N168B", "host"]
  bg <- pmax(gm[, "IL64B"] * s64, gm[, "N168B"] * s168)
  flag <- bg > as.numeric(P$background_max_share) * gm[, "N269B"]
  refs <- c("C2B", TUM)
  cpm_on <- function(v, genes) v[genes] / sum(v[genes]) * 1e6
  cr_pc <- sapply(refs, function(r) cpm_on(gm[, r], pc))
  expressed <- pc[rowMeans(cr_pc) >= as.numeric(P$universe_min_mean_cpm)]
  U <- expressed[!flag[expressed]]
  kv("B", "protein_coding_genes", length(pc)); kv("B", "expressed_mean_cpm_ge_threshold", length(expressed))
  kv("B", "background_flagged_among_expressed", sum(flag[expressed])); kv("B", "universe_genes", length(U))
  kv("B", "N269B_share_of_pc_counts_in_flagged_genes", sprintf("%.4f", sum(gm[pc, "N269B"][flag[pc]]) / sum(gm[pc, "N269B"])))
  wt(data.frame(gene_id = pc, gene_name = gname[pc], N269B = gm[pc, "N269B"], IL64B = gm[pc, "IL64B"], N168B = gm[pc, "N168B"],
                background_expected_in_N269B = bg[pc], flagged = flag[pc], expressed = pc %in% expressed, in_universe = pc %in% U),
     file.path(Bd, "background_flags.tsv.gz"))

  # N269B's deduplicated depth (Picard MarkDuplicates on the graft-only run)
  pic <- read_picard(file.path(P$graft_picard_dir, "N269B_graft.markdup.sorted.MarkDuplicates.metrics.txt"))
  dup <- pic$PERCENT_DUPLICATION
  nU <- sum(gm[U, "N269B"]); D <- round(nU * (1 - dup))
  kv("B", "N269B_picard", sprintf("pairs examined %s, pair duplicates %s, PERCENT_DUPLICATION %.4f, unique pairs %s",
                                   fint(pic$READ_PAIRS_EXAMINED), fint(pic$READ_PAIR_DUPLICATES), dup,
                                   fint(pic$READ_PAIRS_EXAMINED - pic$READ_PAIR_DUPLICATES)))
  kv("B", "N269B_universe_counts", round(nU)); kv("B", "thinning_target_D", D)
  cpmU <- function(v) v / sum(v) * 1e6
  X <- sapply(c(refs, "N269B"), function(r) cpmU(gm[U, r]))
  thin <- function(v) { v <- round(v); p <- min(1, D / sum(v)); matrix(rbinom(length(v) * ND, rep(v, ND), p), nrow = length(v)) }
  tI <- thin(gm[U, "IL69B"]); tC <- thin(gm[U, "C2B"])
  others5 <- setdiff(TUM, "IL69B")
  tO <- lapply(setNames(others5, others5), function(k) thin(gm[U, k]))
  cpmM <- function(m) sweep(m, 2, colSums(m), "/") * 1e6
  tIc <- cpmM(tI); tCc <- cpmM(tC); tOc <- lapply(tO, cpmM)
  kv("B", "thin_draws", ND)

  # --- Spearman: N269B and the depth-matched references against every reference profile
  sp_rows <- list()
  for (r in refs) {
    rho_n <- cor(X[, "N269B"], X[, r], method = "spearman")
    rI <- apply(tIc, 2, function(v) cor(v, X[, r], method = "spearman"))
    rC <- apply(tCc, 2, function(v) cor(v, X[, r], method = "spearman"))
    sp_rows[[r]] <- data.frame(reference = r, rho_N269B = rho_n,
                               thinIL69B_median = median(rI), thinIL69B_q025 = quantile(rI, .025), thinIL69B_q975 = quantile(rI, .975),
                               thinC2B_median = median(rC), thinC2B_q025 = quantile(rC, .025), thinC2B_q975 = quantile(rC, .975),
                               note = ifelse(r == "IL69B", "thinIL69B vs IL69B is its own depth ceiling",
                                             ifelse(r == "C2B", "thinC2B vs C2B is its own depth ceiling", "")),
                               stringsAsFactors = FALSE)
  }
  spr <- do.call(rbind, sp_rows); rownames(spr) <- NULL
  wt(spr, file.path(Bd, "spearman.tsv"))

  # --- culture-to-in-vivo axis: genes where the five other tumours all sit on one side of C2B by >= 2-fold
  L <- function(x) log2(x + 1)
  axis_def <- function(tums) {
    Lc <- L(X[, "C2B"]); Lt <- L(rowMeans(X[, tums, drop = FALSE])); d <- Lt - Lc
    sd_ <- sign(L(X[, tums, drop = FALSE]) - Lc)
    same <- apply(sd_, 1, function(s) all(s == s[1]) && s[1] != 0)
    ax <- abs(d) >= as.numeric(P$axis_min_abs_log2) & same
    list(ax = ax, Lc = Lc[ax], d = d[ax])
  }
  tpos <- function(a, v) median((L(v[a$ax]) - a$Lc) / a$d)
  a5 <- axis_def(others5)
  kv("B", "axis_genes_C2B_vs_5_tumours", sum(a5$ax))
  t_n <- tpos(a5, X[, "N269B"]); t_i <- tpos(a5, X[, "IL69B"])
  t_tI <- apply(tIc, 2, function(v) tpos(a5, v)); t_tC <- apply(tCc, 2, function(v) tpos(a5, v))
  bio <- list(); t_tO <- list()
  for (k in others5) {
    ak <- axis_def(setdiff(others5, k))
    bio[[k]] <- data.frame(tumour = k, axis_genes = sum(ak$ax), t_full_depth = tpos(ak, X[, k]),
                           t_thinned_median = NA, t_thinned_q025 = NA, t_thinned_q975 = NA)
    tt <- apply(tOc[[k]], 2, function(v) tpos(ak, v)); t_tO[[k]] <- tt
    bio[[k]][, c("t_thinned_median", "t_thinned_q025", "t_thinned_q975")] <- c(median(tt), quantile(tt, c(.025, .975)))
  }
  bio <- do.call(rbind, bio)
  invivo <- c(t_tI, unlist(t_tO)); invivo_iv <- quantile(invivo, c(.025, .975)); cult_iv <- quantile(t_tC, c(.025, .975))
  in_vivo_ok <- t_n >= invivo_iv[1] && t_n <= invivo_iv[2]
  in_cult <- t_n >= cult_iv[1] && t_n <= cult_iv[2]
  verdict <- if (in_vivo_ok && !in_cult) "in-vivo U251 state, not culture-like" else
             if (in_cult && !in_vivo_ok) "culture-like" else "not resolved (inside both or neither null)"
  axis_tab <- data.frame(profile = c("N269B", "IL69B (full depth)", "thinned IL69B", "thinned other tumours (pooled, own leave-out axes)", "thinned C2B"),
                         t_median_or_value = c(t_n, t_i, median(t_tI), median(unlist(t_tO)), median(t_tC)),
                         q025 = c(NA, NA, quantile(t_tI, .025), quantile(unlist(t_tO), .025), quantile(t_tC, .025)),
                         q975 = c(NA, NA, quantile(t_tI, .975), quantile(unlist(t_tO), .975), quantile(t_tC, .975)),
                         draws = c(1, 1, ND, ND * length(others5), ND))
  wt(axis_tab, file.path(Bd, "axis_position.tsv")); wt(bio, file.path(Bd, "axis_tumour_spread.tsv"))
  kv("B", "axis_t_N269B", sprintf("%.3f", t_n))
  kv("B", "axis_invivo_null_95", sprintf("%.3f-%.3f", invivo_iv[1], invivo_iv[2]))
  kv("B", "axis_culture_null_95", sprintf("%.3f-%.3f", cult_iv[1], cult_iv[2]))
  kv("B", "verdict_rule", "in-vivo if N269B's axis position is inside the 95% interval of thinned in-vivo profiles and outside that of thinned C2B")
  kv("B", "verdict", verdict)

  # --- NNLS on C2B plus the six tumours, leave each donor out. Rows are divided by the square root of the reference
  # mean CPM (Pearson scaling: counting noise has variance close to the mean), so neither the few very abundant genes
  # nor the many noisy low ones dominate the fit.
  m_g <- sqrt(rowMeans(X[, refs]))
  y_n <- X[, "N269B"] / m_g
  nn <- list()
  addnn <- function(target, refset, y, lab, draws = 1) {
    f <- nnls_enum(y, X[, refset, drop = FALSE] / m_g)
    nn[[length(nn) + 1]] <<- data.frame(target = target, fit = lab, t(setNames(f$share[refs], refs)), sum_w = sum(f$w), r2 = f$r2, check.names = FALSE)
  }
  addnn("N269B", refs, y_n, "all 7")
  for (r in refs) addnn("N269B", setdiff(refs, r), y_n, paste0("without ", r))
  for (k in TUM) addnn(k, setdiff(refs, k), X[, k] / m_g, paste0("tumour vs C2B + other 5 (", k, " left out)"))
  nnt <- do.call(rbind, nn); nnt[is.na(nnt)] <- 0
  wt(nnt, file.path(Bd, "nnls_weights.tsv"))
  refs_noI <- setdiff(refs, "IL69B")
  cw <- apply(tIc, 2, function(v) { f <- nnls_enum(v / m_g, X[, refs_noI] / m_g); f$share["C2B"] })
  n_noI <- nnls_enum(y_n, X[, refs_noI] / m_g)
  kv("B", "nnls_N269B_C2B_share_without_IL69B", sprintf("%.3f", n_noI$share["C2B"]))
  kv("B", "nnls_thinIL69B_C2B_share_null_95", sprintf("%.3f-%.3f (median %.3f)", quantile(cw, .025), quantile(cw, .975), median(cw)))
  wt(data.frame(draw = seq_along(cw), c2b_share_thinIL69B = cw, axis_t_thinIL69B = t_tI, axis_t_thinC2B = t_tC),
     file.path(Bd, "null_draws.tsv"))

  # --- index pairs: libraries sharing an i7 or i5 with N269B (index hopping cannot be excluded for them)
  # prep_reads wrote M/reads/index_seq/<lib>.txt: key, value, note; the "modal_index" row is the library's index pair
  idx <- vapply(LIBS, function(l) {
    f <- file.path(P$index_seq_dir, paste0(l, ".txt"))
    if (!file.exists(f)) return(NA_character_)
    x <- readLines(f, warn = FALSE)
    mo <- grep("^modal_index\t", x, value = TRUE)
    if (length(mo)) return(strsplit(mo[1], "\t")[[1]][2])
    m <- regmatches(x, regexpr("[ACGTN]{6,}\\+[ACGTN]{6,}", x))
    if (!length(m)) NA_character_ else m[1]
  }, "")
  if (!is.na(idx["N269B"])) {
    p7 <- sub("\\+.*", "", idx); p5 <- sub(".*\\+", "", idx)
    sh <- LIBS[LIBS != "N269B" & !is.na(idx) & (p7 == p7["N269B"] | p5 == p5["N269B"])]
    kv("B", "libraries_sharing_an_index_with_N269B", paste(sh, collapse = ","))
  } else kv("B", "libraries_sharing_an_index_with_N269B", "index files not found")
  wt(data.frame(library = LIBS, index_pair = idx), file.path(Bd, "index_pairs.tsv"))

  # --- gene by gene: N269B against thinned IL69B and the six-tumour spread
  q <- t(apply(tIc, 1, quantile, c(0, .025, .975, 1)))
  qc <- t(apply(tCc, 1, quantile, c(0, 1)))
  six <- X[, TUM]
  nb <- X[, "N269B"]
  out_thin <- nb < q[, 1] | nb > q[, 4]
  out_six <- nb < apply(six, 1, min) | nb > apply(six, 1, max)
  out_c2b <- nb < qc[, 1] | nb > qc[, 2]
  # calibration of the rule on the null itself: each of 20 draws against the other draws' min and max
  cal_rate <- mean(vapply(1:min(20, ND), function(j) { o <- tIc[, -j, drop = FALSE]
    mean(tIc[, j] < apply(o, 1, min) | tIc[, j] > apply(o, 1, max)) }, 0))
  gl <- data.frame(gene_id = U, gene_name = gname[U], cpm_N269B = nb, cpm_IL69B = X[, "IL69B"], cpm_C2B = X[, "C2B"],
                   thinIL69B_min = q[, 1], thinIL69B_q025 = q[, 2], thinIL69B_q975 = q[, 3], thinIL69B_max = q[, 4],
                   six_min = apply(six, 1, min), six_max = apply(six, 1, max),
                   log2_N269B_over_IL69B = log2((nb + 1) / (X[, "IL69B"] + 1)),
                   outside_thinned_IL69B = out_thin, outside_six_tumour_spread = out_six, counts_as_difference = out_thin & out_six,
                   outside_thinned_C2B = out_c2b, stringsAsFactors = FALSE)
  gl <- gl[order(!gl$counts_as_difference, -abs(gl$log2_N269B_over_IL69B)), ]
  wt(gl, file.path(Bd, "gene_level_N269B_vs_IL69B.tsv.gz"))
  kv("B", "genes_outside_thinned_IL69B", sum(out_thin)); kv("B", "genes_outside_six_tumour_spread", sum(out_six))
  kv("B", "genes_counting_as_difference", sprintf("%d (%d up, %d down)", sum(out_thin & out_six), sum(out_thin & out_six & nb > X[, "IL69B"]), sum(out_thin & out_six & nb < X[, "IL69B"])))
  kv("B", "null_rate_outside_thinned_minmax", sprintf("%.4f (expected %.0f genes if N269B were thinned IL69B)", cal_rate, cal_rate * length(U)))
  kv("B", "genes_outside_thinned_C2B", sum(out_c2b))
})

## ============================================================================================== PART D
run_part("D", {
  Dd <- file.path(OUT, "D_host_niche")
  rg <- read_gene_table(file.path(P$gene_tables_dir, "mRatBN7.2_110_genes.tsv"))
  HC <- read_counts(P$host_counts)
  hm <- HC$m
  hl <- colnames(hm)
  stopifnot(all(TUM %in% hl), all(c("N168B", "N269B", "IL64B") %in% hl))
  hcpm <- sweep(hm, 2, colSums(hm), "/") * 1e6
  rid <- function(nm) rg$gene_id[rg$gene_name == nm & rg$gene_id %in% rownames(hm)]
  panels <- list(myeloid = splitv(P$panel_myeloid), reactive_glia = splitv(P$panel_reactive_glia),
                 fibroblast_collagen = splitv(P$panel_fibroblast_collagen), meninges = splitv(P$panel_meninges),
                 neuron = splitv(P$panel_neuron), white_matter = splitv(P$panel_white_matter))
  getc <- function(nm) { id <- rid(nm); if (!length(id)) return(setNames(rep(NA_real_, length(hl)), hl)); colSums(hcpm[id, , drop = FALSE]) }
  pt <- list(); missing <- character(0)
  for (pn in names(panels)) for (g in panels[[pn]]) {
    id <- rid(g); if (!length(id)) missing <- c(missing, g)
    v <- getc(g)
    pt[[length(pt) + 1]] <- data.frame(panel = pn, gene = g, rat_gene_ids = paste(id, collapse = ","), t(v),
                                       tumour_min = min(v[TUM]), tumour_median = median(v[TUM]), tumour_max = max(v[TUM]),
                                       log2_N269B_over_N168B = log2((v["N269B"] + 0.5) / (v["N168B"] + 0.5)),
                                       log2_N269B_over_IL64B = log2((v["N269B"] + 0.5) / (v["IL64B"] + 0.5)),
                                       log2_N269B_over_tumour_median = log2((v["N269B"] + 0.5) / (median(v[TUM]) + 0.5)),
                                       check.names = FALSE, stringsAsFactors = FALSE)
  }
  pt <- do.call(rbind, pt)
  wt(pt, file.path(Dd, "marker_panels_cpm.tsv"))
  ps <- do.call(rbind, lapply(split(pt, pt$panel), function(d) data.frame(panel = d$panel[1], genes_found = sum(nzchar(d$rat_gene_ids)),
         median_log2_N269B_over_N168B = median(d$log2_N269B_over_N168B, na.rm = TRUE),
         median_log2_N269B_over_IL64B = median(d$log2_N269B_over_IL64B, na.rm = TRUE),
         median_log2_N269B_over_tumour_median = median(d$log2_N269B_over_tumour_median, na.rm = TRUE))))
  ps <- ps[match(names(panels), ps$panel), ]
  wt(ps, file.path(Dd, "marker_panels_summary.tsv"))
  kv("D", "panel_genes_missing_from_rat_gtf", if (length(missing)) paste(missing, collapse = ",") else "none")
  for (i in seq_len(nrow(ps))) kv("D", paste0("panel_", ps$panel[i]), sprintf("N269B/N168B median log2 %.2f; N269B/IL64B %.2f; N269B/tumour median %.2f",
                                 ps$median_log2_N269B_over_N168B[i], ps$median_log2_N269B_over_IL64B[i], ps$median_log2_N269B_over_tumour_median[i]))

  # female host check, outside the human side: rat chrY genes and Xist in the host reads
  ychr <- rg$gene_id[rg$chr %in% c("Y", "chrY") & rg$gene_id %in% rownames(hm)]
  xist <- rid("Xist")
  fh <- data.frame(library = hl, rat_chrY_genes = length(ychr), rat_chrY_counts = if (length(ychr)) colSums(hm[ychr, , drop = FALSE]) else NA,
                   rat_chrY_cpm = if (length(ychr)) colSums(hcpm[ychr, , drop = FALSE]) else NA,
                   Xist_cpm = if (length(xist)) colSums(hcpm[xist, , drop = FALSE]) else NA)
  wt(fh, file.path(Dd, "female_host_check.tsv"))
  kv("D", "rat_chrY_cpm_by_library", paste(sprintf("%s=%.2f", fh$library, fh$rat_chrY_cpm), collapse = ";"))
  kv("D", "rat_Xist_cpm_by_library", if (length(xist)) paste(sprintf("%s=%.1f", fh$library, fh$Xist_cpm), collapse = ";") else "Xist not named in the rat GTF")

  # Mmp13: rank of the true split among the ten host relabellings (host_verify fits)
  mid <- rid("Mmp13"); if (length(mid) != 1) stop("Mmp13 does not resolve to exactly one rat gene id: ", paste(mid, collapse = ","))
  man <- read.delim(P$host_fits_manifest, stringsAsFactors = FALSE)
  man <- man[man$type == "relabel_3v3", ]
  if (nrow(man) != 10 || !setequal(man$label, names(SP))) stop("host fits manifest does not hold the same ten splits")
  mm <- do.call(rbind, lapply(man$label, function(l) {
    f <- read.delim(gzfile(file.path(P$host_fits_dir, paste0(l, ".tsv.gz"))), stringsAsFactors = FALSE)
    r <- f[f$gene_id == mid, ]
    data.frame(label = l, true_split = as.logical(man$true_split[man$label == l]), stat = if (nrow(r)) r$stat else NA,
               lfc_mle = if (nrow(r)) r$lfc_mle else NA, padj = if (nrow(r)) r$padj else NA,
               genes_padj05_in_split = sum(!is.na(f$padj) & f$padj < 0.05))
  }))
  stopifnot(sum(mm$true_split) == 1, mm$label[mm$true_split] == TL)
  st <- setNames(mm$stat, mm$label); rk <- rank_true(st, TL)
  mm$abs_stat_rank <- rank(-abs(mm$stat), ties.method = "min")
  wt(mm, file.path(Dd, "mmp13_relabel_ranks.tsv"))
  kv("D", "mmp13_true_split_rank_by_abs_wald", sprintf("%s (Wald %.2f, lfc %.2f, padj %s)", fmt_rank(rk), st[TL], mm$lfc_mle[mm$true_split], fsig(mm$padj[mm$true_split])))
  file.copy(P$host_r1_mmp13, file.path(Dd, "r1_mmp13_host_verify.tsv"), overwrite = TRUE)
  file.copy(P$host_r2_splits, file.path(Dd, "r2_relabel_splits_host_verify.tsv"), overwrite = TRUE)

  # Mmp13 per unit of stroma and of myeloid cells
  lm13 <- log2(getc("Mmp13") + 0.5)
  st_g <- splitv(P$mmp13_stroma); my_g <- splitv(P$mmp13_myeloid)
  lg <- function(g) log2(getc(g) + 0.5)
  norms <- c(list(raw = rep(0, length(hl))), setNames(lapply(c(st_g, my_g), lg), paste0("per_", c(st_g, my_g))),
             list(per_stroma_geomean = Reduce(`+`, lapply(st_g, lg)) / length(st_g),
                  per_myeloid_geomean = Reduce(`+`, lapply(my_g, lg)) / length(my_g)))
  nr <- do.call(rbind, lapply(names(norms), function(n) {
    v <- setNames(lm13 - norms[[n]], hl); s <- split_stats(v, SP); r <- rank_true(s, TL)
    data.frame(measure = paste0("log2 Mmp13 ", n), mean_primary = mean(v[PRIM]), mean_recurrent = mean(v[RECUR]),
               diff_recurrent_minus_primary = mean(v[RECUR]) - mean(v[PRIM]), welch_t_true = s[TL], rank = fmt_rank(r),
               rank_min = r[1], N168B = v["N168B"], N269B = v["N269B"], IL64B = v["IL64B"], stringsAsFactors = FALSE)
  }))
  wt(nr, file.path(Dd, "mmp13_normalised.tsv"))
  for (i in seq_len(nrow(nr))) kv("D", gsub(" ", "_", nr$measure[i]), sprintf("recurrent-primary %.2f log2; Welch t %.2f; rank %s", nr$diff_recurrent_minus_primary[i], nr$welch_t_true[i], nr$rank[i]))
  cr <- data.frame(marker = c(st_g, my_g), pearson_log_with_Mmp13_six_tumours = vapply(c(st_g, my_g), function(g) cor(lm13[TUM], lg(g)[TUM]), 0))
  wt(cr, file.path(Dd, "mmp13_marker_correlation_tumours.tsv"))

  # the transcript is host: human MMP13 in the human counts beside rat Mmp13 in the host counts
  hid <- hg$gene_id[hg$gene_name == "MMP13" & hg$chr == "chr11"]
  if (length(hid) != 1) stop("human MMP13 does not resolve to one GENCODE id")
  hum <- data.frame(library = LIBS, human_MMP13_graft_only = gm[hid, LIBS], human_MMP13_graft_plus_both = pm[hid, LIBS],
                    human_MMP13_cpm_graft_only = gm[hid, LIBS] / colSums(gm[, LIBS]) * 1e6,
                    rat_Mmp13_counts = ifelse(LIBS %in% hl, hm[mid, match(LIBS, hl)], NA),
                    rat_Mmp13_cpm = ifelse(LIBS %in% hl, hcpm[mid, match(LIBS, hl)], NA))
  wt(hum, file.path(Dd, "mmp13_human_vs_rat.tsv"))
  kv("D", "human_MMP13_graft_only_counts", paste(sprintf("%s=%.0f", hum$library, hum$human_MMP13_graft_only), collapse = ";"))
  kv("D", "rat_Mmp13_counts", paste(sprintf("%s=%s", hum$library, fint(hum$rat_Mmp13_counts)), collapse = ";"))
})

## ============================================================================================== DESeq2 inputs for part C
run_part("C_inputs", {
  Cd <- file.path(OUT, "C_programmes")
  suppressPackageStartupMessages({ library(DESeq2); library(BiocParallel) })
  kv("C", "DESeq2_version_stepA", as.character(packageVersion("DESeq2")))
  ncpu <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  BP <- if (ncpu > 1) MulticoreParam(ncpu) else SerialParam()
  srcs <- list(graft_only = gm, graft_plus_both = pm)
  vst_of <- function(M, libs) {
    cnt <- round(M[, libs]); cnt <- cnt[rowSums(cnt) > 0, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(cnt, data.frame(row.names = libs, lib = libs), ~ 1)
    sft <- tryCatch({ dds <- estimateSizeFactors(dds); "ratio" },
                    error = function(e) { dds <<- estimateSizeFactors(dds, type = "poscounts"); "poscounts" })
    v <- assay(vst(dds, blind = TRUE))
    list(v = v, sf = sizeFactors(dds), type = sft)
  }
  for (s in names(srcs)) for (set in c("all10", "six")) {
    libs <- if (set == "all10") LIBS else TUM
    r <- vst_of(srcs[[s]], libs)
    wt(data.frame(gene_id = rownames(r$v), gene_name = gname[rownames(r$v)], r$v, check.names = FALSE),
       file.path(Cd, "expr", sprintf("vst_%s_%s.tsv.gz", s, set)))
    kv("C", sprintf("vst_%s_%s_size_factors_%s", s, set, r$type), paste(sprintf("%s=%.3f", names(r$sf), r$sf), collapse = ";"))
  }
  for (s in names(srcs)) {
    M <- srcs[[s]][, LIBS]; lc <- log2(sweep(M, 2, colSums(M), "/") * 1e6 + 1); lc <- lc[rowSums(M) > 0, ]
    wt(data.frame(gene_id = rownames(lc), gene_name = gname[rownames(lc)], lc, check.names = FALSE),
       file.path(Cd, "expr", sprintf("logcpm_%s_all10.tsv.gz", s)))
  }
  # the ten relabelling fits on the six tumours: this job's own, used only when human_de_null's are absent
  mc <- as.numeric(P$fit_min_count); ms <- as.integer(P$fit_min_samples)
  for (s in names(srcs)) {
    cnt <- round(srcs[[s]][, TUM]); keep <- rowSums(cnt >= mc) >= ms; cnt <- cnt[keep, ]
    kv("C", paste0("own_fits_universe_", s), sprintf("%d genes (>= %g counts in >= %d of 6)", nrow(cnt), mc, ms))
    man <- list()
    for (l in names(SP)) {
      f <- file.path(Cd, "fits_own", s, paste0(l, ".tsv.gz"))
      grp <- factor(ifelse(TUM %in% SP[[l]]$B, "B", "A"), levels = c("A", "B"))
      if (!file.exists(f)) {
        dds <- DESeqDataSetFromMatrix(cnt, data.frame(row.names = TUM, group = grp), ~ group)
        dds <- DESeq(dds, test = "Wald", fitType = "parametric", minReplicatesForReplace = 99, quiet = TRUE,
                     parallel = ncpu > 1, BPPARAM = BP)
        r <- results(dds, contrast = c("group", "B", "A"), alpha = 0.1)
        wt(data.frame(gene_id = rownames(r), gene_name = gname[rownames(r)], baseMean = r$baseMean, lfc_mle = r$log2FoldChange,
                      stat = r$stat, pvalue = r$pvalue, padj = r$padj), f)
      }
      man[[l]] <- data.frame(label = l, type = "relabel_3v3", reference = paste(SP[[l]]$A, collapse = ","),
                             target = paste(SP[[l]]$B, collapse = ","), true_split = SP[[l]]$true)
    }
    wt(do.call(rbind, man), file.path(Cd, "fits_own", s, "fits.tsv"))
  }
  writeLines(capture.output(sessionInfo()), file.path(OUT, "inputs", "sessionInfo_stepA.txt"))
})
wt(status, file.path(OUT, "STATUS_A.tsv"))
cat("step A done\n"); print(status)
quit(status = if (all(status$status == "ok")) 0 else 2)
