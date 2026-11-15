#!/usr/bin/env Rscript
# integration_rna_prep (SELECTION S8 / EVAL_integration P0): the RNA half of the U251N methylation integration.
# Reads nothing from the methylation run. Base R only, so it runs unchanged in the pinned DESeq2 1.34.0 image (R 4.1).
# README.md in this folder says why every setting is what it is.
#
#   Rscript --vanilla prep_rna_refs.R loci <params.yaml> <samplesheet.csv> <workdir>
#       picks the MGMT / CDKN2A-CDKN2B positive-control genes (same chromosome, C2B TPM 1-5) and writes the loci and
#       their merged exons for run.sbatch's samtools pass
#   Rscript --vanilla prep_rna_refs.R main <params.yaml> <samplesheet.csv> <workdir> <outdir> <published|graft>
#       items 1-6 of the spec; mode graft reruns items 1, 2, 5 and 6 on the rnaseq_graft counts (sensitivity analysis)
#
# Exit codes: 0 complete (failed non-fatal checks are listed in summary.json); 2 a fatal check failed;
#             5 identity STOP (C2B is not nearest to the U-251 family); 1 unexpected error.
options(stringsAsFactors = FALSE, warn = 1, scipen = 12)
A <- commandArgs(trailingOnly = TRUE)
if (length(A) < 4) stop("usage: prep_rna_refs.R loci|main <params.yaml> <samplesheet.csv> <workdir> [<outdir> <mode>]")
STAGE <- A[1]; PF <- A[2]; SS <- A[3]; WD <- A[4]
OUT <- if (length(A) >= 5) A[5] else NA_character_
MODE <- if (length(A) >= 6) A[6] else "published"
HERE <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  dirname(normalizePath(f))
})
source(file.path(HERE, "irp_json.R"))

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
say <- function(...) cat(sprintf("[%s] %s\n", ts(), paste0(...)))
f4 <- function(x, d = 4) formatC(x, digits = d, format = "f")

read_params <- function(f) {
  L <- readLines(f, warn = FALSE)
  L <- L[!grepl("^[[:space:]]*(#|$)", L)]
  L <- sub("[[:space:]]+#.*$", "", L)
  k <- trimws(sub(":.*$", "", L))
  v <- trimws(sub("^[^:]*:", "", L))
  setNames(as.list(v), k)
}
P <- read_params(PF)
pv <- function(k) { v <- P[[k]]; if (is.null(v) || !nzchar(v)) stop("params.yaml: ", k, " missing"); v }
ps <- function(k) trimws(strsplit(pv(k), ",", fixed = TRUE)[[1]])
pn <- function(k) as.numeric(ps(k))

LIBS <- ps("libraries")
HCAL <- pv("human_calibrator")
RCAL <- ps("rat_calibrators")
PRIM <- ps("group_primary")
RECU <- ps("group_recurrent")
HIST_RX <- pv("histone_regex")
group_of <- function(l) ifelse(l %in% PRIM, "primary", ifelse(l %in% RECU, "recurrent",
                        ifelse(l == HCAL, "culture", ifelse(l %in% RCAL, "control_no_tumour", "contralateral"))))

read_gene_tab <- function(path, cols) {
  d <- read.delim(path, check.names = FALSE, quote = "", comment.char = "", colClasses = "character")
  if (!identical(names(d)[1:2], c("gene_id", "gene_name"))) stop(path, ": first columns are not gene_id, gene_name")
  miss <- setdiff(cols, names(d))
  if (length(miss)) stop(path, ": columns missing: ", paste(miss, collapse = ", "))
  if (anyDuplicated(d$gene_id)) stop(path, ": duplicated gene_id")
  m <- sapply(cols, function(cn) as.numeric(d[[cn]]))
  if (!is.matrix(m)) m <- matrix(m, ncol = length(cols), dimnames = list(NULL, cols))
  rownames(m) <- d$gene_id
  if (anyNA(m)) stop(path, ": non-numeric values")
  list(id = d$gene_id, name = d$gene_name, m = m)
}
read_genes <- function() {
  g <- read.delim(file.path(WD, "genes.tsv"), quote = "", comment.char = "",
                  colClasses = c("character", "character", "character", "character", "numeric", "numeric", "character"))
  if (anyDuplicated(g$gene_id)) stop("genes.tsv: duplicated gene_id")
  g
}
merge_iv <- function(s, e) {
  o <- order(s); s <- s[o]; e <- e[o]
  rs <- s[1]; re <- e[1]; out_s <- c(); out_e <- c()
  if (length(s) > 1) for (i in 2:length(s)) {
    if (s[i] <= re + 1) re <- max(re, e[i]) else { out_s <- c(out_s, rs); out_e <- c(out_e, re); rs <- s[i]; re <- e[i] }
  }
  data.frame(start = c(out_s, rs), end = c(out_e, re))
}

# ====================================================================================================== stage: loci
if (STAGE == "loci") {
  genes <- read_genes()
  tp <- read_gene_tab(pv("tpm_published"), LIBS)
  pc <- genes[genes$gene_type == "protein_coding" & !grepl("_PAR_Y$", genes$gene_id), ]
  tg <- ps("locus_targets")
  one <- function(nm) {
    r <- pc[pc$gene_name == nm, ]
    if (nrow(r) != 1) stop(sprintf("%s: %d protein-coding records in the GTF, expected 1", nm, nrow(r)))
    r
  }
  TG <- do.call(rbind, lapply(tg, one))
  if (!all(TG$gene_id %in% tp$id)) stop("a target gene id from the GTF is not in the salmon TPM table")
  c2b <- setNames(tp$m[, HCAL], tp$id)
  lo <- pn("posctrl_tpm_lo"); hi <- pn("posctrl_tpm_hi"); k <- pn("posctrl_n")
  TG$c2b_tpm <- c2b[TG$gene_id]; TG$distance_kb <- 0; TG$role <- "target"
  keep <- c("gene_name", "gene_id", "chr", "start", "end", "strand", "role", "c2b_tpm", "distance_kb")
  L <- TG[, keep]
  for (ch in unique(TG$chr)) {
    g <- TG[TG$chr == ch, ]
    s0 <- min(g$start); e0 <- max(g$end); mid <- (s0 + e0) / 2
    cand <- pc[pc$chr == ch & !(pc$gene_id %in% TG$gene_id) & !grepl(HIST_RX, pc$gene_name) &
               (pc$end < s0 | pc$start > e0), ]
    cand$c2b_tpm <- c2b[cand$gene_id]
    cand <- cand[!is.na(cand$c2b_tpm) & cand$c2b_tpm >= lo & cand$c2b_tpm <= hi, ]
    cand$distance_kb <- abs((cand$start + cand$end) / 2 - mid) / 1000
    cand <- cand[order(cand$distance_kb), ]
    if (nrow(cand) < k) stop(sprintf("%s: only %d candidate controls with C2B TPM %s-%s", ch, nrow(cand), lo, hi))
    cand <- cand[seq_len(k), ]
    cand$role <- paste0("control_for_", paste(g$gene_name, collapse = "_"))
    L <- rbind(L, cand[, keep])
  }
  ex <- read.delim(file.path(WD, "exons_chr9_chr10.tsv"), quote = "", comment.char = "",
                   colClasses = c("character", "character", "numeric", "numeric"))
  L$exon_intervals <- NA; L$exonic_bp <- NA
  for (i in seq_len(nrow(L))) {
    e <- ex[ex$gene_id == L$gene_id[i], ]
    if (nrow(e) == 0) stop(L$gene_name[i], ": no exons in the chr9/chr10 extract")
    m <- merge_iv(e$start, e$end)
    write.table(m, file.path(WD, paste0("exons_", L$gene_name[i], ".txt")), sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE)
    L$exon_intervals[i] <- nrow(m); L$exonic_bp[i] <- sum(m$end - m$start + 1)
  }
  write.table(L, file.path(WD, "loci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  say("loci: ", paste(sprintf("%s (%s, %s, C2B TPM %s)", L$gene_name, L$role, L$chr, f4(L$c2b_tpm, 2)), collapse = "; "))
  quit(save = "no", status = 0)
}
if (STAGE != "main") stop("unknown stage ", STAGE)

# ====================================================================================================== stage: main
if (!MODE %in% c("published", "graft")) stop("mode must be published or graft")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
unlink(file.path(OUT, c("summary.json", "REPORT.md")))
REFS <- pv("refs_dir")
CHECKS <- list(); FAILED <- character(0); SUMMARY <- list(); OUTFILES <- character(0); INPUTS <- character(0)
reg <- function(label, path) INPUTS[[label]] <<- path

wt <- function(d, f) {
  tmp <- file.path(OUT, paste0(f, ".part"))
  write.table(d, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  file.rename(tmp, file.path(OUT, f))
  OUTFILES <<- unique(c(OUTFILES, f))
  say("wrote ", f, " (", nrow(d), " rows)")
}
md_table <- function(d, digits = 4) {
  if (is.null(d) || nrow(d) == 0) return("(none)")
  cv <- lapply(d, function(x) if (is.numeric(x)) ifelse(is.na(x), "NA", trimws(formatC(x, digits = digits, format = "fg"))) else as.character(x))
  rows <- do.call(paste, c(cv, sep = " | "))
  c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
    paste0("| ", rows, " |"))
}
input_md5 <- function() {
  out <- list()
  for (k in names(INPUTS)) {
    p <- INPUTS[[k]]
    out[[k]] <- list(path = p, md5 = if (file.exists(p)) unname(tools::md5sum(p)) else NA)
  }
  bm <- file.path(WD, "bam_md5.txt")
  if (MODE == "published" && file.exists(bm)) {
    b <- read.table(bm, header = FALSE, col.names = c("md5", "path"))
    for (i in seq_len(nrow(b))) out[[paste0("bam:", basename(b$path[i]))]] <- list(path = b$path[i], md5 = b$md5[i])
  }
  out
}
env_info <- function() {
  e <- list(R = R.version.string, node = unname(Sys.info()[["nodename"]]),
            base_packages = paste(sessionInfo()$basePkgs, collapse = ","))
  f <- file.path(WD, "env.tsv")
  if (file.exists(f)) {
    d <- read.delim(f, header = FALSE, col.names = c("key", "value"), quote = "")
    for (i in seq_len(nrow(d))) e[[d$key[i]]] <- d$value[i]
  }
  e
}
REPORT_EXTRA <- character(0)
write_report <- function(S) {
  L <- c(sprintf("# integration_rna_prep: %s counts", MODE), "",
         sprintf("Status **%s**, finished %s. Failed checks: %s.", S$status, S$finished,
                 if (length(FAILED)) paste(FAILED, collapse = ", ") else "none"), "", "## Checks", "")
  for (k in names(CHECKS)) L <- c(L, sprintf("- %s **%s** (%s): %s", if (CHECKS[[k]]$pass) "PASS" else "FAIL", k,
                                            if (CHECKS[[k]]$hard) "fatal" else "flag", CHECKS[[k]]$detail))
  writeLines(c(L, "", REPORT_EXTRA), file.path(OUT, "REPORT.md"))
}
finish <- function(status, code = 0) {
  S <- list(status = status, mode = MODE, finished = ts(), failed_checks = as.list(FAILED), checks = CHECKS,
            results = SUMMARY, outputs = as.list(OUTFILES), inputs = input_md5(), environment = env_info())
  write_json_file(S, file.path(OUT, "summary.json"))
  tryCatch(write_report(S), error = function(e) say("REPORT.md not written: ", conditionMessage(e)))
  say("status: ", status)
  quit(save = "no", status = code)
}
rec <- function(id, pass, hard, detail, data = NULL) {
  pass <- isTRUE(pass)
  CHECKS[[id]] <<- list(pass = pass, hard = hard, detail = detail, data = data)
  say(if (pass) "PASS " else if (hard) "FAIL (fatal) " else "FAIL (flag) ", id, ": ", detail)
  if (!pass) FAILED <<- c(FAILED, id)
  if (!pass && hard) finish(paste0("failed:", id), 2)
}

main <- function() {
  reg("params.yaml", PF); reg("samplesheet.csv", SS)
  for (f in c("prep_rna_refs.R", "irp_json.R")) reg(f, file.path(HERE, f))

  # ------------------------------------------------------------------------------------------ GTF identity
  genes <- read_genes()
  hdr <- readLines(file.path(WD, "gtf_header.txt"), warn = FALSE)
  reg("gtf", pv("gtf")); reg("gtf_MD5SUMS", pv("gtf_md5sums"))
  rec("gtf_is_gencode_v44", any(grepl(pv("gtf_version_text"), hdr, fixed = TRUE)), TRUE,
      paste0("GTF header: ", paste(grep("description", hdr, value = TRUE), collapse = " ")))
  ms <- read.table(pv("gtf_md5sums"), header = FALSE, col.names = c("md5", "file"))
  want <- ms$md5[ms$file == pv("gtf_md5sums_name")]
  got <- unname(tools::md5sum(pv("gtf")))
  rec("gtf_md5_matches_gencode_list", length(want) == 1 && identical(got, want), TRUE,
      sprintf("md5 %s; GENCODE MD5SUMS lists %s for %s", got, paste(want, collapse = ","), pv("gtf_md5sums_name")))

  # ------------------------------------------------------------------------------------------ item 1: human share
  ss <- read.csv(SS, check.names = FALSE, colClasses = "character")
  rec("samplesheet_10_libraries", identical(ss$sample, LIBS) && length(LIBS) == 10, TRUE,
      paste0("samplesheet: ", paste(ss$sample, collapse = ",")))
  parse_xs <- function(f, lib) {
    L <- readLines(f, warn = FALSE)
    i <- grep("^prefix[[:space:]]", L)
    if (length(i) != 1) stop(f, ": ", length(i), " 'prefix' header lines")
    h <- strsplit(trimws(L[i]), "[[:space:]]+")[[1]]
    v <- strsplit(trimws(L[i + 1]), "[[:space:]]+")[[1]]
    if (v[1] != lib) stop(f, ": counts row is for ", v[1], ", expected ", lib)
    x <- setNames(as.numeric(v[-1]), h[-1])
    x[c("host", "graft", "ambiguous", "both", "neither")]
  }
  B <- t(sapply(LIBS, function(l) { f <- ss$xengsort_log[ss$sample == l]; reg(paste0("xengsort:", l), f); parse_xs(f, l) }))
  rec("xengsort_bins_parsed", !anyNA(B) && all(B >= 0), TRUE, "five bins parsed for every library")
  tot <- rowSums(B)
  rec("xengsort_total_pairs", sum(tot) == pn("expected_total_pairs"), TRUE,
      sprintf("sum of the five bins over the ten libraries = %.0f (expected %.0f)", sum(tot), pn("expected_total_pairs")))
  reg("fastqc_raw", pv("fastqc_raw"))
  fq <- read.delim(pv("fastqc_raw"), check.names = FALSE, quote = "", colClasses = "character")
  fqtot <- function(l, r) { x <- fq[["Total Sequences"]][fq$Sample == paste0(l, "_", r)]; if (length(x) != 1) NA else as.numeric(x) }
  hf <- B[, "graft"] + B[, "both"]
  fqt <- cbind(R1 = sapply(LIBS, fqtot, 1), R2 = sapply(LIBS, fqtot, 2))
  bad <- LIBS[is.na(fqt[, 1]) | is.na(fqt[, 2]) | fqt[, 1] != hf | fqt[, 2] != hf]
  rec("graft_plus_both_equals_fastqc_total", length(bad) == 0, TRUE,
      if (length(bad)) paste0("mismatch: ", paste(bad, collapse = ",")) else
        sprintf("graft + both = FastQC Total Sequences (R1 and R2) for all ten; sum %.0f", sum(hf)))
  fr <- B / tot
  hum <- fr[HCAL, ]; ratm <- colMeans(fr[RCAL, , drop = FALSE])
  near <- function(x, s) abs(x - s) <= 0.5 * 10^-(nchar(sub("^[^.]*\\.", "", format(s, scientific = FALSE))))
  stated <- list(c(hum["graft"], pn("stated_human_graft")), c(hum["both"], pn("stated_human_both")),
                 c(hum["host"], pn("stated_human_host")))
  rh <- pn("stated_rat_host"); rg <- pn("stated_rat_graft"); rb <- pn("stated_rat_both")
  for (i in seq_along(RCAL)) stated <- c(stated, list(c(fr[RCAL[i], "host"], rh[i]), c(fr[RCAL[i], "graft"], rg[i]),
                                                      c(fr[RCAL[i], "both"], rb[i])))
  ok <- all(sapply(stated, function(z) near(z[1], z[2])))
  rec("calibrator_rates_match_stated", ok, TRUE,
      sprintf("C2B graft %s both %s host %s; %s host %s/%s graft %s/%s both %s/%s", f4(hum["graft"]), f4(hum["both"]),
              f4(hum["host"]), paste(RCAL, collapse = "/"), f4(fr[RCAL[1], "host"]), f4(fr[RCAL[2], "host"]),
              f4(fr[RCAL[1], "graft"]), f4(fr[RCAL[2], "graft"]), f4(fr[RCAL[1], "both"]), f4(fr[RCAL[2], "both"])))
  h_graft <- (fr[, "graft"] - ratm["graft"]) / (hum["graft"] - ratm["graft"])
  h_host <- (ratm["host"] - fr[, "host"]) / (ratm["host"] - hum["host"])
  dh <- h_host - h_graft
  tol <- pn("h_agree_tol")
  bad <- LIBS[abs(dh) > tol]
  # Not fatal: the preview already puts IL66B at 0.031. A disagreement is a finding about the two-rate model, and both
  # estimates are carried forward as an interval; nothing downstream picks the one that suits it.
  rec("h_graft_h_host_agree_within_tol", length(bad) == 0, FALSE,
      if (length(bad)) sprintf("|h_host - h_graft| > %s in %s (%s)", tol, paste(bad, collapse = ","),
                               paste(f4(dh[bad]), collapse = ",")) else sprintf("max |h_host - h_graft| = %s", f4(max(abs(dh)))),
      data = as.list(setNames(round(dh, 5), LIBS)))
  both_pred <- h_host * hum["both"] + (1 - h_host) * ratm["both"]

  # ------------------------------------------------------------------------------------------ counts for this mode
  sfx <- if (MODE == "graft") pv("graft_suffix") else ""
  cpath <- if (MODE == "graft") pv("counts_graft") else pv("counts_published")
  tpath <- if (MODE == "graft") pv("tpm_graft") else pv("tpm_published")
  reg("counts", cpath); reg("tpm", tpath)
  cnt <- read_gene_tab(cpath, paste0(LIBS, sfx)); colnames(cnt$m) <- LIBS
  tpm <- read_gene_tab(tpath, paste0(LIBS, sfx)); colnames(tpm$m) <- LIBS
  rec("counts_tpm_same_genes", identical(cnt$id, tpm$id), TRUE, sprintf("%d genes in counts and TPM, same order", length(cnt$id)))
  miss <- setdiff(cnt$id, genes$gene_id)
  rec("count_genes_in_gtf", length(miss) == 0, TRUE, sprintf("%d of %d count-table gene ids missing from the GTF", length(miss), length(cnt$id)))
  gi <- genes[match(cnt$id, genes$gene_id), ]

  # chrY estimate of the human share: the rats are female and U251 is male, so chrY reads come only from U251
  cpm_c2b <- cnt$m[, HCAL] / sum(cnt$m[, HCAL]) * 1e6
  ysel <- cnt$id[gi$chr == "chrY" & gi$gene_type == "protein_coding" & !grepl("_PAR_Y$", cnt$id) & cpm_c2b > 10]
  ysum <- colSums(cnt$m[ysel, , drop = FALSE])
  h_Y <- (ysum / tot) / (ysum[HCAL] / tot[HCAL])
  H <- data.frame(library = LIBS, group = group_of(LIBS), total_pairs = tot, host = B[, "host"], graft = B[, "graft"],
                  ambiguous = B[, "ambiguous"], both = B[, "both"], neither = B[, "neither"], human_file_pairs = hf,
                  fastqc_total_R1 = fqt[, 1], raw_graft_frac = round(fr[, "graft"], 6), host_frac = round(fr[, "host"], 6),
                  both_frac = round(fr[, "both"], 6), h_graft = round(h_graft, 5), h_host = round(h_host, 5),
                  h_diff = round(dh, 5), h_agree = abs(dh) <= tol, h_lo = round(pmin(h_graft, h_host), 5),
                  h_hi = round(pmax(h_graft, h_host), 5), both_frac_predicted = round(both_pred, 5),
                  both_obs_minus_pred = round(fr[, "both"] - both_pred, 5), chrY_reads = ysum, h_chrY = round(h_Y, 4))
  wt(H, "h_share.tsv")
  wt(data.frame(bin = c("host", "graft", "ambiguous", "both", "neither"), human_rate_C2B = round(hum, 6),
                rat_rate_mean = round(ratm, 6), rat_rate_IL64B = round(fr[RCAL[1], ], 6), rat_rate_N168B = round(fr[RCAL[2], ], 6)),
     "calibrator_rates.tsv")
  SUMMARY$human_share <<- list(rates_human = as.list(round(hum, 6)), rates_rat_mean = as.list(round(ratm, 6)),
                              chrY_genes = length(ysel), note_both_bin = "both_obs_minus_pred: the both bin is not used by either estimate, so it is a closure check on the two-rate model",
                              by_library = H[, c("library", "group", "raw_graft_frac", "h_graft", "h_host", "h_diff",
                                                 "both_obs_minus_pred", "h_chrY")])

  # ------------------------------------------------------------------------------------------ item 2: rat leak per gene
  Fc <- if (MODE == "graft") B[, "graft"] else hf
  m_rate <- if (MODE == "graft") ratm[["graft"]] else ratm[["graft"]] + ratm[["both"]]
  pc_each <- sapply(RCAL, function(c) cnt$m[, c] / Fc[c])    # per-pair yield of each control alone
  p_ctrl <- rowMeans(pc_each)
  R_model <- (1 - h_host) * tot * m_rate
  R_hg <- (1 - h_graft) * tot * m_rate
  R_resid <- if (MODE == "published") pmax(0, B[, "both"] - h_host * tot * hum[["both"]]) + (1 - h_host) * tot * ratm[["graft"]] else rep(NA, length(LIBS))
  names(R_resid) <- LIBS
  obs <- cnt$m
  leak <- function(R, p = p_ctrl) { e <- outer(p, R); f <- e / obs; f[obs <= 0] <- NA; list(e = e, f = f) }
  L0 <- leak(R_model); Lg <- leak(R_hg)
  Lr <- if (MODE == "published") leak(R_resid) else NULL
  L1 <- leak(R_model, pc_each[, 1]); L2 <- leak(R_model, pc_each[, 2])   # each control's profile alone: the spread
  thr <- pn("rat_leak_mark")
  mark <- !is.na(L0$f) & L0$f > thr
  mark_hg <- !is.na(Lg$f) & Lg$f > thr
  mark_rs <- if (!is.null(Lr)) (!is.na(Lr$f) & Lr$f > thr) else NULL
  mark_cx <- (!is.na(L1$f) & L1$f > thr) | (!is.na(L2$f) & L2$f > thr)
  expressed <- tpm$m >= pn("tpm_min")
  minc <- pn("leak_min_count_selfcheck")
  selfc <- sapply(RCAL, function(c) median(L0$f[obs[, c] >= minc, c], na.rm = TRUE))
  selft <- sapply(RCAL, function(c) sum(L0$e[, c]) / sum(obs[, c]))
  lr12 <- log2(pc_each[, 1] / pc_each[, 2])
  lr12 <- lr12[is.finite(lr12) & obs[, RCAL[1]] >= minc & obs[, RCAL[2]] >= minc]
  rec("leak_model_reproduces_controls", all(selfc > 0.8 & selfc < 1.25), FALSE,
      sprintf("calibrators run through their own model (expect about 1): median per gene (>= %d counts) %s; total %s. Per-gene log2(%s/%s) yield: median %s, IQR %s to %s over %d genes",
              minc, paste(sprintf("%s %s", RCAL, f4(selfc, 3)), collapse = ", "), paste(sprintf("%s %s", RCAL, f4(selft, 3)), collapse = ", "),
              RCAL[1], RCAL[2], f4(median(lr12), 3), f4(quantile(lr12, 0.25), 3), f4(quantile(lr12, 0.75), 3), length(lr12)))
  keep <- rowSums(obs > 0) > 0
  G <- data.frame(gene_id = cnt$id, gene_name = cnt$name, gene_type = gi$gene_type,
                  p_ctrl_per_pair = signif(p_ctrl, 5))[keep, ]
  for (l in LIBS) G[[paste0("frac_", l)]] <- signif(L0$f[keep, l], 4)
  for (l in LIBS) G[[paste0("mark_", l)]] <- as.integer(mark[keep, l])
  wt(G, "rat_leak_genes.tsv")
  LB <- data.frame(library = LIBS, group = group_of(LIBS), h_host = round(h_host, 5), file_pairs = Fc,
                   rat_pairs_in_file_model = round(R_model), rat_pairs_in_file_hgraft = round(R_hg),
                   rat_pairs_in_file_bothresid = round(R_resid),
                   rat_share_of_file_model = round(R_model / Fc, 5),
                   rat_share_of_counts_model = round(colSums(L0$e) / colSums(obs), 5),
                   rat_share_of_counts_hgraft = round(colSums(Lg$e) / colSums(obs), 5),
                   rat_share_of_counts_bothresid = if (!is.null(Lr)) round(colSums(Lr$e) / colSums(obs), 5) else NA,
                   rat_share_of_counts_ctrl1_profile = round(colSums(L1$e) / colSums(obs), 5),
                   rat_share_of_counts_ctrl2_profile = round(colSums(L2$e) / colSums(obs), 5),
                   genes_tpm_ge_min = colSums(expressed), genes_tpm_ge_min_marked = colSums(expressed & mark),
                   marks_changed_if_hgraft = colSums(expressed & (mark != mark_hg)),
                   marks_changed_if_bothresid = if (!is.null(mark_rs)) colSums(expressed & (mark != mark_rs)) else NA,
                   marks_added_if_either_ctrl = colSums(expressed & mark_cx & !mark),
                   median_frac_selfcheck = ifelse(LIBS %in% RCAL, round(selfc[LIBS], 4), NA),
                   total_frac_selfcheck = ifelse(LIBS %in% RCAL, round(selft[LIBS], 4), NA))
  names(LB)[names(LB) == "rat_share_of_counts_ctrl1_profile"] <- paste0("rat_share_of_counts_", RCAL[1], "_profile")
  names(LB)[names(LB) == "rat_share_of_counts_ctrl2_profile"] <- paste0("rat_share_of_counts_", RCAL[2], "_profile")
  wt(LB, "rat_leak_by_library.tsv")
  SUMMARY$rat_leak <<- list(
    formula = if (MODE == "graft") "expected_rat(g,t) = (1 - h_host,t) x total_t x g_R x p_ctrl(g); p_ctrl = mean over IL64B, N168B of graft-only count / graft pairs"
              else "expected_rat(g,t) = (1 - h_host,t) x total_t x (g_R + b_R) x p_ctrl(g); p_ctrl = mean over IL64B, N168B of human-file count / human-file pairs",
    mark_threshold = thr, prior_rat_share_of_counts = "EVAL_transcripts: primary 1.9-5.1 %, recurrent 5.1-8.0 % (graft+both counts)",
    by_library = LB[, c("library", "group", "rat_share_of_counts_model", "rat_share_of_counts_bothresid",
                        paste0("rat_share_of_counts_", RCAL, "_profile"), "genes_tpm_ge_min", "genes_tpm_ge_min_marked",
                        "marks_added_if_either_ctrl")])

  # ------------------------------------------------------------------------------------------ item 3: loci (published)
  if (MODE == "published") {
    loci <- read.delim(file.path(WD, "loci.tsv"), quote = "", comment.char = "")
    lc <- read.delim(file.path(WD, "locus_counts.tsv"), header = FALSE, quote = "",
                     col.names = c("library", "gene_name", "frag_all", "frag_exonic_all", "frag_intronic_all", "frag_dedup",
                                   "frag_exonic_dedup", "frag_intronic_dedup", "sense_exonic_dedup", "antisense_exonic_dedup",
                                   "sense_all", "antisense_all", "reads_all", "reads_dedup"))
    rec("locus_counts_complete", nrow(lc) == nrow(loci) * length(LIBS) && all(LIBS %in% lc$library), TRUE,
        sprintf("%d locus x library rows (expected %d)", nrow(lc), nrow(loci) * length(LIBS)))
    nmh <- read.delim(file.path(WD, "locus_nm.tsv"), header = FALSE, quote = "", colClasses = c("character", "character", "character", "character", "numeric"),
                      col.names = c("library", "gene_name", "duplicate", "NM", "reads"))
    star_u <- sapply(LIBS, function(l) {
      f <- ss$star_log[ss$sample == l]; reg(paste0("star_log:", l), f)
      x <- grep("Uniquely mapped reads number", readLines(f), value = TRUE)
      as.numeric(trimws(sub(".*\\|", "", x)))
    })
    pic <- t(sapply(LIBS, function(l) {
      f <- ss$picard_metrics[ss$sample == l]; reg(paste0("picard:", l), f)
      x <- readLines(f); i <- grep("^## METRICS CLASS", x)
      h <- strsplit(x[i + 1], "\t")[[1]]; v <- strsplit(x[i + 2], "\t")[[1]]
      c(pairs = as.numeric(v[h == "READ_PAIRS_EXAMINED"]), dup_pairs = as.numeric(v[h == "READ_PAIR_DUPLICATES"]))
    }))
    dedup_pairs <- pic[, "pairs"] - pic[, "dup_pairs"]
    lvl <- pn("poisson_ub_level")
    ub <- function(k) qchisq(lvl, 2 * (round(k) + 1)) / 2        # one-sided Poisson; k = 0 gives 2.996 (rule of three)
    csum <- colSums(obs)
    lc$group <- group_of(lc$library)
    lc <- merge(lc, loci[, c("gene_name", "gene_id", "chr", "start", "end", "strand", "role", "exonic_bp")], by = "gene_name")
    lc$salmon_count <- mapply(function(g, l) obs[g, l], lc$gene_id, lc$library)
    lc$salmon_tpm <- mapply(function(g, l) tpm$m[g, l], lc$gene_id, lc$library)
    lc$salmon_cpm <- lc$salmon_count / csum[lc$library] * 1e6
    lc$salmon_cpm_ub95 <- ub(lc$salmon_count) / csum[lc$library] * 1e6
    lc$bam_sense_exonic_dedup_per_M_dedup_pairs <- lc$sense_exonic_dedup / dedup_pairs[lc$library] * 1e6
    lc$bam_cpm_ub95_dedup <- ub(lc$sense_exonic_dedup) / dedup_pairs[lc$library] * 1e6
    lc$rat_leak_frac <- mapply(function(g, l) L0$f[g, l], lc$gene_id, lc$library)
    nd <- nmh[nmh$duplicate == "nondup", ]
    lc$nondup_reads_NM0 <- mapply(function(g, l) sum(nd$reads[nd$library == l & nd$gene_name == g & nd$NM == "0"]), lc$gene_name, lc$library)
    lc$nondup_reads_NMge3 <- mapply(function(g, l) sum(nd$reads[nd$library == l & nd$gene_name == g & nd$NM %in% c("3", "4", "5", "6+")]), lc$gene_name, lc$library)
    lc <- lc[order(match(lc$gene_name, loci$gene_name), match(lc$library, LIBS)), ]
    wt(lc, "mgmt_cdkn2a_loci.tsv")
    wt(nmh[order(nmh$gene_name, match(nmh$library, LIBS), nmh$duplicate, nmh$NM), ], "locus_nm_hist.tsv")
    wt(loci, "locus_definitions.tsv")
    wt(data.frame(library = LIBS, salmon_assigned_counts = round(csum), star_unique_pairs = star_u,
                  picard_pairs = pic[, "pairs"], picard_dup_pairs = pic[, "dup_pairs"], dedup_pairs = dedup_pairs,
                  dup_fraction = round(pic[, "dup_pairs"] / pic[, "pairs"], 4)), "library_denominators.tsv")
    ctl <- lc[grepl("^control", lc$role) & lc$library == HCAL, ]
    rec("positive_controls_detected_in_C2B", nrow(ctl) > 0 && all(ctl$frag_exonic_dedup >= 10), FALSE,
        paste(sprintf("%s %d exonic dedup fragments (salmon %s)", ctl$gene_name, ctl$frag_exonic_dedup, f4(ctl$salmon_count, 1)), collapse = "; "))
    mg <- lc[lc$gene_name == "MGMT", ]
    SUMMARY$loci <<- list(
      mgmt_salmon_counts = as.list(setNames(round(mg$salmon_count, 3), mg$library)),
      mgmt_bam_sense_exonic_dedup = as.list(setNames(mg$sense_exonic_dedup, mg$library)),
      mgmt_bam_intronic_dedup = as.list(setNames(mg$frag_intronic_dedup, mg$library)),
      mgmt_salmon_cpm_ub95 = as.list(setNames(signif(mg$salmon_cpm_ub95, 3), mg$library)),
      preview = "0 MGMT reads in C2B and in 5 of 6 tumours, 2 in IL69B (salmon, EVAL_integration section 5)",
      controls = unique(loci$gene_name[grepl("^control", loci$role)]))
    REPORT_EXTRA <<- c(REPORT_EXTRA, "## MGMT, CDKN2A, CDKN2B and controls (BAM fragments, MAPQ 255)", "",
                       md_table(lc[, c("gene_name", "role", "library", "salmon_count", "frag_exonic_dedup", "sense_exonic_dedup",
                                       "frag_intronic_dedup", "frag_exonic_all", "salmon_cpm_ub95")], 3), "")
  }

  # ------------------------------------------------------------------------------------------ item 4: references
  st <- read.delim(file.path(REFS, "refs_status.tsv"), quote = "", colClasses = "character")
  for (i in seq_len(nrow(st))) reg(paste0("ref:", st$name[i]), file.path(REFS, st$name[i]))
  rec("references_present_md5_ok", all(st$status == "ok"), FALSE,
      paste(sprintf("%s %s", st$name, st$status), collapse = "; "), data = st)
  need <- c("OmicsExpressionProteinCodingGenesTPMLogp1.csv", "Model.csv", "OmicsDefaultModelProfiles.csv", "cytoBand.hg38.txt.gz")
  rec("references_needed_here_ok", all(need %in% st$name[st$status == "ok"]), TRUE, paste(need, collapse = ", "))
  SUMMARY$references <<- st

  # ------------------------------------------------------------------------------------------ DepMap expression
  say("reading DepMap expression")
  ex_path <- file.path(REFS, "OmicsExpressionProteinCodingGenesTPMLogp1.csv")
  # scan() into one character vector is about 3x faster than read.csv on 19,194 columns (1.8 against 5.8 min locally)
  con <- file(ex_path, "r")
  h1 <- gsub('^"|"$', "", strsplit(readLines(con, n = 1), ",", fixed = TRUE)[[1]])
  v <- scan(con, what = "", sep = ",", quiet = TRUE, quote = "\"", na.strings = character(0))
  close(con)
  nc <- length(h1)
  rec("depmap_matrix_rectangular", length(v) %% nc == 0, TRUE, sprintf("%d fields, %d columns", length(v), nc))
  first <- seq(1, length(v), by = nc)
  ids <- v[first]
  X <- matrix(suppressWarnings(as.numeric(v[-first])), nrow = length(ids), byrow = TRUE, dimnames = list(ids, h1[-1]))
  rm(v, first); invisible(gc())
  rec("depmap_rows_are_models", all(grepl("^ACH-", ids)) && !anyDuplicated(ids) && !anyNA(X), TRUE,
      sprintf("%d rows (ModelID, unique), %d gene columns, %d non-numeric values", length(ids), ncol(X), sum(is.na(X))))
  dsym <- sub(" \\([0-9]+\\)$", "", colnames(X))
  dupd <- dsym %in% dsym[duplicated(dsym)]
  X <- X[, !dupd]; colnames(X) <- dsym[!dupd]
  med <- apply(X, 2, median, na.rm = TRUE)
  model <- read.csv(file.path(REFS, "Model.csv"), check.names = FALSE, colClasses = "character", encoding = "UTF-8")
  fam <- ps("depmap_family"); famn <- ps("depmap_family_names")
  fam_names <- model$StrippedCellLineName[match(fam, model$ModelID)]
  rec("depmap_family_ids_match_names", identical(fam_names, famn), TRUE,
      sprintf("Model.csv: %s", paste(sprintf("%s = %s", fam, fam_names), collapse = ", ")))
  fam_in <- fam[fam %in% ids]
  fam_out <- setdiff(fam, ids)
  rec("depmap_family_expression_present", length(fam_in) >= 1, TRUE,
      sprintf("in the 24Q4 expression matrix: %s; absent: %s", paste(fam_in, collapse = ","),
              if (length(fam_out)) paste(sprintf("%s (%s)", fam_out, famn[match(fam_out, fam)]), collapse = ",") else "none"))
  nm_of <- function(id) model$CellLineName[match(id, model$ModelID)]

  # gene matching: GENCODE v44 protein-coding gene_name against DepMap's symbol
  chrs <- paste0("chr", c(1:22, "X", "Y", "M"))
  pc <- genes[genes$gene_type == "protein_coding" & genes$chr %in% chrs & !grepl("_PAR_Y$", genes$gene_id), ]
  un <- pc$gene_name[!(pc$gene_name %in% pc$gene_name[duplicated(pc$gene_name)])]
  pcu <- pc[pc$gene_name %in% un, ]
  matched <- intersect(pcu$gene_name, colnames(X))
  histone <- matched[grepl(HIST_RX, matched)]
  gset <- setdiff(matched, histone)
  gid <- setNames(pcu$gene_id, pcu$gene_name)
  rec("matched_genes_in_count_table", all(gid[gset] %in% cnt$id), TRUE, "every matched gene id is in the count table")
  GM <- list(depmap_gene_columns = length(dsym), depmap_symbols_unique = ncol(X), depmap_symbols_duplicated_dropped = sum(dupd),
             gencode_pc_genes = nrow(pc), gencode_pc_names_unique = nrow(pcu), matched_by_name = length(matched),
             histone_excluded = length(histone), used = length(gset))
  SUMMARY$gene_matching <<- GM
  say("gene matching: ", paste(sprintf("%s %s", names(GM), unlist(GM)), collapse = ", "))

  # DepMap anchors for the loci
  anc_ids <- unique(c(fam_in, model$ModelID[model$StrippedCellLineName == pv("depmap_mgmt_posctrl")]))
  anc_ids <- anc_ids[anc_ids %in% ids]
  ag <- intersect(c(ps("locus_targets"), "GAPDH"), colnames(X))
  AN <- do.call(rbind, lapply(ag, function(g) data.frame(gene = g, model = c(anc_ids, "pan-cancer median"),
                                                          name = c(nm_of(anc_ids), "all rows"),
                                                          log2_tpm1 = round(c(X[anc_ids, g], med[g]), 4),
                                                          tpm = round(2^c(X[anc_ids, g], med[g]) - 1, 3))))
  wt(AN, "depmap_anchor_genes.tsv")

  # ------------------------------------------------------------------------------------------ item 6: identity
  cns <- intersect(model$ModelID[model$OncotreeLineage == pv("depmap_lineage")], ids)
  Vc <- X[cns, gset, drop = FALSE]
  vv <- apply(Vc, 2, var, na.rm = TRUE)
  ntop <- min(pn("identity_top_var_genes"), length(vv))
  topg <- names(sort(vv, decreasing = TRUE))[seq_len(ntop)]
  Vcc <- sweep(Vc, 2, med[gset])
  ident <- function(lib) {
    ya <- log2(tpm$m[gid[topg], lib] + 1)
    yb <- log2(tpm$m[gid[gset], lib] + 1) - med[gset]
    ra <- apply(Vc[, topg, drop = FALSE], 1, function(r) suppressWarnings(cor(r, ya, method = "spearman", use = "complete.obs")))
    rb <- apply(Vcc, 1, function(r) suppressWarnings(cor(r, yb, use = "complete.obs")))
    data.frame(library = lib, model = cns, name = nm_of(cns), u251_family = cns %in% fam, rho_topvar = round(ra, 5),
               rank_topvar = rank(-ra, ties.method = "min"), r_centred = round(rb, 5), rank_centred = rank(-rb, ties.method = "min"))
  }
  IDT <- do.call(rbind, lapply(LIBS, ident))
  wt(IDT[order(match(IDT$library, LIBS), IDT$rank_topvar), ], "depmap_identity.tsv")
  byl <- do.call(rbind, lapply(LIBS, function(l) {
    d <- IDT[IDT$library == l, ]
    a <- d[order(d$rank_topvar), ]; b <- d[order(d$rank_centred), ]
    data.frame(library = l, group = group_of(l), top_topvar = a$name[1], top_topvar_in_family = a$u251_family[1],
               best_family_rank_topvar = min(d$rank_topvar[d$u251_family]),
               margin_topvar = round(max(d$rho_topvar[d$u251_family]) - max(d$rho_topvar[!d$u251_family]), 5),
               top_centred = b$name[1], top_centred_in_family = b$u251_family[1],
               best_family_rank_centred = min(d$rank_centred[d$u251_family]),
               margin_centred = round(max(d$r_centred[d$u251_family]) - max(d$r_centred[!d$u251_family]), 5))
  }))
  wt(byl, "depmap_identity_by_library.tsv")
  c2 <- byl[byl$library == HCAL, ]
  SUMMARY$identity <<- list(cns_lines = length(cns), top_var_genes = ntop, family_in_matrix = fam_in, family_absent = fam_out,
                           by_library = byl)
  REPORT_EXTRA <<- c(REPORT_EXTRA, "## Identity against DepMap CNS lines", "",
                     sprintf("%d CNS/Brain lines with 24Q4 expression; U-251 family in the matrix: %s; absent: %s.", length(cns),
                             paste(sprintf("%s (%s)", fam_in, nm_of(fam_in)), collapse = ", "), paste(fam_out, collapse = ", ")),
                     "", md_table(byl, 4), "")
  rec("identity_C2B_nearest_U251_family", isTRUE(c2$top_topvar_in_family), FALSE,
      sprintf("C2B nearest CNS line (Spearman, %d most variable genes): %s; family margin %s", ntop, c2$top_topvar, f4(c2$margin_topvar)))
  rec("identity_C2B_nearest_U251_family_centred", isTRUE(c2$top_centred_in_family), FALSE,
      sprintf("C2B nearest CNS line (Pearson, centred on the pan-cancer median, %d genes): %s; margin %s", length(gset),
              c2$top_centred, f4(c2$margin_centred)))
  # Both metrics must agree: in the dry run the top-variance Spearman also put the two mostly-rat controls nearest
  # U-251 MG (margins 0.01 and 0.04), so on its own it is not a strict enough gate.
  if (!isTRUE(c2$top_topvar_in_family) || !isTRUE(c2$top_centred_in_family)) {
    SUMMARY$identity$stop <<- "C2B is not nearest to the U-251 family on both metrics: arm dosage not computed; look at sample identity first"
    finish("STOP_identity", 5)
  }

  # ------------------------------------------------------------------------------------------ item 5: arm dosage
  cb <- read.delim(gzfile(file.path(REFS, "cytoBand.hg38.txt.gz")), header = FALSE,
                   col.names = c("chr", "start", "end", "band", "stain"), colClasses = c("character", "numeric", "numeric", "character", "character"))
  achr <- paste0("chr", c(1:22, "X", "Y"))
  cen <- sapply(achr, function(ch) max(cb$end[cb$chr == ch & cb$stain == "acen" & substr(cb$band, 1, 1) == "p"]))
  clen <- sapply(achr, function(ch) max(cb$end[cb$chr == ch]))
  rec("cytoband_centromeres", all(is.finite(cen)) && all(cen > 0 & cen < clen), TRUE, sprintf("%d chromosomes with a p-arm acen band", sum(is.finite(cen))))
  ARMS <- as.vector(t(outer(c(1:22, "X", "Y"), c("p", "q"), paste0)))
  ga <- pcu[match(gset, pcu$gene_name), ]
  gmid <- (ga$start + ga$end) / 2
  garm <- ifelse(ga$chr %in% achr, paste0(sub("^chr", "", ga$chr), ifelse(gmid <= cen[ga$chr], "p", "q")), NA)
  names(garm) <- gset
  amin <- pn("arm_min_genes")
  arm_scores <- function(x, keep) {
    ok <- keep & !is.na(garm) & !is.na(x)
    s <- sapply(ARMS, function(a) { v <- x[ok & garm == a]; if (length(v) >= amin) median(v) else NA })
    n <- sapply(ARMS, function(a) sum(ok & garm == a))
    list(score = s - median(s, na.rm = TRUE), raw = s, n = n)
  }
  tmin <- pn("tpm_min")
  ours <- function(lib, mk) {
    tv <- tpm$m[gid[gset], lib]
    arm_scores(log2(tv + 1) - med[gset], tv >= tmin & !mk[gid[gset], lib])
  }
  AS <- list(); AN_ <- list()
  for (l in LIBS) { r <- ours(l, mark); AS[[l]] <- r$score; AN_[[l]] <- r$n }
  dm_lab <- function(id) paste0("DepMap_", model$StrippedCellLineName[match(id, model$ModelID)], "_", id)
  for (id in fam_in) {
    xv <- X[id, gset]
    r <- arm_scores(xv - med[gset], xv >= log2(tmin + 1))
    AS[[dm_lab(id)]] <- r$score; AN_[[dm_lab(id)]] <- r$n
  }
  ASd <- data.frame(arm = ARMS, sapply(AS, function(v) round(v, 5)), check.names = FALSE)
  ANd <- data.frame(arm = ARMS, sapply(AN_, identity), check.names = FALSE)
  wt(ASd, "arm_dosage.tsv"); wt(ANd, "arm_dosage_ngenes.tsv")
  sens <- NULL
  if (MODE == "published") {
    mx <- function(d) if (all(is.na(d))) NA else round(max(abs(d), na.rm = TRUE), 5)
    sens <- do.call(rbind, lapply(LIBS, function(l) {
      a0 <- AS[[l]]
      data.frame(library = l, arms_scored = sum(!is.na(a0)),
                 max_abs_change_hgraft_marks = mx(ours(l, mark_hg)$score - a0),
                 max_abs_change_bothresid_marks = mx(ours(l, mark_rs)$score - a0),
                 max_abs_change_ctrlmax_marks = mx(ours(l, mark_cx)$score - a0))
    }))
    wt(sens, "arm_dosage_mark_sensitivity.tsv")
  }

  # DepMap arm copy number: the positive control for the dosage method (DepMap's own expression against its own CN)
  dmp <- read.csv(file.path(REFS, "OmicsDefaultModelProfiles.csv"), colClasses = "character")
  cnpath <- file.path(REFS, "OmicsCNSegmentsProfile.csv")
  cnres <- NULL
  if (file.exists(cnpath) && "OmicsCNSegmentsProfile.csv" %in% st$name[st$status == "ok"]) {
    seg <- read.csv(cnpath, colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "character"))
    seg$chr <- paste0("chr", sub("^chr", "", seg$Chromosome))
    dna <- dmp[dmp$ProfileType == "dna", ]
    arm_cn <- function(pid) {
      s <- seg[seg$ProfileID == pid, ]
      if (nrow(s) == 0) return(rep(NA, length(ARMS)))
      v <- sapply(ARMS, function(a) {
        ch <- paste0("chr", sub("[pq]$", "", a)); p <- substring(a, nchar(a))
        a0 <- if (p == "p") 0 else cen[ch]; a1 <- if (p == "p") cen[ch] else clen[ch]
        t <- s[s$chr == ch, ]
        w <- pmax(0, pmin(t$End, a1) - pmax(t$Start, a0))
        if (sum(w) <= 0) return(NA)
        log2(max(sum(w * t$SegmentMean) / sum(w), 1e-3))
      })
      v - median(v, na.rm = TRUE)
    }
    cn_ids <- intersect(cns, dna$ModelID)
    CN <- sapply(cn_ids, function(id) arm_cn(dna$ProfileID[dna$ModelID == id][1]))
    rownames(CN) <- ARMS
    fam_cn <- intersect(fam, cn_ids)
    wt(data.frame(arm = ARMS, round(CN[, fam_cn, drop = FALSE], 5), check.names = FALSE), "depmap_arm_cn_family.tsv")
    sp <- function(a, b) { ok <- !is.na(a) & !is.na(b); if (sum(ok) < 10) return(c(NA, sum(ok))); c(cor(a[ok], b[ok], method = "spearman"), sum(ok)) }
    cnres <- do.call(rbind, lapply(names(AS), function(s) do.call(rbind, lapply(fam_cn, function(f) {
      z <- sp(AS[[s]], CN[, f]); data.frame(sample = s, cn_of = paste0(nm_of(f), " ", f), spearman = round(z[1], 4), arms = z[2])
    }))))
    wt(cnres, "arm_dosage_vs_depmap_cn.tsv")
    ic <- sapply(cn_ids, function(f) sp(AS[[HCAL]], CN[, f])[1])
    icd <- data.frame(model = cn_ids, name = nm_of(cn_ids), u251_family = cn_ids %in% fam, spearman_C2B_dosage_vs_cn = round(ic, 4),
                      rank = rank(-ic, ties.method = "min", na.last = "keep"))
    wt(icd[order(icd$rank), ], "identity_c2b_dosage_vs_cns_cn.tsv")
    pos <- cnres[cnres$sample == dm_lab(fam[1]) & grepl(fam[1], cnres$cn_of), ]
    SUMMARY$arm_dosage_cn <<- list(positive_control = pos, c2b_best_cn_family_rank = min(icd$rank[icd$u251_family], na.rm = TRUE),
                                  cns_lines_with_cn = length(cn_ids))
    rec("dosage_method_positive_control", nrow(pos) == 1 && isTRUE(pos$spearman > 0), FALSE,
        sprintf("DepMap %s: its arm expression dosage against its own arm copy number, Spearman %s over %s arms", fam[1],
                if (nrow(pos)) f4(pos$spearman) else "NA", if (nrow(pos)) pos$arms else "NA"))
  } else {
    rec("depmap_cn_segments_available", FALSE, FALSE, "OmicsCNSegmentsProfile.csv missing: no copy-number positive control")
  }
  SUMMARY$arm_dosage <<- list(arms_min_genes = amin, tpm_min = tmin, arms_scored = as.list(sapply(AS, function(v) sum(!is.na(v)))),
                             mark_sensitivity = sens)
  REPORT_EXTRA <<- c(REPORT_EXTRA, "## Arm dosage (median-arm centred log2 expression relative to the DepMap pan-cancer median)", "",
                     md_table(ASd, 3), "", "## Arm dosage against DepMap arm copy number (Spearman over arms)", "",
                     md_table(cnres, 3), "")

  # ------------------------------------------------------------------------------------------ CCLE RRBS and EPIC checks
  if (MODE == "published") {
    rr <- file.path(REFS, "CCLE_RRBS_TSS1kb_20181022.txt.gz")
    if ("CCLE_RRBS_TSS1kb_20181022.txt.gz" %in% st$name[st$status == "ok"]) {
      LL <- readLines(gzfile(rr), warn = FALSE)
      hh <- strsplit(LL[1], "\t", fixed = TRUE)[[1]]
      pref <- ps("rrbs_loci_prefixes")
      rows <- LL[grepl(paste0("^(", paste(pref, collapse = "|"), ")"), LL)]
      rm(LL)
      M <- t(sapply(strsplit(rows, "\t", fixed = TRUE), function(p) { length(p) <- length(hh); p }))
      cc <- c(1:3, grep("_CENTRAL_NERVOUS_SYSTEM$", hh))
      RR <- as.data.frame(M[, cc, drop = FALSE]); names(RR) <- hh[cc]
      RR[] <- lapply(RR, trimws)
      wt(RR, "ccle_rrbs_loci_cns.tsv")
      uc <- grep("^U251MG_", names(RR), value = TRUE)
      SUMMARY$ccle_rrbs <<- list(cns_columns = length(cc) - 3, u251_column = uc,
                                u251_values = if (length(uc)) as.list(setNames(RR[[uc[1]]], RR[[1]])) else NULL)
    }
    ep <- file.path(REFS, "EPIC.hg38.manifest.gencode.v36.tsv.gz")
    if ("EPIC.hg38.manifest.gencode.v36.tsv.gz" %in% st$name[st$status == "ok"]) {
      eh <- readLines(gzfile(ep), n = 1)
      rec("epic_manifest_header", grepl("probeID", eh, fixed = TRUE), FALSE, substr(eh, 1, 160))
    }
  }

  # ------------------------------------------------------------------------------------------ graft mode: two leak estimates
  if (MODE == "graft") {
    reg("counts_published", pv("counts_published"))
    pub <- read_gene_tab(pv("counts_published"), LIBS)
    Pm <- pub$m[intersect(pub$id, cnt$id), ]; Gm <- cnt$m[rownames(Pm), ]
    rg_ <- Gm[, HCAL] / Pm[, HCAL]
    elig <- rownames(Pm)[which(Pm[, HCAL] >= 50 & rg_ >= 0.2)]   # the rnaseq_graft check 7 rule
    pctrl_pub <- rowMeans(sapply(RCAL, function(c) Pm[, c] / hf[c]))
    Rpub <- (1 - h_host) * tot * (ratm[["graft"]] + ratm[["both"]])
    two <- do.call(rbind, lapply(LIBS, function(t) {
      rat_counts <- Pm[elig, t] - Gm[elig, t] / rg_[elig]
      e_model <- pctrl_pub[elig] * Rpub[t]
      sP <- sum(Pm[elig, t])
      data.frame(library = t, group = group_of(t), eligible_genes = length(elig),
                 rat_share_counts_estimate = round(sum(rat_counts) / sP, 5),
                 rat_share_rate_model = round(sum(e_model) / sP, 5),
                 spearman_gene_level = if (t %in% c(PRIM, RECU)) round(suppressWarnings(cor(rat_counts, e_model, method = "spearman")), 4) else NA)
    }))
    ck <- pv("graft_checks"); reg("rnaseq_graft_checks.json", ck)
    jt <- paste(readLines(ck, warn = FALSE), collapse = "\n")
    two$check7_rat_share_signed <- sapply(two$library, function(t) {
      m <- regmatches(jt, regexpr(paste0('"', t, '": \\{[^}]*"rat_share_signed": [-0-9.eE+]+'), jt))
      if (length(m) == 0) NA else as.numeric(sub('.*"rat_share_signed": ', "", m))
    })
    wt(two, "rat_leak_two_estimates.tsv")
    ok7 <- !all(is.na(two$check7_rat_share_signed)) &&
      max(abs(two$check7_rat_share_signed - two$rat_share_counts_estimate), na.rm = TRUE) < 1e-4
    rec("check7_reproduced", ok7, FALSE, "count-based rat share recomputed here equals rnaseq_graft check 7 (to 1e-4)")
    SUMMARY$two_leak_estimates <<- two
    pa <- file.path(dirname(OUT), "p0", "arm_dosage.tsv")
    if (file.exists(pa)) {
      a0 <- read.delim(pa, check.names = FALSE)
      cmp <- do.call(rbind, lapply(LIBS, function(l) {
        x <- a0[[l]][match(ARMS, a0$arm)]; y <- AS[[l]]; ok <- !is.na(x) & !is.na(y)
        data.frame(library = l, arms_both = sum(ok), spearman = if (sum(ok) >= 10) round(cor(x[ok], y[ok], method = "spearman"), 4) else NA,
                   max_abs_diff = if (any(ok)) round(max(abs(x[ok] - y[ok])), 4) else NA)
      }))
      wt(cmp, "arm_dosage_graft_vs_published.tsv")
      SUMMARY$arm_dosage_graft_vs_published <<- cmp
    }
    REPORT_EXTRA <<- c(REPORT_EXTRA, "## Rat leak: rate model against the count-based estimate (rnaseq_graft check 7)", "",
                       md_table(two, 4), "")
  }

  REPORT_EXTRA <<- c("## Human share (calibrated) and rat leak", "", md_table(H[, c("library", "group", "raw_graft_frac", "h_graft", "h_host",
                     "h_diff", "both_obs_minus_pred", "h_chrY")], 4), "",
                     md_table(LB[, c("library", "rat_share_of_counts_model", "rat_share_of_counts_bothresid",
                                     paste0("rat_share_of_counts_", RCAL, "_profile"), "genes_tpm_ge_min",
                                     "genes_tpm_ge_min_marked", "marks_changed_if_bothresid", "marks_added_if_either_ctrl")], 4),
                     "", REPORT_EXTRA)
  finish(if (length(FAILED)) "complete_with_failed_checks" else "complete", 0)
}

tryCatch(main(), error = function(e) {
  msg <- conditionMessage(e)
  say("ERROR: ", msg)
  CHECKS[["unexpected_error"]] <<- list(pass = FALSE, hard = TRUE, detail = msg)
  FAILED <<- c(FAILED, "unexpected_error")
  finish(paste0("error: ", msg), 1)
})
