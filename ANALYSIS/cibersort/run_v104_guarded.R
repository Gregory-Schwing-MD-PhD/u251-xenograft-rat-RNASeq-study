# CIBERSORT R script v1.04 (vendor/CIBERSORT.R, unchanged) with the permutation null drawn outside CIBERSORT(), so that
# a draw no nu-SVR can fit is recorded instead of ending the run.
#   Rscript run_v104_guarded.R <signature.txt> <mixture.txt> <out.tsv> <perm>        (CIBERSORT_R=<path> overrides vendor/)
#
# Why: CoreAlg fits nu = 0.25, 0.5, 0.75, clips negative weights to 0 and divides by their sum. When all three fits
# have no positive weight, w = 0/0, all three RMSEs are NaN, which.min() returns integer(0) and out[[mn]] stops with
# "attempt to select less than one element in get1index", killing every sample. With 22 columns (LM22) that almost
# never happens; with 4 (Bowman 2016 myeloid) about one random draw in 1,500 does, and seed 42 met one within the
# first minute of the host run (job 40473363_3).
#
# What is the same as CIBERSORT(): the input lines (read, de-duplicate names, order, anti-log when max < 50, no quantile
# normalisation, intersect, global standardisation of the signature) are copied from CIBERSORT() lines 134-178; the
# draws are doPerm()'s, from the same seed (42) in the same order; every fit is CoreAlg() itself; P is CIBERSORT()'s
# formula. What differs: a draw where CoreAlg() stops is left out of the null (the CIBERSORT null is sort()ed, which
# drops NA, so this is where such a draw would land if CoreAlg() returned NA), and the count is written out. Leaving
# them out is the conservative choice: a draw with no fit cannot beat a sample that has one. A sample where CoreAlg()
# stops is written with NA fractions rather than stopping the others.
# Check: with no degenerate draw the output must equal run_v104.R's to rounding (< 1e-12). The job log compares the
# fractions, r and RMSE with CIBERSORT(perm = 0); the host LM22 and Zhang runs are repeated to compare P.
args <- commandArgs(trailingOnly = TRUE)
sig <- normalizePath(args[1]); mix <- normalizePath(args[2]); out <- args[3]; perm <- as.integer(args[4])
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
cs <- Sys.getenv("CIBERSORT_R", file.path(here, "vendor", "CIBERSORT.R"))
source(cs)
cat("CIBERSORT header:", readLines(cs, n = 1), "\n")
cat("e1071", as.character(packageVersion("e1071")), " preprocessCore", as.character(packageVersion("preprocessCore")), "\n")

# --- CIBERSORT() lines 134-178, QN = FALSE ---
X <- read.table(sig, header = T, sep = "\t", row.names = 1, check.names = F)
Y <- read.table(mix, header = T, sep = "\t", check.names = F)
dups <- dim(Y)[1] - length(unique(Y[, 1]))
if (dups > 0) {
  warning(paste(dups, " duplicated gene symbol(s) found in mixture file!", sep = ""))
  rownames(Y) <- make.names(Y[, 1], unique = TRUE)
} else {rownames(Y) <- Y[, 1]}
Y <- Y[, -1]
X <- data.matrix(X); Y <- data.matrix(Y)
X <- X[order(rownames(X)), ]; Y <- Y[order(rownames(Y)), ]
if (max(Y) < 50) {Y <- 2^Y}
Xgns <- row.names(X); Ygns <- row.names(Y)
Y <- Y[Ygns %in% Xgns, ]
X <- X[Xgns %in% row.names(Y), ]
X <- (X - mean(X)) / sd(as.vector(X))
cat("genes used:", nrow(X), " columns:", ncol(X), " samples:", ncol(Y), "\n")

core <- function(y) tryCatch(CoreAlg(X, y, FALSE, "sig.score"), error = function(e) e)
failed <- function(r) inherits(r, "error")

# --- doPerm(), same draws ---
set.seed(42)
t0 <- Sys.time()
Ylist <- as.list(data.matrix(Y))
nullr <- rep(NA_real_, perm); degenerate <- integer(0)
for (i in seq_len(perm)) {
  yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
  yr <- (yr - mean(yr)) / sd(yr)
  r <- core(yr)
  if (failed(r)) degenerate <- c(degenerate, i) else nullr[i] <- r$mix_r
}
nulldist <- sort(nullr)
cat("null draws:", perm, " fitted:", length(nulldist), " degenerate (no positive weight at any nu):", length(degenerate),
    if (length(degenerate)) paste0("[draw ", paste(degenerate, collapse = ", "), "]") else "", "\n")

# --- CIBERSORT() sample loop ---
rows <- list()
for (j in seq_len(ncol(Y))) {
  y <- Y[, j]; y <- (y - mean(y)) / sd(y)
  r <- core(y)
  if (failed(r)) {
    cat("sample", colnames(Y)[j], "could not be fitted:", conditionMessage(r), "\n")
    w <- rep(NA_real_, ncol(X)); p <- NA_real_; rr <- NA_real_; rmse <- NA_real_
  } else {
    w <- as.numeric(r$w); rr <- r$mix_r; rmse <- r$mix_rmse
    p <- if (perm > 0) 1 - (which.min(abs(nulldist - rr)) / length(nulldist)) else 9999
  }
  rows[[j]] <- c(w, p, rr, rmse)
}
res <- do.call(rbind, rows)
colnames(res) <- c(colnames(X), "P-value", "Correlation", "RMSE")
rownames(res) <- colnames(Y)
cat("seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")

# --- check the fractions against CIBERSORT() itself (perm = 0: no draws, so it can only stop on a sample) ---
tmp <- tempfile("cs_"); dir.create(tmp); old <- setwd(tmp)
ref <- tryCatch(CIBERSORT(sig, mix, perm = 0, QN = FALSE), error = function(e) e)
setwd(old)
if (inherits(ref, "error")) {
  cat("check: CIBERSORT(perm = 0) stopped:", conditionMessage(ref), "\n")
} else {
  k <- c(colnames(X), "Correlation", "RMSE")
  cat("check: max |guarded - CIBERSORT(perm = 0)| over fractions, r, RMSE =",
      format(max(abs(res[rownames(ref), k] - ref[, k]), na.rm = TRUE), digits = 3), "\n")
}

out_df <- data.frame(sample = rownames(res), res, null_fitted = length(nulldist), null_degenerate = length(degenerate),
                     check.names = FALSE)
write.table(out_df, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(out_df, digits = 3)
