# inv_common.R: helpers shared by the invasion_niche steps. Base R only, so it runs in both images
# (DESeq2 1.34.0 image, R 4.1.3; go2432-bioconductor.sif, R 4.3.3).

read_params <- function(f) {
  l <- readLines(f, warn = FALSE)
  l <- sub("[[:space:]]+#.*$", "", l)
  l <- l[grepl("^[A-Za-z0-9_]+:", l)]
  k <- sub(":.*$", "", l)
  v <- trimws(sub("^[^:]+:", "", l))
  v <- gsub('^"|"$', "", v)
  stopifnot(!any(duplicated(k)))
  setNames(as.list(v), k)
}

splitv <- function(x, sep = ",") {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, sep, fixed = TRUE)[[1]])
}

wt <- function(df, f) {
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  if (grepl("\\.gz$", f)) {
    con <- gzfile(f, "w"); on.exit(close(con))
    write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  } else {
    write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  invisible(f)
}

# key/value lines appended to results/summary.tsv (part, key, value)
kv_file <- NULL
kv <- function(part, key, value) {
  if (is.null(kv_file)) return(invisible())
  v <- if (length(value) > 1) paste(value, collapse = ",") else as.character(value)
  cat(part, "\t", key, "\t", gsub("[\t\n]", " ", v), "\n", sep = "", file = kv_file, append = TRUE)
}

read_counts <- function(f) {
  d <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(identical(colnames(d)[1:2], c("gene_id", "gene_name")))
  m <- as.matrix(d[, -(1:2), drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- d$gene_id
  stopifnot(!any(duplicated(rownames(m))))
  list(m = m, name = setNames(d$gene_name, d$gene_id))
}

read_gene_table <- function(f) {
  g <- read.delim(f, stringsAsFactors = FALSE, quote = "", comment.char = "")
  stopifnot(all(c("gene_id", "chr", "gene_type", "gene_name") %in% colnames(g)), !any(duplicated(g$gene_id)))
  g
}

# xengsort classify logs: the tab-separated line under "prefix host graft ambiguous both neither"
read_xeng <- function(dir, libs) {
  do.call(rbind, lapply(libs, function(l) {
    x <- readLines(file.path(dir, paste0(l, ".txt")), warn = FALSE)
    i <- grep("^prefix\thost\tgraft\tambiguous\tboth\tneither", x)
    if (length(i) != 1) stop("xengsort log for ", l, ": header line not found exactly once")
    f <- strsplit(x[i + 1], "\t")[[1]]
    if (length(f) != 6 || f[1] != l) stop("xengsort log for ", l, ": unexpected line ", x[i + 1])
    v <- as.numeric(f[2:6])
    data.frame(lib = l, host = v[1], graft = v[2], ambiguous = v[3], both = v[4], neither = v[5],
               total = sum(v), stringsAsFactors = FALSE)
  }))
}

read_picard <- function(f) {
  x <- readLines(f, warn = FALSE)
  i <- grep("^LIBRARY\t", x)
  if (length(i) != 1) stop("Picard metrics header not found in ", f)
  h <- strsplit(x[i], "\t")[[1]]; v <- strsplit(x[i + 1], "\t")[[1]]
  v <- c(v, rep("", length(h) - length(v)))
  out <- setNames(as.list(v), h)
  for (k in setdiff(h, "LIBRARY")) out[[k]] <- suppressWarnings(as.numeric(out[[k]]))
  out
}

# exact (Garwood) Poisson interval for an integer count
pois_ci <- function(y, conf = 0.95) {
  a <- (1 - conf) / 2
  lo <- ifelse(y == 0, 0, qchisq(a, 2 * y) / 2)
  hi <- qchisq(1 - a, 2 * (y + 1)) / 2
  cbind(lo = lo, hi = hi)
}

# the ten balanced 3-v-3 partitions of the six tumours; IL67B always in group A so no split appears as its mirror.
# Labels follow host_verify_de.R (split_<group B in tumour order>), so host fits join by label.
make_splits <- function(prim, recur, anchor = "IL67B") {
  tum <- c(prim, recur)
  stopifnot(anchor %in% tum, length(prim) == 3, length(recur) == 3)
  sp <- lapply(combn(setdiff(tum, anchor), 2, simplify = FALSE), function(p) {
    A <- tum[tum %in% c(anchor, p)]; B <- setdiff(tum, A)
    list(A = A, B = B, label = paste0("split_", paste(B, collapse = "_")), true = setequal(B, recur))
  })
  names(sp) <- vapply(sp, function(s) s$label, "")
  stopifnot(length(sp) == 10, sum(vapply(sp, function(s) s$true, TRUE)) == 1)
  sp
}
true_label <- function(sp) names(sp)[vapply(sp, function(s) s$true, TRUE)]

welch_t <- function(x, A, B) {
  a <- x[A]; b <- x[B]
  if (any(is.na(c(a, b)))) return(NA_real_)
  se <- sqrt(var(a) / length(a) + var(b) / length(b))
  if (!is.finite(se) || se == 0) return(NA_real_)
  (mean(b) - mean(a)) / se
}
split_stats <- function(x, sp) vapply(sp, function(s) welch_t(x, s$A, s$B), 0)

# rank of the true split's |statistic| among the ten (1 = most extreme); ties give a range
rank_true <- function(s, lab) {
  a <- abs(s); t <- a[lab]
  if (is.na(t)) return(c(rank_min = NA_real_, rank_max = NA_real_))
  c(rank_min = 1 + sum(a > t, na.rm = TRUE), rank_max = sum(a >= t, na.rm = TRUE))
}
fmt_rank <- function(r) {
  if (is.na(r[1])) return("NA")
  if (r[1] == r[2]) sprintf("%d of 10", as.integer(r[1])) else sprintf("%d-%d of 10", as.integer(r[1]), as.integer(r[2]))
}

# exact non-negative least squares for a handful of columns: the optimum is the least-squares fit on its own support
# with positive coefficients, so the feasible support with the smallest residual is the NNLS solution.
nnls_enum <- function(y, X) {
  k <- ncol(X); stopifnot(k <= 10)
  best_rss <- sum(y^2); best_w <- setNames(rep(0, k), colnames(X))
  for (m in seq_len(2^k - 1)) {
    s <- which(bitwAnd(m, 2^(0:(k - 1))) > 0)
    f <- lm.fit(X[, s, drop = FALSE], y)
    b <- f$coefficients
    if (any(is.na(b)) || any(b < 0)) next
    rss <- sum(f$residuals^2)
    if (rss < best_rss) { best_rss <- rss; best_w <- setNames(rep(0, k), colnames(X)); best_w[s] <- b }
  }
  list(w = best_w, share = if (sum(best_w) > 0) best_w / sum(best_w) else best_w, rss = best_rss,
       r2 = 1 - best_rss / sum((y - mean(y))^2))
}

fnum <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, digits = d, format = "f", big.mark = ","))
fint <- function(x) ifelse(is.na(x), "NA", formatC(round(x), format = "d", big.mark = ","))
fsig <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(signif(x, d), format = "g"))

md5_of <- function(files) {
  files <- files[nzchar(files)]
  e <- file.exists(files)
  out <- rep(NA_character_, length(files))
  if (any(e)) out[e] <- unname(tools::md5sum(files[e]))
  data.frame(file = files, exists = e, bytes = ifelse(e, file.size(files), NA), md5 = out, stringsAsFactors = FALSE)
}

md_table <- function(df, digits = 3) {
  if (is.null(df) || !nrow(df)) return("(no rows)")
  f <- function(v) if (is.numeric(v)) ifelse(is.na(v), "NA", formatC(v, digits = digits, format = "g")) else as.character(v)
  cells <- lapply(df, f)
  c(paste0("| ", paste(colnames(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    vapply(seq_len(nrow(df)), function(i) paste0("| ", paste(vapply(cells, function(c) c[i], ""), collapse = " | "), " |"), ""))
}
