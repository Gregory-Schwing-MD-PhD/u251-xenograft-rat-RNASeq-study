# Run CIBERSORT R script v1.04 (vendor/CIBERSORT.R, kept out of git under Stanford's licence) on one signature.
#   Rscript run_v104.R <signature.txt> <mixture.txt> <out.tsv> <perm>
# QN = FALSE (RNA-seq), relative mode, seed 42. The script writes CIBERSORT-Results.txt into the working directory,
# so each call runs in its own temporary directory.
args <- commandArgs(trailingOnly = TRUE)
sig <- normalizePath(args[1]); mix <- normalizePath(args[2]); out <- args[3]; perm <- as.integer(args[4])
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
source(file.path(here, "vendor", "CIBERSORT.R"))
cat("CIBERSORT header:", readLines(file.path(here, "vendor", "CIBERSORT.R"), n = 1), "\n")
cat("e1071", as.character(packageVersion("e1071")), " preprocessCore", as.character(packageVersion("preprocessCore")), "\n")
tmp <- tempfile("cs_"); dir.create(tmp); old <- setwd(tmp)
set.seed(42)
t0 <- Sys.time()
res <- CIBERSORT(sig, mix, perm = perm, QN = FALSE)
setwd(old)
cat("seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
res <- data.frame(sample = rownames(res), res, check.names = FALSE)
write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
print(res, digits = 3)
