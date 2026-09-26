lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
suppressMessages({library(edgeR); library(RNASeqPower); library(ssizeRNA); library(PROPER); library(jsonlite)})
outdir <- commandArgs(TRUE)[1]
d <- readRDS(file.path(outdir, "disp_all6.rds"))
mu <- as.numeric(d$mu); disp <- as.numeric(d$disp)
G <- length(mu); cat("genes", G, "median mu", median(mu), "median BCV", median(sqrt(disp)), "\n")
res <- list(G = G)

# ---- Jung (2005) per-test alpha that yields FDR f with r1 = power * m1 expected true rejections
jung <- function(m, pi0, power = 0.8, f = 0.05) { m1 <- round(m * (1 - pi0)); m0 <- m - m1; r1 <- power * m1; r1 * f / (m0 * (1 - f)) }
alphas <- c(`pi0=0.99` = jung(G, 0.99), `pi0=0.95` = jung(G, 0.95), `pi0=0.90` = jung(G, 0.90), `single gene 0.05` = 0.05)
print(alphas); res$jung_alpha <- as.list(alphas)

# ---- Hart et al. 2013 (RNASeqPower): n = 2 (z_{1-a/2} + z_{power})^2 (1/depth + cv^2) / (ln fc)^2
depths <- c(5, 10, 20, 50, 100, 463)
cvs <- c(0.10, 0.20, 0.269, 0.40, 0.53)
hart <- list()
for (a in names(alphas)) {
  tab <- rnapower(depth = depths, cv = cvs, effect = 1.5, alpha = alphas[[a]], power = 0.8)
  tab <- drop(tab); print(a); print(ceiling(tab))
  hart[[a]] <- ceiling(tab)
}
res$hart <- lapply(hart, function(t) { t <- as.data.frame(t); t$depth <- rownames(t); t })

# ---- ssizeRNA (Bi & Liu 2016): gene-wise mean and dispersion from the six tumours
set.seed(20260926)
pdf(file.path(outdir, "ssizeRNA_curves.pdf"))
ss <- list()
for (p0 in c(0.99, 0.95, 0.90)) {
  r <- ssizeRNA_vary(nGenes = G, pi0 = p0, m = 200, mu = mu, disp = disp, fc = 1.5, up = 0.5, replace = TRUE,
                     fdr = 0.05, power = 0.8, maxN = 60)
  cat("ssizeRNA_vary pi0", p0, "n/group", r$ssize, "\n")
  ss[[paste0("vary_pi0_", p0)]] <- list(ssize = r$ssize, power = r$power)
}
r1 <- ssizeRNA_single(nGenes = G, pi0 = 0.95, m = 200, mu = median(mu), disp = 0.0725, fc = 1.5, fdr = 0.05, power = 0.8, maxN = 60)
cat("ssizeRNA_single median mu, common disp, pi0 0.95:", r1$ssize, "\n")
ss$single_median_mu_common_disp_pi0_0.95 <- list(ssize = r1$ssize, power = r1$power)
dev.off()
res$ssizeRNA <- ss

write_json(res, file.path(outdir, "power_rnaseq.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
