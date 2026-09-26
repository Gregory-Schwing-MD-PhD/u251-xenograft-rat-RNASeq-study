lib <- Sys.getenv("R_LIBS_USER"); .libPaths(c(lib, .libPaths()))
suppressMessages({library(edgeR); library(PROPER); library(jsonlite)})
outdir <- commandArgs(TRUE)[1]
d <- readRDS(file.path(outdir, "disp_all6.rds")); mu <- as.numeric(d$mu); disp <- as.numeric(d$disp); G <- length(mu)
opts <- RNAseq.SimOptions.2grp(ngenes = G, lBaselineExpr = log(mu), lOD = log(disp), p.DE = 0.05,
                               lfc = function(n) sample(c(-1, 1), n, replace = TRUE) * log(1.5), sim.seed = 777)
Nreps <- c(3, 5, 8, 10, 12, 15, 20, 25)
t0 <- Sys.time()
sims <- runSims(Nreps = Nreps, sim.opts = opts, DEmethod = "edgeR", nsims = 5, verbose = TRUE)
print(Sys.time() - t0)
pw <- comparePower(sims, alpha.type = "fdr", alpha.nominal = 0.05, stratify.by = "expr", filter.by = "expr", strata.filtered = 1, target.by = "lfc", delta = log(1.5) / 2)
sp <- summaryPower(pw); print(sp)
strat <- apply(pw$power, c(1, 2), mean, na.rm = TRUE); colnames(strat) <- Nreps; print(round(strat, 3))
write_json(list(summary = as.data.frame(sp), by_expr = as.data.frame(strat), FDR = rowMeans(pw$FDR.marginal)), file.path(outdir, "proper_quick.json"), auto_unbox = TRUE, pretty = TRUE, digits = 6)
