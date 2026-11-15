#!/usr/bin/env Rscript
# invasion_niche: assemble results/REPORT.md from what steps A and C wrote. Base R. A part that did not run is reported
# as not run, with its error, never filled in.
#   Rscript inv_report.R <params.yaml>
args <- commandArgs(TRUE)
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "inv_common.R"))
P <- read_params(args[1]); OUT <- P$outdir
rd <- function(f) if (file.exists(f)) read.delim(if (grepl("\\.gz$", f)) gzfile(f) else f, check.names = FALSE, stringsAsFactors = FALSE) else NULL
S <- rbind(rd(file.path(OUT, "summary_A.tsv")), rd(file.path(OUT, "summary_C.tsv")))
g <- function(part, key) { if (is.null(S)) return("NA"); v <- S$value[S$part == part & S$key == key]; if (length(v)) v[length(v)] else "NA" }
st <- rd(file.path(OUT, "STATUS_A.tsv"))
part_ok <- function(p) !is.null(st) && any(st$part == p & st$status == "ok")
part_msg <- function(p) if (is.null(st) || !any(st$part == p)) "step A did not reach this part" else st$message[st$part == p][1]
c_ok <- file.exists(file.path(OUT, "C_programmes", "gsva_split_summary.tsv"))
L <- character(0); add <- function(...) L <<- c(L, ...)

add("# invasion_niche: U251 in the contralateral hemisphere, its state, invasion programmes by arm, host niche", "",
    sprintf("Written %s by `inv_report.R`. Results: `%s`. Scripts, parameters and the rules below: `ANALYSIS/multiomics/invasion_niche/` (README.md).",
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), OUT), "")
gate <- if (file.exists(file.path(OUT, "inputs", "gate.txt"))) readLines(file.path(OUT, "inputs", "gate.txt"), warn = FALSE) else "gate record missing"
add(paste0("**Gate.** ", paste(gate, collapse = " ")), "",
    sprintf("**Parts run.** A (Y dosage): %s. B (contralateral state): %s. C (programmes, GSVA/fgsea): %s. D (host niche): %s.",
            if (part_ok("A")) "ok" else paste("NOT RUN:", part_msg("A")), if (part_ok("B")) "ok" else paste("NOT RUN:", part_msg("B")),
            if (c_ok) "ok" else "NOT RUN (see logs/stepC)", if (part_ok("D")) "ok" else paste("NOT RUN:", part_msg("D"))), "")

add("## Claim ceiling", "",
    "- U251 RNA is present in the contralateral hemisphere of rats 68 and 69, both primary arm.",
    "- Migration cannot be separated from dissection carry-over.",
    "- Nothing can be said about LITT and invasion (no ablated rat has contralateral tissue).",
    "- With 3 primaries against 3 recurrences there are exactly 10 balanced relabellings: no primary-versus-recurrent result here can go below p = 0.10, and rank 1 of 10 is not a finding.", "")

## ---------------------------------------------------------------- A
add("## A. How much U251 is in each sample (human chrY dosage; the rats are female)", "")
yd <- rd(file.path(OUT, "A_y_dosage", "y_dosage.tsv"))
if (!is.null(yd)) {
  add(sprintf("Genes: the 11 X-degenerate Y genes, ids resolved from the GENCODE v44 GTF by name on chrY outside the PAR (`A_y_dosage/gene_ids_from_gtf.tsv`): %s.", g("A", "y_gene_ids")), "",
      sprintf("Calibration, Y counts per million graft pairs in the six tumours: graft+both %s; graft-only %s.",
              g("A", "graft_plus_both_calibration_range"), g("A", "graft_only_calibration_range")),
      sprintf("Denominator check against something outside the computation: the ten xengsort totals equal the ten SRA spot counts (SRR39547363-72): **%s**.", g("inputs", "xengsort_totals_equal_SRA_spot_counts")), "")
  t <- yd[yd$library %in% splitv(P$y_targets), c("source", "library", "y_count", "poisson95_lo", "poisson95_hi", "u251_pct_of_library",
                                                        "u251_pct_lo_poisson_and_cal", "u251_pct_hi_poisson_and_cal", "u251_pct_of_graft_bin", "preliminary_pct_range", "matches_preliminary")]
  add(md_table(t, 3), "",
      "`u251_pct_of_library` is U251 graft-bin-equivalent pairs over all pairs of the library (the calibration counts only graft pairs); the interval combines the exact Poisson interval on the Y count with the range of the six tumour calibrations. Multiplying by C2B's (graft+both)/graft converts to all U251 pairs (`u251_pct_all_pairs_c2b_conversion` in the table file).", "")
  rk <- rd(file.path(OUT, "A_y_dosage", "y_x_split_ranks.tsv"))
  if (!is.null(rk)) add("Y/X gametolog ratios (clonal marker; `y_x_ratios.tsv`), recurrent minus primary, rank of the true split among the 10:", "",
                        md_table(rk[, c("source", "pair", "mean_log2_primary", "mean_log2_recurrent", "log2_C2B", "welch_t_true", "rank")], 3), "")
  add("**Can show:** that human Y RNA, which a female host cannot make, is present in N269B and N168B at the stated shares, on counts that do not depend on xengsort's call for these reads beyond the human bin; that the published graft+both estimate is or is not reproduced.",
      "**Cannot show:** where the U251 RNA came from (migration, carry-over at dissection, index hopping); cell numbers (RNA per cell differs between states); the IL64B figure is a floor-level count and says nothing about cells. The calibration assumes the Y share of U251 RNA is the same in contralateral cells as in the tumours; a Y-negative subclone would bias it down.", "")
} else add("Not run: ", part_msg("A"), "")

## ---------------------------------------------------------------- B
add("## B. Are the contralateral U251 cells (N269B) in an in-vivo or a culture-like state?", "")
if (part_ok("B")) {
  add(sprintf("Universe: %s protein-coding genes with mean CPM >= %s over C2B and the six tumours, minus %s genes where the controls' graft background (scaled to N269B's host depth) exceeds %s of N269B's count.",
              g("B", "universe_genes"), P$universe_min_mean_cpm, g("B", "background_flagged_among_expressed"), P$background_max_share),
      sprintf("Depth-matched null: IL69B, C2B and the other tumours binomially thinned to N269B's deduplicated depth D = %s fragments (%s), %s draws.",
              g("B", "thinning_target_D"), g("B", "N269B_picard"), g("B", "thin_draws")), "")
  ax <- rd(file.path(OUT, "B_contralateral_state", "axis_position.tsv"))
  add("**Culture-to-in-vivo axis** (position 0 = C2B, 1 = the tumours; median over the axis genes):", "", md_table(ax, 3), "",
      sprintf("Rule fixed before the run: %s. N269B = %s; in-vivo null %s; culture null %s. **Verdict: %s.**",
              g("B", "verdict_rule"), g("B", "axis_t_N269B"), g("B", "axis_invivo_null_95"), g("B", "axis_culture_null_95"), g("B", "verdict")), "")
  sp <- rd(file.path(OUT, "B_contralateral_state", "spearman.tsv"))
  add("Spearman over the universe (N269B, and the thinned references' 95 % ranges):", "", md_table(sp, 3), "")
  nn <- rd(file.path(OUT, "B_contralateral_state", "nnls_weights.tsv"))
  add("NNLS weights (shares of the fitted mixture), N269B on C2B plus the six tumours, and each donor left out; each tumour against the others as the in-vivo spread:", "", md_table(nn, 3), "",
      sprintf("N269B's C2B share with IL69B left out: %s; the same fit on thinned IL69B: %s.", g("B", "nnls_N269B_C2B_share_without_IL69B"), g("B", "nnls_thinIL69B_C2B_share_null_95")),
      sprintf("Libraries sharing an i7 or i5 index with N269B (index hopping from them cannot be excluded): %s.", g("B", "libraries_sharing_an_index_with_N269B")), "",
      sprintf("Gene by gene (N269B against thinned IL69B): %s genes outside the thinned min-max, %s outside the six-tumour spread, **%s count as differences** (both). The rule's own null rate: %s.",
              g("B", "genes_outside_thinned_IL69B"), g("B", "genes_outside_six_tumour_spread"), g("B", "genes_counting_as_difference"), g("B", "null_rate_outside_thinned_minmax")),
      "PCR duplicates make N269B noisier than binomial thinning at D, so the gene count is an upper bound on real differences.", "",
      "**Can show:** whether the U251 RNA in N269B looks like U251 grown in brain or like U251 in culture, with the depth of N269B taken into account.",
      "**Cannot show:** that the cells migrated (carry-over of tumour tissue from the same animal gives the same in-vivo profile); anything about a population (one animal); hopped reads from index-sharing libraries would carry their own libraries' state into N269B. The claim stops at 'in-vivo U251 state, not culture-like' (or its negation); never 'migrated'.", "")
} else add("Not run: ", part_msg("B"), "")

## ---------------------------------------------------------------- C
add("## C. Do LITT recurrences carry more or less of published glioma infiltration programmes?", "")
if (c_ok) {
  su <- rd(file.path(OUT, "genesets", "sets_used.tsv"))
  add(sprintf("Gene sets from files only (`genesets/`, `genesets/SOURCES.tsv` holds each URL and md5). Pre-declared primary sets: %s. Missing primary sets: %s. Optional sets not available: %s.",
              g("C", "primary_sets"), g("C", "primary_sets_missing"), g("C", "optional_sets_missing")),
      sprintf("GSVA %s, fgsea %s, %s; msigdbr %s.", g("C", "GSVA_version"), g("C", "fgsea_version"), g("C", "R_version_stepC"), g("C", "msigdbr_version")), "",
      md_table(su[, intersect(c("set", "primary", "genes", "matched_graft_only", "source_file", "rule"), colnames(su))], 3), "")
  lk <- rd(file.path(OUT, "C_programmes", "leak_checks.tsv"))
  add("GSVA (Gaussian kernel on DESeq2 vst, all ten libraries scored), Welch t recurrent minus primary on the six tumours, rank of the true split among the 10; graft-only against graft+both, with the leak flags:", "",
      md_table(lk[, c("set", "primary_set", "matched_graft_only", "rank_graft_only", "t_graft_only", "rank_graft_plus_both", "t_graft_plus_both",
                      "high_in_controls_graft_only", "tracks_m_graft_only", "moves_only_with_both_bin", "leak_verdict")], 3), "",
      sprintf("Sets at rank 1 of 10 on graft-only: %s. With this many sets some reach rank 1 by chance.", g("C", "gsva_sets_at_rank1_graft_only")), "")
  fr <- rd(file.path(OUT, "C_programmes", "fgsea_split_ranks.tsv"))
  if (!is.null(fr)) add(sprintf("fgsea on each relabelling's Wald statistic (%s):", g("C", "fgsea_fits_graft_only")), "",
                        md_table(fr[fr$counts == "graft_only", c("set", "primary_set", "NES_true", "rank")], 3), "")
  bt <- rd(file.path(OUT, "C_programmes", "scores_beside_tumour_content.tsv"))
  if (!is.null(bt)) add("Primary-set GSVA scores (graft-only vst) beside tumour content. Graft % separates the arms completely; the DNA tumour fraction column is empty until the methylation integration fills it, and nothing is adjusted away:", "", md_table(bt, 3), "")
  add("**Can show:** whether any published infiltration programme separates the arms more than a random relabelling of the same six tumours does, and whether a separation seen on graft+both counts survives removal of the both bin.",
      "**Cannot show:** a primary-versus-recurrent difference at better than p = 0.10; whether a difference is LITT or tumour content (arm and graft % are confounded completely); anything about invasion in the animal (these are expression programmes in bulk tumour RNA).", "")
} else add("Not run (step C did not finish; see the step C log).", "")

## ---------------------------------------------------------------- D
add("## D. Host niche: is Mmp13 the one host change, and what tissue is the N269B piece?", "")
if (part_ok("D")) {
  nr <- rd(file.path(OUT, "D_host_niche", "mmp13_normalised.tsv"))
  add(sprintf("Mmp13, rank of the true split among the 10 host relabellings: %s.", g("D", "mmp13_true_split_rank_by_abs_wald")), "",
      "Mmp13 per unit of stroma and of myeloid markers (log2 CPM differences), recurrent minus primary, 10-split rank:", "", md_table(nr[, c("measure", "diff_recurrent_minus_primary", "welch_t_true", "rank", "N168B", "N269B")], 3), "",
      sprintf("Human MMP13 in the graft-only human counts: %s. Rat Mmp13 host counts: %s.", g("D", "human_MMP13_graft_only_counts"), g("D", "rat_Mmp13_counts")), "")
  ps <- rd(file.path(OUT, "D_host_niche", "marker_panels_summary.tsv"))
  add("Marker panels (rat host reads), median log2 CPM ratio per panel:", "", md_table(ps, 3), "",
      sprintf("Panel genes missing from the rat GTF: %s. Female-host check (rat chrY CPM): %s; Xist CPM: %s.",
              g("D", "panel_genes_missing_from_rat_gtf"), g("D", "rat_chrY_cpm_by_library"), g("D", "rat_Xist_cpm_by_library")), "",
      "**Can show:** whether the Mmp13 fall in recurrences outranks every relabelling, whether it survives division by stromal or myeloid content, whether the transcript is rat, and which compartments (meninges, white matter, grey matter, myeloid, reactive glia) N269B's piece is richer or poorer in than N168B's.",
      "**Cannot show:** where Mmp13 is expressed (needs ISH); Mmp13 stays 'exploratory; needs ISH'. The panel ratios describe the dissected piece, not the tumour's microenvironment, and N269B's host reads include any tumour-associated host cells that came with it.", "")
} else add("Not run: ", part_msg("D"), "")

add("## Not known", "",
    "- Whether the contralateral U251 cells migrated or were carried over; whether any N269B reads are index-hopped (the index pairs are shared with human-rich libraries).",
    "- Whether the Y share of U251 RNA is the same in contralateral cells as in the tumours (the dosage calibration assumes it).",
    "- The DNA tumour fraction per library (to be joined when the methylation integration runs).",
    "- The direction of Venkataramani 2022's two gene groups is read from the file's own legends (Table S2e), not from the paper's text.",
    "- Darmanis 2017, Hoelzinger 2005, Gill 2014 and Mariani 2001 lists could not be obtained as files (see genesets/SOURCES.tsv); they were not transcribed.", "")
im <- rd(file.path(OUT, "inputs", "images.tsv")); md <- rd(file.path(OUT, "inputs", "inputs_md5.tsv"))
add("## Inputs and images", "", if (!is.null(im)) md_table(im) else "images.tsv missing", "", if (!is.null(md)) md_table(md) else "inputs_md5.tsv missing", "")
writeLines(L, file.path(OUT, "REPORT.md.part"))
file.rename(file.path(OUT, "REPORT.md.part"), file.path(OUT, "REPORT.md"))
cat("REPORT.md written:", length(L), "lines\n")
