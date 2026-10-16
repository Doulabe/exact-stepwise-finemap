# =============================================================================
# File: make_geuvadis_cis.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Block-bootstrap confidence intervals for the LD-aware replication rates of
# the GEUVADIS application. The hundred cis windows were merged into 66
# disjoint genomic regions (geuvadis_regions.txt); the bootstrap resamples the
# regions, so that genes sharing overlapping windows (hence LD and possibly
# individuals' local haplotypes) enter or leave the resample together.
# Outputs a supplementary table (fixed-K rates with 95% CIs per prior, plus
# the per-gene TSS - Uniform difference) and prints the BIC analogue.
#
# Usage (from code/):
#   Rscript make_geuvadis_cis.R \
#     --results ../results_geuvadis_real_ld.rds --rds_dir ../geuvadis_rds \
#     --regions ../geuvadis_data/geuvadis_regions.txt \
#     --out ../tab_geuvadis_ci.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            rds_dir = "../geuvadis_rds",
            regions = "../geuvadis_data/geuvadis_regions.txt",
            cache   = "../geuvadis_data/geuvadis_gene_regions.csv",
            B = 2000L, seed = 2026L,
            out = "../tab_geuvadis_ci.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

res <- readRDS(cfg$results)
tab <- res$table
genes <- res$genes

# ---- gene -> merged region (cached: reading the gene RDS files is slow) -----
if (file.exists(cfg$cache)) {
  map <- read.csv(cfg$cache, stringsAsFactors = FALSE)
} else {
  reg <- read.table(cfg$regions, col.names = c("chr", "start", "end"))
  reg$region <- seq_len(nrow(reg))
  map <- do.call(rbind, lapply(genes, function(g) {
    gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
    r <- reg$region[reg$chr == gx$chr & reg$start <= gx$tss & gx$tss <= reg$end]
    data.frame(gene = g, chr = gx$chr, tss = gx$tss,
               region = if (length(r)) r[1L] else NA_integer_)
  }))
  stopifnot(!anyNA(map$region))
  write.csv(map, cfg$cache, row.names = FALSE)
}
cat(sprintf("%d genes in %d regions (sizes: %s)\n", nrow(map),
            length(unique(map$region)),
            paste(range(table(map$region)), collapse = "-")))

# ---- bootstrap over regions -------------------------------------------------
METHODS <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM", "SuSiE")
PRIORS  <- c("uniform", "polyfun", "tss", "discordant")
regions_of <- split(map$gene, map$region)
R <- length(regions_of)

boot_means <- function(d, value = "rep_rate_ld") {
  # d: data frame with columns gene and the value; returns B bootstrap means
  v <- setNames(d[[value]], d$gene)
  sapply(seq_len(cfg$B), function(b) {
    gs <- unlist(regions_of[sample.int(R, R, replace = TRUE)], use.names = FALSE)
    mean(v[gs], na.rm = TRUE)
  })
}

set.seed(cfg$seed)
ci_cell <- function(sel, pr, me) {
  d <- tab %>% filter(selection == sel, prior == pr, method == me) %>%
    select(gene, rep_rate_ld)
  bs <- boot_means(d)
  c(est = mean(d$rep_rate_ld, na.rm = TRUE),
    lo = unname(quantile(bs, .025, na.rm = TRUE)),
    hi = unname(quantile(bs, .975, na.rm = TRUE)))
}
ci_diff <- function(sel, me, pr1 = "tss", pr0 = "uniform") {
  d <- tab %>% filter(selection == sel, method == me, prior %in% c(pr1, pr0)) %>%
    select(gene, prior, rep_rate_ld) %>%
    pivot_wider(names_from = prior, values_from = rep_rate_ld)
  if (!all(c(pr1, pr0) %in% names(d)) || !nrow(d))
    return(c(est = NA_real_, lo = NA_real_, hi = NA_real_))
  d <- d %>% mutate(diff = .data[[pr1]] - .data[[pr0]])
  bs <- boot_means(d, "diff")
  c(est = mean(d$diff, na.rm = TRUE),
    lo = unname(quantile(bs, .025, na.rm = TRUE)),
    hi = unname(quantile(bs, .975, na.rm = TRUE)))
}

fmt <- function(x) sprintf("%.2f [%.2f, %.2f]", x["est"], x["lo"], x["hi"])
build <- function(sel) {
  rows <- vapply(METHODS, function(me) {
    cells <- vapply(PRIORS, function(pr) fmt(ci_cell(sel, pr, me)), "")
    dd <- fmt(ci_diff(sel, me))
    paste0(me, " & ", paste(cells, collapse = " & "), " & ", dd, " \\\\")
  }, "")
  rows
}

rows_fixed <- build("fixedK")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Block-bootstrap confidence intervals for the LD-aware ",
         "replication rates.} Mean LD-aware replication rate at fixed $K=5$ with $95\\%$ ",
         "percentile intervals from $", cfg$B, "$ bootstrap resamples of the $66$ disjoint ",
         "genomic regions formed by merging the hundred overlapping cis windows, so that ",
         "genes sharing local LD enter the resample together. Last column: per-gene ",
         "difference between the TSS-proximity and uniform priors. The BIC-stopping ",
         "analogue gives the same picture.}}"),
  "\\label{tab:geuvadis_ci}",
  "\\centering\\scriptsize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "\\textbf{Method} & Uniform & PolyFun & TSS & Discordant & TSS $-$ Uniform \\\\",
  "\\midrule",
  rows_fixed,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n\n=== BIC analogue (console only) ===\n")
for (me in METHODS) {
  if (!nrow(tab %>% filter(selection == "BIC", method == me))) next
  cells <- vapply(PRIORS, function(pr) fmt(ci_cell("BIC", pr, me)), "")
  cat(sprintf("%-6s %s | diff %s\n", me, paste(cells, collapse = "  "),
              fmt(ci_diff("BIC", me))))
}
