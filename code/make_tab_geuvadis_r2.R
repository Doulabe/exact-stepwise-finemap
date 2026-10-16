# =============================================================================
# File: make_tab_geuvadis_r2.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Out-of-sample adjusted R2 of the selected sets in the YRI replication panel,
# the second replication metric announced in Section 4.3: fixed-K block and
# BIC block (with the mean selected size), by prior specification.
#
# Usage (from code/):
#   Rscript make_tab_geuvadis_r2.R \
#     --results ../results_geuvadis_real_ld.rds --out ../tab_geuvadis_r2.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            out = "../tab_geuvadis_r2.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}

tab <- readRDS(cfg$results)$table
METHODS <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM", "SuSiE")
PRIORS  <- c(uniform = "Uniform", polyfun = "PolyFun", tss = "TSS",
             discordant = "Discordant")
fmt <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))

block <- function(sel, with_k = FALSE) {
  w <- tab %>% filter(selection == sel, method %in% METHODS, prior %in% names(PRIORS)) %>%
    group_by(method, prior) %>%
    summarise(v = mean(R2_yri, na.rm = TRUE), K = mean(K, na.rm = TRUE), .groups = "drop")
  vapply(METHODS[METHODS %in% w$method], function(m) {
    cells <- vapply(names(PRIORS), function(p) {
      z <- w %>% filter(method == m, prior == p)
      if (!nrow(z)) return("--")
      if (with_k) sprintf("%s (%.1f)", fmt(z$v), z$K) else fmt(z$v)
    }, "")
    paste0(m, " & ", paste(cells, collapse = " & "), " \\\\")
  }, "")
}
lasso <- tab %>% filter(method == "LASSO") %>%
  summarise(v = mean(R2_yri, na.rm = TRUE), K = mean(K, na.rm = TRUE))

writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Out-of-sample adjusted $R^2$ in the YRI replication panel.} ",
         "Mean adjusted $R^2$ of the selected set over the $100$ genes, by prior specification, at ",
         "fixed $K=5$ (top block) and under BIC stopping (bottom block, mean selected size $\\hat K$ ",
         "in parentheses). The larger single-SNP selections under BIC buy a larger out-of-sample ",
         "$R^2$ at a lower per-variant replication rate (Table~\\ref{tab:geuvadis_adaptive}). ",
         "The prior-free LASSO attains $", fmt(lasso$v), "$ with $\\hat K=", sprintf("%.1f", lasso$K), "$.}}"),
  "\\label{tab:geuvadis_r2}",
  "\\centering\\small",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"),
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{Fixed $K=5$}} \\\\",
  block("fixedK"),
  "\\midrule",
  "\\multicolumn{5}{l}{\\emph{BIC stopping}} \\\\",
  block("BIC", with_k = TRUE),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
