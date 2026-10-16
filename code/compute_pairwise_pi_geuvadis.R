# =============================================================================
# File: compute_pairwise_pi_geuvadis.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Pairwise probability of
# superiority between methods ACROSS THE GEUVADIS GENES, using the per-gene
# LD-aware replication rate as the per-replicate performance measure:
#     pi_hat(A > B) = (1/G) sum_g 1{ r_A(g) > r_B(g) }
#                   + (1/(2G)) sum_g 1{ r_A(g) = r_B(g) }
# (same tie convention as the simulation version, Section 2.4). Ties are
# frequent here because the fixed-K rates live on a discrete grid.
#
# Usage (from code/):
#   Rscript compute_pairwise_pi_geuvadis.R --prior tss \
#       --out ../tab_pi_geuvadis.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(results  = "../results_geuvadis_real_ld.rds",
            prior    = "tss",
            selection= "fixedK",
            metric   = "rep_rate_ld",
            digits   = 2L,
            out      = "../tab_pi_geuvadis.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

METHODS <- c("NaW","AsW","ExW","AsS","NaM","ExM","SuSiE","LASSO")
tab <- readRDS(cfg$results)$table %>%
  filter(selection == cfg$selection,
         (prior == cfg$prior & method != "LASSO") | method == "LASSO") %>%
  select(gene, method, r = all_of(cfg$metric))
W <- tab %>% pivot_wider(names_from = method, values_from = r)
W <- W[complete.cases(W[, METHODS]), ]
G <- nrow(W)
cat(sprintf("prior=%s, metric=%s: %d genes with complete method sets\n",
            cfg$prior, cfg$metric, G))

P <- matrix(NA_real_, length(METHODS), length(METHODS),
            dimnames = list(METHODS, METHODS))
for (A in METHODS) for (B in METHODS) if (A != B)
  P[A, B] <- mean(W[[A]] > W[[B]]) + 0.5 * mean(W[[A]] == W[[B]])
pbar <- rowMeans(P, na.rm = TRUE)
ord  <- order(-pbar)

fmt <- function(x) ifelse(is.na(x), "--", sprintf(paste0("%.", cfg$digits, "f"), x))
rows <- vapply(METHODS[ord], function(A) {
  paste0(A, " & ", paste(fmt(P[A, METHODS[ord]]), collapse = " & "),
         " & ", fmt(pbar[A]), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Pairwise probability of superiority on the real GEUVADIS genes.} ",
         "Entries are $\\hat\\pi(A > B)$, the proportion of the ", G, " genes on which the row method ",
         "attains a strictly higher LD-aware replication rate than the column method (ties counted ",
         "$1/2$), fixed $K=5$, ",
         if (cfg$prior == "tss") "TSS-proximity prior" else paste0(cfg$prior, " prior"),
         " (LASSO is prior-free). $\\bar\\pi$: row mean. Rows and columns sorted by $\\bar\\pi$.}}"),
  "\\label{tab:pi_geuvadis}",
  "\\centering\\small",
  paste0("\\begin{tabular}{l", paste(rep("c", length(METHODS) + 1L), collapse = ""), "}"),
  "\\toprule",
  paste0(" & ", paste(METHODS[ord], collapse = " & "), " & $\\bar\\pi$ \\\\"),
  "\\midrule",
  rows,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n\nMatrix (sorted):\n")
print(round(P[ord, ord], 2))
cat("\nRow means:\n"); print(round(sort(pbar, decreasing = TRUE), 3))
