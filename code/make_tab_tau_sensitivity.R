# =============================================================================
# File: make_tab_tau_sensitivity.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Sensitivity of the fixed-K recovery to the prior scale tau, from the n = 80
# full-grid evaluation (informative prior, gamma = 0.6, mean causal variants
# recovered out of 15 over rho and replications).
#
# Usage (from code/):
#   Rscript make_tab_tau_sensitivity.R \
#     --in ../../server_rerun_fix_20260928/results_exact/results_bench_synth_informative_n80.rds \
#     --out ../tab_tau_sensitivity.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(`in` = "../results_bench_synth_informative_n80.rds",
            gamma = 0.6, out = "../tab_tau_sensitivity.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

PRETTY <- c(Naive_Wald = "NaW", Asymptotic_Wald = "AsW", Exact_Wald = "ExW",
            Asymptotic_Score = "AsS", Naive_Multiple = "NaM", Exact_Multiple = "ExM")
d <- readRDS(cfg$`in`)$all_results %>%
  filter(method %in% names(PRETTY), gamma == cfg$gamma) %>%
  group_by(method, tau) %>% summarise(M = mean(M), .groups = "drop") %>%
  pivot_wider(names_from = tau, values_from = M)
d <- d[match(names(PRETTY), d$method), ]
taus <- setdiff(names(d), "method")

rows <- vapply(seq_len(nrow(d)), function(i)
  paste0(PRETTY[d$method[i]], " & ",
         paste(sprintf("%.2f", as.numeric(d[i, taus])), collapse = " & "), " \\\\"), "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Sensitivity to the prior scale $\\tau$.} Mean number of causal ",
         "variants recovered (out of $15$) at fixed $K=15$ in the controlled AR(1) design at $n=80$ ",
         "(informative prior, $\\gamma=0.6$, averaged over $\\rho$ and replications), for ",
         "$\\tau \\in \\{0.2, 0.5, 1\\}$. The ordering of the methods is unchanged across $\\tau$, and ",
         "the exact scores of the joint methods are in addition stable in $\\tau$ on the real-LD ",
         "selection experiment of Table~\\ref{tab:exact_select}.}}"),
  "\\label{tab:tau_sensitivity}",
  "\\centering\\small",
  paste0("\\begin{tabular}{l", paste(rep("c", length(taus)), collapse = ""), "}"),
  "\\toprule",
  paste0("Method & ", paste(sprintf("$\\tau=%s$", taus), collapse = " & "), " \\\\"),
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
print(as.data.frame(d), row.names = FALSE)
