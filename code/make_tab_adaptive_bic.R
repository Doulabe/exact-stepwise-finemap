# =============================================================================
# File: make_tab_adaptive_bic.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Adaptive-K tables of the controlled AR(1) design (informative prior,
# n = 1000) with a PRESPECIFIED stopping rule for the in-house methods:
#   - main table: in-house methods stopped by BIC (the default rule of the
#     paper), external comparators by their native rules;
#   - supplementary table: the per-method best data-driven rule (the previous
#     presentation), kept as a sensitivity envelope.
#
# Usage (from code/):
#   Rscript make_tab_adaptive_bic.R \
#     --in ../../server_rerun_fix_20260928/results_exact/results_bench_synth_informative_adaptive.rds \
#     --out_main ../tab_synthetic_informative_adaptive.tex \
#     --out_best ../tab_synth_adaptive_best.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(`in` = "../results_bench_synth_informative_adaptive.rds",
            out_main = "../tab_synthetic_informative_adaptive.tex",
            out_best = "../tab_synth_adaptive_best.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}

PRETTY <- c(Naive_Wald = "NaW", Asymptotic_Wald = "AsW", Exact_Wald = "ExW",
            Asymptotic_Score = "AsS", Naive_Multiple = "NaM", Exact_Multiple = "ExM",
            SuSiE_unif = "SuSiE\\textsubscript{unif}", SuSiE_prior = "SuSiE\\textsubscript{prior}",
            LASSO_1se = "LASSO", BayesC = "BayesC")
IH <- names(PRETTY)[1:6]
d <- readRDS(cfg$`in`)$all_results %>% filter(method %in% names(PRETTY))
d$F1 <- 2 * d$M / (d$K + 15)          # F1 convention of the paper (empty = failure)

agg <- function(dd) dd %>% group_by(method, criterion, rho_intra) %>%
  summarise(F1 = mean(F1), K = mean(K), .groups = "drop")
A <- agg(d)
rhos <- sort(unique(A$rho_intra))

RULE <- c(R2 = "$R^2$", BIC = "BIC", eBIC = "eBIC",
          Post_0.5 = "$q>.5$", Post_0.95 = "$q>.95$", Post_0.99 = "$q>.99$",
          natural_adaptive = "native", PIP_0.5 = "PIP$>$.5")
mk_rows <- function(tab, show_rule = FALSE) {
  vapply(seq_len(nrow(tab)), function(i) {
    r <- tab[i, ]
    cells <- vapply(rhos, function(rr) {
      z <- A %>% filter(method == r$method, criterion == r$criterion, rho_intra == rr)
      rl <- if (r$criterion %in% names(RULE)) RULE[[r$criterion]] else gsub("_", "", r$criterion)
      sprintf("%.3f%s", z$F1,
              if (show_rule) sprintf(" [%s; %.1f]", rl, z$K) else
                sprintf(" (%.1f)", z$K))
    }, "")
    paste0(PRETTY[r$method], " & ", paste(cells, collapse = " & "), " \\\\")
  }, "")
}
header <- function(size = "\\small", colsep = "4pt") c(
  paste0("\\centering", size),
  paste0("\\setlength{\\tabcolsep}{", colsep, "}"),
  "\\begin{tabular}{lccc}", "\\toprule",
  paste0("\\textbf{Method} & ", paste(sprintf("$\\rho=%.2f$", rhos), collapse = " & "), " \\\\"),
  "\\midrule")

# ---- main: prespecified rules (in-house BIC, externals native) --------------
main_tab <- tibble(method = IH, criterion = "BIC")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{As Table~\\ref{tab:synth_informative_fixed} but in the adaptive-$K$ regime ",
         "with the prespecified BIC rule, the default of this paper, for the in-house methods ",
         "(mean selected size $\\hat K$ in parentheses). The external comparators are compared at ",
         "fixed $K$ in Table~\\ref{tab:synth_informative_fixed}, their native adaptive selections ",
         "being degenerate on this design. The per-method best data-driven ",
         "rule, a post hoc envelope, is reported as a sensitivity in ",
         "Table~\\ref{tab:synth_adaptive_best}.}}"),
  "\\label{tab:synth_informative_adaptive}",
  header(), mk_rows(main_tab), "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out_main)
cat("Wrote", cfg$out_main, "\n")

# ---- supplement: best rule per method ---------------------------------------
best_tab <- A %>% filter(rho_intra == rhos[2] | length(rhos) == 1L) %>%
  group_by(method) %>% slice_max(F1, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(method, criterion)
best_tab <- best_tab[match(names(PRETTY), best_tab$method), ] %>% filter(!is.na(method))
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{Adaptive-$K$ sensitivity envelope, giving for each method the best data-driven ",
         "stopping rule at $\\rho=0.95$ (rule and mean selected size in brackets). This post hoc ",
         "selection of the rule favors every method equally and is reported only as an upper ",
         "envelope of Table~\\ref{tab:synth_informative_adaptive}.}}"),
  "\\label{tab:synth_adaptive_best}",
  header("\\scriptsize", "2.5pt"), mk_rows(best_tab, show_rule = TRUE), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  cfg$out_best)
cat("Wrote", cfg$out_best, "\n")
print(as.data.frame(A %>% filter(criterion == "BIC", method %in% IH) %>%
  pivot_wider(names_from = rho_intra, values_from = c(F1, K))), row.names = FALSE)
