# =============================================================================
# File: fill_geuvadis_tables.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Generate the LaTeX tables of the GEUVADIS application from the benchmark
# outputs. Main text: tab_geuvadis_fixedK (LD-aware replication, fixed-K,
# nominal / permutation) and tab_geuvadis_adaptive (BIC stopping, same dual
# reading). Supplement: tab_geuvadis_sensitivity (exact-variant and per-window
# Bonferroni readings merged) and tab_geuvadis_semisynth (semi-synthetic arm).
# Run FROM code/. Tables written to the repo root (../).
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

ld    <- readRDS("../results_geuvadis_real_ld.rds")$table
bonf  <- readRDS("../results_geuvadis_real_ld_bonf.rds")$table
semis <- bind_rows(lapply(Sys.glob("../results_geuvadis_semisynth50_*.rds"),
                          function(f) readRDS(f)$summary))

METHODS <- c("NaW","AsW","ExW","AsS","NaM","ExM","SuSiE")
PRIORS  <- c(uniform = "Uniform", polyfun = "PolyFun", tss = "TSS",
             discordant = "Discordant")
fmt <- function(x, d = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", d, "f"), x))

wide_rate <- function(tab, sel, value_col) {
  tab %>% filter(selection == sel, method %in% METHODS, prior %in% names(PRIORS)) %>%
    group_by(method, prior) %>%
    summarise(v = mean(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = prior, values_from = v) %>%
    slice(match(METHODS, method))
}

rate_table <- function(tab, sel, value_col, caption, label, file, khat = NULL) {
  w <- wide_rate(tab, sel, value_col)
  rows <- vapply(seq_len(nrow(w)), function(i) {
    cells <- vapply(names(PRIORS), function(p) {
      s <- fmt(w[[p]][i])
      if (!is.null(khat)) s <- sprintf("%s (%.1f)", s, khat[[p]][i])
      s
    }, "")
    paste0(w$method[i], " & ", paste(cells, collapse = " & "), " \\\\")
  }, "")
  writeLines(c(
    "\\begin{table}[htbp]",
    paste0("\\caption{", caption, "}"),
    paste0("\\label{", label, "}"),
    "\\centering\\small",
    "\\begin{tabular}{lcccc}", "\\toprule",
    paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
    rows,
    "\\bottomrule", "\\end{tabular}", "\\end{table}"), file)
  cat("Wrote", file, "\n")
}

# ---- 1. main: fixed-K, LD-aware, nominal / permutation (+ LASSO row) --------
permt <- if (file.exists("../results_geuvadis_real_ld_perm.rds"))
  readRDS("../results_geuvadis_real_ld_perm.rds")$table else NULL
w_n <- wide_rate(ld, "fixedK", "rep_rate_ld")
w_p <- if (!is.null(permt)) wide_rate(permt, "fixedK", "rep_rate_ld") else NULL
rows1 <- vapply(seq_len(nrow(w_n)), function(i) {
  cells <- vapply(names(PRIORS), function(p) {
    s <- fmt(w_n[[p]][i])
    if (!is.null(w_p)) s <- paste0(s, " / ", fmt(w_p[[p]][i]))
    s
  }, "")
  paste0(w_n$method[i], " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
lstat <- ld %>% filter(method == "LASSO") %>%
  summarise(v = mean(rep_rate_ld, na.rm = TRUE), K = mean(K, na.rm = TRUE))
lstat_p <- if (!is.null(permt)) permt %>% filter(method == "LASSO") %>%
  summarise(v = mean(rep_rate_ld, na.rm = TRUE)) else NULL
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\textbf{Cross-population replication of fixed-$K$ selections in the GEUVADIS application.} ",
         "\\rev{Each cell reports the mean LD-aware replication rate over the 100 genes under the nominal ",
         "replication test, followed by the rate under a per-window permutation threshold (the $95$th ",
         "percentile of the maximal $|t|$ statistic over $500$ permutations of the residualized YRI ",
         "phenotype, a window-wide family-wise criterion that adapts to the local LD):} the fraction of ",
         "the $K=5$ selected variants in replication-panel LD ($|r| \\geq 0.7$) with a variant that ",
         "replicates in the YRI set \\rev{(sign-consistent, significant at the respective threshold)}, by ",
         "prior specification. \\rev{The nominal reading is the more sensitive and the permutation reading ",
         "the more specific, and they give the same ordering.} The last rows report the prior-free LASSO, ",
         "which ignores the prior arms and retains $\\hat K$ variants at the one-standard-error penalty",
         "\\rev{, and the native conditional cis scan of QTLtools \\parencite{delaneau2017}, run on the ",
         "same windows, phenotype and covariates ($1{,}000$-permutation gene-level calibration at ",
         "$\\alpha = 0.05$, forward-backward conditional pass, $\\hat K$ being its number of declared ",
         "independent signals)}."),
  "}",
  "\\label{tab:geuvadis_fixedK}",
  "\\centering\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
  rows1,
  "\\midrule",
  sprintf("LASSO (prior-free) & \\multicolumn{4}{c}{%s%s \\quad ($\\hat K = %.1f$)} \\\\",
          fmt(lstat$v),
          if (!is.null(lstat_p)) paste0(" / ", fmt(lstat_p$v)) else "", lstat$K),
  if (file.exists("../results_qtltools_eval.rds")) {
    q <- readRDS("../results_qtltools_eval.rds")$table
    sprintf("\\rev{QTLtools (prior-free)} & \\multicolumn{4}{c}{%s / %s \\quad ($\\hat K = %.1f$)} \\\\",
            fmt(mean(q$rep_rate_ld, na.rm = TRUE)),
            fmt(mean(q$rep_rate_ld_perm, na.rm = TRUE)), mean(q$K))
  } else character(0),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "../tab_geuvadis_fixedK.tex")
cat("Wrote ../tab_geuvadis_fixedK.tex \n")

# ---- 2. main: adaptive BIC, nominal / permutation (+ mean selected size) -----
# The permutation reading is merged into the main adaptive table (second value
# of each cell); no separate perm_bic supplement table is emitted any more.
w_k     <- wide_rate(ld, "BIC", "K") %>% filter(method != "SuSiE")
w_n_bic <- wide_rate(ld %>% filter(method != "SuSiE"), "BIC", "rep_rate_ld")
w_p_bic <- if (!is.null(permt))
  wide_rate(permt %>% filter(method != "SuSiE"), "BIC", "rep_rate_ld") else NULL
rows2 <- vapply(seq_len(nrow(w_n_bic)), function(i) {
  cells <- vapply(names(PRIORS), function(p) {
    s <- fmt(w_n_bic[[p]][i])
    if (!is.null(w_p_bic)) s <- paste0(s, " / ", fmt(w_p_bic[[p]][i]))
    sprintf("%s (%.1f)", s, w_k[[p]][i])
  }, "")
  paste0(w_n_bic$method[i], " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\textbf{Adaptive BIC-stopped selections in the GEUVADIS application.} ",
         "\\rev{Each cell reports the mean LD-aware replication rate over the 100 genes under the ",
         "nominal replication test, followed by the rate under the per-window permutation threshold ",
         "of Table~\\ref{tab:geuvadis_fixedK}} (mean selected size $\\hat K$ in parentheses, ",
         "iteration cap 30), by prior specification. The joint methods stop at compact selections ",
         "\\rev{($\\hat K \\approx 1.5$--$2$) whose replication rates exceed, under both readings, ",
         "the field-standard QTLtools conditional scan, the natural comparator at this operating ",
         "point ($\\hat K = 1.5$, $0.79$ nominal and $0.52$ permutation-controlled, ",
         "Table~\\ref{tab:geuvadis_fixedK}), and the fixed-$K=5$ SuSiE selection of the same table}.}"),
  "\\label{tab:geuvadis_adaptive}",
  "\\centering\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
  rows2,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "../tab_geuvadis_adaptive.tex")
cat("Wrote ../tab_geuvadis_adaptive.tex \n")

# ---- 3+4. supplement: exact-variant and Bonferroni sensitivities (merged) ----
w_ex <- wide_rate(ld, "fixedK", "rep_rate")
w_bf <- wide_rate(bonf, "fixedK", "rep_rate_ld")
rows34 <- vapply(seq_len(nrow(w_ex)), function(i) {
  cells <- vapply(names(PRIORS), function(p)
    sprintf("%s / %s", fmt(w_ex[[p]][i]), fmt(w_bf[[p]][i])), "")
  paste0(w_ex$method[i], " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Replication sensitivities of the fixed-$K=5$ selections.} ",
         "Each cell reports, by prior specification, the exact-variant replication rate (crediting ",
         "a selected variant only when it replicates itself, no LD credit), followed by the ",
         "LD-aware rate under a per-window Bonferroni replication threshold $\\alpha/m$, at which ",
         "the $n=89$ replication set is underpowered. Both variants preserve the ordering of ",
         "priors and methods of Table~\\ref{tab:geuvadis_fixedK} of the main text.}}"),
  "\\label{tab:geuvadis_sensitivity}",
  "\\centering\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
  rows34, "\\bottomrule", "\\end{tabular}", "\\end{table}"), "../tab_geuvadis_sensitivity.tex")
cat("Wrote ../tab_geuvadis_sensitivity.tex \n")

# ---- 5. supplement: semi-synthetic benchmark --------------------------------
sm <- semis %>% filter(method %in% METHODS, prior %in% c("polyfun", "uniform")) %>%
  group_by(n, prior, method) %>% summarise(M_clust = mean(M_clust), .groups = "drop") %>%
  pivot_wider(names_from = c(prior, n), values_from = M_clust)
rows <- vapply(METHODS, function(m) {
  r <- sm %>% filter(method == m)
  paste0(m, " & ",
         paste(vapply(c("polyfun_80","uniform_80","polyfun_373","uniform_373"),
                      function(cc) fmt(r[[cc]], 2), ""), collapse = " & "), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\textbf{Semi-synthetic evaluation on real GEUVADIS cis windows.} ",
         "Mean number of causal LD clusters tagged (out of $K_{true}=5$, $|r| \\geq 0.7$) over thirty gene ",
         "windows spanning the eGene list ($B=50$ replications per window and sample size), by subsample ",
         "size $n$ and prior arm. True causal variants are drawn with probability proportional to the ",
         "PolyFun per-SNP weights, so the PolyFun arm is the truth-aligned prior in this design and the ",
         "uniform arm is neutral.}"),
  "\\label{tab:geuvadis_semisynth}",
  "\\centering\\small",
  "\\begin{tabular}{lcccc}", "\\toprule",
  " & \\multicolumn{2}{c}{$n=80$} & \\multicolumn{2}{c}{$n=373$} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Method & PolyFun & Uniform & PolyFun & Uniform \\\\", "\\midrule",
  rows,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "../tab_geuvadis_semisynth.tex")
cat("Wrote ../tab_geuvadis_semisynth.tex\n")

# ---- 6. main: computational cost (per GEUVADIS window) ----------------------
tm <- ld %>% filter(method %in% c(METHODS, "LASSO")) %>%
  group_by(method, selection) %>% summarise(s = mean(runtime, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = selection, values_from = s)
rows <- vapply(c(METHODS, "LASSO"), function(m) {
  r <- tm %>% filter(method == m)
  paste0(m, " & ", fmt(r$fixedK, 2), " & ",
         if (!is.null(r$BIC) && !is.na(r$BIC)) fmt(r$BIC, 2) else "--", " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Mean wall-clock time per gene window (s)} in the GEUVADIS application ",
         "(median $m = 4{,}424$ candidate variants, $n = 373$, single-threaded ",
         "within each window), for the fixed-$K=5$ and the adaptive BIC-stopped runs. SuSiE and LASSO are ",
         "run at fixed $K$ only.}}"),
  "\\label{tab:timing}",
  "\\centering\\small",
  "\\begin{tabular}{lcc}", "\\toprule",
  "Method & Fixed-$K$ & Adaptive (BIC) \\\\", "\\midrule",
  rows,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), "../tab_timing_geuvadis.tex")
cat("Wrote ../tab_timing_geuvadis.tex\n")

# ---- headline numbers for the manuscript text -------------------------------
cat("\n--- headline numbers ---\n")
g <- ld %>% filter(selection == "fixedK", method %in% METHODS[1:6]) %>% group_by(prior) %>%
  summarise(ld = round(mean(rep_rate_ld, na.rm = TRUE), 3))
print(as.data.frame(g), row.names = FALSE)
h <- ld %>% filter((method == "SuSiE" & selection == "fixedK") | (method %in% c("ExM","NaM") & selection == "BIC")) %>%
  group_by(method, selection) %>%
  summarise(K = round(mean(K), 1), ld = round(mean(rep_rate_ld, na.rm = TRUE), 3),
            R2 = round(mean(R2_yri, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(h), row.names = FALSE)
p <- semis %>% filter(prior %in% c("uniform","polyfun","tss")) %>% group_by(n, method) %>%
  summarise(M_clust = mean(M_clust), .groups = "drop") %>% pivot_wider(names_from = method, values_from = M_clust) %>%
  mutate(dExW = round(ExW - AsW, 2), dExM = round(ExM - NaM, 2),
         best_ih = round(pmax(NaW,AsW,ExW,AsS,NaM,ExM), 2), SuSiE = round(SuSiE, 2))
print(as.data.frame(p %>% select(n, dExW, dExM, best_ih, SuSiE)), row.names = FALSE)
sg <- semis %>% filter(n == 80, method %in% METHODS[1:6]) %>% group_by(prior) %>%
  summarise(M_clust = round(mean(M_clust), 2))
print(as.data.frame(sg), row.names = FALSE)
