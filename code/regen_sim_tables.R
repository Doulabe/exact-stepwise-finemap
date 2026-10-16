# =============================================================================
# regen_sim_tables.R  (Article 1, SiM version)
# Regenerates the SiM tables and the figure that depend on AsS, NaM and ExM
# from a directory of result files, into a staging directory (the manuscript
# directory is never written here).
#
#   main text : tab_geuvadis_fixedK, tab_geuvadis_adaptive,
#               tab_synthetic_informative, Figure_premium_vs_n.pdf
#   supplement: tab_pi_geuvadis, tab_geuvadis_exact_snp, tab_geuvadis_bonferroni,
#               tab_locus_sim_informative_fixed, tab_locus_sim_informative_adaptive,
#               tab_carma_prior_sensitivity
#   headline_numbers.txt: the numbers quoted in the text.
#
# The GEUVADIS tables and the figure are produced by the repository scripts
# (fill_geuvadis_tables.R, compute_pairwise_pi_geuvadis.R, figure_premium_vs_n.R)
# run in a temporary workspace; the other tables reproduce the layout of the
# hand-built SiM tables. Files absent from --in are skipped.
#
# Usage (from code/; point --in/--orig at a results directory, e.g.
#        ../../server_rerun_fix_20260928/results_exact):
#   Rscript regen_sim_tables.R --in results_orig  --out staging_orig
#   Rscript regen_sim_tables.R --in results_fixed --out staging_fixed
#   --semisynth DIR : semi-synthetic GEUVADIS files (default: the original ones
#                     in the repository root, so tab_geuvadis_semisynth is only
#                     copied when DIR is given explicitly)
#   --orig DIR      : original result files; the summary rows of the methods
#                     that were not rerun are taken from them. The merged
#                     summaries hold the same means up to 1e-15 (summation
#                     order), which is enough to flip the last printed digit of
#                     exact halves such as 0.965.
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(`in` = "results_fixed", out = "staging_fixed", semisynth = "", orig = "",
            repo = "..", prior_only = "../results_prior_only.rds")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}
IN   <- normalizePath(cfg$`in`, mustWork = TRUE)
REPO <- normalizePath(cfg$repo, mustWork = TRUE)
dir.create(cfg$out, showWarnings = FALSE, recursive = TRUE)
OUT  <- normalizePath(cfg$out)
f_in <- function(x) file.path(IN, x)
has  <- function(x) file.exists(f_in(x))
HEAD <- file.path(OUT, "headline_numbers.txt")
if (file.exists(HEAD)) invisible(file.remove(HEAD))
say <- function(...) cat(paste(c(...), collapse = "\n"), "\n", file = HEAD, append = TRUE, sep = "")

fmt  <- function(x, d) ifelse(is.na(x) | !is.finite(x), "{--}", sprintf(paste0("%.", d, "f"), x))
pretty <- c(Naive_Wald = "NaW", Asymptotic_Wald = "AsW", Exact_Wald = "ExW",
            Asymptotic_Score = "AsS", Naive_Multiple = "NaM", Exact_Multiple = "ExM",
            SuSiE_unif = "SuSiE\\textsubscript{unif}", SuSiE_prior = "SuSiE\\textsubscript{prior}",
            varbvs = "varbvs", LASSO_topK = "LASSO\\textsubscript{topK}",
            LASSO_1se = "LASSO\\textsubscript{1se}", BayesC = "BayesC",
            CARMA = "CARMA", FINEMAP = "FINEMAP")
INHOUSE <- c("Naive_Wald", "Asymptotic_Wald", "Exact_Wald",
             "Asymptotic_Score", "Naive_Multiple", "Exact_Multiple")
RERUN <- c("Asymptotic_Score", "Naive_Multiple", "Exact_Multiple")
ORIG  <- if (nzchar(cfg$orig)) normalizePath(cfg$orig, mustWork = TRUE) else ""
# F1 is recomputed from the per-replication results as 2M/(K+q), q the number of
# causal variants. It equals 0 when no causal variant is selected, whereas the
# stored F1 is NA in that case and the stored means simply drop those
# replications, which inflates the F1 of methods that often select nothing.
`%||%` <- function(a, b) if (is.null(a)) b else a
fixed_f1_summary <- function(obj) {
  s <- as.data.frame(obj$summary)
  r <- obj$all_results %||% obj$raw
  q <- obj$config$n_causal %||% obj$config$q_locus
  if (is.null(r) || is.null(q)) stop("cannot recompute F1: per-replication results or q missing")
  r <- as.data.frame(r)
  for (k in intersect(c("n", "rho_intra", "gamma", "tau", "M", "K"), names(r))) r[[k]] <- as.numeric(r[[k]])
  keys <- intersect(c("n", "rho_intra", "gamma", "tau", "method", "criterion"), names(s))
  f <- r %>% mutate(F1c = ifelse(M > 0, 2 * M / (K + q), 0)) %>%
    group_by(across(all_of(keys))) %>% summarise(F1c = mean(F1c), .groups = "drop")
  s <- s %>% left_join(f, by = keys)
  stopifnot(!any(is.na(s$F1c)))
  s %>% mutate(F1 = F1c) %>% select(-F1c)
}
read_summary <- function(f) {
  s <- fixed_f1_summary(readRDS(f_in(f)))
  if (nzchar(ORIG) && file.exists(file.path(ORIG, f))) {
    so <- fixed_f1_summary(readRDS(file.path(ORIG, f)))
    s <- bind_rows(so[!(so$method %in% RERUN), ], s[s$method %in% RERUN, ])
  }
  s
}

# ---- 1. Synthetic AR(1), informative, fixed K (main text) -------------------
if (has("results_bench_synth_informative.rds")) {
  sm <- read_summary("results_bench_synth_informative.rds") %>%
    filter(n == 1000, criterion == "K_fixed", gamma == 0.6, tau == 0.2)
  rhos <- c(0.8, 0.95, 0.98)
  meths <- c(INHOUSE, "SuSiE_unif", "SuSiE_prior", "LASSO_topK", "LASSO_1se", "BayesC")
  rows <- character(0)
  # Prior-only reference row (rank by p_j, ties by marginal |r|): under this
  # oracle prior it is essentially perfect by construction, so the informative
  # column reads as prior tracking (compute_prior_only.R, same data sets).
  if (file.exists(cfg$prior_only)) {
    po <- readRDS(cfg$prior_only)$all_results %>%
      filter(file == "results_bench_synth_informative.rds", n == 1000)
    v <- vapply(rhos, function(r) {
      x <- po$F1[abs(po$rho_intra - r) < 1e-6]
      if (length(x)) mean(x) else NA_real_ }, 0)
    rows <- c(rows, sprintf("\\rev{Prior-only ranking} & %s \\\\",
                            paste(fmt(v, 3), collapse = " & ")), "\\midrule")
  }
  for (m in meths) {
    v <- vapply(rhos, function(r) {
      x <- sm$F1[sm$method == m & abs(sm$rho_intra - r) < 1e-6]
      if (length(x) == 1L) x else NA_real_ }, 0)
    if (all(is.na(v))) next
    rows <- c(rows, sprintf("%s & %s \\\\", pretty[m], paste(fmt(v, 3), collapse = " & ")))
  }
  writeLines(c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{\\textbf{Synthetic AR(1) LD sensitivity}: F1 across $\\rho\\in\\{0.80,0.95,0.98\\}$ on a synthetic design ($m=5{,}000$, $100$ AR(1) blocks of $50$), \\textit{informative} prior, fixed-$K$, $\\gamma=0.6$, $\\tau_j=0.2$. varbvs omitted \\rev{because its runtime was prohibitive at $m=5{,}000$ in this simulation grid}. \\rev{The prior-only reference ranks by $p_j$ with ties broken by the marginal $|r|$ on the same data sets. This oracle prior identifies the causal set by construction, so the column reads as prior tracking, and the margin of a method over the reference is the contribution of the data.}}",
    "\\label{tab:synth_informative_fixed}",
    "\\setlength{\\tabcolsep}{4pt}", "\\renewcommand{\\arraystretch}{1.0}", "\\small",
    "\\begin{tabular}{lccc}", "\\toprule",
    "\\textbf{Method} & \\multicolumn{3}{c}{$n = 1000$} \\\\", "\\cmidrule(lr){2-4}",
    " & $\\rho=0.80$ & $\\rho=0.95$ & $\\rho=0.98$ \\\\", "\\midrule",
    rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(OUT, "tab_synthetic_informative.tex"))
  say("== Synthetic AR(1), informative, fixed K, n = 1000 (F1; M = causal SNPs recovered /15)")
  w <- sm %>% filter(method %in% meths) %>% select(method, rho_intra, F1, mean_M) %>%
    arrange(match(method, meths), rho_intra)
  capture.output(print(as.data.frame(w), digits = 4, row.names = FALSE)) %>% say()
}

if (has("results_bench_synth_informative_adaptive.rds")) {
  # The adaptive tables have a single source of truth: make_tab_adaptive_bic.R
  # writes the main table (prespecified BIC rule, in-house methods) and the
  # supplementary per-method best-rule envelope (tab_synth_adaptive_best).
  system(paste("Rscript make_tab_adaptive_bic.R --in",
               f_in("results_bench_synth_informative_adaptive.rds"),
               "--out_main", file.path(OUT, "tab_synthetic_informative_adaptive.tex"),
               "--out_best", file.path(OUT, "tab_synth_adaptive_best.tex")))
}

# ---- 2. Single-locus tables (supplement) --------------------------------
locus_row_cells <- function(r) c(fmt(r$mean_M, 2), fmt(r$mean_K, 1), fmt(r$precision, 3),
                                 fmt(r$recall, 3), fmt(r$F1, 3), fmt(r$FDR, 3))
bold <- function(cells) ifelse(cells == "{--}", cells, paste0("\\textbf{", cells, "}"))

if (has("results_locus_sim_informative_fixed.rds")) {
  s <- read_summary("results_locus_sim_informative_fixed.rds") %>%
    filter(n == 1000, criterion == "K_fixed")
  meths <- c(INHOUSE, "SuSiE_unif", "SuSiE_prior", "varbvs", "LASSO_topK", "LASSO_1se",
             "CARMA", "FINEMAP")
  meths <- meths[meths %in% s$method]
  best <- meths[which.max(vapply(meths, function(m) s$F1[s$method == m], 0))]
  body <- vapply(seq_along(meths), function(i) {
    m <- meths[i]; cells <- locus_row_cells(s[s$method == m, ])
    if (m == best) cells <- bold(cells)
    lead <- if (i == 1L) sprintf("\\multirow{%d}{*}{$1000$}", length(meths)) else ""
    paste0(lead, " & ", pretty[m], " & ", paste(cells, collapse = " & "), " \\\\")
  }, "")
  writeLines(c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{\\textbf{\\rev{Single-locus} assessment vs.\\ external comparators}, \\textit{informative} prior, fixed-$K$: single AR(1) block ($m_{\\text{locus}}=200$, $\\rho=0.9$), $\\gamma=0.6$, $\\tau_j=0.2$, top-$K=3$. \\textbf{Bold} best F1.}",
    "\\label{tab:locus_sim_informative_fixed}",
    "\\setlength{\\tabcolsep}{3.5pt}", "\\renewcommand{\\arraystretch}{1.0}", "\\scriptsize",
    "\\begin{tabular}{c l c c c c c c}", "\\toprule",
    "$n$ & \\textbf{Method} & $\\hat M$ & $\\hat K$ & Precision & Recall & F1 & FDR \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(OUT, "tab_locus_sim_informative_fixed.tex"))
}

if (has("results_locus_sim_informative_adaptive.rds")) {
  s <- read_summary("results_locus_sim_informative_adaptive.rds") %>% filter(n == 1000)
  crit_ih  <- c(R2 = "$R^2$", BIC = "$BIC$", eBIC = "$eBIC$", Post_0.5 = "$Post(0.5)$",
                Post_0.95 = "$Post(0.95)$", Post_0.99 = "$Post(0.99)$")
  crit_pip <- c(PIP_0.5 = "$PIP{>}0.5$", PIP_0.95 = "$PIP{>}0.95$", PIP_0.99 = "$PIP{>}0.99$")
  blocks <- list()
  for (m in INHOUSE[INHOUSE %in% s$method])
    blocks[[length(blocks) + 1L]] <- list(m = m, crit = crit_ih, grouped = TRUE)
  for (m in c("SuSiE_unif", "SuSiE_prior", "varbvs", "LASSO_1se", "CARMA", "FINEMAP")) {
    if (!m %in% s$method) next
    cr <- if (m %in% c("varbvs", "CARMA", "FINEMAP")) crit_pip else c(natural_adaptive = "$Adaptive$")
    blocks[[length(blocks) + 1L]] <- list(m = m, crit = cr, grouped = FALSE)
  }
  rows_all <- bind_rows(lapply(blocks, function(b)
    tibble(method = b$m, criterion = names(b$crit), label = unname(b$crit), grouped = b$grouped)))
  rows_all <- rows_all %>% left_join(s, by = c("method", "criterion"))
  ib <- which.max(rows_all$F1)
  n_rows <- nrow(rows_all)
  body <- character(0); first <- TRUE
  for (b in blocks) {
    idx <- which(rows_all$method == b$m & rows_all$criterion %in% names(b$crit))
    for (j in seq_along(idx)) {
      r <- rows_all[idx[j], ]; cells <- locus_row_cells(r)
      if (idx[j] == ib) cells <- bold(cells)
      lead <- if (first) sprintf("\\multirow{%d}{*}{$1000$}", n_rows) else ""
      first <- FALSE
      mcol <- if (b$grouped) { if (j == 1L) sprintf("\\multirow{%d}{*}{%s}", length(idx), pretty[b$m]) else "" }
              else pretty[b$m]
      sep <- if (b$grouped && j > 1L) " &  & " else paste0(" & ", mcol, " & ")
      body <- c(body, paste0(lead, sep, r$label, " & ", paste(cells, collapse = " & "), " \\\\"))
    }
    if (b$grouped) body <- c(body, "\\cmidrule(lr){2-9}")
  }
  writeLines(c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{As Table~\\ref{tab:locus_sim_informative_fixed} but adaptive-$K$ ($\\hat K$ by the \\emph{Criterion} column), while CARMA/FINEMAP run prior-free, varbvs via \\texttt{logodds}. \\rev{Precision and FDR are averaged over the replications with a non-empty selection.}}",
    "\\label{tab:locus_sim_informative_adaptive}",
    "\\setlength{\\tabcolsep}{3.5pt}", "\\renewcommand{\\arraystretch}{1.0}", "\\scriptsize",
    "\\begin{tabular}{c l l c c c c c c}", "\\toprule",
    "$n$ & \\textbf{Method} & \\textbf{Criterion} & $\\hat M$ & $\\hat K$ & Precision & Recall & F1 & FDR \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(OUT, "tab_locus_sim_informative_adaptive.tex"))
}

# CARMA sensitivity: n = 1000, best F1 over the stopping criteria of each method.
carma_files <- c("informative_fixed", "informative_adaptive", "uniform_fixed", "uniform_adaptive")
if (all(has(sprintf("results_locus_sim_%s.rds", carma_files)))) {
  best_f1 <- function(sc, md, m) {
    s <- read_summary(sprintf("results_locus_sim_%s_%s.rds", sc, md))
    x <- s$F1[s$n == 1000 & s$method == m]
    if (length(x)) max(x, na.rm = TRUE) else NA_real_
  }
  body <- character(0)
  say("== CARMA sensitivity (n = 1000, best criterion): ExM, SuSiE_prior, CARMA, CARMA_prior")
  for (sc in c("informative", "uniform")) for (md in c("fixed", "adaptive")) {
    v <- vapply(c("Exact_Multiple", "SuSiE_prior", "CARMA", "CARMA_prior"),
                function(m) best_f1(sc, md, m), 0)
    body <- c(body, sprintf("%-11s & %-8s & %s \\\\",
                            tools::toTitleCase(sc), tools::toTitleCase(md),
                            paste(fmt(v + 1e-9, 2), collapse = " & ")))  # half-up, as typeset
    say(sprintf("  %s %s: %s", sc, md, paste(sprintf("%.4f", v), collapse = "  ")))
  }
  writeLines(c(
    "\\begin{table}[htbp]", "\\centering",
    "\\caption{\\textbf{Prior-informed CARMA sensitivity analysis.} Entries are F1 scores in the \\rev{single-locus} assessment. $\\mathrm{CARMA}_{p}$ denotes CARMA run with the same prior-weight information used by the prior-informed competitors, when applicable.}",
    "\\label{tab:carma_prior_sensitivity}", "\\small",
    "\\setlength{\\tabcolsep}{6pt}", "\\renewcommand{\\arraystretch}{1.12}",
    "\\begin{tabular}{llcccc}", "\\toprule",
    "\\textbf{Scenario} & \\textbf{Regime} & \\textbf{ExM} & \\textbf{SuSiE$_p$} & \\textbf{CARMA} & \\textbf{CARMA$_p$} \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(OUT, "tab_carma_prior_sensitivity.tex"))
}
if (has("results_locus_sim_informative_adaptive.rds")) {
  s <- read_summary("results_locus_sim_informative_adaptive.rds") %>%
    filter(n == 1000, method %in% INHOUSE) %>% select(method, criterion, mean_M, mean_K, F1)
  say("== Locus informative adaptive, n = 1000, in-house rows")
  capture.output(print(as.data.frame(s), digits = 3, row.names = FALSE)) %>% say()
}

# ---- 2b. Small n (n = 80), informative prior, tau = 0.2, averaged over rho ----
signed <- function(x, d) ifelse(is.na(x), "{--}", sprintf(paste0("%+.", d, "f"), x))
if (has("results_bench_synth_informative_n80.rds")) {
  s <- read_summary("results_bench_synth_informative_n80.rds") %>%
    filter(n == 80, criterion == "K_fixed", abs(tau - 0.2) < 1e-9) %>%
    group_by(method, gamma) %>% summarise(M = mean(mean_M), F1 = mean(F1), .groups = "drop")
  gams <- c(0.2, 0.6, 0.9)
  cell <- function(m, g) { r <- s[s$method == m & abs(s$gamma - g) < 1e-9, ]; c(r$M, r$F1) }
  line <- function(lab, m) paste0(lab, " & ", paste(vapply(gams, function(g) {
    v <- cell(m, g); sprintf("%.2f (%.3f)", v[1], v[2]) }, ""), collapse = " & "), " \\\\")
  delta <- function(lab, a, b) paste0(lab, " & ", paste(vapply(gams, function(g) {
    v <- cell(a, g) - cell(b, g); sprintf("%s (%s)", signed(v[1], 2), signed(v[2], 3)) }, ""),
    collapse = " & "), " \\\\")
  writeLines(c(
    "\\begin{table}[H]", "\\centering",
    "\\caption{\\textbf{Small-sample evaluation at $n=80$ (fixed $K$).} Synthetic AR(1) design, informative prior, $\\tau_j=0.2$, averaged over $\\rho\\in\\{0.8,0.95,0.98\\}$ ($B=50$). Cells: mean causal SNPs recovered out of $15$ (F1 in parentheses). \\rev{The two $\\Delta$ rows give the exact-minus-approximate gain, a finite-sample effect, below $0.2$ by $n=1000$.}}",
    "\\label{tab:smalln}", "\\scriptsize", "\\setlength{\\tabcolsep}{5pt}",
    "\\renewcommand{\\arraystretch}{1.12}", "\\begin{tabular}{lccc}", "\\toprule",
    "\\textbf{Method / contrast} & $\\gamma=0.2$ & $\\gamma=0.6$ & $\\gamma=0.9$ \\\\", "\\midrule",
    line("AsW", "Asymptotic_Wald"), line("ExW", "Exact_Wald"),
    delta("$\\Delta$(ExW-AsW)", "Exact_Wald", "Asymptotic_Wald"), "\\addlinespace[2pt]",
    line("AsS", "Asymptotic_Score"), "\\addlinespace[2pt]",
    line("NaM", "Naive_Multiple"), line("ExM", "Exact_Multiple"),
    delta("$\\Delta$(ExM-NaM)", "Exact_Multiple", "Naive_Multiple"), "\\addlinespace[2pt]",
    line("SuSiE\\textsubscript{prior}", "SuSiE_prior"), line("SuSiE\\textsubscript{unif}", "SuSiE_unif"),
    line("varbvs", "varbvs"), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(OUT, "tab_smalln_premium.tex"))
  say("== n = 80, fixed K, tau = 0.2 (M and F1 averaged over rho)")
  capture.output(print(as.data.frame(s), digits = 4, row.names = FALSE)) %>% say()
}
if (has("results_bench_synth_informative_n80_adaptive.rds")) {
  s <- read_summary("results_bench_synth_informative_n80_adaptive.rds") %>%
    filter(n == 80, abs(tau - 0.2) < 1e-9, method %in% INHOUSE) %>%
    group_by(method, gamma, criterion) %>%
    summarise(F1 = mean(F1), K = mean(mean_K), .groups = "drop") %>%
    group_by(method, gamma) %>% slice_max(F1, n = 1, with_ties = FALSE) %>% ungroup()
  rule <- c(R2 = "adj.~$R^2$", BIC = "BIC", eBIC = "eBIC", Post_0.5 = "$q>0.5$",
            Post_0.95 = "$q>0.95$", Post_0.99 = "$q>0.99$")
  gams <- c(0.2, 0.6, 0.9)
  get <- function(m, g) s[s$method == m & abs(s$gamma - g) < 1e-9, ]
  line <- function(lab, m) paste0(lab, " & ", paste(vapply(gams, function(g) {
    r <- get(m, g); sprintf("%.3f [%s; $\\hat K$=%.1f]", r$F1, rule[[r$criterion]], r$K) }, ""),
    collapse = " & "), " \\\\")
  delta <- function(lab, a, b) paste0(lab, " & ", paste(vapply(gams, function(g)
    signed(get(a, g)$F1 - get(b, g)$F1, 3), ""), collapse = " & "), " \\\\")
  writeLines(c(
    "\\begin{table}[H]", "\\centering",
    "\\caption{\\textbf{Small-sample evaluation at $n=80$, adaptive stopping (robustness check).} Same design as Table~\\ref{tab:smalln}, each in-house method stopping by its best data-driven rule; we report the best F1 (averaged over $\\rho$), the selected rule and mean size $\\hat K$ (iteration cap $30$). \\rev{The joint variants stop near the true size ($15$) under the posterior rule $q>0.5$.}}",
    "\\label{tab:smalln_adaptive}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.12}", "\\begin{tabularx}{\\textwidth}{lXXX}", "\\toprule",
    "\\textbf{Method / contrast} & $\\gamma=0.2$ & $\\gamma=0.6$ & $\\gamma=0.9$ \\\\", "\\midrule",
    line("AsW", "Asymptotic_Wald"), line("ExW", "Exact_Wald"),
    delta("$\\Delta$F1(ExW - AsW)", "Exact_Wald", "Asymptotic_Wald"), "\\addlinespace[2pt]",
    line("AsS", "Asymptotic_Score"), "\\addlinespace[2pt]",
    line("NaM", "Naive_Multiple"), line("ExM", "Exact_Multiple"),
    delta("$\\Delta$F1(ExM - NaM)", "Exact_Multiple", "Naive_Multiple"),
    "\\bottomrule", "\\end{tabularx}", "\\end{table}"),
    file.path(OUT, "tab_smalln_adaptive.tex"))
  say("== n = 80, adaptive, tau = 0.2: best rule per method (F1 and K averaged over rho)")
  capture.output(print(as.data.frame(s), digits = 4, row.names = FALSE)) %>% say()
}

# ---- 3. Repository scripts in a temporary workspace -------------------------
run_repo <- function(script, inputs, outputs, args = character(0)) {
  ws <- tempfile("regen_"); dir.create(file.path(ws, "code"), recursive = TRUE)
  file.copy(file.path(REPO, "code", script), file.path(ws, "code", script))
  for (x in inputs) file.symlink(x, file.path(ws, basename(x)))
  owd <- setwd(file.path(ws, "code")); on.exit(setwd(owd))
  log <- system2("Rscript", c(script, args), stdout = TRUE, stderr = TRUE)
  st <- attr(log, "status")
  if (!is.null(st) && st != 0) stop(script, " failed:\n", paste(log, collapse = "\n"))
  for (o in outputs) file.copy(file.path(ws, o), file.path(OUT, o), overwrite = TRUE)
  log
}

geu <- c("results_geuvadis_real_ld.rds", "results_geuvadis_real_ld_bonf.rds")
if (all(has(geu))) {
  semis_dir <- if (nzchar(cfg$semisynth)) normalizePath(cfg$semisynth) else REPO
  semis <- Sys.glob(file.path(semis_dir, "results_geuvadis_semisynth50_*.rds"))
  outs <- c("tab_geuvadis_fixedK.tex", "tab_geuvadis_adaptive.tex",
            "tab_geuvadis_exact_snp.tex", "tab_geuvadis_bonferroni.tex")
  if (nzchar(cfg$semisynth)) outs <- c(outs, "tab_geuvadis_semisynth.tex")
  log <- run_repo("fill_geuvadis_tables.R", c(f_in(geu), semis), outs)
  say("== GEUVADIS (fill_geuvadis_tables.R headline numbers)")
  say(paste(log[seq(which(grepl("headline numbers", log))[1], length(log))], collapse = "\n"))
  log <- run_repo("compute_pairwise_pi_geuvadis.R", f_in(geu[1]), "tab_pi_geuvadis.tex")
  say("== GEUVADIS pairwise pi (TSS prior, fixed K)")
  say(paste(log, collapse = "\n"))
  ld <- readRDS(f_in(geu[1]))$table
  say("== GEUVADIS LD-aware rate and size by method x prior x selection")
  tb <- ld %>% filter(selection %in% c("fixedK", "BIC")) %>%
    group_by(selection, method, prior) %>%
    summarise(rate = round(mean(rep_rate_ld, na.rm = TRUE), 3), K = round(mean(K, na.rm = TRUE), 1),
              R2_yri = round(mean(R2_yri, na.rm = TRUE), 3), .groups = "drop")
  capture.output(print(as.data.frame(tb), row.names = FALSE)) %>% say()
}

prem <- c("results_synth_premium_vs_n.rds", "results_synth_premium_vs_n_band.rds")
if (has(prem[1])) {
  run_repo("figure_premium_vs_n.R", f_in(prem[has(prem)]), "Figure_premium_vs_n.pdf")
  s <- bind_rows(lapply(f_in(prem[has(prem)]), function(f) readRDS(f)$summary)) %>%
    filter(tau == 0.2, gamma == 0.6, criterion == "K_fixed") %>%
    select(n, rho_intra, method, M = mean_M) %>%
    pivot_wider(names_from = method, values_from = M) %>%
    transmute(n, rho_intra, ExW_AsW = Exact_Wald - Asymptotic_Wald,
              ExM_NaM = Exact_Multiple - Naive_Multiple) %>%
    group_by(n) %>%
    summarise(across(c(ExW_AsW, ExM_NaM), list(mean = mean, min = min, max = max)),
              n_rho = dplyr::n(), .groups = "drop")
  say("== Gain vs n (extra causal SNPs /15; mean and range over rho)")
  capture.output(print(as.data.frame(s), digits = 3, row.names = FALSE)) %>% say()
}
cat("Staging directory:", OUT, "\n")
if (file.exists(HEAD)) cat(readLines(HEAD), sep = "\n")
