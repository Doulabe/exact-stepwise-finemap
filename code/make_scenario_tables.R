# =============================================================================
# make_scenario_tables.R  (Article 1, SiM supplement)
# Builds three supplementary tables from the exact-run result files:
#   tab_synth_scenarios.tex  AR(1) m = 5000, fixed K = 15, four prior scenarios,
#                            n in {80, 1000}, all methods (Stage A files)
#   tab_locus_scenarios.tex  locus m = 200, fixed K = 3, four scenarios, n = 1000,
#                            with CARMA, CARMA_p and FINEMAP (merged locus files)
#   tab_stop_joint.tex       stopping rules for the joint scores (ExM; NaM agrees)
# Usage (from code/):
#   Rscript make_scenario_tables.R --in results_exact --out staging_exact
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
get  <- function(k, d) { i <- match(k, args); if (is.na(i)) d else args[i + 1L] }
IN   <- get("--in", "results_exact"); OUT <- get("--out", "staging_exact")
PO   <- get("--prior_only", "../results_prior_only.rds")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SC   <- c(informative = "Informative", uniform = "Uniform", misspecified = "Misspecified", mixed = "Mixed")
LAB  <- c(Naive_Wald = "NaW", Asymptotic_Wald = "AsW", Exact_Wald = "ExW", Asymptotic_Score = "AsS",
          Naive_Multiple = "NaM", Exact_Multiple = "ExM",
          SuSiE_unif = "SuSiE\\textsubscript{unif}", SuSiE_prior = "SuSiE\\textsubscript{prior}",
          varbvs = "varbvs", CARMA = "CARMA", CARMA_prior = "CARMA\\textsubscript{p}", FINEMAP = "FINEMAP",
          LASSO_topK = "LASSO\\textsubscript{topK}", LASSO_1se = "LASSO\\textsubscript{1se}")
fmt  <- function(x, d = 2) ifelse(is.na(x) | !is.finite(x), "{--}", sprintf(paste0("%.", d, "f"), x + 1e-9))
summ <- function(f) as.data.frame(readRDS(file.path(IN, f))$summary)

# ---- 1. AR(1), four scenarios, all methods ----------------------------------
sa <- bind_rows(lapply(names(SC), function(s) summ(sprintf("results_synth_scen_%s.rds", s)) %>% mutate(scen = s))) %>%
  filter(criterion == "K_fixed") %>%
  group_by(scen, n, method) %>% summarise(M = mean(mean_M), .groups = "drop")
meths <- c("Naive_Wald", "Asymptotic_Wald", "Exact_Wald", "Asymptotic_Score", "Naive_Multiple", "Exact_Multiple",
           "SuSiE_unif", "SuSiE_prior", "LASSO_topK", "LASSO_1se")
row1 <- function(m) {
  v <- unlist(lapply(c(80, 1000), function(nn) sapply(names(SC), function(s) {
    x <- sa$M[sa$method == m & sa$n == nn & sa$scen == s]; if (length(x)) fmt(x) else "{--}" })))
  paste0(LAB[[m]], " & ", paste(v, collapse = " & "), " \\\\")
}
# Prior-only reference row (rank by p_j, ties by marginal |r|), computed on
# the same data sets by compute_prior_only.R.
po_row <- character(0)
if (file.exists(PO)) {
  po <- readRDS(PO)$all_results %>%
    filter(grepl("^results_synth_scen_", file)) %>%
    group_by(scenario, n) %>% summarise(M = mean(M), .groups = "drop")
  v <- unlist(lapply(c(80, 1000), function(nn) sapply(names(SC), function(s) {
    x <- po$M[po$scenario == s & po$n == nn]; if (length(x)) fmt(x) else "{--}" })))
  po_row <- c(paste0("\\rev{Prior-only ranking} & ", paste(v, collapse = " & "), " \\\\"),
              "\\midrule")
}
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{\\rev{\\textbf{Synthetic AR(1) design under the four prior scenarios.} Mean number of causal SNPs recovered (out of $15$) at fixed $K=15$, averaged over $\\rho\\in\\{0.8,0.95,0.98\\}$ ($m=5{,}000$, $\\gamma=0.6$, $\\tau_j=0.2$, $B=50$). Scenarios as in Table~\\ref{tab:parameters}. The first row is the prior-only reference, which ranks by $p_j$ with ties broken by the marginal $|r|$ on the same data sets. It bounds what prior tracking alone achieves in each scenario, and the margin of a method over it is the contribution of the data.}}",
  "\\label{tab:synth_scenarios}", "\\scriptsize", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lcccccccc}", "\\toprule",
  " & \\multicolumn{4}{c}{$n=80$} & \\multicolumn{4}{c}{$n=1000$} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}",
  paste0("Method & ", paste(rep(c("Inf.", "Unif.", "Missp.", "Mixed"), 2), collapse = " & "), " \\\\"),
  "\\midrule", po_row, vapply(meths, row1, ""), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(OUT, "tab_synth_scenarios.tex"))

# ---- 2. Locus, four scenarios, with CARMA / FINEMAP -------------------------
lo <- bind_rows(lapply(names(SC), function(s) summ(sprintf("results_locus_sim_%s_fixed.rds", s)) %>% mutate(scen = s))) %>%
  filter(n == 1000, criterion == "K_fixed")
lm <- c("Naive_Wald", "Asymptotic_Wald", "Exact_Wald", "Asymptotic_Score", "Naive_Multiple", "Exact_Multiple",
        "SuSiE_unif", "SuSiE_prior", "varbvs", "CARMA", "CARMA_prior", "FINEMAP", "LASSO_1se")
row2 <- function(m) {
  v <- sapply(names(SC), function(s) { x <- lo$mean_M[lo$method == m & lo$scen == s]; if (length(x)) fmt(x) else "{--}" })
  paste0(LAB[[m]], " & ", paste(v, collapse = " & "), " \\\\")
}
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  "\\caption{\\rev{\\textbf{Single-locus design under the four prior scenarios.} Mean number of causal variants recovered (out of $3$) at fixed $K=3$ ($m_{\\text{locus}}=200$, $\\rho=0.9$, $n=1000$, $B=200$). CARMA, FINEMAP and SuSiE\\textsubscript{unif} ignore the prior, while CARMA\\textsubscript{p} receives $\\mathrm{logit}(p_j)$ as an annotation whose weight it estimates.}}",
  "\\label{tab:locus_scenarios}", "\\small", "\\setlength{\\tabcolsep}{6pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(SC, collapse = " & "), " \\\\"), "\\midrule",
  vapply(lm, row2, ""), "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(OUT, "tab_locus_scenarios.tex"))

# ---- 3. Stopping rules for the joint scores ----------------------------------
rules <- c(R2 = "adj.~$R^2$", BIC = "BIC", eBIC = "eBIC", Post_0.5 = "$q>0.5$", Post_0.95 = "$q>0.95$", Post_0.99 = "$q>0.99$")
crit_tab <- function(d) d %>% filter(criterion %in% names(rules)) %>%
  group_by(method, criterion) %>% summarise(M = mean(mean_M), K = mean(mean_K), .groups = "drop")
blocks <- list()
for (s in names(SC)) blocks[[paste("ar", s)]] <- list(lab = SC[[s]], d = crit_tab(summ(sprintf("results_bench_synth_%s_adaptive.rds", s)) %>% filter(n == 1000)))
blocks[["ar80"]] <- list(lab = "Informative, $n=80$",
  d = crit_tab(summ("results_bench_synth_informative_n80_adaptive.rds") %>% filter(abs(gamma - 0.6) < 1e-9, abs(tau - 0.2) < 1e-9)))
for (s in names(SC)) blocks[[paste("lo", s)]] <- list(lab = SC[[s]], d = crit_tab(summ(sprintf("results_locus_sim_%s_adaptive.rds", s)) %>% filter(n == 1000)))
row3 <- function(b) {
  d <- b$d %>% filter(method == "Exact_Multiple")
  v <- sapply(names(rules), function(r) { x <- d[d$criterion == r, ]; if (nrow(x)) sprintf("%s [%.1f]", fmt(x$M), x$K) else "{--}" })
  paste0("\\quad ", b$lab, " & ", paste(v, collapse = " & "), " \\\\")
}
gaps <- sapply(blocks, function(b) { w <- b$d %>% select(method, criterion, M) %>%
  pivot_wider(names_from = method, values_from = M); max(abs(w$Exact_Multiple - w$Naive_Multiple), na.rm = TRUE) })
gap <- gaps[["ar80"]]; gap_other <- max(gaps[names(gaps) != "ar80"])
cat(sprintf("max |M(ExM) - M(NaM)|: n = 80 %.3f, other designs %.3f\n", gap, gap_other))
nb <- names(blocks)
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  sprintf("\\caption{\\rev{\\textbf{Stopping rules for the joint scores.} Mean number of causal variants recovered by ExM and mean selected size $\\hat K$ (in brackets) under each rule, for $15$ causal SNPs (AR(1), averaged over $\\rho$, cap $30$) and $3$ causal variants (locus), $n=1000$ unless stated. NaM differs by at most $%.2f$ at $n=1000$ and $%.1f$ at $n=80$.}}", max(0.01, ceiling(gap_other * 100) / 100), ceiling(gap * 10) / 10),
  "\\label{tab:stop_joint}", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
  paste0("\\begin{tabular}{l", strrep("c", length(rules)), "}"), "\\toprule",
  paste0("Design and prior & ", paste(rules, collapse = " & "), " \\\\"), "\\midrule",
  "\\multicolumn{7}{l}{\\textit{AR(1), $m=5{,}000$}} \\\\",
  vapply(blocks[grep("^ar", nb)], row3, ""),
  "\\addlinespace[2pt]", "\\multicolumn{7}{l}{\\textit{Locus, $m_{\\text{locus}}=200$}} \\\\",
  vapply(blocks[grep("^lo", nb)], row3, ""),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(OUT, "tab_stop_joint.tex"))
cat("written: tab_synth_scenarios.tex, tab_locus_scenarios.tex, tab_stop_joint.tex\n")
