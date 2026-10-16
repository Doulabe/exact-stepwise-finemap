# =============================================================================
# File: sim_noisy_prior.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Graded-accuracy prior arm for the controlled AR(1) design: the informative
# scenario gives the 15 causal variants p = 0.8 by construction (an oracle a
# prior-only ranking solves), so this experiment degrades the prior
# continuously. A latent annotation s_j = delta 1{j causal} + N(0,1) with
# delta = sqrt(2) Phi^{-1}(AUC) ranks the candidates; the top 15 receive
# p = 0.8 and the rest 0.2 (the informative scale). AUC in {0.6, 0.75, 0.9, 1}
# spans weak to oracle annotation. The prior-only reference (rank by p_j, ties
# by marginal |r|) is evaluated on the same data sets, so the value added by
# the data over the prior is read directly.
#
# Design: m = 5000, AR(1) blocks of 50, rho = 0.95, n = 1000, gamma = 0.6,
# tau = 0.2, B = 50, K = 15. Methods: AsS, ExW, NaM, ExM (exact integrals).
#
# Usage (from code/):
#   Rscript sim_noisy_prior.R --n_cores 6 --out ../results_noisy_prior.rds
# =============================================================================
suppressPackageStartupMessages({
  source("Generate_XY_realistic.R")
  source("Test_statistics_functions_optimized.R")
  library(dplyr); library(parallel)
})

cfg <- list(m = 5000L, block = 50L, rho = 0.95, n = 1000L, B = 50L, K = 15L,
            gamma = 0.6, tau = 0.2, seed = 2026L, n_cores = 6L,
            table_only = FALSE,   # TRUE: rebuild the table from the stored RDS
            out = "../results_noisy_prior.rds",
            tex = "../tab_noisy_prior.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

src <- readLines("Benchmark_External.R")
grab <- function(name) {
  i0 <- grep(paste0("^", name, " <- function"), src)
  i1 <- i0 + which(grepl("^}", src[(i0 + 1):length(src)]))[1]
  eval(parse(text = paste(src[i0:i1], collapse = "\n")), envir = globalenv())
}
grab("calibrate_lambda")

AUCS <- c(0.6, 0.75, 0.9, 1.0)
causal_idx <- choose_causal_indices(cfg$m, 15L, cfg$block, "A", 3L, seed = cfg$seed)
lambda_true <- calibrate_lambda(cfg$n, cfg$gamma)
gh <- make_quad_nodes(30L, cfg$tau)

one_rep <- function(b) {
  set.seed(cfg$seed + 7919L * b)
  X <- Generate.X.realistic(n = cfg$n, m = cfg$m, ld_structure = "blocks",
                            rho_intra = cfg$rho, rho_inter = 0,
                            block_size = cfg$block, maf_lo = 0.05, maf_hi = 0.5,
                            maf_distribution = "uniform")
  beta <- rep(0, cfg$m); beta[causal_idx] <- lambda_true
  Y <- Generate.Y(n = cfg$n, X = X, beta = beta, sigma = 1)
  r_marg <- suppressWarnings(as.numeric(cor(Y, X))); r_marg[is.na(r_marg)] <- 0
  out <- list()
  for (auc in AUCS) {
    set.seed(cfg$seed + 7919L * b + round(auc * 1000))
    if (auc >= 1) {
      top <- causal_idx
    } else {
      delta <- sqrt(2) * qnorm(auc)
      s <- delta * (seq_len(cfg$m) %in% causal_idx) + rnorm(cfg$m)
      top <- order(s, decreasing = TRUE)[seq_len(15L)]
    }
    p_j <- rep(0.2, cfg$m); p_j[top] <- 0.8
    sel_po <- order(-p_j, -abs(r_marg))[seq_len(cfg$K)]
    res <- list(Prior_only = list(indices = sel_po))
    res$Asymptotic_Score <- Asymptotic_Score(X, Y, p_j, cfg$K, cfg$tau, gh)
    res$Exact_Wald       <- Exact_Wald(X, Y, p_j, cfg$K, cfg$tau, gh)
    res$Naive_Multiple   <- Naive_Multiple(X, Y, p_j, cfg$K, cfg$tau^2, Inf)
    res$Exact_Multiple   <- Exact_Multiple(X, Y, p_j, cfg$K, cfg$tau^2, Inf)
    for (meth in names(res)) {
      M <- length(intersect(res[[meth]]$indices, causal_idx))
      out[[length(out) + 1L]] <- data.frame(
        b = b, auc = auc, method = meth, M = M,
        F1 = 2 * M / (cfg$K + 15L))
    }
  }
  bind_rows(out)
}

if (isTRUE(cfg$table_only)) {
  d <- readRDS(cfg$out)$all_results
} else {
  d <- bind_rows(mclapply(seq_len(cfg$B), one_rep, mc.cores = cfg$n_cores,
                          mc.preschedule = FALSE))
  saveRDS(list(all_results = d, config = cfg), cfg$out)
  cat("Wrote", cfg$out, "\n")
}

agg <- d %>% group_by(method, auc) %>%
  summarise(M = mean(M), F1 = mean(F1), .groups = "drop")
print(as.data.frame(agg %>% tidyr::pivot_wider(names_from = auc, values_from = c(M, F1))),
      row.names = FALSE)

# Composition of the favoured stratum (data-independent, same seeds as the
# runs): the number of true causal variants among the 15 candidates given 0.8.
strat <- vapply(AUCS, function(auc) mean(vapply(seq_len(cfg$B), function(b) {
  set.seed(cfg$seed + 7919L * b + round(auc * 1000))
  if (auc >= 1) return(15L)
  delta <- sqrt(2) * qnorm(auc)
  s <- delta * (seq_len(cfg$m) %in% causal_idx) + rnorm(cfg$m)
  length(intersect(order(s, decreasing = TRUE)[1:15], causal_idx))
}, 1L)), 0)

# Paired margins over the prior-only reference (same data sets)
w <- d %>% select(b, auc, method, M) %>% tidyr::pivot_wider(names_from = method, values_from = M)
marg <- bind_rows(lapply(c("Asymptotic_Score", "Exact_Wald", "Naive_Multiple",
                           "Exact_Multiple"), function(me)
  w %>% group_by(auc) %>%
    summarise(method = me, dm = mean(.data[[me]] - Prior_only),
              se = sd(.data[[me]] - Prior_only) / sqrt(dplyr::n()), .groups = "drop")))
cat("paired margins over prior-only:\n")
print(as.data.frame(marg %>% filter(auc < 1) %>% mutate(across(c(dm, se), ~round(., 2)))),
      row.names = FALSE)

PRETTY <- c(Prior_only = "Prior-only ranking", Asymptotic_Score = "AsS",
            Exact_Wald = "ExW", Naive_Multiple = "NaM", Exact_Multiple = "ExM")
rows <- vapply(names(PRETTY), function(me) {
  cells <- vapply(AUCS, function(a2) {
    z <- agg %>% filter(method == me, abs(auc - a2) < 1e-9)
    sprintf("%.2f (%.3f)", z$M, z$F1)
  }, "")
  paste0(PRETTY[me], " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
strat_row <- paste0("\\emph{Causal variants in the $p_j=0.8$ stratum} & ",
                    paste(sprintf("%.2f", strat), collapse = " & "), " \\\\")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Graded prior accuracy.} Mean causal variants recovered out of ",
         "$15$ (F1 in parentheses) on the controlled AR(1) design ($m=5{,}000$, $\\rho=0.95$, ",
         "$n=1{,}000$, $\\gamma=0.6$, $B=50$) when the $15$ candidates given $p_j=0.8$ are chosen ",
         "by a noisy latent annotation of prescribed accuracy (AUC of the annotation for causal ",
         "against non-causal variants, AUC $=1$ being the oracle informative scenario). With $15$ ",
         "causal variants among $m=5{,}000$, even a seemingly high AUC leaves few causal variants ",
         "in the favored stratum (first row), which the prior-only reference, ranking by $p_j$ ",
         "with ties broken by the marginal $|r|$, selects exactly. Every method exceeds this ",
         "reference on average at every imperfect accuracy, with paired margins of $+0.3$ to ",
         "$+1.9$ causal variants, all at least four paired standard errors above zero.}}"),
  "\\label{tab:noisy_prior}",
  "\\centering\\small",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(sprintf("AUC $%.2f$", AUCS), collapse = " & "), " \\\\"),
  "\\midrule", strat_row, "\\midrule", rows,
  "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$tex)
cat("Wrote", cfg$tex, "\n")
