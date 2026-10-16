# =============================================================================
# Robustness_Checks.R
# Greedy selected-set check on controlled AR(1) high-LD loci: averaged over B
# replicates, do NaM and ExM (greedy) recover the causal clusters as often as
# SuSiE (non-greedy)? Recovery is counted within +/-2 SNPs of each causal.
#
# NaM and ExM use the exact prior integrals (S = Inf); the Monte-Carlo version
# (S = 100 draws) is reported alongside on the same data sets. Each replicate
# has its own seed, so that all scores are compared on identical data. The
# prior-free baseline uses a uniform prior (p = 0.5) for the stepwise scores.
#
# Usage (from code/):  Rscript Robustness_Checks.R
# =============================================================================
suppressPackageStartupMessages({ library(MASS); library(susieR) })
source("Test_statistics_functions_optimized.R")

m <- 100L; n <- 1000L; tau2 <- 0.2^2; K <- 3L; B <- 50L
causal   <- c(20L, 50L, 80L)
beta_str <- 8 / sqrt(n)                                    # strong marginal |z| ~ 8
p_info   <- rep(0.2, m); p_info[causal] <- 0.8             # informative prior
p_unif   <- rep(0.5, m)
chol_cache <- list()
gen_X <- function(rho) {
  key <- as.character(rho)
  if (is.null(chol_cache[[key]]))
    chol_cache[[key]] <<- chol(rho^abs(outer(seq_len(m), seq_len(m), "-")))
  matrix(rnorm(n * m), n, m) %*% chol_cache[[key]]
}
within <- function(sel, tol = 2L) sum(vapply(causal, function(c0) any(abs(sel - c0) <= tol), NA))

run_arm <- function(rho, arm_id) {
  rec <- matrix(NA_real_, B, 7L, dimnames = list(NULL,
    c("ExM", "NaM", "SuSiE", "ExM_MC100", "NaM_MC100", "ExM_unif", "NaM_unif")))
  for (b in seq_len(B)) {
    set.seed(2026L + 1000L * arm_id + b)
    Xb <- gen_X(rho)
    Yb <- as.numeric(scale(Xb, scale = FALSE)[, causal] %*% rep(beta_str, 3) + rnorm(n))
    rec[b, "ExM"]   <- within(Exact_Multiple(Xb, Yb, p_info, K, tau2, Inf)$indices)
    rec[b, "NaM"]   <- within(Naive_Multiple(Xb, Yb, p_info, K, tau2, Inf)$indices)
    rec[b, "SuSiE"] <- within(order(susie(Xb, Yb, L = K, prior_weights = p_info / sum(p_info))$pip,
                                    decreasing = TRUE)[seq_len(K)])
    rec[b, "ExM_unif"] <- within(Exact_Multiple(Xb, Yb, p_unif, K, tau2, Inf)$indices)
    rec[b, "NaM_unif"] <- within(Naive_Multiple(Xb, Yb, p_unif, K, tau2, Inf)$indices)
    set.seed(9000L + 1000L * arm_id + b)                   # Monte-Carlo draws only
    rec[b, "ExM_MC100"] <- within(Exact_Multiple(Xb, Yb, p_info, K, tau2, 100L)$indices)
    rec[b, "NaM_MC100"] <- within(Naive_Multiple(Xb, Yb, p_info, K, tau2, 100L)$indices)
  }
  colMeans(rec)
}
arms <- lapply(seq_along(c(0.90, 0.95)), function(i) run_arm(c(0.90, 0.95)[i], i))
names(arms) <- c("0.90", "0.95")
cat(sprintf("\nMean within-cluster causal recovery (out of 3) over B = %d replicates, m = %d, n = %d:\n", B, m, n))
print(round(do.call(rbind, arms), 2))
saveRDS(list(arms = arms, B = B, m = m, n = n, tau2 = tau2, causal = causal),
        "../results_robustness_checks.rds")
cat("\nSaved ../results_robustness_checks.rds\n")
