# =============================================================================
# File: compute_prior_only.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Prior-only reference ranking for the synthetic AR(1) designs: candidates are
# ordered by the prior weight p_j alone, with ties broken by the marginal |r|,
# and the first K are selected. This reference separates what the prior alone
# achieves from what the data add: under the oracle informative scenario it is
# essentially perfect by construction, so any method matching it there is
# tracking the prior, and the value added by the data appears in the partial
# scenarios. The datasets are regenerated with the exact per-replicate seeds of
# Benchmark_External.R, so the reference is evaluated on the same data as the
# stored results.
#
# Usage (from code/):
#   Rscript compute_prior_only.R \
#     --results_dir ../../server_rerun_fix_20260928/results_exact \
#     --out ../results_prior_only.rds
# =============================================================================
suppressPackageStartupMessages({source("Generate_XY_realistic.R"); library(dplyr)})

cfg <- list(results_dir = "../../server_rerun_fix_20260928/results_exact",
            K = 15L, out = "../results_prior_only.rds")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}
cfg$K <- as.integer(cfg$K)

# calibrate_lambda and build_prior exactly as in Benchmark_External.R
src <- readLines("Benchmark_External.R")
grab <- function(name) {
  i0 <- grep(paste0("^", name, " <- function"), src)
  i1 <- i0 + which(grepl("^}", src[(i0 + 1):length(src)]))[1]
  eval(parse(text = paste(src[i0:i1], collapse = "\n")), envir = globalenv())
}
grab("calibrate_lambda"); grab("build_prior")

FILES <- c("results_bench_synth_informative.rds",
           "results_synth_scen_informative.rds", "results_synth_scen_mixed.rds",
           "results_synth_scen_uniform.rds", "results_synth_scen_misspecified.rds")

rows <- list()
for (f in FILES) {
  path <- file.path(cfg$results_dir, f)
  if (!file.exists(path)) next
  rc <- readRDS(path)$config
  n_grid   <- as.integer(strsplit(rc$n_grid, ",")[[1L]])
  rho_grid <- as.numeric(strsplit(rc$rho_grid, ",")[[1L]])
  gamma_v  <- if (nzchar(rc$gamma_grid)) as.numeric(strsplit(rc$gamma_grid, ",")[[1L]])[1L] else rc$gamma
  tau_v    <- if (nzchar(rc$tau_grid))   as.numeric(strsplit(rc$tau_grid,   ",")[[1L]])[1L] else rc$tau
  m        <- rc$m_synth
  # pre-warm the AR(1) factors (hits the AR1_FACTOR_DIR disk cache); without
  # this, every Generate.X.realistic call at m = 5000 redoes a full eigen()
  ar1_prewarm_factors(m, rho_grid, rc$rho_inter, rc$block_size, cores = 1L)
  causal_idx <- choose_causal_indices(m, rc$n_causal, rc$block_size,
                                      rc$placement, rc$causals_per_window,
                                      seed = rc$seed)
  p_j <- build_prior(rc$scenario, m, causal_idx)
  cat(sprintf("%s: scen=%s n={%s} rho={%s} B=%d m=%d\n", f, rc$scenario,
              rc$n_grid, rc$rho_grid, rc$B, m))
  for (n in n_grid) for (rho in rho_grid) for (b in seq_len(rc$B)) {
    set.seed(rc$seed + b * 1000L + n * 7L + as.integer(rho * 1e4) +
             as.integer(gamma_v * 1e4) + as.integer(tau_v * 1e4) * 13L)
    lambda_true <- calibrate_lambda(n, gamma_v)
    X <- Generate.X.realistic(n = n, m = m, ld_structure = "blocks",
                              rho_intra = rho, rho_inter = rc$rho_inter,
                              block_size = rc$block_size, maf_lo = rc$maf_lo,
                              maf_hi = rc$maf_hi,
                              maf_distribution = rc$maf_distribution)
    beta <- rep(0, m); beta[causal_idx] <- lambda_true
    Y <- Generate.Y(n = n, X = X, beta = beta, sigma = 1)
    r <- suppressWarnings(as.numeric(cor(Y, X))); r[is.na(r)] <- 0
    sel <- order(-p_j, -abs(r))[seq_len(cfg$K)]
    M <- length(intersect(sel, causal_idx))
    rows[[length(rows) + 1L]] <- data.frame(
      file = f, scenario = rc$scenario, n = n, rho_intra = rho, b = b,
      M = M, K = cfg$K, F1 = 2 * M / (cfg$K + rc$n_causal))
  }
}
d <- bind_rows(rows)
saveRDS(list(all_results = d), cfg$out)
cat("Wrote", cfg$out, "\n")
print(as.data.frame(d %>% group_by(file, scenario, n) %>%
  summarise(M = round(mean(M), 2), F1 = round(mean(F1), 3), .groups = "drop")),
  row.names = FALSE)
