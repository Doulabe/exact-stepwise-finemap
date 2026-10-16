# =============================================================================
# File: Benchmark_Locus_Sim.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe 
#
# Description:
#   Single-locus simulation benchmark. Complements Benchmark_External.R
#   (genome-wide) by testing all comparators in the regime for which CARMA,
#   FINEMAP, SuSiE and varbvs are actually designed: a single LD block of
#   m ~ 200 SNPs with a small number of causal variants.
#
#   For each replication we generate
#     - a fresh AR(1) LD block of m_locus SNPs (rho_intra controls strength)
#     - q causal SNPs placed at spread positions within the locus
#     - a continuous Y with effect size lambda calibrated to power gamma
#
#   We then run six in-house methods + SuSiE (no-prior + with-prior) +
#   varbvs + LASSO + CARMA + FINEMAP. The wall-time on m=200 is seconds.
#
#   Two K modes are reported:
#     - K_mode = "fixed"    -> top-K_fixed selection per method.
#     - K_mode = "adaptive" -> 4 criteria for in-house (R^2, BIC, eBIC,
#                              Posterior @ 0.5/0.95/0.99) and the natural
#                              adaptive setting of each external method
#                              (SuSiE 95% CS, varbvs PIP>0.5/0.95/0.99,
#                              LASSO lambda.1se, CARMA PIP>0.5/0.95/0.99,
#                              FINEMAP PIP>0.5/0.95/0.99).
#
#   Metrics: M_hat, K_hat, precision, recall, F1, FDR, runtime (one row
#   per method x criterion x n x replication). Output format matches that
#   of Benchmark_External.R so the table generators can render them directly.
#
# Usage:
#   Rscript Benchmark_Locus_Sim.R \
#     --B 50 --n_cores 12 \
#     --scenario informative --K_mode fixed \
#     --m_locus 200 --rho_intra 0.9 --q_locus 3 \
#     --K_fixed 3 --gamma 0.6 --tau 0.2 \
#     --outfile ../results_locus_sim_informative.rds
# =============================================================================

suppressPackageStartupMessages({
  library(MASS); library(stats); library(Matrix); library(mvtnorm)
  library(foreach); library(doParallel); library(dplyr); library(statmod)
})
source("Generate_XY_realistic.R")
source("Test_statistics_functions_optimized.R")
source("Test_statistics_functions_Criteria_optimized.R")
source("External_Comparators.R")

# ---------- 1. Configuration -------------------------------------------------
config_default <- list(
  B            = 50L,
  n_grid       = "500,1000",
  m_locus      = 200L,
  rho_intra    = 0.9,
  block_size   = 200L,   # one block = the entire locus
  q_locus      = 3L,
  K_fixed      = 3L,
  K_mode       = "fixed",
  K_threshold  = 30L,
  target_prob  = 0.95,
  gamma        = 0.6,
  tau          = 0.2,
  S_mc         = 100L,
  multiple_exact = FALSE,   # TRUE: NaM/ExM prior integrals computed exactly (S_mc unused)
  scenario     = "informative",
  seed         = 2026L,
  n_cores      = parallel::detectCores() - 1L,
  maf_lo       = 0.05,
  maf_hi       = 0.50,
  maf_distribution = "uniform",
  outfile      = "results_locus_sim.rds",
  inhouse_only     = "",     # comma-separated subset of in-house methods (empty = all)
  skip_externals   = FALSE,  # if TRUE, run the in-house methods only
  only_carma_prior = FALSE   # if TRUE, run ONLY the prior-aware CARMA arm
)
args <- commandArgs(trailingOnly = TRUE)
parse_cli <- function(args, defaults) {
  i <- 1L
  while (i <= length(args) - 1L) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i]); val <- args[i + 1L]
      if (!is.null(defaults[[key]])) {
        defaults[[key]] <- methods::as(val, class(defaults[[key]]))
      }
      i <- i + 2L
    } else i <- i + 1L
  }
  defaults
}
cfg <- parse_cli(args, config_default)
stopifnot(cfg$K_mode %in% c("fixed", "adaptive"))
n_values <- as.integer(strsplit(cfg$n_grid, ",", fixed = TRUE)[[1L]])
cat("Configuration:\n"); str(cfg); cat("n grid:", n_values, "\n")

# ---------- 2. Causal placement and prior -----------------------------------
# Causals spread within the locus
causal_idx <- as.integer(round(seq(0.10, 0.90, length.out = cfg$q_locus)
                                * cfg$m_locus))
causal_idx <- sort(unique(causal_idx))
cat("Causal SNPs (locus indices):", causal_idx, "\n")

build_prior <- function(scenario, m, causal_idx, target_set = integer(0)) {
  switch(
    scenario,
    informative = {
      v <- rep(0.2, m); v[causal_idx] <- 0.8; v
    },
    uniform = rep(0.5, m),
    misspecified = {
      v <- rep(0.5, m); v[causal_idx] <- 0.2
      decoy <- setdiff(seq_len(m), causal_idx)
      n_decoy <- min(length(causal_idx), length(decoy))
      v[sample(decoy, n_decoy)] <- 0.8; v
    },
    mixed = {
      v <- rep(0.5, m); q <- length(causal_idx)
      v[causal_idx[seq_len(min(1L, q))]] <- 0.2
      if (q >= 2L) v[causal_idx[seq(2L, min(2L, q))]] <- 0.5
      if (q >= 3L) v[causal_idx[seq(3L, q)]] <- 0.8
      v
    },
    stop("Unknown scenario")
  )
}
set.seed(cfg$seed)   # make the misspecified decoy draw reproducible
p_j <- build_prior(cfg$scenario, cfg$m_locus, causal_idx)

# ---------- 3. Effect size calibration ---------------------------------------
calibrate_lambda <- function(n, gamma, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_gamma <- qnorm(gamma)
  (z_alpha + z_gamma) / sqrt(n)
}

# ---------- 4. Helpers shared with Benchmark_External ----------------------
build_inhouse_args <- function(fn_name, X, Y, K, gh) {
  base <- list(X = X, Y = Y, a.priori = p_j, K = K)
  extras <- switch(fn_name,
    Naive_Wald       = list(tau = cfg$tau),
    Asymptotic_Wald  = list(tau = cfg$tau, gh = gh),
    Exact_Wald       = list(tau = cfg$tau, gh = gh),
    Asymptotic_Score = list(tau = cfg$tau, gh = gh),
    Naive_Multiple   = list(tau2 = cfg$tau^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
    Exact_Multiple   = list(tau2 = cfg$tau^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
    list(tau = cfg$tau))
  c(base, extras)
}
build_criteria_args <- function(fn_name, X, Y, crit, gh) {
  base <- list(X = X, Y = Y, a.priori = p_j,
               threshold = cfg$K_threshold,
               target_prob = cfg$target_prob, criteria = crit)
  extras <- switch(fn_name,
    Naive_Wald       = list(tau = cfg$tau),
    Asymptotic_Wald  = list(tau = cfg$tau, gh = gh),
    Exact_Wald       = list(tau = cfg$tau, gh = gh),
    Asymptotic_Score = list(tau = cfg$tau, gh = gh),
    Naive_Multiple   = list(tau2 = cfg$tau^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
    Exact_Multiple   = list(tau2 = cfg$tau^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
    list(tau = cfg$tau))
  c(base, extras)
}
derive_stop_k <- function(hist, want_increasing) {
  if (length(hist) <= 1L) return(length(hist))
  diffs <- diff(hist)
  bad <- if (want_increasing) which(diffs < 0) else which(diffs > 0)
  if (length(bad) == 0L) length(hist) else bad[1L]
}
# Stringify an integer vector for storage in raw rows. Empty -> "".
indices_to_csv <- function(idx) {
  if (length(idx) == 0L) return("")
  paste(as.integer(idx), collapse = ",")
}

# ---------- 5. One replication -----------------------------------------------
one_replication <- function(b, n) {
  set.seed(cfg$seed + b * 1000L + n * 7L)
  lambda_true <- calibrate_lambda(n, cfg$gamma)
  X <- Generate.X.realistic(
    n = n, m = cfg$m_locus,
    ld_structure = "blocks", rho_intra = cfg$rho_intra,
    rho_inter = 0, block_size = cfg$block_size,
    maf_lo = cfg$maf_lo, maf_hi = cfg$maf_hi,
    maf_distribution = cfg$maf_distribution
  )
  beta <- rep(0, cfg$m_locus); beta[causal_idx] <- lambda_true
  Y <- Generate.Y(n = n, X = X, beta = beta, sigma = 1)
  m <- cfg$m_locus
  K <- cfg$K_fixed
  gh <- make_quad_nodes(30L, cfg$tau)
  results <- list()

  in_house <- c("Naive_Wald","Asymptotic_Wald","Exact_Wald",
                "Asymptotic_Score","Naive_Multiple","Exact_Multiple")
  if (nzchar(cfg$inhouse_only)) {
    short <- c(NaW = "Naive_Wald", AsW = "Asymptotic_Wald", ExW = "Exact_Wald",
               AsS = "Asymptotic_Score", NaM = "Naive_Multiple", ExM = "Exact_Multiple")
    sel <- strsplit(cfg$inhouse_only, ",", fixed = TRUE)[[1L]]
    sel <- ifelse(sel %in% names(short), short[sel], sel)
    in_house <- intersect(in_house, sel)
  }
  shape_criteria <- c("R2","BIC","eBIC")
  post_thr <- c(0.5, 0.95, 0.99)

  if (!isTRUE(cfg$only_carma_prior)) {   # skip all other methods in CARMA-only mode
  if (cfg$K_mode == "fixed") {
    for (fn_name in in_house) {
      if (!exists(fn_name)) next
      t0 <- Sys.time()
      out <- tryCatch(do.call(fn_name, build_inhouse_args(fn_name, X, Y, K, gh)),
                       error = function(e) NULL)
      rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(out)) next
      met <- evaluate_selection(out$indices, causal_idx, m)
      results[[fn_name]] <- c(method = fn_name, criterion = "K_fixed",
                              b = b, n = n, rho_intra = cfg$rho_intra,
                              runtime = rt, met,
                              indices_csv = indices_to_csv(out$indices))
    }
  } else {
    # Shape criteria via single ALL fit
    for (fn_name in in_house) {
      fn_adapt <- paste0(fn_name, "_Criteria")
      if (!exists(fn_adapt)) next
      t0 <- Sys.time()
      out <- tryCatch(do.call(fn_adapt, build_criteria_args(fn_name, X, Y, "ALL", gh)),
                       error = function(e) NULL)
      rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(out)) next
      stops <- list(
        R2   = derive_stop_k(out$r_squared_adj, TRUE),
        BIC  = derive_stop_k(out$bic_values,    FALSE),
        eBIC = derive_stop_k(out$ebic_values,   FALSE))
      for (crit in shape_criteria) {
        Ks <- stops[[crit]]
        idx <- if (Ks > 0L) out$indices[seq_len(min(Ks, length(out$indices)))]
               else integer(0L)
        results[[paste0(fn_name, "_", crit)]] <- c(
          method = fn_name, criterion = crit,
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = rt / 3,
          evaluate_selection(as.integer(idx), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(idx)))
      }
      # Posterior thresholds derived from post.probs
      for (thr in post_thr) {
        Kp <- if (length(out$post.probs)) sum(out$post.probs >= thr) else 0L
        idx <- if (Kp > 0L) out$indices[seq_len(min(Kp, length(out$indices)))]
               else integer(0L)
        results[[paste0(fn_name, "_Post_", thr)]] <- c(
          method = fn_name, criterion = paste0("Post_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = rt / length(post_thr),
          evaluate_selection(as.integer(idx), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(idx)))
      }
    }
  }

  if (!isTRUE(cfg$skip_externals)) {
  # External comparators
  ext_label <- if (cfg$K_mode == "fixed") "K_fixed" else "natural_adaptive"
  K_susie <- if (cfg$K_mode == "fixed") K else NULL
  s_unif <- tryCatch(run_susie(X, Y, K = K_susie, prior_p = NULL), error = function(e) NULL)
  if (!is.null(s_unif))
    results[["SuSiE_unif"]] <- c(method="SuSiE_unif", criterion=ext_label,
       b=b, n=n, rho_intra=cfg$rho_intra, runtime=s_unif$runtime,
       evaluate_selection(s_unif$indices, causal_idx, m),
       indices_csv = indices_to_csv(s_unif$indices))
  s_pr <- tryCatch(run_susie(X, Y, K = K_susie, prior_p = p_j), error = function(e) NULL)
  if (!is.null(s_pr))
    results[["SuSiE_prior"]] <- c(method="SuSiE_prior", criterion=ext_label,
       b=b, n=n, rho_intra=cfg$rho_intra, runtime=s_pr$runtime,
       evaluate_selection(s_pr$indices, causal_idx, m),
       indices_csv = indices_to_csv(s_pr$indices))

  # varbvs (always with prior; 3 PIP thresholds in adaptive)
  if (cfg$K_mode == "fixed") {
    v <- tryCatch(run_varbvs(X, Y, K = K, prior_p = p_j), error = function(e) NULL)
    if (!is.null(v))
      results[["varbvs"]] <- c(method="varbvs", criterion="K_fixed",
        b=b, n=n, rho_intra=cfg$rho_intra, runtime=v$runtime,
        evaluate_selection(v$indices, causal_idx, m),
        indices_csv = indices_to_csv(v$indices))
  } else {
    vf <- tryCatch(run_varbvs(X, Y, K = NULL, prior_p = p_j, pip_threshold = 1e-12),
                    error = function(e) NULL)
    if (!is.null(vf)) {
      pips <- vf$post.probs
      for (thr in post_thr) {
        sel <- which(pips > thr)
        results[[paste0("varbvs_PIP_", thr)]] <- c(
          indices_csv = indices_to_csv(as.integer(sel)),
          method = "varbvs", criterion = paste0("PIP_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = vf$runtime / length(post_thr),
          evaluate_selection(as.integer(sel), causal_idx, m))
      }
    }
  }

  # LASSO
  if (cfg$K_mode == "fixed") {
    lK <- tryCatch(run_lasso(X, Y, K = K, lambda_rule = "1se"), error = function(e) NULL)
    if (!is.null(lK))
      results[["LASSO_topK"]] <- c(method="LASSO_topK", criterion="K_fixed",
        b=b, n=n, rho_intra=cfg$rho_intra, runtime=lK$runtime,
        evaluate_selection(lK$indices, causal_idx, m),
        indices_csv = indices_to_csv(lK$indices))
  }
  lf <- tryCatch(run_lasso(X, Y, K = NULL, lambda_rule = "1se"), error = function(e) NULL)
  if (!is.null(lf))
    results[["LASSO_1se"]] <- c(method="LASSO_1se", criterion=ext_label,
      b=b, n=n, rho_intra=cfg$rho_intra, runtime=lf$runtime,
      evaluate_selection(lf$indices, causal_idx, m),
      indices_csv = indices_to_csv(lf$indices))

  # CARMA -- locus regime natural
  ca <- tryCatch(run_carma(X, Y, K = if (cfg$K_mode == "fixed") K else NULL,
                            pip_threshold = 1e-12),
                  error = function(e) NULL)
  if (!is.null(ca)) {
    if (cfg$K_mode == "fixed") {
      results[["CARMA"]] <- c(method="CARMA", criterion="K_fixed",
         b=b, n=n, rho_intra=cfg$rho_intra, runtime=ca$runtime,
         evaluate_selection(ca$indices, causal_idx, m),
         indices_csv = indices_to_csv(ca$indices))
    } else {
      for (thr in post_thr) {
        sel <- which(ca$post.probs > thr)
        results[[paste0("CARMA_PIP_", thr)]] <- c(
          method = "CARMA", criterion = paste0("PIP_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = ca$runtime / length(post_thr),
          evaluate_selection(as.integer(sel), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(sel)))
      }
    }
  }

  # BayesC (BGLR) -- spike-and-slab Bayesian MCMC counterpart of NaM/ExM and
  # of varbvs (variational). Returns true PIPs via $ETA[[1]]$d.
  bc <- tryCatch(run_bayesc(X, Y, K = if (cfg$K_mode == "fixed") K else NULL,
                            pip_threshold = 1e-12,
                            nIter = 1500L, burnIn = 500L),
                 error = function(e) NULL)
  if (!is.null(bc)) {
    if (cfg$K_mode == "fixed") {
      results[["BayesC"]] <- c(method = "BayesC", criterion = "K_fixed",
         b = b, n = n, rho_intra = cfg$rho_intra, runtime = bc$runtime,
         evaluate_selection(bc$indices, causal_idx, m),
         indices_csv = indices_to_csv(bc$indices))
    } else {
      for (thr in post_thr) {
        sel <- which(bc$post.probs > thr)
        results[[paste0("BayesC_PIP_", thr)]] <- c(
          method = "BayesC", criterion = paste0("PIP_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = bc$runtime / length(post_thr),
          evaluate_selection(as.integer(sel), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(sel)))
      }
    }
  }

  # FINEMAP -- locus regime natural
  fm <- tryCatch(run_finemap(X, Y, K = if (cfg$K_mode == "fixed") K else NULL,
                              max_causal = 5L, pip_threshold = 1e-12),
                  error = function(e) NULL)
  if (!is.null(fm)) {
    if (cfg$K_mode == "fixed") {
      results[["FINEMAP"]] <- c(method="FINEMAP", criterion="K_fixed",
         b=b, n=n, rho_intra=cfg$rho_intra, runtime=fm$runtime,
         evaluate_selection(fm$indices, causal_idx, m),
         indices_csv = indices_to_csv(fm$indices))
    } else {
      for (thr in post_thr) {
        sel <- which(fm$post.probs > thr)
        results[[paste0("FINEMAP_PIP_", thr)]] <- c(
          method = "FINEMAP", criterion = paste0("PIP_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = fm$runtime / length(post_thr),
          evaluate_selection(as.integer(sel), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(sel)))
      }
    }
  }

  }   # end if (!skip_externals)
  }   # end if (!only_carma_prior)

  # CARMA prior-aware -- same locus regime, but supplied the scenario annotation
  # (binary indicator of prior-boosted SNPs) through CARMA's native w.list.
  # Placed last so it does not perturb the RNG stream of the methods above
  # (their results stay bit-identical to the prior-free run).
  cap <- if (isTRUE(cfg$skip_externals)) NULL else
         tryCatch(run_carma(X, Y, K = if (cfg$K_mode == "fixed") K else NULL,
                            pip_threshold = 1e-12, prior_p = p_j),
                  error = function(e) NULL)
  if (!is.null(cap)) {
    if (cfg$K_mode == "fixed") {
      results[["CARMA_prior"]] <- c(method="CARMA_prior", criterion="K_fixed",
         b=b, n=n, rho_intra=cfg$rho_intra, runtime=cap$runtime,
         evaluate_selection(cap$indices, causal_idx, m),
         indices_csv = indices_to_csv(cap$indices))
    } else {
      for (thr in post_thr) {
        sel <- which(cap$post.probs > thr)
        results[[paste0("CARMA_prior_PIP_", thr)]] <- c(
          method = "CARMA_prior", criterion = paste0("PIP_", thr),
          b = b, n = n, rho_intra = cfg$rho_intra,
          runtime = cap$runtime / length(post_thr),
          evaluate_selection(as.integer(sel), causal_idx, m),
          indices_csv = indices_to_csv(as.integer(sel)))
      }
    }
  }

  do.call(rbind, lapply(results, function(r) {
    fields <- c("method","criterion","b","n","rho_intra",
                "M","K","precision","recall","F1","FDR","runtime",
                "indices_csv")
    # Tolerate older code paths that don't set indices_csv.
    if (is.null(r[["indices_csv"]])) r[["indices_csv"]] <- ""
    as.data.frame(r[fields], stringsAsFactors = FALSE)
  }))
}

# ---------- 6. Parallel benchmark -------------------------------------------
use_fork <- (.Platform$OS.type == "unix")
cl <- if (use_fork) {
  parallel::makeForkCluster(cfg$n_cores)
} else {
  parallel::makeCluster(cfg$n_cores, type = "PSOCK")
}
registerDoParallel(cl); on.exit(stopCluster(cl), add = TRUE)

grid <- expand.grid(b = seq_len(cfg$B), n = n_values)
cat(sprintf("\nRunning %d reps x %d n values on %d cores (locus m=%d)\n",
            cfg$B, length(n_values), cfg$n_cores, cfg$m_locus))
t_start <- Sys.time()

all_results <- foreach(
  i = seq_len(nrow(grid)),
  .combine = rbind, .errorhandling = "remove",
  .packages = c("MASS","stats","Matrix","mvtnorm","susieR","varbvs","glmnet","statmod")
) %dopar% {
  invisible(gc(verbose = FALSE))
  source("Generate_XY_realistic.R", local = TRUE)
  source("External_Comparators.R", local = TRUE)
  source("Test_statistics_functions_optimized.R", local = TRUE)
  source("Test_statistics_functions_Criteria_optimized.R", local = TRUE)
  one_replication(grid$b[i], grid$n[i])
}
cat(sprintf("Done in %.1f seconds.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

# ---------- 7. Aggregate ----------------------------------------------------
summary_df <- all_results %>%
  mutate(M = as.numeric(M), K = as.numeric(K),
         precision = as.numeric(precision), recall = as.numeric(recall),
         F1 = as.numeric(F1), FDR = as.numeric(FDR),
         runtime = as.numeric(runtime), n = as.integer(n)) %>%
  group_by(n, method, criterion) %>%
  summarise(
    n_reps    = dplyr::n(),
    mean_M    = mean(M, na.rm = TRUE),
    mean_K    = mean(K, na.rm = TRUE),
    precision = mean(precision, na.rm = TRUE),
    recall    = mean(recall, na.rm = TRUE),
    F1        = mean(F1, na.rm = TRUE),
    FDR       = mean(FDR, na.rm = TRUE),
    runtime_s = mean(runtime, na.rm = TRUE),
    .groups   = "drop"
  ) %>% arrange(n, method, criterion)

cat("\n=========== AGGREGATED RESULTS (locus sim) ===========\n")
print(as.data.frame(summary_df), row.names = FALSE, digits = 3)

saveRDS(list(config = cfg, summary = summary_df, raw = all_results),
        cfg$outfile)
cat(sprintf("\nResults written to: %s\n", cfg$outfile))
