# ==============================================================================
# File: Benchmark_External.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Description:
#   Reproducible evaluation of the six proposed in-house methods
#       (NaW, AsW, ExW, AsS, NaM, ExM)
#   against external comparators (SuSiE, varbvs, LASSO, BayesC) on a synthetic
#   genotype matrix with a block-AR(1) LD structure (Generate_XY_realistic.R).
#
#   For each replication:
#     1. Generate an n x m genotype matrix with intra-block correlation rho.
#     2. Place q = 15 causal SNPs (placement A: one at the centre of each of
#        the first 15 blocks; see choose_causal_indices for the variants).
#     3. Simulate Y = X beta + eps with effect sizes lambda calibrated to
#        match a marginal power gamma.
#     4. Apply the selection methods at top-K (or with adaptive stopping) and
#        compute the metrics.
#
#   Metrics: M, K, precision, recall, F1, FDR, runtime (seconds).
#
# Usage:
#   Rscript Benchmark_External.R \
#       --B 50 --scenario informative --n_grid "500,1000" \
#       --rho_grid "0.8,0.95,0.98" --multiple_exact TRUE
# ==============================================================================

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(methods)
})

# ------------------------------------------------------------------------------
# 0. Project sources
# ------------------------------------------------------------------------------
source("Generate_XY_realistic.R")                           # synthetic AR(1) X
source("Test_statistics_functions_optimized.R")             # in-house methods (fixed-K)
source("Test_statistics_functions_Criteria_optimized.R")    # in-house methods (adaptive K via stopping rule)
source("External_Comparators.R")                            # SuSiE / varbvs / LASSO

# ------------------------------------------------------------------------------
# 1. Configuration (override via CLI args)
# ------------------------------------------------------------------------------
config_default <- list(
  # ----- Design (kept in the output for compatibility; synthetic AR(1) only) -----
  design             = "synthetic",

  n_grid             = "500,1000",                # subsample / individuals sizes
  gamma_grid         = "",                          # comma-separated; overrides --gamma
  tau_grid           = "",                          # comma-separated; overrides --tau
  n_causal           = 15L,
  placement          = "A",                        # "A" / "B" / "random"
  gamma              = 0.6,                        # target marginal power
  tau                = 0.20,
  K_fixed            = 15L,                        # top-K used in fixed mode
  K_mode             = "fixed",                    # "fixed" or "adaptive"
  K_threshold        = 30L,                        # max iterations cap for adaptive in-house.
                                                   # Set to 30 (= 2*q with q=15 causals) which
                                                   # is sufficient for BIC/eBIC/Posterior to
                                                   # stop naturally while keeping wall time
                                                   # ~40% shorter than the previous default 50.
  target_prob        = 0.95,                       # threshold for "Posterior" stopping rule
  scenario           = "informative",              # informative / uniform / misspecified / mixed
  B                  = 50L,                        # set >= 500 for production
  S_mc               = 100L,
  multiple_exact     = FALSE,                      # TRUE: NaM/ExM prior integrals computed exactly (S_mc unused)
  seed               = 2026L,
  # Memory-conservative default: cap at 4 cores. Each worker keeps its own
  # copy of the n x m design matrix plus the allocations made by
  # susie/varbvs/glmnet (~2-3x X), so many workers can require several GB.
  # Override with --n_cores when you have abundant RAM.
  n_cores            = min(parallel::detectCores() - 1L, 4L),
  outfile            = "results_benchmark_external.rds",

  # ----- Method selection (run subset of comparators) -----
  # At m ~ 5000, varbvs's variational EM can take hours per fit (the
  # mean-field landscape has many quasi-equivalent local maxima). The flags
  # below let the user split work:
  #   - --skip_varbvs   : skip varbvs only (run in-house + SuSiE + LASSO)
  #   - --only_varbvs   : run varbvs only (skip in-house + SuSiE + LASSO)
  #   - --skip_susie    : skip SuSiE
  #   - --skip_lasso    : skip LASSO
  # The expected usage on a memory-bound machine: run --skip_varbvs on the
  # main server and --only_varbvs separately for the same scenario; merge the
  # resulting RDS at table-generation time.
  skip_varbvs        = FALSE,
  only_varbvs        = FALSE,
  only_bayesc        = FALSE,                      # run ONLY BayesC (skip in-house + SuSiE +
                                                   # varbvs + LASSO); BayesC is prior-free, so a
                                                   # single run fills every prior scenario's table.
  skip_susie         = FALSE,
  susie_L            = 10L,                        # susieR L (number of effects); must be >= the
                                                   # number of simulated causal variants for a fair
                                                   # fixed-K comparison (q = 15 -> use L >= 15)
  skip_lasso         = FALSE,
  skip_bayesc        = FALSE,                      # BayesC (BGLR) is the Bayesian MCMC sibling
                                                   # of NaM/ExM and varbvs. ~3-5 min per fit at
                                                   # m ~ 3000; set TRUE to bypass when compute-bound.

  # Comma-separated subset of in-house methods to run (long or short names,
  # e.g. "AsS,NaM,ExM"); empty = all six.
  inhouse_only       = "",

  # ----- Synthetic AR(1) design -----
  causals_per_window = 3L,                        # causals per block (placement B)
  m_synth            = 5000L,                     # # SNPs for the synthetic design
  rho_grid           = "0.8,0.95,0.98",            # intra-block LD levels (synthetic only)
  rho_inter          = 0.0,
  block_size         = 50L,
  maf_lo             = 0.05,
  maf_hi             = 0.50,
  maf_distribution   = "uniform"
)

args <- commandArgs(trailingOnly = TRUE)
parse_cli <- function(args, defaults) {
  i <- 1L
  while (i <= length(args) - 1L) {
    if (startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      val <- args[i + 1L]
      if (!is.null(defaults[[key]])) {
        defaults[[key]] <- methods::as(val, class(defaults[[key]]))
      }
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  defaults
}
cfg <- parse_cli(args, config_default)
if (cfg$design != "synthetic") stop("Only the synthetic AR(1) design is supported (--design synthetic).")
stopifnot(cfg$K_mode %in% c("fixed", "adaptive"))
n_values <- as.integer(strsplit(cfg$n_grid, ",", fixed = TRUE)[[1L]])
stopifnot(all(is.finite(n_values)), all(n_values >= 50L))

# Parse gamma/tau grids; fall back to the scalar cfg$gamma / cfg$tau when
# the corresponding grid is not provided. This keeps the single-(γ,τ) usage
# of the script backward-compatible.
gamma_values <- if (is.null(cfg[["gamma_grid"]]) || !nzchar(cfg[["gamma_grid"]])) {
  as.numeric(cfg[["gamma"]])
} else {
  as.numeric(strsplit(cfg[["gamma_grid"]], ",", fixed = TRUE)[[1L]])
}
tau_values <- if (is.null(cfg[["tau_grid"]]) || !nzchar(cfg[["tau_grid"]])) {
  as.numeric(cfg[["tau"]])
} else {
  as.numeric(strsplit(cfg[["tau_grid"]], ",", fixed = TRUE)[[1L]])
}
stopifnot(all(is.finite(gamma_values)), all(gamma_values > 0 & gamma_values < 1))
stopifnot(all(is.finite(tau_values)), all(tau_values > 0))

# Intra-block LD levels
rho_values <- as.numeric(strsplit(cfg$rho_grid, ",", fixed = TRUE)[[1L]])
stopifnot(all(is.finite(rho_values)), all(abs(rho_values) < 1))

cat("Configuration:\n"); str(cfg)
cat("n grid:", n_values, "\n")
cat("rho_intra grid:", rho_values, "\n")

# ------------------------------------------------------------------------------
# 2. Causal SNPs of the synthetic block layout
# ------------------------------------------------------------------------------
m_eff <- cfg$m_synth
causal_idx <- choose_causal_indices(
  m                 = m_eff,
  n_causal          = cfg$n_causal,
  block_size        = cfg$block_size,
  placement         = cfg$placement,
  causals_per_block = cfg$causals_per_window,
  seed              = cfg$seed
)
cat(sprintf("\nCausal SNP indices (synthetic, placement %s, %d SNPs):\n",
            cfg$placement, length(causal_idx)))
print(causal_idx)

# ------------------------------------------------------------------------------
# 4. Prior inclusion probability vector — scenario-dependent
# ------------------------------------------------------------------------------
# In simulation, true causals are known so we may use them in the "informative"
# prior. The values are kept at the simulation-oracle calibration (0.2/0.5/0.8);
# the real-data application uses the milder scale 0.3/0.5/0.7. target_set
# (optional) lists SNPs given 0.5 under the informative prior and used as decoys
# under the misspecified prior; it is empty for the synthetic design.
build_prior <- function(scenario, m, causal_idx, target_set = integer(0)) {
  p_j <- switch(
    scenario,
    informative   = {
      v <- rep(0.2, m)
      if (length(target_set) > 0L) v[target_set] <- 0.5
      v[causal_idx] <- 0.8
      v
    },
    uniform       = rep(0.5, m),
    misspecified  = {
      v <- rep(0.5, m)
      v[causal_idx] <- 0.2
      decoy <- setdiff(target_set, causal_idx)
      n_decoy <- min(length(causal_idx), length(decoy))
      if (n_decoy > 0L) v[sample(decoy, n_decoy)] <- 0.8
      v
    },
    mixed         = {
      v <- rep(0.5, m)
      n_c <- length(causal_idx)
      lo  <- causal_idx[seq_len(min(5L, n_c))]
      mid <- causal_idx[seq(6L, min(10L, n_c))]
      hi  <- causal_idx[seq(11L, n_c)]
      v[lo]  <- 0.2
      if (length(mid) > 0L) v[mid] <- 0.5
      if (length(hi)  > 0L) v[hi]  <- 0.8
      v
    },
    stop("Unknown scenario: ", scenario)
  )
  p_j
}
p_j <- build_prior(cfg$scenario, m_eff, causal_idx)

# ------------------------------------------------------------------------------
# 5. Effect-size calibration to target power gamma
# ------------------------------------------------------------------------------
calibrate_lambda <- function(n, gamma, alpha = 0.05) {
  (qnorm(1 - alpha / 2) + qnorm(gamma)) / sqrt(n)
}
cat("gamma grid:", gamma_values, "\n")
cat("tau   grid:", tau_values,   "\n")
for (g_v in gamma_values) for (n_v in n_values) {
  cat(sprintf("Calibrated lambda for n=%d, gamma=%.2f: %.4f\n",
              n_v, g_v, calibrate_lambda(n_v, g_v)))
}

# ------------------------------------------------------------------------------
# 6. Single replication function (parameterised by n and rho_intra)
# ------------------------------------------------------------------------------
one_replication <- function(b, n, rho_intra, gamma_v, tau_v,
                            cfg, p_j, causal_idx) {
  set.seed(cfg$seed + b * 1000L + n * 7L +
             ifelse(is.finite(rho_intra), as.integer(rho_intra * 1e4), 0L) +
             as.integer(gamma_v * 1e4) + as.integer(tau_v * 1e4) * 13L)
  lambda_true <- calibrate_lambda(n, gamma_v)

  # Synthetic AR(1) genotype matrix at this rho_intra
  X_sub <- Generate.X.realistic(
    n                = n,
    m                = cfg$m_synth,
    ld_structure     = "blocks",
    rho_intra        = rho_intra,
    rho_inter        = cfg$rho_inter,
    block_size       = cfg$block_size,
    maf_lo           = cfg$maf_lo,
    maf_hi           = cfg$maf_hi,
    maf_distribution = cfg$maf_distribution
  )
  beta <- rep(0, cfg$m_synth); beta[causal_idx] <- lambda_true
  Y <- Generate.Y(n = n, X = X_sub, beta = beta, sigma = 1)
  m <- cfg$m_synth

  K <- cfg$K_fixed
  results <- list()

  # 3. In-house methods — skip entirely if --only_varbvs (the varbvs-only run
  # exists for the sole purpose of producing the varbvs rows of the table).
  in_house <- if (cfg$only_varbvs || cfg$only_bayesc) character(0L) else
              c("Naive_Wald", "Asymptotic_Wald", "Exact_Wald",
                "Asymptotic_Score", "Naive_Multiple", "Exact_Multiple")
  if (nzchar(cfg$inhouse_only)) {
    short <- c(NaW = "Naive_Wald", AsW = "Asymptotic_Wald", ExW = "Exact_Wald",
               AsS = "Asymptotic_Score", NaM = "Naive_Multiple", ExM = "Exact_Multiple")
    sel <- strsplit(cfg$inhouse_only, ",", fixed = TRUE)[[1L]]
    sel <- ifelse(sel %in% names(short), short[sel], sel)
    in_house <- intersect(in_house, sel)
  }

  # Pre-compute Gauss-Hermite nodes (used by Asymptotic/Exact Wald, Score).
  # Computing once per replication is much cheaper than recomputing for each
  # method (~50 us).
  gh <- make_quad_nodes(30L, tau_v)

  # Build the argument list specific to each method signature.
  # - Naive_Wald: (X, Y, a.priori, K, tau)
  # - Asymptotic_Wald, Exact_Wald, Asymptotic_Score: (..., tau, gh)
  # - Naive_Multiple, Exact_Multiple: (..., tau2, S)
  build_inhouse_args <- function(fn_name) {
    base <- list(X = X_sub, Y = Y, a.priori = p_j, K = K)
    extras <- switch(
      fn_name,
      Naive_Wald        = list(tau = tau_v),
      Asymptotic_Wald   = list(tau = tau_v, gh = gh),
      Exact_Wald        = list(tau = tau_v, gh = gh),
      Asymptotic_Score  = list(tau = tau_v, gh = gh),
      Naive_Multiple    = list(tau2 = tau_v^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
      Exact_Multiple    = list(tau2 = tau_v^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
      list(tau = tau_v)
    )
    c(base, extras)
  }

  if (cfg$K_mode == "fixed") {
    # ---- Fixed K=K_fixed: 6 in-house rows --------------------------------
    for (fn_name in in_house) {
      if (!exists(fn_name)) next
      t0 <- Sys.time()
      out <- tryCatch(
        do.call(fn_name, build_inhouse_args(fn_name)),
        error = function(e) NULL
      )
      rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(out)) next
      met <- evaluate_selection(out$indices, causal_idx, m)
      results[[fn_name]] <- c(method = fn_name, criterion = "K_fixed",
                              b = b, n = n, rho_intra = rho_intra,
                              runtime = rt, met)
    }
  } else {
    # ---- Adaptive K via stopping criterion: 6 methods x (3 + 3) criteria ----
    # Three "shape" criteria + three Posterior thresholds (0.5, 0.95, 0.99)
    posterior_thresholds <- c(0.5, 0.95, 0.99)
    shape_criteria <- c("R2", "BIC", "eBIC")

    # Build the adaptive (Criteria) argument list per method signature
    build_criteria_args <- function(fn_name, crit, target_prob = cfg$target_prob) {
      base <- list(X = X_sub, Y = Y, a.priori = p_j,
                   threshold = cfg$K_threshold,
                   target_prob = target_prob,
                   criteria = crit)
      extras <- switch(
        fn_name,
        Naive_Wald        = list(tau = tau_v),
        Asymptotic_Wald   = list(tau = tau_v, gh = gh),
        Exact_Wald        = list(tau = tau_v, gh = gh),
        Asymptotic_Score  = list(tau = tau_v, gh = gh),
        Naive_Multiple    = list(tau2 = tau_v^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
        Exact_Multiple    = list(tau2 = tau_v^2, S = (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
        list(tau = tau_v)
      )
      c(base, extras)
    }

    # Shape criteria (R2, BIC, eBIC) --- OPTIMISED:
    # for each in-house method we make ONE adaptive call with criteria="ALL"
    # (iterates up to K_threshold without early-stop) and then derive the three
    # stopping points post-hoc from the returned histories. Saves ~3x on this
    # block since the three shape criteria differ only in when they break.
    derive_stop_k <- function(hist, want_increasing) {
      # First k where the metric stops improving. For R^2 (want_increasing=TRUE)
      # we stop at the last k before a decrease; for BIC/eBIC (FALSE), at the
      # last k before an increase. Mirrors check_stop()'s "first non-improvement"
      # rule and yields indices.inclus[1:K_stop] as the selected set.
      if (length(hist) <= 1L) return(length(hist))
      diffs <- diff(hist)
      bad <- if (want_increasing) which(diffs < 0) else which(diffs > 0)
      if (length(bad) == 0L) length(hist) else bad[1L]
    }
    for (fn_name in in_house) {
      fn_adaptive <- paste0(fn_name, "_Criteria")
      if (!exists(fn_adaptive)) next
      t0 <- Sys.time()
      out <- tryCatch(
        do.call(fn_adaptive, build_criteria_args(fn_name, "ALL")),
        error = function(e) NULL
      )
      rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(out)) next
      rt_share <- rt / 3
      stops <- list(
        R2   = derive_stop_k(out$r_squared_adj, want_increasing = TRUE),
        BIC  = derive_stop_k(out$bic_values,    want_increasing = FALSE),
        eBIC = derive_stop_k(out$ebic_values,   want_increasing = FALSE)
      )
      for (crit in shape_criteria) {
        K_stop <- stops[[crit]]
        idx <- if (K_stop > 0L) out$indices[seq_len(min(K_stop, length(out$indices)))]
               else integer(0L)
        results[[paste0(fn_name, "_", crit)]] <- c(
          method = fn_name, criterion = crit,
          b = b, n = n, rho_intra = rho_intra,
          runtime = rt_share,
          evaluate_selection(as.integer(idx), causal_idx, m)
        )
      }
    }

    # Posterior threshold criterion at three values --- OPTIMISED:
    # for each in-house method we make ONE fixed-K call up to K_threshold and
    # then derive the three Posterior stopping points post-hoc from the
    # trajectory of posterior probabilities. Saves ~3x on this block.
    for (fn_name in in_house) {
      if (!exists(fn_name)) next
      t0 <- Sys.time()
      # Reuse build_inhouse_args but override K to K_threshold
      args_post <- build_inhouse_args(fn_name)
      args_post$K <- cfg$K_threshold
      full_run <- tryCatch(do.call(fn_name, args_post),
                           error = function(e) NULL)
      rt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (is.null(full_run)) next
      # full_run$indices : K_threshold indices in order of selection
      # full_run$post.probs : posterior prob of the selected SNP at each step
      idx_traj <- as.integer(full_run$indices)
      prob_traj <- as.numeric(full_run$post.probs)
      # Per-threshold rows share the cost of the single full_run call
      rt_share <- rt / length(posterior_thresholds)
      for (thr in posterior_thresholds) {
        first_below <- which(prob_traj < thr)
        stop_step <- if (length(first_below) == 0L) {
          length(idx_traj)
        } else {
          first_below[1L] - 1L
        }
        sel <- if (stop_step >= 1L) idx_traj[seq_len(stop_step)] else integer(0)
        crit_label <- sprintf("Post_%g", thr)
        met <- evaluate_selection(sel, causal_idx, m)
        results[[paste0(fn_name, "_", crit_label)]] <- c(
          method = fn_name, criterion = crit_label,
          b = b, n = n, rho_intra = rho_intra,
          runtime = rt_share, met
        )
      }
    }
  }

  # 4. External comparators — adaptive setting determined by K_mode
  #    In fixed mode: top-K=K_fixed across the board.
  #    In adaptive mode: SuSiE -> 95% credible sets, varbvs -> PIP > 0.5,
  #                      LASSO -> lambda.1se nonzero set.
  ext_label <- if (cfg$K_mode == "fixed") "K_fixed" else "natural_adaptive"

  # SuSiE (uniform + informative prior) -- skip if cfg$only_varbvs or cfg$skip_susie
  K_susie <- if (cfg$K_mode == "fixed") K else NULL
  if (!cfg$skip_susie && !cfg$only_varbvs && !cfg$only_bayesc) {
    s_unif <- tryCatch(run_susie(X_sub, Y, K = K_susie, L = cfg$susie_L, prior_p = NULL),
                       error = function(e) NULL)
    if (!is.null(s_unif)) {
      results[["SuSiE_unif"]] <- c(
        method = "SuSiE_unif", criterion = ext_label,
        b = b, n = n, rho_intra = rho_intra,
        runtime = s_unif$runtime,
        evaluate_selection(s_unif$indices, causal_idx, m)
      )
    }
    s_info <- tryCatch(run_susie(X_sub, Y, K = K_susie, L = cfg$susie_L, prior_p = p_j),
                       error = function(e) NULL)
    if (!is.null(s_info)) {
      results[["SuSiE_prior"]] <- c(
        method = "SuSiE_prior", criterion = ext_label,
        b = b, n = n, rho_intra = rho_intra,
        runtime = s_info$runtime,
        evaluate_selection(s_info$indices, causal_idx, m)
      )
    }
  }

  # varbvs — skip if cfg$skip_varbvs (may run separately). By design, varbvs is a Bayesian method whose
  # prior inclusion probabilities are first-class arguments; evaluating it
  # under an artificial uniform prior would handicap it relative to its
  # intended use, so we report only the scenario-prior variant.
  if (!cfg$skip_varbvs && !cfg$only_bayesc) {
  if (cfg$K_mode == "fixed") {
    v <- tryCatch(run_varbvs(X_sub, Y, K = K, prior_p = p_j),
                  error = function(e) NULL)
    if (!is.null(v)) {
      results[["varbvs"]] <- c(
        method = "varbvs", criterion = "K_fixed",
        b = b, n = n, rho_intra = rho_intra,
        runtime = v$runtime,
        evaluate_selection(v$indices, causal_idx, m)
      )
    }
  } else {
    # OPTIMISED: a single varbvs fit produces the full PIP vector; we then
    # apply the three thresholds post-hoc rather than refit 3 times.
    v_full <- tryCatch(
      run_varbvs(X_sub, Y, K = NULL, prior_p = p_j, pip_threshold = 1e-12),
      error = function(e) NULL)
    if (!is.null(v_full)) {
      pips <- v_full$post.probs
      rt_share <- v_full$runtime / 3
      for (thr in c(0.5, 0.95, 0.99)) {
        crit_label <- sprintf("PIP_%g", thr)
        sel <- which(pips > thr)
        results[[paste0("varbvs_", crit_label)]] <- c(
          method = "varbvs", criterion = crit_label,
          b = b, n = n, rho_intra = rho_intra,
          runtime = rt_share,
          evaluate_selection(as.integer(sel), causal_idx, m)
        )
      }
    }
  }
  }  # end if (!cfg$skip_varbvs)

  # LASSO — skip if cfg$only_varbvs or cfg$skip_lasso.
  # Top-K is only reported in fixed mode (in adaptive we use 1se only).
  if (!cfg$skip_lasso && !cfg$only_varbvs && !cfg$only_bayesc) {
    if (cfg$K_mode == "fixed") {
      l_K <- tryCatch(run_lasso(X_sub, Y, K = K, lambda_rule = "1se"),
                      error = function(e) NULL)
      if (!is.null(l_K)) {
        results[["LASSO_topK"]] <- c(
          method = "LASSO_topK", criterion = "K_fixed",
          b = b, n = n, rho_intra = rho_intra,
          runtime = l_K$runtime,
          evaluate_selection(l_K$indices, causal_idx, m)
        )
      }
    }
    l_full <- tryCatch(run_lasso(X_sub, Y, K = NULL, lambda_rule = "1se"),
                       error = function(e) NULL)
    if (!is.null(l_full)) {
      results[["LASSO_1se"]] <- c(
        method = "LASSO_1se", criterion = ext_label,
        b = b, n = n, rho_intra = rho_intra,
        runtime = l_full$runtime,
        evaluate_selection(l_full$indices, causal_idx, m)
      )
    }
  }

  # BayesC (BGLR spike-and-slab MCMC) -- Bayesian sibling of NaM/ExM and
  # varbvs. Skipped if cfg$skip_bayesc or cfg$only_varbvs. In fixed mode
  # we keep the top-K by PIP; in adaptive mode we report the three PIP
  # thresholds {0.5, 0.95, 0.99} matched to the in-house Posterior ladder.
  if (!cfg$skip_bayesc && !cfg$only_varbvs) {
    bc <- tryCatch(
      run_bayesc(X_sub, Y, K = NULL,
                 pip_threshold = 1e-12,
                 nIter = 1500L, burnIn = 500L),
      error = function(e) NULL
    )
    if (!is.null(bc)) {
      if (cfg$K_mode == "fixed") {
        sel_topK <- order(bc$post.probs, decreasing = TRUE)[seq_len(min(K, m))]
        results[["BayesC"]] <- c(
          method = "BayesC", criterion = "K_fixed",
          b = b, n = n, rho_intra = rho_intra,
          runtime = bc$runtime,
          evaluate_selection(as.integer(sel_topK), causal_idx, m)
        )
      } else {
        for (thr in c(0.5, 0.95, 0.99)) {
          sel <- which(bc$post.probs > thr)
          results[[paste0("BayesC_PIP_", thr)]] <- c(
            method = "BayesC", criterion = paste0("PIP_", thr),
            b = b, n = n, rho_intra = rho_intra,
            runtime = bc$runtime / 3,
            evaluate_selection(as.integer(sel), causal_idx, m)
          )
        }
      }
    }
  }

  do.call(rbind, lapply(results, function(r) {
    as.data.frame(r[c("method", "criterion", "b", "n", "rho_intra",
                      "M", "K", "precision", "recall", "F1", "FDR",
                      "runtime")],
                  stringsAsFactors = FALSE)
  }))
}

# ------------------------------------------------------------------------------
# 7. Parallel benchmark
# ------------------------------------------------------------------------------
cat(sprintf("\nRunning %d reps x %d n values x %d rho values on %d cores (design = %s)...\n",
            cfg$B, length(n_values), length(rho_values), cfg$n_cores, cfg$design))

# Cluster type: on Linux/macOS use forking (copy-on-write, MUCH cheaper
# memory-wise; the parent objects are shared rather than serialised per worker).
# On Windows (no fork available), fall back to PSOCK with explicit env export.
use_fork <- (.Platform$OS.type == "unix")
if (use_fork) {
  cat(sprintf("Using FORK cluster (%d workers, copy-on-write shared memory).\n",
              cfg$n_cores))
  # Must run before forking: workers then share the AR(1) factors
  # copy-on-write (Generate_XY_realistic.R, section 4b).
  ar1_prewarm_factors(cfg$m_synth, rho_values, cfg$rho_inter, cfg$block_size,
                      cores = cfg$n_cores)
  cl <- parallel::makeForkCluster(cfg$n_cores)
} else {
  cat(sprintf("Using PSOCK cluster (%d workers, per-worker memory).\n",
              cfg$n_cores))
  cl <- parallel::makeCluster(cfg$n_cores, type = "PSOCK")
}
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

t_start <- Sys.time()
grid <- expand.grid(b = seq_len(cfg$B), n = n_values, rho_intra = rho_values,
                    gamma = gamma_values, tau = tau_values)

all_results <- foreach(
  i = seq_len(nrow(grid)),
  .combine = rbind,
  .errorhandling = "remove",
  .packages = c("MASS", "stats", "Matrix", "BGLR", "mvtnorm",
                "susieR", "varbvs", "glmnet", "statmod")
) %dopar% {
  # Aggressive GC at start of each worker iteration to prevent memory
  # creep across the (n x m) genotype-matrix copies allocated by SuSiE
  # and varbvs during the cell run.
  invisible(gc(verbose = FALSE))
  source("Generate_XY_realistic.R",                  local = TRUE)
  source("Test_statistics_functions_optimized.R",    local = TRUE)
  source("External_Comparators.R",                   local = TRUE)
  res <- one_replication(grid$b[i], grid$n[i], grid$rho_intra[i],
                          grid$gamma[i], grid$tau[i],
                          cfg, p_j, causal_idx)
  if (!is.null(res) && nrow(res) > 0L) {
    res$gamma <- grid$gamma[i]
    res$tau   <- grid$tau[i]
  }
  res
}
# stopCluster is handled by the on.exit() registered above.

cat(sprintf("Done in %.1f seconds.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "secs"))))

# ------------------------------------------------------------------------------
# 8. Aggregate and report
# ------------------------------------------------------------------------------
group_cols <- c("n", "rho_intra", "gamma", "tau", "method", "criterion")

summary_df <- all_results %>%
  mutate(
    M = as.numeric(M), K = as.numeric(K),
    precision = as.numeric(precision), recall = as.numeric(recall),
    F1 = as.numeric(F1), FDR = as.numeric(FDR),
    runtime = as.numeric(runtime),
    n = as.integer(n),
    rho_intra = as.numeric(rho_intra),
    gamma = as.numeric(gamma),
    tau = as.numeric(tau)
  ) %>%
  group_by(across(all_of(group_cols))) %>%
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
  ) %>%
  arrange(across(all_of(group_cols)), desc(F1))

cat(sprintf("\n=========== AGGREGATED RESULTS (design = %s) ===========\n", cfg$design))
print(as.data.frame(summary_df), row.names = FALSE, digits = 3)

saveRDS(list(config       = cfg,
             design       = cfg$design,
             n_grid       = n_values,
             gamma_grid   = gamma_values,
             tau_grid     = tau_values,
             rho_grid     = rho_values,
             causal_idx   = causal_idx,
             causal_map   = NULL,
             all_results  = all_results,
             summary      = summary_df),
        cfg$outfile)
cat(sprintf("\nResults saved to: %s\n", cfg$outfile))
