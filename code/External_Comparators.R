# ==============================================================================
# File: External_Comparators.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe 
#
# Description:
#   Wrapper functions for three external fine-mapping / variable-selection
#   methods used as benchmarks in the revised manuscript:
#       * SuSiE   (susieR::susie)
#       * varbvs  (varbvs::varbvs, Carbonetto & Stephens 2012)
#       * LASSO   (glmnet::cv.glmnet, Tibshirani 1996)
#
#   Each wrapper returns a list with the same structure as the in-house
#   methods (NaW, AsW, ExW, AsS, NaM, ExM):
#
#       list(
#         indices    = integer vector of selected SNP indices (length <= K),
#         post.probs = numeric vector of PIPs (or proxy) over all m SNPs,
#         runtime    = elapsed wall-clock time in seconds
#       )
#
#   Plus a helper `evaluate_selection()` that returns standard metrics
#   (precision, recall, F1, FDR) for downstream comparison.
#
# Dependencies:
#   install.packages(c("susieR", "varbvs", "glmnet"))
# ==============================================================================

suppressPackageStartupMessages({
  library(susieR)
  library(varbvs)
  library(glmnet)
})

# ------------------------------------------------------------------------------
# 1. SuSiE wrapper
# ------------------------------------------------------------------------------

#' Run SuSiE on standardized X and continuous Y
#'
#' @param X         Matrix (n x m), centered & scaled (as produced by Generate.X).
#' @param Y         Numeric vector of length n.
#' @param K         Maximum number of selected SNPs (top-K by PIP). If NULL,
#'                  takes union of credible sets.
#' @param L         Maximum number of single-effect components (default 10).
#' @param prior_p  Optional prior inclusion probability vector (length m).
#'                  Passed to susie via `prior_weights`.
#' @param coverage  Credible-set coverage level (default 0.95).
#' @param ...       Passed to susieR::susie.
#'
#' @return list(indices, post.probs, runtime, cs).
run_susie <- function(X, Y, K = NULL, L = 10L,
                      prior_p = NULL, coverage = 0.95, ...) {
  t0 <- Sys.time()
  fit <- susieR::susie(
    X, as.numeric(Y),
    L              = L,
    prior_weights  = prior_p,
    coverage       = coverage,
    standardize    = FALSE,   # we already scaled in Generate.X
    intercept      = TRUE,
    ...
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  pip <- as.numeric(fit$pip)
  m   <- ncol(X)
  if (length(pip) != m) pip <- c(pip, rep(0, m - length(pip)))

  if (is.null(K)) {
    cs_list <- fit$sets$cs
    sel <- if (is.null(cs_list)) integer(0) else sort(unique(unlist(cs_list)))
  } else {
    sel <- order(pip, decreasing = TRUE)[seq_len(min(K, m))]
  }

  list(indices = as.integer(sel),
       post.probs = pip,
       runtime = runtime,
       cs = fit$sets)
}

# ------------------------------------------------------------------------------
# 2. varbvs wrapper
# ------------------------------------------------------------------------------

#' Run varbvs (variational Bayesian variable selection)
#'
#' @param X        Matrix (n x m).
#' @param Y        Numeric vector of length n.
#' @param K        Optional cap on selected SNPs (top-K by PIP).
#' @param prior_p Optional vector of prior inclusion probability $p_j$ in (0,1) for each SNP; the
#'                vector is internally converted to base-10 log-odds for varbvs.
#' @param ...      Passed to varbvs::varbvs.
#'
#' @param pip_threshold Posterior inclusion probability cut-off for selection
#'                     when K is NULL (default 0.5). Allowed to vary in
#'                     {0.5, 0.95, 0.99} for fairness with the in-house
#'                     Posterior stopping criterion.
#' @return list(indices, post.probs, runtime).
run_varbvs <- function(X, Y, K = NULL, prior_p = NULL,
                       pip_threshold = 0.5, ...) {
  t0 <- Sys.time()

  # varbvs uses logodds; convert if a prior is provided. CRITICAL: passing
  # `logodds = NULL` explicitly crashes varbvs (it expects either a numeric
  # vector or to be left unset to use the internal default grid). We therefore
  # build the argument list dynamically and only include `logodds` when we
  # actually have one.
  # Build args list. Any parameter passed via ... should take precedence
  # over our defaults (in particular `verbose`); we therefore inject
  # defaults only when not already supplied by the caller.
  extra <- list(...)
  varbvs_args <- list(
    X       = X,
    Z       = NULL,
    y       = as.numeric(Y),
    family  = "gaussian"
  )
  if (!"verbose" %in% names(extra)) varbvs_args$verbose <- FALSE
  varbvs_args <- c(varbvs_args, extra)
  if (!is.null(prior_p)) {
    prior_p <- pmin(pmax(prior_p, 1e-6), 1 - 1e-6)
    varbvs_args$logodds <- log(prior_p / (1 - prior_p)) / log(10)
  }

  fit <- do.call(varbvs::varbvs, varbvs_args)
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # PIPs averaged across logodds grid (or single value)
  pip <- if (is.matrix(fit$alpha)) {
    as.numeric(rowMeans(fit$alpha))
  } else {
    as.numeric(fit$alpha)
  }
  m <- ncol(X)
  if (length(pip) != m) pip <- c(pip, rep(0, m - length(pip)))

  if (is.null(K)) {
    stopifnot(pip_threshold > 0, pip_threshold < 1)
    sel <- which(pip > pip_threshold)
  } else {
    sel <- order(pip, decreasing = TRUE)[seq_len(min(K, m))]
  }

  list(indices = as.integer(sel),
       post.probs = pip,
       runtime = runtime)
}

# ------------------------------------------------------------------------------
# 3. LASSO wrapper (cv.glmnet)
# ------------------------------------------------------------------------------

#' Run LASSO with cross-validated lambda
#'
#' @param X            Matrix (n x m).
#' @param Y            Numeric vector of length n.
#' @param K            Optional cap on selected SNPs (top-K by |coef|). If NULL,
#'                     returns all non-zero coefficients at lambda.1se.
#' @param lambda_rule  One of "min" or "1se" (default).
#' @param nfolds       CV folds (default 10).
#' @param ...          Passed to glmnet::cv.glmnet.
#'
#' @return list(indices, post.probs, runtime). Note: post.probs is the
#'         absolute coefficient (proxy for "evidence"), not a true PIP.
run_lasso <- function(X, Y, K = NULL,
                      lambda_rule = c("1se", "min"),
                      nfolds = 10L, ...) {
  lambda_rule <- match.arg(lambda_rule)
  t0 <- Sys.time()
  fit <- glmnet::cv.glmnet(
    X, as.numeric(Y),
    family = "gaussian",
    alpha  = 1,           # LASSO
    nfolds = nfolds,
    standardize = FALSE,
    ...
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  lam <- if (lambda_rule == "1se") fit$lambda.1se else fit$lambda.min
  coefs <- as.numeric(coef(fit, s = lam))[-1L]  # drop intercept
  m <- length(coefs)
  ev <- abs(coefs) / max(abs(coefs) + 1e-12)    # normalized magnitude in [0,1]

  if (is.null(K)) {
    sel <- which(coefs != 0)
  } else {
    sel <- order(ev, decreasing = TRUE)[seq_len(min(K, m))]
    sel <- sel[ev[sel] > 0]
  }

  list(indices = as.integer(sel),
       post.probs = ev,
       runtime = runtime)
}

# ------------------------------------------------------------------------------
# 4. CARMA wrapper
# ------------------------------------------------------------------------------

#' Run CARMA (Causal-Aware Robust Multi-Ancestry, Yang et al. 2023)
#'
#' Single-locus Bayesian fine-mapping using marginal z-scores and the LD
#' correlation matrix. Returns posterior inclusion probabilities (PIPs) over
#' all SNPs in the locus. CARMA is designed for loci of m ~ 200-2000 SNPs.
#'
#' @param X      Matrix (n x m_locus). Used to compute marginal z and LD.
#' @param Y      Numeric vector of length n.
#' @param K      Maximum number of selected SNPs (top-K by PIP). If NULL,
#'               returns the credible set at default coverage (0.95).
#' @param pip_threshold  PIP cut-off when K is NULL (default 0.5).
#' @return list(indices, post.probs, runtime).
run_carma <- function(X, Y, K = NULL, pip_threshold = 0.5, prior_p = NULL) {
  if (!requireNamespace("CARMA", quietly = TRUE))
    stop("Package CARMA required. install.packages('CARMA')")
  t0  <- Sys.time()
  n   <- nrow(X); m <- ncol(X)
  # Marginal z-scores (standardized X, centered Y)
  Y_c <- as.numeric(Y - mean(Y, na.rm = TRUE))
  Xs  <- scale(X); Xs[is.nan(Xs)] <- 0
  # Per-SNP simple regression z = beta_hat / se(beta_hat)
  beta_hat <- as.numeric(crossprod(Xs, Y_c)) / pmax(colSums(Xs^2), 1e-12)
  resid_var <- pmax((sum(Y_c^2) - beta_hat * crossprod(Xs, Y_c)) / (n - 2L), 1e-12)
  se_beta   <- sqrt(resid_var / pmax(colSums(Xs^2), 1e-12))
  z_vec     <- beta_hat / pmax(se_beta, 1e-12)
  z_vec[!is.finite(z_vec)] <- 0
  # LD = correlation matrix of standardized columns
  ld_mat    <- crossprod(Xs) / (n - 1L)
  # Optional prior-aware mode: pass the same per-SNP prior consumed by
  # SuSiE_prior (prior_weights = p_j) and varbvs (logodds = log10 p_j/(1-p_j))
  # to CARMA as a continuous logit(p_j) functional annotation through its native
  # annotation model (w.list). CARMA *learns* the annotation weight, so a
  # misspecified prior down-weights the true causals and is (correctly) harmful.
  # A constant prior (e.g. the uniform scenario) gives a degenerate annotation
  # column -> fall back to prior-free CARMA.
  carma_args <- list(z.list = list(z_vec), ld.list = list(ld_mat),
                     lambda.list = list(1), output.labels = NULL,
                     label.list = list("L"), outlier.switch = FALSE)
  if (!is.null(prior_p)) {
    lp <- qlogis(pmin(pmax(prior_p, 1e-6), 1 - 1e-6))
    if (length(unique(lp)) > 1L)
      carma_args$w.list <- list(cbind(1, lp))   # intercept + logit-prior annotation
  }
  fit <- tryCatch(do.call(CARMA::CARMA, carma_args),
    error = function(e) { warning("CARMA failed: ", conditionMessage(e)); NULL }
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) {
    return(list(indices = integer(0L),
                post.probs = rep(0, m), runtime = runtime))
  }
  pip <- as.numeric(fit[[1L]]$PIPs)
  if (length(pip) != m) pip <- c(pip, rep(0, m - length(pip)))
  if (is.null(K)) {
    sel <- which(pip > pip_threshold)
  } else {
    sel <- order(pip, decreasing = TRUE)[seq_len(min(K, m))]
  }
  list(indices = as.integer(sel),
       post.probs = pip,
       runtime = runtime)
}

# ------------------------------------------------------------------------------
# 5. FINEMAP wrapper
# ------------------------------------------------------------------------------

#' Run FINEMAP (Benner et al. 2016) via direct system2() call
#'
#' Single-locus Bayesian fine-mapping using marginal z-scores + LD matrix.
#' Calls the FINEMAP binary directly (no finemapr dependency). Uses the v1.3.1
#' file format (8-column .z, master with n_samples column).
#'
#' Binary path resolution order:
#'   1. FINEMAP_BIN env var if set;
#'   2. "~/bin/finemap/finemap_v1.3.1" if it exists (local Mac install with
#'      libzstd shipped under ~/bin/finemap/lib);
#'   3. "finemap" assumed on PATH (server Linux default).
#'
#' On macOS x86_64 the binary is dynamically linked against libzstd.1.dylib;
#' if FINEMAP_DYLD_PATH is set, it is exported to the child process as
#' DYLD_LIBRARY_PATH. Defaults to "~/bin/finemap/lib" if that directory exists.
#'
#' @param X      Matrix (n x m_locus).
#' @param Y      Numeric vector of length n.
#' @param K      Top-K by PIP cap (if NULL, returns PIP > pip_threshold).
#' @param max_causal  Maximum number of causal SNPs FINEMAP considers
#'                    (default 5; small values keep wall time bounded).
#' @param pip_threshold Cut-off when K is NULL.
#' @return list(indices, post.probs, runtime).
run_finemap <- function(X, Y, K = NULL, max_causal = 5L,
                        pip_threshold = 0.5) {
  # 1. Resolve binary path
  local_default <- path.expand("~/bin/finemap/finemap_v1.3.1")
  bin_path <- Sys.getenv("FINEMAP_BIN", unset = "")
  if (!nzchar(bin_path)) {
    bin_path <- if (file.exists(local_default)) local_default else "finemap"
  }
  # 2. Resolve DYLD path (macOS x86_64 needs libzstd)
  local_lib <- path.expand("~/bin/finemap/lib")
  dyld_path <- Sys.getenv("FINEMAP_DYLD_PATH", unset = "")
  if (!nzchar(dyld_path) && dir.exists(local_lib)) dyld_path <- local_lib
  env_vec <- if (nzchar(dyld_path))
    c(paste0("DYLD_LIBRARY_PATH=", dyld_path),
      paste0("DYLD_FALLBACK_LIBRARY_PATH=", dyld_path)) else character(0)

  t0  <- Sys.time()
  n   <- nrow(X); m <- ncol(X)
  # 3. Compute marginal z-scores and LD from (X, Y)
  Y_c <- as.numeric(Y - mean(Y, na.rm = TRUE))
  Xs  <- scale(X); Xs[is.nan(Xs)] <- 0
  beta_hat  <- as.numeric(crossprod(Xs, Y_c)) / pmax(colSums(Xs^2), 1e-12)
  resid_var <- pmax((sum(Y_c^2) - beta_hat * crossprod(Xs, Y_c)) / (n - 2L), 1e-12)
  se_beta   <- sqrt(resid_var / pmax(colSums(Xs^2), 1e-12))
  z_vec     <- beta_hat / pmax(se_beta, 1e-12)
  z_vec[!is.finite(z_vec)] <- 0
  ld_mat    <- crossprod(Xs) / (n - 1L)

  # 4. Write FINEMAP v1.3.1 input files
  td <- tempfile("fm_"); dir.create(td, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  zf <- file.path(td, "reg.z");      lf <- file.path(td, "reg.ld")
  sf <- file.path(td, "reg.snp");    cf <- file.path(td, "reg.config")
  rf <- file.path(td, "reg.cred");   mf <- file.path(td, "master")
  rsids <- paste0("v", seq_len(m))
  # 8-column header expected by v1.3.1+: rsid chromosome position allele1
  # allele2 maf beta se. Only beta and se feed the SSS algorithm.
  zdf <- data.frame(rsid = rsids, chromosome = "01", position = seq_len(m),
                    allele1 = "A", allele2 = "G", maf = 0.5,
                    beta = beta_hat, se = pmax(se_beta, 1e-8),
                    stringsAsFactors = FALSE)
  utils::write.table(zdf, zf, col.names = TRUE, row.names = FALSE,
                     quote = FALSE, sep = " ")
  utils::write.table(format(ld_mat, nsmall = 4, trim = TRUE), lf,
                     col.names = FALSE, row.names = FALSE, quote = FALSE)
  writeLines(c("z;ld;snp;config;cred;n_samples",
               paste(zf, lf, sf, cf, rf, round(n), sep = ";")), mf)

  # 5. Run FINEMAP
  ok <- tryCatch(
    system2(bin_path,
            args = c("--sss", "--in-files", mf,
                     "--n-causal-snps", max_causal),
            env = env_vec, stdout = FALSE, stderr = FALSE) == 0L,
    error = function(e) FALSE
  )
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!ok || !file.exists(sf)) {
    return(list(indices = integer(0L), post.probs = rep(0, m),
                runtime = runtime))
  }

  # 6. Parse .snp output -- column 'prob' is the per-SNP marginal PIP in v1.3.1+
  snp_out <- tryCatch(
    utils::read.table(sf, header = TRUE, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(snp_out) || !"prob" %in% names(snp_out) ||
      !"rsid" %in% names(snp_out)) {
    return(list(indices = integer(0L), post.probs = rep(0, m),
                runtime = runtime))
  }
  pip <- stats::setNames(rep(0, m), rsids)
  match_idx <- match(snp_out$rsid, rsids)
  good <- !is.na(match_idx)
  pip[match_idx[good]] <- snp_out$prob[good]

  if (is.null(K)) {
    sel <- which(pip > pip_threshold)
  } else {
    sel <- order(pip, decreasing = TRUE)[seq_len(min(K, m))]
  }
  list(indices = as.integer(sel),
       post.probs = as.numeric(pip),
       runtime = runtime)
}

# ------------------------------------------------------------------------------
# 6. Helper — evaluation metrics
# ------------------------------------------------------------------------------

#' Compute precision / recall / F1 / FDR for a selection vs ground truth
#'
#' @param selected     integer vector of selected SNP indices.
#' @param true_causal  integer vector of indices of true causal SNPs.
#' @param m            total number of SNPs.
#'
#' @return named list with precision, recall, F1, FDR, M (TP), K (set size).
evaluate_selection <- function(selected, true_causal, m) {
  K  <- length(selected)
  TP <- length(intersect(selected, true_causal))
  FP <- K - TP
  FN <- length(setdiff(true_causal, selected))

  precision <- if (K == 0) NA_real_ else TP / K
  recall    <- if (length(true_causal) == 0) NA_real_ else TP / length(true_causal)
  # F1 = 2TP / (K + |C|): the harmonic mean of precision and recall, set to 0
  # when no causal variant is selected (including an empty selection), so that
  # averages over replications do not drop the failures.
  F1        <- if (length(true_causal) == 0) NA_real_ else 2 * TP / (K + length(true_causal))
  FDR <- if (K == 0) NA_real_ else FP / K

  list(M = TP, K = K, FP = FP, FN = FN,
       precision = precision, recall = recall, F1 = F1, FDR = FDR)
}

# ------------------------------------------------------------------------------
# 5. Quick smoke test (disabled in the public release; this file is always
#    sourced as a function library, never executed standalone).
# ------------------------------------------------------------------------------
if (FALSE) {
  source("Generate_XY_realistic.R")
  set.seed(2026)
  n <- 200; m <- 300
  X_t <- Generate.X(n = n, r = 0.2, m = m, f = 0.2)
  beta_t <- c(0, rep(0.6, 5), rep(0, m - 5))
  Y_t <- Generate.Y(n = n, X = X_t, beta = beta_t, sigma = 1)

  cat("\n--- Smoke test: SuSiE ---\n")
  s <- run_susie(X_t, Y_t, K = 5)
  print(s$indices); cat("runtime:", s$runtime, "s\n")

  cat("\n--- Smoke test: varbvs ---\n")
  v <- run_varbvs(X_t, Y_t, K = 5)
  print(v$indices); cat("runtime:", v$runtime, "s\n")

  cat("\n--- Smoke test: LASSO ---\n")
  l <- run_lasso(X_t, Y_t, K = 5)
  print(l$indices); cat("runtime:", l$runtime, "s\n")

  cat("\nMetrics vs truth (1:5):\n")
  print(evaluate_selection(s$indices, 1:5, m))
}

# ------------------------------------------------------------------------------
# 7. BayesC (BGLR) wrapper -- Bayesian spike-and-slab with Normal slab
# ------------------------------------------------------------------------------

#' Run BayesC via BGLR (Habier et al. 2011; Pérez & de los Campos 2014)
#'
#' Fits the spike-and-slab model
#'   $\beta_j \sim (1 - \pi)\,\delta_0 + \pi\,\mathcal{N}(0, \sigma_\beta^2)$
#' via Gibbs sampling. This is the direct MCMC sibling of the in-house
#' framework (same conjugate spike-and-slab prior used by NaM, ExM) and a
#' strict superset of Bayesian Lasso for the purpose of fine-mapping, in
#' that it produces \emph{true} posterior inclusion probabilities
#' (PIPs $= P(\beta_j \neq 0 \mid Y, X)$) rather than a continuous shrinkage
#' score. BGLR exposes the PIP via \code{fit$ETA[[1]]$d} (proportion of MCMC
#' iterations in which $\beta_j$ was sampled from the slab).
#'
#' BayesC is also the closest Bayesian counterpart of varbvs (variational
#' spike-and-slab) and a methodologically more standard choice than
#' Bayesian Lasso for genome-wide / locus fine-mapping, because the spike
#' component directly encodes the "is this SNP causal?" question.
#'
#' Run time is dominated by the MCMC sampler. The default
#' \code{nIter = 1500, burnIn = 500} matches BGLR's package defaults and is
#' adequate for $m \le 3{,}000$.
#'
#' @param X      Matrix (n x m).
#' @param Y      Numeric vector of length n.
#' @param K      Top-K by PIP. If NULL, threshold on \code{pip_threshold}.
#' @param pip_threshold  PIP cut-off when K is NULL (default 0.5).
#' @param nIter,burnIn,thin  MCMC controls (BGLR defaults: 1500, 500, 5).
#' @param verbose  Whether to forward BGLR's chain-progress output.
#'
#' @return list(indices, post.probs, runtime). \code{post.probs} are the
#'   BayesC PIPs (proper posterior inclusion probabilities).
run_bayesc <- function(X, Y, K = NULL,
                       pip_threshold = 0.5,
                       nIter = 1500L, burnIn = 500L, thin = 5L,
                       verbose = FALSE) {
  if (!requireNamespace("BGLR", quietly = TRUE))
    stop("Package BGLR required. install.packages('BGLR')")
  t0 <- Sys.time()
  n  <- nrow(X); m <- ncol(X)
  Y_c <- as.numeric(Y - mean(Y, na.rm = TRUE))

  # BGLR writes intermediate files (mu.dat, varE.dat, ETA_1_lambda.dat, ...).
  # Isolate them in a per-worker tempdir to avoid concurrent collisions.
  td <- tempfile("bayesc_"); dir.create(td, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE, force = TRUE), add = TRUE)
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(normalizePath(td))

  fit <- tryCatch({
    ETA <- list(list(X = as.matrix(X), model = "BayesC"))
    if (verbose) {
      BGLR::BGLR(y = Y_c, ETA = ETA, nIter = nIter, burnIn = burnIn,
                 thin = thin, saveAt = "", verbose = TRUE)
    } else {
      utils::capture.output(
        out <- BGLR::BGLR(y = Y_c, ETA = ETA, nIter = nIter, burnIn = burnIn,
                          thin = thin, saveAt = "", verbose = FALSE),
        type = "output", file = NULL
      )
      out
    }
  }, error = function(e) {
    warning("BGLR::BayesC failed: ", conditionMessage(e))
    NULL
  })
  runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (is.null(fit) || is.null(fit$ETA) || is.null(fit$ETA[[1]]$d))
    return(list(indices = integer(0L), post.probs = rep(0, m), runtime = runtime))

  pip <- as.numeric(fit$ETA[[1]]$d)            # posterior inclusion probability

  if (is.null(K)) {
    sel <- which(pip > pip_threshold)
  } else {
    sel <- order(pip, decreasing = TRUE)[seq_len(min(K, m))]
  }

  list(indices    = as.integer(sel),
       post.probs = pip,
       runtime    = runtime)
}

