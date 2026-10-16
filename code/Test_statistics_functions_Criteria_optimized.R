# ==============================================================================
# File: Test_statistics_functions_Criteria_optimized.R
# Author: N. K. Doulabe
# Description: Variable selection with automatic stopping criteria (BIC, eBIC,
#              Adj-R2, Posterior Probability).
#              Section A : functions for simulation studies.
#              Section B : functions for real data analysis (+ threshold_K).
#
# Depends on: Test_statistics_functions_optimized.R
#   Inherits : compute_qr_stats(), compute_posterior(), make_quad_nodes(),
#              naive_wald_marginal(), asymptotic_wald_marginal_gh(),
#              exact_wald_marginal_gh(), asymptotic_score_marginal_gh(),
#              Fisher_MC_approximation(), Chisq_MC_approximation(),
#              ExM_posterior(), NaM_posterior()
#
# Corrections vs original:
#   - source() points to the optimized version (without integrands)
#   - mapply(integrate(...)) replaced by vectorized marginals everywhere
#   - compute_qr_stats() used in all functions (Section A and B)
#   - Unified stopping logic: length(hist) >= 1 in A and B (consistent)
#   - cbind() inside the loop (Section B) removed — X_inclus rebuilt
#   - All vector initializations: numeric(0) instead of c()
#   - clusterExport() cleaned up (no more integrands)
# ==============================================================================

source("Test_statistics_functions_optimized.R")


# ==============================================================================
# Internal utility: computation of the stopping metrics
# ==============================================================================

#' Compute RSS, adjusted R2, BIC and eBIC after selecting the best candidate
#' @param RSS_new  RSS of the model with the candidate added.
#' @param RSS_current RSS before the addition (unused here, kept for clarity).
#' @param TSS  Total sum of squares.
#' @param n    Sample size.
#' @param k    Current step (df_model = k + 1).
#' @param m    Total number of candidate variables (for the eBIC penalty).
#' @param gamma_ebic  Extended BIC parameter (Chen & Chen 2008), in [0, 1].
#'                   Default 1 (most conservative, recommended for $m \gg n$).
#' @return Named list: r2_adj, ebic, bic.
compute_criteria <- function(RSS_new, TSS, n, k, m = NA_integer_,
                              gamma_ebic = 1.0) {
  df_model <- k + 1L
  r2_adj   <- 1 - (RSS_new / max(1L, n - df_model)) / (TSS / max(1L, n - 1L))
  bic      <- n * log(RSS_new / n) + log(n) * df_model
  # Extended BIC penalty: 2 * gamma * log(choose(m, k)). Uses lchoose for
  # numerical stability. If m is not provided (NA) we fall back to standard
  # BIC by setting the eBIC penalty to zero.
  ebic_penalty <- if (is.finite(m) && k >= 1L && k <= m) {
    2 * gamma_ebic * lchoose(m, k)
  } else {
    0
  }
  ebic <- bic + ebic_penalty
  list(r2_adj = r2_adj, bic = bic, ebic = ebic)
}

#' Evaluation of the stopping criterion
#' @description Unified convention: checks as soon as the history contains
#'   at least one value (length >= 1), identical in Section A and B.
#' @return TRUE if the algorithm must stop.
check_stop <- function(criteria, r2_adj, ebic, bic, best_prob,
                       r2_hist, ebic_hist, bic_hist, target_prob) {
  if (criteria == "R2" && length(r2_hist) >= 1L) {
    prev <- r2_hist[length(r2_hist)]
    if (is.finite(prev) && is.finite(r2_adj) && r2_adj < prev) return(TRUE)
  }
  if (criteria == "eBIC" && length(ebic_hist) >= 1L) {
    prev <- ebic_hist[length(ebic_hist)]
    if (is.finite(prev) && is.finite(ebic) && ebic > prev) return(TRUE)
  }
  if (criteria == "BIC" && length(bic_hist) >= 1L) {
    prev <- bic_hist[length(bic_hist)]
    if (is.finite(prev) && is.finite(bic) && bic > prev) return(TRUE)
  }
  if (criteria == "Posterior" && is.finite(best_prob) && best_prob < target_prob) {
    return(TRUE)
  }
  FALSE
}

#' Compute the RSS after selecting a candidate (common utility)
#' @param beta_best  Estimated coefficient of the best candidate.
#' @param ss_best    Residual sum of squares of the best candidate.
#' @param RSS_current Current RSS before addition.
#' @return Updated RSS, clamped at 1e-9.
update_rss <- function(beta_best, ss_best, RSS_current) {
  RSS_new <- RSS_current - beta_best^2 * ss_best
  if (!is.finite(RSS_new) || RSS_new <= 0) RSS_new <- 1e-9
  RSS_new
}

#' Standard return for all Criteria functions
make_criteria_result <- function(indices.inclus, post.probs,
                                 r2_hist, ebic_hist, bic_hist,
                                 threshold = NULL) {
  out <- list(
    indices       = indices.inclus,
    post.probs    = post.probs,
    r_squared_adj = r2_hist,
    ebic_values   = ebic_hist,
    bic_values    = bic_hist,
    K             = length(indices.inclus)
  )
  if (!is.null(threshold))
    out$indices_under_threshold <- sum(indices.inclus < threshold)
  out
}


# ==============================================================================
# SECTION A : Functions for simulation study
# ==============================================================================

#' Naive Wald with stopping criteria (simulation)
#' @param X        Genotype matrix (n x m), raw (will be centered-scaled).
#' @param Y        Phenotype vector (n x 1).
#' @param a.priori Vector of prior probabilities (length m).
#' @param tau      Standard deviation of the prior.
#' @param threshold Threshold for counting true positives (default 16).
#' @param target_prob Posterior threshold for the "Posterior" criterion (default 0.95).
#' @param criteria  "BIC", "eBIC", "R2", or "Posterior".
Naive_Wald_Criteria <- function(X, Y, a.priori, tau,
                                threshold = 16L, target_prob = 0.95,
                                criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  k_max <- min(m, n - 2L)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(k_max)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      r_vec      <- nums / sqrt(TSS * ss_cand)
      r_vec      <- pmin(pmax(r_vec, -1 + 1e-12), 1 - 1e-12)
      Z_vector   <- r_vec * sqrt(n - 2L) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test))
      beta_cand  <- nums / ss_cand
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "wald"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- naive_wald_marginal(Z_vector, rho_vector, n, tau)
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new    <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals  <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


#' Asymptotic Wald with stopping criteria (simulation)
#' @param gh List of Gauss-Hermite nodes from make_quad_nodes(tau = tau).
Asymptotic_Wald_Criteria <- function(X, Y, a.priori, tau, gh,
                                     threshold = 16L, target_prob = 0.95,
                                     criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  k_max <- min(m, n - 2L)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(k_max)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      r_vec      <- pmin(pmax(nums / sqrt(TSS * ss_cand), -1 + 1e-12), 1 - 1e-12)
      Z_vector   <- r_vec * sqrt(n - 2L) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test))
      beta_cand  <- nums / ss_cand
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "wald"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- asymptotic_wald_marginal_gh(Z_vector, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


#' Exact Wald with stopping criteria (simulation)
#' @param gh List of Gauss-Hermite nodes from make_quad_nodes(tau = tau).
Exact_Wald_Criteria <- function(X, Y, a.priori, tau, gh,
                                threshold = 16L, target_prob = 0.95,
                                criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  k_max <- min(m, n - 2L)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(k_max)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums         <- as.vector(crossprod(Y, X_cand))
      ss_cand      <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand    <- nums / ss_cand
      RSS_prov     <- pmax(TSS - beta_cand^2 * ss_cand, 1e-9)
      sigma2       <- RSS_prov / max(1L, n - 2L)
      Z_vector     <- beta_cand / sqrt(pmax(sigma2 / ss_cand, 1e-12))
      rho_vector   <- numeric(length(indices_remain_to_test))
      Z_hat_vector <- Z_vector
      ss_res       <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats     <- compute_qr_stats(X, Y,
                                       X[, indices.inclus, drop = FALSE],
                                       X_cand, n, k, "wald")
      rho_vector   <- qr_stats$rho_vector
      Z_hat_vector <- qr_stats$Z_vector * sqrt(1 - rho_vector)

      X_inclus     <- X[, indices.inclus, drop = FALSE]
      qr_tmp       <- qr(X_inclus)
      Q_mat        <- qr.Q(qr_tmp)
      X_cand_res   <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res        <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res       <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand    <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_hat_vector[!is.finite(Z_hat_vector)] <- 0
    integrals    <- exact_wald_marginal_gh(Z_hat_vector, rho_vector, n, k, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dt(Z_hat_vector, df = n - k - 1L),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


#' Asymptotic Score with stopping criteria (simulation)
#' @param gh List of Gauss-Hermite nodes from make_quad_nodes(tau = tau).
Asymptotic_Score_Criteria <- function(X, Y, a.priori, tau, gh,
                                      threshold = 16L, target_prob = 0.95,
                                      criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- as.numeric(scale(Y, center = TRUE, scale = FALSE))
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- pmax(sum(Y^2), 1e-9)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    if (n - k - 1L <= 1L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand  <- nums / ss_cand
      Q_hat      <- nums / sqrt(ss_cand * TSS / max(1L, n - 1L))
      rho_vector <- numeric(length(indices_remain_to_test))
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "score"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Q_hat      <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Q_hat[!is.finite(Q_hat)] <- 0
    integrals    <- asymptotic_score_marginal_gh(Q_hat, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Q_hat),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


#' Naive Multiple Test with stopping criteria (simulation)
Naive_Multiple_Criteria <- function(X, Y, a.priori, tau2, S,
                                    threshold = 16L, target_prob = 0.95,
                                    criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    df1 <- k; df2 <- n - k - 1L
    if (df2 <= 0L) break

    if (k == 1L) {
      denom     <- pmax(sqrt(TSS * colSums(X_cand^2)), 1e-12)
      r_sq      <- pmin((as.vector(crossprod(Y, X_cand)) / denom)^2, 1 - 1e-12)
      r_sq[!is.finite(r_sq)] <- 0
      F_stat    <- r_sq * (n - 2L) / (1 - r_sq)
      nums      <- as.vector(crossprod(Y, X_cand))
      ss_cand   <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand <- nums / ss_cand
      ss_res    <- ss_cand
      X_inclus  <- NULL
    } else {
      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_inclus  <- qr(X_inclus)
      Y_res      <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
      RSS_prov   <- sum(Y_res^2)
      RSS_new_v  <- pmax(RSS_prov - beta_cand^2 * ss_res, 1e-9)
      SSR_new_v  <- TSS - RSS_new_v
      F_stat     <- (SSR_new_v / df1) / (RSS_new_v / df2)
    }

    Q_chi <- k * F_stat
    Q_chi[!is.finite(Q_chi) | Q_chi < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Q_chi, X_inclus, X_cand, p_vec, tau2, df2 = Inf)
    } else mapply(function(q, j, p) {
      if (q < 0) return(NA_real_)
      X_mod <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
      NaM_posterior(Q_chi = q, X_tilde_jk = X_mod, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
    }, Q_chi, indices_remain_to_test, p_vec)

    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


#' Exact Multiple Test with stopping criteria (simulation)
Exact_Multiple_Criteria <- function(X, Y, a.priori, tau2, S,
                                    threshold = 16L, target_prob = 0.95,
                                    criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    df1 <- k; df2 <- n - k - 1L
    if (df2 <= 0L) break

    if (k == 1L) {
      denom     <- pmax(sqrt(TSS * colSums(X_cand^2)), 1e-12)
      r_sq      <- pmin((as.vector(crossprod(Y, X_cand)) / denom)^2, 1 - 1e-12)
      r_sq[!is.finite(r_sq)] <- 0
      Z_vector  <- r_sq * (n - 2L) / (1 - r_sq)
      nums      <- as.vector(crossprod(Y, X_cand))
      ss_cand   <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand <- nums / ss_cand
      ss_res    <- ss_cand
      X_inclus  <- NULL
    } else {
      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_inclus  <- qr(X_inclus)
      Y_res      <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
      RSS_prov   <- sum(Y_res^2)
      RSS_new_v  <- pmax(RSS_prov - beta_cand^2 * ss_res, 1e-9)
      SSR_new_v  <- TSS - RSS_new_v
      Z_vector   <- (SSR_new_v / df1) / (RSS_new_v / df2)
    }

    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Z_vector, X_inclus, X_cand, p_vec, tau2, df2 = df2)
    } else mapply(function(z, j, p) {
      if (z < 0) return(NA_real_)
      X_mod <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
      ExM_posterior(Q_W = z, X_tilde_jk = X_mod, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
    }, Z_vector, indices_remain_to_test, p_vec)

    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold) && k >= threshold) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist, threshold)
}


# ==============================================================================
# Simulation Wrapper (Section A)
# ==============================================================================

#' Simulation with stopping criteria
#' @description Runs parallel simulations with the 6 methods under a given
#'   stopping criterion. Returns a summary table.
#' @param B           Number of simulation iterations.
#' @param S           Number of MC samples for the multiple tests.
#' @param target_prob Posterior threshold for the "Posterior" criterion.
#' @param threshold   True positive threshold (default 16).
#' @param criteria    One criterion among "BIC", "eBIC", "R2", "Posterior".
Criteria_Simulation_Table <- function(B, X, sigma, a.priori, p_values, tau_values,
                                      S = 1000L, target_prob = 0.95,
                                      threshold = 16L, criteria = "BIC") {
  m <- ncol(X); n <- nrow(X)

  num_cores <- detectCores() - 1L
  cl        <- makeCluster(num_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  clusterEvalQ(cl, { library(MASS); library(mvtnorm); library(statmod) })

  clusterExport(cl, c(
    "Generate.Y", "Power_Fun",
    "make_quad_nodes", "compute_qr_stats", "compute_posterior",
    "compute_criteria", "check_stop", "update_rss", "make_criteria_result",
    "naive_wald_marginal", "asymptotic_wald_marginal_gh",
    "exact_wald_marginal_gh", "asymptotic_score_marginal_gh",
    "Naive_Wald_Criteria", "Asymptotic_Wald_Criteria",
    "Exact_Wald_Criteria", "Asymptotic_Score_Criteria",
    "Naive_Multiple_Criteria", "Exact_Multiple_Criteria",
    "Fisher_MC_approximation", "Chisq_MC_approximation",
    "ExM_posterior", "NaM_posterior",
    "multiple_theta", "nb_sum_logweights", "multiple_exact_logdens",
    "multiple_posterior_exact", "multiple_posterior_exact_vec", ".multiple_exact_memo"
  ))

  final_summary <- data.frame()

  for (p in p_values) {
    for (tau in tau_values) {
      b_root <- try(uniroot(Power_Fun, lower = 0, upper = 5,
                            n = n, alpha = 0.05, p = p), silent = TRUE)
      b.val  <- if (inherits(b_root, "try-error")) 0.5 else b_root$root
      beta   <- c(1, rep(b.val, 15L), rep(0, m - 15L))

      cat(sprintf("Running: criteria=%s p=%.2f tau=%.2f\n", criteria, p, tau))

      raw_results <- foreach(b = seq_len(B), .combine = rbind,
                             .packages = c("stats", "statmod")) %dopar% {

        gh <- make_quad_nodes(30L, tau)
        Y  <- Generate.Y(n, X, beta, sigma)

        safe_get <- function(vec) if (length(vec) > 0L) tail(vec, 1L) else NA_real_

        runs <- list(
          NaW = function() Naive_Wald_Criteria(X, Y, a.priori, tau, threshold, target_prob, criteria),
          AsW = function() Asymptotic_Wald_Criteria(X, Y, a.priori, tau, gh, threshold, target_prob, criteria),
          ExW = function() Exact_Wald_Criteria(X, Y, a.priori, tau, gh, threshold, target_prob, criteria),
          AsS = function() Asymptotic_Score_Criteria(X, Y, a.priori, tau, gh, threshold, target_prob, criteria),
          NaM = function() Naive_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, criteria),
          ExM = function() Exact_Multiple_Criteria(X, Y, a.priori, tau^2, S, threshold, target_prob, criteria)
        )

        iter_mat <- matrix(NA_real_, nrow = 6L, ncol = 7L)
        colnames(iter_mat) <- c("Method", "K", "Posterior", "R2_adj", "eBIC", "BIC", "indices_under_threshold")

        for (i in seq_along(runs)) {
          res <- tryCatch(runs[[i]](), error = function(e) NULL)
          iter_mat[i, 1L] <- names(runs)[i]
          if (!is.null(res)) {
            iter_mat[i, 2L] <- length(res$indices)
            iter_mat[i, 3L] <- safe_get(res$post.probs)
            iter_mat[i, 4L] <- safe_get(res$r_squared_adj)
            iter_mat[i, 5L] <- safe_get(res$ebic_values)
            iter_mat[i, 6L] <- safe_get(res$bic_values)
            iter_mat[i, 7L] <- res$indices_under_threshold
          }
        }
        iter_mat
      }

      raw_df <- as.data.frame(raw_results, stringsAsFactors = FALSE)
      for (j in 2L:7L) raw_df[, j] <- suppressWarnings(as.numeric(raw_df[, j]))

      for (mn in c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM")) {
        sub  <- raw_df[raw_df[, 1L] == mn, ]
        K_v  <- sub[, 2L]
        avgK <- mean(K_v, na.rm = TRUE)
        cvK  <- if (is.na(avgK) || avgK == 0) 0 else sd(K_v, na.rm = TRUE) / avgK

        final_summary <- rbind(final_summary, data.frame(
          Method                = mn,
          criteria              = criteria,
          p                     = p,
          tau                   = tau,
          Average_K             = avgK,
          Median_K              = median(K_v, na.rm = TRUE),
          CV_K                  = cvK,
          indices_under_threshold = mean(sub[, 7L], na.rm = TRUE),
          Average_prob          = mean(sub[, 3L], na.rm = TRUE),
          Average_R2_adj        = mean(sub[, 4L], na.rm = TRUE),
          Average_eBIC          = mean(sub[, 5L], na.rm = TRUE),
          Average_BIC           = mean(sub[, 6L], na.rm = TRUE)
        ))
      }
    }
  }

  final_summary
}


# ==============================================================================
# SECTION B : Functions for real data (+ threshold_K)
# ==============================================================================
# Identical to Section A but with:
#   - threshold_K : forced stop if too many variables are selected
#   - No return of indices_under_threshold (not relevant outside simulation)

#' Naive Wald for real data
Naive_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau,
                                         threshold_K = NULL, target_prob = 0.95,
                                         criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      r_vec      <- pmin(pmax(nums / sqrt(TSS * ss_cand), -1 + 1e-12), 1 - 1e-12)
      Z_vector   <- r_vec * sqrt(n - 2L) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test))
      beta_cand  <- nums / ss_cand
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "wald"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- naive_wald_marginal(Z_vector, rho_vector, n, tau)
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}


#' Asymptotic Wald for real data
Asymptotic_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, gh,
                                               threshold_K = NULL, target_prob = 0.95,
                                               criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      r_vec      <- pmin(pmax(nums / sqrt(TSS * ss_cand), -1 + 1e-12), 1 - 1e-12)
      Z_vector   <- r_vec * sqrt(n - 2L) / sqrt(1 - r_vec^2)
      rho_vector <- numeric(length(indices_remain_to_test))
      beta_cand  <- nums / ss_cand
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "wald"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- asymptotic_wald_marginal_gh(Z_vector, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}


#' Exact Wald for real data
Exact_Wald_Criteria_For_Data <- function(X, Y, a.priori, tau, gh,
                                          threshold_K = NULL, target_prob = 0.95,
                                          criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums         <- as.vector(crossprod(Y, X_cand))
      ss_cand      <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand    <- nums / ss_cand
      RSS_prov     <- pmax(TSS - beta_cand^2 * ss_cand, 1e-9)
      sigma2       <- RSS_prov / max(1L, n - 2L)
      Z_vector     <- beta_cand / sqrt(pmax(sigma2 / ss_cand, 1e-12))
      rho_vector   <- numeric(length(indices_remain_to_test))
      Z_hat_vector <- Z_vector
      ss_res       <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats     <- compute_qr_stats(X, Y,
                                       X[, indices.inclus, drop = FALSE],
                                       X_cand, n, k, "wald")
      rho_vector   <- qr_stats$rho_vector
      Z_hat_vector <- qr_stats$Z_vector * sqrt(1 - rho_vector)

      X_inclus     <- X[, indices.inclus, drop = FALSE]
      qr_tmp       <- qr(X_inclus)
      Q_mat        <- qr.Q(qr_tmp)
      X_cand_res   <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res        <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res       <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand    <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Z_hat_vector[!is.finite(Z_hat_vector)] <- 0
    integrals    <- exact_wald_marginal_gh(Z_hat_vector, rho_vector, n, k, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dt(Z_hat_vector, df = n - k - 1L),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}


#' Asymptotic Score for real data
Asymptotic_Score_Criteria_For_Data <- function(X, Y, a.priori, tau, gh,
                                                threshold_K = NULL, target_prob = 0.95,
                                                criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- as.numeric(scale(Y, center = TRUE, scale = FALSE))
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- pmax(sum(Y^2), 1e-9)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    if (n - k - 1L <= 1L) break

    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      nums       <- as.vector(crossprod(Y, X_cand))
      ss_cand    <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand  <- nums / ss_cand
      Q_hat      <- nums / sqrt(ss_cand * TSS / max(1L, n - 1L))
      rho_vector <- numeric(length(indices_remain_to_test))
      ss_res     <- ss_cand
    } else {
      if (length(indices.inclus) >= n - 1L) break
      qr_stats_try <- tryCatch(compute_qr_stats(X, Y,
                                     X[, indices.inclus, drop = FALSE],
                                     X_cand, n, k, "score"),
                                 error = function(e) NULL)
      if (is.null(qr_stats_try)) break  # singular matrix safety
      qr_stats   <- qr_stats_try
      Q_hat      <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector

      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_tmp_try <- tryCatch(qr(X_inclus), error = function(e) NULL)
      if (is.null(qr_tmp_try)) break  # singular matrix safety
      qr_tmp     <- qr_tmp_try
      Q_mat      <- qr.Q(qr_tmp)
      X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)
      Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
    }

    Q_hat[!is.finite(Q_hat)] <- 0
    integrals    <- asymptotic_score_marginal_gh(Q_hat, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Q_hat),
                                      a.priori[indices_remain_to_test])
    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}


#' Naive Multiple for real data
Naive_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S,
                                              threshold_K = NULL, target_prob = 0.95,
                                              criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    df1 <- k; df2 <- n - k - 1L
    if (df2 <= 0L) break

    if (k == 1L) {
      denom     <- pmax(sqrt(TSS * colSums(X_cand^2)), 1e-12)
      r_sq      <- pmin((as.vector(crossprod(Y, X_cand)) / denom)^2, 1 - 1e-12)
      r_sq[!is.finite(r_sq)] <- 0
      F_stat    <- r_sq * (n - 2L) / (1 - r_sq)
      nums      <- as.vector(crossprod(Y, X_cand))
      ss_cand   <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand <- nums / ss_cand
      ss_res    <- ss_cand
      X_inclus  <- NULL
    } else {
      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_inclus  <- qr(X_inclus)
      Y_res      <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
      RSS_prov   <- sum(Y_res^2)
      RSS_new_v  <- pmax(RSS_prov - beta_cand^2 * ss_res, 1e-9)
      SSR_new_v  <- TSS - RSS_new_v
      F_stat     <- (SSR_new_v / df1) / (RSS_new_v / df2)
    }

    Q_chi <- k * F_stat
    Q_chi[!is.finite(Q_chi) | Q_chi < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Q_chi, X_inclus, X_cand, p_vec, tau2, df2 = Inf)
    } else mapply(function(q, j, p) {
      if (q < 0) return(NA_real_)
      X_mod <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
      NaM_posterior(Q_chi = q, X_tilde_jk = X_mod, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
    }, Q_chi, indices_remain_to_test, p_vec)

    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}


#' Exact Multiple for real data
Exact_Multiple_Criteria_For_Data <- function(X, Y, a.priori, tau2, S,
                                              threshold_K = NULL, target_prob = 0.95,
                                              criteria = "BIC") {
  if (!criteria %in% c("R2", "eBIC", "BIC", "Posterior", "ALL"))
    stop("Criteria invalid. Must be one of: 'R2', 'eBIC', 'BIC', 'Posterior', 'ALL'.")

  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  r2_hist <- ebic_hist <- bic_hist <- post.probs <- numeric(0)
  TSS         <- sum(Y^2)
  RSS_current <- TSS

  for (k in seq_len(m)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    df1 <- k; df2 <- n - k - 1L
    if (df2 <= 0L) break

    if (k == 1L) {
      denom     <- pmax(sqrt(TSS * colSums(X_cand^2)), 1e-12)
      r_sq      <- pmin((as.vector(crossprod(Y, X_cand)) / denom)^2, 1 - 1e-12)
      r_sq[!is.finite(r_sq)] <- 0
      Z_vector  <- r_sq * (n - 2L) / (1 - r_sq)
      nums      <- as.vector(crossprod(Y, X_cand))
      ss_cand   <- pmax(colSums(X_cand^2), 1e-12)
      beta_cand <- nums / ss_cand
      ss_res    <- ss_cand
      X_inclus  <- NULL
    } else {
      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_inclus  <- qr(X_inclus)
      Y_res      <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_res     <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(crossprod(Y_res, X_cand_res)) / ss_res
      RSS_prov   <- sum(Y_res^2)
      RSS_new_v  <- pmax(RSS_prov - beta_cand^2 * ss_res, 1e-9)
      SSR_new_v  <- TSS - RSS_new_v
      Z_vector   <- (SSR_new_v / df1) / (RSS_new_v / df2)
    }

    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vec <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Z_vector, X_inclus, X_cand, p_vec, tau2, df2 = df2)
    } else mapply(function(z, j, p) {
      if (z < 0) return(NA_real_)
      X_mod <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
      ExM_posterior(Q_W = z, X_tilde_jk = X_mod, X_k = X_inclus, n = n, p = p, tau2 = tau2, S = S)
    }, Z_vector, indices_remain_to_test, p_vec)

    a.posteriori[!is.finite(a.posteriori)] <- -Inf
    if (all(!is.finite(a.posteriori))) break

    best_local <- which.best(a.posteriori)
    best_prob  <- a.posteriori[best_local]
    best_var   <- indices_remain_to_test[best_local]

    RSS_new   <- update_rss(beta_cand[best_local], ss_res[best_local], RSS_current)
    crit_vals <- compute_criteria(RSS_new, TSS, n, k, m = m)

    if (!is.null(threshold_K) && k >= threshold_K) break  # safety cap
    if (criteria != "ALL" &&
        check_stop(criteria, crit_vals$r2_adj, crit_vals$ebic, crit_vals$bic,
                   best_prob, r2_hist, ebic_hist, bic_hist, target_prob)) break

    indices.inclus         <- c(indices.inclus, best_var)
    indices_remain_to_test <- setdiff(indices_remain_to_test, best_var)
    r2_hist    <- c(r2_hist,    crit_vals$r2_adj)
    ebic_hist   <- c(ebic_hist,   crit_vals$ebic)
    bic_hist   <- c(bic_hist,   crit_vals$bic)
    post.probs <- c(post.probs, best_prob)
    RSS_current <- RSS_new

    if (!is.null(threshold_K) && length(indices.inclus) >= threshold_K) break
  }

  make_criteria_result(indices.inclus, post.probs, r2_hist, ebic_hist, bic_hist)
}
