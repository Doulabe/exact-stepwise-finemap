# ==============================================================================
# File: Test_statistics_functions_optimized.R
# Author: N. K. Doulabe
# Description: Implements the finite-sample exact, asymptotic and naive Wald, Score
#              and Multiple-SNP test-statistic distributions and their stepwise
#              posterior scores for fine-mapping.
#
# Optimizations vs original:
#   - compute_qr_stats() : shared QR block extracted (removes ~60 duplicated lines)
#   - Naive_Wald         : exact closed-form solution (replaces mapply + integrate)
#   - Asymptotic_Wald    : vectorized Gauss-Hermite quadrature
#   - Exact_Wald         : GH quadrature with non-central dt (loop over nodes)
#   - Asymptotic_Score   : vectorized GH quadrature (mean.Q / var.Q as matrices)
#   - make_quad_nodes()  : GH nodes pre-computed only once outside the simulation
# ==============================================================================

library(MASS)
library(stats)
library(Matrix)
library(doParallel)
library(mvtnorm)
library(dplyr)
library(statmod)   # gauss.quad.prob()


# ==============================================================================
# Section 1 : Common utilities
# ==============================================================================

#' Pre-computation of the Gauss-Hermite nodes and weights
#' @description To be called ONLY ONCE before the simulation. The result is
#'   shared between Asymptotic_Wald, Exact_Wald and Asymptotic_Score.
#' @param n_nodes Number of nodes (30 recommended — good accuracy/speed trade-off).
#' @param tau     Standard deviation of the prior (must match the one used in the functions).
#' @return List with $nodes (vector) and $weights (vector).
make_quad_nodes <- function(n_nodes = 30L, tau) {
  gh <- statmod::gauss.quad.prob(n_nodes, dist = "normal", mu = 0, sigma = tau)
  list(nodes = gh$nodes, weights = gh$weights)
}

#' Shared QR block — projection and computation of rho
#' @description Extracts the duplicated code present in the 4 main functions.
#'   Computes the residuals, beta, Z-score and the rho correction factor for
#'   all candidate SNPs simultaneously (vectorized).
#' @param X        Full genotype matrix (n x m), centered and scaled.
#' @param Y        Phenotype vector (n x 1), centered.
#' @param X_inclus Sub-matrix of the already included variables (n x k).
#' @param X_cand   Sub-matrix of the candidate variables (n x m_cand).
#' @param n        Sample size.
#' @param k        Current step (determines df = n - k - 1).
#' @param type     "wald" : Z = beta/se  |  "score" : Q = <Y, X_res> / sqrt(D * sigma0^2),
#'                 with sigma0^2 the residual variance of the null fit.
#' @return List with $Z_vector (m_cand,) and $rho_vector (m_cand,).
compute_qr_stats <- function(X, Y, X_inclus, X_cand, n, k,
                             type = c("wald", "score")) {
  type      <- match.arg(type)
  qr_inclus <- qr(X_inclus)
  Q_mat     <- qr.Q(qr_inclus)

  Y_res      <- Y - Q_mat %*% crossprod(Q_mat, Y)
  X_cand_res <- X_cand - Q_mat %*% crossprod(Q_mat, X_cand)

  ss_cand   <- pmax(colSums(X_cand_res^2), 1e-9)

  # Numerator : Y_res for Wald, Y for Score (efficient statistic)
  num_vec   <- as.vector(crossprod(
    if (type == "score") Y else Y_res,
    X_cand_res
  ))

  beta_cand   <- num_vec / ss_cand
  RSS_current <- sum(Y_res^2)
  RSS_new     <- RSS_current - beta_cand^2 * ss_cand
  sigma2_new  <- RSS_new / (n - k - 1L)

  Z_vector <- if (type == "wald") {
    beta_cand / sqrt(sigma2_new / ss_cand)
  } else {
    # Score: residual variance of the null fit (intercept + k-1 selected SNPs).
    sigma2_null <- RSS_current / (n - k)
    num_vec / sqrt(ss_cand * pmax(sigma2_null, 1e-12))
  }

  # rho correction factor (identical for Wald and Score)
  XtX_inv    <- chol2inv(qr.R(qr_inclus))
  R_matrix   <- crossprod(X_inclus, X_cand) / n
  M          <- (n * XtX_inv) %*% R_matrix
  rho_vector <- pmax(0, pmin(colSums(R_matrix * M), 1 - 1e-9))

  list(Z_vector = Z_vector, rho_vector = rho_vector)
}

#' Posterior update (common to all single-SNP methods)
#' @description Factored out here to avoid repetition in each function.
compute_posterior <- function(integrals, densities, p_vector) {
  numerator    <- p_vector * integrals
  denominator  <- numerator + (1 - p_vector) * densities
  q <- numerator / (denominator + 1e-20)
  # Log posterior odds, carried as an attribute for the selection step: the
  # posterior itself rounds to exactly 1 when (1-p) B < eps * p A (|Z| above
  # roughly 9 on strong real-data signals), and which.max would then break the
  # ties by column order, i.e. by genomic position. Selection uses which.best.
  lo <- log(p_vector) - log1p(-p_vector) + log(integrals) - log(densities)
  structure(q, logodds = lo)
}

#' Selection by log posterior odds (saturation-safe argmax)
#' @description Strictly monotone in the posterior, so the selected candidate
#'   is unchanged wherever the posterior does not saturate; where it does, the
#'   log odds still discriminate. Falls back to the posterior when the log odds
#'   are unavailable (e.g. Monte-Carlo scores) or undefined.
which.best <- function(post) {
  lo <- attr(post, "logodds")
  if (is.null(lo)) return(which.max(post))
  bad <- is.na(lo)
  if (any(bad)) lo[bad] <- ifelse(is.na(post[bad]), NA_real_,
                                  qlogis(pmin(pmax(post[bad], 1e-300), 1 - 1e-16)))
  which.max(lo)
}


# ==============================================================================
# Section 2 : Vectorized marginal likelihoods
# ==============================================================================

#' Exact marginal likelihood — Naive Wald
#' @description Closed-form analytical solution :
#'   ∫ dnorm(Z|λ√(n(1-ρ)),1) × dnorm(λ|0,τ) dλ = dnorm(Z|0, √(1 + n(1-ρ)τ²))
#' @param Z_vector   Z vector (m_cand,).
#' @param rho_vector rho vector (m_cand,).
naive_wald_marginal <- function(Z_vector, rho_vector, n, tau) {
  sd_marginal <- sqrt(1 + n * (1 - rho_vector) * tau^2)
  dnorm(Z_vector, mean = 0, sd = sd_marginal)
}

#' Marginal likelihood — Asymptotic Wald (Gauss-Hermite)
#' @description ∫ dnorm(Z|mean.Z, √var.Z) × dnorm(λ|0,τ) dλ
#'   Approximated by GH quadrature : O(n_nodes × m_cand) matrix operations.
#' @param gh List produced by make_quad_nodes().
asymptotic_wald_marginal_gh <- function(Z_vector, rho_vector, n, gh) {
  nodes   <- gh$nodes
  weights <- gh$weights
  m_cand  <- length(Z_vector)

  mean_mat <- outer(nodes, sqrt(n * (1 - rho_vector)))
  var_mat  <- 1 + 0.5 * outer(nodes^2, (1 - rho_vector))
  Z_mat    <- matrix(Z_vector, nrow = length(nodes), ncol = m_cand, byrow = TRUE)

  lik_mat  <- dnorm(Z_mat, mean = mean_mat, sd = sqrt(var_mat))
  lik_mat[!is.finite(lik_mat)] <- 0

  as.vector(weights %*% lik_mat)
}

#' Marginal likelihood — Exact Wald (Gauss-Hermite with non-central t)
#' @description ∫ dt(Z|df, ncp=λ√(n(1-ρ))) × dnorm(λ|0,τ) dλ
#'   Note : dt() does not accept a vectorized ncp — loop over n_nodes (not over m_cand).
#'   Gain : m_cand / n_nodes calls vs m_cand in the original.
#' @param k Current step (df = n - k - 1).
exact_wald_marginal_gh <- function(Z_vector, rho_vector, n, k, gh) {
  nodes   <- gh$nodes
  weights <- gh$weights
  m_cand  <- length(Z_vector)
  df      <- n - k - 1L

  ncp_mat <- outer(nodes, sqrt(n * (1 - rho_vector)))  # (n_nodes × m_cand)
  lik_mat <- matrix(0, nrow = length(nodes), ncol = m_cand)

  for (i in seq_along(nodes)) {
    vals <- dt(Z_vector, df = df, ncp = ncp_mat[i, ])
    vals[!is.finite(vals)] <- 0
    lik_mat[i, ] <- vals
  }

  as.vector(weights %*% lik_mat)
}

#' Marginal likelihood — Asymptotic Score (Gauss-Hermite)
#' @description ∫ dnorm(Q|mean.Q, √var.Q) × dnorm(λ|0,τ) dλ
#'   mean.Q and var.Q built via outer() — fully vectorized.
asymptotic_score_marginal_gh <- function(Q_vector, rho_vector, n, gh) {
  nodes   <- gh$nodes
  weights <- gh$weights
  m_cand  <- length(Q_vector)

  rho1_vec  <- 1 - rho_vector
  lam2_mat  <- outer(nodes^2, rho1_vec)
  denom_mat <- 1 + lam2_mat

  mean_mat  <- outer(nodes, sqrt(n * rho1_vec)) / sqrt(denom_mat)
  var_mat   <- pmax((2 + lam2_mat) / (2 * denom_mat^3), 1e-12)

  Q_mat     <- matrix(Q_vector, nrow = length(nodes), ncol = m_cand, byrow = TRUE)
  lik_mat   <- dnorm(Q_mat, mean = mean_mat, sd = sqrt(var_mat))
  lik_mat[!is.finite(lik_mat)] <- 0

  as.vector(weights %*% lik_mat)
}


# ==============================================================================
# Section 3 : Main single-SNP functions
# ==============================================================================

#' Naive Wald Test (NaW) — optimized
#' @description Step-by-step selection by naive Wald with exact marginal likelihood.
#' @param X        Genotype matrix (n x m).
#' @param Y        Phenotype vector (n x 1).
#' @param a.priori Vector of prior probabilities (length m).
#' @param K        Maximum number of variables to select.
#' @param tau      Standard deviation of the prior on the effect.
Naive_Wald <- function(X, Y, a.priori, K, tau) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- numeric(K)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      r          <- cor(Y, X_cand)[1L, ]; r[is.na(r)] <- 0
      Z_vector   <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2L) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      qr_stats   <- compute_qr_stats(X, Y, X[, indices.inclus, drop = FALSE], X_cand, n, k, "wald")
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- naive_wald_marginal(Z_vector, rho_vector, n, tau)
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector), a.priori[indices_remain_to_test])

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}

#' Asymptotic Wald Test (AsW) — optimized
#' @description Selection by asymptotic Wald with variance correction.
#' @param gh List produced by make_quad_nodes(tau = tau).
Asymptotic_Wald <- function(X, Y, a.priori, K, tau, gh) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- numeric(K)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      r          <- cor(Y, X_cand)[1L, ]; r[is.na(r)] <- 0
      Z_vector   <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2L) / sqrt(1 - r^2))
      rho_vector <- numeric(length(indices_remain_to_test))
    } else {
      qr_stats   <- compute_qr_stats(X, Y, X[, indices.inclus, drop = FALSE], X_cand, n, k, "wald")
      Z_vector   <- qr_stats$Z_vector
      rho_vector <- qr_stats$rho_vector
    }

    Z_vector[!is.finite(Z_vector)] <- 0
    integrals    <- asymptotic_wald_marginal_gh(Z_vector, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Z_vector), a.priori[indices_remain_to_test])

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}

#' Exact Wald Test (ExW) — optimized
#' @description Exact selection by non-central t distribution (robust to high LD).
#' @param gh List produced by make_quad_nodes(tau = tau).
Exact_Wald <- function(X, Y, a.priori, K, tau, gh) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- numeric(K)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      r            <- cor(Y, X_cand)[1L, ]; r[is.na(r)] <- 0
      Z_hat_vector <- ifelse(abs(r) == 1, sign(r) * Inf, r * sqrt(n - 2L) / sqrt(1 - r^2))
      rho_vector   <- numeric(length(indices_remain_to_test))
    } else {
      qr_stats     <- compute_qr_stats(X, Y, X[, indices.inclus, drop = FALSE], X_cand, n, k, "wald")
      # Exact correction : Z_hat = Z × √(1-ρ) to align with the t assumptions
      Z_hat_vector <- qr_stats$Z_vector * sqrt(1 - qr_stats$rho_vector)
      rho_vector   <- qr_stats$rho_vector
    }

    Z_hat_vector[!is.finite(Z_hat_vector)] <- 0
    integrals    <- exact_wald_marginal_gh(Z_hat_vector, rho_vector, n, k, gh)
    integrals[!is.finite(integrals)] <- 0
    # H0 : t with n-k-1 degrees of freedom
    a.posteriori <- compute_posterior(integrals, dt(Z_hat_vector, df = n - k - 1L),
                                      a.priori[indices_remain_to_test])

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}

#' Asymptotic Score Test (AsS) — optimized
#' @description Selection by the efficient Score statistic.
#' @param gh List produced by make_quad_nodes(tau = tau).
Asymptotic_Score <- function(X, Y, a.priori, K, tau, gh) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- numeric(K)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]

    if (k == 1L) {
      num_vec      <- as.vector(crossprod(Y, X_cand))
      sigma2_null  <- sum(Y^2) / (n - 1L)
      Q_hat_vector <- num_vec / sqrt(colSums(X_cand^2) * sigma2_null)
      rho_vector   <- numeric(length(indices_remain_to_test))
    } else {
      qr_stats     <- compute_qr_stats(X, Y, X[, indices.inclus, drop = FALSE], X_cand, n, k, "score")
      Q_hat_vector <- qr_stats$Z_vector
      rho_vector   <- qr_stats$rho_vector
    }

    Q_hat_vector[!is.finite(Q_hat_vector)] <- 0
    integrals    <- asymptotic_score_marginal_gh(Q_hat_vector, rho_vector, n, gh)
    integrals[!is.finite(integrals)] <- 0
    a.posteriori <- compute_posterior(integrals, dnorm(Q_hat_vector),
                                      a.priori[indices_remain_to_test])

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}


# ==============================================================================
# Section 4 : Monte Carlo approximations for multi-SNP tests
# ==============================================================================

#' Fisher MC Approximation (Exact Multiple Test)
#' @description Approximates the posterior probability for the joint tests via
#'   Monte Carlo simulation of the non-central F distribution.
#' @param Q_W       Observed F statistic.
#' @param XtX_tilde Gram matrix of the candidate model.
#' @param XtX_k     Gram matrix of the null model (already included variables).
#' @param S         Number of Monte Carlo samples.
Fisher_MC_approximation <- function(Q_W, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA  <- ncol(XtX_tilde)
  df1 <- rA
  df2 <- n - (rA + 1)

  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  # ncp in R's convention (squared norm of the standardized mean), Lambda ~ N(0, tau2 I)
  delta_A  <- pmax(colSums(Lambda_A * (XtX_tilde %*% Lambda_A)), 0)
  densA    <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat    <- mean(densA)

  if (is.null(XtX_k) || ncol(XtX_k) == 0) {
    densB <- df(Q_W, df1 = df1, df2 = df2, ncp = 0)
    B_hat <- if (is.finite(densB)) densB else 0
  } else {
    rB       <- ncol(XtX_k)
    Lambda_B <- matrix(rnorm(rB * S, sd = sqrt(tau2)), nrow = rB)
    delta_B  <- pmax(colSums(Lambda_B * (XtX_k %*% Lambda_B)), 0)
    densB    <- df(Q_W, df1 = df1, df2 = df2, ncp = delta_B)
    densB[!is.finite(densB)] <- 0
    B_hat    <- mean(densB)
  }

  num <- p * A_hat
  den <- num + (1 - p) * B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

#' Posterior Wrapper — Exact Multiple Test
#' @param S Number of Monte Carlo samples; S = Inf computes the prior integrals
#'   exactly (Section 4b), which is the limit of the Monte Carlo as S grows.
ExM_posterior <- function(Q_W, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k)
  if (is.infinite(S)) {
    return(multiple_posterior_exact(Q_W, XtX_tilde, XtX_k, p, tau2,
                                    df2 = n - (ncol(XtX_tilde) + 1)))
  }
  Fisher_MC_approximation(Q_W, XtX_tilde, XtX_k, n, p, tau2, S)
}

#' Chi-Square MC Approximation (Naive Multiple Test)
#' @description Approximates the posterior probability via non-central Chi-2.
Chisq_MC_approximation <- function(Q_chi, XtX_tilde, XtX_k, n, p, tau2, S) {
  stopifnot(is.matrix(XtX_tilde))
  rA <- ncol(XtX_tilde)
  df <- rA

  Lambda_A <- matrix(rnorm(rA * S, sd = sqrt(tau2)), nrow = rA)
  # ncp in R's convention (squared norm of the standardized mean), Lambda ~ N(0, tau2 I)
  delta_A  <- pmax(colSums(Lambda_A * (XtX_tilde %*% Lambda_A)), 0)
  densA    <- dchisq(Q_chi, df = df, ncp = delta_A)
  densA[!is.finite(densA)] <- 0
  A_hat    <- mean(densA)

  if (is.null(XtX_k) || ncol(XtX_k) == 0L) {
    B_hat <- dchisq(Q_chi, df = df, ncp = 0)
    B_hat <- if (is.finite(B_hat)) B_hat else 0
  } else {
    rB       <- ncol(XtX_k)
    Lambda_B <- matrix(rnorm(rB * S, sd = sqrt(tau2)), nrow = rB)
    delta_B  <- pmax(colSums(Lambda_B * (XtX_k %*% Lambda_B)), 0)
    densB    <- dchisq(Q_chi, df = df, ncp = delta_B)
    densB[!is.finite(densB)] <- 0
    B_hat    <- mean(densB)
  }

  num <- p * A_hat
  den <- num + (1 - p) * B_hat
  if (!is.finite(den) || den <= 0) return(NA_real_)
  num / den
}

#' Posterior Wrapper — Naive Multiple Test
#' @param S Number of Monte Carlo samples; S = Inf computes the prior integrals exactly.
NaM_posterior <- function(Q_chi, X_tilde_jk, X_k, n, p, tau2, S) {
  XtX_tilde <- crossprod(X_tilde_jk)
  XtX_k     <- if (is.null(X_k)) NULL else crossprod(X_k)
  if (is.infinite(S)) {
    return(multiple_posterior_exact(Q_chi, XtX_tilde, XtX_k, p, tau2, df2 = Inf))
  }
  Chisq_MC_approximation(Q_chi, XtX_tilde, XtX_k, n, p, tau2, S)
}


# ==============================================================================
# Section 4b : Exact prior integrals for multi-SNP tests (S = Inf)
# ==============================================================================
# With Lambda ~ N(0, tau2 I) and delta = Lambda' G Lambda, averaging the Poisson
# mixture of the non-central law over the prior gives
#   E[f(x; delta)] = sum_r pi_r h_r(x),
# where pi is the law of N = sum_i NB(1/2, theta_i), theta_i = tau2 mu_i / (1 + tau2 mu_i),
# mu_i are the eigenvalues of G, and h_r is the central F (chi-square for NaM)
# density with k + 2r numerator degrees of freedom, rescaled as the non-central
# law requires.

.multiple_exact_memo <- new.env(parent = emptyenv())

#' theta_i = tau2 mu_i / (1 + tau2 mu_i) for the positive eigenvalues mu_i of G.
multiple_theta <- function(G, tau2) {
  if (is.null(G) || length(G) == 0L) return(numeric(0))
  mu <- pmax(eigen(G, symmetric = TRUE, only.values = TRUE)$values, 0)
  th <- tau2 * mu / (1 + tau2 * mu)
  th[th > 0]
}

#' log pi_r for r = 0..R, through r pi_r = (1/2) sum_{l <= r} p_l pi_{r-l} with
#' p_l = sum_i theta_i^l. All terms are positive, so there is no cancellation.
#' memo = TRUE caches the weights of the null model, which all candidates of a
#' step share.
nb_sum_logweights <- function(theta, R, memo = FALSE) {
  if (memo && identical(theta, .multiple_exact_memo$theta) &&
      length(.multiple_exact_memo$lw) > R) {
    return(.multiple_exact_memo$lw[seq_len(R + 1L)])
  }
  if (!length(theta)) return(c(0, rep(-Inf, R)))
  p <- colSums(outer(theta, seq_len(R), `^`))
  g <- numeric(R + 1L); g[1L] <- 1; off <- 0
  for (r in seq_len(R)) {
    g[r + 1L] <- sum(p[seq_len(r)] * g[r:1L]) / (2 * r)
    if (g[r + 1L] > 1e250) { g <- g * 1e-250; off <- off + 250 * log(10) }
  }
  lw <- 0.5 * sum(log1p(-theta)) + log(g) + off
  if (memo) {
    .multiple_exact_memo$theta <- theta
    .multiple_exact_memo$lw    <- lw
  }
  lw
}

#' log of the prior-integrated density E[f(x; delta)], with df2 = Inf for the
#' chi-square version. r -> log h_r(x) is concave, so once past its mode the
#' neglected tail is at most h_{R+1}(x) P(N > R); R grows until that bound is
#' below eps times the partial sum.
multiple_exact_logdens <- function(x, theta, df1, df2 = Inf, memo = FALSE, eps = 1e-12) {
  x <- max(x, .Machine$double.xmin)
  log_h <- function(r) {
    a <- df1 + 2 * r
    if (is.finite(df2)) log(df1 / a) + stats::df(df1 * x / a, a, df2, log = TRUE)
    else stats::dchisq(x, a, log = TRUE)
  }
  if (!length(theta)) return(log_h(0))
  half <- max(0, (if (is.finite(df2)) df1 * x else x) - df1) / 2
  R <- as.integer(ceiling(half + 8 * sqrt(half + 1) + 16))
  repeat {
    lw <- nb_sum_logweights(theta, R, memo)
    lh <- log_h(0:(R + 1L))
    lt <- lw + lh[seq_len(R + 1L)]
    mx <- max(lt)
    if (!is.finite(mx)) return(-Inf)
    ls <- mx + log(sum(exp(lt - mx)))
    tail_mass <- max(0, 1 - sum(exp(lw)))
    if (lh[R + 2L] < lh[R + 1L] && lh[R + 2L] + log(tail_mass) <= log(eps) + ls) return(ls)
    if (R >= 65536L) {
      warning("multiple_exact_logdens: truncation bound not reached")
      return(ls)
    }
    R <- 2L * R
  }
}

#' Exact posterior scores of all the candidates of a step at once (S = Inf).
#' @param stat     Joint statistics of the candidates: F (ExM, df2 finite) or
#'                 chi-square k F (NaM, df2 = Inf). Negative or non-finite
#'                 entries give NA, as in the Monte Carlo path.
#' @param X_inclus Columns already selected (NULL at the first step).
#' @param X_cand   Candidate columns, in the order of stat and p.
#' @details Under the prior the statistic of the candidate model is a sum of
#'   chi2_1 weighted by w_i = 1 + tau2 mu_i. For A, Ruben's representation is
#'   used with the scale beta = min_i w_i (theta_i = 1 - beta / w_i, kernels
#'   rescaled by beta): the number of terms is then about q / (2 beta) instead
#'   of q / 2, which matters for strong signals. The null model keeps beta = 1
#'   (its weights include the candidate coordinate, w = 1); its weights are
#'   shared by all candidates of the step and computed once.
multiple_posterior_exact_vec <- function(stat, X_inclus, X_cand, p, tau2, df2,
                                         eps = 1e-12, chunk = 2000L) {
  out    <- rep(NA_real_, length(stat))
  lo_out <- rep(NA_real_, length(stat))
  idx <- which(is.finite(stat) & stat >= 0)
  if (!length(idx)) return(structure(out, logodds = lo_out))
  k   <- if (is.null(X_inclus)) 1L else ncol(X_inclus) + 1L
  Gss <- if (k > 1L) crossprod(X_inclus) else NULL
  thB <- multiple_theta(Gss, tau2)

  # log of the density of beta * (central law with k + 2r degrees of freedom)
  # at x, for r = 0..Rmax+1 (rows) and the columns of x and beta.
  log_h_mat <- function(x, beta, Rmax) {
    a <- k + 2 * (0:(Rmax + 1L)); y <- x / beta
    if (is.finite(df2)) {
      u <- k * y / df2
      outer(lgamma((a + df2) / 2) - lgamma(a / 2), rep(1, length(y))) +
        outer(a / 2, log(u) - log1p(u)) +
        rep(-lgamma(df2 / 2) - log(y) - (df2 / 2) * log1p(u) - log(beta), each = length(a))
    } else {
      outer(-(a / 2) * log(2) - lgamma(a / 2), rep(1, length(y))) +
        outer(a / 2 - 1, log(y)) + rep(-y / 2 - log(beta), each = length(a))
    }
  }
  col_lse <- function(L) {
    mx <- apply(L, 2L, max); res <- rep(-Inf, ncol(L)); ok <- is.finite(mx)
    res[ok] <- mx[ok] + log(colSums(exp(L[, ok, drop = FALSE] - rep(mx[ok], each = nrow(L)))))
    res
  }
  n_terms <- function(q) {
    h <- pmax(0, q - k) / 2
    as.integer(ceiling(h + 8 * sqrt(h + 1) + 16))
  }
  # Stop once past the mode of r -> h_r (log-concave) with the neglected tail,
  # at most h_{R+1} P(N > R), below eps times the partial sum.
  certified <- function(LH, Rj, ls, tail) {
    hR1 <- LH[cbind(Rj + 2L, seq_along(Rj))]; hR <- LH[cbind(Rj + 1L, seq_along(Rj))]
    (hR1 < hR & hR1 + log(tail) <= log(eps) + ls) | !is.finite(ls) | Rj >= 65536L
  }

  # A: candidate model, one weight vector per column (columns sorted by R so
  # that those still active at step r of the recursion form a prefix).
  part_A <- function(Th, lp0, x, beta, R) {
    res <- rep(NA_real_, length(x)); todo <- seq_along(x)
    while (length(todo)) {
      o <- todo[order(R[todo], decreasing = TRUE)]; Rj <- R[o]; Rmax <- Rj[1L]
      Tj <- Th[, o, drop = FALSE]; J <- length(o)
      P  <- matrix(0, Rmax, J); Tl <- Tj
      for (l in seq_len(Rmax)) { P[l, ] <- colSums(Tl); Tl <- Tl * Tj }
      G <- matrix(0, Rmax + 1L, J); G[1L, ] <- 1; off <- numeric(J)
      for (r in seq_len(Rmax)) {
        s <- seq_len(sum(Rj >= r))
        G[r + 1L, s] <- colSums(P[seq_len(r), s, drop = FALSE] * G[r:1L, s, drop = FALSE]) / (2 * r)
        big <- which(G[r + 1L, s] > 1e250)
        if (length(big)) { G[, big] <- G[, big] * 1e-250; off[big] <- off[big] + 250 * log(10) }
      }
      LW <- log(G) + rep(lp0[o] + off, each = Rmax + 1L)
      LW[outer(0:Rmax, Rj, `>`)] <- -Inf
      LH <- log_h_mat(x[o], beta[o], Rmax)
      ls <- col_lse(LW + LH[seq_len(Rmax + 1L), , drop = FALSE])
      ok <- certified(LH, Rj, ls, pmax(0, 1 - colSums(exp(LW))))
      res[o[ok]] <- ls[ok]
      todo <- o[!ok]; R[todo] <- 2L * R[todo]
    }
    res
  }
  # B: null model, weights shared by the columns (beta = 1)
  part_B <- function(x, R) {
    res <- rep(NA_real_, length(x)); todo <- seq_along(x)
    while (length(todo)) {
      Rj <- R[todo]; Rmax <- max(Rj)
      lw <- nb_sum_logweights(thB, Rmax, memo = TRUE)
      LH <- log_h_mat(x[todo], rep(1, length(todo)), Rmax)
      L  <- lw + LH[seq_len(Rmax + 1L), , drop = FALSE]
      L[outer(0:Rmax, Rj, `>`)] <- -Inf
      ls <- col_lse(L)
      ok <- certified(LH, Rj, ls, pmax(0, 1 - cumsum(exp(lw))[Rj + 1L]))
      res[todo[ok]] <- ls[ok]
      nt <- todo[!ok]; R[nt] <- 2L * R[nt]; todo <- nt
    }
    res
  }

  for (ch in split(idx, ceiling(seq_along(idx) / chunk))) {
    x  <- pmax(stat[ch], .Machine$double.xmin)
    Xc <- X_cand[, ch, drop = FALSE]
    d  <- colSums(Xc^2)
    if (k == 1L) {
      W <- matrix(1 + tau2 * d, nrow = 1L)
    } else {
      Bm <- crossprod(X_inclus, Xc)
      W  <- vapply(seq_along(ch), function(j) {
        mu <- eigen(rbind(cbind(Gss, Bm[, j]), c(Bm[, j], d[j])),
                    symmetric = TRUE, only.values = TRUE)$values
        1 + tau2 * pmax(mu, 0)
      }, numeric(k))
    }
    beta <- apply(W, 2L, min)
    ThA  <- 1 - rep(beta, each = k) / W
    ThA[ThA < 0] <- 0
    lp0A <- 0.5 * colSums(log(rep(beta, each = k) / W))
    q    <- if (is.finite(df2)) k * x else x
    lA <- part_A(ThA, lp0A, x, beta, n_terms(q / beta))
    lB <- part_B(x, n_terms(q))
    pp <- p[ch]
    post <- 1 / (1 + (1 - pp) / pp * exp(lB - lA))
    post[pp <= 0] <- 0; post[pp >= 1] <- 1
    post[is.nan(post) | (lA == -Inf & lB == -Inf)] <- NA_real_
    lo <- log(pp) - log1p(-pp) + lA - lB
    lo[pp <= 0] <- -Inf; lo[pp >= 1] <- Inf
    lo[is.na(post)] <- NA_real_
    out[ch] <- post
    lo_out[ch] <- lo
  }
  structure(out, logodds = lo_out)
}

#' Posterior score p A / (p A + (1 - p) B) with A and B computed exactly.
multiple_posterior_exact <- function(stat, XtX_tilde, XtX_k, p, tau2, df2) {
  if (p <= 0) return(0)
  if (p >= 1) return(1)
  k  <- ncol(XtX_tilde)
  lA <- multiple_exact_logdens(stat, multiple_theta(XtX_tilde, tau2), k, df2)
  lB <- multiple_exact_logdens(stat, multiple_theta(XtX_k, tau2), k, df2, memo = TRUE)
  if (is.nan(lA) || is.nan(lB) || (lA == -Inf && lB == -Inf)) return(NA_real_)
  1 / (1 + (1 - p) / p * exp(lB - lA))
}


# ==============================================================================
# Section 5 : Multi-SNP (joint) tests
# ==============================================================================

#' Naive Multiple Test (NaM)
#' @description Joint test by Chi-2 approximation.
Naive_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- rep(NA_real_, K)
  TSS                    <- sum(Y^2)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1L
    if (df2 <= 0) break

    if (k == 1L) {
      r_sq          <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      F_block_vector <- r_sq * (n - 2) / (1 - r_sq)
      X_inclus       <- NULL
    } else {
      X_inclus      <- X[, indices.inclus, drop = FALSE]
      qr_inclus     <- qr(X_inclus)
      Y_res         <- qr.resid(qr_inclus, Y)
      X_cand_res    <- qr.resid(qr_inclus, X_cand)
      ss_cand       <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand     <- as.vector(t(Y_res) %*% X_cand_res) / ss_cand
      RSS_current   <- sum(Y_res^2)
      RSS_new       <- RSS_current - beta_cand^2 * ss_cand
      SSR_new       <- TSS - RSS_new
      F_block_vector <- (SSR_new / df1) / (RSS_new / df2)
    }

    Q_chi_vector <- rA * F_block_vector
    Q_chi_vector[!is.finite(Q_chi_vector) | Q_chi_vector < 0] <- -1
    p_vector     <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Q_chi_vector, X_inclus, X_cand, p_vector, tau2, df2 = Inf)
    } else mapply(
      function(q_chi, j, p) {
        if (q_chi < 0) return(NA_real_)
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        NaM_posterior(q_chi, X.modele, X_inclus, n, p, tau2, S)
      },
      Q_chi_vector, indices_remain_to_test, p_vector
    )

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}

#' Exact Multiple Test (ExM)
#' @description Joint test by exact F distribution (MC approximation).
Exact_Multiple <- function(X, Y, a.priori, K, tau2, S) {
  m <- ncol(X); n <- nrow(X)
  Y <- scale(Y, center = TRUE, scale = FALSE)
  X <- scale(X)

  indices.inclus         <- integer(0)
  indices_remain_to_test <- seq_len(m)
  post.probs             <- rep(NA_real_, K)
  TSS                    <- sum(Y^2)

  for (k in seq_len(K)) {
    if (length(indices_remain_to_test) == 0L) break
    X_cand <- X[, indices_remain_to_test, drop = FALSE]
    rA <- k; df1 <- rA; df2 <- n - rA - 1L
    if (df2 <= 0) break

    if (k == 1L) {
      r_sq     <- (as.vector(t(Y) %*% X_cand) / sqrt(TSS * colSums(X_cand^2)))^2
      r_sq[is.na(r_sq)] <- 0
      Z_vector <- r_sq * (n - 2) / (1 - r_sq)
      X_inclus <- NULL
    } else {
      X_inclus   <- X[, indices.inclus, drop = FALSE]
      qr_inclus  <- qr(X_inclus)
      Y_res      <- qr.resid(qr_inclus, Y)
      X_cand_res <- qr.resid(qr_inclus, X_cand)
      ss_cand    <- pmax(colSums(X_cand_res^2), 1e-9)
      beta_cand  <- as.vector(t(Y_res) %*% X_cand_res) / ss_cand
      RSS_current <- sum(Y_res^2)
      RSS_new    <- RSS_current - beta_cand^2 * ss_cand
      SSR_new    <- TSS - RSS_new
      Z_vector   <- (SSR_new / df1) / (RSS_new / df2)
    }

    Z_vector[!is.finite(Z_vector) | Z_vector < 0] <- -1
    p_vector <- a.priori[indices_remain_to_test]

    a.posteriori <- if (is.infinite(S)) {
      multiple_posterior_exact_vec(Z_vector, X_inclus, X_cand, p_vector, tau2, df2 = df2)
    } else mapply(
      function(z, j, p) {
        if (z < 0) return(NA_real_)
        X.modele <- if (is.null(X_inclus)) X[, j, drop = FALSE] else cbind(X_inclus, X[, j, drop = FALSE])
        ExM_posterior(z, X.modele, X_inclus, n, p, tau2, S)
      },
      Z_vector, indices_remain_to_test, p_vector
    )

    if (all(!is.finite(a.posteriori))) break
    best_local             <- which.best(a.posteriori)
    post.probs[k]          <- a.posteriori[best_local]
    a.inclure              <- indices_remain_to_test[best_local]
    indices.inclus         <- c(indices.inclus, a.inclure)
    indices_remain_to_test <- setdiff(indices_remain_to_test, a.inclure)
  }

  list(post.probs = post.probs, indices = indices.inclus)
}


# ==============================================================================
# Section 6 : Simulation and evaluation
# ==============================================================================

#' Simulation Wrapper
#' @description Runs parallel simulations to compare the type I error
#'   and the power between methods.
#' @param B         Number of simulation iterations per configuration.
#' @param S         Number of Monte Carlo samples for the multiple tests.
#' @param threshold Threshold for the selected variables (default 16).
Simulation <- function(X, p_values, tau_values, a.priori, alpha, K, B, S,
                       threshold = 16) {
  m <- ncol(X); n <- nrow(X)
  num_cores <- detectCores() - 1
  cl        <- makeCluster(num_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    library(doParallel); library(MASS); library(mvtnorm)
    library(stats);      library(statmod)
  })

  clusterExport(cl, c(
    "Generate.Y", "Power_Fun",
    "make_quad_nodes", "compute_qr_stats", "compute_posterior",
    "naive_wald_marginal", "asymptotic_wald_marginal_gh",
    "exact_wald_marginal_gh", "asymptotic_score_marginal_gh",
    "Naive_Wald", "Asymptotic_Wald", "Exact_Wald", "Asymptotic_Score",
    "Naive_Multiple", "Exact_Multiple",
    "Fisher_MC_approximation", "Chisq_MC_approximation",
    "ExM_posterior", "NaM_posterior",
    "multiple_theta", "nb_sum_logweights", "multiple_exact_logdens",
    "multiple_posterior_exact", "multiple_posterior_exact_vec", ".multiple_exact_memo"
  ))

  param_grid           <- expand.grid(p = p_values, tau = tau_values, stringsAsFactors = FALSE)
  final_main_results   <- list()
  final_sel_summaries  <- list()
  method_names         <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM")

  for (i in seq_len(nrow(param_grid))) {
    p   <- param_grid$p[i]
    tau <- param_grid$tau[i]

    # GH nodes pre-computed once per configuration (tau fixed within the loop)
    gh <- make_quad_nodes(n_nodes = 30L, tau = tau)

    b_val <- uniroot(Power_Fun, lower = 0, upper = 1, n = n, alpha = alpha, p = p)$root
    beta  <- c(1, rep(b_val, 15), rep(0, m - 15))

    mc_results <- foreach(b = seq_len(B), .combine = rbind,
                          .packages = c("stats", "statmod")) %dopar% {
      Y <- Generate.Y(n, X, beta, sigma = 1)
      gh_worker <- make_quad_nodes(30L, tau)

      methods_results <- list(
        NaW = Naive_Wald(X, Y, a.priori, K, tau),
        AsW = Asymptotic_Wald(X, Y, a.priori, K, tau, gh_worker),
        ExW = Exact_Wald(X, Y, a.priori, K, tau, gh_worker),
        AsS = Asymptotic_Score(X, Y, a.priori, K, tau, gh_worker),
        NaM = Naive_Multiple(X, Y, a.priori, K, tau^2, S),
        ExM = Exact_Multiple(X, Y, a.priori, K, tau^2, S)
      )

      c(
        sapply(methods_results, function(res) sum(res$indices < threshold)),
        sapply(methods_results, function(res) sum(res$indices %in% 1:15))
      )
    }

    # Aggregation of the main results
    mean_counts  <- colMeans(mc_results[, seq_along(method_names)])
    mean_results <- data.frame(p = p, tau = tau)
    for (j in seq_along(method_names))
      mean_results[[paste0("mean_", method_names[j])]] <- mean_counts[j]
    final_main_results[[i]] <- mean_results

    # Aggregation of the selection summaries
    for (j in seq_along(method_names)) {
      method_col  <- mc_results[, length(method_names) + j]
      sel_summary <- tabulate(method_col + 1L, nbins = 16L) / B
      summary_df  <- as.data.frame(t(sel_summary))
      colnames(summary_df) <- paste0("choix.", 0:15)
      summary_df$p <- p; summary_df$tau <- tau; summary_df$method <- method_names[j]
      final_sel_summaries[[length(final_sel_summaries) + 1L]] <- summary_df
    }
  }

  stopCluster(cl)

  main_results      <- do.call(rbind, final_main_results)
  selection_summary <- do.call(rbind, final_sel_summaries)
  row.names(main_results) <- row.names(selection_summary) <- NULL

  list(main_results = main_results, selection_summary = selection_summary)
}


# ==============================================================================
# Section 7 : Comparison of the methods
# ==============================================================================

#' Superiority probability
#' @description Computes P(method X selects more true variables than Y).
Get_Superior_Probability <- function(x, y) {
  sum(sapply(0:15, function(i) {
    xi <- x[i + 1]
    sum(sapply(i:15, function(j) {
      yj     <- y[j + 1]
      weight <- if (i == j) 0.5 else 1
      xi * yj * weight
    }))
  }))
}

#' Pairwise comparison matrix
#' @description Generates the comparison matrices for all the methods.
Compare_Methods_Matrix <- function(selection_summary) {
  methods      <- unique(selection_summary$method)
  param_groups <- unique(selection_summary[, c("p", "tau")])
  all_matrices <- list()

  for (i in seq_len(nrow(param_groups))) {
    p_val   <- param_groups$p[i]
    tau_val <- param_groups$tau[i]

    subset_data   <- selection_summary %>% filter(p == p_val, tau == tau_val)
    matrix_result <- matrix(NA, nrow = length(methods), ncol = length(methods),
                            dimnames = list(methods, methods))

    for (m1 in methods) {
      for (m2 in methods) {
        if (m1 != m2) {
          x <- as.numeric(subset_data %>% filter(method == m1) %>% select(starts_with("choix.")))
          y <- as.numeric(subset_data %>% filter(method == m2) %>% select(starts_with("choix.")))
          matrix_result[m1, m2] <- round(Get_Superior_Probability(x, y), 3)
        }
      }
    }
    all_matrices[[paste0("p=", p_val, "_tau=", tau_val)]] <- matrix_result
  }
  all_matrices
}
