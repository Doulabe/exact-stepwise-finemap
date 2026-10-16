# Prototype: exact prior-predictive densities of the joint Multiple-SNP statistic.
#   A = E_Lambda[ f_{k,nu}(z; Lambda' G Lambda) ],  Lambda ~ N(0, tau^2 I_k)
# Representation: A = sum_r pi_r h_r(z), pi = law of N = sum_i NB(1/2, theta_i),
#   theta_i = tau^2 mu_i / (1 + tau^2 mu_i), mu_i = eigenvalues of G,
#   h_r(z) = (k/(k+2r)) f_{k+2r,nu}(k z/(k+2r))  (ExM, central F)
#          = f_{chi2_{k+2r}}(q)                  (NaM, q = k z)
suppressPackageStartupMessages(library(CompQuadForm))
# Run FROM code/ (sources the repository implementations below).
source("Generate_XY_realistic.R")

nb_logpmf <- function(theta, R) {
  theta <- theta[theta > 0]
  lp0 <- 0.5 * sum(log1p(-theta))
  out <- rep(-Inf, R + 1L); out[1L] <- lp0
  if (!length(theta) || R == 0L) return(out)
  p <- vapply(seq_len(R), function(l) sum(theta^l), 0)
  g <- numeric(R + 1L); g[1L] <- 1
  for (r in seq_len(R)) g[r + 1L] <- sum(p[seq_len(r)] * g[r:1]) / (2 * r)
  lp0 + log(g)
}
log_h <- function(x, k, r, nu) {
  a <- k + 2 * r
  if (is.finite(nu)) log(k / a) + df(k * x / a, a, nu, log = TRUE) else dchisq(x, a, log = TRUE)
}
# log of the exact prior-predictive density, relative truncation error <= eps
exact_logdens <- function(x, theta, k, nu = Inf, eps = 1e-12) {
  R <- 64L
  repeat {
    lh <- log_h(x, k, 0:R, nu)
    lt <- nb_logpmf(theta, R) + lh
    mx <- max(lt); ls <- mx + log(sum(exp(lt - mx)))
    if (lh[R + 1L] < lh[R] && log_h(x, k, R + 1L, nu) <= log(eps) + ls)
      return(c(logdens = ls, R = R))
    R <- 2L * R
    if (R > 2^20) stop("no convergence")
  }
}
mc_dens <- function(x, G, tau, k, nu, S, chunk = 2e5) {
  vals <- numeric(0)
  while (length(vals) < S) {
    s <- min(chunk, S - length(vals))
    L <- matrix(rnorm(nrow(G) * s, sd = tau), nrow(G))
    d <- pmax(colSums(L * (G %*% L)), 0)
    vals <- c(vals, if (is.finite(nu)) df(x, k, nu, ncp = d) else dchisq(x, k, ncp = d))
  }
  c(mean = mean(vals), relse1 = sd(vals) / mean(vals))   # relse1: relative SE for S = 1
}

