# Validation table for the note: exact series vs independent inversion and vs Monte Carlo,
# on one set of simulated designs (AR(1) block, rho = 0.95, half of the k effects at power 0.6).
suppressPackageStartupMessages(library(CompQuadForm))
source("exact_multiple_funs.R")   # run FROM code/
gp_dens_F <- function(x, w, k, nu) {
  cc <- k * x / nu
  f <- function(u) vapply(u, function(uu) Re(exp(-0.5 * sum(log(1 - 2i * uu * w)) - (nu / 2 + 1) * log(1 + 2i * uu * cc))), 0)
  k * integrate(f, 0, Inf, subdivisions = 5000L, rel.tol = 1e-10)$value / pi
}
davies_dens <- function(q, w) {
  if (length(w) == 1L) return(dchisq(q / w, 1) / w)
  h <- 1e-4 * q; a1 <- davies(q - h, w, acc = 1e-11, lim = 1e7); a2 <- davies(q + h, w, acc = 1e-11, lim = 1e7)
  if (a1$ifault != 0 || a2$ifault != 0) NA_real_ else (a1$Qq - a2$Qq) / (2 * h)
}
set.seed(20260930)
cases <- expand.grid(k = c(1, 3, 8, 15), n = c(80, 1000), tau = c(0.2, 0.5))
rows <- list()
for (i in seq_len(nrow(cases))) {
  n <- cases$n[i]; k <- cases$k[i]; tau <- cases$tau[i]
  X <- scale(Generate.X.realistic(n = n, m = 100, rho_intra = 0.95, block_size = 50)); X[is.nan(X)] <- 0
  XA <- X[, sort(sample(ncol(X), k)), drop = FALSE]
  lam <- rep(0, k); lam[seq_len(max(1, k %/% 2))] <- (qnorm(0.975) + qnorm(0.6)) / sqrt(n)
  Y <- as.vector(XA %*% lam + rnorm(n)); Y <- Y - mean(Y)
  fit <- lm.fit(XA, Y); RSS <- sum(fit$residuals^2); SSR <- sum(Y^2) - RSS
  nu <- n - k - 1; z <- (SSR / k) / (RSS / nu); q <- k * z
  G <- crossprod(XA); mu <- pmax(eigen(G, symmetric = TRUE, only.values = TRUE)$values, 0)
  th <- tau^2 * mu / (1 + tau^2 * mu); w <- 1 + tau^2 * mu
  eF <- exact_logdens(z, th, k, nu); eC <- exact_logdens(q, th, k, Inf)
  AF <- exp(eF[["logdens"]]); AC <- exp(eC[["logdens"]])
  mF <- mc_dens(z, G, tau, k, nu, 1e6)
  rows[[i]] <- data.frame(n = n, k = k, tau = tau, ntau2 = n * tau^2, z = z, R = eF[["R"]], ExM_exact = AF,
    ExM_vs_inv = gp_dens_F(z, w, k, nu) / AF - 1, NaM_vs_davies = davies_dens(q, w) / AC - 1,
    MC1e6_vs_exact = mF[["mean"]] / AF - 1, MC_relSE_S100 = mF[["relse1"]] / 10)
}
res <- do.call(rbind, rows)
options(width = 200); print(res, digits = 3, row.names = FALSE)
write.csv(res, "../exact_multiple_validation.csv", row.names = FALSE)
