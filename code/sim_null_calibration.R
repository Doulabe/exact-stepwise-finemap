# =============================================================================
# File: sim_null_calibration.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Calibration of the exact conditional laws under the global null, along a
# data-driven selection path. At each step k, the exact conditional p-value of
# the SELECTED candidate (the best of the m - k + 1 remaining ones) is compared
# with that of a data-independent reference candidate evaluated given the same
# selected set, under the exact t law and under the standard normal
# approximation shared by the naive and asymptotic variants. The exact
# reference p-values behave like Uniform(0,1), confirming the conditional
# exactness; the normal-approximation reference p-values exceed the nominal
# level at small n (the miscalibration the exact laws remove); the selected
# p-values concentrate near zero, quantifying the selection effect that no
# candidate-level law corrects.
#
# Design: AR(1) blocks (m = 1000, blocks of 50, rho = 0.95), beta = 0,
# B = 500 data sets, K = 10 steps, n in {50, 80, 1000}. The reference candidate
# is the centre of a fixed block and is excluded from selection.
#
# Usage (from code/):
#   Rscript sim_null_calibration.R --out ../tab_null_calibration.tex
# =============================================================================
suppressPackageStartupMessages({source("Generate_XY_realistic.R")})

cfg <- list(m = 1000L, block = 50L, rho = 0.95, B = 500L, K = 10L,
            ref = 525L, seed = 2026L, out = "../tab_null_calibration.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

run_n <- function(n) {
  p_sel <- p_ref <- p_ref_norm <- matrix(NA_real_, cfg$B, cfg$K)
  for (b in seq_len(cfg$B)) {
    set.seed(cfg$seed + 7919L * b + n)
    X <- Generate.X.realistic(n = n, m = cfg$m, rho_intra = cfg$rho,
                              rho_inter = 0, block_size = cfg$block)
    X <- scale(X); X[is.nan(X)] <- 0
    Y <- rnorm(n)
    sel <- integer(0)
    Yr <- Y - mean(Y); Xr <- X
    for (k in seq_len(cfg$K)) {
      df <- n - k - 1L
      ss <- pmax(colSums(Xr^2), 1e-12)
      tt <- as.vector(crossprod(Xr, Yr)) / sqrt(ss)
      tt <- tt / sqrt(pmax((sum(Yr^2) - tt^2) / df, 1e-300) )
      tt[c(sel, cfg$ref)] <- NA
      j <- which.max(abs(tt))
      p_sel[b, k] <- 2 * pt(-abs(tt[j]), df)
      tr <- as.vector(crossprod(Xr[, cfg$ref], Yr)) / sqrt(ss[cfg$ref])
      tr <- tr / sqrt(pmax((sum(Yr^2) - tr^2) / df, 1e-300))
      p_ref[b, k] <- 2 * pt(-abs(tr), df)
      # the naive and asymptotic variants share the standard normal null law
      p_ref_norm[b, k] <- 2 * pnorm(-abs(tr))
      v <- Xr[, j] / sqrt(ss[j])
      Yr <- Yr - v * sum(v * Yr)
      Xr <- Xr - outer(v, as.vector(crossprod(Xr, v)))
      sel <- c(sel, j)
    }
  }
  list(sel = p_sel, ref = p_ref, refn = p_ref_norm)
}

NS <- c(50L, 80L, 1000L)
res <- lapply(NS, run_n)
names(res) <- as.character(NS)

steps <- c(1L, 2L, 5L, 10L)
fmt <- function(x) sprintf("%.3f", x)
rows <- vapply(steps, function(k) {
  cells <- unlist(lapply(as.character(NS), function(nn) c(
    fmt(mean(res[[nn]]$sel[, k]  < 0.05)),
    fmt(mean(res[[nn]]$ref[, k]  < 0.05)),
    fmt(mean(res[[nn]]$refn[, k] < 0.05)))))
  paste0(k, " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Calibration under the global null along a data-driven path.} ",
         "Empirical frequency of conditional $p$-values below $0.05$ at step $k$ of a null ",
         "stepwise scan (AR(1) design, $m=1{,}000$, $\\rho=0.95$, $\\beta = 0$, $B=500$ data sets), ",
         "for the candidate selected at that step under the exact $t$ law (sel.), and for a fixed, ",
         "data-independent reference candidate evaluated given the same selected set, under the ",
         "exact $t$ law (ref.\\ $t$) and under the standard normal approximation shared by the ",
         "naive and asymptotic variants (ref.\\ $\\mathcal{N}$). The exact reference frequencies ",
         "stay at the nominal level, confirming the conditional exactness of the laws. The ",
         "normal-approximation frequencies slightly exceed it at small $n$, a mild null-side ",
         "miscalibration at these residual degrees of freedom, the consequential difference ",
         "between the kernels appearing under the alternative (Supplementary ",
         "Figure~\\ref{fig:premium_by_prior}). The selected frequencies quantify the ",
         "selection effect, which no candidate-level law corrects. The scores are for ",
         "prioritization, not post-selection testing.}}"),
  "\\label{tab:null_calibration}",
  "\\centering\\scriptsize",
  "\\setlength{\\tabcolsep}{3.5pt}",
  "\\begin{tabular}{lccccccccc}", "\\toprule",
  paste0(" & ", paste(sprintf("\\multicolumn{3}{c}{$n=%d$}", NS), collapse = " & "), " \\\\"),
  "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
  paste0("Step $k$", strrep(" & sel. & ref.\\ $t$ & ref.\\ $\\mathcal{N}$", length(NS)), " \\\\"),
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
for (nn in names(res)) cat(sprintf(
  "n=%s: ref-t p<.05 = %.3f | ref-norm p<.05 = %.3f | mean ref-t p = %.3f\n",
  nn, mean(res[[nn]]$ref < 0.05), mean(res[[nn]]$refn < 0.05), mean(res[[nn]]$ref)))
