# =============================================================================
# File: figure_budget_curve.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Replication-versus-budget curve for the GEUVADIS application: mean LD-aware
# replication rate of the first K' selected variants, K' = 1..5, for the joint
# stepwise scores (ExM, NaM), SuSiE and LASSO, under the TSS-proximity prior
# (LASSO is prior-free). Selections are stored in path / PIP / coefficient
# order, so the curve needs no re-run; LASSO retains K_hat <= 5 variants at the
# one-standard-error penalty, so its curve is computed on the first
# min(K', K_hat) variants and plateaus at its own selection size.
#
# Usage (from code/):
#   Rscript figure_budget_curve.R \
#     --results ../../server_rerun_fix_20260928/results_exact/results_geuvadis_real_ld.rds \
#     --rds_dir ../geuvadis_rds --prior tss --out ../Figure_budget_curve.pdf
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            rds_dir = "../geuvadis_rds",
            prior   = "tss",
            r_clust = 0.7, rep_alpha = 0.05, Kmax = 5L,
            out     = "../Figure_budget_curve.pdf")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

res   <- readRDS(cfg$results)
genes <- res$genes
METH  <- c(ExM = paste0(cfg$prior, ".ExM.fixedK"), NaM = paste0(cfg$prior, ".NaM.fixedK"),
           SuSiE = paste0(cfg$prior, ".SuSiE.fixedK"), LASSO = "none.LASSO.fixedK")

rate <- array(NA_real_, dim = c(length(genes), length(METH), cfg$Kmax),
              dimnames = list(genes, names(METH), NULL))
for (i in seq_along(genes)) {
  mg <- res$marginals[[i]]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", genes[i])))
  Gv <- gx$G_yri[, match(mg$pos, gx$pos), drop = FALSE]
  n_v <- nrow(Gv)
  t_crit <- qt(1 - cfg$rep_alpha / 2, df = n_v - 2)
  t_all  <- mg$r_yri * sqrt((n_v - 2) / pmax(1 - mg$r_yri^2, 1e-12))
  rep_set <- which((sign(mg$r_yri) == sign(mg$r_eur)) & (abs(t_all) > t_crit))
  for (mn in names(METH)) {
    sel <- match(res$selections[[i]][[METH[[mn]]]], mg$pos)
    sel <- sel[!is.na(sel)]
    if (!length(sel)) next
    cred <- if (length(rep_set)) {
      ldC <- suppressWarnings(abs(cor(Gv[, sel, drop = FALSE], Gv[, rep_set, drop = FALSE])))
      ldC[is.na(ldC)] <- 0
      rowSums(ldC >= cfg$r_clust) > 0
    } else rep(FALSE, length(sel))
    for (K in seq_len(cfg$Kmax)) {
      kk <- min(K, length(sel))
      rate[i, mn, K] <- mean(cred[seq_len(kk)])
    }
  }
}
curve <- apply(rate, c(2, 3), mean, na.rm = TRUE)
cat("Mean LD-aware replication rate by budget K' (prior:", cfg$prior, ")\n")
print(round(curve, 3))

pdf(cfg$out, width = 5.6, height = 4.2, pointsize = 10)
par(mar = c(4, 4, 1, 1) + 0.1)
cols <- c(ExM = "#D55E00", NaM = "#E69F00", SuSiE = "#0072B2", LASSO = "#009E73")
pchs <- c(ExM = 16, NaM = 17, SuSiE = 15, LASSO = 18)
matplot(seq_len(cfg$Kmax), t(curve), type = "o", lty = 1, lwd = 1.6,
        pch = pchs[rownames(curve)], col = cols[rownames(curve)],
        xlab = "Selection budget K' (first K' selected variants)",
        ylab = "Mean LD-aware replication rate", ylim = c(0.55, 0.9), xaxt = "n")
axis(1, at = seq_len(cfg$Kmax))
legend("topright", legend = rownames(curve), col = cols[rownames(curve)],
       pch = pchs[rownames(curve)], lty = 1, lwd = 1.6, bty = "n")
invisible(dev.off())
cat("Wrote", cfg$out, "\n")
