# =============================================================================
# File: make_qtltools_eval.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Evaluates the QTLtools conditional selections (run_qtltools_geuvadis100.R)
# with the replication machinery of the main analysis: LD-aware credit in YRI
# under the nominal test and under the per-window permutation threshold (same
# seeds as add_ldaware_replication.R, so the thresholds are identical to those
# of Table 9), plus the redundancy diagnostics (distinct LD clusters, mean
# within-selection |r|, cluster-level replication).
#
# Usage (from code/):
#   Rscript make_qtltools_eval.R --out ../results_qtltools_eval.rds
# =============================================================================
suppressPackageStartupMessages({library(dplyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            qtl = "../results_qtltools_geuvadis.rds",
            rds_dir = "../geuvadis_rds", r_clust = 0.7, rep_alpha = 0.05,
            perm = 500L, perm_seed = 2026L,
            out = "../results_qtltools_eval.rds")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

residualize <- function(M, C) {
  q <- qr(cbind(1, as.matrix(C)))
  qr.resid(q, as.matrix(M))
}
greedy_clusters <- function(C, thr) {
  K <- nrow(C); id <- integer(K)
  for (i in seq_len(K)) {
    if (id[i] > 0L) next
    id[i] <- max(id) + 1L
    for (j in seq_len(K)) if (id[j] == 0L && C[i, j] >= thr) id[j] <- id[i]
  }
  id
}

res <- readRDS(cfg$results)
qtl <- readRDS(cfg$qtl)$results
genes <- res$genes

rows <- list()
for (i in seq_along(genes)) {
  g  <- genes[i]
  qr_ <- qtl[[g]]
  if (is.null(qr_)) next
  mg <- res$marginals[[i]]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
  colsel <- match(mg$pos, gx$pos)
  Gv <- gx$G_yri[, colsel, drop = FALSE]
  Ge <- gx$G_eur[, colsel, drop = FALSE]
  n_v <- nrow(Gv)
  df_v <- n_v - qr(cbind(1, as.matrix(gx$covar_yri)))$rank - 1L
  t_all <- mg$r_yri * sqrt(df_v / pmax(1 - mg$r_yri^2, 1e-12))
  sgn <- sign(mg$r_yri) == sign(mg$r_eur)
  # nominal threshold
  t_nom <- qt(1 - cfg$rep_alpha / 2, df = df_v)
  # permutation threshold, identical seeds to add_ldaware_replication.R
  yv <- as.numeric(residualize(gx$y_yri, gx$covar_yri))
  Xv <- residualize(gx$G_yri, gx$covar_yri)[, colsel, drop = FALSE]
  Xs <- scale(Xv); Xs[is.na(Xs)] <- 0
  set.seed(cfg$perm_seed + i)
  Yp <- scale(vapply(seq_len(cfg$perm), function(p) sample(yv), numeric(n_v)))
  rmax <- apply(abs(crossprod(Xs, Yp)) / (n_v - 1), 2L, max, na.rm = TRUE)
  rmax <- pmin(rmax, 1 - 1e-12)
  t_perm <- quantile(rmax * sqrt(df_v / (1 - rmax^2)), 0.95, names = FALSE)

  idx <- match(qr_$sel_pos, mg$pos)
  idx <- idx[!is.na(idx)]
  K <- length(idx)
  cred <- function(t_crit) {
    rep_set <- which(sgn & (abs(t_all) > t_crit))
    if (!K) return(c(rate = NA_real_, clrep = NA_real_))
    if (!length(rep_set)) return(c(rate = 0, clrep = 0))
    ldC <- suppressWarnings(abs(cor(Gv[, idx, drop = FALSE],
                                    Gv[, rep_set, drop = FALSE])))
    ldC[is.na(ldC)] <- 0
    ok <- rowSums(ldC >= cfg$r_clust) > 0
    if (K >= 2L) {
      Ce <- suppressWarnings(abs(cor(Ge[, idx, drop = FALSE]))); Ce[is.na(Ce)] <- 0
      cl <- greedy_clusters(Ce, cfg$r_clust)
    } else cl <- rep(1L, K)
    c(rate = mean(ok),
      clrep = mean(vapply(unique(cl), function(cc) any(ok[cl == cc]), TRUE)))
  }
  cn <- cred(t_nom); cp <- cred(t_perm)
  n_cl <- if (K >= 2L) {
    Ce <- suppressWarnings(abs(cor(Ge[, idx, drop = FALSE]))); Ce[is.na(Ce)] <- 0
    length(unique(greedy_clusters(Ce, cfg$r_clust)))
  } else K
  rows[[length(rows) + 1L]] <- data.frame(
    gene = g, K = K, rep_rate_ld = cn["rate"], rep_rate_ld_perm = cp["rate"],
    n_clusters = n_cl, cluster_rep = cn["clrep"])
  if (i %% 25 == 0) cat(i, "genes\n")
}
d <- bind_rows(rows)
saveRDS(list(table = d), cfg$out)
cat("Wrote", cfg$out, "\n")
cat(sprintf("QTLtools: K^ = %.2f | LD-aware nominal = %.3f | perm = %.3f | clusters = %.2f | cluster-rep = %.3f | K=0 on %d genes\n",
            mean(d$K), mean(d$rep_rate_ld, na.rm = TRUE),
            mean(d$rep_rate_ld_perm, na.rm = TRUE), mean(d$n_clusters),
            mean(d$cluster_rep, na.rm = TRUE), sum(d$K == 0)))
