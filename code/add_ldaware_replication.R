# =============================================================================
# File: add_ldaware_replication.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Post-hoc computation of the LD-aware replication rate on an existing
# Benchmark_GEUVADIS.R real-mode result, without re-running any selection:
# a selected SNP receives credit when it is in YRI genotype LD
# (|r| >= r_clust) with at least one SNP that itself replicates
# (sign-consistent with EUR and nominally significant in YRI). The stored
# `selections` (positions) and `marginals` (r_eur / r_yri) are reused; YRI
# genotypes come from the per-gene RDS.
#
# Usage (from code/):
#   Rscript add_ldaware_replication.R \
#       --results ../results_geuvadis_real.rds --rds_dir ../geuvadis_rds \
#       --out ../results_geuvadis_real_ld.rds
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

cfg <- list(results = "../results_geuvadis_real.rds",
            rds_dir = "../geuvadis_rds",
            r_clust = 0.7, rep_alpha = 0.05,
            bonferroni = FALSE,   # TRUE: per-gene threshold rep_alpha / m (window-wide Bonferroni)
            perm = 0L,            # >0: per-locus permutation threshold (this many
                                  # permutations of the residualized YRI phenotype;
                                  # threshold = 95th pct of the max |t| null)
            perm_seed = 2026L,
            out = "../results_geuvadis_real_ld.rds")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

res <- readRDS(cfg$results)
stopifnot(res$mode == "real")
genes <- res$genes

# FWL residualization (same helper as Benchmark_GEUVADIS.R), needed for the
# permutation null, which must be built on the same residualized YRI data as
# the stored marginal correlations r_yri.
residualize <- function(M, C) {
  q <- qr(cbind(1, C))
  qr.resid(q, as.matrix(M))
}

rows <- list()
for (i in seq_along(genes)) {
  g  <- genes[i]
  mg <- res$marginals[[i]]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
  Gv <- gx$G_yri[, match(mg$pos, gx$pos), drop = FALSE]   # align to filtered pos
  n_v <- nrow(Gv)
  # Exact residual degrees of freedom of the marginal YRI test: the stored
  # r_yri are computed on data residualized on [1, covar_yri], so the t test
  # has n - rank([1, C]) - 1 degrees of freedom, not n - 2.
  df_v <- n_v - qr(cbind(1, as.matrix(gx$covar_yri)))$rank - 1L
  if (cfg$perm > 0L) {
    # Per-locus permutation threshold: 95th percentile of the max |t| over the
    # window under permutations of the residualized YRI phenotype.
    yv <- as.numeric(residualize(gx$y_yri, gx$covar_yri))
    Xv <- residualize(gx$G_yri, gx$covar_yri)[, match(mg$pos, gx$pos), drop = FALSE]
    Xs <- scale(Xv); Xs[is.na(Xs)] <- 0
    set.seed(cfg$perm_seed + i)
    Yp <- scale(vapply(seq_len(cfg$perm), function(p) sample(yv), numeric(n_v)))
    rmax <- apply(abs(crossprod(Xs, Yp)) / (n_v - 1), 2L, max, na.rm = TRUE)
    rmax <- pmin(rmax, 1 - 1e-12)
    t_crit <- quantile(rmax * sqrt(df_v / (1 - rmax^2)), 0.95, names = FALSE)
  } else {
    alpha_g <- if (isTRUE(cfg$bonferroni)) cfg$rep_alpha / length(mg$pos) else cfg$rep_alpha
    t_crit <- qt(1 - alpha_g / 2, df = df_v)
  }
  t_all  <- mg$r_yri * sqrt(df_v / pmax(1 - mg$r_yri^2, 1e-12))
  rep_set <- which((sign(mg$r_yri) == sign(mg$r_eur)) & (abs(t_all) > t_crit))

  for (key in names(res$selections[[i]])) {
    sel_pos <- res$selections[[i]][[key]]
    parts <- strsplit(key, ".", fixed = TRUE)[[1L]]
    idx <- match(sel_pos, mg$pos)
    idx <- idx[!is.na(idx)]
    rr_ld <- if (!length(idx)) NA_real_ else if (!length(rep_set)) 0 else {
      ldC <- suppressWarnings(abs(cor(Gv[, idx, drop = FALSE],
                                      Gv[, rep_set, drop = FALSE])))
      ldC[is.na(ldC)] <- 0
      mean(rowSums(ldC >= cfg$r_clust) > 0)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, prior = parts[1], method = parts[2], selection = parts[3],
      rep_rate_ld = rr_ld, stringsAsFactors = FALSE)
  }
}
ld <- do.call(rbind, rows)
# Recent Benchmark_GEUVADIS.R runs store an inline rep_rate_ld; drop it so the
# join does not produce .x/.y suffixes (this script recomputes it, and is the
# only path that supports the Bonferroni variant).
tab <- res$table %>% select(-any_of("rep_rate_ld")) %>%
  left_join(ld, by = c("gene", "prior", "method", "selection"))
res$table <- tab
res$config$r_clust_ld <- cfg$r_clust
saveRDS(res, cfg$out)
cat(sprintf("Wrote %s (%d rows, rep_rate_ld non-NA: %d)\n\n",
            cfg$out, nrow(tab), sum(!is.na(tab$rep_rate_ld))))

s <- tab %>% filter(prior != "none") %>% group_by(prior, method, selection) %>%
  summarise(K = round(mean(K), 1),
            rep_rate = round(mean(rep_rate, na.rm = TRUE), 3),
            rep_rate_ld = round(mean(rep_rate_ld, na.rm = TRUE), 3),
            R2_yri = round(mean(R2_yri, na.rm = TRUE), 3), .groups = "drop")
cat("=== fixed-K (K=5), mean over genes ===\n")
print(as.data.frame(s %>% filter(selection == "fixedK") %>% arrange(prior, desc(rep_rate_ld))), row.names = FALSE)
cat("\n=== adaptive BIC ===\n")
print(as.data.frame(s %>% filter(selection == "BIC") %>% arrange(prior, desc(rep_rate_ld))), row.names = FALSE)
cat("\n=== gradient (in-house methods, fixed-K) ===\n")
g <- tab %>% filter(selection == "fixedK", method %in% c("NaW","AsW","ExW","AsS","NaM","ExM")) %>%
  group_by(prior) %>%
  summarise(rep_rate = round(mean(rep_rate, na.rm = TRUE), 3),
            rep_rate_ld = round(mean(rep_rate_ld, na.rm = TRUE), 3), .groups = "drop")
print(as.data.frame(g), row.names = FALSE)
