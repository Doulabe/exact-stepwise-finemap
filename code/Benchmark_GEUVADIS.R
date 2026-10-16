#!/usr/bin/env Rscript
# =============================================================================
# File: Benchmark_GEUVADIS.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# GEUVADIS cis-eQTL driver. Consumes the per-gene RDS built by
# prep_geuvadis_gene.R (EUR discovery / YRI replication split, four priors
# spanning a fidelity gradient: uniform / polyfun / tss / discordant).
#
# Two modes:
#
#   --mode real (default) — real expression phenotype.
#     Discovery on EUR under each prior: the 6 in-house methods (fixed-K and,
#     optionally, adaptive stopping via *_Criteria_For_Data) + SuSiE (flat for
#     the uniform arm, prior-weighted otherwise) + LASSO (prior-free).
#     Evaluation = cross-population replication in YRI, whose different LD
#     penalizes selections that merely tag an EUR-specific proxy:
#       * rep_rate       : share of selected SNPs that are sign-consistent with
#                          EUR and nominally significant in YRI (|t| test at
#                          level --rep_alpha),
#       * mean_abs_r_yri : mean |marginal correlation| of the selection in YRI,
#       * R2_yri / R2_eur: adjusted R^2 of the selected set out-of-population
#                          (YRI) and in-sample (EUR).
#     Covariates (population dummies + sex) are removed from phenotype AND
#     genotypes by FWL residualization before selection; the induced df offset
#     (q <= 4 at n >= 87) is a documented, negligible approximation.
#
#   --mode semisynth — spike-in phenotypes on the real EUR genotypes (real cis
#     LD), truth known. K_true causals drawn proportional to PolyFun snpvar
#     (--truth polyfun, default) or uniformly; per-SNP effect calibrated to
#     marginal power --gamma; B replications; n-gradient by row subsampling
#     (--n_grid). Metrics: exact-SNP M/F1 (evaluate_selection) + LD-cluster
#     credit at |r| > --r_clust.
#
# Heavy comparators (varbvs / CARMA / BayesC) are deliberately left out here.
#
# Usage (from code/):
#   Rscript Benchmark_GEUVADIS.R --rds ../geuvadis_data/geuvadis_ERAP2_pilot.rds \
#       --mode real --K_fixed 5 --adaptive TRUE --outfile ../results_geuvadis_real_ERAP2.rds
#   Rscript Benchmark_GEUVADIS.R --rds_dir ../geuvadis_rds --mode real
#   Rscript Benchmark_GEUVADIS.R --rds ../geuvadis_data/geuvadis_ERAP2_pilot.rds \
#       --mode semisynth --B 200 --n_grid 80,200,358 --gamma 0.6 --K_true 5 \
#       --outfile ../results_geuvadis_semisynth_ERAP2.rds
# =============================================================================

suppressPackageStartupMessages({
  library(parallel); library(dplyr); library(methods)
})
source("Generate_XY_realistic.R")                       # Generate.Y
source("Test_statistics_functions_optimized.R")         # in-house fixed-K + make_quad_nodes
source("Test_statistics_functions_Criteria_optimized.R")# in-house adaptive (*_Criteria_For_Data)
source("External_Comparators.R")                        # run_susie / run_lasso / evaluate_selection

# ---------- 1. Configuration -------------------------------------------------
cfg <- list(
  rds        = "",                 # single per-gene RDS
  rds_dir    = "",                 # or a directory of geuvadis_*.rds (batch)
  mode       = "real",             # "real" or "semisynth"
  # shared
  K_fixed    = 5L,
  tau        = 0.2,
  S_mc       = 100L,               # MC draws for NaM/ExM (raise for production if needed)
  multiple_exact = FALSE,         # TRUE: NaM/ExM prior integrals computed exactly (S_mc unused)
  seed       = 2026L,
  n_cores    = max(1L, parallel::detectCores() - 1L),
  outfile    = "../results_geuvadis.rds",
  skip_susie = FALSE,
  skip_lasso = FALSE,
  inhouse_only = "",               # comma-separated subset, e.g. "AsS,NaM,ExM" (empty = all)
  # real mode
  prior_scale = "primary",         # amplitude of the three-level priors:
                                   # "primary" {0.3,0.5,0.7} (as stored),
                                   # "weak" {0.4,0.5,0.6}, "strong" {0.2,0.5,0.8}
  adaptive   = TRUE,               # also run the adaptive-stopping variants
  criteria   = "BIC",
  K_threshold= 30L,                # iteration cap for the adaptive variants
  target_prob= 0.95,
  rep_alpha  = 0.05,               # nominal level of the YRI replication test
  # semisynth mode
  B          = 200L,
  n_grid     = "80,200,358",
  gamma      = 0.6,
  K_true     = 5L,
  truth      = "polyfun",          # causal sampling: "polyfun" (~ snpvar) or "uniform"
  r_clust    = 0.7
)
args <- commandArgs(trailingOnly = TRUE)
i <- 1L
while (i <= length(args) - 1L) {
  if (startsWith(args[i], "--")) {
    key <- sub("^--", "", args[i])
    if (!is.null(cfg[[key]])) cfg[[key]] <- methods::as(args[i + 1L], class(cfg[[key]]))
    i <- i + 2L
  } else i <- i + 1L
}
stopifnot(cfg$mode %in% c("real", "semisynth"))
cat("Configuration:\n"); str(cfg)

rds_files <- if (nzchar(cfg$rds_dir)) {
  list.files(cfg$rds_dir, pattern = "^geuvadis_.*\\.rds$", full.names = TRUE)
} else cfg$rds
if (!length(rds_files)) stop("Provide --rds or --rds_dir.")

in_house <- c(NaW = "Naive_Wald",       AsW = "Asymptotic_Wald",
              ExW = "Exact_Wald",       AsS = "Asymptotic_Score",
              NaM = "Naive_Multiple",   ExM = "Exact_Multiple")
if (nzchar(cfg$inhouse_only)) {
  sel <- strsplit(cfg$inhouse_only, ",", fixed = TRUE)[[1L]]
  in_house <- in_house[names(in_house) %in% sel | in_house %in% sel]
}

fixedK_call <- function(fn, X, Y, pj, K, gh) {
  switch(fn,
    Naive_Wald       = Naive_Wald(X, Y, pj, K, cfg$tau),
    Asymptotic_Wald  = Asymptotic_Wald(X, Y, pj, K, cfg$tau, gh),
    Exact_Wald       = Exact_Wald(X, Y, pj, K, cfg$tau, gh),
    Asymptotic_Score = Asymptotic_Score(X, Y, pj, K, cfg$tau, gh),
    Naive_Multiple   = Naive_Multiple(X, Y, pj, K, cfg$tau^2, (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)),
    Exact_Multiple   = Exact_Multiple(X, Y, pj, K, cfg$tau^2, (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc)))
}
adaptive_call <- function(fn, X, Y, pj, gh) {
  switch(fn,
    Naive_Wald       = Naive_Wald_Criteria_For_Data(X, Y, pj, cfg$tau,
                         cfg$K_threshold, cfg$target_prob, cfg$criteria),
    Asymptotic_Wald  = Asymptotic_Wald_Criteria_For_Data(X, Y, pj, cfg$tau, gh,
                         cfg$K_threshold, cfg$target_prob, cfg$criteria),
    Exact_Wald       = Exact_Wald_Criteria_For_Data(X, Y, pj, cfg$tau, gh,
                         cfg$K_threshold, cfg$target_prob, cfg$criteria),
    Asymptotic_Score = Asymptotic_Score_Criteria_For_Data(X, Y, pj, cfg$tau, gh,
                         cfg$K_threshold, cfg$target_prob, cfg$criteria),
    Naive_Multiple   = Naive_Multiple_Criteria_For_Data(X, Y, pj, cfg$tau^2, (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc),
                         cfg$K_threshold, cfg$target_prob, cfg$criteria),
    Exact_Multiple   = Exact_Multiple_Criteria_For_Data(X, Y, pj, cfg$tau^2, (if (isTRUE(cfg$multiple_exact)) Inf else cfg$S_mc),
                         cfg$K_threshold, cfg$target_prob, cfg$criteria))
}

# FWL residualization: remove [1, C] from the columns of M
residualize <- function(M, C) {
  Q <- qr.Q(qr(cbind(1, as.matrix(C))))
  M <- as.matrix(M)
  M - Q %*% crossprod(Q, M)
}

# Exact residual-dimension reduction. The FWL-residualized data live in the
# (n - q)-dimensional orthogonal complement of [1, C] (q = rank of [1, C]).
# This maps them onto n - q + 1 rows through V W'M, where the columns of W are
# an orthonormal basis of that complement and the columns of V an orthonormal
# basis of the centered subspace of R^(n-q+1) (V'V = I, 1'V = 0), so inner
# products and column norms are preserved and the image is exactly centered.
# Every downstream method then computes its residual degrees of freedom from
# the row dimension: n_rows - k - 1 = n - q - k, the exact residual dimension
# after projecting out the covariates, with no approximation.
reduce_residual <- function(M, C) {
  M <- as.matrix(M); n <- nrow(M)
  qrC <- qr(cbind(1, as.matrix(C)))
  q <- qrC$rank
  W <- qr.Q(qrC, complete = TRUE)[, (q + 1L):n, drop = FALSE]
  np <- n - q + 1L
  V <- qr.Q(qr(matrix(1, np, 1L)), complete = TRUE)[, -1L, drop = FALSE]
  V %*% crossprod(W, M)
}

# =============================================================================
# 2. REAL mode: EUR discovery -> YRI replication
# =============================================================================
run_real_gene <- function(gx) {
  YX <- reduce_residual(cbind(as.numeric(gx$y_eur), as.matrix(gx$G_eur)), gx$covar_eur)
  y  <- YX[, 1L]
  X  <- YX[, -1L, drop = FALSE]
  keep <- apply(X, 2, sd) > 1e-8
  X  <- scale(X[, keep, drop = FALSE]); X[is.nan(X)] <- 0
  YXv <- reduce_residual(cbind(as.numeric(gx$y_yri), as.matrix(gx$G_yri)[, keep, drop = FALSE]),
                         gx$covar_yri)
  yv <- YXv[, 1L]
  Xv <- scale(YXv[, -1L, drop = FALSE]); Xv[is.nan(Xv)] <- 0
  pos    <- gx$pos[keep]
  priors <- lapply(gx$priors, function(p) p[keep])
  # Prior-amplitude sensitivity: remap the stored three-level scale
  # {0.3, 0.5, 0.7} onto a weaker or stronger contrast (level by level; the
  # uniform arm, constant at 0.5, is invariant).
  if (cfg$prior_scale != "primary") {
    lv <- switch(cfg$prior_scale,
                 weak   = c(`0.3` = 0.4, `0.5` = 0.5, `0.7` = 0.6),
                 strong = c(`0.3` = 0.2, `0.5` = 0.5, `0.7` = 0.8),
                 stop("prior_scale must be primary, weak or strong"))
    priors <- lapply(priors, function(p) unname(lv[sprintf("%.1f", p)]))
    stopifnot(!anyNA(unlist(priors)))
  }
  n <- nrow(X); n_v <- length(yv)   # exact residual dimensions (+1 for centering)
  gh <- make_quad_nodes(30L, cfg$tau)
  cat(sprintf("[%s] real mode: m=%d, n_EUR_res=%d, n_YRI_res=%d\n", gx$gene, ncol(X), n, n_v))

  r_eur <- as.numeric(cor(y, X))
  r_yri <- suppressWarnings(as.numeric(cor(yv, Xv))); r_yri[is.na(r_yri)] <- 0
  t_crit <- qt(1 - cfg$rep_alpha / 2, df = n_v - 2)
  # LD-aware replication credit: a selected SNP also counts as replicated when
  # it is in YRI genotype LD (|r| >= r_clust) with a SNP that itself replicates
  # (sign-consistent, nominally significant in YRI). YRI LD is the right panel:
  # an EUR-specific proxy that decorrelates from the causal in YRI gets no credit.
  Gv_raw  <- gx$G_yri[, keep, drop = FALSE]
  t_all   <- r_yri * sqrt((n_v - 2) / pmax(1 - r_yri^2, 1e-12))
  rep_set <- which((sign(r_yri) == sign(r_eur)) & (abs(t_all) > t_crit))

  eval_sel <- function(idx) {
    idx <- as.integer(idx); K <- length(idx)
    if (!K) return(data.frame(K = 0L, rep_rate = NA_real_, rep_rate_ld = NA_real_,
                              mean_abs_r_yri = NA_real_,
                              R2_yri = NA_real_, R2_eur = NA_real_))
    r_y <- r_yri[idx]
    t_y <- r_y * sqrt((n_v - 2) / pmax(1 - r_y^2, 1e-12))
    repl <- (sign(r_y) == sign(r_eur[idx])) & (abs(t_y) > t_crit)
    if (length(rep_set)) {
      ldC <- suppressWarnings(abs(cor(Gv_raw[, idx, drop = FALSE],
                                      Gv_raw[, rep_set, drop = FALSE])))
      ldC[is.na(ldC)] <- 0
      rep_ld <- rowSums(ldC >= cfg$r_clust) > 0
    } else rep_ld <- rep(FALSE, K)
    R2v <- tryCatch(summary(lm(yv ~ Xv[, idx, drop = FALSE]))$adj.r.squared,
                    error = function(e) NA_real_)
    R2e <- tryCatch(summary(lm(y ~ X[, idx, drop = FALSE]))$adj.r.squared,
                    error = function(e) NA_real_)
    data.frame(K = K, rep_rate = mean(repl), rep_rate_ld = mean(rep_ld),
               mean_abs_r_yri = mean(abs(r_y)),
               R2_yri = R2v, R2_eur = R2e)
  }

  rows <- list(); sels <- list()
  rec <- function(prior, method, selection, idx, rt) {
    rows[[length(rows) + 1L]] <<- cbind(
      data.frame(gene = gx$gene, prior = prior, method = method,
                 selection = selection, stringsAsFactors = FALSE),
      eval_sel(idx), data.frame(runtime = round(rt, 2)))
    sels[[paste(prior, method, selection, sep = ".")]] <<- pos[as.integer(idx)]
  }

  for (pr in names(priors)) {
    pj <- priors[[pr]]
    for (mn in names(in_house)) {
      t0 <- Sys.time()
      out <- tryCatch(fixedK_call(in_house[[mn]], X, y, pj, cfg$K_fixed, gh),
                      error = function(e) { warning(mn, ": ", conditionMessage(e)); NULL })
      if (!is.null(out)) rec(pr, mn, "fixedK", out$indices,
                             as.numeric(difftime(Sys.time(), t0, units = "secs")))
      if (isTRUE(cfg$adaptive)) {
        t0 <- Sys.time()
        out <- tryCatch(adaptive_call(in_house[[mn]], X, y, pj, gh),
                        error = function(e) { warning(mn, ": ", conditionMessage(e)); NULL })
        if (!is.null(out)) rec(pr, mn, cfg$criteria, out$indices,
                               as.numeric(difftime(Sys.time(), t0, units = "secs")))
      }
    }
    if (!isTRUE(cfg$skip_susie)) {
      s <- tryCatch(run_susie(X, y, K = cfg$K_fixed,
                              prior_p = if (pr == "uniform") NULL else pj),
                    error = function(e) NULL)
      if (!is.null(s)) rec(pr, "SuSiE", "fixedK", s$indices, s$runtime)
    }
  }
  if (!isTRUE(cfg$skip_lasso)) {
    l <- tryCatch(run_lasso(X, y, K = cfg$K_fixed, lambda_rule = "1se"),
                  error = function(e) NULL)
    if (!is.null(l)) rec("none", "LASSO", "fixedK", l$indices, l$runtime)
  }
  list(table = do.call(rbind, rows), selections = sels,
       marginals = list(pos = pos, r_eur = r_eur, r_yri = r_yri))
}

# =============================================================================
# 3. SEMISYNTH mode: spike-in on real EUR genotypes, known truth, n-gradient
# =============================================================================
run_semisynth_gene <- function(gx) {
  n_values <- as.integer(strsplit(cfg$n_grid, ",", fixed = TRUE)[[1L]])
  n_values <- pmin(n_values, nrow(gx$G_eur))
  sv <- gx$snpvar
  sv[is.na(sv) | sv <= 0] <- min(sv[sv > 0], na.rm = TRUE) * 1e-3
  grid <- expand.grid(b = seq_len(cfg$B), n = n_values)
  cat(sprintf("[%s] semisynth: m=%d, n in {%s}, B=%d -> %d cells on %d cores\n",
              gx$gene, ncol(gx$G_eur), paste(n_values, collapse = ","),
              cfg$B, nrow(grid), cfg$n_cores))

  one_cell <- function(g) {
    b <- grid$b[g]; n <- grid$n[g]
    set.seed(cfg$seed + 7919L * g)
    rows_sub <- if (n < nrow(gx$G_eur)) sample(nrow(gx$G_eur), n) else seq_len(nrow(gx$G_eur))
    Xraw <- gx$G_eur[rows_sub, , drop = FALSE]
    ok <- apply(Xraw, 2, sd) > 1e-6
    X <- scale(Xraw[, ok, drop = FALSE]); X[is.nan(X)] <- 0
    mk <- ncol(X)
    truth <- sort(sample(mk, cfg$K_true,
                         prob = if (cfg$truth == "polyfun") sv[ok] else NULL))
    lambda <- (qnorm(0.975) + qnorm(cfg$gamma)) / sqrt(n)
    beta <- numeric(mk); beta[truth] <- lambda
    Y <- Generate.Y(n = n, X = X, beta = beta, sigma = 1)
    priors_ok <- lapply(gx$priors, function(p) p[ok])
    gh <- make_quad_nodes(30L, cfg$tau)

    rec <- function(method, prior, idx, rt) {
      idx <- as.integer(idx)
      e <- evaluate_selection(idx, truth, mk)
      if (length(idx)) {          # LD-cluster credit on the subsample correlations
        Ct <- suppressWarnings(abs(cor(X[, truth, drop = FALSE], X[, idx, drop = FALSE])))
        Ct[is.na(Ct)] <- 0
        Mc    <- sum(apply(Ct, 1, max) > cfg$r_clust)
        precC <- mean(apply(Ct, 2, max) > cfg$r_clust)
        recC  <- Mc / length(truth)
        F1c   <- if (precC + recC > 0) 2 * precC * recC / (precC + recC) else 0
      } else { Mc <- 0L; F1c <- 0 }
      data.frame(gene = gx$gene, n = n, rep = b, method = method, prior = prior,
                 M = e$M, K = e$K, precision = e$precision, recall = e$recall,
                 F1 = e$F1, FDR = e$FDR, M_clust = Mc, F1_clust = F1c,
                 runtime = round(rt, 2), stringsAsFactors = FALSE)
    }

    out <- list()
    for (pr in names(priors_ok)) {
      pj <- priors_ok[[pr]]
      for (mn in names(in_house)) {
        t0 <- Sys.time()
        r <- tryCatch(fixedK_call(in_house[[mn]], X, Y, pj, cfg$K_true, gh),
                      error = function(e) NULL)
        if (!is.null(r)) out[[length(out) + 1L]] <-
          rec(mn, pr, r$indices, as.numeric(difftime(Sys.time(), t0, units = "secs")))
      }
      if (!isTRUE(cfg$skip_susie)) {
        s <- tryCatch(run_susie(X, Y, K = cfg$K_true,
                                prior_p = if (pr == "uniform") NULL else pj),
                      error = function(e) NULL)
        if (!is.null(s)) out[[length(out) + 1L]] <- rec("SuSiE", pr, s$indices, s$runtime)
      }
    }
    if (!isTRUE(cfg$skip_lasso)) {
      l <- tryCatch(run_lasso(X, Y, K = cfg$K_true, lambda_rule = "1se"),
                    error = function(e) NULL)
      if (!is.null(l)) out[[length(out) + 1L]] <- rec("LASSO", "none", l$indices, l$runtime)
    }
    do.call(rbind, out)
  }

  res <- mclapply(seq_len(nrow(grid)),
                  function(g) tryCatch(one_cell(g), error = function(e) NULL),
                  mc.cores = cfg$n_cores)
  do.call(rbind, Filter(Negate(is.null), res))
}

# =============================================================================
# 4. Run over genes + aggregate + save
# =============================================================================
t0 <- Sys.time()
if (cfg$mode == "real") {
  per_gene <- mclapply(rds_files, function(f)
    tryCatch(run_real_gene(readRDS(f)),
             error = function(e) { message("SKIP ", f, ": ", conditionMessage(e)); NULL }),
    mc.cores = max(1L, min(cfg$n_cores, length(rds_files))))
  per_gene <- Filter(Negate(is.null), per_gene)
  if (!length(per_gene)) stop("real mode: no gene produced results.")
  tab <- do.call(rbind, lapply(per_gene, `[[`, "table"))
  cat("\n=========== GEUVADIS real cis-eQTL: EUR discovery -> YRI replication ===========\n")
  print(tab %>% arrange(gene, prior, selection, desc(rep_rate)), row.names = FALSE)
  saveRDS(list(config = cfg, mode = "real", table = tab,
               selections = lapply(per_gene, `[[`, "selections"),
               marginals  = lapply(per_gene, `[[`, "marginals"),
               genes = vapply(per_gene, function(z) z$table$gene[1], "")),
          cfg$outfile)
} else {
  raw <- do.call(rbind, lapply(rds_files, function(f) run_semisynth_gene(readRDS(f))))
  summary_df <- raw %>%
    group_by(gene, n, prior, method) %>%
    summarise(n_reps = dplyr::n(),
              M = round(mean(M, na.rm = TRUE), 2), F1 = round(mean(F1, na.rm = TRUE), 3),
              M_clust = round(mean(M_clust, na.rm = TRUE), 2),
              F1_clust = round(mean(F1_clust, na.rm = TRUE), 3),
              precision = round(mean(precision, na.rm = TRUE), 3),
              recall = round(mean(recall, na.rm = TRUE), 3),
              runtime_s = round(mean(runtime, na.rm = TRUE), 2), .groups = "drop") %>%
    arrange(gene, n, prior, desc(F1))
  cat("\n=========== GEUVADIS semi-synthetic (real cis LD, known truth) ===========\n")
  print(as.data.frame(summary_df), row.names = FALSE)
  saveRDS(list(config = cfg, mode = "semisynth", summary = summary_df, raw = raw),
          cfg$outfile)
}
cat(sprintf("\nDone in %.1f s. Results written to: %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")), cfg$outfile))
