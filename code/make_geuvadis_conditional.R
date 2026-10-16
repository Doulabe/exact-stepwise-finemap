# =============================================================================
# File: make_geuvadis_conditional.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Conditional (path-wise) replication sensitivity for the GEUVADIS
# application. The stepwise selections are conditional by construction, so
# this analysis replicates them conditionally as well: for each ordered
# selection (fixed K), variant k is tested in the YRI panel GIVEN the k - 1
# variants selected before it (t statistic of the k-th coefficient in the
# incremental regression on residualized data, exact residual degrees of
# freedom n_YRI - rank([1, C]) - k), and counts as conditionally replicating
# when this effect is sign-consistent with the corresponding EUR conditional
# effect and significant at alpha. Collinear or undefined coefficients count
# as non-replicating.
#
# Usage (from code/):
#   Rscript make_geuvadis_conditional.R \
#     --results ../results_geuvadis_real_ld.rds --rds_dir ../geuvadis_rds \
#     --out ../tab_geuvadis_cond.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            rds_dir = "../geuvadis_rds", alpha = 0.05,
            out = "../tab_geuvadis_cond.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

residualize <- function(M, C) {
  q <- qr(cbind(1, as.matrix(C)))
  qr.resid(q, as.matrix(M))
}

res <- readRDS(cfg$results)
stopifnot(res$mode == "real")
genes <- res$genes

# Conditional replication of one ordered selection: share of its variants
# whose YRI effect, given the variants selected before them, is
# sign-consistent with the EUR conditional effect and significant at alpha.
cond_rate <- function(idx, y, X, yv, Xv, q_eur, q_yri) {
  K <- length(idx)
  if (!K) return(NA_real_)
  n  <- nrow(X); n_v <- nrow(Xv)
  ok <- logical(K)
  for (k in seq_len(K)) {
    sub  <- idx[seq_len(k)]
    fe <- tryCatch(lm.fit(cbind(1, X[, sub, drop = FALSE]), y),  error = function(e) NULL)
    fv <- tryCatch(lm.fit(cbind(1, Xv[, sub, drop = FALSE]), yv), error = function(e) NULL)
    if (is.null(fe) || is.null(fv)) next
    be <- fe$coefficients[k + 1L]; bv <- fv$coefficients[k + 1L]
    if (is.na(be) || is.na(bv)) next
    df_v <- n_v - q_yri - k                 # exact residual df in YRI
    R    <- qr.R(fv$qr)
    XtXi <- tryCatch(chol2inv(R), error = function(e) NULL)
    if (is.null(XtXi) || df_v <= 0) next
    s2 <- sum(fv$residuals^2) / df_v
    se <- sqrt(pmax(s2 * XtXi[k + 1L, k + 1L], 1e-300))
    tv <- bv / se
    ok[k] <- (sign(bv) == sign(be)) && (abs(tv) > qt(1 - cfg$alpha / 2, df_v))
  }
  mean(ok)
}

rows <- list()
for (i in seq_along(genes)) {
  g  <- genes[i]
  mg <- res$marginals[[i]]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
  colsel <- match(mg$pos, gx$pos)
  y  <- as.numeric(residualize(gx$y_eur, gx$covar_eur))
  X  <- residualize(gx$G_eur, gx$covar_eur)[, colsel, drop = FALSE]
  yv <- as.numeric(residualize(gx$y_yri, gx$covar_yri))
  Xv <- residualize(gx$G_yri, gx$covar_yri)[, colsel, drop = FALSE]
  q_eur <- qr(cbind(1, as.matrix(gx$covar_eur)))$rank
  q_yri <- qr(cbind(1, as.matrix(gx$covar_yri)))$rank

  for (key in names(res$selections[[i]])) {
    parts <- strsplit(key, ".", fixed = TRUE)[[1L]]
    if (parts[3] != "fixedK") next
    idx <- match(res$selections[[i]][[key]], mg$pos)
    idx <- idx[!is.na(idx)]
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, prior = parts[1], method = parts[2],
      cond_rate = cond_rate(idx, y, X, yv, Xv, q_eur, q_yri),
      stringsAsFactors = FALSE)
  }
  if (i %% 20 == 0) cat(i, "genes\n")
}
d <- bind_rows(rows)

METHODS <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM", "SuSiE")
PRIORS  <- c(uniform = "Uniform", polyfun = "PolyFun", tss = "TSS",
             discordant = "Discordant")
w <- d %>% filter(method %in% METHODS, prior %in% names(PRIORS)) %>%
  group_by(method, prior) %>%
  summarise(v = mean(cond_rate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = prior, values_from = v) %>%
  slice(match(METHODS, method))
lasso <- d %>% filter(method == "LASSO") %>%
  summarise(v = mean(cond_rate, na.rm = TRUE)) %>% pull(v)

fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
trows <- vapply(seq_len(nrow(w)), function(i)
  paste0(w$method[i], " & ",
         paste(vapply(names(PRIORS), function(p) fmt(w[[p]][i]), ""), collapse = " & "),
         " \\\\"), "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Conditional replication sensitivity (fixed $K=5$).} ",
         "The stepwise selections are conditional by construction, so this variant replicates ",
         "them conditionally as well. Variant $k$ of each ordered selection is tested in YRI ",
         "given the $k-1$ variants selected before it (incremental regression on the ",
         "residualized data, exact residual degrees of freedom), and counts as replicating ",
         "when this conditional effect is sign-consistent with the corresponding EUR ",
         "conditional effect and significant at $\\alpha=0.05$. Entries: mean share of ",
         "conditionally replicating variants over the $100$ genes. The prior-free LASSO ",
         "(ordered by $|$coefficient$|$) attains $", sprintf("%.3f", lasso), "$.}}"),
  "\\label{tab:geuvadis_cond}",
  "\\centering\\small",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
  trows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
print(as.data.frame(w), row.names = FALSE)
cat("LASSO:", lasso, "\n")
