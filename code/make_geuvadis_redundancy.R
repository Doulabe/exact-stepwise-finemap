# =============================================================================
# File: make_geuvadis_redundancy.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Redundancy diagnostics of the fixed-K GEUVADIS selections, addressing the
# question of whether the LD-aware replication metric rewards selections of
# redundant proxies of a single signal. For each selection: (i) the mean
# pairwise |r| between the selected variants (EUR raw dosages); (ii) the
# number of distinct LD clusters covered (greedy clustering at |r| >= 0.7,
# the credit threshold); (iii) the share of those clusters containing at
# least one variant with LD-aware replication credit in YRI (cluster-level
# replication, which counts each tagged signal once however many proxies
# are selected).
#
# Usage (from code/):
#   Rscript make_geuvadis_redundancy.R \
#     --results ../results_geuvadis_real_ld.rds --rds_dir ../geuvadis_rds \
#     --out ../tab_geuvadis_clusters.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(results = "../results_geuvadis_real_ld.rds",
            rds_dir = "../geuvadis_rds", r_clust = 0.7, rep_alpha = 0.05,
            out = "../tab_geuvadis_clusters.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}

res <- readRDS(cfg$results)
genes <- res$genes

greedy_clusters <- function(C, thr) {
  # C: |r| matrix between selected variants; returns cluster id per variant
  K <- nrow(C); id <- integer(K)
  for (i in seq_len(K)) {
    if (id[i] > 0L) next
    id[i] <- max(id) + 1L
    for (j in seq_len(K)) if (id[j] == 0L && C[i, j] >= thr) id[j] <- id[i]
  }
  id
}

rows <- list()
for (i in seq_along(genes)) {
  g  <- genes[i]
  mg <- res$marginals[[i]]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
  colsel <- match(mg$pos, gx$pos)
  Ge <- gx$G_eur[, colsel, drop = FALSE]     # raw dosages, discovery panel
  Gv <- gx$G_yri[, colsel, drop = FALSE]
  n_v <- nrow(Gv)
  df_v <- n_v - qr(cbind(1, as.matrix(gx$covar_yri)))$rank - 1L
  t_crit <- qt(1 - cfg$rep_alpha / 2, df = df_v)
  t_all  <- mg$r_yri * sqrt(df_v / pmax(1 - mg$r_yri^2, 1e-12))
  rep_set <- which((sign(mg$r_yri) == sign(mg$r_eur)) & (abs(t_all) > t_crit))

  for (key in names(res$selections[[i]])) {
    parts <- strsplit(key, ".", fixed = TRUE)[[1L]]
    if (parts[3] != "fixedK") next
    idx <- match(res$selections[[i]][[key]], mg$pos)
    idx <- idx[!is.na(idx)]
    K <- length(idx)
    if (K < 2L) next
    Ce <- suppressWarnings(abs(cor(Ge[, idx, drop = FALSE]))); Ce[is.na(Ce)] <- 0
    mean_r <- mean(Ce[upper.tri(Ce)])
    cl <- greedy_clusters(Ce, cfg$r_clust)
    # LD-aware credit per selected variant (YRI panel, as in the main metric)
    if (length(rep_set)) {
      ldC <- suppressWarnings(abs(cor(Gv[, idx, drop = FALSE],
                                      Gv[, rep_set, drop = FALSE])))
      ldC[is.na(ldC)] <- 0
      cred <- rowSums(ldC >= cfg$r_clust) > 0
    } else cred <- rep(FALSE, K)
    cl_rep <- mean(vapply(unique(cl), function(cc) any(cred[cl == cc]), TRUE))
    rows[[length(rows) + 1L]] <- data.frame(
      gene = g, prior = parts[1], method = parts[2],
      mean_r = mean_r, n_clusters = length(unique(cl)), cluster_rep = cl_rep)
  }
  if (i %% 25 == 0) cat(i, "genes\n")
}
d <- bind_rows(rows)

METHODS <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM", "SuSiE")
PRIORS  <- c(uniform = "Uniform", polyfun = "PolyFun", tss = "TSS",
             discordant = "Discordant")
w <- d %>% filter(method %in% METHODS, prior %in% names(PRIORS)) %>%
  group_by(method, prior) %>%
  summarise(r = mean(mean_r), C = mean(n_clusters), cr = mean(cluster_rep),
            .groups = "drop")
lass <- d %>% filter(method == "LASSO") %>%
  summarise(r = mean(mean_r), C = mean(n_clusters), cr = mean(cluster_rep))

trows <- vapply(METHODS, function(m) {
  cells <- vapply(names(PRIORS), function(p) {
    z <- w %>% filter(method == m, prior == p)
    if (!nrow(z)) return("--")
    sprintf("%.1f / %.2f / %.2f", z$C, z$r, z$cr)
  }, "")
  paste0(m, " & ", paste(cells, collapse = " & "), " \\\\")
}, "")
writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Redundancy diagnostics of the fixed-$K=5$ selections.} ",
         "Each cell reports, over the $100$ genes, the mean number of distinct LD clusters ",
         "covered by the five selected variants (greedy clustering of the selection at ",
         "$|r| \\geq 0.7$ on EUR raw dosages, maximum $5$), the mean pairwise $|r|$ within the ",
         "selection, and the mean share of the covered clusters containing at least one variant ",
         "with LD-aware replication credit, which counts each tagged signal once however many ",
         "proxies are selected. The prior-free LASSO ($\\hat K = 3.8$) attains ",
         sprintf("$%.1f$ clusters, within-selection $|r|$ of $%.2f$ and a cluster-level rate of $%.2f$.", lass$C, lass$r, lass$cr),
         "}}"),
  "\\label{tab:geuvadis_clusters}",
  "\\centering\\scriptsize",
  "\\setlength{\\tabcolsep}{3pt}",
  "\\begin{tabular}{lcccc}", "\\toprule",
  paste0("Method & ", paste(PRIORS, collapse = " & "), " \\\\"),
  "\\midrule", trows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
print(as.data.frame(w %>% mutate(across(where(is.numeric), ~round(., 3)))), row.names = FALSE)
cat("LASSO:", round(lass$C, 2), round(lass$r, 3), round(lass$cr, 3), "\n")
