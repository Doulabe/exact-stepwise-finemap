# =============================================================================
# File: merge_fixed_methods.R
# Project: exact-stepwise-finemap
# Author: N.K. Doulabe
#
# Replaces the rows of a subset of in-house methods (typically AsS, NaM, ExM
# after the 2026-09-28 correction of the score statistic and of the Multiple
# non-centrality) in an existing result file by the rows of a rerun restricted
# to those methods (--inhouse_only), then recomputes the summary with the same
# aggregation as the benchmark script that produced the file. Rows of all other
# methods are kept unchanged.
#
# Usage (from code/):
#   Rscript merge_fixed_methods.R --type external --orig ../results_orig/X.rds \
#       --new ../results_rerun/X.rds --out ../results_fixed/X.rds
#   --type : external (Benchmark_External.R), locus (Benchmark_Locus_Sim.R),
#            geuvadis_real, geuvadis_semisynth (Benchmark_GEUVADIS.R)
#   --keep_runtime TRUE keeps the original runtimes of the replaced rows
#            (used when the rerun hardware differs from the original one).
#   --note TEXT : description of the correction, stored with the merged rows.
#   --canary M : a method rerun alongside but NOT replaced (e.g. Naive_Wald or
#            NaW, deterministic); its rows are compared with the original ones
#            to check that the rerun regenerated exactly the same datasets.
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

cfg <- list(type = "", orig = "", new = "", out = "", keep_runtime = FALSE, canary = "", note = "")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}
stopifnot(cfg$type %in% c("external", "locus", "geuvadis_real", "geuvadis_semisynth"),
          file.exists(cfg$orig), file.exists(cfg$new), nzchar(cfg$out))

orig <- readRDS(cfg$orig)
new  <- readRDS(cfg$new)
canary_status <- "NA"

report <- function(o, n, m, key, what) {
  for (mm in m) {
    f1o <- suppressWarnings(mean(as.numeric(o$F1[o$method == mm]), na.rm = TRUE))
    f1n <- suppressWarnings(mean(as.numeric(n$F1[n$method == mm]), na.rm = TRUE))
    cat(sprintf("  %-18s rows %5d -> %5d   mean F1 %s -> %s\n", mm,
                sum(o$method == mm), sum(n$method == mm),
                if (is.finite(f1o)) sprintf("%.3f", f1o) else "NA",
                if (is.finite(f1n)) sprintf("%.3f", f1n) else "NA"))
  }
}

# Replace the rows of the rerun methods; optionally carry over the old runtimes.
replace_rows <- function(o, n, key_cols) {
  if (nzchar(cfg$canary)) {
    ko <- do.call(paste, c(o[o$method == cfg$canary, key_cols, drop = FALSE], sep = "|"))
    kn <- do.call(paste, c(n[n$method == cfg$canary, key_cols, drop = FALSE], sep = "|"))
    vals <- intersect(c("M", "K", "F1", "rep_rate", "M_clust"), names(o))
    vo <- do.call(paste, c(o[o$method == cfg$canary, vals, drop = FALSE], sep = "|"))
    vn <- do.call(paste, c(n[n$method == cfg$canary, vals, drop = FALSE], sep = "|"))
    same <- sum(vo[match(kn, ko)] == vn, na.rm = TRUE)
    # With real genotypes subsampled to small n, several SNPs can be identical
    # columns; which of them is selected then depends on the last bit of the
    # floating-point arithmetic, which differs between machines. Such rows keep
    # the same number of selected SNPs and the same LD-cluster metrics.
    ties <- 0L
    if (same < length(kn) && all(c("M_clust", "F1_clust") %in% names(o))) {
      cl <- c("K", "M_clust", "F1_clust")
      co <- do.call(paste, c(o[o$method == cfg$canary, cl, drop = FALSE], sep = "|"))
      cn <- do.call(paste, c(n[n$method == cfg$canary, cl, drop = FALSE], sep = "|"))
      ties <- sum(co[match(kn, ko)] == cn, na.rm = TRUE) - same
    }
    status <- if (!length(kn)) "MISMATCH" else if (same == length(kn)) "OK" else
              if (same + ties == length(kn)) "TIES" else "MISMATCH"
    cat(sprintf("  canary %s: %d / %d rerun rows identical to the original%s\n",
                cfg$canary, same, length(kn),
                switch(status, OK = "",
                       TIES = sprintf("; the %d others differ only in which of several identical genotype columns was selected (LD-cluster metrics identical)", ties),
                       MISMATCH = "   CANARY MISMATCH")))
    canary_status <<- status
    n <- n[n$method != cfg$canary, , drop = FALSE]
  }
  meths <- unique(n$method)
  n <- n[, intersect(names(o), names(n)), drop = FALSE]
  for (cc in setdiff(names(o), names(n))) n[[cc]] <- NA
  n <- n[, names(o), drop = FALSE]
  if (isTRUE(cfg$keep_runtime) && "runtime" %in% names(o)) {
    ko <- do.call(paste, c(o[o$method %in% meths, key_cols, drop = FALSE], sep = "|"))
    kn <- do.call(paste, c(n[, key_cols, drop = FALSE], sep = "|"))
    rt <- o$runtime[o$method %in% meths][match(kn, ko)]
    n$runtime[!is.na(rt)] <- rt[!is.na(rt)]
  }
  ok <- do.call(paste, c(o[o$method %in% meths, key_cols, drop = FALSE], sep = "|"))
  nk <- do.call(paste, c(n[, key_cols, drop = FALSE], sep = "|"))
  miss <- setdiff(ok, nk)
  if (length(miss)) warning(sprintf("%d original key(s) have no rerun row (e.g. %s)",
                                    length(miss), miss[1L]))
  out <- rbind(o[!(o$method %in% meths), , drop = FALSE], n)
  rownames(out) <- NULL
  list(rows = out, methods = meths)
}

stamp <- function(obj, meths) {
  obj$fix_20260928 <- list(
    note = if (nzchar(cfg$note)) paste("Rows of", paste(meths, collapse = ", "), "recomputed:", cfg$note)
           else paste("Rows of", paste(meths, collapse = ", "),
                      "recomputed after the 2026-09-28 code fix (AsS score statistic and",
                      "law; NaM/ExM non-centrality in R's convention)."),
    rerun_file = normalizePath(cfg$new), keep_runtime = cfg$keep_runtime,
    canary = cfg$canary, canary_status = canary_status)
  obj
}

if (cfg$type == "external") {
  design <- if (!is.null(orig$design)) orig$design else "synthetic"
  key <- c("b", "n", "rho_intra", "gamma", "tau", "method", "criterion")
  r <- replace_rows(orig$all_results, new$all_results, key)
  cat("external:", basename(cfg$orig), "\n"); report(orig$all_results, r$rows, r$methods)
  group_cols <- if (design == "real") c("n", "gamma", "tau", "method", "criterion")
                else c("n", "rho_intra", "gamma", "tau", "method", "criterion")
  orig$all_results <- r$rows
  orig$summary <- r$rows %>%
    mutate(M = as.numeric(M), K = as.numeric(K),
           precision = as.numeric(precision), recall = as.numeric(recall),
           F1 = as.numeric(F1), FDR = as.numeric(FDR),
           runtime = as.numeric(runtime), n = as.integer(n),
           rho_intra = as.numeric(rho_intra),
           gamma = as.numeric(gamma), tau = as.numeric(tau)) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(n_reps = dplyr::n(),
              mean_M = mean(M, na.rm = TRUE), mean_K = mean(K, na.rm = TRUE),
              precision = mean(precision, na.rm = TRUE),
              recall = mean(recall, na.rm = TRUE),
              F1 = mean(F1, na.rm = TRUE), FDR = mean(FDR, na.rm = TRUE),
              runtime_s = mean(runtime, na.rm = TRUE), .groups = "drop") %>%
    arrange(across(all_of(group_cols)), desc(F1))
  orig <- stamp(orig, r$methods)

} else if (cfg$type == "locus") {
  r <- replace_rows(orig$raw, new$raw, c("b", "n", "method", "criterion"))
  cat("locus:", basename(cfg$orig), "\n"); report(orig$raw, r$rows, r$methods)
  orig$raw <- r$rows
  orig$summary <- r$rows %>%
    mutate(M = as.numeric(M), K = as.numeric(K),
           precision = as.numeric(precision), recall = as.numeric(recall),
           F1 = as.numeric(F1), FDR = as.numeric(FDR),
           runtime = as.numeric(runtime), n = as.integer(n)) %>%
    group_by(n, method, criterion) %>%
    summarise(n_reps = dplyr::n(),
              mean_M = mean(M, na.rm = TRUE), mean_K = mean(K, na.rm = TRUE),
              precision = mean(precision, na.rm = TRUE),
              recall = mean(recall, na.rm = TRUE),
              F1 = mean(F1, na.rm = TRUE), FDR = mean(FDR, na.rm = TRUE),
              runtime_s = mean(runtime, na.rm = TRUE), .groups = "drop") %>%
    arrange(n, method, criterion)
  orig <- stamp(orig, r$methods)

} else if (cfg$type == "geuvadis_semisynth") {
  r <- replace_rows(orig$raw, new$raw, c("gene", "n", "rep", "method", "prior"))
  cat("semisynth:", basename(cfg$orig), "\n"); report(orig$raw, r$rows, r$methods)
  orig$raw <- r$rows
  orig$summary <- r$rows %>%
    group_by(gene, n, prior, method) %>%
    summarise(n_reps = dplyr::n(),
              M = round(mean(M, na.rm = TRUE), 2), F1 = round(mean(F1, na.rm = TRUE), 3),
              M_clust = round(mean(M_clust, na.rm = TRUE), 2),
              F1_clust = round(mean(F1_clust, na.rm = TRUE), 3),
              precision = round(mean(precision, na.rm = TRUE), 3),
              recall = round(mean(recall, na.rm = TRUE), 3),
              runtime_s = round(mean(runtime, na.rm = TRUE), 2), .groups = "drop") %>%
    arrange(gene, n, prior, desc(F1))
  orig <- stamp(orig, r$methods)

} else {  # geuvadis_real: table + per-gene selections
  r <- replace_rows(orig$table, new$table, c("gene", "prior", "method", "selection"))
  cat("geuvadis real:", basename(cfg$orig), "\n")
  report(transform(orig$table, F1 = rep_rate), transform(r$rows, F1 = rep_rate), r$methods)
  cat("  (column reported as 'F1' above is rep_rate)\n")
  orig$table <- r$rows
  pat <- paste0("\\.(", paste(r$methods, collapse = "|"), ")\\.")
  for (i in seq_along(orig$genes)) {
    j <- match(orig$genes[i], new$genes)
    keep <- orig$selections[[i]][!grepl(pat, names(orig$selections[[i]]))]
    add  <- if (is.na(j)) list() else new$selections[[j]][grepl(pat, names(new$selections[[j]]))]
    orig$selections[[i]] <- c(keep, add)
  }
  orig <- stamp(orig, r$methods)
}

dir.create(dirname(cfg$out), showWarnings = FALSE, recursive = TRUE)
saveRDS(orig, cfg$out)
cat("written:", cfg$out, "\n")
cat(sprintf("CANARY_STATUS %s %s\n", basename(cfg$out), canary_status))
