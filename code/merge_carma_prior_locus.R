#!/usr/bin/env Rscript
# merge_carma_prior_locus.R
# Merge the CARMA-only arm into the existing
# locus-sim RDS for the deterministic-prior scenarios (informative/uniform/mixed),
# whose other-method values are unchanged. Misspecified is NOT merged here: it
# was fully re-run (so all methods share the seeded decoy). Run from project root.
suppressPackageStartupMessages(library(dplyr))

scns <- c("informative", "uniform", "misspecified", "mixed")
kms  <- c("fixed", "adaptive")
for (km in kms) for (scn in scns) {
  base_f <- sprintf("results_locus_sim_%s_%s.rds", scn, km)
  co_f   <- sprintf("carma_only_%s_%s.rds", scn, km)
  if (!file.exists(base_f) || !file.exists(co_f)) {
    cat("SKIP (missing):", base_f, "or", co_f, "\n"); next
  }
  base <- readRDS(base_f); co <- readRDS(co_f)
  # drop any prior CARMA_prior rows (idempotent re-merge)
  base$raw     <- base$raw[base$raw$method != "CARMA_prior", , drop = FALSE]
  base$summary <- base$summary[base$summary$method != "CARMA_prior", , drop = FALSE]
  # align columns and append
  base$raw     <- rbind(base$raw,     co$raw[, names(base$raw), drop = FALSE])
  base$summary <- rbind(base$summary, co$summary[, names(base$summary), drop = FALSE])
  saveRDS(base, base_f)
  np <- sum(base$summary$method == "CARMA_prior")
  cat(sprintf("merged CARMA_prior into %-40s (+%d summary rows)\n", base_f, np))
  file.remove(co_f)
}
cat("done.\n")
