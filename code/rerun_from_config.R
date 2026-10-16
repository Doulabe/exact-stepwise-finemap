# =============================================================================
# File: rerun_from_config.R
# Project: exact-stepwise-finemap
# Author: N.K. Doulabe
#
# Re-runs a benchmark with the exact configuration stored in an existing result
# file (its $config element), overriding only the output file, the number of
# cores and any extra flags (e.g. --inhouse_only AsS,NaM,ExM plus the flags that
# skip the external comparators). Used to recompute a subset of methods after a
# code correction, before merge_fixed_methods.R.
#
# Usage (from code/):
#   Rscript rerun_from_config.R --orig ../results_orig/X.rds \
#       --script Benchmark_External.R --out ../results_rerun/X.rds \
#       --n_cores 28 --extra "--inhouse_only AsS,NaM,ExM --skip_susie TRUE"
# =============================================================================
cfg <- list(orig = "", script = "", out = "", n_cores = 1L, extra = "")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}
stopifnot(file.exists(cfg$orig), file.exists(cfg$script), nzchar(cfg$out))

old <- readRDS(cfg$orig)$config
skip_keys <- c("outfile", "n_cores")
extra <- if (nzchar(cfg$extra)) strsplit(trimws(cfg$extra), "[[:space:]]+")[[1L]] else character(0)
extra_keys <- sub("^--", "", extra[startsWith(extra, "--")])

args <- character(0)
for (k in setdiff(names(old), c(skip_keys, extra_keys))) {
  v <- old[[k]]
  if (is.null(v) || length(v) != 1L) next
  v <- if (is.logical(v)) (if (isTRUE(v)) "TRUE" else "FALSE") else as.character(v)
  if (!nzchar(v)) next              # empty grids fall back to the script defaults
  args <- c(args, paste0("--", k), v)
}
args <- c(args, extra, "--n_cores", as.character(cfg$n_cores), "--outfile", cfg$out)
cat("Rscript", cfg$script, paste(args, collapse = " "), "\n")
status <- system2("Rscript", c(cfg$script, shQuote(args)))
quit(status = status)
