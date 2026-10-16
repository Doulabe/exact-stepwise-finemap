# =============================================================================
# make_tab_timing.R
# Mean wall-clock time per GEUVADIS gene window (s), by method, for the fixed-K
# and the adaptive BIC-stopped runs, from the `table` of a Benchmark_GEUVADIS.R
# real-mode result (runtime column, averaged over genes and prior arms).
#
#   Rscript make_tab_timing.R --in results_timing/results_geuvadis_real_timing.rds \
#       --out staging_exact/tab_timing.tex [--rows_only TRUE]
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
get <- function(k, d) { i <- match(k, args); if (is.na(i)) d else args[i + 1L] }
fin <- get("--in", "results_timing/results_geuvadis_real_timing.rds")
out <- get("--out", "")
tb  <- readRDS(fin)$table
ord <- c("NaW", "AsW", "ExW", "AsS", "NaM", "ExM", "SuSiE", "LASSO")
mt  <- function(meth, sel) {
  v <- tb$runtime[tb$method == meth & tb$selection == sel]
  if (!length(v)) "--" else sprintf("%.2f", mean(v, na.rm = TRUE))
}
rows <- vapply(ord, function(m) sprintf("%s & %s & %s \\\\", m, mt(m, "fixedK"), mt(m, "BIC")), "")
if (nzchar(out)) writeLines(rows, out) else cat(rows, sep = "\n")
