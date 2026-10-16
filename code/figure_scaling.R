# =============================================================================
# File: figure_scaling.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Empirical computational scaling: per-window wall-clock time of the fixed-K
# runs against the number of candidate variants m, across the hundred GEUVADIS
# windows (log-log), with the fitted power log t = a + b log m per method.
#
# Usage (from code/):
#   Rscript figure_scaling.R \
#     --timing ../../server_rerun_fix_20260928/results_timing/results_geuvadis_real_timing.rds \
#     --rds_dir ../geuvadis_rds --out ../Figure_scaling_m.pdf
# =============================================================================
suppressPackageStartupMessages(library(dplyr))

cfg <- list(timing = "../results_geuvadis_real_timing.rds",
            rds_dir = "../geuvadis_rds",
            out = "../Figure_scaling_m.pdf")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}

tb <- readRDS(cfg$timing)$table %>% filter(selection == "fixedK")
genes <- unique(tb$gene)
mlen <- vapply(genes, function(g)
  length(readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))$pos), 1L)
tb$m <- mlen[tb$gene]
METH <- c("NaW", "AsS", "NaM", "ExM", "SuSiE", "LASSO")
tb <- tb %>% filter(method %in% METH) %>%
  group_by(gene, method, m) %>% summarise(s = mean(runtime, na.rm = TRUE), .groups = "drop")

slopes <- tb %>% group_by(method) %>%
  summarise(b = round(coef(lm(log(s) ~ log(m)))[2], 2), .groups = "drop")
cat("Pentes log-log (temps ~ m^b):\n"); print(as.data.frame(slopes), row.names = FALSE)

cols <- c(NaW = "#999999", AsS = "#56B4E9", NaM = "#E69F00", ExM = "#D55E00",
          SuSiE = "#0072B2", LASSO = "#009E73")
pdf(cfg$out, width = 5.8, height = 4.4, pointsize = 10)
par(mar = c(4, 4, 1, 1) + 0.1)
plot(NA, xlim = range(tb$m), ylim = range(tb$s), log = "xy",
     xlab = "Candidate variants per window (m)",
     ylab = "Wall-clock time per window (s)")
for (mm in METH) {
  d <- tb %>% filter(method == mm) %>% arrange(m)
  points(d$m, d$s, col = adjustcolor(cols[mm], 0.45), pch = 16, cex = 0.55)
  f <- lm(log(s) ~ log(m), data = d)
  xx <- range(d$m)
  lines(xx, exp(predict(f, data.frame(m = xx))), col = cols[mm], lwd = 1.8)
}
legend("topleft", legend = sprintf("%s  (b = %.2f)", slopes$method,
       slopes$b)[match(METH, slopes$method)],
       col = cols[METH], lwd = 1.8, bty = "n", cex = 0.82)
invisible(dev.off())
cat("Wrote", cfg$out, "\n")
