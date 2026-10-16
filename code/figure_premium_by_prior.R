# =============================================================================
# File: figure_premium_by_prior.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Companion to Figure 1 (gain of the exact Wald law against sample size):
# the same sweep under the uniform and misspecified priors, isolating what the
# exact calibration contributes without an informative prior. The gain is
# specific to the informative-prior setting: without it, selection reduces to
# (near-)monotone functions of the same statistics and the exact and
# asymptotic variants select essentially the same candidates.
#
# Usage (from code/):
#   Rscript figure_premium_by_prior.R \
#     --informative ../../server_rerun_fix_20260928/results_exact/results_synth_premium_vs_n.rds \
#     --uniform ../results_synth_premium_vs_n_uniform.rds \
#     --misspecified ../results_synth_premium_vs_n_misspecified.rds \
#     --out ../Figure_premium_by_prior.pdf
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(informative = "../../server_rerun_fix_20260928/results_exact/results_synth_premium_vs_n.rds",
            uniform = "../results_synth_premium_vs_n_uniform.rds",
            misspecified = "../results_synth_premium_vs_n_misspecified.rds",
            out = "../Figure_premium_by_prior.pdf")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}

prem <- function(f, lab) {
  readRDS(f)$all_results %>%
    filter(method %in% c("Asymptotic_Wald", "Exact_Wald")) %>%
    group_by(method, n) %>% summarise(M = mean(M), .groups = "drop") %>%
    pivot_wider(names_from = method, values_from = M) %>%
    transmute(n = n, d = Exact_Wald - Asymptotic_Wald, prior = lab) %>% arrange(n)
}
d <- bind_rows(prem(cfg$informative, "Informative"),
               prem(cfg$uniform, "Uniform"),
               prem(cfg$misspecified, "Misspecified"))
print(as.data.frame(pivot_wider(d, names_from = prior, values_from = d)), row.names = FALSE)

cols <- c(Informative = "#D55E00", Uniform = "#0072B2", Misspecified = "#009E73")
pch  <- c(Informative = 16, Uniform = 17, Misspecified = 15)
pdf(cfg$out, width = 5.8, height = 4.2, pointsize = 10)
par(mar = c(4, 4, 1, 1) + 0.1)
plot(NA, xlim = range(d$n), ylim = range(c(0, d$d)), log = "x",
     xlab = "Sample size n", ylab = "Extra causal SNPs recovered (ExW - AsW)")
abline(h = 0, col = "grey70", lty = 3)
for (p in names(cols)) {
  dd <- d %>% filter(prior == p)
  lines(dd$n, dd$d, col = cols[p], lwd = 1.8)
  points(dd$n, dd$d, col = cols[p], pch = pch[p], cex = 0.9)
}
legend("topright", legend = names(cols), col = cols, pch = pch, lwd = 1.8, bty = "n")
invisible(dev.off())
cat("Wrote", cfg$out, "\n")
