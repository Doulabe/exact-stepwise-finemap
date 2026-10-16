#!/usr/bin/env Rscript
# =============================================================================
# figure_premium_vs_n.R
# Builds the main-text "finite-sample gain vs sample size" figure: the extra
# number of causal SNPs recovered by each EXACT variant relative to its
# asymptotic/naive counterpart, on the synthetic AR(1) design (informative prior,
# gamma=0.6, tau=0.2), as a function of n. The gain is large at small n and
# decays to zero as n grows.
#
# Inputs (produced by the in-house-only synthetic sweep, e.g.
#   Rscript Benchmark_External.R --design synthetic --m_synth 1000 --S_mc 30 \
#     --rho_grid 0.95 --n_grid 50,80,150,300,600,1000 --gamma_grid 0.6 --tau_grid 0.2 \
#     --scenario informative --K_mode fixed --B 20 \
#     --skip_susie TRUE --skip_lasso TRUE --skip_bayesc TRUE --skip_varbvs TRUE \
#     --outfile ../results_synth_premium_vs_n.rds
# and optionally a second file over the other rho values for a min--max band):
#   ../results_synth_premium_vs_n.rds        (required)
#   ../results_synth_premium_vs_n_band.rds   (optional; adds the rho band)
#
# Run FROM code/. Output: ../Figure_premium_vs_n.pdf
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

files <- c("../results_synth_premium_vs_n.rds", "../results_synth_premium_vs_n_band.rds")
files <- files[file.exists(files)]
if (!length(files)) stop("No premium-sweep RDS found in ../ (see header for how to produce it).")
s <- bind_rows(lapply(files, function(f) readRDS(f)$summary))

# common n grid across all rho values present
ns <- s %>% filter(tau == 0.2, gamma == 0.6) %>%
  count(n, rho_intra) %>% count(n) %>% filter(nn == max(nn)) %>% pull(n)

prem <- s %>% filter(tau == 0.2, gamma == 0.6, n %in% ns) %>%
  select(n, rho_intra, method, M = mean_M) %>%
  pivot_wider(names_from = method, values_from = M) %>%
  transmute(n, rho_intra,
            `Exact vs. asymptotic Wald`    = Exact_Wald - Asymptotic_Wald,
            `Exact vs. naive Multiple-SNP` = Exact_Multiple - Naive_Multiple) %>%
  pivot_longer(c(-n, -rho_intra), names_to = "Contrast", values_to = "prem") %>%
  group_by(n, Contrast) %>%
  summarise(mean = mean(prem), lo = min(prem), hi = max(prem), .groups = "drop")

p <- ggplot(prem, aes(n, mean, colour = Contrast, fill = Contrast)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey55") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_x_log10(breaks = sort(unique(prem$n))) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) +
  labs(x = "Sample size  n  (log scale)",
       y = "Extra causal SNPs recovered  (/15)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
        plot.margin = margin(t = 6, r = 10, b = 4, l = 8),
        axis.title.y = element_text(margin = margin(r = 6)))

ggsave("../Figure_premium_vs_n.pdf", p, width = 6.6, height = 4.2)
cat("Wrote ../Figure_premium_vs_n.pdf\n")
