# =============================================================================
# File: make_tab_prior_scale.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Prior-amplitude sensitivity of the GEUVADIS application: the three-level
# scale {0.3, 0.5, 0.7} of the real-data priors is remapped onto a weaker
# {0.4, 0.5, 0.6} and a stronger {0.2, 0.5, 0.8} contrast (Benchmark_GEUVADIS
# --prior_scale) and the fixed-K selections are rerun. This table reports the
# nominal LD-aware replication rates per method, scale and prior arm. The
# uniform arm (constant 0.5) is invariant by construction and serves as an
# internal check.
#
# Usage (from code/):
#   Rscript make_tab_prior_scale.R --out ../tab_geuvadis_prior_scale.tex
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(tidyr)})

cfg <- list(primary = "../results_geuvadis_real_ld.rds",
            weak    = "../results_geuvadis_real_scale_weak_ld.rds",
            strong  = "../results_geuvadis_real_scale_strong_ld.rds",
            out     = "../tab_geuvadis_prior_scale.tex")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- a[i + 1L]
}

METHODS <- c("AsS", "ExW", "NaM", "ExM")
PRIORS  <- c(uniform = "Uniform", polyfun = "PolyFun", tss = "TSS",
             discordant = "Discordant")
SCALES  <- c(weak    = "$\\{0.4, 0.5, 0.6\\}$",
             primary = "$\\{0.3, 0.5, 0.7\\}$",
             strong  = "$\\{0.2, 0.5, 0.8\\}$")

agg <- function(f, sc) readRDS(f)$table %>%
  filter(selection == "fixedK", method %in% METHODS, prior %in% names(PRIORS)) %>%
  group_by(method, prior) %>%
  summarise(v = mean(rep_rate_ld, na.rm = TRUE), .groups = "drop") %>%
  mutate(scale = sc)
d <- bind_rows(agg(cfg$weak, "weak"), agg(cfg$primary, "primary"),
               agg(cfg$strong, "strong"))

fmt <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
rows <- unlist(lapply(METHODS, function(m) {
  r <- vapply(names(SCALES), function(sc) {
    cells <- vapply(names(PRIORS), function(p) {
      z <- d %>% filter(method == m, scale == sc, prior == p)
      if (!nrow(z)) "--" else fmt(z$v)
    }, "")
    paste0(if (sc == "weak") paste0("\\multirow{3}{*}{", m, "} ") else "",
           "& ", SCALES[sc], " & ", paste(cells, collapse = " & "), " \\\\")
  }, "")
  c(r, "\\addlinespace")
}))
rows <- rows[-length(rows)]

writeLines(c(
  "\\begin{table}[htbp]",
  paste0("\\caption{\\rev{\\textbf{Sensitivity to the amplitude of the prior weights.} ",
         "Mean LD-aware replication rate of the fixed-$K=5$ selections (nominal test) when the ",
         "three-level prior scale of the main analysis, $\\{0.3, 0.5, 0.7\\}$, is replaced by a ",
         "weaker $\\{0.4, 0.5, 0.6\\}$ or a stronger $\\{0.2, 0.5, 0.8\\}$ contrast, level by ",
         "level. The uniform arm is constant at $0.5$ and hence invariant by construction. Its ",
         "repetition across scales is an internal consistency check.}}"),
  "\\label{tab:geuvadis_prior_scale}",
  "\\centering\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{llcccc}", "\\toprule",
  paste0("Method & Scale & ", paste(PRIORS, collapse = " & "), " \\\\"), "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"), cfg$out)
cat("Wrote", cfg$out, "\n")
print(as.data.frame(d %>% pivot_wider(names_from = prior, values_from = v) %>%
  arrange(match(method, METHODS), match(scale, names(SCALES)))), row.names = FALSE)
