# =============================================================================
# File: install_dependencies.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe (review additions)
#
# Description:
#   Install all R packages required to reproduce the manuscript from a clean
#   R installation. The script:
#     1. Lists the canonical set of CRAN dependencies with minimum versions.
#     2. Checks which packages are missing or out-of-date.
#     3. Installs them from CRAN using a configurable mirror.
#     4. Prints a summary of what was installed and what is now available.
#
# Usage:
#   Interactive:
#     setwd("path/to/Projet1_donnees_brutes")
#     source("install_dependencies.R")
#     install_dependencies()                            # install missing/out-of-date
#     install_dependencies(verify_only = TRUE)          # report without installing
#     install_dependencies(mirror = "https://cran.r-project.org")
#
#   Command line:
#     Rscript install_dependencies.R                    # install
#     Rscript install_dependencies.R --verify           # report only
#     Rscript install_dependencies.R --mirror <URL>
# =============================================================================

# ----------------------------------------------------------------------------
# 0. Pretty logging (re-used from run_all.R)
# ----------------------------------------------------------------------------
.term <- isatty(stdout())
.col  <- function(code, x) if (.term) paste0("\033[", code, "m", x, "\033[0m") else x
.stage <- function(...) cat(.col("1;34", "[stage]"), ..., "\n")
.info  <- function(...) cat(.col("32", "[info] "), ..., "\n")
.warn  <- function(...) cat(.col("33", "[warn] "), ..., "\n", file = stderr())
.err   <- function(...) cat(.col("1;31", "[err]  "), ..., "\n", file = stderr())
.ok    <- function(...) cat(.col("32", "  ok   "), ..., "\n")
.miss  <- function(...) cat(.col("31", "  miss "), ..., "\n")
.old   <- function(...) cat(.col("33", "  old  "), ..., "\n")

# ----------------------------------------------------------------------------
# 1. Canonical dependency list (CRAN)
# ----------------------------------------------------------------------------
# Format: list(pkg = "name", min_version = "x.y.z" or NULL, reason = "...")
.dependencies <- list(
  # ----- Core base/recommended packages (usually pre-installed) -----
  list(pkg = "MASS",        min_version = NULL,
       reason = "Linear algebra utilities used in Test_statistics_functions"),
  list(pkg = "Matrix",      min_version = NULL,
       reason = "Sparse / dense matrix operations"),
  list(pkg = "stats",       min_version = NULL,
       reason = "Distributions, base statistics"),
  list(pkg = "methods",     min_version = NULL,
       reason = "Type coercion in CLI parsers"),
  list(pkg = "utils",       min_version = NULL,
       reason = "install.packages, write.csv, etc."),

  # ----- Tidy data manipulation & plotting -----
  list(pkg = "dplyr",       min_version = "1.0.0",
       reason = "Aggregation pipelines in fill_table / Benchmark_External"),
  list(pkg = "tidyr",       min_version = NULL,
       reason = "Wide/long reshaping in the table generators"),
  list(pkg = "ggplot2",     min_version = "3.3.0",
       reason = "Figures (e.g. Figure_premium_vs_n)"),
  list(pkg = "scales",      min_version = NULL,
       reason = "Colour/axis scaling in the LD-structure figures"),
  list(pkg = "data.table",  min_version = NULL,
       reason = "Fast I/O helper (Load_packages)"),

  # ----- Probability / numerical primitives -----
  list(pkg = "mvtnorm",     min_version = NULL,
       reason = "Multivariate normal sampling in Generate_XY_realistic"),
  list(pkg = "statmod",     min_version = NULL,
       reason = "Gauss-Hermite quadrature in optimised test statistics"),

  # ----- Parallel backend -----
  list(pkg = "doParallel",  min_version = NULL,
       reason = "Parallel foreach backend in Benchmark_External"),
  list(pkg = "foreach",     min_version = NULL,
       reason = "Companion to doParallel"),

  # ----- BayesC comparator -----
  list(pkg = "BGLR",        min_version = NULL,
       reason = "BayesC spike-and-slab MCMC comparator"),

  # ----- External comparators -----
  list(pkg = "susieR",      min_version = "0.12.0",
       reason = "External fine-mapping benchmark (Wang et al. 2020)"),
  list(pkg = "varbvs",      min_version = NULL,
       reason = "External variational Bayesian VS (Carbonetto & Stephens 2012)"),
  list(pkg = "glmnet",      min_version = "4.0",
       reason = "LASSO baseline (Tibshirani 1996)"),
  list(pkg = "CARMA",       min_version = NULL,
       reason = "Locus-restricted fine-mapping comparator (Yang et al. 2023)"),
  list(pkg = "duckdb",      min_version = NULL,
       reason = "Reads the PolyFun snpvar Parquet (prep_geuvadis_gene.R)")
)
# Note: FINEMAP is an external command-line binary, not an R package. Install it
# separately (v1.3.1) and put it on the PATH or set the FINEMAP_BIN environment
# variable. It is only needed for the locus-restricted simulation and the
# real-data locus analysis.

# ----------------------------------------------------------------------------
# 2. Inspection helpers
# ----------------------------------------------------------------------------
.installed_packages_cached <- function() {
  rownames(utils::installed.packages())
}

.pkg_status <- function(spec, installed) {
  pkg <- spec$pkg
  if (!(pkg %in% installed)) {
    return(list(status = "missing", current = NA_character_))
  }
  v <- as.character(utils::packageVersion(pkg))
  if (!is.null(spec$min_version) &&
      utils::compareVersion(v, spec$min_version) < 0L) {
    return(list(status = "outdated", current = v))
  }
  list(status = "ok", current = v)
}

# ----------------------------------------------------------------------------
# 3. Main entry point
# ----------------------------------------------------------------------------
install_dependencies <- function(verify_only = FALSE,
                                 mirror      = "https://cloud.r-project.org",
                                 dependencies = c("Depends", "Imports"),
                                 quiet       = FALSE) {
  .stage("Inspecting installed R packages")
  installed <- .installed_packages_cached()

  status_tbl <- data.frame(
    pkg         = vapply(.dependencies, `[[`, character(1L), "pkg"),
    min_version = vapply(.dependencies, function(d) {
      if (is.null(d$min_version)) NA_character_ else d$min_version
    }, character(1L)),
    current     = NA_character_,
    status      = NA_character_,
    reason      = vapply(.dependencies, `[[`, character(1L), "reason"),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(.dependencies)) {
    s <- .pkg_status(.dependencies[[i]], installed)
    status_tbl$status[i]  <- s$status
    status_tbl$current[i] <- s$current
    pkg <- status_tbl$pkg[i]
    msg <- sprintf("%-12s  %s%s",
                   pkg,
                   if (is.na(s$current)) "-" else s$current,
                   if (!is.na(status_tbl$min_version[i]))
                     sprintf("  (min %s)", status_tbl$min_version[i]) else "")
    switch(s$status,
      ok       = .ok(msg),
      outdated = .old(msg),
      missing  = .miss(msg))
  }

  to_install <- status_tbl$pkg[status_tbl$status %in% c("missing", "outdated")]
  ok_count   <- sum(status_tbl$status == "ok")
  .info(sprintf("Summary: %d ok / %d need (re)install (out of %d total)",
                ok_count, length(to_install), nrow(status_tbl)))

  if (length(to_install) == 0L) {
    .info("All required packages are present and up-to-date.")
  } else if (verify_only) {
    .warn("--verify mode: not installing. Missing/outdated packages:")
    cat("  ", paste(to_install, collapse = ", "), "\n")
  } else {
    .stage(sprintf("Installing %d packages from %s", length(to_install), mirror))
    tryCatch({
      options(repos = c(CRAN = mirror))
      utils::install.packages(to_install,
                              dependencies = dependencies,
                              quiet = quiet)
    }, error = function(e) {
      .err(sprintf("install.packages failed: %s", conditionMessage(e)))
      .warn("Try a different mirror via install_dependencies(mirror = ...).")
    })

    # Re-inspect to verify
    installed <- .installed_packages_cached()
    failed <- character(0L)
    for (p in to_install) {
      v <- if (p %in% installed) as.character(utils::packageVersion(p)) else NA
      if (is.na(v)) failed <- c(failed, p)
    }
    if (length(failed) > 0L) {
      .err(sprintf("Failed to install: %s",
                   paste(failed, collapse = ", ")))
    } else {
      .info("All requested packages installed successfully.")
    }
  }

  # --- Final summary table ---
  cat("\n", .col("1", "Final package status:\n"), sep = "")
  print(status_tbl, row.names = FALSE)
  cat("\n")

  invisible(list(table = status_tbl,
                 installed = to_install))
}

# ----------------------------------------------------------------------------
# 4. CLI entry point
# ----------------------------------------------------------------------------
.parse_cli <- function(args) {
  cfg <- list(verify_only = FALSE,
              mirror      = "https://cloud.r-project.org",
              quiet       = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      cat(
        "Usage: Rscript install_dependencies.R [options]\n",
        "  --verify            Only report missing/outdated, do not install\n",
        "  --mirror <URL>      CRAN mirror (default: https://cloud.r-project.org)\n",
        "  --quiet             Quiet install (suppress download progress)\n",
        "  -h, --help          Show this message\n",
        sep = "")
      quit(save = "no", status = 0L)
    }
    if (a == "--verify") {
      cfg$verify_only <- TRUE; i <- i + 1L
    } else if (a == "--quiet") {
      cfg$quiet <- TRUE; i <- i + 1L
    } else if (a == "--mirror" && i < length(args)) {
      cfg$mirror <- args[i + 1L]; i <- i + 2L
    } else {
      .warn(sprintf("Unknown option ignored: %s", a))
      i <- i + 1L
    }
  }
  cfg
}

if (!interactive() && sys.nframe() == 0L) {
  cfg <- .parse_cli(commandArgs(trailingOnly = TRUE))
  do.call(install_dependencies, cfg)
}
