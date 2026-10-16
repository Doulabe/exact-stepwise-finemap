# =============================================================================
# File: run_qtltools_geuvadis100.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Native QTLtools comparator for the GEUVADIS application: the field-standard
# conditional cis scan (Delaneau et al. 2017), run with the QTLtools 1.3.1
# binary on the same inputs as the proposed methods. For each of the 100 cis
# windows: a 1,000-permutation pass calibrates the gene-level 0.05 nominal
# cutoff through the beta approximation (same degrees-of-freedom conversion as
# the official qtltools_runFDR_cis.R), then the conditional pass
# (--mapping) returns the ordered independent signals (forward-backward scan).
# Genotypes are passed as dosage BED (no artificial alleles), the phenotype is
# the per-population rank-INT quantification used throughout Section 4, and
# sex plus population indicators are passed through --cov, so QTLtools applies
# its own covariate handling to the same covariates.
#
# Usage (from code/):
#   Rscript run_qtltools_geuvadis100.R --n_cores 3 \
#     --out ../results_qtltools_geuvadis.rds
# =============================================================================
suppressPackageStartupMessages({library(parallel)})

cfg <- list(rds_dir = "../geuvadis_rds",
            bin = file.path("/Users/kossi/Desktop/Dossier_these/Redaction/LMM_Gaussian_Bayesian",
                            "cbf-lmm/research/variance_uncertainty/software_comparison",
                            "vendor/qtltools-master/bin/QTLtools"),
            work = "", n_cores = 3L, seed = 26000L,
            out = "../results_qtltools_geuvadis.rds")
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}
stopifnot(file.exists(cfg$bin))
if (!nzchar(cfg$work)) cfg$work <- file.path(dirname(cfg$out), "qtltools_work")
dir.create(cfg$work, showWarnings = FALSE, recursive = TRUE)

residualize <- function(M, C) {
  q <- qr(cbind(1, as.matrix(C)))
  qr.resid(q, as.matrix(M))
}

genes <- sort(sub("^geuvadis_(ENSG[0-9]+)\\.rds$", "\\1",
                  grep("^geuvadis_ENSG", list.files(cfg$rds_dir), value = TRUE)))
stopifnot(length(genes) >= 1)

one_gene <- function(gi) {
  g <- genes[gi]
  gx <- readRDS(file.path(cfg$rds_dir, sprintf("geuvadis_%s.rds", g)))
  dest <- file.path(cfg$work, g); dir.create(dest, showWarnings = FALSE)
  # same candidate set as run_real_gene (drop columns degenerate after FWL)
  keep <- apply(residualize(gx$G_eur, gx$covar_eur), 2, sd) > 1e-8
  G   <- gx$G_eur[, keep, drop = FALSE]
  pos <- gx$pos[keep]
  smp <- gx$samples_eur
  ids <- sprintf("v%06d_%d", seq_along(pos), pos)
  o <- order(pos)
  bed <- data.frame(chr = gx$chr, start = pos[o] - 1L, end = pos[o],
                    id = ids[o], gid = ids[o], strand = "+",
                    t(G[, o, drop = FALSE]), check.names = FALSE)
  colnames(bed) <- c("#chr", "start", "end", "id", "gid", "strand", smp)
  data.table::fwrite(bed, file.path(dest, "genotypes.bed"), sep = "\t", quote = FALSE)
  ph <- data.frame(chr = gx$chr, start = gx$tss - 1L, end = gx$tss, id = g,
                   gid = g, strand = "+", t(as.numeric(gx$y_eur)), check.names = FALSE)
  colnames(ph) <- colnames(bed)
  data.table::fwrite(ph, file.path(dest, "expression.bed"), sep = "\t", quote = FALSE)
  for (nm in c("genotypes", "expression")) {
    fn <- file.path(dest, paste0(nm, ".bed"))
    stopifnot(system2("bgzip", c("-f", shQuote(fn))) == 0,
              system2("tabix", c("-f", "-p", "bed", shQuote(paste0(fn, ".gz")))) == 0)
  }
  cv <- t(as.matrix(gx$covar_eur))
  cov_file <- file.path(dest, "covariates.txt")
  writeLines(c(paste(c("id", smp), collapse = "\t"),
               vapply(seq_len(nrow(cv)), function(r)
                 paste(c(rownames(cv)[r], cv[r, ]), collapse = "\t"), "")), cov_file)

  W <- max(abs(pos - gx$tss)) + 1L
  common <- c("cis", "--vcf", shQuote(file.path(dest, "genotypes.bed.gz")),
              "--bed", shQuote(file.path(dest, "expression.bed.gz")),
              "--cov", shQuote(cov_file),
              "--window", format(W, scientific = FALSE), "--seed", cfg$seed + gi)
  per <- file.path(dest, "permutation.txt")
  t0 <- proc.time()[["elapsed"]]
  # resume-friendly: reuse existing outputs; retry once on a transient failure
  run2 <- function(args, log) {
    s <- system2(cfg$bin, args, stdout = log, stderr = log)
    if (s != 0) { Sys.sleep(5); s <- system2(cfg$bin, args, stdout = log, stderr = log) }
    s
  }
  if (!file.exists(per) || file.info(per)$size == 0)
    stopifnot(run2(c(common, "--permute", "1000", "--out", shQuote(per)),
                   file.path(dest, "perm.log")) == 0)
  p <- read.table(per)
  stopifnot(nrow(p) == 1, ncol(p) == 20)
  # gene-level alpha = 0.05 cutoff; dof conversion as in qtltools_runFDR_cis.R
  pb <- qbeta(.05, p$V14, p$V15)
  tt <- qf(pb, 1, p$V13, lower.tail = FALSE)
  corr2 <- tt / (p$V13 + tt)
  cutoff <- pf(p$V12 * corr2 / (1 - corr2), 1, p$V12, lower.tail = FALSE)
  thr <- file.path(dest, "threshold.txt")
  write.table(data.frame(g, cutoff), thr, row.names = FALSE, col.names = FALSE, quote = FALSE)
  mapping <- file.path(dest, "mapping.txt")
  if (!file.exists(mapping))
    stopifnot(run2(c(common, "--mapping", shQuote(thr), "--out", shQuote(mapping)),
                   file.path(dest, "map.log")) == 0)
  tab <- if (file.info(mapping)$size > 0) read.table(mapping) else data.frame()
  hits <- if (nrow(tab)) tab[tab$V21 == 1, , drop = FALSE] else tab
  sel_pos <- integer(0)
  if (nrow(hits)) {
    hits <- hits[order(hits$V12), , drop = FALSE]      # by signal rank
    sel_pos <- as.integer(sub("^v[0-9]+_", "", hits$V8))
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  cat(sprintf("[%s] QTLtools K=%d cutoff=%.3g (%.1f s)\n", g, length(sel_pos), cutoff, elapsed))
  list(gene = g, sel_pos = sel_pos, K = length(sel_pos), cutoff = cutoff,
       beta = c(p$V14, p$V15), elapsed = elapsed)
}

res <- mclapply(seq_along(genes), one_gene, mc.cores = cfg$n_cores, mc.preschedule = FALSE)
names(res) <- genes
bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r$gene), TRUE)
if (any(bad)) stop("failed genes: ", paste(genes[bad], collapse = ", "))
saveRDS(list(results = res, genes = genes,
             config = cfg[c("seed", "bin")]), cfg$out)
cat("Wrote", cfg$out, "\n")
cat("K distribution:\n"); print(table(vapply(res, `[[`, 1L, "K")))
