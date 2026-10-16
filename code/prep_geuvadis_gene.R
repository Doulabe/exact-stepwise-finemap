# =============================================================================
# File: prep_geuvadis_gene.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Build one RDS per gene for the GEUVADIS cis-eQTL application: cis-window
# genotypes split into an EUR discovery set and a YRI replication set, the
# normalized expression phenotype, and four prior specifications spanning a
# fidelity gradient (uniform / PolyFun off-the-shelf / TSS-distance / discordant).
#
# Inputs (local files, see README):
#   GEUVADIS_DIR (env), containing
#     * per-chromosome genotype VCFs (E-GEUV-1, GRCh37, 465 samples):
#         GEUVADIS.chr<CHR>.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz (+ .tbi)
#     * PEER-residualized gene quantifications (462 samples):
#         GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
#       (leading columns: TargetID, Gene_Symbol, Chr, Coord; Coord = TSS proxy;
#        the ".norm" per-gene standard-normal step is NOT distributed by EBI --
#        it is applied here as a rank-based INT, see --rank_normal)
#     * 1000G population panel (sample / pop / super_pop / gender):
#         integrated_call_samples_v3.20130502.ALL.panel
#   POLYFUN_DIR (env): PolyFun snpvar parquet(s), read with duckdb
#     (snpvar_meta.chr1_7.parquet and, for chr >= 8, snpvar_meta.chr8_22.parquet).
#
# Output RDS (one per gene): list(gene, symbol, chr, tss, window, pos, maf_eur,
#   G_eur, G_yri, y_eur, y_yri, snpvar, priors = list(uniform, polyfun, tss,
#   discordant), covar_eur, covar_yri, samples_eur, samples_yri, cfg).
#   NOTE: the LD matrix R is NOT stored (one dense
#   m x m block per gene is prohibitive across genes); compute it at analysis
#   time from G_eur / G_yri when a comparator needs it.
#
# Usage (from code/):
#   Rscript prep_geuvadis_gene.R --gene ENSG00000198502 --out ../geuvadis_ENSG00000198502.rds
#   Rscript prep_geuvadis_gene.R --genes_file ../geuvadis_genes.txt --out_dir ../geuvadis_rds
# =============================================================================
suppressPackageStartupMessages({ library(duckdb) })

DG <- Sys.getenv("GEUVADIS_DIR", "../geuvadis_data")
D4 <- Sys.getenv("POLYFUN_DIR",  "../polyfun_data")

cfg <- list(
  gene       = "",                     # ENSG id (version suffix ok) or gene symbol
  genes_file = "",                     # one gene per line (batch mode)
  rank_normal= TRUE,                   # per-gene, per-population rank-based inverse
                                       # normal transform of expression (the ".norm"
                                       # step of the GEUVADIS pipeline; the EBI
                                       # distribution ships only the PEER-residualized
                                       # resk10 file, so we apply INT here)
  window     = 1000000L,               # cis half-window around the TSS (bp)
  maf_lo     = 0.05,                   # MAF filter, computed on the EUR set
  tss_near   = 10000L,                 # |pos-TSS| <  tss_near -> p = 0.7
  tss_mid    = 100000L,                # |pos-TSS| <  tss_mid  -> p = 0.5, else 0.3
  vcf_pattern= file.path(DG, "GEUVADIS.chr%s.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"),
  expr_file  = file.path(DG, "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"),
  panel_file = file.path(DG, "integrated_call_samples_v3.20130502.ALL.panel"),
  out        = "",                     # single-gene output RDS
  out_dir    = "../geuvadis_rds"       # batch-mode output directory
)
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}
EUR_POPS <- c("CEU", "FIN", "GBR", "TSI")

# ---- 0. gene list -----------------------------------------------------------
genes <- if (nzchar(cfg$genes_file)) readLines(cfg$genes_file) else cfg$gene
genes <- trimws(genes[nzchar(trimws(genes))])
if (!length(genes)) stop("Provide --gene or --genes_file.")

# ---- 1. expression matrix + panel (loaded once) -----------------------------
cat("Reading expression matrix:", cfg$expr_file, "\n")
expr <- read.delim(gzfile(cfg$expr_file), check.names = FALSE, stringsAsFactors = FALSE)
ann_cols <- match(c("TargetID", "Gene_Symbol", "Chr", "Coord"), colnames(expr))
if (anyNA(ann_cols)) stop("Expression file: expected leading columns TargetID, Gene_Symbol, Chr, Coord.")
expr_samples <- setdiff(colnames(expr), colnames(expr)[ann_cols])
cat(sprintf("  %d genes x %d samples\n", nrow(expr), length(expr_samples)))

panel <- read.delim(cfg$panel_file, stringsAsFactors = FALSE)
panel <- panel[, c("sample", "pop", "super_pop", "gender")]

strip_ver <- function(x) sub("\\..*$", "", x)

# rank-based inverse normal transform (per gene, applied per population)
maybe_int <- function(y) {
  if (!isTRUE(cfg$rank_normal)) return(y)
  qnorm((rank(y, ties.method = "average") - 0.5) / length(y))
}

# ---- 2. per-gene builder ----------------------------------------------------
build_gene <- function(g) {
  row <- which(strip_ver(expr$TargetID) == strip_ver(g) | expr$Gene_Symbol == g)
  if (length(row) != 1L) stop(sprintf("Gene '%s': %d matches in expression matrix.", g, length(row)))
  chr <- as.character(expr$Chr[row]); tss <- as.integer(expr$Coord[row])
  lo  <- max(1L, tss - cfg$window); hi <- tss + cfg$window
  cat(sprintf("[%s|%s] chr%s TSS=%d, cis window %d-%d\n",
              strip_ver(expr$TargetID[row]), expr$Gene_Symbol[row], chr, tss, lo, hi))

  # -- genotypes from the E-GEUV-1 VCF (local file or remote indexed URL; a
  #    1000G phase-3 VCF also works, since the GEUVADIS individuals are 1000G
  #    samples, and the sample intersection below handles the panel
  #    difference). Samples are restricted at bcftools level to
  #    those present in the expression matrix. ---------------------------------
  vcf <- sprintf(cfg$vcf_pattern, chr)
  is_url <- grepl("^(https?|ftp)://", vcf)
  if (!is_url && !file.exists(vcf)) stop("VCF not found: ", vcf)
  vcf_samples <- system(sprintf("bcftools query -l %s 2>/dev/null", shQuote(vcf)), intern = TRUE)
  if (!length(intersect(vcf_samples, expr_samples)))
    stop("No overlap between VCF samples and expression samples.")
  # NOTE: no bcftools -S sample subsetting here. The E-GEUV-1 headers do not
  # declare the FORMAT tags PP/BD, and bcftools refuses sample-subset mode on
  # undeclared tags ("Undefined tags in the header"); all samples are queried
  # and the EUR/YRI subsetting is done by name in R below.
  gfile <- tempfile(fileext = ".gt")
  run_extract <- function(region) {
    cmd <- sprintf(
      "bcftools view -r %s -m2 -M2 -v snps %s 2>/dev/null | bcftools query -f '%%POS[\\t%%GT]\\n' > %s",
      region, shQuote(vcf), shQuote(gfile))
    system(cmd); readLines(gfile)
  }
  raw <- run_extract(sprintf("%s:%d-%d", chr, lo, hi))
  if (!length(raw))                       # chr-prefixed contig names
    raw <- run_extract(sprintf("chr%s:%d-%d", chr, lo, hi))
  unlink(gfile)
  if (!length(raw)) stop("bcftools returned no variants for the cis window.")
  parts <- strsplit(raw, "\t", fixed = TRUE)
  pos_g <- as.integer(vapply(parts, `[`, "", 1L))
  gt2dose <- function(v) {
    d <- integer(length(v))
    d[v %in% c("0|1", "1|0", "0/1", "1/0")] <- 1L
    d[v %in% c("1|1", "1/1")]               <- 2L
    d[v %in% c(".|.", "./.")]               <- NA_integer_
    d
  }
  G <- t(vapply(parts, function(p) gt2dose(p[-1L]), integer(length(vcf_samples))))
  G <- t(G)                                        # samples x SNPs
  rownames(G) <- vcf_samples; colnames(G) <- pos_g

  # -- sample split: EUR discovery / YRI replication, matched to expression ---
  common  <- intersect(vcf_samples, expr_samples)
  pop     <- panel$pop[match(common, panel$sample)]
  s_eur   <- common[pop %in% EUR_POPS & !is.na(pop)]
  s_yri   <- common[pop == "YRI" & !is.na(pop)]
  cat(sprintf("  samples: %d EUR (discovery), %d YRI (replication)\n", length(s_eur), length(s_yri)))

  G_eur <- G[s_eur, , drop = FALSE]; G_yri <- G[s_yri, , drop = FALSE]

  # -- QC: EUR MAF filter; residual missing -> column-mean imputation ---------
  maf <- pmin(colMeans(G_eur, na.rm = TRUE) / 2, 1 - colMeans(G_eur, na.rm = TRUE) / 2)
  keep <- which(is.finite(maf) & maf >= cfg$maf_lo)
  G_eur <- G_eur[, keep, drop = FALSE]; G_yri <- G_yri[, keep, drop = FALSE]
  pos <- pos_g[keep]; maf <- maf[keep]
  impute_cm <- function(M) { if (anyNA(M)) for (j in which(colSums(is.na(M)) > 0))
    M[is.na(M[, j]), j] <- round(mean(M[, j], na.rm = TRUE)); M }
  G_eur <- impute_cm(G_eur); G_yri <- impute_cm(G_yri)
  cat(sprintf("  after EUR MAF>=%.2f: m=%d SNPs\n", cfg$maf_lo, length(pos)))

  # -- PolyFun snpvar (off-the-shelf functional prior), matched by position ---
  pq <- file.path(D4, if (as.integer(chr) <= 7L) "snpvar_meta.chr1_7.parquet" else "snpvar_meta.chr8_22.parquet")
  snpvar <- rep(NA_real_, length(pos))
  if (file.exists(pq)) {
    con <- dbConnect(duckdb::duckdb()); on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)
    sv <- dbGetQuery(con, sprintf(
      "SELECT BP AS pos, MAX(snpvar_bin) AS snpvar FROM read_parquet(%s)
       WHERE CHR=%s AND BP>=%d AND BP<%d GROUP BY BP", shQuote(pq), chr, lo, hi))
    snpvar <- sv$snpvar[match(pos, sv$pos)]
    cat(sprintf("  PolyFun snpvar matched: %d/%d SNPs\n", sum(!is.na(snpvar)), length(pos)))
  } else {
    warning("PolyFun parquet not found (", pq, "): polyfun prior falls back to 0.5 everywhere.")
  }

  # -- priors: the real-data fidelity gradient, on the manuscript's
  #    chromosome-level convention p_j in {0.3, 0.5, 0.7} ----------------------
  p_unif <- rep(0.5, length(pos))
  p_poly <- rep(0.5, length(pos))                  # NA snpvar -> neutral 0.5
  ok <- which(!is.na(snpvar))
  if (length(ok) >= 3L) {
    sv <- snpvar[ok]
    if (any(sv == 0) && any(sv > 0)) {
      # UKB snpvar_bin has a point mass at 0 (no annotation signal): plain
      # terciles would put the zeros in the middle bin. Explicit rule instead:
      # zero -> 0.3; positive below/above the median of positives -> 0.5/0.7.
      med <- median(sv[sv > 0])
      p_poly[ok] <- ifelse(sv == 0, 0.3, ifelse(sv <= med, 0.5, 0.7))
    } else {
      q <- quantile(sv, c(1/3, 2/3), names = FALSE)
      p_poly[ok] <- 0.3 + 0.2 * findInterval(sv, q)         # terciles -> 0.3/0.5/0.7
    }
  }
  d      <- abs(pos - tss)
  p_tss  <- ifelse(d < cfg$tss_near, 0.7, ifelse(d < cfg$tss_mid, 0.5, 0.3))
  p_disc <- 1 - p_tss                              # inverted tissue-matched prior

  # -- covariates (FWL block): EUR population dummies + sex; YRI sex only -----
  mk_cov <- function(s) {
    p  <- panel[match(s, panel$sample), ]
    cv <- cbind(sex = as.integer(p$gender == "male"))
    if (any(p$pop %in% EUR_POPS))
      for (pp in setdiff(EUR_POPS, "CEU")) cv <- cbind(cv, as.integer(p$pop == pp))
    colnames(cv) <- c("sex", setdiff(EUR_POPS, "CEU"))[seq_len(ncol(cv))]
    rownames(cv) <- s; cv
  }

  list(gene = strip_ver(expr$TargetID[row]), symbol = expr$Gene_Symbol[row],
       chr = chr, tss = tss, window = cfg$window, pos = pos, maf_eur = maf,
       G_eur = G_eur, G_yri = G_yri,
       y_eur = maybe_int(as.numeric(unlist(expr[row, s_eur], use.names = FALSE))),
       y_yri = maybe_int(as.numeric(unlist(expr[row, s_yri], use.names = FALSE))),
       snpvar = snpvar,
       priors = list(uniform = p_unif, polyfun = p_poly, tss = p_tss, discordant = p_disc),
       covar_eur = mk_cov(s_eur), covar_yri = mk_cov(s_yri),
       samples_eur = s_eur, samples_yri = s_yri, cfg = cfg)
}

# ---- 3. run -----------------------------------------------------------------
if (length(genes) == 1L && nzchar(cfg$out)) {
  saveRDS(build_gene(genes), cfg$out); cat("Saved:", cfg$out, "\n")
} else {
  dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
  for (g in genes) {
    out <- file.path(cfg$out_dir, sprintf("geuvadis_%s.rds", strip_ver(g)))
    if (file.exists(out)) { cat("SKIP (exists):", out, "\n"); next }
    res <- tryCatch(build_gene(g), error = function(e) { message("SKIP ", g, ": ", conditionMessage(e)); NULL })
    if (!is.null(res)) { saveRDS(res, out); cat("Saved:", out, "\n") }
  }
}
