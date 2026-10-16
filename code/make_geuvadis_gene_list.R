# =============================================================================
# File: make_geuvadis_gene_list.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe
#
# Build the eGene list for the GEUVADIS cis-eQTL application from the official
# GEUVADIS EUR373 best-association-per-gene file (FDR 5%): the --n_top
# strongest eGenes (by -log10 p) on chromosomes 1-7 (the chromosomes covered
# by the locally available PolyFun snpvar parquet), with ERAP2 forced in.
# Also writes the small subset used for the semi-synthetic benchmark
# (top eGene + a mid-list eGene + ERAP2).
#
# Usage (from code/):
#   Rscript make_geuvadis_gene_list.R --n_top 100
# =============================================================================
cfg <- list(
  fdr5          = "../geuvadis_data/EUR373.gene.cis.FDR5.best.rs137.txt.gz",
  n_top         = 100L,
  chrs          = "1,2,3,4,5,6,7",
  margin        = 1100000L,            # regions half-window (> prep's 1 Mb, safety margin)
  out           = "../geuvadis_data/geuvadis_genes.txt",
  out_semisynth = "../geuvadis_data/geuvadis_genes_semisynth.txt",
  out_regions   = "../geuvadis_data/geuvadis_regions.txt"
)
a <- commandArgs(trailingOnly = TRUE)
for (k in names(cfg)) {
  i <- match(paste0("--", k), a)
  if (!is.na(i)) cfg[[k]] <- methods::as(a[i + 1L], class(cfg[[k]]))
}
ERAP2 <- "ENSG00000164308"

# Columns (no header): 1 snp_id, 3 gene_id, 5 chr_snp, 6 chr_gene, 7 snp_pos,
# 8 gene_tss, 9 distance, 10 slope, 11 pvalue, 12 -log10(pvalue)
d <- read.delim(gzfile(cfg$fdr5), header = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(d) >= 12L)
d$gene   <- sub("\\..*$", "", d$V3)
d$log10p <- as.numeric(d$V12)
keep_chr <- trimws(strsplit(cfg$chrs, ",", fixed = TRUE)[[1L]])
d <- d[as.character(d$V6) %in% keep_chr & !duplicated(d$gene), ]
d <- d[order(-d$log10p), ]
cat(sprintf("eGenes at FDR5 on chr {%s}: %d\n", cfg$chrs, nrow(d)))

genes <- head(d$gene, cfg$n_top)
if (!ERAP2 %in% genes) genes <- c(genes, ERAP2)
writeLines(genes, cfg$out)
cat(sprintf("Wrote %d genes to %s (top -log10p: %.1f ... %.1f)\n",
            length(genes), cfg$out, d$log10p[1], d$log10p[min(cfg$n_top, nrow(d))]))

semis <- unique(c(genes[1], genes[ceiling(length(genes) / 2)], ERAP2))
writeLines(semis, cfg$out_semisynth)
cat("Semisynth subset:", paste(semis, collapse = ", "), "\n")

# Merged cis regions (chr <TAB> start <TAB> end) for the disk-frugal per-chr
# VCF window extraction (fetch_geuvadis_windows.sh / bcftools view -R).
sel <- d[d$gene %in% genes, c("V6", "V8")]
if (!all(genes %in% d$gene))
  warning("genes without an FDR5 row (no region emitted): ",
          paste(setdiff(genes, d$gene), collapse = ", "))
reg <- data.frame(chr = as.character(sel$V6),
                  start = pmax(1L, as.integer(sel$V8) - cfg$margin),
                  end   = as.integer(sel$V8) + cfg$margin)
reg <- reg[order(reg$chr, reg$start), ]
merged <- list()
for (cc in unique(reg$chr)) {
  r <- reg[reg$chr == cc, ]
  cur <- r[1, ]
  for (j in seq_len(nrow(r))[-1]) {
    if (r$start[j] <= cur$end) cur$end <- max(cur$end, r$end[j])
    else { merged[[length(merged) + 1L]] <- cur; cur <- r[j, ] }
  }
  merged[[length(merged) + 1L]] <- cur
}
merged <- do.call(rbind, merged)
write.table(merged, cfg$out_regions, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
cat(sprintf("Wrote %d merged regions (%.1f Mb total) to %s\n",
            nrow(merged), sum(merged$end - merged$start) / 1e6, cfg$out_regions))
