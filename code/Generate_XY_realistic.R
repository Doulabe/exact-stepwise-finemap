# ==============================================================================
# File: Generate_XY_realistic.R
# Project: exact-stepwise-finemap (Exact association-test distributions for prior-informed stepwise SNP prioritization)
# Author: N.K. Doulabe 
#
# Description:
#   Synthetic GWAS-like data generator for the controlled simulations. It
#   gives parametric control over the LD strength (rho) and the MAF
#   distribution, which the evaluations on real human LD (GEUVADIS windows)
#   complement.
#
#   Key features:
#     (1) m SNPs partitioned into LD blocks of fixed size (default 50 SNPs);
#     (2) Within each block, the SNP correlation follows an AR(1) structure:
#             R[a, b] = rho^|a - b|                  for a, b in the same block
#             R[a, b] = rho_inter (~ 0)              for different blocks
#     (3) Minor allele frequencies are independent Uniform(0.05, 0.50) by
#         default, mimicking the realistic spread of real GWAS data;
#     (4) Three causal-placement plans (A / B / random).
#
# Usage:
#   X <- Generate.X.realistic(n = 1000, m = 5000)
#   X <- Generate.X.realistic(n = 500,  m = 5000, rho_intra = 0.95)
# ==============================================================================

suppressPackageStartupMessages({
  library(mvtnorm)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# 1. MAF sampling
# ------------------------------------------------------------------------------
sample_maf <- function(m, maf_lo = 0.05, maf_hi = 0.50,
                       distribution = c("uniform", "fixed")) {
  distribution <- match.arg(distribution)
  if (distribution == "fixed") return(rep(maf_lo, m))
  stopifnot(maf_lo > 0, maf_hi < 1, maf_lo < maf_hi)
  runif(m, min = maf_lo, max = maf_hi)
}

# ------------------------------------------------------------------------------
# 2. AR(1) correlation block
# ------------------------------------------------------------------------------
make_ar1_block <- function(size, rho) {
  stopifnot(abs(rho) < 1, size >= 1)
  idx <- 0:(size - 1)
  rho ^ abs(outer(idx, idx, "-"))
}

# ------------------------------------------------------------------------------
# 3. Block-diagonal AR(1) correlation matrix
# ------------------------------------------------------------------------------
make_corr_blocks_ar1 <- function(m, block_size = 50L,
                                  rho_intra = 0.9, rho_inter = 0.0) {
  stopifnot(abs(rho_intra) < 1, abs(rho_inter) < 1, block_size >= 1)
  block_id <- rep(seq_len(ceiling(m / block_size)),
                  each = block_size, length.out = m)
  R <- matrix(rho_inter, nrow = m, ncol = m)
  for (b in unique(block_id)) {
    idx <- which(block_id == b)
    R[idx, idx] <- make_ar1_block(length(idx), rho_intra)
  }
  diag(R) <- 1
  R
}

# ------------------------------------------------------------------------------
# 4. Compound-symmetry baseline (weak-LD reference)
# ------------------------------------------------------------------------------
make_corr_compound <- function(m, rho) {
  stopifnot(abs(rho) < 1)
  R <- matrix(rho, nrow = m, ncol = m); diag(R) <- 1
  R
}

build_correlation <- function(m,
                              ld_structure = c("blocks", "compound"),
                              rho_intra = 0.9,
                              rho_inter = 0.0,
                              block_size = 50L) {
  ld_structure <- match.arg(ld_structure)
  if (ld_structure == "compound") return(make_corr_compound(m, rho_intra))
  make_corr_blocks_ar1(m, block_size = block_size,
                       rho_intra = rho_intra, rho_inter = rho_inter)
}

# ------------------------------------------------------------------------------
# 4b. Pre-warmed square-root factors of R (speed only, identical draws)
# ------------------------------------------------------------------------------
# The cache lives in globalenv() because every foreach task re-sources this
# file with local = TRUE: an existing cache must survive, and forked workers
# read the factors pre-warmed by the parent copy-on-write. Access it only via
# its name as a string (get0, assign), never as a variable: foreach exports the
# globals of exported functions and would serialise the factors to every worker.
if (!exists(".ar1_factor_cache", envir = globalenv(), inherits = FALSE)) {
  assign(".ar1_factor_cache", new.env(parent = emptyenv()), envir = globalenv())
}

# %.17g makes the key identify the doubles, hence the matrix R, exactly.
ar1_factor_key <- function(m, ld_structure, rho_intra, rho_inter, block_size) {
  paste(c(ld_structure,
          sprintf("%.17g", as.double(c(m, rho_intra, rho_inter, block_size)))),
        collapse = "|")
}

# These two functions must reproduce mvtnorm::rmvnorm(n, sigma = R) (method
# "eigen", mvtnorm 1.3.6) bit for bit: rows simulated through the cache are
# merged with rows simulated by the uncached path.
ar1_sqrt_factor <- function(sigma) {
  ev <- eigen(sigma, symmetric = TRUE)
  if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
    warning("sigma is numerically not positive semidefinite")
  }
  t(ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
}

ar1_rmvnorm_factor <- function(n, fac) {
  mean <- rep(0, nrow(fac))
  # byrow = TRUE is rmvnorm's default (byrow = !pre0.9_9994).
  retval <- matrix(stats::rnorm(n * ncol(fac)), nrow = n, byrow = TRUE) %*% fac
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}

# NULL unless the parent pre-warmed this key. Workers never fill the cache, so
# a process without pre-warm keeps the original path and its memory use.
ar1_cached_factor <- function(m, ld_structure, rho_intra, rho_inter, block_size) {
  cache <- get0(".ar1_factor_cache", envir = globalenv(), inherits = FALSE)
  if (!is.environment(cache)) return(NULL)
  get0(ar1_factor_key(m, ld_structure, rho_intra, rho_inter, block_size),
       envir = cache, inherits = FALSE)
}

# Call in the parent before makeForkCluster(). The caller's RNG state is left
# untouched, and nothing is cached unless the self-check passes.
ar1_prewarm_factors <- function(m, rho_values, rho_inter, block_size,
                                ld_structure = "blocks", cores = 1L) {
  t0 <- Sys.time()
  elapsed <- function() as.numeric(difftime(Sys.time(), t0, units = "secs"))
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })

  S  <- build_correlation(60L, "blocks", rho_intra = 0.9, rho_inter = 0,
                          block_size = 20L)
  ok <- tryCatch({
    set.seed(1L); a <- mvtnorm::rmvnorm(7L, sigma = S)
    set.seed(1L); b <- ar1_rmvnorm_factor(7L, ar1_sqrt_factor(S))
    identical(dim(a), dim(b)) && isTRUE(all(a == b))
  }, error = function(e) FALSE)
  if (!ok) {
    cat(sprintf("AR(1) factor cache disabled: self-check against mvtnorm::rmvnorm failed (%.1f s).\n",
                elapsed()))
    return(invisible(FALSE))
  }

  cache <- get(".ar1_factor_cache", envir = genv, inherits = FALSE)
  keys  <- vapply(rho_values, function(r)
    ar1_factor_key(m, ld_structure, r, rho_inter, block_size), "")
  todo  <- which(!duplicated(keys) &
                   !vapply(keys, exists, NA, envir = cache, inherits = FALSE))
  # Optional disk cache shared by successive runs on the same machine
  # (AR1_FACTOR_DIR): a factor saved by a previous run is the same eigen()
  # result, stored exactly by saveRDS, so reloading it keeps the draws bit for bit.
  fdir  <- Sys.getenv("AR1_FACTOR_DIR")
  fpath <- function(key) file.path(fdir, paste0("ar1_factor_", gsub("[^A-Za-z0-9.]", "_", key), ".rds"))
  if (nzchar(fdir) && length(todo)) {
    for (i in todo) {
      if (!file.exists(fpath(keys[i]))) next
      f <- tryCatch(readRDS(fpath(keys[i])), error = function(e) NULL)
      if (is.matrix(f) && is.double(f) && all(dim(f) == m)) assign(keys[i], f, envir = cache)
    }
    todo <- todo[!vapply(keys[todo], exists, NA, envir = cache, inherits = FALSE)]
  }
  facs  <- parallel::mclapply(rho_values[todo], function(r)
    ar1_sqrt_factor(build_correlation(m, ld_structure = ld_structure,
                                      rho_intra = r, rho_inter = rho_inter,
                                      block_size = block_size)),
    mc.cores = max(1L, min(length(todo), as.integer(cores))),
    mc.set.seed = FALSE)
  # A failed or killed child leaves its rho uncached (rmvnorm path).
  for (i in seq_along(todo)) {
    f <- facs[[i]]
    if (is.matrix(f) && is.double(f) && all(dim(f) == m)) {
      assign(keys[todo[i]], f, envir = cache)
      if (nzchar(fdir) && dir.exists(fdir)) {
        tmp <- paste0(fpath(keys[todo[i]]), ".tmp", Sys.getpid())
        saveRDS(f, tmp, compress = FALSE)
        file.rename(tmp, fpath(keys[todo[i]]))
      }
    }
  }
  invisible(gc(verbose = FALSE))
  n_ok <- sum(vapply(unique(keys), exists, NA, envir = cache, inherits = FALSE))
  cat(sprintf("AR(1) factor cache %s: %d/%d factor(s) for m = %s pre-warmed in %.1f s.\n",
              if (n_ok > 0L) "active" else "disabled", n_ok,
              length(unique(keys)), format(m), elapsed()))
  invisible(n_ok > 0L)
}

# ------------------------------------------------------------------------------
# 5. Realistic genotype matrix generator
# ------------------------------------------------------------------------------
Generate.X.realistic <- function(n = 1000L, m = 5000L,
                                 ld_structure = c("blocks", "compound"),
                                 rho_intra = 0.9,
                                 rho_inter = 0.0,
                                 block_size = 50L,
                                 maf_lo = 0.05,
                                 maf_hi = 0.50,
                                 maf_distribution = c("uniform", "fixed"),
                                 return_maf = FALSE) {
  ld_structure <- match.arg(ld_structure)
  maf_distribution <- match.arg(maf_distribution)

  fac <- ar1_cached_factor(m, ld_structure, rho_intra, rho_inter, block_size)
  if (is.null(fac)) {
    R   <- build_correlation(m, ld_structure = ld_structure,
                             rho_intra = rho_intra, rho_inter = rho_inter,
                             block_size = block_size)
    maf <- sample_maf(m, maf_lo, maf_hi, maf_distribution)

    Z1 <- mvtnorm::rmvnorm(n, sigma = R)
    Z2 <- mvtnorm::rmvnorm(n, sigma = R)
  } else {
    # Same RNG draws, in the same order, as the branch above.
    maf <- sample_maf(m, maf_lo, maf_hi, maf_distribution)
    Z1  <- ar1_rmvnorm_factor(n, fac)
    Z2  <- ar1_rmvnorm_factor(n, fac)
  }
  U1 <- pnorm(Z1); U2 <- pnorm(Z2)
  maf_mat <- matrix(maf, nrow = n, ncol = m, byrow = TRUE)
  G  <- (U1 < maf_mat) * 1L + (U2 < maf_mat) * 1L

  X <- scale(G); X[is.nan(X)] <- 0
  if (return_maf) return(list(X = X, maf = maf))
  X
}

# ------------------------------------------------------------------------------
# 6. Causal indices
# ------------------------------------------------------------------------------
choose_causal_indices <- function(m, n_causal = 15L, block_size = 50L,
                                  placement = c("A", "B", "random"),
                                  causals_per_block = 3L,
                                  seed = 1L) {
  placement <- match.arg(placement)
  set.seed(seed)
  n_blocks <- ceiling(m / block_size)

  if (placement == "A") {
    stopifnot(n_causal <= n_blocks)
    block_offsets <- (seq_len(n_causal) - 1L) * block_size
    centre <- block_size %/% 2L
    return(as.integer(block_offsets + centre))
  }

  if (placement == "B") {
    k_blocks <- ceiling(n_causal / causals_per_block)
    stopifnot(k_blocks <= n_blocks)
    chosen_blocks <- sample.int(n_blocks, k_blocks)
    idx <- integer(0)
    for (b in chosen_blocks) {
      n_pick <- min(causals_per_block, n_causal - length(idx))
      pos_in_block <- sample.int(block_size, n_pick)
      idx <- c(idx, (b - 1L) * block_size + pos_in_block)
    }
    return(sort(as.integer(idx))[seq_len(n_causal)])
  }

  # placement == "random"
  blocks_pick <- sample.int(n_blocks, n_causal, replace = (n_causal > n_blocks))
  centre <- block_size %/% 2L
  jitter <- sample.int(block_size, n_causal, replace = TRUE) - centre
  sort(as.integer((blocks_pick - 1L) * block_size + centre + jitter))
}

# ------------------------------------------------------------------------------
# 7. Linear-model phenotype simulator
# ------------------------------------------------------------------------------
Generate.Y <- function(n, X, beta, sigma) {
  m <- ncol(X)
  if (length(beta) == m + 1L) {
    Y <- beta[1L] + X %*% beta[-1L] + rnorm(n, sd = sigma)
  } else if (length(beta) == m) {
    Y <- X %*% beta + rnorm(n, sd = sigma)
  } else {
    stop("beta must have length m or m+1.")
  }
  as.numeric(Y)
}

# ------------------------------------------------------------------------------
# 8. Smoke test
# ------------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  set.seed(2026)
  cat("--- Generate.X.realistic n=200, m=500, AR(1) rho=0.9 ---\n")
  out <- Generate.X.realistic(n = 200L, m = 500L,
                              block_size = 50L, rho_intra = 0.9,
                              return_maf = TRUE)
  cat("MAF summary:\n"); print(summary(out$maf))
  cat(sprintf("Empirical correlation |a-b|=1: %.3f (theoretical %.2f)\n",
              cor(out$X[,1], out$X[,2]), 0.9))
  cat(sprintf("Empirical correlation between blocks: %.3f (theoretical 0)\n",
              cor(out$X[,1], out$X[,100])))
}
