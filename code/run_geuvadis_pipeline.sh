#!/bin/bash
# =============================================================================
# run_geuvadis_pipeline.sh
# Unattended end-to-end GEUVADIS chain for the Statistics in Medicine version:
#
#   Stage 1  eGene list + merged cis regions (top-100 FDR5 eGenes chr1-7
#            + ERAP2)                                    [make_geuvadis_gene_list.R]
#   Stage 2  disk-frugal genotype fetch: per chromosome download -> extract
#            the gene windows -> delete the 3-4 GB full VCF
#                                                        [fetch_geuvadis_windows.sh]
#   Stage 3  per-gene RDS (462 samples, EUR/YRI split, INT, 4 priors)
#                                                        [prep_geuvadis_gene.R]
#   Stage 4  real cis-eQTL application, EUR discovery -> YRI replication,
#            then LD-aware replication rates (and their Bonferroni variant)
#                                                        [Benchmark_GEUVADIS.R --mode real,
#                                                         add_ldaware_replication.R]
#   Stage 5  semi-synthetic evaluation on 30 windows, B=50, n in {80,373}
#                                                        [run_semisynth30.sh]
#
# NaM/ExM use the exact prior integrals (--multiple_exact TRUE). POLYFUN_DIR
# must point to the PolyFun snpvar parquet(s) (default ../polyfun_data).
#
# Idempotent: every stage skips outputs that already exist, so the script can
# be re-run after an interruption. All output appended to geuvadis_data/pipeline.log.
#
# Usage (from code/):   bash run_geuvadis_pipeline.sh [N_CORES]
# =============================================================================
set -u
cd "$(dirname "$0")"
export PATH="/opt/homebrew/bin:$PATH"
export POLYFUN_DIR="${POLYFUN_DIR:-../polyfun_data}"
NC=${1:-7}
DG="../geuvadis_data"
VCF_PATTERN="../geuvadis_data/vcf_subsets/GEUVADIS.chr%s.subset.vcf.gz"
LOG="$DG/pipeline.log"
ts() { echo "[$(date '+%F %T')] $*" | tee -a "$LOG"; }
fail() { ts "PIPELINE_FAILED: $*"; exit 1; }

ts "=== GEUVADIS pipeline started (n_cores=$NC) ==="

# ---- Stage 1: eGene list + regions ------------------------------------------
if [ -s "$DG/geuvadis_genes.txt" ] && [ -s "$DG/geuvadis_regions.txt" ]; then
  ts "Stage 1 SKIP: gene list exists ($(wc -l < "$DG/geuvadis_genes.txt" | tr -d ' ') genes)"
else
  ts "Stage 1: building eGene list + regions ..."
  Rscript make_geuvadis_gene_list.R --n_top 100 >> "$LOG" 2>&1 || fail "gene list"
fi

# ---- Stage 2: disk-frugal genotype fetch ------------------------------------
ts "Stage 2: fetching genotype windows (download -> extract -> delete, chr1-7) ..."
bash fetch_geuvadis_windows.sh >> "$LOG" 2>&1 || fail "genotype fetch (see log)"
NS=$(ls "$DG"/vcf_subsets/GEUVADIS.chr*.subset.vcf.gz 2>/dev/null | wc -l | tr -d ' ')
ts "Stage 2 OK: $NS chromosome subsets ($(du -sh "$DG/vcf_subsets" | cut -f1))"

# ---- Stage 3: per-gene RDS ----------------------------------------------------
ts "Stage 3: preparing per-gene RDS (462 samples, INT, 4 priors) ..."
Rscript prep_geuvadis_gene.R --genes_file "$DG/geuvadis_genes.txt" \
    --vcf_pattern "$VCF_PATTERN" \
    --panel_file "$DG/geuvadis_panel.txt" \
    --out_dir ../geuvadis_rds >> "$LOG" 2>&1 || fail "prep"
NR=$(ls ../geuvadis_rds/geuvadis_*.rds 2>/dev/null | wc -l | tr -d ' ')
[ "$NR" -ge 50 ] || fail "only $NR gene RDS produced"
ts "Stage 3 OK: $NR gene RDS in geuvadis_rds/"

# ---- Stage 4: real cis-eQTL benchmark ---------------------------------------
if [ -s ../results_geuvadis_real.rds ]; then
  ts "Stage 4 SKIP: results_geuvadis_real.rds exists"
else
  ts "Stage 4: real benchmark (EUR discovery -> YRI replication, $NR genes) ..."
  Rscript Benchmark_GEUVADIS.R --rds_dir ../geuvadis_rds --mode real \
      --adaptive TRUE --K_fixed 5 --multiple_exact TRUE --n_cores "$NC" \
      --outfile ../results_geuvadis_real.rds >> "$LOG" 2>&1 || fail "real benchmark"
  ts "Stage 4 OK"
fi
if [ -s ../results_geuvadis_real_ld_bonf.rds ]; then
  ts "Stage 4b SKIP: LD-aware replication results exist"
else
  ts "Stage 4b: LD-aware replication rates ..."
  Rscript add_ldaware_replication.R --results ../results_geuvadis_real.rds \
      --rds_dir ../geuvadis_rds --out ../results_geuvadis_real_ld.rds >> "$LOG" 2>&1 \
    && Rscript add_ldaware_replication.R --results ../results_geuvadis_real.rds \
      --rds_dir ../geuvadis_rds --bonferroni TRUE \
      --out ../results_geuvadis_real_ld_bonf.rds >> "$LOG" 2>&1 || fail "LD-aware replication"
  ts "Stage 4b OK"
fi

# ---- Stage 5: semi-synthetic evaluation (30 windows) --------------------------
ts "Stage 5: semi-synthetic evaluation on 30 windows (B=50, n in {80,373}) ..."
bash run_semisynth30.sh "$NC" || fail "semi-synthetic evaluation (see $DG/semisynth30.log)"
ts "Stage 5 OK"

ts "=== PIPELINE_ALL_DONE ==="
