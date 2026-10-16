#!/bin/bash
# =============================================================================
# run_semisynth30.sh
# Semi-synthetic evaluation on real GEUVADIS cis windows: 30 windows spanning
# the eGene list (the genes ranked round(seq(1, 100, length.out = 30)) in
# geuvadis_genes.txt, written by make_geuvadis_gene_list.R), B = 50 replications
# per window, n in {80, 373}, exact NaM/ExM prior integrals.
# Idempotent: a window whose result file exists is skipped.
#
# Usage (from code/):   bash run_semisynth30.sh [N_CORES]
# =============================================================================
set -u
cd "$(dirname "$0")"
NC=${1:-7}
DG=../geuvadis_data
LOG=$DG/semisynth30.log
GENES=$DG/geuvadis_genes_semisynth30.txt
: > "$LOG"
if [ ! -s "$GENES" ]; then
  [ -s "$DG/geuvadis_genes.txt" ] || { echo "missing $DG/geuvadis_genes.txt (run make_geuvadis_gene_list.R)"; exit 1; }
  Rscript -e 'g <- readLines("../geuvadis_data/geuvadis_genes.txt"); writeLines(g[round(seq(1, length(g), length.out = 30))], "../geuvadis_data/geuvadis_genes_semisynth30.txt")'
fi
while read -r g; do
  [ -z "$g" ] && continue
  out="../results_geuvadis_semisynth50_${g}.rds"
  if [ -s "$out" ]; then echo "SKIP $g" >> "$LOG"; continue; fi
  [ -s "../geuvadis_rds/geuvadis_${g}.rds" ] || { echo "WARN no RDS $g" >> "$LOG"; continue; }
  echo "[$(date '+%F %T')] >>> $g" >> "$LOG"
  Rscript Benchmark_GEUVADIS.R --rds "../geuvadis_rds/geuvadis_${g}.rds" \
      --mode semisynth --B 50 --n_grid 80,373 --K_true 5 --gamma 0.6 \
      --multiple_exact TRUE --n_cores "$NC" --outfile "$out" >> "$LOG" 2>&1 || { echo "FAIL $g" >> "$LOG"; exit 1; }
done < "$GENES"
echo "[$(date '+%F %T')] SEMISYNTH30_ALL_DONE" >> "$LOG"
