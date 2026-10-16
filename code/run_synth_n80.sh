#!/bin/bash
# =============================================================================
# run_synth_n80.sh
# Reproduce the synthetic AR(1) benchmark at SMALL n (default n=80), with the
# EXACT same configuration as the existing n in {500,1000} runs, to expose the
# finite-sample gain of the exact laws (Exact Wald >> Asymptotic/Naive Wald)
# that is negligible at large n.
#
# Config: design=synthetic, m_synth=5000, block_size=50, n_causal=15, K_fixed=15,
#   exact NaM/ExM integrals, B=50, rho={0.8,0.95,0.98}, MAF~U(0.05,0.50), seed=2026,
#   placement=A, n_grid=80. Sweeps the full grid
#   gamma in {0.2,0.6,0.9} x tau in {0.2,0.5,1.0} (the existing synthetic n=500/1000
#   runs used only the central cell gamma=0.6, tau=0.2 -> that cell is the apples-to-
#   apples large-n reference; the other cells map the exact-law gain over power x scale).
#
# Runs all 4 prior scenarios x {fixed, adaptive K} = 8 cells, like the existing
# synthetic supplement. Idempotent: skips a cell whose --outfile already exists.
#
# Usage (from code/):
#   bash run_synth_n80.sh [N_CORES] [N]
#   e.g.  bash run_synth_n80.sh 30 80
# =============================================================================
set -e
cd "$(dirname "$0")"
NC=${1:-30}        # cores
N=${2:-80}         # small sample size

common=(--design synthetic --m_synth 5000 --block_size 50 --n_causal 15 --K_fixed 15
        --rho_grid 0.8,0.95,0.98 --n_grid "$N"
        --gamma_grid 0.2,0.6,0.9 --tau_grid 0.2,0.5,1.0
        --multiple_exact TRUE --B 50 --maf_lo 0.05 --maf_hi 0.5 --maf_distribution uniform
        --seed 2026 --placement A --n_cores "$NC")

run() {                                   # run() <outfile> <extra args...>
  local out="$1"; shift
  if [ -f "$out" ] && [ "$(wc -c < "$out")" -gt 1024 ]; then
    echo "SKIP (exists): $out"; return
  fi
  echo ">>> $out"
  Rscript Benchmark_External.R "${common[@]}" "$@" --outfile "$out"
}

echo "Synthetic AR(1) at n=$N on $NC cores (matched to the n=500/1000 config)."
for scn in informative uniform misspecified mixed; do
  run "../results_bench_synth_${scn}_n${N}.rds"          --scenario "$scn" --K_mode fixed
  run "../results_bench_synth_${scn}_n${N}_adaptive.rds" --scenario "$scn" --K_mode adaptive
done
echo "ALL_N${N}_DONE"
