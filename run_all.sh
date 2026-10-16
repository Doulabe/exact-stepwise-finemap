#!/usr/bin/env bash
# =============================================================================
# run_all.sh  --  convenience SMOKE run of the pipeline.
#
# Runs a small, fast end-to-end subset to confirm that the installation works
# and the core scripts execute without error. It does NOT reproduce the full
# manuscript (that uses m = 5,000, four prior scenarios, the GEUVADIS data and
# larger B); see README.md for the complete, per-stage reproduction commands.
#
# Usage:  bash run_all.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/code"

echo "==> [1/2] Evaluation on a small synthetic locus (in-house + SuSiE/LASSO/BayesC) ..."
Rscript Benchmark_External.R --design synthetic --m_synth 300 --rho_grid 0.9 \
    --gamma_grid 0.6 --tau_grid 0.2 --n_grid 500 --K_mode fixed --n_causal 5 \
    --K_fixed 5 --scenario informative --B 3 --skip_varbvs TRUE \
    --multiple_exact TRUE --outfile ../results_bench_synth_smoke.rds

echo "==> [2/2] Robustness checks (exact versus Monte-Carlo scores, greedy penalty) ..."
Rscript Robustness_Checks.R

echo
echo "Smoke run complete. Outputs were written to the repository root (../)."
echo "For the full results, run the per-stage commands documented in README.md."
