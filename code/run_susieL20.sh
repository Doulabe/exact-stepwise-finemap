#!/bin/zsh
# SuSiE L=20 rerun (reviewer M7: L=10 default mechanically caps recall at 10/15
# in the K=15 synthetic designs). Reruns the SuSiE rows only, then merges them
# into copies of the results_exact files under results_exact_v2/.
set -e
cd "$(dirname "$0")"
SC=/private/tmp/claude-501/-Users-kossi-Desktop-Dossier-these-Redaction-Projet2-Summary-data-Claude-FineMap/7ef7ba61-a826-450e-a31f-126fccac561c/scratchpad
RE=../../server_rerun_fix_20260928/results_exact
V2=../../server_rerun_fix_20260928/results_exact_v2
mkdir -p $SC/susieL20 $V2
EXTRA="--inhouse_only none --skip_lasso TRUE --skip_varbvs TRUE --skip_bayesc TRUE --skip_susie FALSE --susie_L 20"

for f in results_bench_synth_informative_adaptive.rds \
         results_bench_synth_informative.rds \
         results_synth_scen_informative.rds results_synth_scen_mixed.rds \
         results_synth_scen_uniform.rds results_synth_scen_misspecified.rds \
         results_bench_synth_informative_n80.rds; do
  echo "=== $(date) rerun SuSiE L=20: $f ==="
  AR1_FACTOR_DIR=$SC/ar1cache Rscript rerun_from_config.R --orig $RE/$f \
    --script Benchmark_External.R --out $SC/susieL20/$f --n_cores 5 --extra "$EXTRA"
  Rscript merge_fixed_methods.R --type external --orig $RE/$f \
    --new $SC/susieL20/$f --out $V2/$f --canary FALSE --keep_runtime FALSE
done
echo "=== $(date) SUSIEL20_DONE ==="
