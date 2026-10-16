#!/bin/zsh
# Post-hoc chain after the exact-df GEUVADIS rerun (results_geuvadis_real_exactdf.rds):
# waits for the RDS, then regenerates every downstream result and table and
# copies them into Article1_SiM. Launched detached; log: geuvadis_data/posthoc_exactdf.log
set -e
cd "$(dirname "$0")"
RDS=../results_geuvadis_real_exactdf.rds

while [ ! -s "$RDS" ]; do sleep 120; done
# let the writer finish flushing
sleep 30

echo "=== $(date) RDS present, starting post-hoc chain ==="
Rscript add_ldaware_replication.R --results $RDS --rds_dir ../geuvadis_rds \
  --out ../results_geuvadis_real_ld.rds
Rscript add_ldaware_replication.R --results $RDS --rds_dir ../geuvadis_rds \
  --bonferroni TRUE --out ../results_geuvadis_real_ld_bonf.rds
Rscript add_ldaware_replication.R --results $RDS --rds_dir ../geuvadis_rds \
  --perm 500 --out ../results_geuvadis_real_ld_perm.rds
echo "=== $(date) LD variants done ==="

Rscript fill_geuvadis_tables.R
Rscript compute_pairwise_pi_geuvadis.R --prior tss --out ../tab_pi_geuvadis.tex
Rscript make_geuvadis_cis.R
Rscript make_tab_geuvadis_r2.R
Rscript make_geuvadis_conditional.R
Rscript figure_budget_curve.R --results ../results_geuvadis_real_ld.rds \
  --rds_dir ../geuvadis_rds --prior tss --out ../Figure_budget_curve.pdf
echo "=== $(date) tables and figures done ==="

for f in tab_geuvadis_fixedK tab_geuvadis_adaptive tab_geuvadis_exact_snp \
         tab_geuvadis_bonferroni tab_geuvadis_perm tab_geuvadis_semisynth \
         tab_geuvadis_ci tab_geuvadis_r2 tab_geuvadis_cond tab_pi_geuvadis; do
  cp ../$f.tex ../../Article1_SiM/
done
cp ../Figure_budget_curve.pdf ../../Article1_SiM/
echo "=== $(date) POSTHOC_EXACTDF_DONE ==="
