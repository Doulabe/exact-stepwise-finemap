# Exact association-test distributions for prior-informed stepwise SNP prioritization

Reproducibility code for the manuscript *"Exact association-test distributions for
prior-informed stepwise SNP prioritization"*.

The method is a forward stepwise variable-selection framework for fine-mapping.
At each step it scores every candidate SNP, conditional on the SNPs already
selected, by a **local posterior inclusion score** that combines a per-SNP prior
weight with the **exact finite-sample distribution** of a Wald, Score or joint
Multiple-SNP statistic, which accounts for the variability of the
residual-variance estimate that the usual normal and chi-square approximations
treat as fixed. For the joint statistic, the prior-integrated score is computed
exactly, as a convergent mixture of central densities truncated with a
guaranteed error bound (`--multiple_exact TRUE`; a Monte-Carlo average is used
otherwise). Six in-house variants are studied:

| Code | Method | Family |
|------|--------|--------|
| NaW | Naive Wald | single-SNP, normal approx. with fixed variance |
| AsW | Asymptotic Wald | single-SNP, large-sample |
| ExW | Exact Wald | single-SNP, exact non-central *t* |
| AsS | Asymptotic Score | single-SNP, score test |
| NaM | Naive Multiple | joint (Multiple-SNP), non-central chi-square |
| ExM | Exact Multiple | joint, exact non-central *F* |

They are compared with **SuSiE**, **varbvs**, **LASSO**, **CARMA**, **FINEMAP**
and a spike-and-slab MCMC sampler (**BayesC**, via BGLR) on a synthetic
block-AR(1) LD design, a single-locus design, a semi-synthetic design built
on real GEUVADIS cis windows, and a cis-eQTL application to GEUVADIS gene
expression (European discovery, Yoruba replication).

---

## Repository structure

```
.
├── README.md
├── LICENSE                     # MIT
├── install_dependencies.R      # installs all required R packages
├── run_all.sh                  # convenience orchestrator (smoke run)
└── code/                       # all scripts (run FROM this directory)
    ├── Load_packages.R         # package-loading helper
    │
    ├── # --- function libraries (sourced by the scripts below) ---
    ├── Test_statistics_functions_optimized.R          # NaW/AsW/ExW/AsS/NaM/ExM (fixed K), exact prior integrals
    ├── Test_statistics_functions_Criteria_optimized.R # the same under stopping criteria
    ├── External_Comparators.R                         # SuSiE / varbvs / LASSO / CARMA / CARMA_p / FINEMAP / BayesC, metrics
    ├── Generate_XY_realistic.R                        # synthetic AR(1) genotypes and phenotypes
    │
    ├── # --- simulations (produce results_*.rds in the repository root) ---
    ├── Benchmark_External.R          # synthetic AR(1) design (m = 5,000), fixed or adaptive K
    ├── run_synth_n80.sh              # the same design at small n (finite-sample gain)
    ├── Benchmark_Locus_Sim.R         # single-locus design (adds CARMA / CARMA_p / FINEMAP)
    ├── merge_carma_prior_locus.R     # merge the prior-aware CARMA_p arm into results_locus_sim_*.rds
    ├── Robustness_Checks.R           # exact versus Monte-Carlo scores, greedy penalty under extreme LD
    │
    ├── # --- GEUVADIS application and semi-synthetic evaluation ---
    ├── run_geuvadis_pipeline.sh      # end-to-end driver (stages below)
    ├── make_geuvadis_gene_list.R     # the 100 strongest EUR eGenes on chromosomes 1-7
    ├── fetch_geuvadis_windows.sh     # disk-frugal genotype download, cis windows only
    ├── prep_geuvadis_gene.R          # per-gene RDS: EUR/YRI split, normalized expression, four priors
    ├── Benchmark_GEUVADIS.R          # --mode real (EUR discovery, YRI replication) or --mode semisynth
    ├── add_ldaware_replication.R     # LD-aware replication rates (and their Bonferroni variant)
    ├── run_semisynth30.sh            # semi-synthetic evaluation on 30 windows
    │
    ├── # --- re-running a subset of methods after a code change ---
    ├── rerun_from_config.R           # re-run with the configuration stored in a result file
    ├── merge_fixed_methods.R         # merge the recomputed methods into the original result file
    │
    ├── # --- figures and LaTeX tables (read results_*.rds -> Figure_*.pdf, tab_*.tex) ---
    ├── figure_premium_vs_n.R         # gain-vs-n figure (file names keep the premium_vs_n stem)
    ├── regen_sim_tables.R            # regenerate every simulation and GEUVADIS table from result files
    ├── make_scenario_tables.R        # four-prior-scenario tables (AR(1) m = 5,000, locus, stopping rules)
    ├── make_tab_timing.R             # timing-table rows from a GEUVADIS real-mode result
    ├── fill_geuvadis_tables.R        # GEUVADIS tables (fixed K, adaptive, exact-variant, Bonferroni, semi-synthetic, timing)
    ├── compute_pairwise_pi_geuvadis.R  # pairwise probability of superiority across the GEUVADIS genes
    │
    ├── # --- independent validation of the exact prior integrals ---
    ├── exact_multiple_funs.R         # standalone series, Gil-Pelaez inversion and Davies comparison
    └── exact_multiple_validation.R   # accuracy comparison behind the supplement's validation table
```

> **Prior-aware CARMA (`CARMA_p`).** The single-locus simulation also runs
> CARMA with the same per-SNP prior given to `SuSiE_prior` and `varbvs`, passed
> as a continuous `logit(p_j)` functional annotation through CARMA's native
> `w.list`. Because CARMA *learns* the weight of the annotation, it can exploit
> an informative prior and also a reversed (misspecified) one, whereas the
> proposed scores take `p_j` as given.

> **Path convention.** Every script is meant to be run *from the `code/`
> directory*; outputs (`results_*.rds`) and generated tables/figures
> (`tab_*.tex`, `Figure_*.pdf`) are written to the repository root (one level up,
> `../`). These generated artefacts are not tracked (see `.gitignore`); the
> pipeline regenerates them.

---

## Requirements

- **R ≥ 4.0**.
- R packages: `BGLR`, `susieR`, `varbvs`, `glmnet`, `CARMA`, `MASS`, `Matrix`,
  `mvtnorm`, `doParallel`, `foreach`, `dplyr`, `tidyr`, `ggplot2`, `scales`,
  `statmod`, `data.table`, `duckdb`. Install them all with:

  ```bash
  Rscript install_dependencies.R
  ```

- **FINEMAP** is an external command-line binary (not an R package). Install
  FINEMAP v1.3.1 and either put it on your `PATH` or set the environment
  variable `FINEMAP_BIN` to its full path. FINEMAP is only required for the
  single-locus simulation (`Benchmark_Locus_Sim.R`).

- The **GEUVADIS** analyses additionally need `curl`, `bcftools`, `bgzip` and
  `tabix` (htslib) on the `PATH`, and the PolyFun per-SNP `snpvar` table (Parquet
  files, read with `duckdb`) in the directory given by the environment variable
  `POLYFUN_DIR` (default `../polyfun_data`).

### Data

No genotype data are shipped. The synthetic designs are generated on the fly by
`Generate_XY_realistic.R`. The GEUVADIS analyses use public data placed in
`geuvadis_data/` (not tracked): the E-GEUV-1 genotypes, which
`fetch_geuvadis_windows.sh` downloads and reduces to the cis windows, and the
gene quantifications (`GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz`), the
EUR373 best cis-association file (`EUR373.gene.cis.FDR5.best.rs137.txt.gz`) and
the 1000 Genomes sample panel, downloaded from the GEUVADIS and 1000 Genomes
repositories (see the header of `prep_geuvadis_gene.R` for the expected files).

---

## Reproducing the results

All commands below are run from the `code/` directory. Use `--n_cores <k>` to
parallelise. The manuscript uses `B = 50` replications for the AR(1) design,
`B = 200` for the single-locus design and `B = 50` per window for the
semi-synthetic evaluation, with the exact prior integrals for NaM and ExM.

### 1. Synthetic AR(1) design (m = 5,000)

For each prior scenario `{informative, uniform, misspecified, mixed}` and each
`--K_mode {fixed, adaptive}`:

```bash
Rscript Benchmark_External.R --design synthetic --m_synth 5000 --block_size 50 \
    --rho_grid "0.8,0.95,0.98" --n_grid "500,1000" --gamma 0.6 --tau 0.2 \
    --K_mode fixed --scenario informative --B 50 --multiple_exact TRUE \
    --outfile ../results_bench_synth_informative.rds
# ... repeat with --K_mode adaptive (outfile ..._adaptive.rds) and the other 3 scenarios
```

The four-scenario table with all the methods at `n` in `{80, 1000}`:

```bash
Rscript Benchmark_External.R --design synthetic --m_synth 5000 --block_size 50 \
    --n_causal 15 --K_fixed 15 --K_mode fixed --placement A \
    --rho_grid 0.8,0.95,0.98 --n_grid 80,1000 --gamma 0.6 --tau 0.2 \
    --multiple_exact TRUE --B 50 --seed 2026 --skip_varbvs TRUE --skip_bayesc TRUE \
    --scenario informative --outfile ../results_synth_scen_informative.rds
# ... repeat for the other 3 scenarios
```

### 2. Finite-sample gain at small n

```bash
# the m = 5,000 design at n = 80, full gamma x tau grid, 4 scenarios x {fixed, adaptive}
bash run_synth_n80.sh 8 80

# gain-vs-n sweep for the figure (m = 1,000, 20 blocks of 50, B = 20 per n)
Rscript Benchmark_External.R --design synthetic --m_synth 1000 --multiple_exact TRUE \
    --rho_grid 0.95 --n_grid 50,80,150,300,600,1000 --gamma_grid 0.6 --tau_grid 0.2 \
    --scenario informative --K_mode fixed --B 20 \
    --skip_susie TRUE --skip_lasso TRUE --skip_bayesc TRUE --skip_varbvs TRUE \
    --outfile ../results_synth_premium_vs_n.rds
# (the band over rho uses the same command with --rho_grid 0.8,0.98 and
#  --n_grid 50,80,150,300,1000, outfile ../results_synth_premium_vs_n_band.rds)
Rscript figure_premium_vs_n.R     # -> ../Figure_premium_vs_n.pdf
```

### 3. Single-locus design (CARMA / CARMA_p / FINEMAP)

```bash
Rscript Benchmark_Locus_Sim.R --scenario informative --K_mode fixed --B 200 \
    --multiple_exact TRUE
# repeat for the four scenarios x {fixed, adaptive} -> results_locus_sim_<scenario>_<mode>.rds
# the prior-aware CARMA_p arm is included automatically; under the uniform prior
# it falls back to plain CARMA (degenerate annotation).
```

If the `CARMA_p` arm is regenerated on its own, merge it back into the eight
`results_locus_sim_*.rds` with `Rscript merge_carma_prior_locus.R`.

### 4. GEUVADIS application and semi-synthetic evaluation

```bash
POLYFUN_DIR=/path/to/polyfun bash run_geuvadis_pipeline.sh 7
```

The driver is idempotent and runs, in order: the eGene list
(`make_geuvadis_gene_list.R`), the genotype download (`fetch_geuvadis_windows.sh`),
the per-gene RDS (`prep_geuvadis_gene.R`), the real analysis with EUR discovery and
YRI replication (`Benchmark_GEUVADIS.R --mode real`), the LD-aware replication rates
(`add_ldaware_replication.R`), and the semi-synthetic evaluation on 30 windows
(`run_semisynth30.sh`). The tables are then written by

```bash
Rscript fill_geuvadis_tables.R
Rscript compute_pairwise_pi_geuvadis.R --prior tss --out ../tab_pi_geuvadis.tex
```

### 5. Robustness checks

```bash
Rscript Robustness_Checks.R   # exact versus Monte-Carlo scores and greedy penalty under extreme LD
```

### 6. Re-running a subset of methods after a code change

`rerun_from_config.R` re-runs a script with the configuration stored in an
existing result file, for a subset of methods, and `merge_fixed_methods.R` puts
the recomputed rows back into the original file. A canary method recomputed
alongside checks that the data sets are regenerated identically.

```bash
Rscript rerun_from_config.R --orig ../results_bench_synth_informative.rds \
    --script Benchmark_External.R --out ../rerun/results_bench_synth_informative.rds \
    --extra "--inhouse_only NaW,NaM,ExM --multiple_exact TRUE --skip_susie TRUE --skip_lasso TRUE --skip_bayesc TRUE --skip_varbvs TRUE"
Rscript merge_fixed_methods.R --type external --orig ../results_bench_synth_informative.rds \
    --new ../rerun/results_bench_synth_informative.rds --out ../merged/results_bench_synth_informative.rds \
    --canary Naive_Wald
```

---

## Script → manuscript output

| Script | Produces |
|--------|----------|
| `Benchmark_External.R` | synthetic AR(1) results: main simulation tables, four-scenario table, stopping rules |
| `run_synth_n80.sh` + `Benchmark_External.R` (gain sweep) + `figure_premium_vs_n.R` | small-`n` tables and the gain-vs-`n` figure |
| `Benchmark_Locus_Sim.R` (+ `merge_carma_prior_locus.R`) | single-locus tables, including CARMA_p |
| `run_geuvadis_pipeline.sh` + `fill_geuvadis_tables.R` | GEUVADIS replication tables, semi-synthetic table, timing table |
| `compute_pairwise_pi_geuvadis.R` | pairwise probability of superiority across the GEUVADIS genes |
| `Robustness_Checks.R` | exact versus Monte-Carlo scores, greedy-penalty check |

---

## A convenience orchestrator

`run_all.sh` runs a small, fast end-to-end **smoke** version of the pipeline to
confirm that the installation works: an evaluation on a small synthetic locus
(in-house methods + SuSiE/LASSO/BayesC) and the robustness checks. It does not
reproduce the manuscript.

```bash
bash run_all.sh
```

---

## License

Released under the MIT License (see `LICENSE`).

## Citation

If you use this code, please cite the accompanying manuscript. (Citation details
will be added upon publication.)
