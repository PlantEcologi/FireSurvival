# FireSurvival Stan Reimplementation

This repository now includes a Stan-based reimplementation of the non-spatial wildfire survival model from Wilson et al. (2010).

## Model overview

The Stan model is in:

- `stan/firemodel.stan`

It fits a hierarchical-ready probit regression for binary fire occurrence:

- \( y_n \sim \text{Bernoulli}(\Phi(\eta_n)) \)
- \( \eta_n = X_n \beta \)

where \(\Phi\) is the standard normal CDF.

### Priors

To match the original R/MCMC implementation, coefficients use diffuse priors:

- \( \beta_k \sim \mathcal{N}(0, 1000) \) in Stan parameterization (SD = 1000)

### Left censoring and legacy exclusions

The legacy non-spatial code excludes:

1. Pre-1980 rows (first 112 seasons per grid cell), and
2. The season immediately after each observed fire.

These exclusions are passed as `include_row` and applied directly in the Stan likelihood.

## Workflow

### 1) Prepare data

Script:

- `R/01_prepare_data.R`

What it does:

- Downloads `conc.csv`, `data.csv`, and `static.csv` from the `AdamWilsonLab/FireSurvival2010` release tag `Data` (if missing)
- Merges covariates (`data.csv` + `static.csv` + `conc.csv`)
- Scales `pptconc` (as in legacy script)
- Fills missing `aao` with the mean (legacy behavior)
- Builds seasonal indicators (spring/summer/fall)
- Builds Stan-ready data list and saves to:
  - `data/processed/fire_stan_data.rds`

Run:

```r
Rscript R/01_prepare_data.R
```

### 2) Fit Stan model

Script:

- `R/02_fit_model.R`

What it does:

- Loads prepared RDS data
- Compiles Stan model
- Samples with 3 chains (`iter_warmup = 2000`, `iter_sampling = 10000` by default)
- Uses `cmdstanr` if available, otherwise `rstan`
- Saves fitted object to:
  - `results/model/fire_stan_fit.rds`

Run:

```r
Rscript R/02_fit_model.R
```

### 3) Examine output

Script:

- `R/03_examine_output.R`

What it does:

- Loads fitted model
- Creates trace plots, forest plot, and posterior histograms
- Creates posterior coefficient summary (mean, sd, 2.5%, 50%, 97.5%)
- Compares coefficient order with legacy output ordering
- Summarizes acceptance statistics and Rhat/ESS diagnostics
- Writes outputs under:
  - `results/summary/`

Run:

```r
Rscript R/03_examine_output.R
```

## Output files

Key outputs include:

- `results/summary/beta_summary.csv`
- `results/summary/coefficient_order_comparison.csv`
- `results/summary/diagnostics_rhat_ess.csv`
- `results/summary/acceptance_summary.csv`
- `results/summary/traceplots_beta.pdf`
- `results/summary/forestplot_beta.pdf`
- `results/summary/posterior_hist_beta.pdf`

## Differences and improvements vs original R sampler

- Uses Stan HMC/NUTS instead of custom random-walk Metropolis updates.
- Provides stronger diagnostics by default (Rhat, ESS, sampler diagnostics).
- Keeps model structure simple and extensible for future additions (e.g., spatial or ecoregion random effects).
- Improves reproducibility with explicit prepare/fit/examine scripts.

## Notes on computational efficiency

- CmdStan (`cmdstanr`) is preferred for compilation and runtime performance.
- The full default run (`3` chains, `2000` warmup, `10000` sampling) can be computationally heavy; reduce iterations for quick checks.
