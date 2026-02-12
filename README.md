# Reproducing the Analyses

## Required R Packages

Install the required packages:

```r
install.packages(c("rstan", "mgcv", "tidyverse", "pROC", "fields", "boot"))
```

## 1. Fit the SCMM Surface

```bash
Rscript scmm_surface_fit.R
```

This produces:

- `scmm_surface_gam.rds`
- `scmm_grid_pi.csv`

## 2. Run the Simulation Study

```bash
Rscript simulate_and_run.R
```

This produces:

- `sim_results_*.rds` for each scenario

## 3. NHANES Analysis

Prepare the processed NHANES dataset (no raw PHI included in this repository), then run:

```bash
Rscript nhanes_pipeline_and_run.R
```

## Output

All results, figures, and tables are saved in the working directory.

---

Synthetic NHANES-style data are provided to allow reproducible testing without access to restricted data.
