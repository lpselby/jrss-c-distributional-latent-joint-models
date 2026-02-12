---
title: "nhanes_pipeline_and_run.R"
output: html_document
---



# nhanes_pipeline_and_run.R
library(tidyverse)
library(haven)   # for SAS/XPT reading if needed
library(survey)
library(rstan)

# Load NHANES data (replace with your local path)
# nhanes_raw <- read_xpt("DEMO_J.XPT") # example
# or read csv you prepared
# df <- read_csv("nhanes_prepared.csv")
df <- read.csv("nhanes_synthetic.csv", stringsAsFactors = FALSE)
# map columns expected by Stan script:
df$yb <- df$diabetes_self
df$yc <- df$fasting_glucose
# create model matrix/etc as in your earlier script


# Example placeholder: prepare df with columns
#yb: self_report_diabetes (0/1), yc: fasting_glucose (continuous)
# covariates: age, sex (0/1), race_eth, education, income_pov, bmi, sbp, dbp, wbc, smoking, physical_activity
# weights: wt_mec, strata, psu

# PREPROCESS - standardize continuous covariates
cont_vars <- c("age","bmi","sbp","dbp","wbc")
df <- df %>% mutate(across(all_of(cont_vars), ~ scale(.)[,1]))

#List variables needed
vars <- c("age", "sex", "race", "education", "income_pov", "bmi", "sbp", "dbp", "wbc", "smoker", "phys_act", "yb", "yc")
#keep only complete cases
df <- df[complete.cases(df[ , vars]), ]

# Build design matrices
X <- model.matrix(~ age + sex + race + education + income_pov + bmi + sbp + dbp + wbc + smoker + phys_act, data=df)

stopifnot(nrow(X)==nrow(df))


# Create Stan data list as in simulation, but include weights by scaling log-likelihood in Stan or use pseudo-likelihood
# For weighted Stan: pass weights and multiply log_lik contributions in model block; simpler: use replicate weights or resampling outside Stan.

# Example: assemble data list for Stan:
idx <- which(!is.na(df$yc))
df2 <- df[idx, , drop = FALSE]
X2 <- model.matrix(~ age + sex + race + education + income_pov + bmi + sbp + dbp + wbc + smoker + phys_act, data=df2)
stan_data <- list(
  N = nrow(df2),
  yb = df2$yb,
  yc = df2$yc,
  p_mu = ncol(X2),
  p_sigma = ncol(X2),
  p_alpha = ncol(X2),
  p_kappa = ncol(X2),
  X_mu = X2,
  X_sigma = X2,
  X_alpha = X2,
  X_kappa = X2,
  sens_prior_a = 18, sens_prior_b = 3,
  spec_prior_a = 20, spec_prior_b = 2
)

fit_nhanes <- stan(file="stan_pmm.stan", data=stan_data, iter=4000, warmup=2000, chains=4, control=list(adapt_delta=0.95, max_treedepth=12))
saveRDS(fit_nhanes, file="nhanes_pmm_fit.rds")

# postprocess: extract betas, compute latent moments for each individual, plot Figure 2 heatmap, Figure 3 densities, compute AUC etc.







