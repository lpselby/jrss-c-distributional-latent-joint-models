---
title: "simulate_and_run.R"
output: html_document
---

# simulate_and_run.R
library(tidyverse)

library(rstan)
library(pROC)
library(boot)
library(mgcv)  # for SCMM surface predict


# set rstan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# paths
stan_file <- "stan_pmm.stan"
scmm_rds <- "scmm_surface_gam.rds"

# load SCMM surface
scmm_fit <- readRDS(scmm_rds)

# functions to simulate Z under different scenarios
simulate_Z_sas <- function(mu, logsig2, alpha, kappa, n) {
  sigma <- sqrt(exp(logsig2))
  U <- rnorm(n)
  Z <- mu + sigma * sinh( (asinh(U) + alpha) / kappa )
  Z
}
simulate_Z_skewt <- function(mu, logsig2, alpha, kappa, n, df=5) {
  # approximate skew-t via location-scale transformation with pearson matching, simpler: use rt
  sigma <- sqrt(exp(logsig2))
  T <- rt(n, df=df)
  # center/scale to match variance 1 -> variance of rt is df/(df-2) for df>2
  Tstd <- T / sqrt(df/(df-2))
  Z <- mu + sigma * (Tstd + alpha) # crude add skew via alpha shift
  Z
}
simulate_Z_gaussian_re <- function(mu, logsig2, alpha, kappa, n, rho=0.6) {
  # simulate shared RE scenario: here Z ~ N(mu, sigma^2) independently, but introduce correlated RE with separate D?
  sigma <- sqrt(exp(logsig2))
  rnorm(n, mean=mu, sd=sigma)
}

# evaluation helper functions
compute_metrics <- function(yb, p_hat, yc, yc_hat, yc_samples) {
  # AUC
  auc <- tryCatch(as.numeric(roc(yb, p_hat)$auc), error=function(e) NA)
  brier <- mean((p_hat - yb)^2)
  calib_slope <- tryCatch(coef(lm(yb ~ p_hat))[2], error=function(e) NA)
  rmse <- sqrt(mean((yc - yc_hat)^2))
  # coverage
  # yc_samples: matrix n x nsamples of predictive draws
  lower <- apply(yc_samples, 1, quantile, 0.025)
  upper <- apply(yc_samples, 1, quantile, 0.975)
  coverage <- mean((yc >= lower) & (yc <= upper))
  list(AUC=auc, Brier=brier, CalibSlope=calib_slope, RMSE=rmse, Coverage=coverage)
}

# main loop: simulate and fit - show structure for one dataset per scenario
set.seed(2025)
n <- 1000
nsim <- 500
scenarios <- c("S1_SAS","S2_SkewT","S3_Copula","S4_Misclass","S5_Hetero")
results <- list()

for (sc in scenarios) {
  message("Running scenario ", sc)
  res_sc <- vector("list", nsim)
  for (rep in 1:nsim) {
    # generate covariates
    X <- matrix(runif(n*4, -1, 1), ncol=4)
    x1 <- X[,1]; x2 <- X[,2]; x3 <- X[,3]; x4 <- X[,4]
    mu_i <- 0.5*x1 - 0.3*x2
    logsig2_i <- 0.2*x3
    alpha_i <- 0.4*x4
    kappa_i <- rep(0.5, n)
    # simulate Z according to scenario
    if (sc == "S1_SAS") {
      Z <- simulate_Z_sas(mu_i, logsig2_i, alpha_i, kappa_i, n)
    } else if (sc == "S2_SkewT") {
      Z <- simulate_Z_skewt(mu_i, logsig2_i, alpha_i, kappa_i, n, df=5)
    } else if (sc == "S3_Copula") {
      Z <- simulate_Z_gaussian_re(mu_i, logsig2_i, alpha_i, kappa_i, n)
    } else if (sc == "S4_Misclass") {
      Z <- simulate_Z_sas(mu_i, logsig2_i, alpha_i, kappa_i, n)
    } else if (sc == "S5_Hetero") {
      Z <- simulate_Z_sas(mu_i, logsig2_i, alpha_i, kappa_i, n)
    }
    # measurement errors
    eps_c <- rnorm(n, 0, sqrt(0.25))
    yc <- 1 + 0.8 * Z + eps_c
    eps_b <- rnorm(n)
    D <- as.integer(Z + eps_b > 0)  # latent diagnosis
    if (sc == "S4_Misclass") {
      # flip 10% random observations (simple misclassification)
      flip <- sample(n, size = round(0.1*n))
      yb <- D
      yb[flip] <- 1 - yb[flip]
    } else {
      yb <- D
    }
    # Fit PMM via Stan
    # Create design matrices
    X_mu <- cbind(1, X) # intercept + covariates; adjust p
    X_sigma <- cbind(1, X)
    X_alpha <- cbind(1, X)
    X_kappa <- cbind(1, X)
    stan_data <- list(N=n, yb=yb, yc=yc,
                      p_mu=ncol(X_mu), p_sigma=ncol(X_sigma), p_alpha=ncol(X_alpha), p_kappa=ncol(X_kappa),
                      X_mu=X_mu, X_sigma=X_sigma, X_alpha=X_alpha, X_kappa=X_kappa,
                      sens_prior_a=18, sens_prior_b=3,
                      spec_prior_a=20, spec_prior_b=2)
    fit <- stan(file=stan_file, data=stan_data, iter=2000, warmup=1000, chains=4, control=list(adapt_delta=0.9, max_treedepth=12))
    samp <- rstan::extract(fit, permuted=TRUE)
    # posterior mean predictions:
    # compute posterior predictive p(D=1) using posterior draws: here approximate using Phi(Z_draw)
    # We need to compute Z draws for posterior samples: reconstruct from posterior betas and U draws if accessible.
    # For simplicity, take posterior mean of beta and reconstruct mean Z.
    beta_mu_hat <- apply(samp$beta_mu, 2, mean)
    beta_sigma_hat <- apply(samp$beta_sigma, 2, mean)
    beta_alpha_hat <- apply(samp$beta_alpha, 2, mean)
    beta_kappa_hat <- apply(samp$beta_kappa, 2, mean)
    U_hat <- apply(samp$U, 2, mean)
    mu_hat <- X_mu %*% beta_mu_hat
    logsig2_hat <- X_sigma %*% beta_sigma_hat
    alpha_hat <- X_alpha %*% beta_alpha_hat
    kappa_hat <- X_kappa %*% beta_kappa_hat
    sigma_hat <- sqrt(exp(logsig2_hat))
    Z_hat <- mu_hat + sigma_hat * sinh( (asinh(U_hat) + alpha_hat) / kappa_hat )
    p_bin_hat <- pnorm(Z_hat)
    # predicted binary probability p_hat
    # take posterior mean of pi_sens and pi_spec
    pi_sens_hat <- mean(samp$pi_sens)
    pi_spec_hat <- mean(samp$pi_spec)
    p_yb_hat <- p_bin_hat * pi_sens_hat + (1 - p_bin_hat) * (1 - pi_spec_hat)
    # continuous prediction:
    mu2_hat <- mean(samp$mu2)
    lambda_hat <- mean(samp$lambda)
    yc_hat <- mu2_hat + lambda_hat * Z_hat
    # posterior predictive samples for yc: use posterior param draws for subset
    nsamp_pp <- 200
    idx <- sample.int(length(samp$mu2), nsamp_pp)
    yc_samples <- matrix(NA, n, nsamp_pp)
    for (k in seq_len(nsamp_pp)) {
      mu2_k <- samp$mu2[idx[k]]
      lambda_k <- samp$lambda[idx[k]]
      # reconstruct Z draw via U draw and betas
      beta_mu_k <- samp$beta_mu[idx[k],]
      beta_sigma_k <- samp$beta_sigma[idx[k],]
      beta_alpha_k <- samp$beta_alpha[idx[k],]
      beta_kappa_k <- samp$beta_kappa[idx[k],]
      U_k <- samp$U[idx[k],]
      mu_k <- X_mu %*% beta_mu_k
      logsig2_k <- X_sigma %*% beta_sigma_k
      alpha_k <- X_alpha %*% beta_alpha_k
      kappa_k <- X_kappa %*% beta_kappa_k
      sigma_k <- sqrt(exp(logsig2_k))
      Z_k <- mu_k + sigma_k * sinh( (asinh(U_k) + alpha_k) / kappa_k )
      yc_samples[,k] <- rnorm(n, mu2_k + lambda_k * Z_k, sd=mean(samp$sigma_eps))
    }
    # compute metrics
    mets <- compute_metrics(yb, p_yb_hat, yc, yc_hat, yc_samples)
    res_sc[[rep]] <- list(mets=mets, fit=fit)
    if (rep %% 20 == 0) message("Rep ", rep, " done")
  } # end rep
  results[[sc]] <- res_sc
  # save intermediate
  saveRDS(results, file=paste0("sim_results_", sc, ".rds"))
} # end scenario loop





