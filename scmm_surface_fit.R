---
title: "scmm_surface_fit.R"
output: html_document
---



# scmm_surface_fit.R
library(tidyverse)
library(mgcv)   # for gam with tensor product smooths
library(fields) # optional alternative smoother

# 1. build grid for moments
mu_grid <- seq(-1.5,1.5,length.out=15)
logsig_grid <- seq(log(0.5), log(2), length.out=12)
alpha_grid <- seq(-1,1,length.out=11)
kappa_grid <- seq(0.6,1.6,length.out=9)

grid <- expand.grid(mu=mu_grid, logsig2=logsig_grid, alpha=alpha_grid, kappa=kappa_grid)
nrep <- 200  # number of replicates per grid point

compute_p <- function(mu, logsig2, alpha, kappa, nrep=200) {
  sigma <- sqrt(exp(logsig2))
  # simulate many U's and compute Z = mu + sigma * sinh((asinh(U) + alpha)/kappa)
  U <- rnorm(nrep*2000) # big draw to reduce Monte Carlo noise; adjust as needed
  # vectorized approach:
  Um <- matrix(U, nrow=nrep, ncol=2000, byrow=TRUE)[,1:nrep] # simpler: do smaller batches
  # simpler implementation: simulate nrep draws for each grid point
  U2 <- rnorm(nrep)
  asu <- asinh(U2)
  Z <- mu + sigma * sinh( (asu + alpha)/kappa )
  # latent binary D ~ 1(Z + eps > 0), eps ~ N(0,1)
  eps1 <- rnorm(nrep)
  D <- as.integer(Z + eps1 > 0)
  mean(D)
}

# for speed, loop & store
set.seed(2025)
grid$pi_hat <- NA_real_
for (i in seq_len(nrow(grid))) {
  gr <- grid[i,]
  # simulate nrep draws for this grid point
  U <- rnorm(nrep)
  asu <- asinh(U)
  sigma <- sqrt(exp(gr$logsig2))
  Z <- gr$mu + sigma * sinh( (asu + gr$alpha)/gr$kappa )
  eps1 <- rnorm(nrep)
  D <- as.integer(Z + eps1 > 0)
  grid$pi_hat[i] <- mean(D)
  if (i %% 200 == 0) message("Computed ", i, " / ", nrow(grid))
}

# fit a flexible surface: gam with tensor product smooth
# Use family=gaussian since pi_hat is continuous between 0 and 1
library(mgcv)
fit_g <- gam(pi_hat ~ te(mu, logsig2, alpha, kappa, k=rep(6,4)), data=grid, method="REML")

# save
saveRDS(fit_g, file="scmm_surface_gam.rds")
write.csv(grid, file="scmm_grid_pi.csv", row.names=FALSE)
message("SCMM surface fit saved.")



