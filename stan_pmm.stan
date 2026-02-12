// stan_pmm.stan
data {
  int<lower=1> N;
  int<lower=0,upper=1> yb[N];
  vector[N] yc;
  int<lower=1> p_mu;
  int<lower=1> p_sigma;
  int<lower=1> p_alpha;
  int<lower=1> p_kappa;
  matrix[N,p_mu] X_mu;
  matrix[N,p_sigma] X_sigma;
  matrix[N,p_alpha] X_alpha;
  matrix[N,p_kappa] X_kappa;
  // if using misclassification (optional)
  real<lower=0> sens_prior_a;
  real<lower=0> sens_prior_b;
  real<lower=0> spec_prior_a;
  real<lower=0> spec_prior_b;
}
parameters {
  vector[p_mu] beta_mu;
  vector[p_sigma] beta_sigma;
  vector[p_alpha] beta_alpha;
  vector[p_kappa] beta_kappa;
  real mu2;
  real lambda;
  real<lower=0> sigma_eps;
  // optional misclassification parameters (on probability scale)
  real<lower=0,upper=1> pi_sens;
  real<lower=0,upper=1> pi_spec;
  // latent U (non-centered approach)
  vector[N] U; // standard normal latent draws
}
transformed parameters {
  vector[N] mu;
  vector[N] logsig2;
  vector[N] alpha;
  vector[N] kappa;
  vector[N] sigma;
  vector[N] Z;
  vector[N] p_bin; // probability of binary=1 under SAS latent

  // vectorized linear predictors
  mu      = X_mu * beta_mu;
  logsig2 = X_sigma * beta_sigma;
  alpha   = X_alpha * beta_alpha;

  // compute (vector) log-kappa and force kappa > 0 by exponentiating;
  // assign to the pre-declared kappa (do not redeclare)
  vector[N] log_kappa = X_kappa * beta_kappa;
  kappa = exp(log_kappa);

  for (i in 1:N) {
    sigma[i] = sqrt(exp(logsig2[i]));
    // SAS transform: Z = mu + sigma * sinh( (asinh(U) + alpha)/kappa )
    real asu = asinh(U[i]);
    Z[i] = mu[i] + sigma[i] * sinh( (asu + alpha[i]) / kappa[i] );
    // p_bin: prob(D_i = 1 conditional on Z_i) under standard-normal error
    p_bin[i] = Phi( Z[i] );
  }
}
model {
  // priors (weakly informative)
  beta_mu ~ normal(0,2.5);
  beta_sigma ~ normal(0,2.5);
  beta_alpha ~ normal(0,2.5);
  beta_kappa ~ normal(0,2.5);
  mu2 ~ normal(0,5);
  lambda ~ normal(0,1);
  sigma_eps ~ cauchy(0,2.5);
  pi_sens ~ beta(sens_prior_a, sens_prior_b);
  pi_spec ~ beta(spec_prior_a, spec_prior_b);
  U ~ normal(0,1);
  // likelihood:
  for (i in 1:N) {
    // continuous measurement: Yc ~ N(mu2 + lambda * Z, sigma_eps^2)
    yc[i] ~ normal( mu2 + lambda * Z[i], sigma_eps );
    // binary measurement with possible misclassification:
    real p_yb = p_bin[i] * pi_sens + (1 - p_bin[i]) * (1 - pi_spec);
    yb[i] ~ bernoulli( p_yb );
  }
}
generated quantities{
  // posterior predictive samples or summary can be added here if desired
}
