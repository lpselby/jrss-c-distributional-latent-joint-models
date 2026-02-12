---
  title: "generate_synthetic_nhanes"
output: html_document
---
  
  install.packages(c("dplyr","tibble","readr","rstan","cmdstanr","survey"))

library(dplyr)
library(tibble)



set.seed(2025)
n <- 2943  # match analytic sample size in paper

# 1. Demographics
age <- round(rnorm(n, mean = 50, sd = 16))           # approximate adult ages
age[age < 20] <- sample(20:25, sum(age < 20), replace=TRUE)
sex <- rbinom(n, 1, 0.48)                            # 1 = male, 0 = female
race_cat <- sample(c("Non-Hispanic White","Non-Hispanic Black","Hispanic","Other"),
                   size=n, replace=TRUE, prob=c(0.62,0.12,0.17,0.09))

# 2. Socioeconomic
education_cat <- sample(c("Less than HS","HS/GED","Some college","College+"),
                        size=n, replace=TRUE, prob=c(0.12,0.28,0.30,0.30))
income_pov <- pmax(0.1, rlnorm(n, mean=0.2, sd=0.6)) # income-to-poverty ratio (>=0.1)

# 3. Clinical and biomarkers
bmi <- rnorm(n, mean=28.5, sd=6.0)                  # body-mass index
sbp <- rnorm(n, mean=125, sd=15)                    # systolic bp
dbp <- rnorm(n, mean=78, sd=10)                     # diastolic bp
wbc <- rnorm(n, mean=6.6, sd=1.6)                   # white blood cell count (10^3/µL)
smoker <- sample(c("Never","Former","Current"), size=n, replace=TRUE, prob=c(0.55,0.28,0.17))
phys_act <- sample(c("Low","Moderate","High"), size=n, replace=TRUE, prob=c(0.35,0.45,0.20))

# 4. Design variables (weights, strata, PSU)
# Simulate survey weights with reasonable variability and sum near sample size * design factor
base_weight <- runif(n, 0.5, 3.0)
# create strata and psu as integers
strata <- sample(1:40, n, replace=TRUE)
psu <- sample(1:200, n, replace=TRUE)
# create MEC weight-like variable and normalized weight for Stan use
wt_mec <- base_weight * (1 + 0.2 * (age - mean(age))/sd(age))
wt_norm <- wt_mec / mean(wt_mec)   # normalized weights (mean 1)

# 5. Create latent variable Z and outcomes
# Build covariate matrix used in manuscript: x1..x4 from simulated covariates
# Use standardized continuous covariates: age, bmi, sbp, wbc
age_std <- as.numeric(scale(age))
bmi_std <- as.numeric(scale(bmi))
sbp_std <- as.numeric(scale(sbp))
wbc_std <- as.numeric(scale(wbc))

# latent moment functions (as in manuscript)
mu_z <- 0.5 * age_std - 0.3 * bmi_std
logsig2_z <- 0.2 * sbp_std
alpha_z <- 0.4 * wbc_std
kappa_z <- rep(0.5, n)   # fixed for simplicity

# simulate latent Z under sinh-arcsinh-like transform approximation
# draw U ~ N(0,1), then SAS transform: Z = mu + sigma * sinh((asinh(U) + alpha)/kappa)
U <- rnorm(n)
sigma_z <- sqrt(exp(logsig2_z))
Z <- mu_z + sigma_z * sinh( (asinh(U) + alpha_z) / kappa_z )

# continuous outcome: fasting glucose (mg/dL) scaled to realistic values
# choose mu2=90 baseline and lambda=10 to give glucose in plausible range after transform
mu2 <- 90
lambda <- 10
eps_c <- rnorm(n, sd = sqrt(0.25)) # small measurement noise
glucose <- mu2 + lambda * Z + eps_c
# clamp glucose to plausible lab values
glucose <- pmax(40, pmin(400, glucose))

# binary outcome: self-reported diabetes (with some misclassification)
# true latent diagnosis D = 1(Z + eps_b > 0)
eps_b <- rnorm(n, sd = 1)
D_true <- as.integer(Z + eps_b > 0)
# add misclassification: sensitivity 0.86, specificity 0.91 (approx)
sens <- 0.86; spec <- 0.91
yb <- sapply(seq_len(n), function(i) {
  if (D_true[i]==1) {
    rbinom(1,1,sens)
  } else {
    rbinom(1,1,1-spec)
  }
})

# 6. Add some realistic covariate noise / missingness to test pipeline
# Introduce 2% missingness at random in WBC and 1% in glucose (simulate typical missing)
wbc[sample(n, size = round(0.02*n))] <- NA
glucose[sample(n, size = round(0.01*n))] <- NA

# 7. Assemble data.frame
df <- tibble(
  id = sprintf("SYN%04d", 1:n),
  age = age,
  sex = sex,
  race = race_cat,
  education = education_cat,
  income_pov = income_pov,
  bmi = bmi,
  sbp = sbp,
  dbp = dbp,
  wbc = wbc,
  smoker = smoker,
  phys_act = phys_act,
  wt_mec = wt_mec,
  wt_norm = wt_norm,
  strata = strata,
  psu = psu,
  # latent & outcomes
  Z = Z,
  fasting_glucose = glucose,
  diabetes_self = yb,
  diabetes_true = D_true
)

# 8. Standardize continuous covariates used in models (as in manuscript)
df <- df %>%
  mutate(across(c(age, bmi, sbp, dbp, wbc, income_pov),
                ~ as.numeric(scale(.)), .names = "std_{col}"))

# 9. Save CSV and small README
write.csv(df, "nhanes_synthetic.csv", row.names = FALSE)

readme_text <- c(
  "Synthetic NHANES-style dataset for reproducible code testing",
  "",
  "File: nhanes_synthetic.csv",
  "Rows: 2943",
  "Columns include:",
  "  id: synthetic ID",
  "  age, sex (0=female,1=male), race, education, income_pov",
  "  bmi, sbp, dbp, wbc, smoker, phys_act",
  "  wt_mec: simulated survey weight",
  "  wt_norm: normalized weight (mean 1)",
  "  strata, psu: simulated survey design variables",
  "  Z: simulated latent variable (for demonstration only)",
  "  fasting_glucose: continuous outcome (mg/dL), contains ~1% missing",
  "  diabetes_self: synthetic self-reported diabetes (binary), contains misclassification",
  "  diabetes_true: underlying synthetic diagnosis (latent)",
  "  std_*: standardized versions of continuous covariates used in modelling",
  "",
  "NOTE: This dataset is synthetic and intended for method testing and review only."
)

writeLines(readme_text, con = "README_synthetic.txt")

cat("Synthetic NHANES-style dataset 'nhanes_synthetic.csv' and README_synthetic.txt created.\n")