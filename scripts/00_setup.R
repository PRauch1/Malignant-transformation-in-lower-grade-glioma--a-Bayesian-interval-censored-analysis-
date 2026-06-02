# =====================================================================
# 00_setup.R — environment setup for Bayesian TTM analysis
# Run this ONCE on your MacBook before any other scripts.
#
# Expected total runtime: 15-40 minutes (most of it is brms/cmdstanr install)
#
# Prerequisites on the Mac:
#   1. R >= 4.3 installed (https://cran.r-project.org)
#   2. Xcode Command Line Tools installed:
#        Open Terminal, run: xcode-select --install
#      If already installed, this command does nothing harmful.
# =====================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- Step 1: install standard CRAN packages ----
needed_cran <- c("dplyr", "readr", "ggplot2", "survival", "posterior",
                 "bayesplot", "tidybayes", "loo")
to_install <- setdiff(needed_cran, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

# ---- Step 2: install cmdstanr (from Stan's R-universe, not CRAN) ----
if (!"cmdstanr" %in% rownames(installed.packages())) {
  install.packages("cmdstanr",
    repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
}
library(cmdstanr)

# ---- Step 3: install CmdStan (the Stan C++ backend) ----
# This downloads ~150MB and compiles. On M4 it takes ~3-6 min.
# If install_cmdstan() complains about make / C++ toolchain, run in Terminal:
#    xcode-select --install
# and retry.
if (is.null(cmdstan_version(error_on_NA = FALSE))) {
  cat("Installing CmdStan...\n")
  install_cmdstan(cores = parallel::detectCores() - 1, overwrite = FALSE)
}
cat("cmdstan version:", cmdstan_version(), "\n")

# ---- Step 4: install brms ----
if (!"brms" %in% rownames(installed.packages())) {
  install.packages("brms")
}
library(brms)

# Direct brms to use cmdstanr backend (faster, no rstan headache)
options(brms.backend = "cmdstanr")

# ---- Step 5: smoke test ----
cat("\n=== Smoke test ===\n")
test_data <- data.frame(y = rnorm(50), x = rnorm(50))
test_fit <- brm(y ~ x, data = test_data, chains = 2, iter = 500,
                refresh = 0, silent = 2, backend = "cmdstanr")
cat("Smoke test fit succeeded.\n")
cat("R-hat for x coefficient:", round(rhat(test_fit)["b_x"], 3), "\n")
cat("\n--- Setup complete. Proceed to 01_load_data.R ---\n")
