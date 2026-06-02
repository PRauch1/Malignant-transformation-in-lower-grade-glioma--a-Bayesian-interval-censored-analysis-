# =====================================================================
# 02_bayes_baseline.R — primary Bayesian analysis
#
# Model: Bayesian interval-censored Weibull PH with regularized horseshoe
#        prior on radiomic coefficients, weakly informative priors on
#        clinical coefficients.
#
# Endpoint: time to malignant transformation (TTM), interval-censored
#           with both pre-surgical (Case 2) and post-surgical (Case 3)
#           events properly bounded.
#
# Cohorts:
#   (1) Full cohort (n=155, 77 events)
#   (2) IDH-mutant subgroup (n=109, 48 events)
#   (3) IDH-wildtype subgroup (n=42, 26 events)
#
# Estimated runtime on M4: 15-25 min total for all three fits
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(brms)
  library(posterior)
  library(bayesplot)
  library(loo)
})
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores() - 1)
set.seed(2026)

baseline <- readRDS("01_baseline_with_rad.rds")

# ---- brms interval-censored survival setup ----
# brms expects: time1 (lower), time2 (upper for cens="interval"), or single time + cens
# brms cens codes: "none" (event), "right", "left", "interval"
# We code right-censored when R == Inf, interval-censored otherwise.
prep_brms <- function(d) {
  d %>%
    mutate(
      time1 = ifelse(is.infinite(R), L, L),
      time2 = ifelse(is.infinite(R), L, R),
      cens = ifelse(is.infinite(R), "right", "interval"),
      time1 = ifelse(time1 == 0, 0.001, time1),  # brms dislikes time=0
      time2 = ifelse(time2 == time1 & cens == "interval", time1 + 0.001, time2)
    )
}

baseline_brms <- prep_brms(baseline)

# ---- Identify radiomic features for the horseshoe-penalized block ----
rad_cols <- grep("^original_", colnames(baseline_brms), value = TRUE)
cat("Radiomic features in horseshoe block:", length(rad_cols), "\n")

# ---- Build formula ----
# Clinical effects: weakly informative Student-t(3, 0, 2.5) priors via default in brms
# Radiomic effects: regularized horseshoe via set_prior(horseshoe(...))
clin_terms <- c("contrast_enhancingI", "age_atDiag45", "vol_diag_centered")
all_terms <- c(clin_terms, rad_cols)

f_full <- bf(
  paste0("time1 | cens(cens, time2) ~ ", paste(all_terms, collapse = " + "))
)

# ---- Priors ----
# Regularized horseshoe applied class-wide to all population-level coefficients
# (clinical + radiomic). brms 2.23 does not permit per-coefficient priors alongside
# a class-wide horseshoe; see RESULTS_REPORT.md §H for the methods-reporting
# language covering this deviation from the original plan.
# par_ratio = (expected_nonzero / total) ~ 10/102 = 0.1; df_global=1, df_slab=4,
# scale_slab=2 keeps strong clinical effects (CE) essentially unshrunk.
priors_full <- c(
  prior(horseshoe(par_ratio = 0.1, df = 1, df_global = 1,
                  scale_global = 2, scale_slab = 2, df_slab = 4),
        class = "b")
)

# ---- Fit on full cohort ----
cat("\n=== Fitting full-cohort model (Weibull PH, horseshoe on radiomics) ===\n")
cat("This will take 5-10 minutes on M4...\n")
fit_full <- brm(
  formula = f_full,
  data = baseline_brms,
  family = weibull(),
  prior = priors_full,
  chains = 4, iter = 8000, warmup = 4000,
  control = list(adapt_delta = 0.999, max_treedepth = 14),
  seed = 2026, refresh = 200, silent = 1
)
saveRDS(fit_full, "02_fit_full.rds")
cat("Saved: 02_fit_full.rds\n")

# ---- Diagnostics ----
cat("\n--- Convergence diagnostics (full cohort) ---\n")
sumr <- summary(fit_full)
cat("  Max R-hat:", round(max(sumr$fixed[, "Rhat"], na.rm = TRUE), 4), "\n")
cat("  Min ESS_bulk:", round(min(sumr$fixed[, "Bulk_ESS"], na.rm = TRUE)), "\n")
cat("  Divergent transitions:", sum(rstan::get_num_divergent(fit_full$fit)), "\n")
cat("  Max treedepth hits:", sum(rstan::get_num_max_treedepth(fit_full$fit)), "\n")

# ---- Fit on IDH-mutant subgroup ----
cat("\n=== Fitting IDH-mutant subgroup model ===\n")
baseline_mut <- baseline_brms %>% filter(idh == 1, !is.na(idh))
cat("N:", nrow(baseline_mut), "events:", sum(baseline_mut$event == 1, na.rm = TRUE), "\n")

fit_mut <- brm(
  formula = f_full,
  data = baseline_mut,
  family = weibull(),
  prior = priors_full,
  chains = 4, iter = 8000, warmup = 4000,
  control = list(adapt_delta = 0.999, max_treedepth = 14),
  seed = 2026, refresh = 200, silent = 1
)
saveRDS(fit_mut, "02_fit_mut.rds")

# ---- Fit on IDH-wildtype subgroup (smaller, expect convergence challenges) ----
cat("\n=== Fitting IDH-wildtype subgroup model ===\n")
baseline_wt <- baseline_brms %>% filter(idh == 0, !is.na(idh))
cat("N:", nrow(baseline_wt), "events:", sum(baseline_wt$event == 1, na.rm = TRUE), "\n")

fit_wt <- try(brm(
  formula = f_full,
  data = baseline_wt,
  family = weibull(),
  prior = priors_full,
  chains = 4, iter = 8000, warmup = 4000,
  control = list(adapt_delta = 0.999, max_treedepth = 14),
  seed = 2026, refresh = 200, silent = 1
), silent = TRUE)
if (!inherits(fit_wt, "try-error")) {
  saveRDS(fit_wt, "02_fit_wt.rds")
} else {
  cat("WARNING: IDH-wildtype model failed (likely n too small).\n")
  cat("  Error:", attr(fit_wt, "condition")$message, "\n")
}

cat("\n--- Primary Bayesian models complete. Proceed to 03_bayes_tvc.R ---\n")
