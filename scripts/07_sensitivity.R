# =====================================================================
# 07_sensitivity.R — sensitivity analyses
#
# 1. Horseshoe prior sensitivity: refit primary model with par_ratio
#    in {0.05, 0.1, 0.2}. Confirm conclusions are robust.
# 2. Adult-only cohort: exclude n=6 pediatric patients (age < 18).
# 3. Era adjustment: include diagnosis era as covariate to address
#    scanner-generation confounding.
#
# Run AFTER 02_bayes_baseline.R has produced 02_fit_full.rds.
# Estimated runtime on M4: 30-60 minutes (3 sensitivity refits)
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(posterior)
})
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores() - 1)
set.seed(2026)

baseline <- readRDS("01_baseline_with_rad.rds")

# Same prep as in 02_bayes_baseline.R
prep_brms <- function(d) {
  d %>%
    mutate(
      time1 = ifelse(is.infinite(R), L, L),
      time2 = ifelse(is.infinite(R), L, R),
      cens  = ifelse(is.infinite(R), "right", "interval"),
      time1 = ifelse(time1 == 0, 0.001, time1),
      time2 = ifelse(time2 == time1 & cens == "interval", time1 + 0.001, time2)
    )
}

baseline_brms <- prep_brms(baseline)
rad_cols <- grep("^original_", colnames(baseline_brms), value = TRUE)
clin_terms <- c("contrast_enhancingI", "age_atDiag45", "vol_diag_centered")
all_terms  <- c(clin_terms, rad_cols)
f_full <- bf(paste0("time1 | cens(cens, time2) ~ ", paste(all_terms, collapse = " + ")))

# Per-coefficient Student-t priors removed: brms 2.23 does not permit per-coef
# priors alongside a class-wide horseshoe. All population-level coefficients
# (clinical + radiomic + any era covariate) receive the class-wide horseshoe
# specified per loop iteration. See RESULTS_REPORT.md §H.
# (No leading `base_priors` object is used — c(NULL, prior(...)) drops the
# brmsprior class, so priors are assigned directly below.)

# ---- Sensitivity 1: horseshoe par_ratio ----
cat("\n=== SENSITIVITY 1: horseshoe par_ratio ===\n")
sens_results <- list()
for (pr in c(0.05, 0.1, 0.2)) {
  cat("\nFitting with par_ratio =", pr, "...\n")
  prior_pr <- prior_string(
    paste0("horseshoe(par_ratio = ", pr,
           ", df = 1, df_global = 1, scale_global = 1, scale_slab = 2, df_slab = 4)"),
    class = "b")
  fit_pr <- brm(f_full, data = baseline_brms, family = weibull(),
                prior = prior_pr,
                chains = 4, iter = 4000, warmup = 2000,
                control = list(adapt_delta = 0.995, max_treedepth = 14),
                seed = 2026, refresh = 0, silent = 2)
  draws <- as_draws_df(fit_pr)
  alpha <- draws[["shape"]]   # AFT->PH conversion requires joint shape α draws
  # Extract clinical posteriors and count of "active" radiomic features
  ce_beta <- draws[["b_contrast_enhancingI1"]]
  ce_hr   <- exp(-alpha * ce_beta)  # draw-by-draw PH hazard ratio
  n_active <- 0
  for (rc in rad_cols) {
    coef_name <- paste0("b_", rc)
    if (coef_name %in% colnames(draws)) {
      x <- draws[[coef_name]]
      # ci_excludes_0 is invariant under multiplication by α > 0, so the
      # β-based check gives the same count as a log-HR-based check.
      if (quantile(x, 0.025) > 0 || quantile(x, 0.975) < 0) n_active <- n_active + 1
    }
  }
  sens_results[[as.character(pr)]] <- list(
    par_ratio = pr,
    HR_CE_median = median(ce_hr),
    HR_CE_q025   = quantile(ce_hr, 0.025),
    HR_CE_q975   = quantile(ce_hr, 0.975),
    n_active_radiomic = n_active
  )
  saveRDS(fit_pr, paste0("07_sens_horseshoe_pr", pr, ".rds"))
}
sens_df <- do.call(rbind, lapply(sens_results, as.data.frame))
print(sens_df, row.names = FALSE)
saveRDS(sens_df, "07_sens_horseshoe_summary.rds")
write.csv(sens_df, "07_sens_horseshoe_summary.csv", row.names = FALSE)

# ---- Sensitivity 2: adult-only cohort ----
cat("\n=== SENSITIVITY 2: adult-only cohort (exclude age < 18) ===\n")
# age_atDiag45 is age - 45, so age < 18 corresponds to age_atDiag45 < -27
adult <- baseline_brms %>% filter(age_atDiag45 >= -27)
cat("Adult-only N:", nrow(adult), "events:", sum(adult$event == 1, na.rm = TRUE), "\n")

prior_default <- prior(horseshoe(par_ratio = 0.1, df = 1, df_global = 1,
                                 scale_global = 1, scale_slab = 2, df_slab = 4),
                       class = "b")

fit_adult <- brm(f_full, data = adult, family = weibull(),
                 prior = prior_default,
                 chains = 4, iter = 4000, warmup = 2000,
                 control = list(adapt_delta = 0.995, max_treedepth = 14),
                 seed = 2026, refresh = 0, silent = 2)
saveRDS(fit_adult, "07_sens_adult.rds")

# Compare adult-only clinical effects to full-cohort
draws_adult <- as_draws_df(fit_adult)
fit_full <- readRDS("02_fit_full.rds")
draws_full <- as_draws_df(fit_full)

alpha_full  <- draws_full[["shape"]]
alpha_adult <- draws_adult[["shape"]]
cat("\nClinical effect comparison (full vs adult-only), HR = exp(-α · β) draw-by-draw:\n")
for (par in c("b_contrast_enhancingI1", "b_age_atDiag45", "b_vol_diag_centered")) {
  if (par %in% colnames(draws_full) && par %in% colnames(draws_adult)) {
    hr_full  <- exp(-alpha_full  * draws_full[[par]])
    hr_adult <- exp(-alpha_adult * draws_adult[[par]])
    cat(sprintf("  %s\n", par))
    cat(sprintf("    Full:        HR = %.3f (%.3f, %.3f)\n",
                median(hr_full),
                quantile(hr_full, 0.025),
                quantile(hr_full, 0.975)))
    cat(sprintf("    Adult-only:  HR = %.3f (%.3f, %.3f)\n",
                median(hr_adult),
                quantile(hr_adult, 0.025),
                quantile(hr_adult, 0.975)))
  }
}

# ---- Sensitivity 3: era adjustment ----
# Era info is not currently in baseline_brms; need to add from raw clinical data
cat("\n=== SENSITIVITY 3: era adjustment ===\n")
# Try to find the source xlsx to derive era
xlsx_candidates <- c(
  "../source_data/Data_LGG_Clinical_update.xlsx",
  Sys.getenv("LGG_CLINICAL_XLSX", unset = NA)
)
xlsx_path <- xlsx_candidates[file.exists(xlsx_candidates)][1]
if (!is.na(xlsx_path)) {
  if (!"readxl" %in% rownames(installed.packages())) install.packages("readxl")
  raw <- readxl::read_excel(xlsx_path, sheet = "Datensaetzte")
  raw$date_diagnosis <- as.Date(raw$date_diagnosis)
  raw$era <- ifelse(format(raw$date_diagnosis, "%Y") < "2010", "early", "late")
  baseline_era <- baseline_brms %>%
    left_join(raw %>% select(patient_number, era), by = "patient_number") %>%
    mutate(era = factor(era, levels = c("early", "late")))
  cat("Era distribution:\n"); print(table(baseline_era$era))
  
  f_era <- bf(paste0("time1 | cens(cens, time2) ~ era + ",
                     paste(all_terms, collapse = " + ")))
  # Per-coef prior on `eralate` removed for the same brms-2.23 constraint;
  # era now absorbs the class-wide horseshoe. See RESULTS_REPORT.md §H.
  prior_era <- prior_default
  fit_era <- brm(f_era, data = baseline_era, family = weibull(),
                 prior = prior_era,
                 chains = 4, iter = 4000, warmup = 2000,
                 control = list(adapt_delta = 0.995, max_treedepth = 14),
                 seed = 2026, refresh = 0, silent = 2)
  saveRDS(fit_era, "07_sens_era.rds")
  draws_era <- as_draws_df(fit_era)
  alpha_era <- draws_era[["shape"]]
  era_beta  <- draws_era[["b_eralate"]]
  era_hr    <- exp(-alpha_era * era_beta)  # AFT->PH draw-by-draw
  cat("\nEra effect (late vs early): HR =",
      sprintf("%.3f (%.3f, %.3f)\n",
              median(era_hr),
              quantile(era_hr, 0.025),
              quantile(era_hr, 0.975)))
} else {
  cat("Source xlsx not found; era sensitivity SKIPPED.\n")
  cat("To run: place Data_LGG_Clinical_update.xlsx in ../source_data/ or set\n")
  cat("LGG_CLINICAL_XLSX env variable, then re-run this script.\n")
}

cat("\n--- Sensitivity analyses complete. ---\n")
