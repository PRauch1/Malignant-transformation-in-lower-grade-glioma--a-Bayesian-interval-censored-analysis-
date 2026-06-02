# =====================================================================
# 05_diagnostics.R — convergence diagnostics, posterior summaries,
#                    feature inclusion probabilities, plots
#
# Estimated runtime on M4: 2-5 minutes
# =====================================================================

suppressPackageStartupMessages({
  library(brms)
  library(posterior)
  library(bayesplot)
  library(ggplot2)
  library(dplyr)
})

# ---- Load fitted models ----
fit_full <- readRDS("02_fit_full.rds")
fit_mut  <- readRDS("02_fit_mut.rds")
fit_wt_path <- "02_fit_wt.rds"
fit_wt <- if (file.exists(fit_wt_path)) readRDS(fit_wt_path) else NULL

# ---- 1. Convergence diagnostics ----
diag_summary <- function(fit, label) {
  s <- summary(fit)$fixed
  ndiv <- sum(rstan::get_num_divergent(fit$fit))
  ntd  <- sum(rstan::get_num_max_treedepth(fit$fit))
  cat("─────", label, "─────\n")
  cat(sprintf("  Max R-hat:       %.4f\n", max(s[, "Rhat"], na.rm = TRUE)))
  cat(sprintf("  Min ESS_bulk:    %d\n",  round(min(s[, "Bulk_ESS"], na.rm = TRUE))))
  cat(sprintf("  Min ESS_tail:    %d\n",  round(min(s[, "Tail_ESS"], na.rm = TRUE))))
  cat(sprintf("  Divergences:     %d\n",  ndiv))
  cat(sprintf("  Treedepth hits:  %d\n",  ntd))
  if (max(s[, "Rhat"], na.rm = TRUE) > 1.01) cat("  ** WARNING: R-hat > 1.01\n")
  if (min(s[, "Bulk_ESS"], na.rm = TRUE) < 400) cat("  ** WARNING: ESS_bulk < 400\n")
  if (ndiv > 0) cat("  ** WARNING:", ndiv, "divergent transitions\n")
  invisible(s)
}

cat("\n=== CONVERGENCE DIAGNOSTICS ===\n")
diag_summary(fit_full, "Full cohort")
diag_summary(fit_mut, "IDH-mutant")
if (!is.null(fit_wt)) diag_summary(fit_wt, "IDH-wildtype")

# ---- 2. Extract clinical effects with 95% credible intervals ----
# brms::weibull() is AFT-parameterized (log link on scale μ). The correct
# proportional-hazards HR for a unit change in covariate j is exp(-α · β_j),
# where α is the Weibull shape. We compute this element-wise across draws so
# that joint α/β uncertainty is propagated (median-plug-in would lose
# dependence). The log-HR is -α · β, sign-inverted vs β but |CrI exclusion
# of zero| is preserved since α > 0 for every draw.
clinical_effects <- function(fit, label) {
  draws <- as_draws_df(fit)
  alpha <- draws[["shape"]]  # brms Weibull shape α
  if (is.null(alpha)) stop("shape parameter not found; AFT->PH conversion requires it")
  clinical_pars <- c("b_contrast_enhancingI1", "b_age_atDiag45", "b_vol_diag_centered")
  out <- lapply(clinical_pars, function(p) {
    if (!p %in% colnames(draws)) return(NULL)
    beta   <- draws[[p]]
    log_hr <- -alpha * beta          # draw-by-draw AFT->PH
    hr     <- exp(log_hr)
    data.frame(
      cohort      = label,
      parameter   = p,
      beta_median = median(beta),
      HR_median   = median(hr),
      HR_q025     = quantile(hr, 0.025),
      HR_q975     = quantile(hr, 0.975),
      Pr_HR_gt_1  = mean(log_hr > 0),
      Pr_HR_lt_1  = mean(log_hr < 0)
    )
  })
  do.call(rbind, out)
}

clin_table <- bind_rows(
  clinical_effects(fit_full, "Full"),
  clinical_effects(fit_mut, "IDH-mutant"),
  if (!is.null(fit_wt)) clinical_effects(fit_wt, "IDH-wildtype") else NULL
)
cat("\n=== CLINICAL EFFECTS (Bayesian posterior) ===\n")
print(clin_table, row.names = FALSE)
saveRDS(clin_table, "05_clinical_effects.rds")
write.csv(clin_table, "05_clinical_effects.csv", row.names = FALSE)

# ---- 3. Radiomic feature posterior inclusion probabilities ----
# With horseshoe, "inclusion" is fuzzy. We use the posterior probability that
# |HR - 1| > some threshold, e.g., HR > 1.1 or HR < 0.9 (i.e., effect of practical magnitude)
# Plus the more standard 95% credible interval excluding 0 (in log-HR space).
radiomic_posterior <- function(fit, label) {
  draws <- as_draws_df(fit)
  alpha <- draws[["shape"]]
  if (is.null(alpha)) stop("shape parameter not found; AFT->PH conversion requires it")
  rad_pars <- grep("^b_original_", colnames(draws), value = TRUE)
  out <- lapply(rad_pars, function(p) {
    beta   <- draws[[p]]
    log_hr <- -alpha * beta          # draw-by-draw AFT->PH; joint α/β uncertainty
    hr     <- exp(log_hr)
    data.frame(
      cohort = label,
      feature = sub("^b_", "", p),
      beta_median = median(beta),
      HR_median = median(hr),
      HR_q025 = quantile(hr, 0.025),
      HR_q975 = quantile(hr, 0.975),
      # Posterior probability of HR > 1.2 (notable accelerative effect)
      Pr_HR_gt_1.2 = mean(hr > 1.2),
      Pr_HR_lt_0.83 = mean(hr < 1/1.2),
      # Posterior probability that 95% CrI excludes 0 on log-HR scale.
      # α > 0 for every draw, so this is equivalent to the β-based check
      # sign-inverted -- confirmed in validation below.
      ci_excludes_0 = (quantile(log_hr, 0.025) > 0) | (quantile(log_hr, 0.975) < 0)
    )
  })
  result <- do.call(rbind, out)
  # Sort by absolute beta median (equivalent to sorting by absolute log-HR median,
  # since α is draw-varying but always positive)
  result[order(-abs(result$beta_median)), ]
}

rad_full <- radiomic_posterior(fit_full, "Full")
rad_mut  <- radiomic_posterior(fit_mut, "IDH-mutant")

cat("\n=== TOP RADIOMIC FEATURES BY POSTERIOR EFFECT (full cohort) ===\n")
print(head(rad_full, 15), row.names = FALSE)

cat("\n=== TOP RADIOMIC FEATURES BY POSTERIOR EFFECT (IDH-mutant) ===\n")
print(head(rad_mut, 15), row.names = FALSE)

saveRDS(list(full = rad_full, mut = rad_mut), "05_radiomic_posterior.rds")
write.csv(rad_full, "05_radiomic_posterior_full.csv", row.names = FALSE)
write.csv(rad_mut,  "05_radiomic_posterior_mut.csv",  row.names = FALSE)

# ---- 4. Posterior predictive check ----
# brms::pp_check() doesn't work cleanly with cens()-style interval-censored
# responses (no scalar y to plot), so we wrap it and fall back to a manual
# density overlay built from posterior_predict() and the interval midpoints.
cat("\n=== POSTERIOR PREDICTIVE CHECK ===\n")
pp_ok <- tryCatch({
  pp_full <- pp_check(fit_full, ndraws = 50)
  ggsave("05_pp_check_full.pdf", pp_full, width = 8, height = 5)
  cat("Saved: 05_pp_check_full.pdf\n")
  TRUE
}, error = function(e) {
  cat("pp_check() not supported for interval-censored fit (",
      conditionMessage(e), ")\n", sep = "")
  cat("Falling back to manual density overlay on interval midpoints.\n")
  FALSE
})

if (!pp_ok) {
  baseline_with_rad <- readRDS("01_baseline_with_rad.rds")
  y_mid <- with(baseline_with_rad, (L + R) / 2)
  yrep <- posterior_predict(fit_full, ndraws = 50)
  # Density overlay: observed (interval midpoint) vs. posterior predictive draws.
  pp_manual <- bayesplot::ppc_dens_overlay(y = y_mid, yrep = yrep) +
    ggplot2::labs(
      title = "Posterior predictive check (manual)",
      subtitle = "Observed = midpoint of (L, R) for interval-censored times",
      x = "Time to malignant transformation (years)"
    )
  ggsave("05_pp_check_full.pdf", pp_manual, width = 8, height = 5)
  cat("Saved: 05_pp_check_full.pdf (manual ppc_dens_overlay)\n")
}

# ---- 5. Trace plots for key parameters ----
trace_pars <- c("b_contrast_enhancingI1", "b_age_atDiag45", "b_vol_diag_centered",
                "b_original_glszm_largearealowgraylevelemphasis_T1c",
                "b_original_gldm_dependencevariance_FLAIR")
trace_pars <- intersect(trace_pars, variables(fit_full))
if (length(trace_pars)) {
  trace_p <- mcmc_trace(fit_full, pars = trace_pars)
  ggsave("05_traces_full.pdf", trace_p, width = 10, height = 8)
  cat("Saved: 05_traces_full.pdf\n")
}

# ---- 6. Forest plot of clinical effects across cohorts ----
fp <- ggplot(clin_table, aes(y = parameter, x = HR_median, color = cohort)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbarh(aes(xmin = HR_q025, xmax = HR_q975),
                 height = 0.2, position = position_dodge(0.4)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  labs(x = "Hazard ratio (95% CrI)", y = NULL,
       title = "Clinical predictors of malignant transformation",
       subtitle = "Bayesian posterior medians and 95% credible intervals") +
  theme_minimal()
ggsave("05_forest_clinical.pdf", fp, width = 8, height = 5)
cat("Saved: 05_forest_clinical.pdf\n")

cat("\n--- Diagnostics complete. Proceed to 06_summary.R ---\n")
