# =====================================================================
# 11_intccr_finegray_ic.R — Fine-Gray subdistribution hazards
#                           natively for interval-censored
#                           competing risks
#
# Addresses statistician review (Baran, May 2026): Fine-Gray
# in cmprsk::crr has no native interval-censoring routine; 04 used
# midpoint plug-in; 10 used multiple imputation. This script uses
# the intccr package (Park, Bakoyannis & Yiannoutsos 2021, JSS 100(8))
# which fits the Fine-Gray subdistribution hazards model directly to
# interval-censored competing-risks data via sieve maximum likelihood
# with B-spline approximation of the cumulative incidence functions.
#
# Methodology note: Baran's email suggested SmoothHazard, which fits
# an illness-death model with cause-specific hazards. We tried that
# first; idm() did not converge on this dataset under either Weibull
# or spline baselines (the optimizer stalls at the zero-coefficient
# starting point despite CE having a strong effect). intccr is the
# more direct tool for the same methodological gap: it estimates the
# subdistribution-hazards version of Fine-Gray with proper interval-
# censoring of the primary event, matching what 04/10 estimate but
# without midpoint plug-in or multiple imputation.
#
# Estimated runtime on M4: 1-3 minutes per model (analytic LS variance,
# nboot = 0). Bootstrap (nboot = 200) is supported by the package but
# took > 1 h on this dataset and was abandoned in favour of the LS
# variance estimator that Park & Bakoyannis (2021, §3.4) derived as a
# closed-form alternative.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

if (!"intccr" %in% rownames(installed.packages())) {
  install.packages("intccr", repos = "https://cloud.r-project.org")
}
library(intccr)

set.seed(2026)

baseline <- readRDS("01_baseline_with_rad.rds")

# ---- Load death information ----
death_lookup <- read.csv("../data/death_info.csv", stringsAsFactors = FALSE)

# ---- Build intccr Surv2(v, u, event) inputs ----
# Conventions (Park & Bakoyannis 2021):
#   v = last observation time before the failure event (or last fu if censored)
#   u = first observation time after the failure event (or Inf if censored,
#       or u = v for exact event times)
#   event = 0 (censored), 1 (cause 1 = MT), 2 (cause 2 = death without MT)
#
# Our mapping:
#   - MT confirmed (event = 1 in baseline, R finite): v = L, u = R, e = 1
#   - Censored (event = 0 in baseline, R = Inf):       v = last_fu, u = Inf, e = 0
#   - Dead without MT (dead_event = 1, time known):    v = u = time_to_death, e = 2
#   - Dead without MT (time unknown, 6 patients):      treat as censored at last_fu
#     (We do not know when they died; the alternative of interval-censoring the
#      death time from last_fu to Inf is supported by intccr but the death
#      indicator is then a 0/1 flag without time; we conservatively right-censor.)

dat <- baseline %>%
  left_join(death_lookup, by = "patient_number") %>%
  mutate(
    death_known = as.integer(dead_event == 1 & !is.na(time_to_death_yr)),
    last_fu_known = ifelse(is.na(last_followup_yr_from_diag),
                           ifelse(is.finite(R), R, L),
                           last_followup_yr_from_diag),
    # Primary event coding (subdistribution hazards: 1 = MT, 2 = death-without-MT)
    ic_event = case_when(
      event == 1                ~ 1L,       # MT confirmed
      death_known == 1 & event != 1 ~ 2L,   # death without MT (time known)
      TRUE                      ~ 0L        # censored
    ),
    # intccr::Surv2 requires v < u strictly for events. For our exact
    # death-date events we bracket by 1 day (1/365.25 yr) to satisfy this.
    one_day = 1 / 365.25,
    ic_v = case_when(
      ic_event == 1 ~ L,
      ic_event == 2 ~ pmax(0, pmax(L, time_to_death_yr) - one_day),
      TRUE          ~ last_fu_known
    ),
    ic_u = case_when(
      ic_event == 1 ~ ifelse(is.finite(R), R, last_fu_known),
      ic_event == 2 ~ pmax(L, time_to_death_yr) + one_day,
      TRUE          ~ Inf
    ),
    # 01_load_data.R relabels contrast_enhancingI to a factor with
    # labels "0"/"1" (so the brms dummy name matches "contrast_enhancingI1").
    # Cast back to a numeric 0/1 indicator for intccr's design matrix.
    CE_num = as.integer(as.character(contrast_enhancingI) == "1")
  ) %>%
  # Patient 285 special case: death at 0.914 yr precedes MT lower bound 0.931.
  # Reclassify as cause-2 (death) with the death time, not MT.
  mutate(
    ic_event = ifelse(patient_number == 285 & ic_event == 0 & death_known == 1,
                      2L, ic_event),
    ic_v     = ifelse(patient_number == 285,
                      pmax(0, time_to_death_yr - one_day), ic_v),
    ic_u     = ifelse(patient_number == 285,
                      time_to_death_yr + one_day, ic_u)
  ) %>%
  # Final safety: any event row where v == u (e.g. degenerate MT bracket)
  # gets widened by 1 day. Censored rows (u = Inf) are unaffected.
  mutate(
    ic_u = ifelse(ic_event > 0 & is.finite(ic_u) & ic_u <= ic_v,
                  ic_v + one_day, ic_u)
  ) %>%
  filter(!is.na(ic_v), !is.na(ic_u), ic_v >= 0,
         (ic_u > ic_v) | (!is.finite(ic_u)))

cat("intccr cohort:", nrow(dat), "of 155 patients\n")
cat("  Event distribution: ", "\n")
print(table(dat$ic_event, useNA = "ifany",
            dnn = "0=cens, 1=MT, 2=death"))
cat("\n  v summary (yrs):\n");      print(summary(dat$ic_v))
cat("  u summary (yrs, Inf = censored):\n"); print(summary(dat$ic_u[is.finite(dat$ic_u)]))
cat("  # u = Inf (censored): ", sum(!is.finite(dat$ic_u)), "\n\n")

# ---- Fit Fine-Gray (alpha = c(0, 0): both causes proportional sub-dist hazards) ----
# Clinical-only specification
cat("============================================================\n")
cat("=== A. Fine-Gray IC — clinical covariates only           ===\n")
cat("============================================================\n")
t0 <- Sys.time()
fit_clin <- ciregic(
  formula = Surv2(v = ic_v, u = ic_u, event = ic_event) ~
            CE_num + age_atDiag45 + vol_diag_centered,
  data    = dat,
  alpha   = c(0, 0),
  nboot   = 0,            # analytic least-squares SE (fast); 200-boot
                          # alternative documented but takes hours on M4
  do.par  = FALSE
)
cat("Fit time:", format(round(Sys.time() - t0, 1)), "\n\n")
print(summary(fit_clin))

cat("\n============================================================\n")
cat("=== B. Fine-Gray IC — clinical + 2-feature signature     ===\n")
cat("============================================================\n")
t0 <- Sys.time()
fit_full <- ciregic(
  formula = Surv2(v = ic_v, u = ic_u, event = ic_event) ~
            CE_num + age_atDiag45 + vol_diag_centered +
            original_glszm_largearealowgraylevelemphasis_T1c +
            original_gldm_dependencevariance_FLAIR,
  data    = dat,
  alpha   = c(0, 0),
  nboot   = 0,            # analytic least-squares SE (fast); 200-boot
                          # alternative documented but takes hours on M4
  do.par  = FALSE
)
cat("Fit time:", format(round(Sys.time() - t0, 1)), "\n\n")
print(summary(fit_full))

# ---- Extract MT-specific (event 1) sHR table ----
extract_intccr <- function(fit, label) {
  # ciregic structure (verified empirically on the fitted object):
  #   fit$varnames    : character vector of p covariate names
  #   fit$coefficients: flat numeric vector of length 2*p, stacked as
  #                     (cov1_event1, ..., covp_event1, cov1_event2, ..., covp_event2)
  #   fit$vcov        : (2p x 2p) covariance with dimnames "var,event type N"
  co  <- fit$varnames
  est <- fit$coefficients
  vcv <- fit$vcov
  p   <- length(co)
  est1 <- est[1:p];           se1 <- sqrt(diag(vcv)[1:p])
  est2 <- est[(p + 1):(2*p)]; se2 <- sqrt(diag(vcv)[(p + 1):(2*p)])
  data.frame(
    model = label,
    cause = rep(c("1=MT", "2=Death"), each = p),
    term  = rep(co, 2),
    sHR        = exp(c(est1, est2)),
    sHR_lower  = exp(c(est1 - 1.96 * se1, est2 - 1.96 * se2)),
    sHR_upper  = exp(c(est1 + 1.96 * se1, est2 + 1.96 * se2)),
    se_log_sHR = c(se1, se2),
    z = c(est1 / se1, est2 / se2),
    p = 2 * pnorm(-abs(c(est1 / se1, est2 / se2))),
    row.names = NULL
  )
}

tab_clin <- extract_intccr(fit_clin, "clinical_only")
tab_full <- extract_intccr(fit_full, "with_signature")

cat("\n=== Pooled intccr Fine-Gray IC results ===\n")
cat("\n  Clinical only:\n"); print(tab_clin, row.names = FALSE, digits = 4)
cat("\n  With 2-feature signature:\n"); print(tab_full, row.names = FALSE, digits = 4)

# ---- Three-way comparison: 04 midpoint, 10 MI, 11 IC ----
cat("\n============================================================\n")
cat("=== CE-coefficient comparison across CR specifications  ===\n")
cat("============================================================\n")
fg_04 <- tryCatch(readRDS("04_finegray_results.rds"), error = function(e) NULL)
fg_10 <- tryCatch(readRDS("10_finegray_mi_results.rds"), error = function(e) NULL)

get_ce <- function(df) {
  if (is.null(df)) return(NULL)
  hit <- df[grep("^contrast|^CE", df$term), ]
  if (nrow(hit) > 0) hit[1, ] else NULL
}

ce_04 <- get_ce(fg_04$clinical_only)
ce_10 <- get_ce(fg_10$clinical_only)
ce_11 <- tab_clin[tab_clin$cause == "1=MT" & tab_clin$term == "CE_num", ]

cat("CE sub-distribution HR for MT (cause 1), clinical-only spec:\n")
if (!is.null(ce_04))
  cat(sprintf("  04 cmprsk midpoint:     sHR = %6.3f (%6.3f, %6.3f), p = %.4f\n",
              ce_04$sHR, ce_04$sHR_lower, ce_04$sHR_upper, ce_04$p))
if (!is.null(ce_10))
  cat(sprintf("  10 cmprsk + MI:         sHR = %6.3f (%6.3f, %6.3f), p = %.4f\n",
              ce_10$sHR, ce_10$sHR_lower, ce_10$sHR_upper, ce_10$p))
cat(sprintf("  11 intccr (native IC):  sHR = %6.3f (%6.3f, %6.3f), p = %.4f\n",
            ce_11$sHR, ce_11$sHR_lower, ce_11$sHR_upper, ce_11$p))

# ---- Save ----
saveRDS(list(clinical_only = fit_clin, with_signature = fit_full,
             tab_clin = tab_clin, tab_full = tab_full),
        "11_intccr_results.rds")

if (dir.exists("final_summary")) {
  out_tab <- rbind(tab_clin, tab_full)
  write.csv(out_tab,
            "final_summary/table7_intccr_finegray_ic.csv", row.names = FALSE)
  cat("\nSaved: final_summary/table7_intccr_finegray_ic.csv\n")
}

cat("\n--- intccr Fine-Gray interval-censored sensitivity complete. ---\n")
