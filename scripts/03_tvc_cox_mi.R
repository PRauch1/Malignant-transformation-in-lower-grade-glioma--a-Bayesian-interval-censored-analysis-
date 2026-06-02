# =====================================================================
# 03_tvc_cox_mi.R — time-varying-covariate Cox via multiple imputation
#
# IMPORTANT methodological note: this script is NOT fully Bayesian.
# Bayesian time-varying-covariate survival with interval censoring is
# not well-supported in current R/Stan tooling. We use a frequentist
# Cox proportional hazards model with multiply-imputed event times
# (M=20 uniform draws within each observation interval), pooling
# coefficients via Rubin's rules. The Bayesian primary analysis in
# 02_bayes_baseline.R uses BASELINE covariates only.
#
# This is a deliberate methodological compromise that is honestly
# documented in the analysis plan and should be reflected in the
# manuscript as such.
#
# Estimated runtime on M4: 10-20 minutes
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(brms)
  library(posterior)
  library(survival)
})
options(brms.backend = "cmdstanr")
options(mc.cores = parallel::detectCores() - 1)
set.seed(2026)

baseline <- readRDS("01_baseline_with_rad.rds")
tvc <- readRDS("01_tvc_with_rad.rds")

# ---- Multiple imputation of event times within observation intervals ----
# For each imputation:
#   - For event patients (R < Inf): draw event time uniformly from (L, R]
#   - For censored patients (R == Inf): use L as censoring time
#   - Then update the counting-process dataset to:
#     (a) close out the patient's row sequence at the imputed event time
#     (b) mark only the LAST row as event = 1
M <- 20
cat("Setting up", M, "multiple imputations of event times...\n")

impute_one <- function(seed) {
  set.seed(seed)
  # Per-patient imputed event time
  pids <- baseline$patient_number
  imp <- baseline %>%
    mutate(
      imputed_time = ifelse(is.infinite(R),
                             L,
                             runif(n(), pmax(L, 1e-4), R)),
      imputed_event = ifelse(is.infinite(R), 0, 1)
    )

  # Build a counting-process dataset whose rows END at imputed_time
  # for each patient
  cp <- list()
  for (i in seq_len(nrow(imp))) {
    pid <- imp$patient_number[i]
    et  <- imp$imputed_time[i]
    ev  <- imp$imputed_event[i]
    pat_rows <- tvc %>% filter(patient_number == pid) %>% arrange(tstart)
    # Truncate any row that extends past et
    pat_rows <- pat_rows %>%
      filter(tstart < et) %>%
      mutate(tstop_new = pmin(tstop, et)) %>%
      mutate(event_in_row = 0)
    # Mark the last row as event if applicable
    if (ev == 1 && nrow(pat_rows) > 0) {
      pat_rows$event_in_row[nrow(pat_rows)] <- 1
    }
    pat_rows <- pat_rows %>% select(-tstop) %>% rename(tstop = tstop_new)
    # Filter zero-length rows
    pat_rows <- pat_rows %>% filter(tstop > tstart)
    cp[[length(cp) + 1]] <- pat_rows
  }
  bind_rows(cp)
}

# ---- Sanity check on a single imputation ----
imp1 <- impute_one(1)
cat("Single imputation: ", nrow(imp1), "rows from",
    length(unique(imp1$patient_number)), "patients,",
    sum(imp1$event_in_row), "events\n")

# ---- Define brms model formula for TVC ----
# Use Cox-style: time2 - time1 with cens marker
# brms supports survival via family = cox() requiring exact event times — perfect for MI
# All time-varying covariates (post_surgery, post_chemo, post_radiation, eor_state)
# get standard normal-ish weakly informative priors
# Radiomic covariates get horseshoe via the unspecified-b default
rad_cols <- grep("^original_", colnames(tvc), value = TRUE)

# Two model variants:
#   M_tvc_basic:    clinical baseline + TVC treatment effects, NO radiomics
#   M_tvc_full:     same + radiomics with horseshoe

clin_terms <- c("contrast_enhancingI", "age_atDiag45", "vol_diag_centered")
tvc_terms  <- c("post_surgery", "eor_state", "post_chemo", "post_radiation")
all_terms_basic <- c(clin_terms, tvc_terms)
all_terms_full  <- c(clin_terms, tvc_terms, rad_cols)

# Use Cox PH via brms (which calls Stan). This requires exact event times per row.
# brms cox model formula: time | cens(censoring_indicator) ~ covariates
# We need to set: censoring = 1 - event_in_row (brms cens convention is "censored = 1 means observation is right-censored")
# Actually brms cens function takes a vector of strings: "none", "right", etc., OR an event indicator with cens="right"

priors_basic <- c(
  prior(student_t(3, 0, 2.5), class = "b", coef = "contrast_enhancingI1"),
  prior(student_t(3, 0, 0.1), class = "b", coef = "age_atDiag45"),
  prior(student_t(3, 0, 0.05), class = "b", coef = "vol_diag_centered"),
  prior(normal(0, 1.5), class = "b", coef = "post_surgery"),
  prior(normal(0, 1.5), class = "b", coef = "eor_stateGTR"),
  prior(normal(0, 1.5), class = "b", coef = "eor_stateSTR"),
  prior(normal(0, 1.5), class = "b", coef = "eor_stateBiopsy"),
  prior(normal(0, 1.5), class = "b", coef = "post_chemo"),
  prior(normal(0, 1.5), class = "b", coef = "post_radiation")
)

priors_full <- c(
  priors_basic,
  prior(horseshoe(par_ratio = 0.1, df = 1, df_global = 1,
                  scale_global = 1, scale_slab = 2, df_slab = 4),
        class = "b")
)

# ---- Fit MI loop ----
fit_one_imputation <- function(imp_data, formula_terms, priors_obj,
                                seed, model_label) {
  imp_data$cens_flag <- ifelse(imp_data$event_in_row == 1, "none", "right")
  imp_data$time_for_brms <- imp_data$tstop  # Cox uses event/censoring time
  # For counting-process Cox in brms, we need entry times via cumulative time approach
  # Actually brms cox doesn't natively handle (tstart, tstop) format. Workaround:
  # use survival::coxph for the actual TVC fit and store the coefs, OR
  # collapse to per-patient summaries using the "last interval" treatment status (less ideal)
  #
  # Best approach for honest TVC: fit a frequentist Cox via survival::coxph with
  # (tstart, tstop, event_in_row) and get coefs + bootstrap SEs. We then report these
  # as a sensitivity analysis to the brms baseline model. The Bayesian TVC piece is
  # a methodological compromise here.
  NULL  # placeholder
}

# ---- Pragmatic compromise: fit TVC via survival::coxph on each imputation ----
# This gives us proper TVC estimation, just without horseshoe shrinkage on the
# treatment terms (which are pre-specified and few, so shrinkage is less critical there).
# We average coefficients across MI replicates per Rubin's rules.

cat("\n=== TVC Cox analysis with multiple imputation (M=", M, ") ===\n")
mi_results <- list()
formula_str <- paste0("Surv(tstart, tstop, event_in_row) ~ ",
                      paste(c(clin_terms, tvc_terms), collapse = " + "))
formula_obj <- as.formula(formula_str)

for (m in 1:M) {
  imp_m <- impute_one(seed = 2026 + m)
  fit_m <- try(coxph(formula_obj, data = imp_m, id = patient_number, robust = TRUE), silent = TRUE)
  if (!inherits(fit_m, "try-error")) {
    co <- coef(fit_m)
    se <- sqrt(diag(vcov(fit_m)))
    mi_results[[m]] <- list(coef = co, se = se, fit = fit_m)
  }
  if (m %% 5 == 0) cat("  Done imputation", m, "/", M, "\n")
}

# ---- Pool via Rubin's rules ----
pool_rubin <- function(mi_results) {
  coefs_mat <- do.call(rbind, lapply(mi_results, function(x) x$coef))
  ses_mat   <- do.call(rbind, lapply(mi_results, function(x) x$se))
  M <- nrow(coefs_mat)
  pooled_est <- colMeans(coefs_mat)
  within_var <- colMeans(ses_mat^2)
  between_var <- apply(coefs_mat, 2, var)
  total_var <- within_var + (1 + 1/M) * between_var
  pooled_se <- sqrt(total_var)
  z <- pooled_est / pooled_se
  p <- 2 * pnorm(-abs(z))
  data.frame(
    term = colnames(coefs_mat),
    estimate = pooled_est,
    HR = exp(pooled_est),
    se = pooled_se,
    HR_lower = exp(pooled_est - 1.96 * pooled_se),
    HR_upper = exp(pooled_est + 1.96 * pooled_se),
    z = z, p = p
  )
}

pooled_full <- pool_rubin(mi_results)
cat("\n=== POOLED TVC Cox results (Full cohort, M=", length(mi_results), ") ===\n")
if (length(mi_results) < M) {
  cat("WARNING:", M - length(mi_results), "imputations failed and were dropped.\n")
  cat("If failure rate > 20%, the model is likely misspecified for some imputed datasets.\n")
}
print(pooled_full, row.names = FALSE)
saveRDS(pooled_full, "03_pooled_tvc_full.rds")

# ---- Extended models with time-since-treatment and time-to-treatment effects ----
# These address whether WHEN a treatment happened (not just whether it did) matters.
#
# Extended model A: time-since-surgery as continuous modifier
#   - Adds `t_since_surgery` = max(0, tstart - t_surgery_per_patient) to each row
#   - Interacted with post_surgery: post_surgery * t_since_surgery
#   - Captures waning / accumulating post-surgical effect over time
#
# Extended model B: time-to-first-treatment as baseline covariate
#   - Adds per-patient `t_to_surgery_at_baseline` and `t_to_chemo_at_baseline`
#     as baseline (time-0) covariates, to capture selection effects
#     (longer delay to treatment indicates clinically stable tumor)

cat("\n=== Extended Model A: time-since-surgery as continuous modifier ===\n")
# For each counting-process row, add t_since_surgery
# Patient's t_surgery (years from diagnosis) is already a per-patient column
impute_with_time_since <- function(seed) {
  cp <- impute_one(seed)
  cp <- cp %>%
    mutate(
      t_since_surgery = ifelse(post_surgery == 1 & !is.na(t_surgery),
                                pmax(tstart - t_surgery, 0), 0),
      t_since_chemo = ifelse(post_chemo == 1 & !is.na(t_chemo_first),
                              pmax(tstart - t_chemo_first, 0), 0),
      t_since_radiation = ifelse(post_radiation == 1 & !is.na(t_radiation_first),
                                  pmax(tstart - t_radiation_first, 0), 0)
    )
  cp
}

mi_results_A <- list()
formula_A <- as.formula(paste0(
  "Surv(tstart, tstop, event_in_row) ~ ",
  paste(c(clin_terms, tvc_terms), collapse = " + "),
  " + t_since_surgery + t_since_chemo + t_since_radiation"
))
for (m in 1:M) {
  imp_m <- impute_with_time_since(seed = 2026 + m)
  fit_m <- try(coxph(formula_A, data = imp_m, id = patient_number, robust = TRUE), silent = TRUE)
  if (!inherits(fit_m, "try-error")) {
    mi_results_A[[m]] <- list(coef = coef(fit_m), se = sqrt(diag(vcov(fit_m))))
  }
}
pooled_A <- pool_rubin(mi_results_A)
cat("\n--- Model A: with time-since-treatment terms ---\n")
print(pooled_A, row.names = FALSE)
saveRDS(pooled_A, "03_pooled_tvc_modelA.rds")

# ---------------------------------------------------------------------
# Model B REMOVED following statistician review (Baran, May 2026).
# Rationale: Model B used `t_to_surgery_at_baseline` (the per-patient
# time from diagnosis to first surgery) as a baseline covariate carried
# into every counting-process interval. For a patient operated at year
# 2, this value is "2" from day 0 onward — i.e., the regressor is not
# yet observable at the time the hazard is being modeled. This is
# anticipation bias of the regressor and violates the principle that
# covariates must be measurable at the time being modeled. The proper
# alternative is a time-varying `t_since_surgery` (Model A above) or a
# landmark analysis; Model B as written is not a valid prognostic
# regression. The reported HR for `t_to_surgery_at_baseline` was 1.03
# (0.94, 1.14), p = 0.53 (null), so its removal changes no substantive
# conclusion.
# ---------------------------------------------------------------------

# ---- Repeat for IDH-mutant subgroup ----
cat("\n=== IDH-mutant subgroup TVC analysis ===\n")
mi_results_mut <- list()
for (m in 1:M) {
  imp_m <- impute_one(seed = 2026 + m)
  imp_m_mut <- imp_m %>% filter(idh == 1, !is.na(idh))
  fit_m <- try(coxph(formula_obj, data = imp_m_mut, id = patient_number, robust = TRUE), silent = TRUE)
  if (!inherits(fit_m, "try-error")) {
    mi_results_mut[[m]] <- list(coef = coef(fit_m), se = sqrt(diag(vcov(fit_m))))
  }
}
pooled_mut <- pool_rubin(mi_results_mut)
cat("\n--- IDH-mutant pooled results ---\n")
print(pooled_mut, row.names = FALSE)
saveRDS(pooled_mut, "03_pooled_tvc_mut.rds")

# ---- Save imputation diagnostics ----
saveRDS(list(M = M,
             full = pooled_full, mut = pooled_mut,
             modelA_full = pooled_A),
        "03_tvc_summary.rds")

cat("\n--- TVC analysis complete. Proceed to 04_competing_risks.R ---\n")
