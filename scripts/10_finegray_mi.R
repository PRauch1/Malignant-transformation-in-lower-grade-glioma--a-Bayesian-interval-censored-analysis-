# =====================================================================
# 10_finegray_mi.R — Fine-Gray competing risks with multiple imputation
#
# Addresses statistician review (Baran, May 2026): script 04
# uses (L + R) / 2 as the transformation event time, which is a
# midpoint plug-in for interval censoring. This script replicates 04
# but draws M = 20 event times t ~ Uniform(L, R) per Rubin's-rules
# imputation, fits Fine-Gray on each, and pools sHR + SE.
#
# Seeds match 03_tvc_cox_mi.R (2026 + m) so the same imputed event
# times are used across the TVC Cox and Fine-Gray analyses.
#
# Estimated runtime on M4: 3-6 minutes
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(cmprsk)
  library(readxl)
})

baseline <- readRDS("01_baseline_with_rad.rds")

# ---- Load death information (same logic as 04_competing_risks.R) ----
death_csv_candidates <- c(
  "../data/death_info.csv",
  "data/death_info.csv",
  Sys.getenv("LGG_DEATH_CSV", unset = NA)
)
death_csv_path <- death_csv_candidates[file.exists(death_csv_candidates)][1]

if (!is.na(death_csv_path)) {
  cat("Loading death info from:", death_csv_path, "\n")
  death_lookup <- read.csv(death_csv_path, stringsAsFactors = FALSE)
} else {
  xlsx_candidates <- c(
    "../source_data/Data_LGG_Clinical_update.xlsx",
    "source_data/Data_LGG_Clinical_update.xlsx",
    "../../source_data/Data_LGG_Clinical_update.xlsx",
    Sys.getenv("LGG_CLINICAL_XLSX", unset = NA)
  )
  xlsx_path <- xlsx_candidates[file.exists(xlsx_candidates)][1]
  if (is.na(xlsx_path)) {
    stop("Cannot find death info. Set LGG_DEATH_CSV or LGG_CLINICAL_XLSX.")
  }
  clin_full <- readxl::read_excel(xlsx_path, sheet = "Datensaetzte") %>%
    mutate(
      last_followup = as.Date(last_followup),
      dod = as.Date(dod),
      date_diagnosis = as.Date(date_diagnosis),
      os_status_n = as.numeric(os_status),
      dead_event = as.integer(os_status_n == 2),
      time_to_death_yr = ifelse(dead_event == 1,
                                as.numeric(difftime(dod, date_diagnosis,
                                                    units = "days")) / 365.25,
                                NA_real_),
      last_followup_yr_from_diag = as.numeric(difftime(last_followup,
                                                       date_diagnosis,
                                                       units = "days")) / 365.25
    )
  death_lookup <- clin_full %>%
    select(patient_number, dead_event, time_to_death_yr,
           last_followup_yr_from_diag)
}

baseline_cr <- baseline %>% left_join(death_lookup, by = "patient_number")

# ---- Build covariate matrices (constant across imputations) ----
X_clin <- model.matrix(~ contrast_enhancingI + age_atDiag45 + vol_diag_centered,
                       data = baseline_cr)[, -1, drop = FALSE]
X_full <- model.matrix(~ contrast_enhancingI + age_atDiag45 + vol_diag_centered +
                         original_glszm_largearealowgraylevelemphasis_T1c +
                         original_gldm_dependencevariance_FLAIR,
                       data = baseline_cr)[, -1, drop = FALSE]

# ---- Single imputation: draw transformation time t ~ Uniform(L, R) ----
# For censored (R = Inf): use L as right-censoring time.
# For competing death (event != 1 & dead_event == 1): use time_to_death_yr.
impute_fg <- function(seed) {
  set.seed(seed)
  baseline_cr %>%
    mutate(
      fg_event = case_when(
        event == 1 ~ 1,
        dead_event == 1 & event != 1 ~ 2,
        TRUE ~ 0
      ),
      fg_time = case_when(
        event == 1 ~ runif(n(), pmax(L, 1e-4), R),
        dead_event == 1 & event != 1 ~ pmax(L, time_to_death_yr, na.rm = TRUE),
        TRUE ~ L
      )
    ) %>%
    filter(!is.na(fg_time) & fg_time > 0)
}

# ---- Fit Fine-Gray on one imputation and return coef + se ----
fit_fg <- function(dat, X) {
  idx <- match(dat$patient_number, baseline_cr$patient_number)
  fit <- crr(ftime = dat$fg_time, fstatus = dat$fg_event,
             cov1 = X[idx, , drop = FALSE], failcode = 1, cencode = 0)
  list(coef = fit$coef, se = sqrt(diag(fit$var)),
       names = colnames(X))
}

# ---- Rubin pooling ----
pool_rubin_fg <- function(mi_list) {
  coefs_mat <- do.call(rbind, lapply(mi_list, function(x) x$coef))
  ses_mat   <- do.call(rbind, lapply(mi_list, function(x) x$se))
  M <- nrow(coefs_mat)
  pooled_est <- colMeans(coefs_mat)
  within_var <- colMeans(ses_mat^2)
  between_var <- apply(coefs_mat, 2, var)
  total_var <- within_var + (1 + 1/M) * between_var
  pooled_se <- sqrt(total_var)
  z <- pooled_est / pooled_se
  p <- 2 * pnorm(-abs(z))
  data.frame(
    term = mi_list[[1]]$names,
    sHR = exp(pooled_est),
    sHR_lower = exp(pooled_est - 1.96 * pooled_se),
    sHR_upper = exp(pooled_est + 1.96 * pooled_se),
    se_log_sHR = pooled_se,
    z = z, p = p
  )
}

# ---- Run M = 20 imputations for both Fine-Gray specifications ----
M <- 20
cat("=== Fine-Gray with multiple imputation (M =", M, ") ===\n")
cat("Imputing transformation times t ~ Uniform(L, R) per patient.\n")
cat("Seeds 2026+m, matching 03_tvc_cox_mi.R for cross-analysis consistency.\n\n")

mi_clin <- list()
mi_full <- list()
for (m in 1:M) {
  dat_m <- impute_fg(seed = 2026 + m)
  mi_clin[[m]] <- fit_fg(dat_m, X_clin)
  mi_full[[m]] <- fit_fg(dat_m, X_full)
  if (m %% 5 == 0) cat("  Done imputation", m, "/", M, "\n")
}

pooled_clin <- pool_rubin_fg(mi_clin)
pooled_full <- pool_rubin_fg(mi_full)

cat("\n=== POOLED Fine-Gray: clinical covariates only ===\n")
print(pooled_clin, row.names = FALSE, digits = 4)

cat("\n=== POOLED Fine-Gray: clinical + 2-feature signature ===\n")
print(pooled_full, row.names = FALSE, digits = 4)

saveRDS(list(clinical_only = pooled_clin, with_signature = pooled_full,
             M = M, seed_base = 2026),
        "10_finegray_mi_results.rds")

# ---- Comparison against midpoint (04_competing_risks.R) ----
cat("\n=== Comparison: MI vs midpoint plug-in (sHR for CE) ===\n")
fg_04 <- readRDS("04_finegray_results.rds")
ce_04_clin <- fg_04$clinical_only[fg_04$clinical_only$term == "contrast_enhancingIY", ]
if (nrow(ce_04_clin) == 0) {
  ce_04_clin <- fg_04$clinical_only[grep("contrast", fg_04$clinical_only$term), ]
}
ce_10_clin <- pooled_clin[grep("contrast", pooled_clin$term), ]
cat(sprintf("  04 midpoint (clin):  sHR = %.3f (%.3f, %.3f), p = %.4f\n",
            ce_04_clin$sHR, ce_04_clin$sHR_lower, ce_04_clin$sHR_upper, ce_04_clin$p))
cat(sprintf("  10 MI (clin):        sHR = %.3f (%.3f, %.3f), p = %.4f\n",
            ce_10_clin$sHR, ce_10_clin$sHR_lower, ce_10_clin$sHR_upper, ce_10_clin$p))

# ---- Final summary table ----
if (dir.exists("final_summary")) {
  write.csv(pooled_clin,
            "final_summary/table6_finegray_mi_clinical.csv", row.names = FALSE)
  write.csv(pooled_full,
            "final_summary/table6_finegray_mi_with_signature.csv", row.names = FALSE)
  cat("\nSaved: final_summary/table6_finegray_mi_*.csv\n")
}

cat("\n--- Fine-Gray MI complete. Proceed to 11_smoothhazard_cr.R for the\n",
    "    interval-censoring-aware cause-specific-hazards sensitivity.\n",
    sep = "")
