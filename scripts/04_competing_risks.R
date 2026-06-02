# =====================================================================
# 04_competing_risks.R — competing-risks analysis
#
# Question: Does the radiomic signature predict cumulative incidence of
#           malignant transformation when accounting for the competing
#           risk of death without transformation?
#
# Approach: Fine–Gray subdistribution hazards model on multiply-imputed
#           event times. Frequentist (cmprsk package) since Bayesian
#           Fine–Gray with interval censoring is not well-supported.
#
# Estimated runtime on M4: 5-10 minutes
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})
set.seed(2026)

# Try to install cmprsk if not present
if (!"cmprsk" %in% rownames(installed.packages())) {
  install.packages("cmprsk")
}
library(cmprsk)

# Load datasets
baseline <- readRDS("01_baseline_with_rad.rds")

# Death information (de-identified): dead_event flag + time-to-death in years
# from diagnosis. Shipped as data/death_info.csv so the pipeline runs without
# the patient-identifying Excel file. To regenerate from the raw clinical
# spreadsheet, see scripts/_make_death_info.R (or the inline recipe below).
death_csv_candidates <- c(
  "../data/death_info.csv",        # canonical: run from scripts/
  "data/death_info.csv",           # run from project root
  Sys.getenv("LGG_DEATH_CSV", unset = NA)
)
death_csv_path <- death_csv_candidates[file.exists(death_csv_candidates)][1]

if (!is.na(death_csv_path)) {
  cat("Loading death info from:", death_csv_path, "\n")
  death_lookup <- read.csv(death_csv_path, stringsAsFactors = FALSE)
} else {
  # Fallback: derive from the original Excel if available.
  if (!"readxl" %in% rownames(installed.packages())) install.packages("readxl")
  xlsx_candidates <- c(
    "../source_data/Data_LGG_Clinical_update.xlsx",
    "source_data/Data_LGG_Clinical_update.xlsx",
    "../../source_data/Data_LGG_Clinical_update.xlsx",
    Sys.getenv("LGG_CLINICAL_XLSX", unset = NA)
  )
  xlsx_path <- xlsx_candidates[file.exists(xlsx_candidates)][1]
  if (is.na(xlsx_path)) {
    stop("Cannot find death info.\n",
         "Expected data/death_info.csv (de-identified, shipped with project)\n",
         "OR source_data/Data_LGG_Clinical_update.xlsx.\n",
         "Set LGG_DEATH_CSV or LGG_CLINICAL_XLSX to point to either.")
  }
  cat("death_info.csv missing — deriving from clinical Excel:", xlsx_path, "\n")
  clin_full <- readxl::read_excel(xlsx_path, sheet = "Datensaetzte") %>%
    mutate(
      last_followup = as.Date(last_followup),
      dod = as.Date(dod),
      date_diagnosis = as.Date(date_diagnosis),
      os_status_n = as.numeric(os_status),
      dead_event = as.integer(os_status_n == 2),
      time_to_death_yr = ifelse(dead_event == 1,
                                as.numeric(difftime(dod, date_diagnosis, units = "days")) / 365.25,
                                NA_real_),
      last_followup_yr_from_diag = as.numeric(difftime(last_followup, date_diagnosis, units = "days")) / 365.25
    )
  death_lookup <- clin_full %>%
    select(patient_number, dead_event, time_to_death_yr, last_followup_yr_from_diag)
}

baseline_cr <- baseline %>%
  left_join(death_lookup, by = "patient_number") %>%
  mutate(
    # Event coding for Fine–Gray:
    #   1 = malignant transformation (primary event)
    #   2 = death without prior transformation (competing event)
    #   0 = censored
    fg_event = case_when(
      event == 1 ~ 1,
      dead_event == 1 & event != 1 ~ 2,
      TRUE ~ 0
    ),
    # Time to event (years from diagnosis): midpoint of (L, R) for
    # transformation, time-to-death for competing event, L for censored.
    fg_time = case_when(
      event == 1 ~ (L + R) / 2,
      dead_event == 1 & event != 1 ~ pmax(L, time_to_death_yr, na.rm = TRUE),
      TRUE ~ L
    )
  )

cat("Competing risk event distribution:\n")
print(table(baseline_cr$fg_event, useNA = "ifany"))
cat("  0 = censored, 1 = malignant transformation, 2 = death w/o MT\n\n")

# ---- Fit Fine–Gray for primary event = transformation ----
clinical_covars <- c("contrast_enhancingI", "age_atDiag45", "vol_diag_centered")
sig_covars <- c("original_glszm_largearealowgraylevelemphasis_T1c",
                "original_gldm_dependencevariance_FLAIR")

# Build covariate matrix
X_clin <- model.matrix(~ contrast_enhancingI + age_atDiag45 + vol_diag_centered,
                       data = baseline_cr)[, -1, drop = FALSE]
X_full <- model.matrix(~ contrast_enhancingI + age_atDiag45 + vol_diag_centered +
                         original_glszm_largearealowgraylevelemphasis_T1c +
                         original_gldm_dependencevariance_FLAIR,
                       data = baseline_cr)[, -1, drop = FALSE]

baseline_cr_complete <- baseline_cr %>%
  filter(!is.na(fg_time) & fg_time > 0)

cat("\n=== Fine-Gray model: clinical covariates only ===\n")
fg_clin <- crr(ftime = baseline_cr_complete$fg_time,
               fstatus = baseline_cr_complete$fg_event,
               cov1 = X_clin[match(baseline_cr_complete$patient_number,
                                   baseline_cr$patient_number), ],
               failcode = 1, cencode = 0)
print(summary(fg_clin))

cat("\n=== Fine-Gray model: clinical + 2-feature signature ===\n")
fg_full <- crr(ftime = baseline_cr_complete$fg_time,
               fstatus = baseline_cr_complete$fg_event,
               cov1 = X_full[match(baseline_cr_complete$patient_number,
                                   baseline_cr$patient_number), ],
               failcode = 1, cencode = 0)
print(summary(fg_full))

# Extract sub-distribution hazard ratios with CIs
extract_fg <- function(fit) {
  s <- summary(fit)
  data.frame(
    term = rownames(s$coef),
    sHR = exp(s$coef[, "coef"]),
    sHR_lower = exp(s$coef[, "coef"] - 1.96 * s$coef[, "se(coef)"]),
    sHR_upper = exp(s$coef[, "coef"] + 1.96 * s$coef[, "se(coef)"]),
    p = s$coef[, "p-value"]
  )
}

fg_results <- list(
  clinical_only = extract_fg(fg_clin),
  with_signature = extract_fg(fg_full)
)
saveRDS(fg_results, "04_finegray_results.rds")

# ---- Cumulative incidence by contrast enhancement status ----
cat("\n=== Cumulative incidence: contrast enhancement strata ===\n")
ci_by_ce <- cuminc(ftime = baseline_cr_complete$fg_time,
                   fstatus = baseline_cr_complete$fg_event,
                   group = baseline_cr_complete$contrast_enhancingI)
print(ci_by_ce$Tests)

saveRDS(ci_by_ce, "04_cuminc_by_ce.rds")

cat("\n--- Competing risks analysis complete. Proceed to 05_diagnostics.R ---\n")
