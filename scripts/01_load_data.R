# =====================================================================
# 01_load_data.R — load and verify all analysis datasets
# Run after 00_setup.R completes successfully.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

DATA_DIR <- "../data"

# ---- Load datasets ----
ttm_baseline <- read_csv(file.path(DATA_DIR, "ttm_baseline.csv"), show_col_types = FALSE)
counting_process <- read_csv(file.path(DATA_DIR, "counting_process.csv"), show_col_types = FALSE)
radiomics_99 <- read_csv(file.path(DATA_DIR, "radiomics_99_zscore.csv"), show_col_types = FALSE)
radiomics_428 <- read_csv(file.path(DATA_DIR, "radiomics_428_zscore.csv"), show_col_types = FALSE)

# ---- Sanity checks ----
stopifnot(nrow(ttm_baseline) == 155)
stopifnot(sum(ttm_baseline$event == 1, na.rm = TRUE) == 77)
stopifnot(nrow(counting_process) == 371)
stopifnot(ncol(radiomics_99) == 100)   # 99 features + patient_number
stopifnot(ncol(radiomics_428) == 429)  # 428 features + patient_number

cat("All datasets loaded and verified.\n")
cat("  ttm_baseline:     ", nrow(ttm_baseline), "patients,",
    sum(ttm_baseline$event == 1, na.rm = TRUE), "events\n")
cat("  counting_process: ", nrow(counting_process), "rows from",
    length(unique(counting_process$patient_number)), "patients\n")
cat("  radiomics_99:     99 z-scored features\n")
cat("  radiomics_428:    428 z-scored features\n\n")

# ---- Build merged datasets used downstream ----
# A. Baseline-covariate dataset for primary horseshoe analysis
baseline_with_rad <- ttm_baseline %>%
  left_join(radiomics_99, by = "patient_number") %>%
  mutate(
    # Raw CE column is coded "N"/"Y"; levels + labels chosen so the resulting
    # brms dummy is `contrast_enhancingI1` (matches the prior spec in 02/03/07).
    contrast_enhancingI = factor(contrast_enhancingI, levels = c("N", "Y"), labels = c("0", "1")),
    EOR = factor(EOR, levels = c("GTR", "STR", "Biopsy")),
    idh_factor = factor(case_when(
      idh == 1 ~ "mutant", idh == 0 ~ "wildtype", TRUE ~ NA_character_),
      levels = c("mutant", "wildtype"))
  )

# B. Time-varying dataset (counting-process) — used for TVC analyses
tvc_with_rad <- counting_process %>%
  left_join(radiomics_99, by = "patient_number") %>%
  mutate(
    contrast_enhancingI = factor(contrast_enhancingI, levels = c("N", "Y"), labels = c("0", "1")),
    eor_state = factor(eor_state, levels = c("pre", "GTR", "STR", "Biopsy", "unknown")),
    idh_factor = factor(case_when(
      idh == 1 ~ "mutant", idh == 0 ~ "wildtype", TRUE ~ NA_character_),
      levels = c("mutant", "wildtype"))
  )

# Save
saveRDS(baseline_with_rad, "01_baseline_with_rad.rds")
saveRDS(tvc_with_rad, "01_tvc_with_rad.rds")

cat("Saved: 01_baseline_with_rad.rds and 01_tvc_with_rad.rds\n")
cat("\n--- Data ready. Proceed to 02_bayes_baseline.R ---\n")
