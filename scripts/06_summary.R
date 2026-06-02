# =====================================================================
# 06_summary.R — final consolidated results for the manuscript
#
# Produces: clean tables ready for paper, summary text describing the
# headline findings, and saves everything to a final_summary directory.
#
# Estimated runtime on M4: < 1 minute
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(brms)
  library(posterior)
})

# Load all results
fit_full <- readRDS("02_fit_full.rds")
fit_mut  <- readRDS("02_fit_mut.rds")
fit_wt   <- if (file.exists("02_fit_wt.rds")) readRDS("02_fit_wt.rds") else NULL
clin_table <- readRDS("05_clinical_effects.rds")
rad_post  <- readRDS("05_radiomic_posterior.rds")
tvc_pooled <- readRDS("03_tvc_summary.rds")
fg_results <- readRDS("04_finegray_results.rds")

dir.create("final_summary", showWarnings = FALSE)

# ===== TABLE 1: clinical effects across cohorts =====
table1 <- clin_table %>%
  mutate(
    HR_CrI = sprintf("%.2f (%.2f, %.2f)", HR_median, HR_q025, HR_q975),
    Pr_HR_gt_1 = sprintf("%.3f", Pr_HR_gt_1)
  ) %>%
  select(cohort, parameter, HR_CrI, Pr_HR_gt_1)
write.csv(table1, "final_summary/table1_clinical_effects.csv", row.names = FALSE)
cat("=== TABLE 1: Clinical predictors (Bayesian) ===\n")
print(table1, row.names = FALSE)

# ===== TABLE 2: top radiomic features by posterior credibility =====
top_features <- function(rad_df, n = 10) {
  rad_df %>%
    mutate(
      HR_CrI = sprintf("%.2f (%.2f, %.2f)", HR_median, HR_q025, HR_q975),
      Pr_active = sprintf("%.2f", pmax(Pr_HR_gt_1.2, Pr_HR_lt_0.83))
    ) %>%
    select(cohort, feature, HR_CrI, Pr_active, ci_excludes_0) %>%
    head(n)
}

table2_full <- top_features(rad_post$full, n = 15)
table2_mut  <- top_features(rad_post$mut,  n = 15)

write.csv(table2_full, "final_summary/table2_radiomic_full.csv", row.names = FALSE)
write.csv(table2_mut,  "final_summary/table2_radiomic_mut.csv",  row.names = FALSE)
cat("\n=== TABLE 2: Top 15 radiomic features (Bayesian, full cohort) ===\n")
print(table2_full, row.names = FALSE)
cat("\n=== TABLE 2b: Top 15 radiomic features (Bayesian, IDH-mutant) ===\n")
print(table2_mut, row.names = FALSE)

# ===== TABLE 3: time-varying treatment effects (frequentist MI Cox) =====
table3 <- tvc_pooled$full %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, HR_lower, HR_upper),
    p_value = sprintf("%.4f", p)
  ) %>%
  select(term, HR_CI, p_value)
write.csv(table3, "final_summary/table3_tvc_treatment_effects.csv", row.names = FALSE)
cat("\n=== TABLE 3: Time-varying treatment effects (full cohort, MI Cox, step-function TVC) ===\n")
print(table3, row.names = FALSE)

# Extended TVC: time-since-treatment
if (!is.null(tvc_pooled$modelA_full)) {
  table3A <- tvc_pooled$modelA_full %>%
    mutate(HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, HR_lower, HR_upper),
           p_value = sprintf("%.4f", p)) %>%
    select(term, HR_CI, p_value)
  write.csv(table3A, "final_summary/table3A_tvc_time_since_treatment.csv", row.names = FALSE)
  cat("\n=== TABLE 3A: TVC extended with time-since-treatment terms ===\n")
  cat("(Tests whether the treatment effect evolves with time since treatment)\n")
  print(table3A, row.names = FALSE)
}

# Extended TVC: time-to-treatment as baseline covariate
if (!is.null(tvc_pooled$modelB_full)) {
  table3B <- tvc_pooled$modelB_full %>%
    mutate(HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, HR_lower, HR_upper),
           p_value = sprintf("%.4f", p)) %>%
    select(term, HR_CI, p_value)
  write.csv(table3B, "final_summary/table3B_tvc_time_to_treatment.csv", row.names = FALSE)
  cat("\n=== TABLE 3B: TVC extended with time-to-treatment baseline covariates ===\n")
  cat("(Tests whether patients treated earlier vs later differ in transformation hazard)\n")
  print(table3B, row.names = FALSE)
}

table3_mut <- tvc_pooled$mut %>%
  mutate(
    HR_CI = sprintf("%.2f (%.2f, %.2f)", HR, HR_lower, HR_upper),
    p_value = sprintf("%.4f", p)
  ) %>%
  select(term, HR_CI, p_value)
write.csv(table3_mut, "final_summary/table3b_tvc_idhmut.csv", row.names = FALSE)
cat("\n=== TABLE 3b: TVC (IDH-mutant) ===\n")
print(table3_mut, row.names = FALSE)

# ===== TABLE 4: Fine-Gray subdistribution hazards =====
table4_clin <- fg_results$clinical_only %>%
  mutate(sHR_CI = sprintf("%.2f (%.2f, %.2f)", sHR, sHR_lower, sHR_upper),
         p = sprintf("%.4f", p)) %>%
  select(term, sHR_CI, p)
table4_full <- fg_results$with_signature %>%
  mutate(sHR_CI = sprintf("%.2f (%.2f, %.2f)", sHR, sHR_lower, sHR_upper),
         p = sprintf("%.4f", p)) %>%
  select(term, sHR_CI, p)
write.csv(table4_clin, "final_summary/table4_finegray_clinical.csv", row.names = FALSE)
write.csv(table4_full, "final_summary/table4_finegray_with_signature.csv", row.names = FALSE)
cat("\n=== TABLE 4: Fine-Gray sHR (clinical only) ===\n")
print(table4_clin, row.names = FALSE)
cat("\n=== TABLE 4b: Fine-Gray sHR (with 2-feature signature) ===\n")
print(table4_full, row.names = FALSE)

# ===== Headline narrative =====
cat("\n\n")
cat("====================================================================\n")
cat("                        HEADLINE FINDINGS\n")
cat("====================================================================\n")

ce_full <- clin_table %>% filter(cohort == "Full", parameter == "b_contrast_enhancingI1")
ce_mut  <- clin_table %>% filter(cohort == "IDH-mutant", parameter == "b_contrast_enhancingI1")
cat(sprintf(
  "1. Contrast enhancement at baseline: HR %.2f (95%% CrI %.2f-%.2f) full cohort,\n   HR %.2f (95%% CrI %.2f-%.2f) within IDH-mutant.\n",
  ce_full$HR_median, ce_full$HR_q025, ce_full$HR_q975,
  ce_mut$HR_median, ce_mut$HR_q025, ce_mut$HR_q975))

n_active_full <- sum(rad_post$full$ci_excludes_0)
n_active_mut  <- sum(rad_post$mut$ci_excludes_0)
cat(sprintf(
  "2. Radiomic features with 95%% CrI excluding zero:\n   Full cohort: %d / %d   IDH-mutant: %d / %d\n",
  n_active_full, nrow(rad_post$full), n_active_mut, nrow(rad_post$mut)))

cat("\n3. Top 5 radiomic features in full cohort by |posterior median log-HR|:\n")
print(head(rad_post$full[, c("feature", "HR_median", "HR_q025", "HR_q975", "ci_excludes_0")], 5),
      row.names = FALSE)

cat("\n4. Time-varying treatment effects (full cohort, pooled MI Cox):\n")
print(table3, row.names = FALSE)

cat("\n5. Fine-Gray (competing risk = death without MT):\n")
print(table4_full, row.names = FALSE)

cat("\n--- All summary tables saved to final_summary/ ---\n")
cat("--- Send these tables back to Claude for the manuscript rewrite. ---\n")
