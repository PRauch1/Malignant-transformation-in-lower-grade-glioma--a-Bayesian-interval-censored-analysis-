# =====================================================================
# 08_eor_binary_sensitivity.R — post-hoc EOR sensitivity
#
# Research question: is the surgical-extent effect preserved when STR is
# pooled with Biopsy under a single "non-GTR" label? If GTR-vs-nonGTR
# resembles the 3-level GTR contrast, that supports a binary interpretation
# (any macroscopic residual confers risk). If GTR-vs-nonGTR is weaker than
# GTR-vs-Biopsy, that argues STR is closer to GTR than to Biopsy.
#
# Structure mirrors 03_tvc_cox_mi.R: 20 multiple imputations of interval-
# censored event times with same seeds (2026 + m), frequentist Cox PH on
# counting-process data with id = patient_number for robust SE, pooled via
# Rubin's rules (extended here to pool the full covariance matrix so the
# linear hypothesis GTR − non-GTR can be tested).
#
# Design choice: post_surgery is NOT included in this model, unlike 03's
# primary spec. In 03, post_surgery + eor_state produced an exact linear
# dependency (see RESULTS_REPORT.md §H.5) that made eor_state contrasts
# relative to Biopsy rather than pre. Including post_surgery here would
# reproduce the same rank deficiency (post_surgery ≡ dummy_GTR + dummy_nonGTR)
# and prevent the direct computation of GTR-vs-pre and nonGTR-vs-pre contrasts
# required by the research question. eor_binary_state already encodes the
# post-surgery indicator via its three non-reference levels.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})
set.seed(2026)

baseline <- readRDS("01_baseline_with_rad.rds")
tvc      <- readRDS("01_tvc_with_rad.rds")

# ---- Build binary-collapsed EOR factor on the TVC dataset ----
tvc <- tvc %>%
  mutate(eor_binary_state = factor(
    case_when(
      eor_state %in% c("STR", "Biopsy") ~ "nonGTR",
      TRUE ~ as.character(eor_state)
    ),
    levels = c("pre", "GTR", "nonGTR", "unknown")
  ))

cat("eor_binary_state row counts:\n"); print(table(tvc$eor_binary_state, useNA = "ifany"))
cat("\neor_binary_state unique patients:\n")
print(tvc %>% group_by(eor_binary_state) %>% summarise(n_patients = n_distinct(patient_number)))

# ---- Imputation function (identical to 03's impute_one) ----
impute_one <- function(seed) {
  set.seed(seed)
  imp <- baseline %>%
    mutate(
      imputed_time = ifelse(is.infinite(R), L, runif(n(), pmax(L, 1e-4), R)),
      imputed_event = ifelse(is.infinite(R), 0, 1)
    )
  cp <- list()
  for (i in seq_len(nrow(imp))) {
    pid <- imp$patient_number[i]
    et  <- imp$imputed_time[i]
    ev  <- imp$imputed_event[i]
    pat_rows <- tvc %>% filter(patient_number == pid) %>% arrange(tstart)
    pat_rows <- pat_rows %>%
      filter(tstart < et) %>%
      mutate(tstop_new = pmin(tstop, et), event_in_row = 0)
    if (ev == 1 && nrow(pat_rows) > 0) pat_rows$event_in_row[nrow(pat_rows)] <- 1
    pat_rows <- pat_rows %>% select(-tstop) %>% rename(tstop = tstop_new)
    pat_rows <- pat_rows %>% filter(tstop > tstart)
    cp[[length(cp) + 1]] <- pat_rows
  }
  bind_rows(cp)
}

# ---- Rubin pooling with full covariance (extended beyond 03's pool_rubin) ----
pool_rubin_vcov <- function(mi_results) {
  mi_results <- mi_results[!sapply(mi_results, is.null)]
  M <- length(mi_results)
  if (M < 2) stop("Need at least 2 successful imputations to pool.")
  all_names <- unique(unlist(lapply(mi_results, function(x) names(x$coef))))
  coefs_mat <- do.call(rbind, lapply(mi_results, function(x) x$coef[all_names]))
  # Within-imputation covariance matrix averaged across M
  vcov_arr <- array(NA_real_,
                    dim = c(length(all_names), length(all_names), M),
                    dimnames = list(all_names, all_names, NULL))
  for (m in seq_len(M)) {
    v <- mi_results[[m]]$vcov
    common <- intersect(rownames(v), all_names)
    vcov_arr[common, common, m] <- v[common, common]
  }
  V_within  <- apply(vcov_arr, c(1, 2), mean, na.rm = TRUE)
  # Between-imputation covariance
  V_between <- cov(coefs_mat, use = "pairwise.complete.obs")
  V_between[is.na(V_between)] <- 0
  V_total   <- V_within + (1 + 1 / M) * V_between
  pooled_est <- colMeans(coefs_mat, na.rm = TRUE)
  pooled_se  <- sqrt(diag(V_total))
  list(est = pooled_est, se = pooled_se, vcov = V_total, M = M,
       M_attempted = NA_integer_)
}

# ---- Model formula: drop post_surgery (see top-of-file note) ----
clin_terms <- c("contrast_enhancingI", "age_atDiag45", "vol_diag_centered")
tvc_terms  <- c("eor_binary_state", "post_chemo", "post_radiation")
f_obj <- as.formula(paste0(
  "Surv(tstart, tstop, event_in_row) ~ ",
  paste(c(clin_terms, tvc_terms), collapse = " + ")
))

M <- 20

fit_cohort_mi <- function(subset_fn, label) {
  results <- vector("list", M)
  for (m in 1:M) {
    imp_m <- impute_one(seed = 2026 + m)
    imp_m <- subset_fn(imp_m)
    if (nrow(imp_m) == 0) next
    fit_m <- try(coxph(f_obj, data = imp_m, id = patient_number, robust = TRUE),
                 silent = TRUE)
    if (!inherits(fit_m, "try-error")) {
      results[[m]] <- list(coef = coef(fit_m), vcov = vcov(fit_m))
    }
    if (m %% 5 == 0) cat("  [", label, "] imputation", m, "/", M, "\n")
  }
  results
}

cat("\n=== Full cohort MI ===\n")
mi_full <- fit_cohort_mi(identity, "full")
n_full_success <- sum(!sapply(mi_full, is.null))
cat("Full cohort successful imputations:", n_full_success, "/", M, "\n")
pooled_full <- pool_rubin_vcov(mi_full)
pooled_full$M_attempted <- M

cat("\n=== IDH-mutant subgroup MI ===\n")
mi_mut <- fit_cohort_mi(function(d) d %>% filter(idh == 1, !is.na(idh)), "IDH-mut")
n_mut_success <- sum(!sapply(mi_mut, is.null))
cat("IDH-mutant successful imputations:", n_mut_success, "/", M, "\n")
pooled_mut <- pool_rubin_vcov(mi_mut)
pooled_mut$M_attempted <- M

# ---- Build per-coefficient pooled table ----
build_table <- function(pooled, label) {
  est <- pooled$est; se <- pooled$se
  z <- est / se
  data.frame(
    cohort   = label,
    term     = names(est),
    estimate = est,
    HR       = exp(est),
    se       = se,
    HR_lower = exp(est - 1.96 * se),
    HR_upper = exp(est + 1.96 * se),
    z        = z,
    p        = 2 * pnorm(-abs(z)),
    row.names = NULL
  )
}

# ---- Linear hypothesis: GTR vs non-GTR ----
gtr_vs_nongtr <- function(pooled, label) {
  est <- pooled$est; V <- pooled$vcov
  i_G  <- "eor_binary_stateGTR"
  i_nG <- "eor_binary_statenonGTR"
  if (!all(c(i_G, i_nG) %in% names(est))) {
    return(data.frame(cohort = label, term = "GTR_vs_nonGTR",
                      estimate = NA_real_, HR = NA_real_, se = NA_real_,
                      HR_lower = NA_real_, HR_upper = NA_real_,
                      z = NA_real_, p = NA_real_))
  }
  d  <- est[i_G] - est[i_nG]
  vd <- V[i_G, i_G] + V[i_nG, i_nG] - 2 * V[i_G, i_nG]
  sd <- sqrt(vd)
  data.frame(
    cohort = label,
    term = "GTR_vs_nonGTR",
    estimate = unname(d),
    HR = exp(unname(d)),
    se = unname(sd),
    HR_lower = exp(unname(d - 1.96 * sd)),
    HR_upper = exp(unname(d + 1.96 * sd)),
    z = unname(d / sd),
    p = unname(2 * pnorm(-abs(d / sd)))
  )
}

out_full <- rbind(build_table(pooled_full, "Full"),
                  gtr_vs_nongtr(pooled_full, "Full"))
out_mut  <- rbind(build_table(pooled_mut,  "IDH-mutant"),
                  gtr_vs_nongtr(pooled_mut, "IDH-mutant"))
all_out <- rbind(out_full, out_mut)

# ---- Print headline rows ----
headline_rows <- c("eor_binary_stateGTR", "eor_binary_statenonGTR",
                   "eor_binary_stateunknown", "GTR_vs_nonGTR")
cat("\n=== EOR BINARY SENSITIVITY: full cohort (headline) ===\n")
print(subset(out_full, term %in% headline_rows), row.names = FALSE, digits = 4)
cat("\n=== EOR BINARY SENSITIVITY: full cohort (full table) ===\n")
print(out_full, row.names = FALSE, digits = 4)

cat("\n=== EOR BINARY SENSITIVITY: IDH-mutant (headline) ===\n")
print(subset(out_mut, term %in% headline_rows), row.names = FALSE, digits = 4)
cat("\n=== EOR BINARY SENSITIVITY: IDH-mutant (full table) ===\n")
print(out_mut, row.names = FALSE, digits = 4)

# ---- Save ----
saveRDS(list(
  full        = pooled_full,
  mut         = pooled_mut,
  table_full  = out_full,
  table_mut   = out_mut,
  M           = M,
  n_success   = list(full = n_full_success, mut = n_mut_success)
), "08_eor_binary_results.rds")

if (!dir.exists("final_summary")) dir.create("final_summary", showWarnings = FALSE)
write.csv(all_out, "final_summary/table5_eor_binary_sensitivity.csv", row.names = FALSE)

cat("\n--- Saved: 08_eor_binary_results.rds and final_summary/table5_eor_binary_sensitivity.csv ---\n")
