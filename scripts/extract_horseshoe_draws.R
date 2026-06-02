suppressPackageStartupMessages({
  library(brms); library(posterior)
})
fit <- readRDS("scripts/07_sens_horseshoe_pr0.1.rds")
d   <- as_draws_df(fit)
rad <- grep("^b_original_", colnames(d), value = TRUE)
stopifnot(length(rad) == 99)
sub <- as.data.frame(d[, rad])
colnames(sub) <- sub("^b_", "", colnames(sub))
set.seed(42)
sub <- sub[sample(seq_len(nrow(sub)), 1000), ]
write.csv(sub, "manuscript/horseshoe_posterior_draws.csv", row.names = FALSE)
cat("Wrote", nrow(sub), "draws x", ncol(sub), "features\n")
