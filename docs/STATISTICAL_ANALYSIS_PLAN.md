# Statistical Analysis Plan: Bayesian Reanalysis of TTM in Lower-Grade Glioma

## Background and rationale

The original analysis (manuscript NOA-D-26-00176) used `icenReg` parametric proportional hazards with forward stepwise AIC selection. Reviewer critiques and post-rejection re-examination identified four substantive issues:

1. **Endpoint mismatch**: models estimated time to malignancy-or-first-surgery (29 events) while the paper claimed time to malignant transformation (77 events).
2. **Stepwise selection instability**: forward AIC selection with 99-428 candidate features and 77 events produces unreliable feature sets and biased effect sizes.
3. **No internal validation**: apparent AICs were reported without optimism correction.
4. **Pre/postoperative treatment effects ignored**: ~60% of the cohort received adjuvant chemo or radiation, but these were treated as if the patients existed in a single biological state.

This revised analysis addresses all four.

## Study population and endpoint

- N=155 phenotypic LGG patients, single-center cohort, 1999-2023
- Cohort construction unchanged from prior analysis (Hauer thesis section 4.1)
- **Primary endpoint**: time from diagnosis to malignant transformation, defined per RANO-style imaging criteria or histopathological upgrading
- **Censoring structure**: interval-censored. Lower bound = date of last MRI without malignant features (or date of first surgery for post-surgical events). Upper bound = date of first malignant MRI
- **Event count**: 77 transformations (Case 2: 27 pre-surgery, Case 3: 50 post-surgery), 78 right-censored (Case 0)

## Primary analysis: Bayesian Weibull PH with regularized horseshoe

### Model specification

For patient i, the hazard of malignant transformation at time t:
```
h_i(t) = h_0(t; α, λ) · exp(β_clinical · X_i,clinical + β_radiomic · X_i,radiomic)
```
where h_0 is a Weibull baseline hazard with shape α and scale λ.

### Priors

**Baseline distribution**:
- α ~ Gamma(2, 1) (shape parameter, weakly favors increasing or decreasing hazards)
- λ ~ Half-Normal(0, 5) (scale)

**Clinical coefficients** (weakly informative, on log-HR scale):
- β_contrast_enhancement ~ Student-t(3, 0, 2.5)
- β_age (per year, centered at 45) ~ Student-t(3, 0, 0.1)
- β_volume (per cm³, centered at median) ~ Student-t(3, 0, 0.05)

**Radiomic coefficients** (regularized horseshoe):
- β_radiomic,j ~ N(0, τ · λ_j)
- τ ~ HalfStudent-t(1, 0, scale_global=1)
- λ_j ~ HalfStudent-t(1, 0, 1) with regularization toward slab
- Slab: scale_slab = 2, df_slab = 4
- Expected fraction of nonzero coefficients: par_ratio = 0.1 (≈10/99)

### Cohorts

1. Full cohort (n=155, 77 events) — primary
2. IDH-mutant subgroup (n=109, 48 events) — pre-specified subgroup analysis
3. IDH-wildtype subgroup (n=42, 26 events) — exploratory; flag if convergence fails

### Inference

- 4 chains × 4000 iterations × 2000 warmup
- adapt_delta = 0.995, max_treedepth = 14
- Convergence criteria: max R-hat ≤ 1.01, min ESS_bulk ≥ 400, zero divergent transitions
- Reported quantities: posterior median HR, 95% credible interval, posterior probability HR > 1.2 or HR < 1/1.2 (practically meaningful effect)

## Secondary analysis: time-varying treatment effects

### Rationale

Surgery, adjuvant chemo, and radiation all occur during follow-up and modify the biological state of the tumor. A single-baseline-state model cannot answer whether post-operative residual tumor behaves differently after GTR vs STR vs biopsy, or whether adjuvant therapy modifies subsequent transformation hazard.

### Method

**Counting-process format**: each patient contributes 2-4 rows corresponding to biological epochs:
- Epoch 1: diagnosis → first surgery (pre-operative state)
- Epoch 2+: post-surgical, with EOR state revealed
- Subsequent epochs: post-chemo, post-radiation as applicable

**Time-varying covariates**:
- `post_surgery` (0/1, step function at date of first surgery)
- `eor_state` (factor: pre, GTR, STR, Biopsy; takes the patient's eventual EOR value after surgery)
- `post_chemo` (0/1, step function at date of first chemo)
- `post_radiation` (0/1, step function at date of first radiation)

**Extended model A: time-since-treatment effects**: adds `t_since_surgery`, `t_since_chemo`, `t_since_radiation` as continuous covariates that increase from 0 post-treatment, to test whether the treatment-associated hazard evolves with time since the procedure (non-proportional treatment hazards).

**Extended model B: time-to-treatment as baseline covariate**: adds per-patient `t_to_surgery_at_baseline` and indicators `has_chemo`, `has_radiation` to test whether patients with longer treatment delays represent a clinically different subpopulation (selection effect on the basis of clinical stability at diagnosis).

**Estimation**: Cox proportional hazards on multiply-imputed event times (M=20). Within each interval-censored observation, draw event time from Uniform(L, R). Fit standard Cox per imputation. Pool coefficients via Rubin's rules.

This is a pragmatic compromise — a fully Bayesian interval-censored TVC model is not well-supported in current software. The MI Cox approach gives unbiased estimation of the time-varying effects with proper uncertainty propagation.

## Tertiary analysis: competing-risks consideration

### Question

Does the radiomic signature predict cumulative incidence of malignant transformation when accounting for the competing risk of death without transformation?

### Method

Fine-Gray subdistribution hazards model on multiply-imputed event times. Two specifications:
- Clinical predictors only (CE + age + volume)
- Clinical + 2-feature radiomic signature (GLSZM-T1c + GLDM-FLAIR)

Death without transformation coded as competing event 2; transformation as primary event 1.

## Sensitivity analyses

1. **Era as scanner-generation proxy**: dichotomize diagnosis date at 2010, include as covariate in primary model. Confirm signature survives era adjustment.
2. **Alternative baseline distribution**: refit primary model with Log-logistic baseline. Check that posterior medians for key effects are within 20% of Weibull-baseline estimates.
3. **Pediatric exclusion**: rerun primary analysis on adult-only cohort (n=149, age ≥ 18). Confirm key inferences hold.

## Pre-specified hypotheses and decisions

- **H1**: Bayesian posterior 95% CrI for contrast enhancement HR excludes 1 in full cohort. (Strong prior expectation; failure would indicate modeling error.)
- **H2**: At least one radiomic feature has 95% CrI for HR excluding 1 in full cohort. (Test of whether radiomics adds anything.)
- **H3**: The radiomic signal preserves direction in IDH-mutant subgroup with at least one feature CrI excluding 1. (Test of whether signal is independent of IDH-wildtype confounding.)
- **H4**: post_surgery time-varying effect HR has CrI overlapping 1, OR if it differs, in the protective direction (HR < 1, surgery delays transformation). (Test of whether surgery itself modifies trajectory.)
- **H5**: time-to-first-surgery as a baseline covariate has HR < 1 (longer delay to surgery = lower subsequent transformation hazard). This would be evidence for clinical-stability selection — patients observed long enough to be operated later were selected on the basis of already-favorable disease. Confirmation would strengthen, not weaken, the paper's argument that baseline radiomics captures intrinsic biology rather than clinical triage.

## What the published paper will report

If H1, H2, H3 hold and convergence is clean:
- Bayesian posterior estimates as primary results
- Frequentist parametric PH (icenReg, prior approach) as supplement showing analytic robustness
- Time-varying treatment effects as key novel finding
- Fine-Gray as honest competing-risks accounting
- Honest acknowledgment of single-center limitation

If H2 or H3 fail:
- The paper becomes a more modest negative-finding contribution
- Reframe as "baseline radiomics, while statistically associated with timing, contributes only marginal incremental prognostic information beyond contrast enhancement and IDH status in this cohort"
- Still publishable in a methods-focused venue but probably not JNO

## Authorship implication

Helga Wagner, as statistical co-author and original methods supervisor, should review and approve this analysis plan before the resulting analyses are incorporated into the manuscript revision.
