# Malignant-transformation risk in lower-grade glioma — an IDH-stratified Bayesian interval-censored analysis

Statistical analysis code for the manuscript *"Malignant-transformation risk in
lower-grade glioma: an IDH-stratified Bayesian interval-censored analysis"*
(submitted to npj Precision Oncology).

The study analyses 155 adult-type WHO grade 2–3 LGG patients (77
malignant-transformation events; 6 non-MT deaths; median imaging surveillance
among non-transformers 4.8 years, IQR 2.0–10.0). Because no single available
estimator simultaneously handles the four structural features of these data —
interval-censored event detection, a competing risk of non-MT death,
time-varying treatment exposures, and ~100 candidate covariates — the contrast-
enhancement (CE) hazard ratio is required to be consistent across four
**complementary survival models**:

1. **Primary** — Bayesian interval-censored Weibull with a horseshoe prior
   over 102 coefficients (3 clinical + 99 PyRadiomics).
2. **Secondary** — Fine-Gray subdistribution-hazards model (midpoint plug-in
   and 20-fold Rubin-pooled multiple imputation).
3. **Secondary** — Cox counting-process time-varying-covariate model with
   20-fold pooled multiple imputation of the interval-censored event time.
4. **Secondary** — Fine-Gray fit with native interval-censored sieve maximum
   likelihood (`intccr`; Park & Bakoyannis 2021), integrating the
   subdistribution-hazards likelihood directly over the bracket [L, R].

Pre-specified subgroups: IDH-mutant (n = 109, 48 events) and IDH-wildtype
(n = 42, 26 events).

## What this repository contains

```
.
├── README.md                          # this file
├── LICENSE.md                         # MIT
├── .gitignore                         # excludes patient data, .rds, build cruft
├── docs/
│   └── STATISTICAL_ANALYSIS_PLAN.md   # pre-specified analysis plan
├── data/                              # de-identified data are NOT included here
│   ├── README.md                      #   (available on request — see "Data availability")
│   └── ttm_baseline_Definition_DE.*   # data dictionary (column definitions)
├── scripts/                           # analysis pipeline (R + Python)
│   ├── 00_setup.R                     # install Stan + dependencies
│   ├── 01_load_data.R                 # load + verify datasets
│   ├── 02_bayes_baseline.R            # primary Bayesian Weibull (horseshoe)
│   ├── 03_tvc_cox_mi.R                # Cox counting-process TVC + multiple imputation
│   ├── 04_competing_risks.R           # Fine-Gray (midpoint) + Gray's test
│   ├── 05_diagnostics.R               # convergence + posterior tables
│   ├── 06_summary.R                   # summary tables
│   ├── 07_sensitivity.R               # horseshoe par_ratio + adult-only + era split
│   ├── 08_eor_binary_sensitivity.R    # binary GTR vs non-GTR sensitivity
│   ├── 08_bootstrap_optimism_cox_tvc.py  # internal validation (Harrell optimism)
│   ├── 10_finegray_mi.R               # Fine-Gray with 20-fold MI
│   ├── 11_intccr_finegray_ic.R        # native interval-censored Fine-Gray (intccr)
│   ├── extract_horseshoe_draws.R      # helper for figure rendering
│   └── era_split/                     # era-split summary outputs
├── figures/                           # final manuscript & supplement figures
│   ├── Figure_1 … Figure_5            # rendered (PNG/JPEG)
│   ├── Supplementary_Figure_S1, S2
│   ├── vector_pdf/                    # vector PDFs (Fig 2, Fig 4, Fig S2)
│   └── code/                          # figure-generation scripts
│       ├── make_figs_corrected.py     # Fig 2 (four-model concordance forest) + Fig 4
│       ├── make_analyticflow_v2.py    # Fig S2 (four-model analytic-flow diagram)
│       ├── make_swimmer_full_cohort.py# Fig 1 (swimmer plot)
│       └── make_horseshoe_ridgeline.py# Fig 3 (per-feature horseshoe ridgeline)
└── submission/                        # final npj Precision Oncology submission package
    ├── 01_Cover_Letter.docx
    ├── 02_Manuscript.docx
    ├── 03_Plain_Language_Summary.docx
    ├── 04_Supplementary_Information.docx
    ├── 05_Supplementary_Table_S6.docx
    ├── figures/                       # figures as submitted
    ├── reporting_checklists/          # STROBE (S7), TRIPOD+AI (S8), Reporting Summary
    └── results_tables/                # numerical result tables (CSV)
```

The folder `source_data/Data_LGG_Clinical_update.xlsx` (raw clinical extract)
and any `final_summary/`, `*.rds`, `*.log` files produced by a script run are
excluded from version control via `.gitignore`.

## Methods overview

| Component | Method | Script |
|---|---|---|
| **Primary inference (baseline → MT)** | Bayesian Weibull AFT (interval-censored) with a horseshoe prior (Piironen–Vehtari) on all 102 population-level coefficients (3 clinical + 99 radiomic). `par_ratio = 0.1`, `df = 1`, `df_global = 1`, `scale_global = 2`, `scale_slab = 2`, `df_slab = 4`. 4 chains × 8 000 iterations × 4 000 warmup, `adapt_delta = 0.999`, `max_treedepth = 14`. Three cohorts: full / IDH-mutant / IDH-wildtype. | `02_bayes_baseline.R` |
| **Time-varying treatment effects** | Frequentist Cox PH on multiply-imputed event times (M = 20 uniform draws within each interval-censored observation), pooled via Rubin's rules with patient-level robust SE clustering. Specification: counting-process with step-function TVCs at EOR / chemo / radiation epochs. Cohorts: full and IDH-mutant. | `03_tvc_cox_mi.R` |
| **Competing risks — midpoint plug-in** | Fine-Gray subdistribution hazards with the midpoint of [L, R] as the event time. Two specifications: clinical-only and clinical + 2-feature radiomic signature. Cumulative incidence by CE strata via Gray's test. | `04_competing_risks.R` |
| **Competing risks — multiple imputation** | Same Fine-Gray model refit with M = 20 multiple imputations of the event time from Uniform(L, R), pooled by Rubin's rules (seeds matched to the Cox-TVC pipeline). | `10_finegray_mi.R` |
| **Competing risks — native interval-censored** | Fine-Gray with native interval-censored sieve maximum likelihood and B-spline approximation of the cumulative incidence functions, via `intccr::ciregic()` (Park & Bakoyannis 2021), integrating the subdistribution-hazards likelihood directly over the bracket [L, R] rather than relying on a midpoint or imputation plug-in. | `11_intccr_finegray_ic.R` |
| **Internal validation** | Harrell's bootstrap-optimism procedure on the primary Cox-TVC counting-process fit (500 patient-level resamples; single mid-interval imputation per replicate). Returns apparent and optimism-corrected C-index and a bias-corrected CE log-hazard ratio. | `08_bootstrap_optimism_cox_tvc.py` |
| **Diagnostics** | R-hat, ESS_bulk / ESS_tail, divergent transitions, max-treedepth hits per fit; clinical-coefficient posteriors with HR derived as `HR = exp(−α · β)` (AFT → PH conversion via joint α/β draws); per-feature radiomic posterior summaries with 95 %-CrI exclusion of the null. | `05_diagnostics.R` |
| **Headline tables** | Consolidated CSV tables per cohort and analysis. | `06_summary.R` |
| **Sensitivity (pre-specified)** | Horseshoe `par_ratio` grid {0.05, 0.10, 0.20}; age ≥ 18 sensitivity (n = 149); era-stratified analysis (diagnosis year dichotomised at 2010). | `07_sensitivity.R` |
| **Sensitivity (EOR)** | Binary EOR collapse: GTR vs non-GTR (STR + Biopsy pooled), full and IDH-mutant cohorts. Linear hypothesis test via pooled covariance. | `08_eor_binary_sensitivity.R` |

The hazard ratios reported throughout the Bayesian analysis use the AFT → PH
conversion `HR = exp(−α · β)` (with α = Weibull shape) computed draw-by-draw
from `as_draws_df(fit)`, propagating joint α / β posterior uncertainty.
`brms::weibull()` is AFT-parameterised; a naive `exp(β)` reading would invert
direction and rescale magnitude.

## Headline results

| Estimand | Estimator | Value (95 % interval) |
|---|---|---|
| CE on MT, full cohort | Bayesian Weibull | HR 22.1 (CrI 10.5–53.2) |
| CE on MT, full cohort | Fine-Gray (midpoint) | sHR 10.5 (CI 5.6–19.5) |
| CE on MT, full cohort | Cox-TVC (Rubin-pooled MI) | HR 13.3 (CI 6.8–26.0) |
| CE on MT, full cohort | Fine-Gray (native IC, intccr) | sHR 2.94 (CI 1.27–6.79) |
| CE on MT, IDH-mutant | Bayesian Weibull | HR 92.2 (CrI 26.9–416) |
| GTR vs Biopsy, full | Cox-TVC | HR 0.32 (CI 0.13–0.77) |
| GTR vs Biopsy, IDH-mutant | Cox-TVC | HR 0.16 (CI 0.05–0.53) |
| Discrimination (apparent → corrected) | Harrell C-index | 0.847 → 0.83 |
| Radiomic features with CrI excluding null | Bayesian horseshoe | 0 / 99 |

See the result tables in `submission/results_tables/` and the manuscript's Supplementary Tables (`submission/04_Supplementary_Information.docx`) for the full numerical output.

## Pre-specified hypotheses

See `docs/STATISTICAL_ANALYSIS_PLAN.md` for the full plan. The pre-specified
contrasts evaluated are:

- **H1**: 95 % CrI for the CE hazard ratio excludes 1 in the full cohort.
- **H2**: At least one radiomic feature has 95 % CrI excluding 1 in the full
  cohort.
- **H3**: The radiomic signal preserves direction in IDH-mutant with at least
  one feature's CrI excluding 1.
- **H4**: The post-surgery time-varying effect HR has CrI overlapping 1, or is
  in the protective direction.
- **H5**: Extent of resection (GTR contrast in the Cox-TVC) is associated with
  reduced post-surgical hazard.

## Dependencies

- **R ≥ 4.3** (developed under R 4.5.3)
- **Python ≥ 3.10** (for `08_bootstrap_optimism_cox_tvc.py`; lifelines + numpy)
- **CmdStan ≥ 2.36** (developed under 2.38.0; installed via
  `cmdstanr::install_cmdstan()` in `00_setup.R`)
- **C++ toolchain** — on macOS, run `xcode-select --install` once if Stan
  compilation fails
- **R packages** (auto-installed by `00_setup.R`):
  `dplyr`, `readr`, `ggplot2`, `survival`, `posterior`, `bayesplot`,
  `tidybayes`, `loo`, `cmdstanr`, `brms`, `cmprsk`, `intccr`, `mice`, `readxl`
- **Python packages**: `lifelines`, `numpy`, `pandas`, `scikit-survival`

Tested on Apple Silicon (M4, 16 GB), R 4.5.3, brms 2.23.0, CmdStan 2.38.0,
Python 3.11. Total wall-clock for the full pipeline: ~120–180 min including
the one-time CmdStan compile.

## Reproducing the analysis

Run from the project root, in order. Each script reports what it produces and
what's next; do not skip ahead unless you already have the prerequisite RDS
files.

```bash
cd scripts/

# One-time: install Stan toolchain (~15-40 min)
Rscript 00_setup.R

# Per-run pipeline:
Rscript 01_load_data.R                         # <1 min
Rscript 02_bayes_baseline.R                    # 15-25 min  (3 brms fits)
Rscript 03_tvc_cox_mi.R                        # 10-20 min  (Cox-TVC + MI)
Rscript 04_competing_risks.R                   # 5-10 min   (Fine-Gray midpoint)
Rscript 05_diagnostics.R                       # 2-5 min    (convergence + tables)
Rscript 06_summary.R                           # <1 min
Rscript 07_sensitivity.R                       # 30-60 min  (5 brms refits)
Rscript 08_eor_binary_sensitivity.R            # 5-10 min   (binary EOR)
python3 08_bootstrap_optimism_cox_tvc.py       # 10-20 min  (Harrell bootstrap)
Rscript 10_finegray_mi.R                       # 5-10 min   (Fine-Gray MI)
Rscript 11_intccr_finegray_ic.R                # 20-40 min  (native IC FG)
```

`04_competing_risks.R`, `10_finegray_mi.R`, and `07_sensitivity.R` (era arm)
require the raw clinical extract `Data_LGG_Clinical_update.xlsx` to compute
death-without-MT and diagnosis era. Expected location:
`source_data/Data_LGG_Clinical_update.xlsx` relative to the project root, or
set the env variable `LGG_CLINICAL_XLSX` to an absolute path. **This file is
NOT included in the repository — it contains identifiable clinical data.**

A successful run produces:
- `scripts/01_*.rds` (prepared datasets)
- `scripts/02_fit_{full,mut,wt}.rds` (Bayesian fits, ~13 MB each)
- `scripts/03_pooled_tvc_*.rds`, `03_tvc_summary.rds`
- `scripts/04_finegray_results.rds`, `04_cuminc_by_ce.rds`
- `scripts/05_clinical_effects.{rds,csv}`, `05_radiomic_posterior.{rds,csv}`
- `scripts/06_summary.R` populates `scripts/final_summary/table*.csv`
- `scripts/07_sens_horseshoe_pr{0.05,0.1,0.2}.rds`, `07_sens_{adult,era}.rds`
- `scripts/08_eor_binary_results.rds`
- `scripts/08_bootstrap_optimism_{summary,draws}.csv`
- `scripts/10_finegray_mi_results.rds`
- `scripts/11_intccr_results.rds`

## Data dictionary

### `data/ttm_baseline.csv` — 155 rows, one per patient

| Column | Type | Description |
|---|---|---|
| `patient_number` | int | Anonymised patient ID |
| `L` | numeric | Lower bound of MT interval, years from diagnosis |
| `R` | numeric | Upper bound of MT interval (`Inf` if right-censored) |
| `event` | 0 / 1 | 1 = transformation observed, 0 = censored |
| `case` | Case_0/2/3 | Case classification: 0 = no MT, 2 = pre-surgical MT, 3 = post-surgical MT |
| `EOR` | factor | Extent of resection: GTR / STR / Biopsy |
| `idh` | 0 / 1 | IDH mutation status (1 = mutant) |
| `p19q` | 0 / 1 | 1p/19q codeletion status |
| `CDKAN2A/B` | text | CDKN2A/B alteration status |
| `age_atDiag45` | numeric | Age at diagnosis, centred at 45 years |
| `vol_diag_centered` | numeric | Tumour volume at diagnosis, cm³, median-centred |
| `contrast_enhancingI` | "N" / "Y" | Contrast enhancement at baseline imaging |

A German-language data dictionary (`data/ttm_baseline_Definition_DE.md`)
provides the same information with extended provenance notes on the variables
sourced from `Data_LGG_Clinical_update.xlsx`.

### `data/counting_process.csv` — 371 rows in (start, stop, status) format

Per-patient row sequences with epoch boundaries at first surgery, first
chemotherapy, first radiation.

| Column | Type | Description |
|---|---|---|
| `patient_number` | int | Anonymised patient ID |
| `tstart`, `tstop` | numeric | Epoch interval, years from diagnosis |
| `post_surgery`, `post_chemo`, `post_radiation` | 0 / 1 | Step-function indicators (current epoch) |
| `eor_state` | factor | `pre` / `GTR` / `STR` / `Biopsy` |
| `L_orig`, `R_orig` | numeric | Original MT interval (constant per patient) |
| `EOR`, `idh`, `p19q`, `CDKN2A_B`, `age_atDiag45`, `vol_diag_centered`, `contrast_enhancingI` | various | Constant per patient |

### `data/event_timing.csv` — 155 rows

Per-patient first-treatment timestamps (years from diagnosis); `NA` if
untreated.

### `data/radiomics_99_zscore.csv` and `data/radiomics_428_zscore.csv`

PyRadiomics-derived features extracted from FLAIR / T1 / T1c / T2 sequences,
z-scored within feature column. The 99-feature panel is a pre-specified
ICC-stable subset; the 428-feature panel is the full PyRadiomics output
retained for sensitivity. Column-naming convention:

```
original_<family>_<descriptor>_<sequence>
```

e.g. `original_glszm_largearealowgraylevelemphasis_T1c`. Families: `shape`,
`firstorder`, `glcm`, `gldm`, `glrlm`, `glszm`, `ngtdm`. Sequences: `FLAIR`,
`T1`, `T1c`, `T2`. The pipeline uses the 99-feature panel by default.

## Known methodological notes

The Methods and Supplementary Information of the manuscript cover these in
detail. Headline items:

- **Prior specification**: brms 2.23 does not allow per-coefficient priors
  alongside a class-wide horseshoe; the horseshoe is therefore applied to all
  102 coefficients (clinical + radiomic). The wide slab (`scale_slab = 2`)
  preserves strong clinical effects; weak effects receive additional
  shrinkage. Prior-sensitivity refit at `scale_global = 1` confirmed headline
  inferences are unchanged.
- **Sampler divergence**: 0.07–0.42 % of post-warmup transitions diverge
  across the three primary fits, after the standard adapt_delta / iter /
  scale_global retry ladder. A divergence-location diagnostic confirms
  divergences concentrate in horseshoe auxiliary-scale hyperparameters
  (`sdb_*`); clinical and radiomic coefficient posteriors are essentially
  divergence-free.
- **AFT → PH conversion**: `brms::weibull()` is AFT-parameterised; the
  proportional-hazards HR is `exp(−α · β)`, computed draw-by-draw to
  propagate α/β joint uncertainty.
- **`post_surgery` × `eor_state` collinearity in 03**: an exact linear
  dependency (`post_surgery ≡ dummy_GTR + dummy_STR + dummy_Biopsy`) makes
  the three-level EOR contrasts relative to Biopsy rather than pre.
  `08_eor_binary_sensitivity.R` provides a clean reparameterisation that
  avoids this dependency.
- **Anticipatory regressor removed**: a previously included "time-to-first-
  surgery as baseline regressor" specification (Cox-TVC "Model B") has been
  removed from the report because it uses a regressor whose value is only
  known after the outcome of interest is partially observed
  (cf. statistical-review note); the effect was null (HR 1.03, p = 0.53) and
  removal does not affect any substantive conclusion.
- **Native interval-censored Fine-Gray**: `intccr::ciregic()` with `alpha =
  c(0, 0)` (proportional subdistribution hazards for both event types) and
  B-spline sieve approximation of the cumulative-incidence functions. Cause 2
  (death without MT, n = 5 with known time-of-death) is too sparse for stable
  joint estimates and is not reported as a separate cause-specific model.
- **No TVC + radiomics joint Bayesian model**: Bayesian time-varying-
  covariate survival with interval censoring is not well-supported by current
  R/Stan tooling; TVC analysis uses frequentist Cox + MI, and the primary
  Bayesian analysis uses baseline covariates only.

## Reproducibility

- All random seeds set to 2026 (or `2026 + m` for per-imputation draws in
  scripts 03, 08, and 10).
- Scripts read from `01_*.rds` (produced by `01_load_data.R`); upstream RDS
  files are deterministic given fixed seeds.
- brms / CmdStan posteriors depend on Stan's compiler: minor numerical
  differences across CmdStan versions are expected at the 4th–5th decimal but
  headline HRs and CrIs are stable.
- The submission version of this repository is tagged `v1.0-npj-submission`.

## Citation

If you use this code, please cite the accompanying manuscript (citation will
be added on acceptance) and the methodological references it relies on:

- Piironen, J. & Vehtari, A. (2017). Sparsity information and regularization
  in the horseshoe and other shrinkage priors. *Electron. J. Statist.* 11(2):
  5018–5051.
- Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models
  Using Stan. *J. Stat. Softw.* 80(1): 1–28.
- Stan Development Team. CmdStan: the command-line interface to Stan.
  https://mc-stan.org/users/interfaces/cmdstan
- Fine, J. P. & Gray, R. J. (1999). A proportional hazards model for the
  subdistribution of a competing risk. *J. Am. Stat. Assoc.* 94(446):
  496–509.
- Park, J. & Bakoyannis, G. (2021). Semiparametric regression on
  cumulative-incidence functions with interval-censored competing-risks
  data. *Stat. Med.* 40(15): 3417–3431.
- Harrell, F. E. Jr., Lee, K. L. & Mark, D. B. (1996). Multivariable
  prognostic models: issues in developing models, evaluating assumptions and
  adequacy, and measuring and reducing errors. *Stat. Med.* 15(4): 361–387.

## License

Code: MIT (see `LICENSE.md`).

Data: code-only repository. No patient-level data is included.

## Data availability

The de-identified analysis-ready datasets (`ttm_baseline.csv`,
`counting_process.csv`, `event_timing.csv`, `radiomics_99_zscore.csv`,
`radiomics_428_zscore.csv`, plus the derived `death_info.csv` used by
`04_competing_risks.R` / `10_finegray_mi.R` / `11_intccr_finegray_ic.R`) and
the original clinical extract (`Data_LGG_Clinical_update.xlsx`) are
**available from the corresponding author on reasonable request**, subject to
institutional and ethical-committee approval.

To reproduce the analysis, place the received CSVs in a local `data/`
directory at the project root (and the Excel, if provided, in `source_data/`).
The pipeline scripts will then run as documented under *Reproducing the
analysis*.

## Contact

Open an issue on this repository or contact the corresponding author of the
manuscript.
</content>
</invoke>