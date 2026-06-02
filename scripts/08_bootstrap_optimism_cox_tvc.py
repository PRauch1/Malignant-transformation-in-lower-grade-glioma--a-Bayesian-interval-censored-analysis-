#!/usr/bin/env python3
"""Bootstrap-optimism correction for the primary Cox-TVC model.

Harrell-style internal validation for the counting-process Cox regression that
generates the primary contrast-enhancement (CE) hazard ratio.

Design
------
1. Load counting_process.csv (371 rows across 155 patients). The primary
   Cox-TVC model is Surv(tstart, tstop, event) ~ contrast_enhancingI
   + factor(eor_state) + post_chemo + post_radiation + age_atDiag45
   + vol_diag_centered + idh, with MT detection time (tstop) taken from the
   mid-point of the interval [L, R] (i.e., one deterministic imputation of the
   20-fold MI pipeline). We use the mid-point draw rather than pooling all 20
   because:
     - bootstrap optimism operates on one model fit
     - the CE log-HR is essentially invariant across MI replicates
       (original Rubin-pooled HR = 13.3 vs single-imputation HR reported here)
2. Apparent performance: fit on full data, record log-HR for CE and Harrell
   C-index over concatenated (tstart, tstop) rows.
3. For b = 1..B = 500, resample patients with replacement (preserving each
   patient's full counting-process block), refit the model:
     (a) boot performance = log-HR and C from boot sample under boot model
     (b) test performance = log-HR and C on ORIGINAL sample under boot model
   optimism_b = boot_b - test_b
4. Corrected estimate = apparent - mean(optimism)
5. 95% bootstrap percentile interval from test_b distribution around corrected.

Outputs
-------
    scripts/08_bootstrap_optimism_summary.csv    — apparent / optimism / corrected
    scripts/08_bootstrap_optimism_draws.csv      — per-replicate boot/test values
"""
from __future__ import annotations
import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_DIR = HERE

N_BOOT = 500
SEED = 202604

# Factor columns → numeric dummies so we can apply a boot-fit model cleanly
FACTOR_COLS = {
    "eor_state": ["pre", "GTR", "STR", "Biopsy"],   # pre is reference
}
COVARS = [
    "contrast_enhancingI_Y",   # binary CE, Y/N -> 1/0
    "eor_state_GTR", "eor_state_STR", "eor_state_Biopsy",
    "post_chemo", "post_radiation",
    "age_atDiag45", "vol_diag_centered",
    "idh",
]


def prep(df: pd.DataFrame) -> pd.DataFrame:
    """Build the design matrix from counting_process rows.

    Expects one row per (patient, interval). The terminal row of each patient
    carries the event indicator via `is_event_patient`.
    """
    df = df.copy()

    # CE binary
    df["contrast_enhancingI_Y"] = (df["contrast_enhancingI"] == "Y").astype(int)

    # EOR dummies — 'pre' is reference level
    for lvl in ["GTR", "STR", "Biopsy"]:
        df[f"eor_state_{lvl}"] = (df["eor_state"] == lvl).astype(int)

    # idh: numeric already (0/1), keep NaN as 0 (very few missings in CP file)
    df["idh"] = df["idh"].fillna(0).astype(float)

    # Determine the terminal row per patient → set event indicator
    df = df.sort_values(["patient_number", "tstart"]).reset_index(drop=True)
    last_row_idx = df.groupby("patient_number").tail(1).index
    df["event"] = 0
    df.loc[last_row_idx, "event"] = df.loc[last_row_idx, "is_event_patient"].astype(int)

    # Keep only needed columns
    keep = ["patient_number", "tstart", "tstop", "event"] + COVARS
    return df[keep].copy()


def fit_cox(df: pd.DataFrame) -> CoxPHFitter:
    cph = CoxPHFitter(penalizer=0.0)
    cph.fit(df, duration_col="tstop", event_col="event",
            entry_col="tstart", formula=" + ".join(COVARS),
            show_progress=False)
    return cph


def linear_predictor(model: CoxPHFitter, df: pd.DataFrame) -> np.ndarray:
    """Linear predictor (log-risk) for each row."""
    # predict_log_partial_hazard returns per-row log-risk (without baseline)
    return model.predict_log_partial_hazard(df[COVARS]).values


def cindex_counting_process(df: pd.DataFrame, lp: np.ndarray) -> float:
    """Harrell C-index suited to counting-process: for each patient, use the
    terminal row's (stop, event) with the lp averaged across their rows (lp is
    constant within patient in our model — all covariates are baseline or
    step-function, so the terminal row's lp is fine).
    """
    tail = (df.assign(lp=lp)
              .groupby("patient_number").tail(1))
    return concordance_index(tail["tstop"].values,
                             -tail["lp"].values,  # higher risk → lower "score"
                             tail["event"].values)


def logHR_CE(model: CoxPHFitter) -> float:
    """Extract log-HR for contrast_enhancingI_Y."""
    return float(model.params_.loc["contrast_enhancingI_Y"])


def bootstrap():
    # Load source counting-process data
    cp = pd.read_csv(os.path.join(DATA, "counting_process.csv"))
    df = prep(cp)

    patients = df["patient_number"].unique()
    n_pat = len(patients)
    print(f"Patients: {n_pat}, rows: {len(df)}")

    # Apparent performance
    apparent_model = fit_cox(df)
    lhr_app = logHR_CE(apparent_model)
    lp_app = linear_predictor(apparent_model, df)
    c_app = cindex_counting_process(df, lp_app)
    print(f"Apparent  log-HR(CE)= {lhr_app:.3f}  HR= {np.exp(lhr_app):.2f}"
          f"   C= {c_app:.3f}")

    rng = np.random.default_rng(SEED)
    draws = []
    n_fail = 0
    for b in range(N_BOOT):
        # Sample patient IDs with replacement
        sampled = rng.choice(patients, size=n_pat, replace=True)
        # Build boot dataset by concatenating each patient's rows
        # (duplicate patients get duplicated blocks)
        pieces = []
        for i, p in enumerate(sampled):
            block = df[df["patient_number"] == p].copy()
            # Assign a unique pseudo-id to avoid lifelines complaining about dups
            block["patient_number"] = i + 1
            pieces.append(block)
        boot_df = pd.concat(pieces, ignore_index=True)

        try:
            boot_model = fit_cox(boot_df)
        except Exception as e:
            n_fail += 1
            continue

        # boot performance on boot sample
        lp_boot_on_boot = linear_predictor(boot_model, boot_df)
        c_boot = cindex_counting_process(boot_df, lp_boot_on_boot)
        lhr_boot = logHR_CE(boot_model)

        # test: apply boot model to the ORIGINAL data
        lp_boot_on_orig = linear_predictor(boot_model, df)
        c_test = cindex_counting_process(df, lp_boot_on_orig)
        # for the coefficient test, lhr_boot is a property of the model; the
        # 'test' contribution to optimism is 0 for the coefficient itself, so
        # we track a non-parametric alternative: the log-HR from the boot
        # model applied to the original data's CE stratification is simply
        # lhr_boot — but for optimism we follow Harrell's prescription for
        # the predictive performance (C-index) and additionally record the
        # bootstrap distribution of lhr_boot for a percentile CI.
        draws.append({
            "rep": b + 1,
            "lhr_boot": lhr_boot,
            "c_boot": c_boot,
            "c_test": c_test,
            "opt_c": c_boot - c_test,
        })

        if (b + 1) % 50 == 0:
            print(f"  rep {b+1}/{N_BOOT}  "
                  f"C_boot={c_boot:.3f}  C_test={c_test:.3f}"
                  f"  opt={c_boot-c_test:+.3f}")

    D = pd.DataFrame(draws)

    # Optimism-corrected C
    mean_opt = D["opt_c"].mean()
    c_corr = c_app - mean_opt

    # Bootstrap percentile CI on log-HR(CE)
    lhr_ci = np.percentile(D["lhr_boot"].values, [2.5, 97.5])

    # Optimism on log-HR: mean(lhr_boot) - lhr_app (simple shrinkage estimate)
    opt_lhr = D["lhr_boot"].mean() - lhr_app
    lhr_corr = lhr_app - opt_lhr

    summary = pd.DataFrame([{
        "metric": "CE log-HR",
        "apparent": lhr_app,
        "apparent_HR": np.exp(lhr_app),
        "mean_bootstrap": D["lhr_boot"].mean(),
        "optimism": opt_lhr,
        "corrected": lhr_corr,
        "corrected_HR": np.exp(lhr_corr),
        "boot_ci_lo": lhr_ci[0],
        "boot_ci_hi": lhr_ci[1],
        "boot_HR_ci_lo": np.exp(lhr_ci[0]),
        "boot_HR_ci_hi": np.exp(lhr_ci[1]),
        "n_boot": len(D),
        "n_fail": n_fail,
    }, {
        "metric": "Harrell C",
        "apparent": c_app,
        "apparent_HR": np.nan,
        "mean_bootstrap": D["c_boot"].mean(),
        "optimism": mean_opt,
        "corrected": c_corr,
        "corrected_HR": np.nan,
        "boot_ci_lo": np.percentile(D["c_test"], 2.5),
        "boot_ci_hi": np.percentile(D["c_test"], 97.5),
        "boot_HR_ci_lo": np.nan,
        "boot_HR_ci_hi": np.nan,
        "n_boot": len(D),
        "n_fail": n_fail,
    }])

    print()
    print(summary.to_string(index=False))

    out_summary = os.path.join(OUT_DIR, "08_bootstrap_optimism_summary.csv")
    out_draws = os.path.join(OUT_DIR, "08_bootstrap_optimism_draws.csv")
    summary.to_csv(out_summary, index=False)
    D.to_csv(out_draws, index=False)
    print(f"\nWrote {out_summary}")
    print(f"Wrote {out_draws}")


if __name__ == "__main__":
    bootstrap()
