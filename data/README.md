# Data directory

The analysis-ready CSVs that the pipeline reads from this directory are not
included in the repository. They contain de-identified patient-level data and
are released only on a controlled basis.

**Available upon reasonable request from the corresponding author**, subject
to institutional and ethical-committee approval:

- `ttm_baseline.csv`
- `counting_process.csv`
- `event_timing.csv`
- `radiomics_99_zscore.csv`
- `radiomics_428_zscore.csv`
- `death_info.csv` (derived; used by `04_competing_risks.R`)

Place the received files directly in this directory, then run the pipeline as
described in the project root `README.md`.

Column-level documentation for each file is in the *Data dictionary* section
of the project root `README.md`.
