# SPINE-TRACT — TBI-TRACT-style Shiny application

## Contents

- `app.R`: Shiny entry point.
- `app_core.R`: source-locked prediction, controls, UI, and validation.
- `test_app_contract.R`: eight-model local inference smoke test.
- `data/deployment_2020_2024/models/`: eight native model JSONs.
- `data/deployment_2020_2024/encoders/`: eight endpoint-specific RDS encoders.
- `data/deployment_2020_2024/artifact_manifest.csv`: model provenance and SHA-256 hashes.

Preserve the existing `www/ohsu_logo.png` in the GitHub repository. Without it the UI displays 'OHSU' text.

## Local RStudio launch

From a fresh R session in this app directory:

```r
install.packages(c("shiny", "bslib", "data.table", "Matrix", "xgboost", "reticulate"))
Sys.setenv(SPINE_TRACT_PY_ENV = "tbi-tract-xgb-gpu")
source("test_app_contract.R")
shiny::runApp(".", launch.browser = TRUE)
```

**Python dependency**: the three multi-quantile models require Python XGBoost >=2 (the models were saved by XGBoost 3.2.0). The app uses `reticulate` to load Python XGBoost. Configure `SPINE_TRACT_PY_ENV` to an EXISTING conda environment, or omit the variable if `reticulate` already selects an appropriate Python. The application never attempts to install an environment at startup.

For instance, in PowerShell:

```powershell
$env:SPINE_TRACT_PY_ENV = "tbi-tract-xgb-gpu"
```

Create the environment separately if necessary; do not assume an existing server has the workstation's conda environment.

**Deployment hosting**: confirm a stable Python runtime with XGBoost >=2 and R reticulate integration. On Posit Connect, configure Python as described in the host's documentation. A plain R-only shinyapps.io deployment has not been verified with these Python-backed quantile models.

## Numerical contract

- Five R XGBoost classifiers: `Matrix::Matrix(..., sparse=TRUE)` with native numeric missing routing.
- Three joint-quantile Python boosters: dense `double` feature matrices; raw log1p predictions are `expm1` transformed, floored to the recorded minimum, then sorted Q10 <= Q50 <= Q90.
- Model feature order, encoder class levels, training years and predictor policy come directly from the Stage 21 encoder RDS objects.
- Context (`race_clean`, `ethnicity_clean`, `payer_clean`) appears **only** in the discharge and categorized hospital-LOS models.
- Quick preview treats hidden inputs as missing/unknown and was not separately evaluated.
- Reported 2024 temporal metrics belong to the 2020–2023 trained evaluation models, **not** to an external validation of these all-2020–2024 deployment models.

**Not tested here**: R Shiny application startup, R 1-row inference, renderUI state preservation and Python/R reticulate interoperability. Run `test_app_contract.R` before deployment and compare predictions against the Stage 21 training environment.

## GitHub and deployment artifacts

This project uses the same Git-tracked deployment structure as TBI-TRACT.
The **eight reviewed** 2020–2024 deployment model JSONs, **eight**
endpoint-specific encoder RDSs and artifact manifest live under
`data/deployment_2020_2024/`. The project's `.gitignore` explicitly permits
these deployment files while excluding other RDS/CSV and patient-level data.

The repository is **public**, so the model and encoder files become publicly
downloadable (including through commit history). Check applicable ACS TQP PUF
license and institutional terms **before committing**. Exclude source PUFs,
derived patient-level records, row-wise predictions, patient identifiers, raw
outcome data, QC folders and temporary training caches. Do not commit any
historical temporal-evaluation model in place of the all-years deployment fit.

A working hosted Shiny app also needs R XGBoost plus Python XGBoost >=2 under
`reticulate`. Git-based code/model distribution alone does not guarantee
that the hosting runtime has the Python environment for the three joint-quantile
models. Smoke-test from a clean checkout and verify that the hosting platform
supports Python before production launch.

## Features mirrored from TBI-TRACT

Complete / Quick mode, accordion clinical inputs, cached values across mode switches, three-column probability cards and risk bars, continuous quantile cards, model performance/context/disclaimer section, responsive OHSU-style visual language. Endpoint-specific differences: thoracolumbar phenotypes, HLOS <=7 / 8–19 / >=20, recorded operation rather than neurosurgical interventions, missing imaging stability/TLICS/AIS data, and the need for Python joint-quantile inference.
