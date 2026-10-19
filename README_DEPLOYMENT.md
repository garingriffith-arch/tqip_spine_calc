# SPINE-TRACT — TBI-TRACT-style Shiny application

## Contents

- `app.R`: Shiny entry point.
- `app_core.R`: source-locked prediction, controls, UI, and validation.
- `test_app_contract.R`: eight-model local inference smoke test.
- `data/deployment_2020_2024/models/`: eight native model JSONs.
- `data/deployment_2020_2024/encoders/`: eight endpoint-specific RDS encoders.
- `data/deployment_2020_2024/artifact_manifest.csv`: model provenance and SHA-256 hashes.

Preserve the existing `www/ohsu_logo.png` in the GitHub repository. Without it the UI displays 'OHSU' text.

## Local RStudio testing and Connect Cloud manifest

From a **fresh** R session at this Git clone's root:

```r
install.packages(c("shiny", "bslib", "data.table", "Matrix",
                   "xgboost", "reticulate", "rsconnect", "jsonlite"))
# Restart R after installing packages if needed.
setwd("C:/Users/garin/Downloads/tqip_spine_calc")
Sys.setenv(RETICULATE_PYTHON =
  reticulate::conda_python("tbi-tract-xgb-gpu"))
source("test_app_contract.R")
shiny::runApp(".", launch.browser = TRUE)
```

To produce the deployable R/Python dependency manifest (without including
machine-specific Python paths), run:

```r
source("prepare_connect_cloud_manifest.R")
```

The helper explicitly includes only the 8 deployment JSONs, 8 encoder
RDS files, model artifact manifest, Shiny sources, the logo and Python
requirements. It verifies that `manifest.json` contains **reticulate**
and a **Python** version with `requirements.txt`. Commit and push only
the generated `manifest.json` (and reviewed deployment source changes).

The three joint-quantile models require Python XGBoost **3.2.0**, pinned in
`requirements.txt`. For Connect Cloud the app relies on host-provided
`RETICULATE_PYTHON` (or reticulate auto-discovery); it never uses
`reticulate::use_condaenv()` at runtime. Do not commit local Conda paths,
`.Renviron` or environment directories.

## Deploy from GitHub to Posit Connect Cloud

1. Open <https://connect.posit.cloud> and select **Publish** > **Shiny**.
2. Choose `garingriffith-arch/tqip_spine_calc`, branch
   `spine-tract-2020-2024-app`, and primary file `app.R`.
3. Make sure the committed `manifest.json` includes R `reticulate` plus
   a Python environment and the repo contains `requirements.txt`.
4. Set automatic republishing as preferred and begin publishing.
5. Check both R and Python environment restoration in the build logs.
   Open the app and test the full model set, including all three
   joint-quantile duration outputs, before merging into `main`.

Connect Cloud supports mixed R/Python content when both dependency
descriptions are present. If the Python environment fails to build, read
the specific build log before altering the models or retraining.
Clinical predictions have **not** been independently externally validated.

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
