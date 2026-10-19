# SPINE-TRACT 2020–2024 deployment artifacts — integration staging

This feature branch stages a replacement Shiny calculator using the repaired-DX
54-predictor clinical and 57-predictor contextual policies. The TBI-TRACT-like
application entrypoint lives at `app.R`, and model loading / UI / prediction
logic lives at `app_core.R`.

## Required runtime artifacts (local/private only; never committed to public GitHub)

Eight XGBoost model JSON files under `data/deployment_2020_2024/models/`
and eight matching encoder RDS files under
`data/deployment_2020_2024/encoders/`. Keep
`data/deployment_2020_2024/artifact_manifest.csv` with them.

Five classifiers use R XGBoost. The three duration models are **joint**
Q10/Q50/Q90 Python XGBoost boosters; run them with Python XGBoost >=2
via `reticulate`. Do not substitute the original legacy
`data/model_bundle.rds` or `data/predictor_metadata.rds`.

## Safety and publishing

- The source repository is **public**. Final JSON models, RDS encoders and
  artifact manifests are permanently excluded from public Git commits.
  Provision them privately to the local application or authorized hosting runtime.
- Git ignore applies only to Git; it does not prevent separately bundled model
  files from being copied to an application host or otherwise served.
- Never commit patient-level ACS PUFs, row-wise predictions, record
  identifiers, or derived patient-level data.
- The main branch and legacy app are intentionally unchanged while
  this integration is tested on a feature branch.
- Preserve `www/ohsu_logo.png` already tracked in this repository.
- Run `source("test_app_contract.R")` locally from the repo root
  before merging. Python XGBoost is a hosting requirement.
- Internal rolling-origin temporal evaluation is not independent external
  validation. These are all-2020–2024 deployment fits.
