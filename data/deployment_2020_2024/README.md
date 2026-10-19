# SPINE-TRACT 2020–2024 deployment artifacts — integration staging

This feature branch stages a replacement Shiny calculator using the repaired-DX
54-predictor clinical and 57-predictor contextual policies. The TBI-TRACT-like
application entrypoint lives at `app.R`, and model loading / UI / prediction
logic lives at `app_core.R`.

## Git-tracked final deployment artifacts

Eight XGBoost model JSON files under `data/deployment_2020_2024/models/`
and eight matching encoder RDS files under
`data/deployment_2020_2024/encoders/`. Keep
`data/deployment_2020_2024/artifact_manifest.csv` with them.

Five classifiers use R XGBoost. The three duration models are **joint**
Q10/Q50/Q90 Python XGBoost boosters; run them with Python XGBoost >=2
via `reticulate`. Do not substitute the original legacy
`data/model_bundle.rds` or `data/predictor_metadata.rds`.

## Safety and publishing

- This repository is **public**. Final model JSONs, encoder RDSs and the
  aggregate artifact manifest are explicitly allowed in Git under
  `data/deployment_2020_2024/` to mirror TBI-TRACT.
- Verify ACS-derived model redistribution and institutional authorization
  **before committing or pushing**. Once published, model files can remain
  accessible through Git history even if removed later.
- Never commit source ACS PUF files, row-wise predictions, patient or
  encounter identifiers, patient-level derivatives, or training/QC extracts.
- Use Git from a local clone for larger JSON objects; GitHub's browser
  uploader may reject those files.
- Preserve `www/ohsu_logo.png` and keep legacy models separate.
- All-years deployment fitting does not constitute independent validation.
- Python XGBoost >=2 must also be supported in the hosting environment
  for the three joint-quantile duration models.
