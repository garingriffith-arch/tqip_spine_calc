# SPINE-TRACT final deployment artifacts

Eight models trained on eligible ACS TQP 2020–2024 records using the repaired-DX predictor lock:

| Model family | Files | Predictors |
| --- | ---: | ---: |
| Discharge disposition and categorical hospital LOS classifiers | 2 | 57 |
| ICU, ventilation, and recorded thoracolumbar operation classifiers | 3 | 54 |
| Hospital LOS, conditional ICU LOS and conditional ventilation duration joint-quantile models | 3 | 54 |

The joint quantile models each output Q10, Q50 and Q90. They are three fitted boosters, not nine separate model objects.

- `models/`: eight native XGBoost JSON files.
- `encoders/`: eight RDS preprocessing metadata files, including feature order and categorical levels.
- `artifact_manifest.csv`: source fit configuration, boosting rounds, eligible training counts, serialization parity, model sizes and SHA-256 checksums.

Run `Rscript analytic_pipeline/verify_deployment_artifacts.R` from the repository root. Deployment fitting source and assumptions are documented in [Stage 21](../../analytic_pipeline/04_deployment/README.md). Five classifiers use R XGBoost; three joint-quantile duration models use Python XGBoost.

The trained objects are **all-years deployment fits**. They are not the models used to generate the 2022–2024 rolling-origin internal temporal performance estimates.

This repository is public. The model and encoder files are downloadable. The source ACS PUFs, derived patient-level datasets, prediction caches and record identifiers are not provided; data and model redistribution remain subject to applicable terms.
