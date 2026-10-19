# Final SPINE-TRACT deployment models

The completed Stage 21 builder is represented by `21_fit_deployment_models.R`, with a corrected CV-parser source under `../implementation/21_build_SPINE_TRACT_deployment_2020_2024_CVFIX_v2.R`. **Source provenance limit:** the deployment artifact manifest records `source_file_md5=64cc76848eeeef02fec044f7345ee464`; the available corrected builder copy does not have that MD5. The source is included as an executable reference, not represented as a verified byte-identical copy of the exact file used for the deployed fit.

The deployment procedure retains each endpoint's repaired-diagnosis predictor policy and the 2024-selected structural configuration, uses three-fold **internal CV only to select boosting rounds**, then fits eligible 2020–2024 data. This is not independent evaluation or an update to the manuscript's rolling-origin performance metrics.

## Resulting model artifacts

`../../data/deployment_2020_2024/models/` contains **eight XGBoost JSON models**, and `../../data/deployment_2020_2024/encoders/` contains the matching endpoint-level RDS preprocessing objects:

- Five classifiers: `discharge3_final`, `hlos3_final`, `icu3_final`, `vent3_final`, `operative_spine_intervention`.
- Three joint Q10/Q50/Q90 boosters: `hospital_los`, `icu_los_conditional`, `ventilator_days_conditional`.

The duration models each produce three quantiles; they are not nine independent fitted models. Two classifier outcomes (discharge and categorical HLOS) have 57 inputs; the remaining classifiers and durations have 54.

Run `Rscript analytic_pipeline/verify_deployment_artifacts.R` from the repository root to compare the JSON/RDS objects to their SHA-256 manifest. The full 2020–2024 deployment builder also records model-serialization prediction parity on its local QC examples (maximum reported difference approximately 1.2e−7).

## Reproducing deployment fitting

Set `SPINE_TRACT_PROJECT_DIR` to the authorized licensed local workspace containing the repaired frozen cohort, predictor lock, the selected final classification configurations (`14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv`) and duration configurations (`04_SELECTED_DURATION_CONFIGS.csv`). Install the packages required in the Stage 21 source. Duration fitting requires Python XGBoost 2+ and uses `reg:quantileerror`.

From the repository root run the numbered entrypoint with `SPINE_TRACT_DEPLOY_STAGE=classification` and then `SPINE_TRACT_DEPLOY_STAGE=duration`. Both stages fit **new deployment models**; they are computationally intensive and are not needed merely to use or redeploy the already-fitted JSON files. Original reported internal CV and serialized output checks are recorded in `artifact_manifest.csv`; patient-level input and QC caches are not public.

## Application

The app loads the eight models from `data/deployment_2020_2024/`. R XGBoost handles classifiers, while Python XGBoost 3.2.0 handles joint-quantile inference through `reticulate`. See [deployment instructions](../../README_DEPLOYMENT.md).

All-years fitting, internal CV for selecting rounds and successful Shiny deployment do **not** provide independent external or prospective model validation.
