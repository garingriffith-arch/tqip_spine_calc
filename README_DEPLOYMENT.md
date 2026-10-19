# SPINE-TRACT deployment

**Live app:** https://01a154f9-e4de-f916-8a37-2663ed8ed9d1.share.connect.posit.cloud

The app is an R Shiny application published through GitHub to Posit Connect Cloud. It loads five R XGBoost classifiers and three Python XGBoost joint-quantile duration models. No licensed ACS source records are required for inference.

## Deployable files

- `app.R`, `app_core.R`: Shiny application and prediction code.
- `data/deployment_2020_2024/models/`: eight native XGBoost JSON model files.
- `data/deployment_2020_2024/encoders/`: eight matching endpoint-specific RDS encoder files.
- `data/deployment_2020_2024/artifact_manifest.csv`: model metadata, provenance and SHA-256 checksums.
- `manifest.json`: R and Python environment metadata for Connect Cloud.
- `requirements.txt`: Python XGBoost 3.2.0.
- `www/ohsu_logo.png`: application logo.

The R manifest specifies the **official CRAN XGBoost 3.2.1.1** for inference. Do not regenerate it in the locally customized GPU-enabled R XGBoost 1.7.2.1 environment: that package has no corresponding CRAN archive tarball and previously caused an HTTP 404 during deployment.

## Check the checkout

From the project root:

```sh
Rscript analytic_pipeline/verify_archive.R
Rscript analytic_pipeline/verify_deployment_artifacts.R
```

Then, with RStudio configured to use Python XGBoost 3.2.0 through `reticulate`, run:

```r
source("test_app_contract.R")
shiny::runApp(".", launch.browser = TRUE)
```

The smoke test checks class probability shape/sums, quantile ordering and per-endpoint feature-column order. It does not establish clinical validity. Separate R and Python dependencies are needed for app startup.

## Connect Cloud

The application has been published from the `spine-tract-2020-2024-app` integration branch. After merging into `main`, update the Connect Cloud source branch if you want future deployments to track `main`. The GitHub repo alone does not verify the active deployed commit.

To republish, select the GitHub repository, correct branch and `app.R`; review R/Python environment restore and the runtime log. If a dependency changes, regenerate the manifest in a **clean R package library** with the official CRAN XGBoost version installed. The helper `prepare_connect_cloud_manifest.R` rejects the custom local GPU build and validates Python requirements.

The app uses dense numeric feature matrices for Python joint-quantile inference and sparse matrices for R classification; both preserve the training encoders' column order, numeric missing-value routing and categorical levels. Duration results use `expm1`, truncate to the outcome's feasible minimum and sort Q10/Q50/Q90.

## Interpretation and restrictions

This tool predicts registry-recorded inpatient trajectories and **recorded operation receipt**, not whether an operation is indicated. It is not a substitute for imaging-based spine stability assessment or neurological examination. The 2022–2024 estimates reflect internal temporal evaluation; external and prospective validation have not been completed.

The deployed models and encoders are publicly downloadable from this repository. The underlying ACS PUFs, derived patient-level records, individual predictions, patient/encounter IDs and caches are not included. Redistribution must conform to applicable ACS/institutional terms.
