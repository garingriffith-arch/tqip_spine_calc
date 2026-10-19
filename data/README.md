# SPINE-TRACT application model artifacts

The active `app.R` loads the **final repaired-diagnosis 2020–2024 deployment models** in `deployment_2020_2024/`. The superseded legacy `model_bundle.rds` and `predictor_metadata.rds` have been removed from the current repository tree.

`deployment_2020_2024/models/`: five classification model JSONs and three joint-quantile duration model JSONs.

`deployment_2020_2024/encoders/`: matching endpoint-specific predictor/encoder metadata RDS files. The training feature order and categorical unknown/unseen-level policy are part of each encoder.

`deployment_2020_2024/artifact_manifest.csv`: training counts, years, boosting rounds, internal-CV settings, serialized prediction parity, file sizes and SHA-256 checksums.

Verify these files from the repository root:

```sh
Rscript analytic_pipeline/verify_deployment_artifacts.R
```

The PUF source data, patient-level derived datasets, row-wise predictions, encounter IDs and QC caches are excluded. The publicly tracked deployment objects must not be mistaken for patient-level ACS records. Model code and redistribution are subject to applicable ACS/institutional terms.

See [deployment instructions](../README_DEPLOYMENT.md) and the [analysis pipeline](../analytic_pipeline/README.md).
