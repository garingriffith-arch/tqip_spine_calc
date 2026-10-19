# Model artifacts and deployment status

The repository-root `app.R` currently reads `data/model_bundle.rds` and `data/predictor_metadata.rds`. These objects belong to a **legacy, superseded** SPINE-TRACT model generation. They must not be presented as the revised 2020–2024 manuscript models or used to interpret the published temporal performance estimates.

The final revised deployment fit has not been completed. It must produce 14 XGBoost model objects: five classifiers (four multiclass trajectories and recorded thoracolumbar operation) and nine duration-quantile models (Q10/Q50/Q90 for each of three duration targets), with endpoint-specific encoders and provenance. The intended destination is `data/deployment_2020_2024/` after training and inference validation. See [deployment status](../analytic_pipeline/04_deployment/README.md).

Do not commit ACS patient-level data, individual predictions, record identifiers, or restricted derivatives to this repository.
