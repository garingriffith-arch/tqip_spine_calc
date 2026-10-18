# Code, data, and model availability

SPINE-TRACT is a prediction-model research project using the American College of Surgeons (ACS) Trauma Quality Programs Participant Use Files (2020–2024).

## Code

The repository now separates the **manuscript-facing rolling temporal study** from the **legacy Shiny calculator**. The current reviewer-facing table builder is archived in `analytic_pipeline/implementation/` and has a numbered reporting entry point.

The frozen 647-code non-patient-level operative-treatment map and the Stage 5 fixed-prediction metric-postprocessing source (together with its pooled-AUROC numeric correction) are also archived under `analytic_pipeline/`.\n\nThe remaining finalized model-development, duration-quantile, sensitivity, and figure scripts have **not yet been archived in this GitHub repository**. They should only be added after their final versions are reconciled with the locked study artifacts. The absence of a file must not be interpreted as evidence that an analysis was not performed; it means the public code archive is incomplete. See `analytic_pipeline/PIPELINE_MANIFEST.md`.

## Participant-level data

No ACS patient-level PUF files, derived patient datasets, prediction caches, or record-level identifiers should be uploaded. The ACS source files are subject to access and data-use restrictions and require independent authorized access.

Aggregate reporting tables and non-patient-level phenotype/procedure code maps may be published after review. The table builder requires local frozen analytic data and aggregate output artifacts; it is **not** executable from a clean clone alone.

## Model objects and calculator

The repository predates the final manuscript model revision and contains `app.R`, `data/model_bundle.rds`, and `data/predictor_metadata.rds` from an **earlier model generation**. These objects must not be represented as the final 54/57-predictor, three-class and duration-quantile SPINE-TRACT models.

The final 2020–2024 deployment model objects have not been archived here and independent external/prospective validation remains pending. The current app has not been migrated or tested against the final architecture.

## Model-reporting and reproducibility status

The manuscript uses rolling-origin internal temporal evaluations for admission years 2022–2024. The 2024 evaluation informed some development choices and must not be described as untouched independent or external validation.

Full reproduction requires the correct ACS source PUFs, the finalized harmonization and model scripts, appropriate computing dependencies, the frozen predictor dictionary, and the evaluated model/prediction artifacts. Reporting code is preserved with its original environment assumptions for provenance, rather than presented as a self-contained push-button reproduction.

Last documentation update: October 2026.
