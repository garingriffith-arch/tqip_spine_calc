# Code, data, and model availability

The manuscript-facing SPINE-TRACT analysis code is organized in `analytic_pipeline/` with numbered entrypoints and original versioned implementation sources. These scripts reproduce analyses **given licensed ACS data and the requisite locally generated intermediate artifacts**; a clean clone does not contain or reconstruct patient-level ACS PUF data. The pipeline manifest and SHA-256 source manifest identify the exact archived scripts.

The original development involved controlled revisions including mandatory correction of diagnosis-derived phenotype variables and targeted duration/pooled-AUROC postprocessing. Scripts before the final repaired-diagnosis policy are not presented as alternative final models.

The final figure-builder revision and final B=1000 TRIPOD+AI script will be added after verifying the executed final local versions. Final full-2020–2024 deployment models **have not yet been fitted or supplied**. The repository-root Shiny calculator and existing `data/model_bundle.rds`/`data/predictor_metadata.rds` remain legacy artifacts, not the final 54/57-predictor models. See `analytic_pipeline/04_deployment/README.md` for the planned model inventory and artifact locations.

No patient-level PUFs, derived records, patient prediction tables, record IDs or licensed extracts are redistributed. Centre-level heterogeneity could not be assessed because facility keys were unavailable.
